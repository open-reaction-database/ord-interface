/**
 * Copyright 2026 Open Reaction Database Project Authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

import React, { useEffect, useState } from 'react';
import { useSearchParams } from 'react-router-dom';
import SearchResults from '../search/SearchResults';
import { useNLQuery } from '../../hooks/useNLQuery';
import type { NLQueryData } from '../../hooks/useNLQuery';
import type { NLInterpretation, ResolvedComponent } from '../../types/search';
import './MainNLSearch.scss';

const EXAMPLES = [
  'reactions for synthesizing ibuprofen',
  'reactions using benzene as an input with yield greater than 70%',
  'reactions that make an aryl boronic acid',
];

const ROLE_LABEL = { OUTPUT: 'product', INPUT: 'reactant/reagent' } as const;

// A "(verbatim)" resolver tag means the model supplied the structure directly (a SMILES
// or SMARTS), so no name resolver was actually called.
const isVerbatim = (resolver: string): boolean => resolver.includes('verbatim');

// The distinct network resolvers that actually ran, ignoring verbatim passthroughs and
// the "(cached)" suffix so a fresh and a cached PubChem hit count as one resolver.
const resolversUsed = (components: ResolvedComponent[]): string[] =>
  Array.from(
    new Set(
      components
        .filter(component => !isVerbatim(component.resolver))
        .map(component => component.resolver.replace(/ \(cached\)$/, '')),
    ),
  );

/** Extracts the model's non-component filters as human-readable strings. */
const filterSummary = (interpretation: NLInterpretation): string[] => {
  const filters: string[] = [];
  if (interpretation.min_yield != null)
    filters.push(`yield ≥ ${interpretation.min_yield}%`);
  if (interpretation.max_yield != null)
    filters.push(`yield ≤ ${interpretation.max_yield}%`);
  if (interpretation.min_conversion != null)
    filters.push(`conversion ≥ ${interpretation.min_conversion}%`);
  if (interpretation.max_conversion != null)
    filters.push(`conversion ≤ ${interpretation.max_conversion}%`);
  if (interpretation.reaction_smarts)
    filters.push(`reaction SMARTS ${interpretation.reaction_smarts}`);
  if (interpretation.similarity_threshold != null)
    filters.push(`similarity threshold ${interpretation.similarity_threshold}`);
  if (interpretation.use_stereochemistry) filters.push('stereochemistry respected');
  if (interpretation.limit != null) filters.push(`limit ${interpretation.limit}`);
  return filters;
};

/**
 * Shows how a question was read, separating the two provenance layers: what the
 * language model extracted (the `build_query` tool call) and what the deterministic
 * resolvers turned each identifier into. A "verbatim" resolver means the model supplied
 * the structure directly (a SMILES or SMARTS), so nothing was looked up.
 */
const Interpretation: React.FC<{ data: NLQueryData }> = ({ data }) => {
  const { interpretation, resolvedComponents } = data;
  const filters = filterSummary(interpretation);
  // Only components whose names actually went through a network resolver; SMILES/SMARTS
  // the model supplied verbatim were never resolved and don't belong in this section.
  const resolved = resolvedComponents.filter(
    component => !isVerbatim(component.resolver),
  );

  return (
    <div className="nl-search__interpretation">
      <section className="nl-search__layer">
        <div className="nl-search__interpretation-title">
          From the model{' '}
          <span className="nl-search__provenance">(build_query tool call)</span>
        </div>
        {interpretation.components.length > 0 ? (
          <ul className="nl-search__interpretation-list">
            {interpretation.components.map((component, index) => (
              <li key={index}>
                <span className="nl-search__role">{ROLE_LABEL[component.target]}</span>{' '}
                <strong>{component.identifier}</strong>{' '}
                <span className="nl-search__mode">
                  ({component.mode.toLowerCase()})
                </span>
              </li>
            ))}
          </ul>
        ) : (
          <div className="nl-search__muted">(no components extracted)</div>
        )}
        {filters.length > 0 && (
          <ul className="nl-search__interpretation-list">
            {filters.map((filter, index) => (
              <li key={`filter-${index}`}>{filter}</li>
            ))}
          </ul>
        )}
        <details className="nl-search__raw">
          <summary>Raw tool call</summary>
          <pre className="nl-search__raw-json">
            {JSON.stringify(interpretation, null, 2)}
          </pre>
        </details>
      </section>

      {resolved.length > 0 && (
        <section className="nl-search__layer">
          <div className="nl-search__interpretation-title">
            Resolved to structures{' '}
            <span className="nl-search__provenance">
              ({resolversUsed(resolved).join(', ')})
            </span>
          </div>
          <ul className="nl-search__resolution-list">
            {resolved.map((component, index) => (
              <li key={index}>
                <span className="nl-search__identifier">{component.identifier}</span>
                <span className="nl-search__arrow">→</span>
                <code>{component.smiles}</code>
                <span className="nl-search__resolver">via {component.resolver}</span>
              </li>
            ))}
          </ul>
        </section>
      )}
    </div>
  );
};

const MainNLSearch: React.FC = () => {
  // The submitted query lives in the URL (?q=…) so searches are shareable and
  // survive reloads, mirroring the structured search page.
  const [searchParams, setSearchParams] = useSearchParams();
  const submittedQuery = searchParams.get('q');
  // Dev mode: translate + resolve but don't run the search; kept in the URL so it's
  // shareable and survives reloads, like the query itself.
  const dryRun = searchParams.get('dry_run') === '1';
  const [input, setInput] = useState(submittedQuery ?? '');
  // Keep the text box in sync with the URL so browser back/forward navigation (which
  // changes ?q=) updates the visible query instead of leaving a stale value.
  useEffect(() => {
    setInput(submittedQuery ?? '');
  }, [submittedQuery]);

  const { data, isFetching, error } = useNLQuery(submittedQuery, true, dryRun);

  const apply = (query: string, dry: boolean) => {
    const next: Record<string, string> = {};
    if (query) next.q = query;
    if (dry) next.dry_run = '1';
    setSearchParams(next);
  };

  const submit = (value: string) => {
    const trimmed = value.trim();
    setInput(trimmed);
    apply(trimmed, dryRun);
  };

  return (
    <div className="nl-search">
      <div
        className="nl-search__banner"
        role="status"
      >
        🚧 This feature is in development — results may be incomplete or change.
      </div>
      <h1 className="nl-search__title">Ask about reactions</h1>
      <p className="nl-search__subtitle">
        Describe what you&apos;re looking for in plain language — compound names are
        resolved to structures and matched against the database.
      </p>

      <form
        className="nl-search__form"
        onSubmit={event => {
          event.preventDefault();
          submit(input);
        }}
      >
        <input
          className="nl-search__input"
          type="text"
          value={input}
          placeholder="e.g. reactions using benzene as an input with yield greater than 70%"
          onChange={event => setInput(event.target.value)}
        />
        <button
          className="nl-search__button"
          type="submit"
          disabled={isFetching}
        >
          {isFetching ? 'Searching…' : 'Search'}
        </button>
      </form>

      <label className="nl-search__dry-run-toggle">
        <input
          type="checkbox"
          checked={dryRun}
          onChange={event =>
            apply(submittedQuery ?? input.trim(), event.target.checked)
          }
        />
        Dry run — translate &amp; resolve only, don&apos;t run the search
      </label>

      <div className="nl-search__examples">
        {EXAMPLES.map(example => (
          <button
            key={example}
            type="button"
            className="nl-search__example"
            onClick={() => submit(example)}
          >
            {example}
          </button>
        ))}
      </div>

      {error && <div className="nl-search__error">{(error as Error).message}</div>}

      {data && submittedQuery && (
        <>
          <Interpretation data={data} />
          {data.dryRun ? (
            <div className="nl-search__dry-run">
              <div className="nl-search__dry-run-title">
                Dry run — search not executed
              </div>
            </div>
          ) : data.results.length > 0 ? (
            <SearchResults searchResults={data.results} />
          ) : (
            <div className="nl-search__empty">
              No reactions matched. Try relaxing a filter, or rephrase the query.
            </div>
          )}
        </>
      )}
    </div>
  );
};

export default MainNLSearch;
