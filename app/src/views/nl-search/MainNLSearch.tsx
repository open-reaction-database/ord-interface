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

import React, { useState } from 'react';
import { useSearchParams } from 'react-router-dom';
import SearchResults from '../search/SearchResults';
import { useNLQuery } from '../../hooks/useNLQuery';
import type { NLQueryData } from '../../hooks/useNLQuery';
import './MainNLSearch.scss';

const EXAMPLES = [
  'reactions for synthesizing ibuprofen',
  'reactions using benzene as an input with yield greater than 70%',
  'reactions that make an aryl boronic acid',
];

/** Renders the model's interpretation so users can see how their question was read. */
const Interpretation: React.FC<{ data: NLQueryData }> = ({ data }) => {
  const { interpretation, resolvedComponents } = data;
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
  if (interpretation.use_stereochemistry) filters.push('stereochemistry respected');

  return (
    <div className="nl-search__interpretation">
      <div className="nl-search__interpretation-title">Interpreted as:</div>
      <ul className="nl-search__interpretation-list">
        {resolvedComponents.map((component, index) => (
          <li key={index}>
            <span className="nl-search__role">
              {component.target === 'OUTPUT' ? 'product' : 'reactant/reagent'}
            </span>{' '}
            <strong>{component.identifier}</strong>{' '}
            <span className="nl-search__mode">({component.mode.toLowerCase()})</span>{' '}
            <code>{component.smiles}</code>{' '}
            <span className="nl-search__resolver">via {component.resolver}</span>
          </li>
        ))}
        {filters.map((filter, index) => (
          <li key={`filter-${index}`}>{filter}</li>
        ))}
        {resolvedComponents.length === 0 && filters.length === 0 && (
          <li>(no constraints extracted)</li>
        )}
      </ul>
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
              <pre className="nl-search__dry-run-query">
                {JSON.stringify(
                  data.queryComponents.map(component => JSON.parse(component)),
                  null,
                  2,
                )}
              </pre>
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
