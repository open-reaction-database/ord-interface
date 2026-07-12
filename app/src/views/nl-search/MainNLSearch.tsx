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
import { useLocation, useSearch } from 'wouter';
import {
  Alert,
  Button,
  Checkbox,
  Code,
  Flex,
  Paper,
  Text,
  TextInput,
  Title,
} from '@mantine/core';
import { IconBarrierBlock, IconSearch } from '@tabler/icons-react';
import SearchResults from '../search/SearchResults';
import { useNLQuery } from '../../hooks/useNLQuery';
import type { NLQueryData } from '../../hooks/useNLQuery';
import type { NLInterpretation, ResolvedComponent } from '../../types/search';
import classes from './MainNLSearch.module.scss';

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
    <Paper
      radius="sm"
      p="lg"
      className={classes.interpretation}
    >
      <section>
        <Text className={classes.layerTitle}>
          From the model{' '}
          <span className={classes.provenance}>(build_query tool call)</span>
        </Text>
        {interpretation.components.length > 0 ? (
          <ul className={classes.list}>
            {interpretation.components.map((component, index) => (
              <li key={index}>
                <span className={classes.role}>{ROLE_LABEL[component.target]}</span>{' '}
                <strong>{component.identifier}</strong>{' '}
                <span className={classes.mode}>({component.mode.toLowerCase()})</span>
              </li>
            ))}
          </ul>
        ) : (
          <Text className={classes.muted}>(no components extracted)</Text>
        )}
        {filters.length > 0 && (
          <ul className={classes.list}>
            {filters.map((filter, index) => (
              <li key={`filter-${index}`}>{filter}</li>
            ))}
          </ul>
        )}
        <details className={classes.raw}>
          <summary>Raw tool call</summary>
          <pre className={classes.rawJson}>
            {JSON.stringify(interpretation, null, 2)}
          </pre>
        </details>
      </section>

      {resolved.length > 0 && (
        <section className={classes.layer}>
          <Text className={classes.layerTitle}>
            Resolved to structures{' '}
            <span className={classes.provenance}>
              ({resolversUsed(resolved).join(', ')})
            </span>
          </Text>
          <ul className={classes.resolutionList}>
            {resolved.map((component, index) => (
              <li key={index}>
                <span>{component.identifier}</span>
                <span className={classes.arrow}>→</span>
                <Code>{component.smiles}</Code>
                <span className={classes.resolver}>via {component.resolver}</span>
              </li>
            ))}
          </ul>
        </section>
      )}
    </Paper>
  );
};

const MainNLSearch: React.FC = () => {
  // The submitted query lives in the URL (?q=…) so searches are shareable and
  // survive reloads, mirroring the structured search page.
  const search = useSearch();
  const [, navigate] = useLocation();
  const searchParams = new URLSearchParams(search);
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
    const next = new URLSearchParams();
    if (query) next.set('q', query);
    if (dry) next.set('dry_run', '1');
    navigate(`/ask?${next.toString()}`);
  };

  const submit = (value: string) => {
    const trimmed = value.trim();
    setInput(trimmed);
    apply(trimmed, dryRun);
  };

  return (
    <div className={classes.nlSearch}>
      <Alert
        className={classes.banner}
        icon={<IconBarrierBlock size={16} />}
        color="orange"
        radius="sm"
      >
        This feature is in development — results may be incomplete or change.
      </Alert>

      <Title order={1}>Ask about reactions</Title>
      <Text className={classes.subtitle}>
        Describe what you&apos;re looking for in plain language — compound names are
        resolved to structures and matched against the database.
      </Text>

      <form
        className={classes.form}
        onSubmit={event => {
          event.preventDefault();
          submit(input);
        }}
      >
        <TextInput
          className={classes.input}
          size="md"
          radius="sm"
          leftSection={<IconSearch size={18} />}
          value={input}
          placeholder="e.g. reactions using benzene as an input with yield greater than 70%"
          onChange={event => setInput(event.currentTarget.value)}
        />
        <Button
          type="submit"
          size="md"
          radius="sm"
          loading={isFetching}
        >
          Search
        </Button>
      </form>

      <Checkbox
        className={classes.dryRunToggle}
        size="sm"
        label="Dry run — translate & resolve only, don't run the search"
        checked={dryRun}
        onChange={event => apply(submittedQuery ?? input.trim(), event.target.checked)}
      />

      <Flex
        className={classes.examples}
        gap="xs"
        wrap="wrap"
      >
        {EXAMPLES.map(example => (
          <Button
            key={example}
            variant="default"
            size="compact-sm"
            radius="xl"
            className={classes.example}
            onClick={() => submit(example)}
          >
            {example}
          </Button>
        ))}
      </Flex>

      {error && (
        <Alert
          color="red"
          radius="sm"
          className={classes.error}
        >
          {(error as Error).message}
        </Alert>
      )}

      {data && submittedQuery && (
        <>
          <Interpretation data={data} />
          {data.dryRun ? (
            <Paper
              radius="sm"
              p="lg"
              className={classes.dryRunNote}
            >
              <Text fw={600}>Dry run — search not executed</Text>
            </Paper>
          ) : data.results.length > 0 ? (
            <SearchResults searchResults={data.results} />
          ) : (
            <Paper
              radius="sm"
              p="lg"
              className={classes.dryRunNote}
            >
              <Text c="secondary.1">
                No reactions matched. Try relaxing a filter, or rephrase the query.
              </Text>
            </Paper>
          )}
        </>
      )}
    </div>
  );
};

export default MainNLSearch;
