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

import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { MemoryRouter } from 'react-router-dom';
import { describe, expect, it, vi } from 'vitest';
import SearchOptions from './SearchOptions';

const renderOptions = (search = '') => {
  const onSearchOptions = vi.fn();
  render(
    <MemoryRouter initialEntries={[`/search${search}`]}>
      <SearchOptions onSearchOptions={onSearchOptions} />
    </MemoryRouter>,
  );
  return { onSearchOptions };
};

// The text inputs in the components section, in reactants-then-products order.
const smilesInputs = (): HTMLInputElement[] =>
  screen
    .getAllByRole('textbox')
    .filter(input => input.closest('.reagent.options') !== null) as HTMLInputElement[];

const sectionIsOpen = (heading: string): boolean =>
  screen.queryByText(heading) !== null;

const search = (params: Record<string, string | string[]>): string => {
  const query = new URLSearchParams();
  for (const [key, value] of Object.entries(params)) {
    for (const item of Array.isArray(value) ? value : [value]) {
      query.append(key, item);
    }
  }
  return `?${query.toString()}`;
};

describe('SearchOptions', () => {
  it('starts with every optional section collapsed', () => {
    renderOptions();
    expect(sectionIsOpen('Reactants & Reagents')).toBe(false);
    expect(sectionIsOpen('Reaction IDs')).toBe(false);
    expect(sectionIsOpen('Dataset IDs')).toBe(false);
    expect(screen.getByLabelText('Result Limit')).toHaveValue(100);
  });

  it('opens a collapsed section on click', async () => {
    const user = userEvent.setup();
    renderOptions();

    await user.click(screen.getByText('Components'));

    expect(sectionIsOpen('Reactants & Reagents')).toBe(true);
  });

  describe('reading components from the URL', () => {
    it('accepts the JSON form', () => {
      renderOptions(
        search({
          component: JSON.stringify({
            pattern: 'CCO',
            target: 'input',
            mode: 'substructure',
          }),
        }),
      );

      expect(smilesInputs().map(input => input.value)).toEqual(['CCO']);
      expect(screen.getByText('substructure')).toHaveClass('selected');
    });

    // Previously shared URLs used "pattern;target;mode".
    it('accepts the legacy delimited form', () => {
      renderOptions(search({ component: 'CCO;input;exact' }));

      expect(smilesInputs().map(input => input.value)).toEqual(['CCO']);
      expect(screen.getByText('exact')).toHaveClass('selected');
    });

    // A SMARTS pattern can contain ";" itself, so the legacy parse has to split
    // from the right.
    it('keeps semicolons inside a legacy SMARTS pattern', () => {
      renderOptions(search({ component: '[C;H2];input;smarts' }));

      expect(smilesInputs().map(input => input.value)).toEqual(['[C;H2]']);
    });

    it('falls back to the legacy parse for JSON that is not a component', () => {
      renderOptions(search({ component: 'CCO;input;exact' }));
      expect(smilesInputs().map(input => input.value)).toEqual(['CCO']);
    });

    it('files inputs as reactants and everything else as products', () => {
      renderOptions(
        search({
          component: [
            JSON.stringify({ pattern: 'CCO', target: 'input', mode: 'exact' }),
            JSON.stringify({ pattern: 'CC=O', target: 'output', mode: 'exact' }),
          ],
        }),
      );

      expect(smilesInputs().map(input => input.value)).toEqual(['CCO', 'CC=O']);
    });

    it('defaults a component with no mode to exact', () => {
      renderOptions(
        search({ component: JSON.stringify({ pattern: 'CCO', target: 'input' }) }),
      );

      expect(screen.getByText('exact')).toHaveClass('selected');
    });

    it('reads the stereochemistry and similarity settings', () => {
      renderOptions(
        search({
          component: JSON.stringify({
            pattern: 'CCO',
            target: 'input',
            mode: 'similar',
          }),
          use_stereochemistry: 'true',
          similarity: '0.75',
        }),
      );

      expect(screen.getByLabelText('Use Stereochemistry')).toBeChecked();
      expect(screen.getByLabelText('Similarity Threshold')).toHaveValue('0.75');
    });

    // The threshold is padded to two decimals for display: 0.5 reads as "0.50".
    it('pads the similarity threshold for display', () => {
      renderOptions(
        search({
          component: JSON.stringify({
            pattern: 'CCO',
            target: 'input',
            mode: 'similar',
          }),
        }),
      );

      expect(screen.getByText('0.50')).toBeInTheDocument();
    });
  });

  describe('reading the other sections from the URL', () => {
    it('opens the datasets section for dataset ids and DOIs', () => {
      renderOptions(
        search({ dataset_id: ['ord_dataset-1', 'ord_dataset-2'], doi: '10.1000/x' }),
      );

      expect(screen.getByText('ord_dataset-1')).toBeInTheDocument();
      expect(screen.getByText('ord_dataset-2')).toBeInTheDocument();
      expect(screen.getByText('10.1000/x')).toBeInTheDocument();
    });

    it('opens the reactions section for reaction ids and SMARTS', () => {
      renderOptions(
        search({ reaction_id: 'ord-abc', reaction_smarts: '[C:1]>>[C:1]O' }),
      );

      expect(screen.getByText('ord-abc')).toBeInTheDocument();
      expect(screen.getByText('[C:1]>>[C:1]O')).toBeInTheDocument();
    });

    it('opens the reactions section for a narrowed yield range', () => {
      renderOptions(search({ min_yield: '20', max_yield: '80' }));

      expect(sectionIsOpen('Reaction IDs')).toBe(true);
      expect(screen.getByText('20% - 80%')).toBeInTheDocument();
    });

    it('opens the reactions section for a narrowed conversion range', () => {
      renderOptions(search({ min_conversion: '10', max_conversion: '90' }));

      expect(sectionIsOpen('Reaction IDs')).toBe(true);
      expect(screen.getByText('10% - 90%')).toBeInTheDocument();
    });

    // Number(undefined) is NaN, and NaN !== 100, so an unguarded comparison
    // against the default bound opened this section on every page load.
    it('leaves the reactions section closed with no reaction params', () => {
      renderOptions(search({ dataset_id: 'ord_dataset-1' }));

      expect(sectionIsOpen('Reaction IDs')).toBe(false);
    });

    it('leaves the reactions section closed for a full-width range', () => {
      renderOptions(search({ max_yield: '100', max_conversion: '100' }));

      expect(sectionIsOpen('Reaction IDs')).toBe(false);
    });

    it('reads the result limit', () => {
      renderOptions(search({ limit: '25' }));
      expect(screen.getByLabelText('Result Limit')).toHaveValue(25);
    });
  });

  describe('emitting options', () => {
    it('reports the components with the shared match mode', async () => {
      const user = userEvent.setup();
      const { onSearchOptions } = renderOptions(
        search({
          component: [
            JSON.stringify({ pattern: 'CCO', target: 'input', mode: 'exact' }),
            JSON.stringify({ pattern: 'CC=O', target: 'output', mode: 'exact' }),
          ],
          similarity: '0.75',
          use_stereochemistry: 'true',
          limit: '25',
        }),
      );

      await user.click(screen.getByRole('button', { name: 'Search' }));

      expect(onSearchOptions).toHaveBeenCalledWith({
        reagent: {
          reagents: [
            { smileSmart: 'CCO', source: 'input', matchMode: 'exact' },
            { smileSmart: 'CC=O', source: 'output', matchMode: 'exact' },
          ],
          useStereochemistry: true,
          similarityThreshold: 0.75,
        },
        reaction: {
          reactionIds: [],
          reactionSmarts: [],
          min_yield: 0,
          max_yield: 100,
          min_conversion: 0,
          max_conversion: 100,
        },
        dataset: { datasetIds: [], DOIs: [] },
        general: { limit: 25 },
      });
    });

    it('applies a match-mode change to every component', async () => {
      const user = userEvent.setup();
      const { onSearchOptions } = renderOptions(
        search({
          component: JSON.stringify({
            pattern: 'CCO',
            target: 'input',
            mode: 'exact',
          }),
        }),
      );

      await user.click(screen.getByText('substructure'));
      await user.click(screen.getByRole('button', { name: 'Search' }));

      expect(onSearchOptions.mock.calls[0][0].reagent.reagents).toEqual([
        { smileSmart: 'CCO', source: 'input', matchMode: 'substructure' },
      ]);
    });

    it('reports components added and edited by hand', async () => {
      const user = userEvent.setup();
      const { onSearchOptions } = renderOptions();

      await user.click(screen.getByText('Components'));
      await user.click(screen.getAllByRole('button', { name: /Add Component/ })[0]);
      await user.type(smilesInputs()[0], 'CCO');
      await user.click(screen.getByRole('button', { name: 'Search' }));

      expect(onSearchOptions.mock.calls[0][0].reagent.reagents).toEqual([
        { smileSmart: 'CCO', source: 'input', matchMode: 'exact' },
      ]);
    });

    it('reports a component removed by hand', async () => {
      const user = userEvent.setup();
      const { onSearchOptions } = renderOptions(
        search({
          component: [
            JSON.stringify({ pattern: 'CCO', target: 'input', mode: 'exact' }),
            JSON.stringify({ pattern: 'CCC', target: 'input', mode: 'exact' }),
          ],
        }),
      );

      const [firstDelete] = screen
        .getAllByRole('button')
        .filter(button => button.closest('.delete') !== null);
      await user.click(firstDelete);
      await user.click(screen.getByRole('button', { name: 'Search' }));

      expect(onSearchOptions.mock.calls[0][0].reagent.reagents).toEqual([
        { smileSmart: 'CCC', source: 'input', matchMode: 'exact' },
      ]);
    });

    it('reports the result limit', async () => {
      const user = userEvent.setup();
      const { onSearchOptions } = renderOptions();

      await user.clear(screen.getByLabelText('Result Limit'));
      await user.type(screen.getByLabelText('Result Limit'), '5');
      await user.click(screen.getByRole('button', { name: 'Search' }));

      expect(onSearchOptions.mock.calls[0][0].general).toEqual({ limit: 5 });
    });
  });

  describe('the structure editor', () => {
    // ModalKetcher embeds an iframe that never boots under jsdom; these tests
    // cover the wiring around it, not the editor itself.
    const drawButtons = () =>
      screen.getAllByRole('button').filter(button => button.closest('.draw') !== null);

    it('opens the editor on the reactant that was clicked', async () => {
      const user = userEvent.setup();
      renderOptions(
        search({
          component: [
            JSON.stringify({ pattern: 'CCO', target: 'input', mode: 'exact' }),
            JSON.stringify({ pattern: 'CCC', target: 'input', mode: 'exact' }),
          ],
        }),
      );

      await user.click(drawButtons()[1]);

      expect(screen.getByTitle('Ketcher Molecular Editor')).toBeInTheDocument();
    });

    it('opens the editor on a product', async () => {
      const user = userEvent.setup();
      renderOptions(
        search({
          component: JSON.stringify({
            pattern: 'CC=O',
            target: 'output',
            mode: 'exact',
          }),
        }),
      );

      await user.click(drawButtons()[0]);

      expect(screen.getByTitle('Ketcher Molecular Editor')).toBeInTheDocument();
    });

    it('closes the editor without changing anything', async () => {
      const user = userEvent.setup();
      renderOptions(
        search({
          component: JSON.stringify({
            pattern: 'CCO',
            target: 'input',
            mode: 'exact',
          }),
        }),
      );

      await user.click(drawButtons()[0]);
      await user.click(screen.getByRole('button', { name: 'Cancel' }));

      expect(screen.queryByTitle('Ketcher Molecular Editor')).not.toBeInTheDocument();
      expect(smilesInputs().map(input => input.value)).toEqual(['CCO']);
    });
  });
});
