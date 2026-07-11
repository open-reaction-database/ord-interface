/**
 * Copyright 2023 Open Reaction Database Project Authors
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

import React, {
  Suspense,
  useCallback,
  useEffect,
  useMemo,
  useState,
  lazy,
} from 'react';
import { useSearch } from 'wouter';
import {
  Accordion,
  ActionIcon,
  Button,
  Checkbox,
  Flex,
  NumberInput,
  RangeSlider,
  SegmentedControl,
  Slider,
  Text,
  TextInput,
} from '@mantine/core';
import { IconPencil, IconPlus, IconSearch, IconTrash } from '@tabler/icons-react';
import SearchItemList from './SearchItemList';
import LoadingSpinner from '../../components/LoadingSpinner';
import classes from './SearchOptions.module.scss';

const ModalKetcher = lazy(() => import('../../components/ModalKetcher'));

interface ReagentComponent {
  smileSmart: string;
  source: string;
  matchMode?: string;
}

interface ParsedComponent {
  pattern: string;
  target: string;
  mode?: string;
}

/**
 * Parses a `component` query-parameter value.
 *
 * New values are JSON objects; the legacy `pattern;target;mode` form is still accepted
 * so previously shared search URLs keep working. The legacy parse splits from the right
 * because a SMARTS pattern may itself contain `;`.
 */
const parseComponent = (comp: string): ParsedComponent => {
  const parseLegacy = (): ParsedComponent => {
    const parts = comp.split(';');
    const mode = parts.pop();
    const target = parts.pop() ?? 'input';
    return { pattern: parts.join(';'), target, mode };
  };
  try {
    const parsed = JSON.parse(comp) as unknown;
    // JSON.parse succeeds for null/numbers/etc.; require the component shape.
    if (
      parsed !== null &&
      typeof parsed === 'object' &&
      typeof (parsed as ParsedComponent).pattern === 'string' &&
      typeof (parsed as ParsedComponent).target === 'string'
    ) {
      return parsed as ParsedComponent;
    }
    return parseLegacy();
  } catch {
    return parseLegacy();
  }
};

interface ReagentOptions {
  reactants: ReagentComponent[];
  products: ReagentComponent[];
  matchMode: string;
  useStereochemistry: boolean;
  similarityThreshold: number;
}

interface ReactionOptions {
  reactionIds: string[];
  reactionSmarts: string[];
  min_yield: number;
  max_yield: number;
  min_conversion: number;
  max_conversion: number;
}

interface DatasetOptions {
  datasetIds: string[];
  DOIs: string[];
}

interface SearchParams {
  limit: number;
}

interface SearchOptionsData {
  reagent: {
    reagents: Array<{ smileSmart: string; source: string; matchMode: string }>;
    useStereochemistry: boolean;
    similarityThreshold: number;
  };
  reaction: ReactionOptions;
  dataset: DatasetOptions;
  general: SearchParams;
}

interface SearchOptionsProps {
  onSearchOptions: (options: SearchOptionsData) => void;
}

const MATCH_MODES = ['exact', 'similar', 'substructure', 'SMARTS'];

interface ComponentListProps {
  title: string;
  components: ReagentComponent[];
  onDraw: (idx: number) => void;
  onUpdate: (idx: number, smiles: string) => void;
  onRemove: (idx: number) => void;
  onAdd: () => void;
}

const ComponentList: React.FC<ComponentListProps> = ({
  title,
  components,
  onDraw,
  onUpdate,
  onRemove,
  onAdd,
}) => (
  <div className={classes.componentList}>
    <Text className={classes.subtitle}>{title}</Text>
    {components.map((component, idx) => (
      <Flex
        key={idx}
        gap={4}
        align="center"
      >
        <ActionIcon
          variant="default"
          size="lg"
          onClick={() => onDraw(idx)}
          aria-label="Draw structure"
        >
          <IconPencil size={16} />
        </ActionIcon>
        <TextInput
          className={classes.componentInput}
          size="xs"
          placeholder="SMILES or SMARTS"
          value={component.smileSmart}
          onChange={event => onUpdate(idx, event.currentTarget.value)}
        />
        <ActionIcon
          variant="transparent"
          color="secondary.1"
          onClick={() => onRemove(idx)}
          aria-label="Remove component"
        >
          <IconTrash size={16} />
        </ActionIcon>
      </Flex>
    ))}
    {!components.length && <Text className={classes.muted}>No components</Text>}
    <Button
      variant="transparent"
      size="compact-sm"
      leftSection={<IconPlus size={14} />}
      className={classes.addButton}
      onClick={onAdd}
    >
      Add component
    </Button>
  </div>
);

const SearchOptions: React.FC<SearchOptionsProps> = ({ onSearchOptions }) => {
  const search = useSearch();
  const [openSections, setOpenSections] = useState<string[]>([]);

  const [reagentOptions, setReagentOptions] = useState<ReagentOptions>({
    reactants: [],
    products: [],
    matchMode: 'exact',
    useStereochemistry: false,
    similarityThreshold: 0.5,
  });

  const [reactionOptions, setReactionOptions] = useState<ReactionOptions>({
    reactionIds: [],
    reactionSmarts: [],
    min_yield: 50,
    max_yield: 100,
    min_conversion: 50,
    max_conversion: 100,
  });

  const [datasetOptions, setDatasetOptions] = useState<DatasetOptions>({
    datasetIds: [],
    DOIs: [],
  });

  const [searchParams, setSearchParams] = useState<SearchParams>({
    limit: 100,
  });

  const [showKetcherModal, setShowKetcherModal] = useState(false);
  const [ketcherModalSmile, setKetcherModalSmile] = useState(0);
  const [ketcherModalSet, setKetcherModalSet] = useState<'reactants' | 'products'>(
    'reactants',
  );

  const defaultQuery = useMemo<Record<string, string | string[]>>(() => {
    const params = new URLSearchParams(search);
    const query: Record<string, string | string[]> = {};
    for (const [key, value] of params.entries()) {
      const existing = query[key];
      if (existing === undefined) {
        query[key] = value;
      } else if (Array.isArray(existing)) {
        existing.push(value);
      } else {
        query[key] = [existing, value];
      }
    }
    return query;
  }, [search]);

  const openKetcherModal = (componentSet: 'reactants' | 'products', idx: number) => {
    setKetcherModalSmile(idx);
    setKetcherModalSet(componentSet);
    setShowKetcherModal(true);
  };

  const emitSearchOptions = () => {
    const allComponents = [...reagentOptions.reactants, ...reagentOptions.products].map(
      component => ({
        ...component,
        matchMode: reagentOptions.matchMode,
      }),
    );

    const searchOptions: SearchOptionsData = {
      reagent: {
        reagents: allComponents,
        useStereochemistry: reagentOptions.useStereochemistry,
        similarityThreshold: reagentOptions.similarityThreshold,
      },
      reaction: reactionOptions,
      dataset: datasetOptions,
      general: searchParams,
    };

    onSearchOptions(searchOptions);
  };

  const setDefaultValues = useCallback(() => {
    const q = defaultQuery;
    const open: string[] = [];

    // Reset reagent options first, then rebuild from URL
    const initialReagentOptions: ReagentOptions = {
      reactants: [],
      products: [],
      matchMode: 'exact',
      useStereochemistry: false,
      similarityThreshold: 0.5,
    };

    // reagent options
    if (q.component?.length) {
      const components = Array.isArray(q.component) ? q.component : [q.component];

      components.forEach((comp: string) => {
        const parsed = parseComponent(comp);
        const compType = parsed.target === 'input' ? 'reactants' : 'products';
        const matchMode = parsed.mode || 'exact';
        const component: ReagentComponent = {
          smileSmart: parsed.pattern,
          source: parsed.target,
          matchMode,
        };

        initialReagentOptions[compType].push(component);
        // Use the last matchMode found (or first if none specified)
        if (matchMode) {
          initialReagentOptions.matchMode = matchMode;
        }
      });

      initialReagentOptions.useStereochemistry =
        q.use_stereochemistry === 'true' || false;
      initialReagentOptions.similarityThreshold = Number(q.similarity) || 0.5;

      open.push('components');
    }
    setReagentOptions(initialReagentOptions);

    // dataset options
    const datasetIds = Array.isArray(q.dataset_id)
      ? q.dataset_id
      : q.dataset_id?.length
        ? [q.dataset_id]
        : [];
    const DOIs = Array.isArray(q.doi) ? q.doi : q.doi?.length ? [q.doi] : [];

    setDatasetOptions({ datasetIds, DOIs });
    if (datasetIds.length || DOIs.length) {
      open.push('datasets');
    }

    // reaction options
    const reactionIds = Array.isArray(q.reaction_id)
      ? q.reaction_id
      : q.reaction_id?.length
        ? [q.reaction_id]
        : [];
    const reactionSmarts = Array.isArray(q.reaction_smarts)
      ? q.reaction_smarts
      : q.reaction_smarts?.length
        ? [q.reaction_smarts]
        : [];

    setReactionOptions(prev => ({
      ...prev,
      reactionIds,
      reactionSmarts,
      min_yield: Number(q.min_yield) || 0,
      max_yield: Number(q.max_yield) || 100,
      min_conversion: Number(q.min_conversion) || 0,
      max_conversion: Number(q.max_conversion) || 100,
    }));

    // Number(q.max_yield) !== 100 used to trip on every page load because
    // Number(undefined) is NaN and NaN !== 100 is true; gate each bound on
    // its URL param being present before comparing.
    if (
      reactionIds.length ||
      reactionSmarts.length ||
      q.min_yield !== undefined ||
      (q.max_yield !== undefined && Number(q.max_yield) !== 100) ||
      q.min_conversion !== undefined ||
      (q.max_conversion !== undefined && Number(q.max_conversion) !== 100)
    ) {
      open.push('reactions');
    }

    // general search params
    setSearchParams({ limit: Number(q.limit) || 100 });
    setOpenSections(open.length ? open : ['components']);
  }, [defaultQuery]);

  const updateComponentSmiles = (
    componentSet: 'reactants' | 'products',
    idx: number,
    smiles: string,
  ) => {
    setReagentOptions(prev => ({
      ...prev,
      [componentSet]: prev[componentSet].map((component, index) =>
        index === idx ? { ...component, smileSmart: smiles } : component,
      ),
    }));
  };

  const addComponent = (componentSet: 'reactants' | 'products') => {
    const source = componentSet === 'reactants' ? 'input' : 'output';
    setReagentOptions(prev => ({
      ...prev,
      [componentSet]: [...prev[componentSet], { smileSmart: '', source }],
    }));
  };

  const removeComponent = (componentSet: 'reactants' | 'products', idx: number) => {
    setReagentOptions(prev => ({
      ...prev,
      [componentSet]: prev[componentSet].filter((_, index) => index !== idx),
    }));
  };

  useEffect(() => {
    setDefaultValues();
  }, [setDefaultValues]);

  return (
    <div className={classes.searchOptions}>
      <Accordion
        multiple
        value={openSections}
        onChange={setOpenSections}
        radius="sm"
      >
        <Accordion.Item value="components">
          <Accordion.Control>Components</Accordion.Control>
          <Accordion.Panel>
            <SegmentedControl
              fullWidth
              size="xs"
              value={reagentOptions.matchMode}
              onChange={mode =>
                setReagentOptions(prev => ({ ...prev, matchMode: mode }))
              }
              data={MATCH_MODES}
            />

            <Checkbox
              className={classes.stereoCheckbox}
              size="sm"
              label="Use stereochemistry"
              checked={reagentOptions.useStereochemistry}
              onChange={event =>
                setReagentOptions(prev => ({
                  ...prev,
                  useStereochemistry: event.currentTarget.checked,
                }))
              }
            />

            {reagentOptions.matchMode === 'similar' && (
              <div className={classes.similarity}>
                <Text className={classes.subtitle}>
                  Similarity threshold: {reagentOptions.similarityThreshold.toFixed(2)}
                </Text>
                <Slider
                  min={0.1}
                  max={1}
                  step={0.01}
                  label={value => value.toFixed(2)}
                  value={reagentOptions.similarityThreshold}
                  onChange={value =>
                    setReagentOptions(prev => ({
                      ...prev,
                      similarityThreshold: value,
                    }))
                  }
                />
              </div>
            )}

            <ComponentList
              title="Reactants & reagents"
              components={reagentOptions.reactants}
              onDraw={idx => openKetcherModal('reactants', idx)}
              onUpdate={(idx, smiles) =>
                updateComponentSmiles('reactants', idx, smiles)
              }
              onRemove={idx => removeComponent('reactants', idx)}
              onAdd={() => addComponent('reactants')}
            />

            <ComponentList
              title="Products"
              components={reagentOptions.products}
              onDraw={idx => openKetcherModal('products', idx)}
              onUpdate={(idx, smiles) => updateComponentSmiles('products', idx, smiles)}
              onRemove={idx => removeComponent('products', idx)}
              onAdd={() => addComponent('products')}
            />
          </Accordion.Panel>
        </Accordion.Item>

        <Accordion.Item value="reactions">
          <Accordion.Control>Reactions</Accordion.Control>
          <Accordion.Panel>
            <SearchItemList
              title="Reaction IDs"
              itemList={reactionOptions.reactionIds}
              onUpdateItemList={newList =>
                setReactionOptions(prev => ({ ...prev, reactionIds: newList }))
              }
            />
            <SearchItemList
              title="Reaction SMARTS"
              itemList={reactionOptions.reactionSmarts}
              onUpdateItemList={newList =>
                setReactionOptions(prev => ({ ...prev, reactionSmarts: newList }))
              }
            />

            <div className={classes.rangeBlock}>
              <Text className={classes.subtitle}>
                Yield: {reactionOptions.min_yield}% – {reactionOptions.max_yield}%
              </Text>
              <RangeSlider
                min={0}
                max={100}
                step={1}
                label={value => `${value}%`}
                value={[reactionOptions.min_yield, reactionOptions.max_yield]}
                onChange={([min, max]) =>
                  setReactionOptions(prev => ({
                    ...prev,
                    min_yield: min,
                    max_yield: max,
                  }))
                }
              />
            </div>

            <div className={classes.rangeBlock}>
              <Text className={classes.subtitle}>
                Conversion: {reactionOptions.min_conversion}% –{' '}
                {reactionOptions.max_conversion}%
              </Text>
              <RangeSlider
                min={0}
                max={100}
                step={1}
                label={value => `${value}%`}
                value={[reactionOptions.min_conversion, reactionOptions.max_conversion]}
                onChange={([min, max]) =>
                  setReactionOptions(prev => ({
                    ...prev,
                    min_conversion: min,
                    max_conversion: max,
                  }))
                }
              />
            </div>
          </Accordion.Panel>
        </Accordion.Item>

        <Accordion.Item value="datasets">
          <Accordion.Control>Datasets</Accordion.Control>
          <Accordion.Panel>
            <SearchItemList
              title="Dataset IDs"
              itemList={datasetOptions.datasetIds}
              onUpdateItemList={newList =>
                setDatasetOptions(prev => ({ ...prev, datasetIds: newList }))
              }
            />
            <SearchItemList
              title="DOIs"
              itemList={datasetOptions.DOIs}
              onUpdateItemList={newList =>
                setDatasetOptions(prev => ({ ...prev, DOIs: newList }))
              }
            />
          </Accordion.Panel>
        </Accordion.Item>
      </Accordion>

      <div className={classes.footer}>
        <NumberInput
          label="Result limit"
          size="xs"
          min={0}
          value={searchParams.limit}
          onChange={value =>
            setSearchParams(prev => ({ ...prev, limit: Number(value) || 0 }))
          }
        />
        <Button
          fullWidth
          leftSection={<IconSearch size={16} />}
          onClick={emitSearchOptions}
        >
          Search
        </Button>
      </div>

      {showKetcherModal && (
        <Suspense fallback={<LoadingSpinner />}>
          <ModalKetcher
            smiles={
              reagentOptions[ketcherModalSet][ketcherModalSmile]?.smileSmart || ''
            }
            onUpdateSmiles={newSmiles =>
              updateComponentSmiles(ketcherModalSet, ketcherModalSmile, newSmiles)
            }
            onCloseModal={() => setShowKetcherModal(false)}
          />
        </Suspense>
      )}
    </div>
  );
};

export default SearchOptions;
