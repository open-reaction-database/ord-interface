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

import React, { useState, useEffect, useMemo, useCallback } from 'react';
import { useLocation } from 'react-router-dom';
import { Range } from 'react-range';
import ModalKetcher from '../../components/ModalKetcher';
import SearchItemList from './SearchItemList';
import './SearchOptions.scss';

interface ReagentComponent {
  smileSmart: string;
  source: string;
  matchMode?: string;
}

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

const SearchOptions: React.FC<SearchOptionsProps> = ({ onSearchOptions }) => {
  const location = useLocation();
  const [showReagentOptions, setShowReagentOptions] = useState(false);
  const [showReactionOptions, setShowReactionOptions] = useState(false);
  const [showDatasetOptions, setShowDatasetOptions] = useState(false);
  
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
  const [ketcherModalSet, setKetcherModalSet] = useState<'reactants' | 'products'>('reactants');
  
  const matchModes = ['exact', 'similar', 'substructure', 'SMARTS'];

  const defaultQuery = useMemo(() => {
    const params = new URLSearchParams(location.search);
    const query: any = {};
    for (const [key, value] of params.entries()) {
      if (query[key]) {
        if (Array.isArray(query[key])) {
          query[key].push(value);
        } else {
          query[key] = [query[key], value];
        }
      } else {
        query[key] = value;
      }
    }
    return query;
  }, [location.search]);

  const simThresholdDisplay = useMemo(() => {
    let trailingZeros = '';
    const simThresh = reagentOptions.similarityThreshold.toString();
    if (simThresh?.length < 2) trailingZeros = '.00';
    else if (simThresh?.length < 4) trailingZeros = '0';
    return simThresh + trailingZeros;
  }, [reagentOptions.similarityThreshold]);

  const openKetcherModal = (componentSet: 'reactants' | 'products', idx: number) => {
    setKetcherModalSmile(idx);
    setKetcherModalSet(componentSet);
    setShowKetcherModal(true);
  };

  const emitSearchOptions = () => {
    const allComponents = [...reagentOptions.reactants, ...reagentOptions.products]
      .map(component => ({ ...component, matchMode: reagentOptions.matchMode }));

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

  const addCompToOptions = useCallback((comp: string) => {
    const compArray = comp.split(';');
    const compType = compArray[1] === 'input' ? 'reactants' : 'products';
    const matchMode = compArray[2];
    const component: ReagentComponent = {
      smileSmart: compArray[0].replaceAll('%3D', '='),
      source: compArray[1],
      matchMode,
    };
    
    setReagentOptions(prev => ({
      ...prev,
      matchMode,
      [compType]: [...prev[compType], component],
    }));
  }, []);

  const setDefaultValues = useCallback(() => {
    const q = defaultQuery;
    
    // reagent options
    if (q.component?.length) {
      if (Array.isArray(q.component)) {
        q.component.forEach((comp: string) => addCompToOptions(comp));
      } else {
        addCompToOptions(q.component);
      }
      setReagentOptions(prev => ({
        ...prev,
        useStereochemistry: q.use_stereochemistry === 'true' || false,
        similarityThreshold: Number(q.similarity) || 0.5,
      }));
      setShowReagentOptions(true);
    }

    // dataset options
    const datasetIds = Array.isArray(q.dataset_id) ? q.dataset_id : q.dataset_id?.length ? [q.dataset_id] : [];
    const DOIs = Array.isArray(q.doi) ? q.doi : q.doi?.length ? [q.doi] : [];
    
    setDatasetOptions({ datasetIds, DOIs });
    if (datasetIds.length || DOIs.length) {
      setShowDatasetOptions(true);
    }

    // reaction options
    const reactionIds = Array.isArray(q.reaction_id) ? q.reaction_id : q.reaction_id?.length ? [q.reaction_id] : [];
    const reactionSmarts = Array.isArray(q.reaction_smarts) ? q.reaction_smarts : q.reaction_smarts?.length ? [q.reaction_smarts] : [];
    
    setReactionOptions(prev => ({
      ...prev,
      reactionIds,
      reactionSmarts,
      min_yield: Number(q.min_yield) || 0,
      max_yield: Number(q.max_yield) || 100,
      min_conversion: Number(q.min_conversion) || 0,
      max_conversion: Number(q.max_conversion) || 100,
    }));
    
    if (reactionIds.length || reactionSmarts.length || Number(q.min_yield) || Number(q.max_yield) !== 100) {
      setShowReactionOptions(true);
    }

    // general search params
    setSearchParams({ limit: Number(q.limit) || 100 });
  }, [defaultQuery, addCompToOptions]);

  const addReactant = () => {
    setReagentOptions(prev => ({
      ...prev,
      reactants: [...prev.reactants, { smileSmart: '', source: 'input' }],
    }));
  };

  const removeReactant = (idx: number) => {
    setReagentOptions(prev => ({
      ...prev,
      reactants: prev.reactants.filter((_, index) => index !== idx),
    }));
  };

  const addProduct = () => {
    setReagentOptions(prev => ({
      ...prev,
      products: [...prev.products, { smileSmart: '', source: 'output' }],
    }));
  };

  const removeProduct = (idx: number) => {
    setReagentOptions(prev => ({
      ...prev,
      products: prev.products.filter((_, index) => index !== idx),
    }));
  };

  const updateReactantSmiles = (idx: number, smiles: string) => {
    setReagentOptions(prev => ({
      ...prev,
      reactants: prev.reactants.map((reactant, index) => 
        index === idx ? { ...reactant, smileSmart: smiles } : reactant
      ),
    }));
  };

  const updateProductSmiles = (idx: number, smiles: string) => {
    setReagentOptions(prev => ({
      ...prev,
      products: prev.products.map((product, index) => 
        index === idx ? { ...product, smileSmart: smiles } : product
      ),
    }));
  };

  const updateKetcherSmiles = (newSmiles: string) => {
    if (ketcherModalSet === 'reactants') {
      updateReactantSmiles(ketcherModalSmile, newSmiles);
    } else {
      updateProductSmiles(ketcherModalSmile, newSmiles);
    }
  };

  useEffect(() => {
    setDefaultValues();
  }, [setDefaultValues]);

  return (
    <div className="search-options">
      {/* Components Section */}
      <div 
        className={`options-title ${showReagentOptions ? '' : 'closed'}`}
        onClick={() => setShowReagentOptions(!showReagentOptions)}
      >
        <span>Components</span>
        <i className="material-icons">expand_less</i>
      </div>
      
      {showReagentOptions && (
        <div id="searchByReagent" className="options-container">
          <div className="section">
            <div className="subtitle">
              <div className="tabs">
                {matchModes.map(mode => (
                  <div
                    key={mode}
                    className={`tab capitalize ${reagentOptions.matchMode === mode ? 'selected' : ''}`}
                    onClick={() => setReagentOptions(prev => ({ ...prev, matchMode: mode }))}
                  >
                    {mode}
                  </div>
                ))}
              </div>
            </div>
          </div>
          
          <div className="section">
            <div className="subtitle">General Options</div>
            <div className="general options">
              <label htmlFor="stereo">Use Stereochemistry</label>
              <input
                id="stereo"
                type="checkbox"
                checked={reagentOptions.useStereochemistry}
                onChange={(e) => setReagentOptions(prev => ({ ...prev, useStereochemistry: e.target.checked }))}
              />
              {reagentOptions.matchMode === 'similar' && (
                <>
                  <label htmlFor="similarity">Similarity Threshold</label>
                  <div className="slider-input">
                    <div className="value">{simThresholdDisplay}</div>
                    <input
                      id="similarity"
                      type="range"
                      min="0.1"
                      max="1.0"
                      step="0.01"
                      value={reagentOptions.similarityThreshold}
                      onChange={(e) => setReagentOptions(prev => ({ ...prev, similarityThreshold: Number(e.target.value) }))}
                    />
                  </div>
                </>
              )}
            </div>
          </div>
          
          <div className="section">
            <div className="subtitle">Reactants & Reagents</div>
            <div className="reagent options">
              {reagentOptions.reactants.map((reactant, idx) => (
                <React.Fragment key={idx}>
                  <div className="draw">
                    <button onClick={() => openKetcherModal('reactants', idx)}>
                      <i className="material-icons">draw</i>
                    </button>
                  </div>
                  <div className="field long">
                    <input
                      type="text"
                      value={reactant.smileSmart}
                      onChange={(e) => updateReactantSmiles(idx, e.target.value)}
                    />
                  </div>
                  <div className="delete">
                    <button onClick={() => removeReactant(idx)}>
                      <i className="material-icons">delete</i>
                    </button>
                  </div>
                </React.Fragment>
              ))}
              {!reagentOptions.reactants?.length && (
                <div className="copy">No components</div>
              )}
              <div id="add-component">
                <button type="button" onClick={addReactant}>
                  <i className="material-icons">add</i>
                  Add Component
                </button>
              </div>
            </div>
          </div>
          
          <div className="section">
            <div className="subtitle">Products</div>
            <div className="reagent options">
              {reagentOptions.products.map((product, idx) => (
                <React.Fragment key={idx}>
                  <div className="draw">
                    <button onClick={() => openKetcherModal('products', idx)}>
                      <i className="material-icons">draw</i>
                    </button>
                  </div>
                  <div className="field long">
                    <input
                      type="text"
                      value={product.smileSmart}
                      onChange={(e) => updateProductSmiles(idx, e.target.value)}
                    />
                  </div>
                  <div className="delete">
                    <button onClick={() => removeProduct(idx)}>
                      <i className="material-icons">delete</i>
                    </button>
                  </div>
                </React.Fragment>
              ))}
              {!reagentOptions.products?.length && (
                <div className="copy">No components</div>
              )}
              <div id="add-component">
                <button type="button" onClick={addProduct}>
                  <i className="material-icons">add</i>
                  Add Component
                </button>
              </div>
            </div>
          </div>
        </div>
      )}

      {/* Reactions Section */}
      <div 
        className={`options-title ${showReactionOptions ? '' : 'closed'}`}
        onClick={() => setShowReactionOptions(!showReactionOptions)}
      >
        <span>Reactions</span>
        <i className="material-icons">expand_less</i>
      </div>
      
      {showReactionOptions && (
        <div id="searchByReaction" className="options-container">
          <SearchItemList
            title="Reaction IDs"
            itemList={reactionOptions.reactionIds}
            onUpdateItemList={(newList) => setReactionOptions(prev => ({ ...prev, reactionIds: newList }))}
          />
          <SearchItemList
            title="Reaction SMARTS"
            itemList={reactionOptions.reactionSmarts}
            onUpdateItemList={(newList) => setReactionOptions(prev => ({ ...prev, reactionSmarts: newList }))}
          />
          
          <div className="slider-input multi">
            <label htmlFor="yield">Yield</label>
            <div className="value">{reactionOptions.min_yield}% - {reactionOptions.max_yield}%</div>
            <div className="range-container">
              <Range
                step={1}
                min={0}
                max={100}
                values={[reactionOptions.min_yield, reactionOptions.max_yield]}
                onChange={(values) => setReactionOptions(prev => ({ 
                  ...prev, 
                  min_yield: values[0], 
                  max_yield: values[1] 
                }))}
                renderTrack={({ props, children }) => (
                  <div
                    {...props}
                    className="range-track"
                    style={props.style}
                  >
                    {children}
                  </div>
                )}
                renderThumb={({ props }) => (
                  <div
                    {...props}
                    className="range-thumb"
                    style={props.style}
                  />
                )}
              />
            </div>
          </div>
          
          <div className="slider-input multi">
            <label htmlFor="conversion">Conversion</label>
            <div className="value">{reactionOptions.min_conversion}% - {reactionOptions.max_conversion}%</div>
            <div className="range-container">
              <Range
                step={1}
                min={0}
                max={100}
                values={[reactionOptions.min_conversion, reactionOptions.max_conversion]}
                onChange={(values) => setReactionOptions(prev => ({ 
                  ...prev, 
                  min_conversion: values[0], 
                  max_conversion: values[1] 
                }))}
                renderTrack={({ props, children }) => (
                  <div
                    {...props}
                    className="range-track"
                    style={props.style}
                  >
                    {children}
                  </div>
                )}
                renderThumb={({ props }) => (
                  <div
                    {...props}
                    className="range-thumb"
                    style={props.style}
                  />
                )}
              />
            </div>
          </div>
        </div>
      )}

      {/* Datasets Section */}
      <div 
        className={`options-title ${showDatasetOptions ? '' : 'closed'}`}
        onClick={() => setShowDatasetOptions(!showDatasetOptions)}
      >
        <span>Datasets</span>
        <i className="material-icons">expand_less</i>
      </div>
      
      {showDatasetOptions && (
        <div id="searchByDataset" className="options-container">
          <SearchItemList
            title="Dataset IDs"
            itemList={datasetOptions.datasetIds}
            onUpdateItemList={(newList) => setDatasetOptions(prev => ({ ...prev, datasetIds: newList }))}
          />
          <SearchItemList
            title="DOIs"
            itemList={datasetOptions.DOIs}
            onUpdateItemList={(newList) => setDatasetOptions(prev => ({ ...prev, DOIs: newList }))}
          />
        </div>
      )}

      {/* Search Parameters Section */}
      <div id="searchParameters" className="options-title">
        Search Parameters
      </div>
      <div className="options-container">
        <div className="section">
          <label htmlFor="limit">Result Limit</label>
          <input
            id="limit"
            type="number"
            min="0"
            value={searchParams.limit}
            onChange={(e) => setSearchParams(prev => ({ ...prev, limit: Number(e.target.value) }))}
          />
        </div>
        <div className="search-button">
          <button onClick={emitSearchOptions}>
            <b>Search</b>
          </button>
        </div>
      </div>

      {/* Ketcher Modal */}
      {showKetcherModal && (
        <ModalKetcher
          smiles={reagentOptions[ketcherModalSet][ketcherModalSmile]?.smileSmart || ''}
          onUpdateSmiles={updateKetcherSmiles}
          onCloseModal={() => setShowKetcherModal(false)}
        />
      )}
    </div>
  );
};

export default SearchOptions;