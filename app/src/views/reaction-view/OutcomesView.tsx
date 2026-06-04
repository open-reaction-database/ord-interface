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

import React, { useState } from 'react';
import reaction_pb from 'ord-schema';
import type { Analysis, ProductMeasurement } from 'ord-schema/proto/reaction_pb';
import CompoundView from './CompoundView';
import FloatingModal from '../../components/FloatingModal';
import { amountObj, amountStr } from '../../utils/amount';
import { enumName } from '../../utils/enum';
import { formatPercentage, formattedTime } from '../../utils/outcomes';
import type { ReactionOutcomeData } from '../../types/search';
import './OutcomesView.scss';

interface OutcomesViewProps {
  outcome: ReactionOutcomeData | undefined;
}

const measurementType = (type: number | undefined): string =>
  enumName(reaction_pb.ProductMeasurement.ProductMeasurementType, type) ?? '';

const analysisType = (type: number | undefined): string => enumName(reaction_pb.Analysis.AnalysisType, type) ?? '';

const measurementValue = (measurement: ProductMeasurement.AsObject): string => {
  if (measurement.percentage) return formatPercentage(measurement.percentage);
  if (measurement.amount) return amountStr(amountObj(measurement.amount));
  if (measurement.floatValue) return String(measurement.floatValue.value);
  if (measurement.stringValue) return measurement.stringValue;
  return '';
};

const measurementWithNamedType = (measurement: ProductMeasurement.AsObject) => ({
  ...measurement,
  type: measurementType(measurement.type),
});

const analysisWithNamedType = (analysis: Analysis.AsObject) => ({
  ...analysis,
  type: analysisType(analysis.type),
});

const OutcomesView: React.FC<OutcomesViewProps> = ({ outcome }) => {
  const [productsIdx, setProductsIdx] = useState(0);
  const [analysesIdx, setAnalysesIdx] = useState(0);
  const [rawMeasurement, setRawMeasurement] = useState<ReturnType<typeof measurementWithNamedType> | null>(null);
  const [rawAnalysis, setRawAnalysis] = useState<ReturnType<typeof analysisWithNamedType> | null>(null);
  const [customDetails, setCustomDetails] = useState<string | null>(null);

  if (!outcome) return null;

  const reactionTime = formattedTime(outcome.reactionTime);
  const conversion = outcome.conversion;
  const showDetails = Boolean(reactionTime || conversion);

  const currentProduct = outcome.productsList?.[productsIdx];
  const currentAnalysis = outcome.analysesMap?.[analysesIdx]?.[1];

  return (
    <div className="outcomes-view">
      {showDetails && (
        <>
          <div className="title">Details</div>
          <div className="details">
            {reactionTime && (
              <>
                <div className="label">Reaction Time</div>
                <div className="value">{reactionTime}</div>
              </>
            )}
            {conversion && (
              <>
                <div className="label">Conversion</div>
                <div className="value">{formatPercentage(conversion)}</div>
              </>
            )}
          </div>
        </>
      )}

      <div className="title">Products</div>
      <div className="sub-section">
        <div className="tabs">
          {outcome.productsList?.map((_product, idx) => (
            <div
              key={idx}
              className={`tab ${productsIdx === idx ? 'selected' : ''}`}
              onClick={() => setProductsIdx(idx)}
            >
              Product {idx + 1}
            </div>
          ))}
        </div>

        {currentProduct && (
          <>
            <div className="compound">
              <CompoundView component={currentProduct} />
            </div>

            <div className="sub-title">Measurements</div>
            <div className="measurements">
              <div className="label">Type</div>
              <div className="label" />
              <div className="label">Value</div>
              <div className="label">Analysis</div>
              <div className="label">Raw</div>
              {currentProduct.measurementsList.map((measurement, idx) => {
                const typeName = measurementType(measurement.type);
                return (
                  <React.Fragment key={idx}>
                    {typeName === 'CUSTOM' ? (
                      <div className="value">
                        <div
                          className="button"
                          onClick={() =>
                            setCustomDetails(
                              measurement.details ||
                                'Please contact the author for details on this custom measurement.',
                            )
                          }
                        >
                          <u>CUSTOM</u>
                        </div>
                      </div>
                    ) : (
                      <div className="value">{typeName}</div>
                    )}
                    <div className="value" />
                    <div className="value">{measurementValue(measurement)}</div>
                    <div className="value">{measurement.analysisKey}</div>
                    <div className="value">
                      <div className="raw">
                        <div
                          className="button"
                          onClick={() => setRawMeasurement(measurementWithNamedType(measurement))}
                        >
                          &lt;&gt;
                        </div>
                      </div>
                    </div>
                  </React.Fragment>
                );
              })}
            </div>
          </>
        )}

        {rawMeasurement && (
          <FloatingModal
            title="Raw Data"
            onCloseModal={() => setRawMeasurement(null)}
          >
            <div className="data">
              <pre>{JSON.stringify(rawMeasurement, null, 2)}</pre>
            </div>
          </FloatingModal>
        )}

        {customDetails !== null && (
          <FloatingModal
            title="Custom Measurement Details"
            onCloseModal={() => setCustomDetails(null)}
          >
            <div className="data">
              <pre>{customDetails}</pre>
            </div>
          </FloatingModal>
        )}
      </div>

      {outcome.analysesMap && outcome.analysesMap.length > 0 && (
        <>
          <div className="title">Analyses</div>
          <div className="sub-section">
            <div className="tabs">
              {outcome.analysesMap.map(([key], idx) => (
                <div
                  key={key}
                  className={`tab ${analysesIdx === idx ? 'selected' : ''}`}
                  onClick={() => setAnalysesIdx(idx)}
                >
                  {key}
                </div>
              ))}
            </div>
            {currentAnalysis && (
              <div className="details">
                <div className="label">Type</div>
                <div className="value">{analysisType(currentAnalysis.type)}</div>
                <div className="label">Details</div>
                <div className="value">{currentAnalysis.details}</div>
                <div className="label">Raw</div>
                <div className="value">
                  <div className="raw">
                    <div
                      className="button"
                      onClick={() => setRawAnalysis(analysisWithNamedType(currentAnalysis))}
                    >
                      &lt;&gt;
                    </div>
                  </div>
                </div>
              </div>
            )}
          </div>

          {rawAnalysis && (
            <FloatingModal
              title="Raw Data"
              onCloseModal={() => setRawAnalysis(null)}
            >
              <div className="data">
                <pre>{JSON.stringify(rawAnalysis, null, 2)}</pre>
              </div>
            </FloatingModal>
          )}
        </>
      )}
    </div>
  );
};

export default OutcomesView;
