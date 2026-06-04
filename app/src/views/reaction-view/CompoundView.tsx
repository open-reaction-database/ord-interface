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

import React, { useState, useEffect, useMemo } from 'react';
import reaction_pb from 'ord-schema';
import type { Compound, ProductCompound } from 'ord-schema/proto/reaction_pb';
import FloatingModal from '../../components/FloatingModal';
import { amountObj, amountStr } from '../../utils/amount';
import { enumName } from '../../utils/enum';
import './CompoundView.scss';

// CompoundView is invoked with both reaction inputs (Compound) and outcome
// products (ProductCompound); the two protobuf messages share most rendered
// fields but each carries a few of its own (Compound: reactionRole / preparationsList,
// ProductCompound: isDesiredProduct / isolatedColor / texture).
type ComponentLike = Partial<Compound.AsObject> & Partial<ProductCompound.AsObject>;

interface CompoundViewProps {
  component: ComponentLike | undefined;
}

interface RawData {
  identifiers?: Array<{ type: string; value: string }>;
  amount?: Record<string, { value: number | undefined; units: string }>;
  preparations?: Array<{ type: string; details: string }>;
  is_desired_product?: boolean;
  isolated_color?: string;
  texture?: { type: string; details?: string };
  reaction_role: string;
}

const CompoundView: React.FC<CompoundViewProps> = ({ component }) => {
  const [compoundSVG, setCompoundSVG] = useState<string | null>(null);
  const [showRawData, setShowRawData] = useState(false);

  const compoundAmountObj = useMemo(() => amountObj(component?.amount), [component?.amount]);
  const compoundAmount = useMemo(() => amountStr(compoundAmountObj), [compoundAmountObj]);

  const compoundRole = useMemo(() => {
    if (!component?.reactionRole) return '';
    return String(enumName(reaction_pb.ReactionRole.ReactionRoleType, component.reactionRole) ?? '');
  }, [component?.reactionRole]);

  const rawData: RawData = useMemo(() => {
    const returnObj: RawData = { reaction_role: compoundRole };

    if (component?.identifiersList?.length) {
      returnObj.identifiers = component.identifiersList.map(identifier => ({
        type: String(enumName(reaction_pb.CompoundIdentifier.CompoundIdentifierType, identifier.type) ?? ''),
        value: identifier.value,
      }));
    }

    if (component?.amount && compoundAmountObj.unitCategory) {
      returnObj.amount = {
        [compoundAmountObj.unitCategory]: {
          value: compoundAmountObj.unitAmount,
          units: compoundAmountObj.unitType ?? '',
        },
      };
    }

    if (component?.preparationsList?.length) {
      returnObj.preparations = component.preparationsList.map(prep => ({
        type: String(enumName(reaction_pb.CompoundPreparation.CompoundPreparationType, prep.type) ?? ''),
        details: prep.details,
      }));
    }

    if (component?.isDesiredProduct) {
      returnObj.is_desired_product = component.isDesiredProduct;
    }

    if (component?.isolatedColor) {
      returnObj.isolated_color = component.isolatedColor;
    }

    if (component?.texture) {
      const textureObj: { type: string; details?: string } = {
        type: component.texture.type !== undefined ? String(component.texture.type) : '',
      };
      if (component.texture.details) {
        textureObj.details = component.texture.details;
      }
      returnObj.texture = textureObj;
    }

    return returnObj;
  }, [component, compoundRole, compoundAmountObj]);

  const gridColumns = useMemo(() => `1fr repeat(${component?.amount ? 3 : 2}, auto)`, [component?.amount]);

  const smilesValue = useMemo(() => {
    if (!component?.identifiersList?.length) return null;
    const smilesType = reaction_pb.CompoundIdentifier.CompoundIdentifierType.SMILES;
    const smilesIdentifier = component.identifiersList.find(identifier => identifier.type === smilesType);
    return smilesIdentifier?.value ?? null;
  }, [component?.identifiersList]);

  useEffect(() => {
    if (!smilesValue) return;

    // prep compound
    const compound = new reaction_pb.Compound();
    const identifier = compound.addIdentifiers();
    identifier.setValue(smilesValue);
    identifier.setType(reaction_pb.CompoundIdentifier.CompoundIdentifierType.SMILES);

    const fetchSVG = async () => {
      try {
        const response = await fetch('/api/compound_svg', {
          method: 'POST',
          headers: {
            'Content-Type': 'application/x-protobuf',
          },
          body: compound.serializeBinary() as BodyInit,
        });

        // Skip the 4xx/5xx body — it's an HTML error page that
        // dangerouslySetInnerHTML would render verbatim in the SVG slot.
        if (!response.ok) {
          console.error(`compound_svg failed (HTTP ${response.status})`);
          return;
        }
        const result = await response.json();
        setCompoundSVG(result);
      } catch (error) {
        console.error('Error fetching compound SVG:', error);
      }
    };

    fetchSVG();
  }, [smilesValue]);

  return (
    <div
      className="compound-view"
      style={{ gridTemplateColumns: gridColumns }}
    >
      <div className="label">Compound</div>
      {component?.amount && <div className="label">Amount</div>}
      <div className="label">Role</div>
      <div className="label">Raw</div>

      <div
        className="svg"
        dangerouslySetInnerHTML={{ __html: compoundSVG || '' }}
      />
      {component?.amount && <div className="amount">{compoundAmount}</div>}
      <div className="role">{compoundRole.toLowerCase()}</div>
      <div className="raw">
        <div
          className="button"
          onClick={() => setShowRawData(true)}
        >
          &lt;&gt;
        </div>
      </div>

      {showRawData && (
        <FloatingModal
          title="Raw Data"
          onCloseModal={() => setShowRawData(false)}
        >
          <div className="data">
            {rawData.identifiers?.map((identifier, index) => (
              <pre key={index}>{`identifiers: ${JSON.stringify(identifier)}`}</pre>
            ))}
            {rawData.amount && <pre>amount: {JSON.stringify(rawData.amount)}</pre>}
            <pre>reaction_role: {rawData.reaction_role}</pre>
            {rawData.is_desired_product && <pre>is_desired_product: {rawData.is_desired_product}</pre>}
            {rawData.isolated_color && <pre>isolated_color: {rawData.isolated_color}</pre>}
            {rawData.texture && <pre>texture: {JSON.stringify(rawData.texture)}</pre>}
            {rawData.preparations?.map((prep, index) => (
              <pre key={index}>preparations: {JSON.stringify(prep)}</pre>
            ))}
          </div>
        </FloatingModal>
      )}
    </div>
  );
};

export default CompoundView;
