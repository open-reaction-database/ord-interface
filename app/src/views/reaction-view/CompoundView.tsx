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
import FloatingModal from '../../components/FloatingModal';
import './CompoundView.scss';

// Inline utility functions for now
const amountObj = (amount: any) => {
  if (!amount) return { unitCategory: '', unitAmount: '', unitType: '' };
  // Simplified implementation - would need full logic from utils/amount.js
  return { unitCategory: 'mass', unitAmount: amount.value || '', unitType: amount.units || '' };
};

const amountStr = (amountObj: any) => {
  if (!amountObj.unitAmount) return '';
  return `${amountObj.unitAmount} ${amountObj.unitType}`;
};

interface CompoundViewProps {
  component: any;
}

interface RawData {
  identifiers?: Array<{ type: string; value: string }>;
  amount?: any;
  preparations?: Array<{ type: string; details: string }>;
  is_desired_product?: boolean;
  isolated_color?: string;
  texture?: { type: string; details?: string };
  reaction_role: string;
}

const CompoundView: React.FC<CompoundViewProps> = ({ component }) => {
  const [compoundSVG, setCompoundSVG] = useState<string | null>(null);
  const [showRawData, setShowRawData] = useState(false);

  const compoundAmountObj = useMemo(() => {
    return amountObj(component?.amount);
  }, [component?.amount]);

  const compoundAmount = useMemo(() => {
    return amountStr(compoundAmountObj);
  }, [compoundAmountObj]);

  const compoundRole = useMemo(() => {
    if (!component?.reactionRole) return '';
    const role = component.reactionRole;
    const types = reaction_pb.ReactionRole.ReactionRoleType;
    return Object.keys(types).find(key => (types as any)[key] === role) || '';
  }, [component?.reactionRole]);

  const rawData: RawData = useMemo(() => {
    const returnObj: RawData = {
      reaction_role: compoundRole
    };
    
    // set identifiers
    if (component?.identifiersList?.length) {
      const idTypes = reaction_pb.CompoundIdentifier.CompoundIdentifierType;
      returnObj.identifiers = component.identifiersList.map((identifier: any) => ({
        type: Object.keys(idTypes).find(key => (idTypes as any)[key] === identifier.type) || '',
        value: identifier.value
      }));
    }
    
    // set amount
    if (component?.amount) {
      returnObj.amount = {
        [compoundAmountObj.unitCategory]: {
          value: compoundAmountObj.unitAmount,
          units: compoundAmountObj.unitType,
        }
      };
    }
    
    // set preparations
    if (component?.preparationsList?.length) {
      const prepTypes = reaction_pb.CompoundPreparation.CompoundPreparationType;
      returnObj.preparations = component.preparationsList.map((prep: any) => ({
        type: Object.keys(prepTypes).find(key => (prepTypes as any)[key] === prep.type) || '',
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
      // Simplified texture handling - would need proper enum mapping
      const textureObj: { type: string; details?: string } = { 
        type: component.texture.type?.toString() || '' 
      };
      if (component.texture.details) {
        textureObj.details = component.texture.details;
      }
      returnObj.texture = textureObj;
    }
    
    return returnObj;
  }, [component, compoundRole, compoundAmountObj]);

  const gridColumns = useMemo(() => {
    return `1fr repeat(${component?.amount ? 3 : 2}, auto)`;
  }, [component?.amount]);

  const getCompoundSVG = async (component: any) => {
    if (!component?.identifiersList?.length) return;
    
    const smilesIdentifier = component.identifiersList.find((identifier: any) => identifier.type === 2); // type 2 is SMILES
    if (!smilesIdentifier) return;
    
    const compoundStr = smilesIdentifier.value;
    
    // prep compound
    const compound = new reaction_pb.Compound();
    const identifier = compound.addIdentifiers();
    identifier.setValue(compoundStr);
    identifier.setType(reaction_pb.CompoundIdentifier.CompoundIdentifierType.SMILES);
    
    try {
      const response = await fetch('/api/compound_svg', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/x-protobuf',
        },
        body: compound.serializeBinary() as BodyInit,
      });
      
      const result = await response.json();
      setCompoundSVG(result);
    } catch (error) {
      console.error('Error fetching compound SVG:', error);
    }
  };

  useEffect(() => {
    if (component) {
      getCompoundSVG(component);
    }
  }, [component]);

  return (
    <div className="compound-view" style={{ gridTemplateColumns: gridColumns }}>
      <div className="label">Compound</div>
      {component?.amount && <div className="label">Amount</div>}
      <div className="label">Role</div>
      <div className="label">Raw</div>
      
      <div className="svg" dangerouslySetInnerHTML={{ __html: compoundSVG || '' }} />
      {component?.amount && <div className="amount">{compoundAmount}</div>}
      <div className="role">{compoundRole.toLowerCase()}</div>
      <div className="raw">
        <div className="button" onClick={() => setShowRawData(true)}>
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