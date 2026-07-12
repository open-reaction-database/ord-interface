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

/* eslint-disable react-refresh/only-export-components -- exports a hook
   alongside the component. */
// Read-only molecule image, mirroring ord-app's ReactionComponentPreview. The
// SVG comes from the backend's RDKit renderer (/api/compound_svg) instead of
// ord-app's Redux preview store; react-query dedupes repeat compounds.
import React from 'react';
import { Loader, Tooltip } from '@mantine/core';
import { IconPhotoOff } from '@tabler/icons-react';
import { useQuery } from '@tanstack/react-query';
import { ord } from 'ord-schema-protobufjs';
import { api } from '../../utils/api';
import { encodeCompound } from '../../utils/proto';
import classes from './reactionPreview.module.scss';

// Identifier types the backend renderer understands.
const STRUCTURAL_TYPES: number[] = [
  ord.CompoundIdentifier.CompoundIdentifierType.SMILES,
  ord.CompoundIdentifier.CompoundIdentifierType.INCHI,
  ord.CompoundIdentifier.CompoundIdentifierType.MOLBLOCK,
];

export type ComponentLike = ord.ICompound & ord.IProductCompound;

interface CompoundPreviewProps {
  component: ComponentLike | null | undefined;
  alt?: string;
}

export function useCompoundSvg(component: ComponentLike | null | undefined) {
  const structural = (component?.identifiers ?? []).filter(
    identifier => identifier.type != null && STRUCTURAL_TYPES.includes(identifier.type),
  );
  const key = structural.map(identifier => identifier.value ?? '').join('|');

  return useQuery<string | null>({
    queryKey: ['compound-svg', key],
    enabled: structural.length > 0,
    staleTime: Infinity,
    retry: false,
    queryFn: async () => {
      const binary = encodeCompound({ identifiers: structural });
      const response = await api.post<string>('/compound_svg', binary, {
        headers: { 'Content-Type': 'application/x-protobuf' },
      });
      return response.data;
    },
  });
}

export const CompoundPreview: React.FC<CompoundPreviewProps> = ({ component, alt }) => {
  const { data: svg, isLoading, isError } = useCompoundSvg(component);

  const hasStructure = (component?.identifiers ?? []).some(
    identifier => identifier.type != null && STRUCTURAL_TYPES.includes(identifier.type),
  );

  if (hasStructure && isLoading) {
    return <Loader size="sm" />;
  }

  if (!hasStructure || isError || !svg) {
    return (
      <Tooltip label="No preview">
        <div className={classes.noPreviewWrapper}>
          <IconPhotoOff className={classes.noPreviewIcon} />
        </div>
      </Tooltip>
    );
  }

  return (
    <div
      className={classes.compoundSvg}
      role="img"
      aria-label={alt ?? 'structure'}
      dangerouslySetInnerHTML={{ __html: svg }}
    />
  );
};

export default CompoundPreview;
