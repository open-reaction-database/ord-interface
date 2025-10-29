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

import React from 'react';
import './ProvenanceView.scss';

interface ProvenanceViewProps {
  provenance: any;
}

const ProvenanceView: React.FC<ProvenanceViewProps> = ({ provenance }) => {
  return (
    <div className="provenance-view">
      {provenance?.experimenter && (
        <>
          <div className="title">Experimenter</div>
          <div className="details experimenter">
            <div className="label">Username</div>
            <div className="label">Name</div>
            <div className="label">ORCID</div>
            <div className="label">Organization</div>
            <div className="label">Email</div>
            <div className="value">{provenance.experimenter.username}</div>
            <div className="value">{provenance.experimenter.name}</div>
            <div className="value">{provenance.experimenter.orcid}</div>
            <div className="value">{provenance.experimenter.organization}</div>
            <div className="value">{provenance.experimenter.email}</div>
          </div>
        </>
      )}
      
      <div className="details">
        {provenance?.city && (
          <>
            <div className="label">City</div>
            <div className="value">{provenance.city}</div>
          </>
        )}
        
        {provenance?.experimentStart && (
          <>
            <div className="label">Experiment Start</div>
            <div className="value">{new Date(provenance.experimentStart).toLocaleDateString()}</div>
          </>
        )}
        
        {provenance?.doi && (
          <>
            <div className="label">DOI</div>
            <div className="value">{provenance.doi}</div>
          </>
        )}
        
        {provenance?.patent && (
          <>
            <div className="label">Patent</div>
            <div className="value">{provenance.patent}</div>
          </>
        )}
        
        {provenance?.publicationUrl && (
          <>
            <div className="label">Publication URL</div>
            <div className="value">
              <a href={provenance.publicationUrl} target="_blank" rel="noopener noreferrer">
                {provenance.publicationUrl}
              </a>
            </div>
          </>
        )}
      </div>
    </div>
  );
};

export default ProvenanceView;