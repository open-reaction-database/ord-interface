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
import './About.scss';

const About: React.FC = () => {
  return (
    <div className="about container">
      <div className="about__section row">
        <div className="col-3 about__section-header">
          <h3 className="about__section-title">About</h3>
        </div>
        <div className="col-9 about__section-content">
          <p className="about__text">
            The Open Reaction Database (ORD) is an open-access schema and infrastructure for structuring and
            sharing organic reaction data, including a centralized data repository. The ORD schema supports
            conventional and emerging technologies, from benchtop reactions to automated high-throughput
            experiments and flow chemistry. Our vision is that a consistent data representation and
            infrastructure to support data sharing will enable downstream applications that will greatly improve
            the state of the art with respect to computer-aided synthesis planning, reaction prediction, and
            other predictive chemistry tasks.
          </p>
          <p className="about__text">
            Since our initial meeting in October 2019, the database has grown to include more than 2M reactions
            (including a large dataset of reactions extracted from USPTO sources) and received contributions
            from academic and industrial users, both from published and unpublished work. Some of our current
            efforts include:
          </p>
          <ul className="about__list">
            <li className="about__list-item">Improving user interfaces and providing support to contributors on GitHub and via email.</li>
            <li className="about__list-item">
              Working with companies to incorporate the ORD schema into their reaction informatics
              infrastructure, including the development of "translators" between the ORD schema and electronic
              lab notebooks (ELNs).
            </li>
            <li className="about__list-item">
              Engaging with journals and other stakeholders to drive adoption of the ORD schema as a{' '}
              <a href="https://en.wikipedia.org/wiki/FAIR_data" className="about__link">FAIR</a> data structure for sharing reaction
              data across academia, government, and industry.
            </li>
          </ul>
          <p className="about__text">
            <a href="https://www.linkedin.com/in/benjamin-deadman-78786242/" className="about__link">Ben Deadman</a>
            {' '}is our Reaction Data Evangelist and ORD Program Manager. Please reach out to him at{' '}
            <a href="mailto:help@open-reaction-database.org" className="about__link">help@open‑reaction‑database.org</a>
            {' '}for help preparing a contribution or to discuss using the ORD in your company or lab.
          </p>
        </div>
      </div>

      <div className="about__section row">
        <div className="col-3 about__section-header">
          <h3 className="about__section-title">Publications and Media</h3>
        </div>
        <div className="col-9 about__section-content">
          <h5 className="about__subsection-title">Journal Articles</h5>
          <ul className="about__list">
            <li className="about__list-item">
              Kearnes SM, Maser MR, Wleklinski M, Kast A, Doyle AG, Dreher SD, Hawkins JM, Jensen KF, Coley
              CW. The Open Reaction Database. <em>J Am Chem Soc</em> 2021, 143(45), 18820-18826. (
              <a href="https://doi.org/10.1021/jacs.1c09820" className="about__link">JACS</a>)
            </li>
            <li className="about__list-item">
              Mercado R, Kearnes SM, Coley C. Data Sharing in Chemistry: Lessons Learned and a Case for
              Mandating Structured Reaction Data. <em>J Chem Inf Model</em> 2023, 63(14), 4253-4265. (
              <a href="https://doi.org/10.1021/acs.jcim.3c00607" className="about__link">JCIM</a>)
            </li>
          </ul>
          <h5 className="about__subsection-title">News</h5>
          <ul className="about__list">
            <li className="about__list-item">
              A new database for machine-learning research (
              <a href="https://cen.acs.org/physical-chemistry/computational-chemistry/new-database-machine-learning-research/99/web/2021/11" className="about__link">C&EN</a>,
              22 November 2021)
            </li>
            <li className="about__list-item">
              Yield-predicting AI needs chemists to stop ignoring failed experiments (
              <a href="https://www.chemistryworld.com/news/yield-predicting-ai-needs-chemists-to-stop-ignoring-failed-experiments/4015662.article" className="about__link">Chemistry World</a>,
              12 May 2022)
            </li>
            <li className="about__list-item">
              Chemists debate machine learning's future in synthesis planning and ask for open data (
              <a href="https://cen.acs.org/physical-chemistry/computational-chemistry/Chemists-debate-machine-learnings-future/100/i18" className="about__link">C&EN</a>,
              18 May 2022)
            </li>
            <li className="about__list-item">
              For chemists, the AI revolution has yet to happen (
              <a href="https://www.nature.com/articles/d41586-023-01612-x" className="about__link">Nature</a>,
              17 May 2023)
            </li>
          </ul>
        </div>
      </div>

      <div className="about__section row">
        <div className="col-3 about__section-header">
          <h3 className="about__section-title">Leadership</h3>
        </div>
        <div className="col-9 about__section-content">
          <h5 className="about__subsection-title">Governing Committee</h5>
          <div className="about__leadership-grid">
            <ul className="about__leadership-column">
              <li className="about__list-item">Connor Coley (MIT, C-CAS)</li>
              <li className="about__list-item">Abby Doyle (UCLA, C-CAS)</li>
              <li className="about__list-item">Spencer Dreher (Merck)</li>
            </ul>
            <ul className="about__leadership-column">
              <li className="about__list-item">Joel Hawkins (Pfizer)</li>
              <li className="about__list-item">Klavs Jensen (MIT)</li>
              <li className="about__list-item">Steven Kearnes (Relay)</li>
            </ul>
          </div>
          <h5 className="about__subsection-title">Advisory Board</h5>
          <div className="about__leadership-grid">
            <ul className="about__leadership-column">
              <li className="about__list-item">Alán Aspuru-Guzik (Toronto, MADNESS)</li>
              <li className="about__list-item">Timothy Cernak (Michigan, Entos)</li>
              <li className="about__list-item">Lucy Colwell (Cambridge, SynTech, Google)</li>
              <li className="about__list-item">Werngard Czechtizky (AstraZeneca)</li>
              <li className="about__list-item">JW Feng</li>
              <li className="about__list-item">Matthew Gaunt (Cambridge, SynTech)</li>
              <li className="about__list-item">Alex Godfrey (NCATS Consultant)</li>
              <li className="about__list-item">Mimi Hii (Imperial, ROAR)</li>
              <li className="about__list-item">Greg Landrum (T5 Informatics)</li>
            </ul>
            <ul className="about__leadership-column">
              <li className="about__list-item">Fabio Lima (Novartis)</li>
              <li className="about__list-item">Christos Nicolaou (Recursion)</li>
              <li className="about__list-item">Sarah Reisman (Caltech, C-CAS)</li>
              <li className="about__list-item">Francesco Rianjongdee (GSK)</li>
              <li className="about__list-item">Marwin Segler (Microsoft)</li>
              <li className="about__list-item">Matthew Sigman (Utah, C-CAS)</li>
              <li className="about__list-item">Jay Stevens (BMS)</li>
              <li className="about__list-item">Sarah Trice (XtalPi)</li>
              <li className="about__list-item">Huimin Zhao (UIUC, MMLI)</li>
            </ul>
          </div>
        </div>
      </div>

      <div className="about__section row">
        <div className="col-3 about__section-header">
          <h3 className="about__section-title">Support</h3>
        </div>
        <div className="col-9 about__section-content">
          <p className="about__text">We gratefully acknowledge support from:</p>
          <ul className="about__list">
            <li className="about__list-item">Google</li>
            <li className="about__list-item">Relay Therapeutics</li>
            <li className="about__list-item">Schmidt Futures</li>
            <li className="about__list-item">University of Notre Dame</li>
          </ul>
        </div>
      </div>
    </div>
  );
};

export default About;