<!--
 Copyright 2023 Open Reaction Database Project Authors

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
-->

<script>
export default {
  name: 'AboutView',
  components: {
  }
}
</script>

<template>
  <div id="about" class="container">
    <div class="row mt-5">
        <div class="col-3 border-top pt-3">
            <h3>About</h3>
        </div>
        <div class="col-9 p-3">
            <p style="text-align: justify; text-justify: inter-word;">
                The Open Reaction Database (ORD) is an open-access schema and infrastructure for structuring and
                sharing organic reaction data, including a centralized data repository. The ORD schema supports
                conventional and emerging technologies, from benchtop reactions to automated high-throughput
                experiments and flow chemistry. Our vision is that a consistent data representation and
                infrastructure to support data sharing will enable downstream applications that will greatly improve
                the state of the art with respect to computer-aided synthesis planning, reaction prediction, and
                other predictive chemistry tasks.
            </p>
            <p style="text-align: justify; text-justify: inter-word;">
                Since our initial meeting in October 2019, the database has grown to include more than 2M reactions
                (including a large dataset of reactions extracted from USPTO sources) and received contributions
                from academic and industrial users, both from published and unpublished work. Some of our current
                efforts include:
            </p>
            <ul>
                <li>Improving user interfaces and providing support to contributors on GitHub and via email.</li>
                <li>
                    Working with companies to incorporate the ORD schema into their reaction informatics
                    infrastructure, including the development of "translators" between the ORD schema and electronic
                    lab notebooks (ELNs).
                </li>
                <li>
                    Engaging with journals and other stakeholders to drive adoption of the ORD schema as a <a
                        href="https://en.wikipedia.org/wiki/FAIR_data">FAIR</a> data structure for sharing reaction
                    data across academia, government, and industry.
                </li>
            </ul>
            <p>
                <a href="https://www.linkedin.com/in/benjamin-deadman-78786242/">Ben Deadman</a>
                is our Reaction Data Evangelist and ORD Program Manager. Please reach out to him at
                <a href="mailto:help@open-reaction-database.org">help@open&#8209;reaction&#8209;database.org</a>
                for help preparing a contribution or to discuss using the ORD in your company or lab.
            </p>
        </div>
    </div>
    <div class="row mt-5">
        <div class="col-3 border-top pt-3">
            <h3>Publications and Media</h3>
        </div>
        <div class="col-9 p-3">
            <h5>Journal Articles</h5>
            <ul>
                <li>
                    Kearnes SM, Maser MR, Wleklinski M, Kast A, Doyle AG, Dreher SD, Hawkins JM, Jensen KF, Coley
                    CW. The Open Reaction Database. <em>J Am Chem Soc</em> 2021, 143(45), 18820-18826. (<a
                        href="https://doi.org/10.1021/jacs.1c09820">JACS</a>)
                </li>
                <li>
                    Mercado R, Kearnes SM, Coley C. Data Sharing in Chemistry: Lessons Learned and a Case for
                    Mandating Structured Reaction Data. <em>J Chem Inf Model</em> 2023, 63(14), 4253-4265. (<a
                        href="https://doi.org/10.1021/acs.jcim.3c00607">JCIM</a>)
                </li>
            </ul>
            <h5>News</h5>
            <ul>
                <li>
                    A new database for machine-learning research (<a
                        href="https://cen.acs.org/physical-chemistry/computational-chemistry/new-database-machine-learning-research/99/web/2021/11">C&EN</a>,
                    22 November 2021)
                </li>
                <li>
                    Yield-predicting AI needs chemists to stop ignoring failed experiments (<a
                        href="https://www.chemistryworld.com/news/yield-predicting-ai-needs-chemists-to-stop-ignoring-failed-experiments/4015662.article">Chemistry
                        World</a>, 12 May 2022)
                </li>
                <li>
                    Chemists debate machine learning&#39;s future in synthesis planning and ask for open data (<a
                        href="https://cen.acs.org/physical-chemistry/computational-chemistry/Chemists-debate-machine-learnings-future/100/i18">C&EN</a>,
                    18 May 2022)
                </li>
                <li>
                    For chemists, the AI revolution has yet to happen (<a
                        href="https://www.nature.com/articles/d41586-023-01612-x">Nature</a>, 17 May 2023)
                </li>
            </ul>
        </div>
    </div>
    <div class="row mt-5">
        <div class="col-3 border-top pt-3">
            <h3>Leadership</h3>
        </div>
        <div class="col-9 p-3">
            <h5>Governing Committee</h5>
            <div class="row">
                <div class="col-6">
                    <ul>
                        <li>Connor Coley (MIT, C-CAS)</li>
                        <li>Abby Doyle (UCLA, C-CAS)</li>
                        <li>Spencer Dreher (Merck)</li>
                    </ul>
                </div>
                <div class="col-6">
                    <ul>
                        <li>Joel Hawkins (Pfizer)</li>
                        <li>Klavs Jensen (MIT)</li>
                        <li>Steven Kearnes (Relay)</li>
                    </ul>
                </div>
            </div>
            <h5>Advisory Board</h5>
            <div class="row">
                <div class="col-6">
                    <ul>
                        <li>Al&aacute;n Aspuru-Guzik (Toronto, MADNESS)</li>
                        <li>Timothy Cernak (Michigan, Entos)</li>
                        <li>Lucy Colwell (Cambridge, SynTech, Google)</li>
                        <li>Werngard Czechtizky (AstraZeneca)</li>
                        <li>JW Feng</li>
                        <li>Matthew Gaunt (Cambridge, SynTech)</li>
                        <li>Alex Godfrey (NCATS Consultant)</li>
                        <li>Mimi Hii (Imperial, ROAR)</li>
                        <li>Greg Landrum (T5 Informatics)</li>
                    </ul>
                </div>
                <div class="col-6">
                    <ul>
                        <li>Fabio Lima (Novartis)</li>
                        <li>Christos Nicolaou (Recursion)</li>
                        <li>Sarah Reisman (Caltech, C-CAS)</li>
                        <li>Francesco Rianjongdee (GSK)</li>
                        <li>Marwin Segler (Microsoft)</li>
                        <li>Matthew Sigman (Utah, C-CAS)</li>
                        <li>Jay Stevens (BMS)</li>
                        <li>Sarah Trice (XtalPi)</li>
                        <li>Huimin Zhao (UIUC, MMLI)</li>
                    </ul>
                </div>
            </div>
        </div>
    </div>
    <div class="row mt-5">
        <div class="col-3 border-top pt-3">
            <h3>Support</h3>
        </div>
        <div class="col-9 p-3">
            <p>We gratefully acknowledge support from:</p>
            <ul>
                <li>Google</li>
                <li>Relay Therapeutics</li>
                <li>Schmidt Futures</li>
                <li>University of Notre Dame</li>
            </ul>
        </div>
    </div>
  </div>
</template>

<style lang="sass" scoped>
#about
  margin-top: 2rem
  background-color: white
  padding: 0 2rem
  border-radius: 0.5rem
</style>
