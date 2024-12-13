<!--
 Copyright 2024 Open Reaction Database Project Authors

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

<style>
.tooltip {
    font: sans-serif 12pt;
    background: #ffffff;
    pointer-events: none;
    border: #000 solid;
    border-radius: 4px;
    padding: 5px;
    position: absolute;
    top: 0px;
    left: 0px;
    z-index: 1;
    visibility: hidden;
    opacity: 1;
}
</style>
<script>

import FloatingModal from "../../components/FloatingModal";
import LoadingSpinner from '@/components/LoadingSpinner';
import reaction_pb from "ord-schema";
import * as d3 from "d3";

export default {
    name: 'ChartView',
    props: {
      uniqueId: String,
      title: String,
      apiCall: String,
      role: String,
      isCollapsed: Boolean
    },
    components: {
      FloatingModal,
      LoadingSpinner
    },
    watch: {
      currentSmiles() {
        this.getMolHtml()
      },
      isCollapsed () {
        this.resize()
      }
   },
    data() {
      return {
        datasetId: "",
        inputsData: [],
        showSmiles: false,
        currentSmiles: "(SMILES string was empty.)",
        currentTimesAppearing: 0,
        molImage: null,
        loading: true,
        molHtml: null,
        tooltipOffsetHorizontal: 0,
        tooltipOffsetVertical: 0,
        showTooltip: "hidden",
        molloading: true
      };
    },
    /*computed: {
      compoundSVG: function () {
        // Add SVG of the molecule to modal
        // prep compound
        const compound = new reaction_pb.Compound()
        const identifier = compound.addIdentifiers()
        identifier.setValue(this.currentSmiles)
        identifier.setType(reaction_pb.CompoundIdentifier.CompoundIdentifierType.SMILES)
        // send http request
        return new Promise(resolve => {
          const xhr = new XMLHttpRequest()
          xhr.open("POST", "/api/compound_svg")
          const binary = compound.serializeBinary()
          xhr.responseType = "json"
          xhr.onload = function () {
            resolve(xhr.response)
          }
          xhr.send(binary)
        }).then(val => {
          this.molImage = val
        })
      }
    },*/
    mounted() {
      this.datasetId = window.location.pathname.split("/")[2]
      fetch("/api/" + this.apiCall + "?dataset_id=" + this.datasetId, {method: "GET"})
        .then(response => response.json())
        .then(data => {
          this.loading = false
          this.inputsData = data

          // Dimensions and margins
          const width = this.isCollapsed ? 180 : 400; // 400
          const height = this.isCollapsed ? 180 : 400; // 400
          const marginTop = 20;
          const marginRight = 20;
          const marginBottom = 30;
          const marginLeft = 40;

          // X-axis
          const x = d3.scaleBand()
              .domain(d3.groupSort(data, ([d]) => -d.times_appearing, (d) => d.smiles)) // descending frequency
              .range([marginLeft, width - marginRight])
              .padding(0.1);
          
          // Y-axis
          const y = d3.scaleLinear()
              .domain([0, d3.max(data, (d) => d.times_appearing)])
              .range([height - marginBottom, marginTop]);

          // SVG
          const svg = d3.select("#" + this.uniqueId)
              .attr("width", width)
              .attr("height", height)
              .attr("viewBox", [0, 0, width, height])
              .attr("style", "max-width: 100%; height: auto;");

          // SVG rects for bars
          svg.append("g")
              .attr("fill", "steelblue")
            .selectAll()
            .data(data)
            .join("rect")
              .attr("x", (d) => x(d.smiles))
              .attr("y", (d) => y(d.times_appearing))
              .attr("height", (d) => y(0) - y(d.times_appearing))
              .attr("width", x.bandwidth())
              .attr("style", "cursor:pointer;")
              .on("mouseover", (d) => {
                this.currentSmiles = d.srcElement.__data__.smiles;
                this.currentTimesAppearing = d.srcElement.__data__.times_appearing;
                this.tooltipOffsetHorizontal = d.clientX;
                this.tooltipOffsetVertical = this.role == "product" ? d.clientY-240 : d.clientY-140;
                this.showTooltip = "visible";
                this.showSmiles = true;
              })
              .on("mouseout", () => {
                this.showTooltip = "hidden";
                this.molHtml = null;
                this.showSmiles = false;
              });

          // X-axis formatting
          svg.append("g")
              .attr("transform", `translate(0,${height - marginBottom})`)
              .call(d3.axisBottom(x).ticks(width / 80).tickSizeOuter(0).tickFormat(function(d) {return ''}))
              .call((g) => g.append("text")
                          .attr("x", width)
                          .attr("y", marginBottom - 4)
                          .attr("fill", "currentColor")
                          .attr("text-anchor", "end")
                          .text("Molecules (Hover to view) →"));

          // Y-axis formatting
          svg.append("g")
              .attr("transform", `translate(${marginLeft},0)`)
              .call(d3.axisLeft(y).tickFormat((y) => (y).toFixed()))
              .call(g => g.select(".domain").remove())
              .call(g => g.append("text")
                  .attr("x", -marginLeft)
                  .attr("y", 10)
                  .attr("fill", "currentColor")
                  .attr("text-anchor", "start")
                  .text("↑ Frequency (no. of occurrences)"));
    })
    },
    methods: {
      getMolHtml () {
        const compoundStr = this.currentSmiles
        // prep compound
        const compound = new reaction_pb.Compound()
        const identifier = compound.addIdentifiers()
        identifier.setValue(compoundStr)
        identifier.setType(reaction_pb.CompoundIdentifier.CompoundIdentifierType.SMILES)
        // send http request
        return new Promise(resolve => {
          const xhr = new XMLHttpRequest()
          xhr.open("POST", "/api/compound_svg")
          const binary = compound.serializeBinary()
          xhr.responseType = "json"
          xhr.onload = function () {
            resolve(xhr.response)
          }
          xhr.send(binary)
        }).then(val => {
          this.molloading = false
          this.molHtml = val
        })
      },
      resize () {
        // Clear before attempting resize - everything needs to be recalculated based on new dimensions
        d3.select("#" + this.uniqueId).selectAll("*").remove()

        this.datasetId = window.location.pathname.split("/")[2]
        // this.inputsData = data

        // Dimensions and margins
        const width = this.isCollapsed ? 180 : 400; // 400
        const height = this.isCollapsed ? 180 : 400; // 400
        const marginTop = 20;
        const marginRight = 20;
        const marginBottom = 30;
        const marginLeft = 40;

        // X-axis
        const x = d3.scaleBand()
            .domain(d3.groupSort(this.inputsData, ([d]) => -d.times_appearing, (d) => d.smiles)) // descending frequency
            .range([marginLeft, width - marginRight])
            .padding(0.1);
        
        // Y-axis
        const y = d3.scaleLinear()
            .domain([0, d3.max(this.inputsData, (d) => d.times_appearing)])
            .range([height - marginBottom, marginTop]);

        // SVG
        const svg = d3.select("#" + this.uniqueId)
            .attr("width", width)
            .attr("height", height)
            .attr("viewBox", [0, 0, width, height])
            .attr("style", "max-width: 100%; height: auto;");

        // SVG rects for bars
        svg.append("g")
            .attr("fill", "steelblue")
          .selectAll()
          .data(this.inputsData)
          .join("rect")
            .attr("x", (d) => x(d.smiles))
            .attr("y", (d) => y(d.times_appearing))
            .attr("height", (d) => y(0) - y(d.times_appearing))
            .attr("width", x.bandwidth())
            .attr("style", "cursor:pointer;")
            .on("mouseover", (d) => {
              this.currentSmiles = d.srcElement.__data__.smiles;
              this.currentTimesAppearing = d.srcElement.__data__.times_appearing;
              this.tooltipOffsetHorizontal = d.clientX;
              this.tooltipOffsetVertical = !this.isCollapsed || this.role == "product" ? d.clientY-240 : d.clientY-140;
              this.showTooltip = "visible";
              this.showSmiles = true;
            })
            .on("mouseout", () => {
              this.showTooltip = "hidden";
              this.molHtml = null;
              this.molloading = true;
              this.showSmiles = false;
            });

        // X-axis formatting
        svg.append("g")
            .attr("transform", `translate(0,${height - marginBottom})`)
            .call(d3.axisBottom(x).ticks(width / 80).tickSizeOuter(0).tickFormat(function(d) {return ''}))
            .call((g) => g.append("text")
                        .attr("x", width)
                        .attr("y", marginBottom - 4)
                        .attr("fill", "currentColor")
                        .attr("text-anchor", "end")
                        .text("Molecules (Hover to view) →"));

        // Y-axis formatting
        svg.append("g")
            .attr("transform", `translate(${marginLeft},0)`)
            .call(d3.axisLeft(y).tickFormat((y) => (y).toFixed()))
            .call(g => g.select(".domain").remove())
            .call(g => g.append("text")
                .attr("x", -marginLeft)
                .attr("y", 10)
                .attr("fill", "currentColor")
                .attr("text-anchor", "start")
                .text("↑ Frequency (no. of occurrences)"));
        
      }
  },
};
</script>
<template lang="pug">
  .div
    .span(:style='isCollapsed ? "font-size: 10pt; width: 150px;" : "font-size: 14pt;"') {{this.title}}
    svg(:id='uniqueId' :style='!loading ? "visibility: visible" : "visibility: hidden"')
    .loading(:style='loading ? "visibility: visible" : "visibility: hidden"')
      LoadingSpinner
    .tooltip(:style='"top: " + tooltipOffsetVertical + "px; left: " + tooltipOffsetHorizontal + "px; visibility: " + showTooltip' v-if='showSmiles')
        pre Count: {{currentTimesAppearing}}
        .svg(
          v-html='molHtml'
        )
        .molloading(:style='this.molloading ? "visibility: visible" : "visibility: hidden"')
          LoadingSpinner
</template>