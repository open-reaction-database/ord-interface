<style>
.tooltip {
    font: sans-serif 12pt;
    background: #eeeeeeee;
    pointer-events: none;
    border-radius: 2px;
    padding: 5px;
    position: absolute;
    top: 0px;
    left: 0px;
    z-index: 1;
}
</style>
<script>

import FloatingModal from "../../components/FloatingModal";
import reaction_pb from "ord-schema"
import * as d3 from "d3";
export default {
    name: 'ChartView',
    props: {
      uniqueId: String,
      title: String,
      apiCall: String,
      role: String
    },
    components: {
      FloatingModal
    },
    data() {
      return {
        datasetId: "",
        inputsData: [],
        showSmiles: false,
        currentSmiles: "(SMILES string was empty.)",
        currentTimesAppearing: 0,
        molImage: null
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
          this.inputsData = data

          // Dimensions and margins
          const width = 400;
          const height = 400;
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
              .on("click", (d) => {
                this.currentSmiles = d.srcElement.__data__.smiles;
                this.currentTimesAppearing = d.srcElement.__data__.times_appearing;
                this.showSmiles = true;
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
                          .text("Molecules (Click each bar to view full molecule SMILES strings) →"));

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
    }
};
</script>
<template lang="pug">
  .div
    .h4 {{this.title}}
    svg(:id='uniqueId')
    FloatingModal(
            v-if='showSmiles'
            title="SMILES"
            @closeModal='showSmiles=false'
    )
      .data
        pre {{currentSmiles}} occurred {{currentTimesAppearing}} {{currentTimesAppearing == 1 ? 'time' : 'times'}} as a {{role}} in this dataset
        .svg(
          v-html='molHtml'
        )
</template>