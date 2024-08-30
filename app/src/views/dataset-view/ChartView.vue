<template>
  <div>
    <h4>Frequency of Reactants</h4>
    <svg></svg>
    <div class="tooltip"></div>
  </div>
</template>
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
import * as d3 from "d3";
export default {
    name: 'ChartView',
    data() {
      return {};
    },
    mounted() {
      const reactants = [
        {id: 'C', timesAppeared: 0},
        {id: 'H2O', timesAppeared: 1},
        {id: 'H2O2', timesAppeared: 2},
        {id: 'HCl', timesAppeared: 0},
        {id: 'CH3', timesAppeared: 3}
      ]

      const testmap = ['C', 'H2O', 'H2O2', 'CH3']
      // Declare the chart dimensions and margins.
      const width = 400;
      const height = 400;
      const marginTop = 20;
      const marginRight = 20;
      const marginBottom = 30;
      const marginLeft = 40;

      // Needs to be updated to handle data that's already been binned.
      // We dont need bins but we need mols on x axis and timesAppeared on y.
      // Bin the data.
      const bins = d3.bin()
          .thresholds(40)
          .value((d) => d.timesAppeared)(reactants);

      // Declare the x (horizontal position) scale.
      // should be scaleOrdinal
      const x = d3.scaleLinear()
          .domain([bins[0].x0, bins[bins.length - 1].x1])
          .range([marginLeft, width - marginRight]);

      // Declare the y (vertical position) scale.
      const y = d3.scaleLinear()
          .domain([0, d3.max(bins, (d) => d.length)])
          .range([height - marginBottom, marginTop]);

      // Add these SVG attrs below
      /*
          .attr("viewBox", [0, 0, width, height])
          .attr("style", "max-width: 100%; height: auto;");*/
      const svg = d3.select("svg").attr("width", width).attr("height", height);

      const tooltip = d3.select(".tooltip");

      // Add a rect for each bin.
      svg.append("g")
          .attr("fill", "steelblue")
        .selectAll()
        .data(bins)
        .join("rect")
          .attr("x", (d) => x(d.x0) + 1)
          .attr("width", (d) => x(d.x1) - x(d.x0) - 1)
          .attr("y", (d) => y(d.length))
          .attr("height", (d) => y(0) - y(d.length))
          .on("mouseenter", (evt, d) => {
            const [mx, my] = d3.pointer(evt);
            const tooltipText = `
              <strong>SMILES</strong>: ${d["id"]} 
              <br> 
              <strong># of Occurrences</strong>: ${d["timesAppeared"]}`;
            tooltip
              .style("top", `${my}px`)
              .style("left", `${mx}px`)
              .html(tooltipText);
          })
          .on("mouseout", () => tooltip.text());

      // Add the x-axis and label.
      svg.append("g")
          .attr("transform", `translate(0,${height - marginBottom})`)
          .call(d3.axisBottom(x).ticks(width / 80).tickSizeOuter(0).tickFormat(function(d) {return testmap[d]}))
          .call((g) => g.append("text")
              .attr("x", width)
              .attr("y", marginBottom - 4)
              .attr("fill", "currentColor")
              .attr("text-anchor", "end")
              .text("Molecules (SMILES) →"));

      // Add the y-axis and label, and remove the domain line.
      svg.append("g")
          .attr("transform", `translate(${marginLeft},0)`)
          .call(d3.axisLeft(y).ticks(height / 40))
          .call((g) => g.select(".domain").remove())
          .call((g) => g.append("text")
              .attr("x", -marginLeft)
              .attr("y", 10)
              .attr("fill", "currentColor")
              .attr("text-anchor", "start")
              .text("↑ Frequency (no. of occurrences)"));
    }
};
</script>