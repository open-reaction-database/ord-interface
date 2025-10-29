import React, { useEffect, useRef, useState, useCallback } from 'react';
import * as d3 from 'd3';
import LoadingSpinner from '../../components/LoadingSpinner';
import reaction_pb from 'ord-schema';
import './ChartView.scss';

interface ChartData {
  smiles: string;
  times_appearing: number;
}

interface ChartViewProps {
  uniqueId: string;
  title: string;
  apiCall: string;
  role: string;
  isCollapsed?: boolean;
}

const ChartView: React.FC<ChartViewProps> = ({
  uniqueId,
  title,
  apiCall,
  role,
  isCollapsed = false
}) => {
  const [loading, setLoading] = useState(true);
  const [inputsData, setInputsData] = useState<ChartData[]>([]);
  const [showTooltip, setShowTooltip] = useState<'visible' | 'hidden'>('hidden');
  const [currentTimesAppearing, setCurrentTimesAppearing] = useState(0);
  const [tooltipOffsetHorizontal, setTooltipOffsetHorizontal] = useState(0);
  const [tooltipOffsetVertical, setTooltipOffsetVertical] = useState(0);
  const [showSmiles, setShowSmiles] = useState(false);
  const [molHtml, setMolHtml] = useState<string | null>(null);
  const [molLoading, setMolLoading] = useState(true);
  
  const svgRef = useRef<SVGSVGElement>(null);

  const getMolHtml = useCallback(async (smiles: string): Promise<string> => {
    const compound = new reaction_pb.Compound();
    const identifier = compound.addIdentifiers();
    identifier.setValue(smiles);
    identifier.setType(reaction_pb.CompoundIdentifier.CompoundIdentifierType.SMILES);

    const binary = compound.serializeBinary();
    
    const response = await fetch('/api/compound_svg', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/x-protobuf',
      },
      body: binary as BodyInit,
    });
    
    return response.json();
  }, []);

  const createChart = useCallback((data: ChartData[], width: number, height: number) => {
    if (!svgRef.current) return;

    // Clear previous chart
    d3.select(svgRef.current).selectAll('*').remove();

    const marginTop = 20;
    const marginRight = 20;
    const marginBottom = 30;
    const marginLeft = 40;

    // X-axis
    const x = d3.scaleBand()
      .domain(d3.groupSort(data, ([d]) => -d.times_appearing, (d) => d.smiles))
      .range([marginLeft, width - marginRight])
      .padding(0.1);

    // Y-axis
    const y = d3.scaleLinear()
      .domain([0, d3.max(data, (d) => d.times_appearing) || 0])
      .range([height - marginBottom, marginTop]);

    // SVG setup
    const svg = d3.select(svgRef.current)
      .attr('width', width)
      .attr('height', height)
      .attr('viewBox', [0, 0, width, height])
      .attr('style', 'max-width: 100%; height: auto;');

    // Create bars
    svg.append('g')
      .attr('fill', 'steelblue')
      .selectAll('rect')
      .data(data)
      .join('rect')
      .attr('x', (d) => x(d.smiles) || 0)
      .attr('y', (d) => y(d.times_appearing))
      .attr('height', (d) => y(0) - y(d.times_appearing))
      .attr('width', x.bandwidth())
      .attr('style', 'cursor: pointer;')
      .on('mouseover', (event: MouseEvent, d: ChartData) => {
        setCurrentTimesAppearing(d.times_appearing);
        setTooltipOffsetHorizontal(event.clientX);
        setTooltipOffsetVertical(role === 'product' ? event.clientY - 240 : event.clientY - 140);
        setShowTooltip('visible');
        setShowSmiles(true);
        setMolLoading(true);
        
        getMolHtml(d.smiles).then((result) => {
          setMolLoading(false);
          setMolHtml(result);
        }).catch((error) => {
          console.error('Error fetching molecule SVG:', error);
          setMolLoading(false);
          setMolHtml(null);
        });
      })
      .on('mouseout', () => {
        setShowTooltip('hidden');
        setMolHtml(null);
        setShowSmiles(false);
        setMolLoading(true);
      });

    // X-axis formatting
    svg.append('g')
      .attr('transform', `translate(0,${height - marginBottom})`)
      .call(d3.axisBottom(x).ticks(width / 80).tickSizeOuter(0).tickFormat(() => ''))
      .call((g) => g.append('text')
        .attr('x', width)
        .attr('y', marginBottom - 4)
        .attr('fill', 'currentColor')
        .attr('text-anchor', 'end')
        .text('Molecules (Hover to view) →'));

    // Y-axis formatting
    svg.append('g')
      .attr('transform', `translate(${marginLeft},0)`)
      .call(d3.axisLeft(y).tickFormat((y) => (y as number).toFixed()))
      .call(g => g.select('.domain').remove())
      .call(g => g.append('text')
        .attr('x', -marginLeft)
        .attr('y', 10)
        .attr('fill', 'currentColor')
        .attr('text-anchor', 'start')
        .text('↑ Frequency (no. of occurrences)'));
  }, [getMolHtml, role]);

  const resize = useCallback(() => {
    const width = isCollapsed ? 180 : 400;
    const height = isCollapsed ? 180 : 400;
    createChart(inputsData, width, height);
  }, [inputsData, isCollapsed, createChart]);

  // Fetch data on mount
  useEffect(() => {
    const datasetId = window.location.pathname.split('/')[2];
    
    fetch(`/api/${apiCall}?dataset_id=${datasetId}`, { method: 'GET' })
      .then(response => response.json())
      .then((data: ChartData[]) => {
        setLoading(false);
        setInputsData(data);
        
        const width = isCollapsed ? 180 : 400;
        const height = isCollapsed ? 180 : 400;
        createChart(data, width, height);
      })
      .catch(() => {
        setLoading(false);
      });
  }, [apiCall, isCollapsed, createChart]);

  // Resize chart when isCollapsed changes
  useEffect(() => {
    if (inputsData.length > 0) {
      resize();
    }
  }, [isCollapsed, resize]);

  return (
    <div className="chart-view">
      <div className="chart-view__title-and-chart">
        <span 
          className="chart-view__title"
          style={isCollapsed 
            ? { fontSize: '10pt', width: '150px' } 
            : { fontSize: '14pt' }
          }
        >
          {title}
        </span>
        
        <svg 
          ref={svgRef}
          id={uniqueId}
          style={{ visibility: loading ? 'hidden' : 'visible' }}
        />
      </div>
      
      <div 
        className="chart-view__loading" 
        style={{ visibility: loading ? 'visible' : 'hidden' }}
      >
        <LoadingSpinner />
      </div>
      
      {showSmiles && (
        <div 
          className="chart-view__tooltip"
          style={{
            top: `${tooltipOffsetVertical}px`,
            left: `${tooltipOffsetHorizontal}px`,
            visibility: showTooltip,
          }}
        >
          <pre>Count: {currentTimesAppearing}</pre>
          <div 
            className="chart-view__svg"
            dangerouslySetInnerHTML={{ __html: molHtml || '' }}
          />
          <div 
            className="chart-view__molloading"
            style={{ visibility: molLoading ? 'visible' : 'hidden' }}
          >
            <LoadingSpinner />
          </div>
        </div>
      )}
    </div>
  );
};

export default ChartView;