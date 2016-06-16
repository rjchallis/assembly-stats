Assembly.prototype.squarePlot = function(parent_div, scale_type, max_count, max_span) {

  // setup plot dimensions
  var size = 600;
  var margin = {top: 80, right: 20, bottom: 60, left: 20},
    width = size - margin.left - margin.right,
    height = size - margin.top - margin.bottom;
  var tick = 10;
  var w = 12; // coloured box size for legend

  if (!scale_type) scale_type = 'linear'
  if (!max_count) max_count = this.scaffold_count //.toLocaleString()
  if (!max_span) max_span = this.assembly

  var x = d3.scale.linear()
    .domain([1, max_count])
    .range([0, width])
    .nice();

  var y = d3.scale.linear()
    .domain([1, max_span])
    .range([height, 0])
    .nice();


  var xAxis = d3.svg.axis()
    .scale(x)
    .orient("bottom")
      .tickFormat(d3.format("s"));
  this.xAxis = xAxis;
  var yAxis = d3.svg.axis()
    .scale(y)
    .orient("left")
      .tickFormat(d3.format("s"));
  this.yAxis = yAxis;


  this.parent_div = parent_div;
  var parent = d3.select('#' + parent_div);
  var svg = parent.append('svg');


  svg.attr('width', '100%')
    .attr('height', '100%')
    .attr('viewBox', '0 0 ' + size + ' ' + size)
    .attr('preserveAspectRatio', 'xMidYMid meet')
  var plot_area = svg.append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  svg.append("g")
    .attr("class", "x axis")
    .attr("transform", "translate(" + margin.left + "," + (height + margin.top) + ")")
    .call(xAxis)
    .selectAll("text")
    .attr("y", 0)
    .attr("x", 9)
    .attr("dy", ".35em")
    .attr("transform", "rotate(90)")
    .style("text-anchor", "start");

  svg.append("g")
    .attr("class", "y axis")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
    .call(yAxis);

  this.svg = svg;
  this.plot_area = plot_area;
  this.xScale = x;
  this.yScale = y;
  this.max_count = max_count;
  this.max_span = max_span;
  this.line_data = [];
  var data = this.prepareLine()
  this.addLine(data,this.name,this);
  this.line_data[0].line.classed('asm-reference-line',true).transition().duration(500).style('opacity',1)



  return this;
}

Assembly.prototype.addKey = function(assemblies) {
  var current = this;
  var plot_area = this.plot_area
  var group = plot_area.append('g').attr('class','asm-square-key');
  var height = assemblies.length * 25;
  var width = 200;
  group.attr('transform','translate(370,'+(430-height)+')')
  var rows = group.selectAll('g').data(assemblies);
  var rows_enter = rows.enter().append('g').attr('class',function(d,i){ var css = 'asm-square-key-row '; if (i == 0){css += 'asm-reference-key'} return css})
                        .attr('transform',function(d,i){return 'translate(0,'+(i*25)+')'})
  var rect = rows_enter.append('rect').attr('height',25).attr('width',150).attr('y',-12.5).style('fill','rgba(255,255,255,0.01)');
  rows_enter.append('text').attr('class',function(d){return 'asm-square-key-text '+d});
  rows.select('text').text(function(d){return d})
  rect.on('mouseover',function(d){
    //d3.select(this).style('fill','rgba(0,0,0,0.3)')
    plot_area.select('text.'+d).classed('asm-square-focus',true)
    plot_area.select('path.'+d).classed('asm-square-focus',true)
  })
  rect.on('click',function(d){
    //d3.select(this).style('fill','rgba(0,0,0,0.3)')
    var obj = findByValue(current.line_data,'name',d)
    current.zoomTo(obj)
    plot_area.select('text.'+d).classed('asm-square-focus',false)
    plot_area.select('path.'+d).classed('asm-square-focus',false)
  })
  rect.on('mouseout',function(d){
    //d3.select(this).style('fill','rgba(0,0,0,0.01)')
    plot_area.select('text.'+d).classed('asm-square-focus',false)
    plot_area.select('path.'+d).classed('asm-square-focus',false)
  })
  return this;
}

function findByValue(source, key, value) {
  for (var i = 0; i < source.length; i++) {
    if (source[i][key] === value) {
      return source[i];
    }
  }
}

Assembly.prototype.zoomTo = function(obj){
  var xScale = this.xScale;
  var yScale = this.yScale;
  xScale.domain([1,obj.count]).nice()
  this.max_count = obj.count;
  yScale.domain([1,obj.span]).nice()
  this.max_span = obj.span;
  this.svg.selectAll("g.y.axis").transition().duration(500)
      .call(this.yAxis);
  this.svg.selectAll("g.x.axis")
      .transition().duration(500).call(this.xAxis)
      .selectAll("text")
      .attr("y", 0)
      .attr("x", 9)
      .attr("dy", ".35em")
      .attr("transform", "rotate(90)")
      .style("text-anchor", "start");
  var line_data = this.line_data;
  var reference = this;
  line_data.forEach(function(line){
    reference.rescaleLine(line);
  })
}

Assembly.prototype.prepareLine = function() {
  var current = this;
  var line_data = [];
  current.npct_count.forEach(function(count,index){
    line_data[index] = {x:count,y:((index+1)/1000*current.assembly)}
  })
  return line_data;
}

Assembly.prototype.addLine = function(data,classname,assembly) {
  var xScale = this.xScale;
  var yScale = this.yScale;
  var plot_area = this.plot_area;
  var current = this;
  var rescale = 0;

  var lineFunction = d3.svg.line()
                            .x(function(d) { return xScale(d.x); })
                            .y(function(d) { return yScale(d.y); })
                            .interpolate("linear");
  this.lineFunction = lineFunction;
  var lineGraph = this.plot_area.append("path")
                              .attr("d", lineFunction(data))
                              .attr("class", 'asm-cumulative-line '+classname)
                              .attr("rel", classname)
                              .style('opacity',0)
  lineGraph.on('mouseover',function(){
    var d = d3.select(this).attr('rel')
    plot_area.select('text.'+d).classed('asm-square-focus',true)
    plot_area.select('path.'+d).classed('asm-square-focus',true)
  })
  lineGraph.on('click',function(){
    var d = d3.select(this).attr('rel')
    var obj = findByValue(current.line_data,'name',d)
    current.zoomTo(obj)
    plot_area.select('text.'+d).classed('asm-square-focus',false)
    plot_area.select('path.'+d).classed('asm-square-focus',false)
  })
  lineGraph.on('mouseout',function(){
    var d = d3.select(this).attr('rel')
    plot_area.select('text.'+d).classed('asm-square-focus',false)
    plot_area.select('path.'+d).classed('asm-square-focus',false)
  })
  this.line_data.push({name:classname,data:data,line:lineGraph,count:assembly.scaffold_count,span:assembly.assembly})

  if (assembly.scaffold_count > this.max_count){
    rescale = 1;
    xScale.domain([1,assembly.scaffold_count]).nice()
    this.max_count = assembly.scaffold_count;
  }
  if (assembly.assembly > this.max_span){
    rescale = 1
    yScale.domain([1,assembly.assembly]).nice()
    this.max_span = assembly.assembly;
  }
  if (rescale){
    this.svg.selectAll("g.y.axis").transition().duration(500)
        .call(this.yAxis);
    this.svg.selectAll("g.x.axis")
        .transition().duration(500).call(this.xAxis)
        .selectAll("text")
        .attr("y", 0)
        .attr("x", 9)
        .attr("dy", ".35em")
        .attr("transform", "rotate(90)")
        .style("text-anchor", "start");
    var line_data = this.line_data;
    var reference = this;
    line_data.forEach(function(line){
      reference.rescaleLine(line);
    })
  }
  else {
    lineGraph.transition().duration(500).style('opacity',1)
  }
}

Assembly.prototype.rescaleLine = function(line_data) {
  line_data.line.transition().duration(500).attr("d", this.lineFunction(line_data.data)).style('opacity',1)
}
