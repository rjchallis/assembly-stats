Assembly.prototype.squarePlot = function(parent_div, scale_type, max_count, max_span) {

  // setup plot dimensions
  var size = 600;
  var margin = {top: 40, right: 40, bottom: 50, left: 50},
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

var mover = document.createEvent('UIEvents');
mover.initUIEvent('mouseover', true, true, window, 1);
var mout = document.createEvent('UIEvents');
mout.initUIEvent('mouseout', true, true, window, 1);
var mclick = document.createEvent('UIEvents');
mclick.initUIEvent('click', true, true, window, 1);
var mmover = document.createEvent('UIEvents');
mmover.initUIEvent('mouseover', true, true, window, 1);
var mmout = document.createEvent('UIEvents');
mmout.initUIEvent('mouseout', true, true, window, 1);
var mmclick = document.createEvent('UIEvents');
mmclick.initUIEvent('click', true, true, window, 1);

Assembly.prototype.addKey = function(assemblies) {
  var current = this;
  var plot_area = this.plot_area
  var group = plot_area.append('g').attr('class','asm-square-key');
  var height = assemblies.length * 25;
  var width = 200;
  group.attr('transform','translate(300,'+(430-height)+')')
  var rows = group.selectAll('g').data(assemblies);
  var rows_enter = rows.enter().append('g').attr('class',function(d,i){ var css = 'asm-square-key-row '; if (i == 0){css += 'asm-reference-key'} return css})
                        .attr('transform',function(d,i){return 'translate(0,'+(i*25)+')'})
  var rect = rows_enter.append('rect').attr('height',25).attr('width',150).attr('y',-12.5).style('fill','rgba(255,255,255,0.01)');
  rows_enter.append('text').attr('class',function(d){return 'asm-square-key-text '+d});
  rows.select('text').text(function(d){return d})
  rect.on('mouseover',function(d){
    d = d.replace('.','_');
    plot_area.select('text.'+d).classed('asm-square-focus',true)
    plot_area.select('path.'+d).node().dispatchEvent(mover);
  })
  rect.on('click',function(d){
    d = d.replace('.','_');
    plot_area.select('text.'+d).classed('asm-square-focus',false)
    plot_area.select('path.'+d).node().dispatchEvent(mclick);
  })
  rect.on('mouseout',function(d){
    d = d.replace('.','_');
    plot_area.select('text.'+d).classed('asm-square-focus',false)
    plot_area.select('path.'+d).node().dispatchEvent(mout);
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
    line.xScale = reference.xScale;
    line.yScale = reference.yScale;
    line.max_span = reference.max_span;
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

Assembly.prototype.span_tip = d3.tip()
  .attr('class', 'd3-tip')
  .offset([-5,2])
  .html(function(d) {
    return '<span class="cu-feat-tip"><strong>Span:</strong> ' + (d.span.toLocaleString()) + '<br/><strong>Count:</strong> ' + (d.count.toLocaleString()) + '<br/></span>';
  })

Assembly.prototype.n50_tip = d3.tip()
  .attr('class', 'd3-tip')
  .offset([-5,2])
  .html(function(d) {
    return '<span class="cu-feat-tip"><strong>N50 length:</strong> ' + (d.n50_length.toLocaleString()) + '<br/><strong>N50 number:</strong> ' + (d.n50_count.toLocaleString()) + '<br/></span>';
  })

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
  var n50_point = this.plot_area.append("circle")
                            .attr("r", 5)
                            .attr("cx",function(d){return xScale(assembly.npct_count[499])})
                            .attr("cy",function(d){return yScale(assembly.assembly/2)})
                            .attr("class", 'asm-cumulative-point '+classname)
                            .attr("rel", classname)
                            .style('opacity',0)
  var span_point = this.plot_area.append("circle")
                            .attr("r", 5)
                            .attr("cx",function(d){return xScale(assembly.scaffold_count)})
                            .attr("cy",function(d){return yScale(assembly.assembly)})
                            .attr("class", 'asm-cumulative-point '+classname)
                            .attr("rel", classname)
                            .style('opacity',0)
  span_point.call(this.span_tip)
  span_point.on('mouseover',function(){
    var d = d3.select(this).attr('rel')
    var obj = findByValue(current.line_data,'name',d)
    current.span_tip.show(obj)
  })
  span_point.on('mouseout',function(){
    var d = d3.select(this).attr('rel')
    current.span_tip.hide()
  })
  n50_point.call(this.n50_tip)
  n50_point.on('mouseover',function(){
    var d = d3.select(this).attr('rel')
    var obj = findByValue(current.line_data,'name',d)
    current.n50_tip.show(obj)
  })
  n50_point.on('mouseout',function(){
    var d = d3.select(this).attr('rel')
    current.n50_tip.hide()
  })
  lineGraph.on('mouseover',function(){
    var d = d3.select(this).attr('rel')
    plot_area.select('text.'+d).classed('asm-square-focus',true)
    plot_area.select('path.'+d).classed('asm-square-focus',true)
    plot_area.selectAll('circle.'+d)[0].forEach(function(c){c.dispatchEvent(mmover)});
  })
  lineGraph.on('click',function(){
    var d = d3.select(this).attr('rel')
    var obj = findByValue(current.line_data,'name',d)
    current.zoomTo(obj)
    plot_area.select('text.'+d).classed('asm-square-focus',false)
    plot_area.select('path.'+d).classed('asm-square-focus',false)
    plot_area.selectAll('circle.'+d)[0].forEach(function(c){c.dispatchEvent(mmout)});
  })
  lineGraph.on('mouseout',function(){
    var d = d3.select(this).attr('rel')
    plot_area.select('text.'+d).classed('asm-square-focus',false)
    plot_area.select('path.'+d).classed('asm-square-focus',false)
    plot_area.selectAll('circle.'+d)[0].forEach(function(c){c.dispatchEvent(mmout)});
  })
  this.line_data.push({name:classname,data:data,line:lineGraph,n50_point:n50_point,span_point:span_point,count:assembly.scaffold_count,span:assembly.assembly,n50_length:assembly.npct_length[499],n50_count:assembly.npct_count[499]})

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
      line.xScale = current.xScale;
      line.yScale = current.yScale;
      reference.rescaleLine(line);
    })
  }
  else {
    lineGraph.transition().duration(500).style('opacity',1)
    n50_point.transition().duration(500).style('opacity',1)
    span_point.transition().duration(500).style('opacity',1)
  }
}

Assembly.prototype.rescaleLine = function(line_data) {
  line_data.line.transition().duration(500).attr("d", this.lineFunction(line_data.data)).style('opacity',1)
  line_data.n50_point.transition().duration(500).attr("cx", this.xScale(line_data.n50_count)).attr("cy", this.yScale(line_data.span/2)).style('opacity',1)
  line_data.span_point.transition().duration(500).attr("cx", this.xScale(line_data.count)).attr("cy", this.yScale(line_data.span)).style('opacity',1)
}
