

 function getReadableSeqSizeString(seqSizeInBases,fixed) {
//http://stackoverflow.com/questions/10420352/converting-file-size-in-bytes-to-human-readable

    var i = -1;
    var baseUnits = [' kB', ' MB', ' GB', ' TB', 'PB', 'EB', 'ZB', 'YB'];
    do {
        seqSizeInBases = seqSizeInBases / 1000;
        i++;
    } while (seqSizeInBases >= 1000);
	fixed = fixed ? fixed : fixed == 0 ? 0 : 1;
    return Math.max(seqSizeInBases, 0.1).toFixed(fixed) + baseUnits[i];
};

function Assembly( stats,scaffolds,contigs ) {
  stats.scaffolds = scaffolds ? scaffolds : stats.scaffolds;
  stats.contigs = contigs ? contigs : stats.contigs;
  var sum = stats.scaffolds.reduce(function(previousValue, currentValue, index, array) {
    return previousValue + currentValue;
  });
  if (stats.contigs){
  	var ctgsum = stats.contigs.reduce(function(previousValue, currentValue, index, array) {
  	  return previousValue + currentValue;
  	});
  }
  this.genome = stats.genome;
  this.assembly = stats.assembly;
  this.N = stats.N ? stats.N <= 100 ? stats.N < 1 ? stats.N * 100 : stats.N : stats.N / this.assembly * 100 : 0;
  this.ATGC = stats.ATGC ? stats.ATGC <= 100 ? stats.ATGC <= 1 ? stats.ATGC * 100 : stats.ATGC : stats.ATGC / this.assembly * 100 : 100 - this.N;
  this.GC = stats.GC ? stats.GC <= 100 ? stats.GC <= 1 ? stats.GC * 100 : stats.GC : stats.GC / this.assembly * 100 : 0;
  this.cegma_complete = stats.cegma_complete;
  this.cegma_partial = stats.cegma_partial;
  this.scaffolds = stats.scaffolds.sort(function(a, b){return b-a});
  var npct_length = {};
  var npct_count = {};
  this.GCs = stats.GCs;
  this.Ns = stats.Ns;
  function getRandomArbitrary(min, max) {
    return Math.random() * (max - min) + min;
  }
  var lsum = 0;
  this.scaffolds.forEach(function(length,index,array){
	var new_sum = lsum + length;
	if (Math.floor(new_sum/sum*1000) > Math.floor(lsum/sum*100)){
		npct_length[Math.floor(new_sum/sum*1000)] = length;
		npct_count[Math.floor(new_sum/sum*1000)] = index+1;
	}
	lsum = new_sum;
  });
  this.seq = Array.apply(0, Array(1000)).map(function (x, y) { return 1000 - y; });
  this.seq.forEach(function(i,index){
  	if (!npct_length[i]) npct_length[i] = npct_length[(i+1)];
  	if (!npct_count[i]) npct_count[i] = npct_count[(i+1)];
  });
  this.npct_length = npct_length;
  this.npct_count = npct_count;

  var nctg_length = {};
  var nctg_count = {};

  if (stats.contigs){
    this.contigs = stats.contigs.sort(function(a, b){return b-a});

  var lsum = 0;
  this.contigs.forEach(function(length,index,array){
	var new_sum = lsum + length;
	if (Math.floor(new_sum/ctgsum*1000) > Math.floor(lsum/ctgsum*100)){
		nctg_length[Math.floor(new_sum/ctgsum*1000)] = length;
		nctg_count[Math.floor(new_sum/ctgsum*1000)] = index+1;
	}
	lsum = new_sum;
  });
  this.seq.forEach(function(i,index){
  	if (!nctg_length[i]) nctg_length[i] = nctg_length[(i+1)];
  	if (!nctg_count[i]) nctg_count[i] = nctg_count[(i+1)];
  });
  this.nctg_length = nctg_length;
  this.nctg_count = nctg_count;
  }

  this.scale = {};
  this.setScale('percent','linear',[0,100],[180* (Math.PI/180),90* (Math.PI/180)]);
  this.setScale('100percent','linear',[0,100],[0,2 * Math.PI]);
  this.setScale('gc','linear',[0,100],[0,100]); // range will be updated when drawing
  this.setScale('count','log',[1,1e6],[100,1]); // range will be updated when drawing
  this.setScale('length','sqrt',[1,1e6],[1,100]); // range will be updated range when drawing
}

Assembly.prototype.setScale = function(element,scaling,domain,range){
  this.scale[element] = scaling == 'log' ? d3.scale.log() : scaling == 'sqrt' ? d3.scale.sqrt() : d3.scale.linear();
  this.scale[element].domain(domain);
  this.scale[element].range(range);
}

Assembly.prototype.drawPlot = function(parent,longest,circle_span){

  // setup plot dimensions
  var size = 600;
  var margin = 100;
  var tick = 10;
  var w = 12; // coloured box size for legend

  parent.attr('width', '100%')
  		.attr('height', '100%')
  		.attr('viewBox','0 0 '+size+' '+size)
  		.attr('preserveAspectRatio','xMidYMid meet')

  // setup radii for circular plots
  var radii = {};
  radii.core = [0,(size-margin*2-tick*2)/2];
  radii.core.majorTick = [radii.core[1],radii.core[1]+tick];
  radii.core.minorTick = [radii.core[1],radii.core[1]+tick/2];

  radii.percent = [radii.core[1]+tick*4,radii.core[1]];;
  radii.percent.majorTick = [radii.percent[0],radii.percent[0]-tick];
  radii.percent.minorTick = [radii.percent[0],radii.percent[0]-tick/2];

  radii.genome = [0,tick*3];;

  radii.ceg = [0,tick*2,tick*4];
  radii.ceg.majorTick = [radii.ceg[2],radii.ceg[2]+tick/1.5];
  radii.ceg.minorTick = [radii.ceg[2],radii.ceg[2]+tick/3];

  this.radii = radii;

  // adjust scales for plot dimensions/data
  if (!longest) longest = this.scaffolds[0]+1
  if (longest < this.scaffolds[0]) longest = this.scaffolds[0]+1
  if (!circle_span) circle_span = this.assembly
  if (circle_span < this.assembly) circle_span = this.assembly
  var span_ratio = this.assembly / circle_span * 100;
  this.scale['length'].domain([1,longest])
  this.scale['length'].range([radii.core[0],radii.core[1]])
  this.scale['count'].range([radii.core[1],radii.core[0]+radii.core[1]/3])
  this.scale['percent'].range([0,(2 * Math.PI * this.assembly / circle_span)])
  this.scale['gc'].range([radii.percent[1],radii.percent[0]])

  var lScale = this.scale['length'];
  var cScale = this.scale['count'];
  var pScale = this.scale['percent'];
  var p100Scale = this.scale['100percent'];
  var gScale = this.scale['gc'];
  var npct_length = this.npct_length;
  var npct_count = this.npct_count;
  var nctg_length = this.nctg_length;
  var nctg_count = this.nctg_count;
  var nctg_GC = this.nctg_GC;
  var nctg_N = this.nctg_N;
  var scaffolds = this.scaffolds;
  var contigs = this.contigs;

  // create a group for the plot
  var g = parent.append('g')
      .attr("transform","translate("+size/2+","+size/2+")")
      .attr("id","asm-g-plot");




  // draw base composition axis fill
  var bcg = g.append('g')
      .attr("id","asm-g-base_composition");
  var bcdg = bcg.append('g')
      .attr("id","asm-g-base_composition_data");
  var n = 100 - this.ATGC;
  var gc_start = n / 100 * this.GC;
  if (this.GCs && this.Ns){
    plot_arc(bcdg,radii.percent[0],radii.percent[1],pScale(0),pScale(100),'asm-ns');
  	var lower = [];
    var upper = [];
    var GCs = this.GCs;
    var Ns = this.Ns;
    Ns.forEach(function(current,i){
      lower.push((current / 100 * GCs[i]));
      upper.push((100 - current + lower[i]));
    });
    var line = d3.svg.line()
    .x(function(d,i) { return Math.cos(pScale(i/10)-Math.PI/2)*(gScale(d)); })
    .y(function(d,i) { return Math.sin(pScale(i/10)-Math.PI/2)*(gScale(d)); });

    var line = d3.svg.line()
      .x(function(d,i) { return Math.cos(pScale(i/10)-Math.PI/2)*(gScale(d)); })
      .y(function(d,i) { return Math.sin(pScale(i/10)-Math.PI/2)*(gScale(d)); });
    var revline = d3.svg.line()
      .x(function(d,i) { return Math.cos(pScale((1000-i)/10)-Math.PI/2)*(gScale(d)); })
      .y(function(d,i) { return Math.sin(pScale((1000-i)/10)-Math.PI/2)*(gScale(d)); });
    var atgc = line([0])+'L'+line(lower).replace(/M[^L]+?/,'')+revline(upper).replace('M','L')

  bcdg.append("path")
  	  .attr("class", "asm-atgc")
      .attr("d", atgc)
      .attr("fill-rule","evenodd");
  var gc = line([0])+'L'+line(lower).replace(/M[^L]+?/,'')+revline(GCs).replace('M','L')
  bcdg.append("path")
  	  .attr("class", "asm-gc")
      .attr("d", gc)
      .attr("fill-rule","evenodd");


  }
  else {
  	plot_arc(bcdg,radii.percent[0],radii.percent[1],pScale(0),pScale(100),'asm-ns');
  	plot_arc(bcdg,radii.percent[0],radii.percent[1],pScale(gc_start),pScale(gc_start+this.ATGC),'asm-atgc');
  	plot_arc(bcdg,radii.percent[0],radii.percent[1],pScale(gc_start),pScale(this.GC),'asm-gc');
  }
  var bcag = bcg.append('g')
      .attr("id","asm-g-base_composition_axis");
  percent_axis(bcag,radii,pScale);

 // plot CEGMA completeness if available
  if (this.cegma_complete){
   var ccg = g.append('g')
       .attr('transform','translate('+(radii.percent[1]+tick*3)+','+(-radii.percent[1]-tick*2)+')')
      .attr("id","asm-cegma_completeness");
  	var ccdg = ccg.append('g')
       .attr("id","asm-cegma_completeness_data");
  	 plot_arc(ccdg,radii.ceg[0],radii.ceg[1],p100Scale(0),p100Scale(this.cegma_complete),'asm-ceg_comp');
     plot_arc(ccdg,radii.ceg[1],radii.ceg[2],p100Scale(0),p100Scale(this.cegma_partial),'asm-ceg_part');
     var ccag = ccg.append('g')
       .attr("id","asm-cegma_completeness_axis");
     ccag.append('circle').attr('r',radii.ceg[1]).attr('class','asm-ceg_line');
     ccag.append('line').attr('y2',-radii.ceg[2]).attr('class','asm-axis');
  	 cegma_axis(ccag,radii,p100Scale);
  }

  var line = d3.svg.line()
    .x(function(d,i) { return Math.cos(pScale(i/10)-Math.PI/2)*(radii.core[1] - cScale(d)); })
    .y(function(d,i) { return Math.sin(pScale(i/10)-Math.PI/2)*(radii.core[1] - cScale(d)); });

  	//plot contig count data if available
  	if(this.contigs){
  	  var ctcg = g.append('g')
          .attr("id","asm-g-contig_count");
      var ctcdg = ctcg.append('g')
          .attr("id","asm-g-contig_count_data");
      var ctg_counts = $.map(this.nctg_count, function(value, index) {
    	  return [value];
		});
	  ctcdg.append("path")
  	      .datum(ctg_counts)
          .attr("class", "asm-contig_count asm-remote")
          .attr("d", line);
    }
  	//plot scaffold count data
  	var scg = g.append('g')
        .attr("id","asm-g-scaffold_count");
    var scdg = scg.append('g')
        .attr("id","asm-g-scaffold_count_data");
    var scaf_counts = $.map(this.npct_count, function(value, index) {
    	 return [value];
	  });
	scdg.append("path")
  	      .datum(scaf_counts)
          .attr("class", "asm-count asm-remote")
          .attr("d", line);


  // plot scaffold lengths
  var slg = g.append('g')
      .attr("id","asm-g-scaffold_length");
  var sldg = slg.append('g')
      .attr("id","asm-g-scaffold_length_data");

  var line = d3.svg.line()
    .x(function(d,i) { return Math.cos(pScale(i/10)-Math.PI/2)*(radii.core[1] - lScale(d)); })
    .y(function(d,i) { return Math.sin(pScale(i/10)-Math.PI/2)*(radii.core[1] - lScale(d)); });
  var revline = d3.svg.line()
    .x(function(d,i) { return Math.cos(pScale((1000-i)/10)-Math.PI/2)*(radii.core[1] - lScale(d)); })
    .y(function(d,i) { return Math.sin(pScale((1000-i)/10)-Math.PI/2)*(radii.core[1] - lScale(d)); });
  var scaf_lengths = $.map(this.npct_length, function(value, index) {
      return [value];
    });
  var zeros = Array.apply(0, Array(1000)).map(function (x, y) { return 0; });
  var hollow = line([0])+'L'+line(scaf_lengths).replace(/M[^L]+?/,'')+revline(zeros).replace('M','L')
  sldg.append("path")
  	  .attr("class", "asm-pie")
      .attr("d", hollow)
      .attr("fill-rule","oddeven");


  // plot contig lengths if available
  if (this.contigs){
  var clg = g.append('g')
      .attr("id","asm-g-contig_length");
  var cldg = slg.append('g')
      .attr("id","asm-g-contig_length_data");
  var ctg_lengths = $.map(this.nctg_length, function(value, index) {
      return [value];
    });
  var hollow = line([0])+'L'+line(ctg_lengths).replace(/M[^L]+?/,'')+revline(zeros).replace('M','L')
  sldg.append("path")
  	  .attr("class", "asm-contig")
      .attr("d", hollow)
      .attr("fill-rule","evenodd");
  }
  // highlight n50, n90 and longest scaffold
  var slhg = slg.append('g')
      .attr("id","asm-g-scaffold_length_highlight");
  var long_pct = scaffolds[0] / this.assembly * 100;
  if (long_pct >= 0.1){
    plot_arc(slhg,radii.core[1] - lScale(scaffolds[0]),radii.core[1],0,pScale(long_pct),'asm-longest_pie');
  }
  plot_arc(slhg,radii.core[1] - lScale(npct_length[500]),radii.core[1],0,pScale(50),'asm-n50_pie');
  plot_arc(slhg,radii.core[1] - lScale(npct_length[900]),radii.core[1],0,pScale(90),'asm-n90_pie');
  plot_arc(slhg,radii.core[1] - lScale(npct_length[500]),radii.core[1],pScale(50),pScale(50),'asm-n50_pie asm-highlight');
  if (long_pct >= 0.1){
    plot_arc(slhg,radii.core[1] - lScale(scaffolds[0]),radii.core[1],pScale(long_pct),pScale(long_pct),'asm-longest_pie asm-highlight');
  }

  // add gridlines at powers of 10
  var length_seq = [];
  var power = 2;
  while (Math.pow(10,power) <= longest){
  	length_seq.push(power)
  	power++;
  }
  var slgg = slg.append('g')
      .attr("id","asm-g-scaffold_length_gridlines");
  length_seq.forEach(function(i,index){
    //if(Math.pow(10,i+4) > longest && Math.pow(10,i+1) > npct_length[1000]){
  if(Math.pow(10,i+4) >= longest && Math.pow(10,i+1) > npct_length[900] && Math.pow(10,i) < npct_length[100]){
    //plot_arc(slgg,radii.core[1]-lScale(Math.pow(10,i)),radii.core[1]-lScale(Math.pow(10,i)),pScale(0),pScale(100),'asm-length_axis asm-dashed');
     slgg.append('circle')
  		.attr('r',radii.core[1]-lScale(Math.pow(10,i)))
  		.attr('cx',0)
  		.attr('cy',0)
  		.attr('stroke-dasharray','10,10')
  		.attr('class', 'asm-length_axis');
        }
  	});

  	// plot scaffold/contig count gridlines
  if (!this.contigs){
  var scgg = scg.append('g')
      .attr("id","asm-g-scaffold_count_gridlines");
  [1,2,3,4,5,6,7].forEach(function(i,index){
  	scgg.append('circle')
  	    .attr('class','asm-count_axis')
  	    .attr('r',(radii.core[1] - cScale(Math.pow(10,i))))
  });
  }
  else {
  var ctcgg = scg.append('g')
      .attr("id","asm-g-contig_count_gridlines");
  [1,2,3,4,5,6,7].forEach(function(i,index){
  	ctcgg.append('circle')
  	    .attr('class','asm-count_axis')
  	    .attr('r',(radii.core[1] - cScale(Math.pow(10,i))))
  });
  }


  // plot radial axis
  var mag = g.append('g')
      .attr("id","asm-g-main_axis");
  var slag = g.append('g')
      .attr("id","asm-g-scaffold_length_axis");

  length_seq.forEach(function(i,index){
        if(Math.pow(10,i+3) > longest && Math.pow(10,i+1) > npct_length[1000]){
  slag.append('text')
  		.attr('transform','translate('+(Math.pow(1.5,i)+2)+','+(-radii.core[1]+lScale(Math.pow(10,i))+4)+')')
  		.text(getReadableSeqSizeString(Math.pow(10,i),0))
  		.attr('class', 'asm-length_label');

        }
  		slag.append('line')
  		.attr('x1',0)
  		.attr('y1',-radii.core[1]+lScale(Math.pow(10,i)))
  		.attr('x2',Math.pow(1.5,i))
  		.attr('y2',-radii.core[1]+lScale(Math.pow(10,i)))
        .attr('class', 'asm-majorTick');
  	});
  	slag.append('line')
  	.attr("class","asm-length asm-axis")
  		.attr('x1',0)
  		.attr('y1',-radii.core[1])
  		.attr('x2',0)
  		.attr('y2',0)

	// draw circumferential axis
	circumference_axis(mag,radii,pScale);

   // draw legends
   var lg = g.append('g')
      .attr("id","asm-g-legend");

  // draw CEGMA legend
  if (this.cegma_complete){
  var lccg = lg.append('g')
      .attr("id","asm-g-cegma_completeness_legend");
  	var txt = lccg.append('text')
        .attr('transform', 'translate('+(size/2-230)+','+(-size/2+20)+')')
        .attr('class','asm-tr_title');
  	txt.append('tspan').text('CEGMA completeness');
  	  var key = lccg.append('g').attr('transform', 'translate('+(size/2-230)+','+(-size/2+28)+')');
  	key.append('rect').attr('height',w).attr('width',w).attr('class','asm-ceg_comp asm-toggle');
  	key.append('text').attr('x',w+3).attr('y',w-1).text('Complete ('+this.cegma_complete.toFixed(1)+'%)').attr('class','asm-key');
  	key.append('rect').attr('y',w*1.5).attr('height',w).attr('width',w).attr('class','asm-ceg_part asm-toggle');
  	key.append('text').attr('x',w+3).attr('y',w*2.5-1).text('Partial ('+this.cegma_partial.toFixed(1)+'%)').attr('class','asm-key');
  }

   //draw base composition legend
   var lbcg = lg.append('g')
      .attr("id","asm-g-base_composition_legend");
   var txt = lbcg.append('text')
        .attr('transform', 'translate('+(size/2-140)+','+(size/2-110)+')')
        .attr('class','asm-br_title');
  	txt.append('tspan').text('Assembly');
  	txt.append('tspan').text('base composition').attr('x',0).attr('dy',18);
  	var key = lbcg.append('g').attr('transform', 'translate('+(size/2-140)+','+(size/2-83)+')');
  	key.append('rect').attr('height',w).attr('width',w).attr('class','asm-gc asm-toggle');
  	key.append('text').attr('x',w+2).attr('y',w-1).text('GC ('+this.GC+'%)').attr('class','asm-key');
  	key.append('rect').attr('y',w*1.5).attr('height',w).attr('width',w).attr('class','asm-atgc asm-toggle');
  	key.append('text').attr('x',w+2).attr('y',w*2.5-1).text('AT ('+(this.ATGC-this.GC).toFixed(1)+'%)').attr('class','asm-key');
  	key.append('rect').attr('y',w*3).attr('height',w).attr('width',w).attr('class','asm-ns asm-toggle');
  	key.append('text').attr('x',w+2).attr('y',w*4-1).text('N ('+n.toFixed(1)+'%)').attr('class','asm-key');


   //draw scaffold legend
   var lsg = lg.append('g')
      .attr("id","asm-g-scaffold_legend");
   var txt = lsg.append('text')
        .attr('transform', 'translate('+(-size/2+10)+','+(-size/2+20)+')')
        .attr('class','asm-tl_title');
  	txt.append('tspan').text('Scaffold statistics');
  	//txt.append('tspan').text('distribution').attr('x',0).attr('dy',20);

  	var key = lsg.append('g').attr('transform', 'translate('+(-size/2+10)+','+(-size/2+28)+')');
  	key.append('rect').attr('height',w).attr('width',w).attr('class','asm-count asm-toggle');
  	var count_txt = key.append('text').attr('x',w+3).attr('y',w-1).attr('class','asm-key')
  		count_txt.append('tspan').text('Log')
  		count_txt.append('tspan').attr('baseline-shift','sub').attr('font-size','75%').text(10)
  		count_txt.append('tspan').text(' scaffold count (total '+this.scaffolds.length.toLocaleString()+')');
  	key.append('rect').attr('y',w*1.5).attr('height',w).attr('width',w).attr('class','asm-pie asm-toggle');
  	key.append('text').attr('x',w+3).attr('y',w*2.5-1).text('Scaffold length (total '+getReadableSeqSizeString(this.assembly,0)+')').attr('class','asm-key');


  	key.append('rect').attr('y',w*3).attr('height',w).attr('width',w).attr('class','asm-longest_pie asm-toggle');
  	key.append('text').attr('x',w+3).attr('y',w*4-1).text('Longest scaffold ('+getReadableSeqSizeString(this.scaffolds[0])+')').attr('class','asm-key');
  	key.append('rect').attr('y',w*4.5).attr('height',w).attr('width',w).attr('class','asm-n50_pie asm-toggle');
  	key.append('text').attr('x',w+3).attr('y',w*5.5-1).text('N50 length ('+getReadableSeqSizeString(this.npct_length[500])+')').attr('class','asm-key');
  	key.append('rect').attr('y',w*6).attr('height',w).attr('width',w).attr('class','asm-n90_pie asm-toggle');
  	key.append('text').attr('x',w+3).attr('y',w*7-1).text('N90 length ('+getReadableSeqSizeString(this.npct_length[900])+')').attr('class','asm-key');

    //draw contig legend if available
	if (this.contigs){
   var lctg = lg.append('g')
      .attr("id","asm-g-contig_legend");
   var txt = lctg.append('text')
        .attr('transform', 'translate('+(-size/2+10)+','+(size/2-70)+')')
        .attr('class','asm-bl_title');
  	txt.append('tspan').text('Contig statistics');

  	var key = lctg.append('g').attr('transform', 'translate('+(-size/2+10)+','+(size/2-62)+')');
  	key.append('rect').attr('height',w).attr('width',w).attr('class','asm-contig_count asm-toggle');
  	var count_txt = key.append('text').attr('x',w+2).attr('y',w-1).attr('class','asm-key')
  		count_txt.append('tspan').text('Log')
  		count_txt.append('tspan').attr('baseline-shift','sub').attr('font-size','75%').text(10)
  		count_txt.append('tspan').text(' contig count (total '+this.contigs.length.toLocaleString()+')');
  	key.append('rect').attr('y',w*1.5).attr('height',w).attr('width',w).attr('class','asm-contig asm-toggle');
  	key.append('text').attr('x',w+3).attr('y',w*2.5-1).text('Contig length').attr('class','asm-key');
  	}


  	// toggle plot features
  	$('.asm-toggle').on('click',function(){
  		var button = this;
  		var classNames = $(this).attr("class").toString().split(' ');
  		if ($(button).css('fill') != "rgb(255, 255, 255)" && $(button).css('fill') != "#ffffff"){
  		  $(button).css({fill: "rgb(255, 255, 255)" });
  		  $.each(classNames, function (i, className) {
              if (className != 'asm-toggle'){
        	    $('.'+className).each(function(){
        	      if (this != button){
            	    $(this).css({visibility: "hidden" })
            	  }
            	});
              }
          });
        }
        else {
          var stroke = $(button).css("stroke");
          $(button).css({fill: stroke });
          $.each(classNames, function (i, className) {
              if (className != 'asm-toggle'){
        	    $('.'+className).each(function(){
            	  if (this != button){
            	    $(this).css({visibility: "visible" })
            	  }
            	});
              }
              if (className == 'asm-count'){
        	    $('.asm-contig_count.asm-toggle').css({fill: "rgb(255, 255, 255)" })
            	$('.asm-contig_count.asm-remote').css({visibility: "hidden" })
              }
              if (className == 'asm-contig_count'){
        	    $('.asm-count.asm-toggle').css({fill: "rgb(255, 255, 255)" })
            	$('.asm-count.asm-remote').css({visibility: "hidden" })
              }
          });
        }
  	});


  	// show stats for any N value on mouseover
  	var overlay = g.append('g');
  	var path = overlay.append('path');
  	var overoverlay = overlay.append('g');
  	var output = overlay.append('g').attr('transform', 'translate('+(size/2-142)+','+(size/2-128)+')');
  	var output_rect = output.append('rect').attr('class', 'asm-live_stats hidden').attr('height',110).attr('width',150);
  	var output_text = output.append('g').attr('transform', 'translate('+(2)+','+(18)+')').attr('class', 'hidden');
  	var stat_circle = overoverlay.append('circle').attr('r',radii.core[1]).attr('fill','white').style('opacity',0);
  	stat_circle.on('mousemove', function () {
  	    output_rect.classed('hidden',false);
   		output_text.classed('hidden',false);
   		path.classed('hidden',false);
   		output_text.selectAll('text').remove();


   		var point = d3.mouse(this);
   		var angle =  (50.5 + 50 / Math.PI * Math.atan2(-point[0],  point[1])).toFixed(0);
      angle = Math.floor(angle * p100Scale(100) / pScale(100) + 0.1)

if (angle <= 100){
	var arc = d3.svg.arc()
      	.innerRadius(radii.core[1])
        .outerRadius(radii.core[0])
        .startAngle(pScale(angle-1))
        .endAngle(pScale(angle));
      path
        .attr('d', arc)
        .attr('class', 'asm-live_segment');


        var txt = output_text.append('text')
            .attr('class','asm-live_title');
  		txt.append('tspan').text('N'+angle);
  		output_text.append('text').attr('y',18).text(npct_count[(angle*10)].toLocaleString() + ' scaffolds').attr('class','asm-key');
  		output_text.append('text').attr('x',120).attr('y',w*1.2+18).text('>= ' + getReadableSeqSizeString(npct_length[(angle*10)])).attr('class','asm-key asm-right');
  		if (nctg_length){
  			output_text.append('text').attr('y',w*3+18).text(nctg_count[(angle*10)].toLocaleString() + ' contigs').attr('class','asm-key');
  			output_text.append('text').attr('x',120).attr('y',w*4.2+18).text('>= ' + getReadableSeqSizeString(nctg_length[(angle*10)])).attr('class','asm-key asm-right');
      }
  		}
      else {
        output_rect.classed('hidden',true);
  	    output_text.classed('hidden',true);
  	    path.classed('hidden',true);
      }
  	});
	stat_circle.on('mouseout', function () {
	    output_rect.classed('hidden',true);
	    output_text.classed('hidden',true);
	    path.classed('hidden',true);
    });

}

function circumference_axis (parent,radii,scale){
	var g = parent.append('g');
  var axis = d3.svg.arc()
      	.innerRadius(radii.core[1])
        .outerRadius(radii.core[1])
        .startAngle(scale(0) )
        .endAngle(scale(100));
      g.append('path')
        .attr('d', axis)
        .attr('class', 'asm-axis');
	var seq = Array.apply(0, Array(50)).map(function (x, y) { return y * 2; });
  	seq.forEach(function(i,index){
  		var tick = d3.svg.arc()
      	.innerRadius(radii.core.minorTick[0])
        .outerRadius(radii.core.minorTick[1])
        .startAngle(scale(i))
        .endAngle(scale(i));
		g.append('path')
        .attr('d', tick)
        .attr('class', 'asm-minorTick');
  	});
  	var seq = Array.apply(0, Array(11)).map(function (x, y) { return y * 10; });
  	seq.forEach(function(i,index){
  		var tick = d3.svg.arc()
      	.innerRadius(radii.core.majorTick[0])
        .outerRadius(radii.core.majorTick[1])
        .startAngle(scale(i))
        .endAngle(scale(i));
		g.append('path')
        .attr('d', tick)
        .attr('class', 'asm-majorTick');
        var x = Math.cos(scale(i)-Math.PI/2)*(radii.core.majorTick[1]+10);
        var y = Math.sin(scale(i)-Math.PI/2)*(radii.core.majorTick[1]+10);
        g.append('text')
        .text(function(){return index > 0 ? (index == 10 && scale(100) < 1.96 * Math.PI) ? '100%' : index < 10 ? index*10 : '' : '0%'})
        .attr('transform', 'translate('+x+','+y+') rotate('+scale(i)/(Math.PI/180)+')');
  	});
}


function percent_axis (parent,radii,scale){
	var g = parent.append('g');
	var axis = d3.svg.arc()
      	.innerRadius(radii.percent[0])
        .outerRadius(radii.percent[0])
        .startAngle(scale(0) )
        .endAngle(scale(100));
      g.append('path')
        .attr('d', axis)
        .attr('class', 'asm-axis');
  var seq = Array.apply(0, Array(11)).map(function (x, y) { return y * 10; });
  seq.forEach(function(d,index){
    var arc = d3.svg.arc()
      			.innerRadius(radii.percent.majorTick[0])
        		.outerRadius(radii.percent.majorTick[1])
        		.startAngle(scale(d) )
        		.endAngle(scale(d));
  	    g.append('path')
  		    .attr('d',arc)
            .attr('class', 'asm-majorTick');

    });
	var seq = Array.apply(0, Array(50)).map(function (x, y) { return y*2; });
  seq.forEach(function(d){
    var arc = d3.svg.arc()
      			.innerRadius(radii.percent.minorTick[0])
        		.outerRadius(radii.percent.minorTick[1])
        		.startAngle(scale(d) )
        		.endAngle(scale(d));
  	g.append('path')
  		.attr('d',arc)
        .attr('class', 'asm-minorTick');
    });
}


function cegma_axis (parent,radii,scale){
	var g = parent.append('g');
	var axis = d3.svg.arc()
      	.innerRadius(radii.ceg[2])
        .outerRadius(radii.ceg[2])
        .startAngle(scale(0) )
        .endAngle(scale(100));
      g.append('path')
        .attr('d', axis)
        .attr('class', 'asm-axis');
  var seq = Array.apply(0, Array(10)).map(function (x, y) { return y * 10; });
  seq.forEach(function(d,index){
    var arc = d3.svg.arc()
      			.innerRadius(radii.ceg.majorTick[0])
        		.outerRadius(radii.ceg.majorTick[1])
        		.startAngle(scale(d) )
        		.endAngle(scale(d));
  	    g.append('path')
  		    .attr('d',arc)
            .attr('class', 'asm-majorTick');
        if (d % 20 == 0){
        var x = Math.cos(scale(d)-Math.PI/2)*(radii.ceg.majorTick[1]+5);
        var y = Math.sin(scale(d)-Math.PI/2)*(radii.ceg.majorTick[1]+5);
        g.append('text')
        .text(function(){return d > 0 ? d : d+'%'})
        .attr('transform', 'translate('+x+','+y+') rotate('+scale(d)/(Math.PI/180)+')');
        }
    });

    var seq = Array.apply(0, Array(21)).map(function (x, y) { return y*5; });
  seq.forEach(function(d){
    var arc = d3.svg.arc()
      			.innerRadius(radii.ceg.minorTick[0])
        		.outerRadius(radii.ceg.minorTick[1])
        		.startAngle(scale(d) )
        		.endAngle(scale(d));
  	g.append('path')
  		.attr('d',arc)
        .attr('class', 'asm-minorTick');
    });

}

function plot_arc (parent,inner,outer,start,end,css){
	var arc = d3.svg.arc()
      	.innerRadius(inner)
        .outerRadius(outer)
        .startAngle(start)
        .endAngle(end);
      parent.append('path')
        .attr('d', arc)
        .attr('class', css);
}

function plot_rect (parent,x1,y1,width,height,css){
	parent.append('rect')
	    .attr('x', x1)
	    .attr('y', y1-height+20)
	    .attr('width', width)
	    .attr('height', height)
        .attr('class', css);
}
