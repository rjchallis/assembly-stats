

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
 
function Assembly( stats,scaffolds ) { 
  var sum = scaffolds.reduce(function(previousValue, currentValue, index, array) {
    return previousValue + currentValue;
  });
  var ctgsum = contigs.reduce(function(previousValue, currentValue, index, array) {
    return previousValue + currentValue;
  });
  this.genome = stats.genome; 
  this.assembly = stats.assembly; 
  this.N = stats.N ? stats.N <= 100 ? stats.N < 1 ? stats.N * 100 : stats.N : stats.N / this.assembly * 100 : 0;
  this.ATGC = stats.ATGC ? stats.ATGC <= 100 ? stats.ATGC <= 1 ? stats.ATGC * 100 : stats.ATGC : stats.ATGC / this.assembly * 100 : 100 - this.N;
  this.GC = stats.GC ? stats.GC <= 100 ? stats.GC <= 1 ? stats.GC * 100 : stats.GC : stats.GC / this.assembly * 100 : 0;
  this.cegma_complete = stats.cegma_complete;
  this.cegma_partial = stats.cegma_partial;
  this.scaffolds = scaffolds.sort(function(a, b){return b-a});
  var npct_length = {};
  var npct_count = {};
  var lsum = 0;
  this.scaffolds.forEach(function(length,index,array){
	var new_sum = lsum + length;
	if (Math.floor(new_sum/sum*1000) > Math.floor(lsum/sum*100)){
		npct_length[Math.floor(new_sum/sum*1000)] = length;
		npct_count[Math.floor(new_sum/sum*1000)] = index;
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
  
  this.contigs = contigs.sort(function(a, b){return b-a});
  var nctg_length = {};
  var nctg_count = {};
  var lsum = 0;
  this.contigs.forEach(function(length,index,array){
	var new_sum = lsum + length;
	if (Math.floor(new_sum/ctgsum*1000) > Math.floor(lsum/ctgsum*100)){
		nctg_length[Math.floor(new_sum/ctgsum*1000)] = length;
		nctg_count[Math.floor(new_sum/ctgsum*1000)] = index;
	}
	lsum = new_sum;
  });
  this.seq.forEach(function(i,index){
  	if (!nctg_length[i]) nctg_length[i] = nctg_length[(i+1)];
  	if (!nctg_count[i]) nctg_count[i] = nctg_count[(i+1)];
  	console.log(i+' '+nctg_count[i]);
  });
  this.nctg_length = nctg_length;
  this.nctg_count = nctg_count;
  
  
  this.scale = {};
  this.setScale('percent','linear',[0,100],[180* (Math.PI/180),90* (Math.PI/180)]);
  this.setScale('count','log',[1,1e6],[100,1]); // range will be updated when drawing
  this.setScale('length','sqrt',[1,1e6],[1,100]); // range will be updated range when drawing
}

Assembly.prototype.setScale = function(element,scaling,domain,range){
  this.scale[element] = scaling == 'log' ? d3.scale.log() : scaling == 'sqrt' ? d3.scale.sqrt() : d3.scale.linear();
  this.scale[element].domain(domain);
  this.scale[element].range(range);
}

Assembly.prototype.drawPlot = function(parent){
  
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
  this.scale['length'].domain([1,this.scaffolds[0]])
  this.scale['length'].range([radii.core[0],radii.core[1]])
  this.scale['count'].range([radii.core[1],radii.core[0]+radii.core[1]/3])
  this.scale['percent'].range([0,(2 * Math.PI)])
  
  var lScale = this.scale['length'];
  var cScale = this.scale['count'];
  var pScale = this.scale['percent'];
  var npct_length = this.npct_length;
  var npct_count = this.npct_count;
  var nctg_length = this.nctg_length;
  var nctg_count = this.nctg_count;
  var scaffolds = this.scaffolds;
  var contigs = this.contigs;
  
  // create a group for the plot
  var g = parent.append('g')
      .attr("transform","translate("+size/2+","+size/2+")")
      .attr("id","asm-g-plot");
  
  /* plot expected genome size if available
  // makes plot too complex as interpretation is unclear
  if (this.genome){
   var egg = g.append('g')
       .attr("id","asm-expected_genome");
  	var egdg = egg.append('g')
       .attr("id","asm-expected_genome_data");
    var pct = this.genome/this.assembly*100;
    plot_arc(egdg,radii.genome[0],radii.genome[1],pScale(0),pScale(pct),'asm-genome');
    while (pct > 100){
    	pct -= 100;
    	plot_arc(egdg,radii.genome[0],radii.genome[1],pScale(0),pScale(pct),'asm-short_genome');
    }
  }
  */
      
  // draw base composition axis fill
  var bcg = g.append('g')
      .attr("id","asm-g-base_composition");
  var bcdg = bcg.append('g')
      .attr("id","asm-g-base_composition_data");
  var atgc = this.ATGC;
  var n = 100 - atgc;
  var gc_start = n / 100 * this.GC;
  plot_arc(bcdg,radii.percent[0],radii.percent[1],pScale(0),pScale(100),'asm-ns');
  plot_arc(bcdg,radii.percent[0],radii.percent[1],pScale(gc_start),pScale(gc_start+atgc),'asm-atgc');
  plot_arc(bcdg,radii.percent[0],radii.percent[1],pScale(gc_start),pScale(this.GC),'asm-gc');
  var bcag = bcg.append('g')
      .attr("id","asm-g-base_composition_axis");
  percent_axis(bcag,radii,pScale);
  
  /* plot expected genome size if available
  if (this.genome){
   var egg = g.append('g')
       .attr('transform','translate('+(radii.percent[1]+tick*4)+','+(radii.percent[1])+')')
      .attr("id","asm-expected_genome_size");
  	var egdg = egg.append('g')
       .attr("id","asm-expected_genome_size");
  	 egdg.append('circle')
  	     .attr('r',(Math.sqrt(this.genome/this.assembly)*radii.genome[1]))
  	     .attr('class','asm-genome');
  	 egdg.append('circle')
  	     .attr('r',radii.genome[1])
  	     .attr('class','asm-assembly');
  	 egdg.append('circle')
  	     .attr('r',(Math.sqrt(this.genome/this.assembly)*radii.genome[1]))
  	     .attr('class','asm-genome')
  	     .attr('style','fill:none;');
  	 egdg.append('circle')
  	     .attr('r',radii.genome[1])
  	     .attr('class','asm-axis');
  	 egdg.append('line')
  	     .attr('y1',-radii.genome[1])
  	     .attr('class','asm-axis');
  	 
  }
  */
  // plot CEGMA completeness if available
  if (this.cegma_complete){
   var ccg = g.append('g')
       .attr('transform','translate('+(radii.percent[1]+tick*3)+','+(-radii.percent[1]-tick*2)+')')
      .attr("id","asm-cegma_completeness");
  	var ccdg = ccg.append('g')
       .attr("id","asm-cegma_completeness_data");
  	 plot_arc(ccdg,radii.ceg[0],radii.ceg[1],pScale(0),pScale(this.cegma_complete),'asm-ceg_comp');
     plot_arc(ccdg,radii.ceg[1],radii.ceg[2],pScale(0),pScale(this.cegma_partial),'asm-ceg_part');
     var ccag = ccg.append('g')
       .attr("id","asm-cegma_completeness_axis");
     ccag.append('circle').attr('r',radii.ceg[1]).attr('class','asm-ceg_line');
     ccag.append('line').attr('y2',-radii.ceg[2]).attr('class','asm-axis');
  	 cegma_axis(ccag,radii,pScale);
  }
  
  var power = 6;
  	
  	//plot contig count data if available
  	if(this.contigs){
  	var ctcg = g.append('g')
      .attr("id","asm-g-contig_count");
  var ctcdg = ctcg.append('g')
      .attr("id","asm-g-contig_count_data");
  while (nctg_count[1000] < Math.pow(10,power)){
  	power--;
  }
  this.seq.forEach(function(i,index){
  	if (i <= 1000){
  		plot_arc(ctcdg,radii.core[0],radii.core[1] - cScale(nctg_count[i]),pScale(i/10),pScale(100),'asm-contig_count asm-remote');
  	  }
  });
  	}
  	//plot scaffold count data
  	var scg = g.append('g')
      .attr("id","asm-g-scaffold_count");
  var scdg = scg.append('g')
      .attr("id","asm-g-scaffold_count_data");
  if (!this.contigs){
  while (npct_count[1000] < Math.pow(10,power)){
  	power--;
  }
  }
  this.seq.forEach(function(i,index){
  	if (i <= 1000){
  		plot_arc(scdg,radii.core[0],radii.core[1] - cScale(npct_count[i]),pScale(i/10),pScale(100),'asm-count');
  	  }
  });
  
  
  
  // plot scaffold lengths 
  var slg = g.append('g')
      .attr("id","asm-g-scaffold_length");
  var sldg = slg.append('g')
      .attr("id","asm-g-scaffold_length_data");
  var long_pct = -1; // thousandths of genome covered by longest scaffold
  this.seq.forEach(function(i,index){
  	if (i <= 1000){
  		if (npct_length[i] == scaffolds[0] && npct_length[(i+1)] < scaffolds[0]){
  		  long_pct = i;
  		}
  		else if (npct_length[i] < scaffolds[0]){
  		  plot_arc(sldg,radii.core[1] - lScale(npct_length[i]),radii.core[1],0,pScale(i/10),'asm-pie');
  		}
  	  }
  });
  
  // plot contig lengths if available
  var clg = g.append('g')
      .attr("id","asm-g-contig_length");
  var cldg = slg.append('g')
      .attr("id","asm-g-contig_length_data");
  this.seq.forEach(function(i,index){
  	if (i <= 1000){
  		plot_arc(cldg,radii.core[1] - lScale(nctg_length[i]),radii.core[1],0,pScale(i/10),'asm-contig');
  	  }
  });
  
  // highlight n50, n90 and longest scaffold 
  var slhg = slg.append('g')
      .attr("id","asm-g-scaffold_length_highlight");
  if (long_pct > -1){
    plot_arc(slhg,radii.core[1] - lScale(npct_length[long_pct]),radii.core[1],0,pScale(long_pct/10),'asm-longest_pie');
  }
  plot_arc(slhg,radii.core[1] - lScale(npct_length[500]),radii.core[1],0,pScale(50),'asm-n50_pie');
  plot_arc(slhg,radii.core[1] - lScale(npct_length[900]),radii.core[1],0,pScale(90),'asm-n90_pie');
  plot_arc(slhg,radii.core[1] - lScale(npct_length[500]),radii.core[1],pScale(50),pScale(50),'asm-n50_pie');
  if (long_pct > -1){
    plot_arc(slhg,radii.core[1] - lScale(npct_length[long_pct]),radii.core[1],pScale(long_pct/10),pScale(long_pct/10),'asm-longest_pie');
  }
  
  // add gridlines at powers of 10
  var length_seq = [];
  var power = 2;
  while (Math.pow(10,power) <= this.scaffolds[0]){
  	length_seq.push(power)
  	power++;
  }
  var slgg = slg.append('g')
      .attr("id","asm-g-scaffold_length_gridlines");
  length_seq.forEach(function(i,index){
  if(Math.pow(10,i+4) > this.scaffolds[0] && Math.pow(10,i+1) > npct_length[900] && Math.pow(10,i) < npct_length[100]){
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
  this.seq.forEach(function(i,index){
  	if (i <= 1000){
  		if (npct_count[i] < Math.pow(10,power)){
			plot_arc(scgg,radii.core[1] - cScale(Math.pow(10,power)),radii.core[1] - cScale(Math.pow(10,power)),pScale(i/10),pScale(100),'asm-count_axis');
			power--;
		}
  	  }
  });
  }
  else {
  var ctcgg = scg.append('g')
      .attr("id","asm-g-contig_count_gridlines");
  this.seq.forEach(function(i,index){
  	if (i <= 1000){
  		if (nctg_count[i] < Math.pow(10,power)){
			plot_arc(ctcgg,radii.core[1] - cScale(Math.pow(10,power)),radii.core[1] - cScale(Math.pow(10,power)),pScale(i/10),pScale(100),'asm-count_axis');
			power--;
		}
  	  }
  });
  }
  
  
  // plot radial axis
  var mag = g.append('g')
      .attr("id","asm-g-main_axis");
  var slag = g.append('g')
      .attr("id","asm-g-scaffold_length_axis");
  
  length_seq.forEach(function(i,index){
        if(Math.pow(10,i+3) > this.scaffolds[0] && Math.pow(10,i+1) > npct_length[1000]){
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
	circumference_axis(mag,radii);
  
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
  	key.append('text').attr('x',w+2).attr('y',w*2.5-1).text('AT ('+(atgc-this.GC).toFixed(1)+'%)').attr('class','asm-key');
  	key.append('rect').attr('y',w*3).attr('height',w).attr('width',w).attr('class','asm-ns asm-toggle');
  	key.append('text').attr('x',w+2).attr('y',w*4-1).text('N ('+n.toFixed(1)+'%)').attr('class','asm-key');
  	

   /*draw genome size  legend if available
   if (this.genome){
   var legg = lg.append('g')
      .attr("id","asm-g-genome_legend");
    var key = legg.append('g').attr('transform', 'translate('+(size/2-190)+','+(size/2-62)+')');
  	key.append('rect').attr('height',w).attr('width',w).attr('class','asm-assembly asm-toggle');
  	key.append('text').attr('x',w+3).attr('y',w-1).text('Assembly length ('+getReadableSeqSizeString(this.assembly,0)+')').attr('class','asm-key');
  	key.append('rect').attr('y',w*1.5).attr('height',w).attr('width',w).attr('class','asm-genome asm-toggle');
  	key.append('text').attr('x',w+3).attr('y',w*2.5-1).text('Expected length ('+getReadableSeqSizeString(this.genome,0)+')').attr('class','asm-key');
  	}*/

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
  	/*
  	key.append('rect').attr('y',w*7.5).attr('height',w).attr('width',w).attr('class','asm-gc asm-toggle');
  	key.append('text').attr('x',w+3).attr('y',w*8.5-1).text('GC ('+this.GC+'%)').attr('class','asm-key');
  	key.append('rect').attr('y',w*9).attr('height',w).attr('width',w).attr('class','asm-atgc asm-toggle');
  	key.append('text').attr('x',w+3).attr('y',w*10-1).text('AT ('+(atgc-this.GC).toFixed(1)+'%)').attr('class','asm-key');
  	key.append('rect').attr('y',w*10.5).attr('height',w).attr('width',w).attr('class','asm-ns asm-toggle');
  	key.append('text').attr('x',w+3).attr('y',w*11.5-1).text('N ('+n.toFixed(1)+'%)').attr('class','asm-key');
  	*/

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
  		if ($(button).css('fill') != "rgb(255, 255, 255)"){
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
  	})
  	
}

function circumference_axis (parent,radii){
	var g = parent.append('g');
	g.append('circle')
		.attr('r',radii.core[1])
		.attr('class','asm-axis');
	var seq = Array.apply(0, Array(50)).map(function (x, y) { return y * 7.2 * (Math.PI/180); });
  	seq.forEach(function(i,index){
  		var tick = d3.svg.arc()
      	.innerRadius(radii.core.minorTick[0])
        .outerRadius(radii.core.minorTick[1])
        .startAngle(i)
        .endAngle(i);
		g.append('path')
        .attr('d', tick)
        .attr('class', 'asm-minorTick');
  	});
  	var seq = Array.apply(0, Array(10)).map(function (x, y) { return y * 36 * (Math.PI/180); });
  	seq.forEach(function(i,index){
  		var tick = d3.svg.arc()
      	.innerRadius(radii.core.majorTick[0])
        .outerRadius(radii.core.majorTick[1])
        .startAngle(i)
        .endAngle(i);
		g.append('path')
        .attr('d', tick)
        .attr('class', 'asm-majorTick');
        var x = Math.cos(i-Math.PI/2)*(radii.core.majorTick[1]+10);
        var y = Math.sin(i-Math.PI/2)*(radii.core.majorTick[1]+10);
        g.append('text')
        .text(function(){return index > 0 ? index*10 : '0%'})
        .attr('transform', 'translate('+x+','+y+') rotate('+i/(Math.PI/180)+')');
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


