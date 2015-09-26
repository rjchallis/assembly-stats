

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
  this.assembly = stats.assembly; 
  this.N = stats.N ? stats.N < 100 ? stats.N < 1 ? stats.N : stats.N / 100 : stats.N / this.assembly : 0;
  this.ATGC = stats.ATGC ? stats.ATGC < 100 ? stats.ATGC < 1 ? stats.ATGC : stats.ATGC / 100 : stats.ATGC / this.assembly : 1 - this.N;
  this.GC = stats.GC;// < 100 ? stats.GC < 1 ? stats.GC : stats.GC / 100 : 50; // TODO: fix last condition 
  this.COMP = stats.CEG_comp;// < 100 ? stats.GC < 1 ? stats.GC : stats.GC / 100 : 50; // TODO: fix last condition 
  this.PART = stats.CEG_part;// < 100 ? stats.GC < 1 ? stats.GC : stats.GC / 100 : 50; // TODO: fix last condition 
  this.scaffolds = scaffolds.sort(function(a, b){return b-a});
  var npct = {};
  var npct_len = {};
  var lsum = 0;
  this.scaffolds.forEach(function(length,index,array){
	var new_sum = lsum + length;
	if (Math.floor(new_sum/sum*1000) > Math.floor(lsum/sum*100)){
		npct[Math.floor(new_sum/sum*1000)] = length;
		npct_len[Math.floor(new_sum/sum*1000)] = index;
	}
	lsum = new_sum;
  });
  this.seq = Array.apply(0, Array(1000)).map(function (x, y) { return 1000 - y; });
  this.seq.forEach(function(i,index){
  	if (!npct[i]) npct[i] = npct[(i+1)];
  	if (!npct_len[i]) npct_len[i] = npct_len[(i+1)];
  });
  this.npct = npct;
  this.npct_len = npct_len;
  this.scale = {};
  this.setScale('percent','linear',[0,100],[180* (Math.PI/180),90* (Math.PI/180)]);
  this.setScale('proportion','log',[1,1e5],[270* (Math.PI/180),360* (Math.PI/180)]);
  this.setScale('count','log',[1,1e6],[100,1]); // TODO: update range when drawing
  this.setScale('length','sqrt',[1,1e6],[1,100]); // TODO: update range when drawing
}

Assembly.prototype.setScale = function(element,scaling,domain,range){
  this.scale[element] = scaling == 'log' ? d3.scale.log() : scaling == 'sqrt' ? d3.scale.sqrt() : d3.scale.linear();
  this.scale[element].domain(domain);
  this.scale[element].range(range);
}

Assembly.prototype.drawPlot = function(parent,size,margin,tick){
  
  size = 600;
  margin = 100;
  tick = 10;
  parent.attr('width', '100%')
  		.attr('height', '100%')
  		.attr('viewBox','0 0 '+size+' '+size)
  		.attr('preserveAspectRatio','xMidYMid meet')
  var radii = {};
  radii.core = [0,(size-margin*2-tick*2)/2];
  radii.core.majorTick = [radii.core[1],radii.core[1]+tick];
  radii.core.minorTick = [radii.core[1],radii.core[1]+tick/2];
  
  radii.proportion = [radii.core[1]+tick*8,radii.core[1]+tick*4];
  radii.proportion.majorTick = [radii.proportion[0],radii.proportion[0]+tick];
  radii.proportion.minorTick = [radii.proportion[0],radii.proportion[0]+tick/2];
  
  radii.percent = [radii.core[1]+tick*4,radii.core[1]];;
  radii.percent.majorTick = [radii.percent[0],radii.percent[0]-tick];
  radii.percent.minorTick = [radii.percent[0],radii.percent[0]-tick/2];
  
  radii.ceg = [0,tick*2,tick*4];
  radii.ceg.majorTick = [radii.ceg[2],radii.ceg[2]+tick/1.5];
  radii.ceg.minorTick = [radii.ceg[2],radii.ceg[2]+tick/3];
  
  this.radii = radii;
  
  this.scale['length'].domain([1,this.scaffolds[0]])
  this.scale['length'].range([radii.core[0],radii.core[1]])
  this.scale['count'].range([radii.core[1],radii.core[0]+radii.core[1]/3])
  
  //square_mod
  this.scale['proportion'].range([radii.core[0],radii.core[1]*2])
  this.scale['percent'].range([radii.core[0],radii.core[1]*2])
  
  //round_mod
  this.scale['percent'].range([0,(2 * Math.PI)])
  
  //this.scale['count'].range([radii.core[0],radii.core[1]])
  var lScale = this.scale['length'];
  var cScale = this.scale['count'];
  var npct = this.npct;
  var npct_len = this.npct_len;
  var scaffolds = this.scaffolds;
  var g = parent.append('g')
      .attr("transform","translate("+size/2+","+size/2+")");
/*  var count_seq = Array.apply(0, Array(6)).map(function (x, y) { return Math.pow(10,y); });
  	count_seq.forEach(function(i,index){
  	  g.append('circle')
		.attr('r',radii.core[1] - cScale(i))
		.attr('class','axis');
    });*/
   
   
   plot_arc(g,radii.percent[0],radii.percent[1],this.scale['percent'](0),this.scale['percent'](100),'asm-ns');
  var atgc = this.ATGC * 100;
  var n = 100 - atgc;
  var gc_start = n / 100 * this.GC;
  plot_arc(g,radii.percent[0],radii.percent[1],this.scale['percent'](gc_start),this.scale['percent'](gc_start+atgc),'asm-atgc');
  plot_arc(g,radii.percent[0],radii.percent[1],this.scale['percent'](gc_start),this.scale['percent'](this.GC),'asm-gc');
  
  if (this.COMP){
   var cg = g.append('g').attr('transform','translate('+(radii.percent[1]+tick*3)+','+(-radii.percent[1]-tick*2)+')');
  	 plot_arc(cg,radii.ceg[0],radii.ceg[1],this.scale['percent'](0),this.scale['percent'](this.COMP),'asm-ceg_comp');
     plot_arc(cg,radii.ceg[1],radii.ceg[2],this.scale['percent'](0),this.scale['percent'](this.PART),'asm-ceg_part');
     cg.append('circle').attr('r',radii.ceg[1]).attr('class','asm-ceg_line');
     cg.append('line').attr('y2',-radii.ceg[2]).attr('class','asm-axis');
  	cegma_axis(cg,radii,this.scale['percent']);
  }
  
  percent_axis(g,radii,this.scale['percent']);
  
   
    var long_pct = -1;
  this.seq.forEach(function(i,index){
  	if (i <= 1000){
  		if (npct[i] == scaffolds[0] && npct[(i+1)] < scaffolds[0]){
  		  long_pct = i;
  		}
  		else if (npct[i] < scaffolds[0]){
  		  plot_arc(g,radii.core[1] - lScale(npct[i]),radii.core[1],0,i * 360 / 1000 * (Math.PI/180),'asm-pie');
  		}
  	  }
  });
  
  plot_arc(g,radii.core[1] - lScale(npct[500]),radii.core[1],0,500 * 360 / 1000 * (Math.PI/180),'asm-n50_pie');
  plot_arc(g,radii.core[1] - lScale(npct[900]),radii.core[1],0,900 * 360 / 1000 * (Math.PI/180),'asm-n90_pie');
  
  if (long_pct > -1){
  plot_arc(g,radii.core[1] - lScale(npct[long_pct]),radii.core[1],0,long_pct * 360 / 1000 * (Math.PI/180),'asm-longest_pie');
  }
  
  
  
  var length_seq = [];
  var power = 2;
  while (Math.pow(10,power) <= this.scaffolds[0]){
  	length_seq.push(power)
  	power++;
  }
  var lg = g.append("g")
		.attr("class","asm-length asm-axis");
  length_seq.forEach(function(i,index){
  if(Math.pow(10,i+4) > this.scaffolds[0] && Math.pow(10,i+1) > npct[900] && Math.pow(10,i) < npct[100]){
     lg.append('circle')
  		.attr('r',radii.core[1]-lScale(Math.pow(10,i)))
  		.attr('cx',0)
  		.attr('cy',0)
  		.attr('stroke-dasharray','10,10')
  		.attr('class', 'asm-length_axis');
  		
        }
  
  	});
  	
  	var power = 6;
  while (npct_len[1000] < Math.pow(10,power)){
  	power--;
  }
  this.seq.forEach(function(i,index){
  	if (i <= 1000){
  		//plot_arc(g,radii.core[1] - cScale(npct_len[i]),radii.core[1] - cScale(npct_len[i]),i * 360 / 1000 * (Math.PI/180),(i+1) * 360 / 1000 * (Math.PI/180),'count');
		plot_arc(g,radii.core[0],radii.core[1] - cScale(npct_len[i]),i * 360 / 1000 * (Math.PI/180),360 * (Math.PI/180),'asm-count');
  	  }
  });
  this.seq.forEach(function(i,index){
  	if (i <= 1000){
  		//plot_arc(g,radii.core[1] - cScale(npct_len[i]),radii.core[1] - cScale(npct_len[i]),i * 360 / 1000 * (Math.PI/180),(i+1) * 360 / 1000 * (Math.PI/180),'count');
		if (npct_len[i] < Math.pow(10,power)){
			plot_arc(g,radii.core[1] - cScale(Math.pow(10,power)),radii.core[1] - cScale(Math.pow(10,power)),i * 360 / 1000 * (Math.PI/180),360 * (Math.PI/180),'asm-count_axis');
			
			g.append('text')
        		.text(Math.pow(10,power))
        		.attr('transform', 'translate(-4,'+(-radii.core[1] + cScale(Math.pow(10,power))+10)+')')
        		.attr('class','asm-count_label');
			power--;
		}
  	  }
  });
  
  length_seq.forEach(function(i,index){
  
        if(Math.pow(10,i+3) > this.scaffolds[0] && Math.pow(10,i+1) > npct[1000]){
  lg.append('text')
  		.attr('transform','translate('+(Math.pow(1.5,i)+2)+','+(-radii.core[1]+lScale(Math.pow(10,i))+4)+')')
  		.text(getReadableSeqSizeString(Math.pow(10,i),0))
  		.attr('class', 'asm-length_label');
  		
        }
  		lg.append('line')
  		.attr('x1',0)
  		.attr('y1',-radii.core[1]+lScale(Math.pow(10,i)))
  		.attr('x2',Math.pow(1.5,i))
  		.attr('y2',-radii.core[1]+lScale(Math.pow(10,i)))
        .attr('class', 'asm-majorTick');
  	});
  	g.append('line')
  	.attr("class","asm-length asm-axis")
  		.attr('x1',0)
  		.attr('y1',-radii.core[1])
  		.attr('x2',0)
  		.attr('y2',0)
		//.attr("transform","translate(25,525)")
	//.call(length_axis);
  
  	main_axis(g,radii);
  
  	/*var x = -Math.pow(1.5,length_seq[length_seq.length-1])
  	var y = -radii.core[1]+lScale(Math.pow(10,length_seq[length_seq.length-1]))
  	g.append('text')
        .text(Math.pow(10,length_seq[length_seq.length-1]))
        .attr('transform', 'translate('+x+','+y+')')
        .attr('class','asm-length_label');
    g.append('text')
        .text(this.scaffolds[0])
        .attr('transform', 'translate('+10+','+-radii.core[1]/2+') rotate(90)')
        .attr('class','asm-length_label');*/
  	
  //square_mod
  //round_mod
  //plot_arc(g,radii.proportion[0],radii.proportion[1],this.scale['proportion'](1),this.scale['proportion'](this.assembly/this.scaffolds[0]),'asm-genome');
  //plot_rect(g,-radii.proportion[0],radii.core[1],Math.abs(radii.proportion[1]-radii.proportion[0]),this.scale['proportion'](this.assembly/this.scaffolds[0]),'asm-genome');
  //proportion_axis(g,radii,this.scale['proportion']);
   var w = 12;
  	/*
  proportion_axis(g,radii,this.scale['proportion']);
  	var x = Math.cos(this.scale['proportion'](25)-90)*(radii.proportion[1]+20);
    var y = Math.sin(this.scale['proportion'](25)-90)*(radii.proportion[1]+20);
	g.append('text')
        .text(this.assembly)
        .attr('transform', 'translate('+x+','+-y+') rotate('+(276)+')')
        .attr('class','asm-assembly_label');
  	*/
  	var txt = g.append('text')
        .attr('transform', 'translate('+(size/2-230)+','+(-size/2+20)+')')
        .attr('class','asm-tr_title');
  	txt.append('tspan').text('CEGMA completeness');
  	//txt.append('tspan').text('size').attr('x',0).attr('dy',18);
  	//txt.append('tspan').text('size').attr('x',0).attr('dy',18);
      var key = g.append('g').attr('transform', 'translate('+(size/2-230)+','+(-size/2+28)+')');
  	key.append('rect').attr('height',w).attr('width',w).attr('class','asm-ceg_comp');
  	key.append('text').attr('x',w+2).attr('y',w-1).text('Complete ('+this.COMP.toFixed(1)+'%)').attr('class','asm-key');
  	key.append('rect').attr('y',w*1.5).attr('height',w).attr('width',w).attr('class','asm-ceg_part');
  	key.append('text').attr('x',w+2).attr('y',w*2.5-1).text('Partial ('+this.PART.toFixed(1)+'%)').attr('class','asm-key');
  	 

  //square_mod
  //round_mod
  //plot_arc(g,radii.percent[0],radii.percent[1],this.scale['percent'](0),this.scale['percent'](100),'asm-ns');
  //plot_arc(g,radii.percent[0],radii.percent[1],this.scale['percent']((1-this.ATGC)/2*100),this.scale['percent'](100*this.ATGC + (1-this.ATGC)/2*100),'asm-atgc');
  //plot_arc(g,radii.percent[0],radii.percent[1],this.scale['percent']((1-this.ATGC)/2*100),this.scale['percent'](this.GC),'asm-gc');
  //plot_rect(g,radii.percent[1],radii.core[1],Math.abs(radii.proportion[1]-radii.proportion[0]),this.scale['percent'](100),'asm-ns');
  //plot_rect(g,radii.percent[1],radii.core[1],Math.abs(radii.proportion[1]-radii.proportion[0]),this.scale['percent'](100*this.ATGC),'asm-atgc');
  //plot_rect(g,radii.percent[1],radii.core[1],Math.abs(radii.proportion[1]-radii.proportion[0]),this.scale['percent'](this.GC),'asm-gc');
  var txt = g.append('text')
        .attr('transform', 'translate('+(size/2-140)+','+(size/2-110)+')')
        .attr('class','asm-br_title');
  	txt.append('tspan').text('Assembly');
  	txt.append('tspan').text('base composition').attr('x',0).attr('dy',18);
  	//txt.append('tspan').text('composition').attr('x',0).attr('dy',18);
  	
  	var key = g.append('g').attr('transform', 'translate('+(size/2-140)+','+(size/2-83)+')');
  	key.append('rect').attr('height',w).attr('width',w).attr('class','asm-gc');
  	key.append('text').attr('x',w+2).attr('y',w-1).text('GC ('+this.GC+'%)').attr('class','asm-key');
  	key.append('rect').attr('y',w*1.5).attr('height',w).attr('width',w).attr('class','asm-atgc');
  	key.append('text').attr('x',w+2).attr('y',w*2.5-1).text('AT ('+(atgc-this.GC).toFixed(1)+'%)').attr('class','asm-key');
  	key.append('rect').attr('y',w*3).attr('height',w).attr('width',w).attr('class','asm-ns');
  	key.append('text').attr('x',w+2).attr('y',w*4-1).text('N ('+n.toFixed(1)+'%)').attr('class','asm-key');
  	

  var txt = g.append('text')
        .attr('transform', 'translate('+(-size/2+10)+','+(-size/2+20)+')')
        .attr('class','asm-tl_title');
  	txt.append('tspan').text('Scaffold length');
  	txt.append('tspan').text('distribution').attr('x',0).attr('dy',20);
  	//txt.append('tspan').text('distribution').attr('x',0).attr('dy',20);
  	
  	var key = g.append('g').attr('transform', 'translate('+(-size/2+10)+','+(-size/2+50)+')');
  	key.append('rect').attr('height',w).attr('width',w).attr('class','asm-pie');
  	//key.append('text').attr('x',w+2).attr('y',w-1).text('(Scaffold length)').attr('class','asm-key').append('tspan').attr('baseline-shift','super').attr('font-size','75%').text(0.5);
  	key.append('text').attr('x',w+2).attr('y',w-1).text('Scaffold length (total '+getReadableSeqSizeString(this.assembly,0)+')').attr('class','asm-key');
  	key.append('rect').attr('y',w*1.5).attr('height',w).attr('width',w).attr('class','asm-pie');
  	key.append('rect').attr('y',w*1.5).attr('height',w).attr('width',w).attr('class','asm-longest_pie');
  	key.append('text').attr('x',w+2).attr('y',w*2.5-1).text('Longest scaffold ('+getReadableSeqSizeString(this.scaffolds[0])+')').attr('class','asm-key');
  	key.append('rect').attr('y',w*3).attr('height',w).attr('width',w).attr('class','asm-pie');
  	key.append('rect').attr('y',w*3).attr('height',w).attr('width',w).attr('class','asm-n50_pie');
  	key.append('text').attr('x',w+2).attr('y',w*4-1).text('N50 length ('+getReadableSeqSizeString(this.npct[500])+')').attr('class','asm-key');
  	key.append('rect').attr('y',w*4.5).attr('height',w).attr('width',w).attr('class','asm-pie');
  	key.append('rect').attr('y',w*4.5).attr('height',w).attr('width',w).attr('class','asm-n90_pie');
  	key.append('text').attr('x',w+2).attr('y',w*5.5-1).text('N90 length ('+getReadableSeqSizeString(this.npct[900])+')').attr('class','asm-key');
  	

    var txt = g.append('text')
        .attr('transform', 'translate('+(-size/2+10)+','+(size/2-70)+')')
        .attr('class','asm-bl_title');
  	txt.append('tspan').text('Cumulative');
  	txt.append('tspan').text('scaffold number').attr('x',0).attr('dy',20);
  	//txt.append('tspan').text('distribution').attr('x',0).attr('dy',20);
  	
  	var key = g.append('g').attr('transform', 'translate('+(-size/2+10)+','+(size/2-43)+')');
  	key.append('rect').attr('height',w).attr('width',w).attr('class','asm-count');
  	var count_txt = key.append('text').attr('x',w+2).attr('y',w-1).attr('class','asm-key')
  		count_txt.append('tspan').text('Log')
  		count_txt.append('tspan').attr('baseline-shift','sub').attr('font-size','75%').text(10)
  		count_txt.append('tspan').text(' scaffold count (total '+this.scaffolds.length.toLocaleString()+')');
  	
  	
}

function main_axis (parent,radii){
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

function proportion_axis (parent,radii,scale){
  var g = parent.append('g');
  g.attr('transform','translate(0,20)')
  var line = g.append('line');
  	line.attr('x1',-radii.proportion[0])
  	    .attr('y1',radii.core[1])
  	    .attr('x2',-radii.proportion[0])
        .attr('y2',radii.core[1]-scale(100000))
        .attr('class', 'asm-axis');
   
   var seq = Array.apply(0, Array(6)).map(function (x, y) { return Math.pow(10,y); });
  seq.forEach(function(d,index){
    var line = g.append('line');
  	line.attr('x1',-radii.proportion.majorTick[0])
  	    .attr('y1',radii.core[1]-scale(d))
  	    .attr('x2',-radii.proportion.majorTick[1])
        .attr('y2',radii.core[1]-scale(d))
        .attr('class', 'asm-axis');
        g.append('text')
          .text(function(){return index < 6 ? Math.pow(10,index) : 10 + '^' +index})
          .attr('transform', 'translate('+(-radii.proportion.majorTick[1]-10)+','+(radii.core[1]-scale(d))+') rotate(270)');
    });

  var minor = [];
    
  seq.forEach(function(d,index){
  	var tmp = Array.apply(0, Array(9)).map(function (x, y) { return d*(y+1) });
  	if (index < seq.length - 1)	minor = minor.concat(tmp);
  });
  minor.forEach(function(d){
   var line = g.append('line');
  	line.attr('x1',-radii.proportion.minorTick[0])
  	    .attr('y1',radii.core[1]-scale(d))
  	    .attr('x2',-radii.proportion.minorTick[1])
        .attr('y2',radii.core[1]-scale(d))
        .attr('class', 'asm-axis');
    });
   //square_mod
   return;
   
   
   
	var g = parent.append('g');
	var axis = d3.svg.arc()
      	.innerRadius(radii.proportion[0])
        .outerRadius(radii.proportion[0])
        .startAngle(scale(1) )
        .endAngle(scale(100000));
      g.append('path')
        .attr('d', axis)
        .attr('class', 'asm-axis');
  var seq = Array.apply(0, Array(6)).map(function (x, y) { return Math.pow(10,y); });
  seq.forEach(function(d,index){
    var arc = d3.svg.arc()
      			.innerRadius(radii.proportion.majorTick[0])
        		.outerRadius(radii.proportion.majorTick[1])
        		.startAngle(scale(d) )
        		.endAngle(scale(d));
  	    g.append('path')
  	  	  .attr('d',arc)
          .attr('class', 'asm-majorTick');
        var x = Math.cos(scale(d)-Math.PI/2)*(radii.proportion.majorTick[1]+10);
        var y = Math.sin(scale(d)-Math.PI/2)*(radii.proportion.majorTick[1]+10);
        g.append('text')
          .text(function(){return index < 6 ? Math.pow(10,index) : 10 + '^' +index})
          .attr('transform', 'translate('+x+','+y+')  rotate('+scale(d)/(Math.PI/180)+')');
    });


	var minor = [];
    
  seq.forEach(function(d,index){
  	var tmp = Array.apply(0, Array(9)).map(function (x, y) { return d*(y+1) });
  	if (index < seq.length - 1)	minor = minor.concat(tmp);
  });
  
  minor.forEach(function(d){
    var arc = d3.svg.arc()
      			.innerRadius(radii.proportion.minorTick[0])
        		.outerRadius(radii.proportion.minorTick[1])
        		.startAngle(scale(d) )
        		.endAngle(scale(d));
  	g.append('path')
  		.attr('d',arc)
        .attr('class', 'asm-minorTick');
    });
    
  
	
}


function percent_axis (parent,radii,scale){

  //round_mod
  /*
  var g = parent.append('g');
  g.attr('transform','translate(0,20)')
  var line = g.append('line');
  	line.attr('x1',radii.percent[0])
  	    .attr('y1',radii.core[1])
  	    .attr('x2',radii.percent[0])
        .attr('y2',radii.core[1]-scale(100))
        .attr('class', 'asm-axis');
   
   var seq = Array.apply(0, Array(11)).map(function (x, y) { return y*10; });
  seq.forEach(function(d,index){
    var line = g.append('line');
  	line.attr('x1',radii.percent.majorTick[0])
  	    .attr('y1',radii.core[1]-scale(d))
  	    .attr('x2',radii.percent.majorTick[1])
        .attr('y2',radii.core[1]-scale(d))
        .attr('class', 'asm-axis');
        g.append('text')
          .text(function(){return d > 0 && d < 100 ? d : d+'%'})
          .attr('transform', 'translate('+(radii.proportion.majorTick[1]+10)+','+(radii.core[1]-scale(d))+') rotate(90)');
    });

  var seq = Array.apply(0, Array(50)).map(function (x, y) { return y*2; });
  seq.forEach(function(d){
   var line = g.append('line');
  	line.attr('x1',radii.proportion.minorTick[0])
  	    .attr('y1',radii.core[1]-scale(d))
  	    .attr('x2',radii.proportion.minorTick[1])
        .attr('y2',radii.core[1]-scale(d))
        .attr('class', 'asm-axis');
    });
   //square_mod
   return;
  */

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
    //round_mod
    /*
    var seq = Array.apply(0, Array(11)).map(function (x, y) { return y * 10; });
  seq.forEach(function(d,index){
  
    var x = Math.cos(scale(d)-Math.PI/2)*(radii.percent.majorTick[1]+10);
    var y = Math.sin(scale(d)-Math.PI/2)*(radii.percent.majorTick[1]+10);
        g.append('text')
          .text(function(){return d > 0 && d < 100 ? d : d+'%'})
          .attr('transform', 'translate('+x+','+y+') rotate('+(180+scale(d)/(Math.PI/180))+')');
	})
    */
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


