Assembly.prototype.tabulate = function() {
  var rowdata = [ {'title':'','val':this.name.replace(/_/g,' ')},
                  {'title':'span (bp)','val':this.assembly},
                  {'title':'N (%)','val':this.N.toFixed(2)},
                  {'title':'GC (%)','val':this.GC.toFixed(2)},
                  {'title':'AT (%)','val':(100-this.GC).toFixed(2)},
                  {'title':'scaffold count','val':this.scaffold_count},
                  {'title':'longest scaffold (bp)','val':this.scaffolds[0]},
                  {'title':'scaffold N50 length (bp)','val':this.npct_length[499]},
                  {'title':'scaffold N50 count','val':this.npct_count[499]},
                  {'title':'scaffold N90 length (bp)','val':this.npct_length[899]},
                  {'title':'scaffold N90 count','val':this.npct_count[899]},
                  {'title':'contig count','val':this.contig_count},
                  {'title':'contig N50 length (bp)','val':this.nctg_length[499]},
                  {'title':'contig N50 count','val':this.nctg_count[499]},
                  {'title':'contig N90 length (bp)','val':this.nctg_length[899]},
                  {'title':'contig N90 count','val':this.nctg_count[899]}
                ]
  this.rowdata = rowdata;
  return this;
}

Assembly.prototype.drawTable = function(parent_div) {
  this.tabulate();
  this.parent_div = parent_div;
  var parent = d3.select('#' + parent_div);
  var table = parent.append('table');
  var rows = table.selectAll('tr').data(this.rowdata);
  rows.enter().append('tr');
  rows.append('td').text(function(d){return d.title});
  rows.append('td').text(function(d){return d.val.toLocaleString()});
  this.table = table;
  return this;
}

Assembly.prototype.addColumn = function(altAssembly) {
  altAssembly.tabulate();
  var rows = this.table.selectAll('tr').data(altAssembly.rowdata);
  rows.append('td').text(function(d){return d.val.toLocaleString()});
  return this;
}
