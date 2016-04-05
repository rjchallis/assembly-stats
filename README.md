# assembly_stats
An assembly metric visualisation to allow rapid visual assessment and comparison of assembly quality.  [![DOI](https://zenodo.org/badge/20772/rjchallis/assembly_stats.svg)](https://zenodo.org/badge/latestdoi/20772/rjchallis/assembly_stats)

![Screenshot](http://rjchallis.github.io/assembly_stats/screenshots/assembly_stats.png "Screenshot")

A _de novo_ genome assembly can be summarised by a number of metrics, including:
- Overall assembly length
- Number of scaffolds/contigs
- Length of longest scaffold/contig
- Scaffold/contig N50 and N90
- Assembly base composition, in particular percentage GC and percentage Ns
- CEGMA completeness
- Scaffold/contig length/count distribution

Simply listing or tabulating these values allows comparison, however the fact that such values can vary by several orders of magnitude can present a barrier to easy interpretation, particularly when the ratio of values is being considered.  Tabulated metrics also necessarily omit distributions.  

Simple plots of cumulative scaffold length can be effective and for comparison of two or more assemblies, plotting both on a single set of axes can reveal differences in quality very clearly.  While metrics such as assembly size, number of scaffolds, N50 and N90 can be derived from these plots, other metrics must still be tabulated and the length of the longest scaffold can be particularly difficult to determine.  The scale for the axes is usually chosen to accommodate the data for a single assembly or set of assemblies, meaning that it is usually necessary to replot the data or consider the relative axis scales carefully to compare assemblies that have been plotted separately.

The goal of assembly_stats is to overcome some of these shortcomings to produce a visualisation that allows rapid assessment of most common assembly metrics. The graphic is essentially scale independent so assemblies of any size with different strengths and weaknesses produce distinct patterns that can be recognised at a glance, allowing comparison of assemblies produced at different times either by looking at the plots side by side or simply remembering the keys features from the plot of a previous assembly when considering a new one.

## examples

### using latest input format
[Danaus plexippus](http://rjchallis.github.io/assembly_stats/Dplex.html)

### previous examples
[Lerema accius v1.1](http://rjchallis.github.io/assembly_stats/index.html)

[Heliconius melpomene assembly comparison](http://rjchallis.github.io/assembly_stats/Hmel1_vs_Hmel2.html)

[Rescaled Heliconius melpomene assembly comparison](http://rjchallis.github.io/assembly_stats/Hmel_rescaled.html)


## plot descritption
- click on any colour tile in the legend to toggle visibility of that feature on/off
- The inner radius of the circular plot represents the length of the longest scaffold in the assembly
- The angle subtended by the first (red) segment within this plot indicates the percentage of the assembly that is in the longest scaffold
- The radial axis originates at the circumference and indicates scaffold length
- Subsequent (grey) segments are plotted from the circumference and the length of segment at a given percentage indicates the cumulative percentage of the assembly that is contained within scaffolds of at least that length
- The N50 and N90 scaffold lengths are indicated respectively by dark and light orange arcs that connect to the radial axis for ease of comparison
- The cumulative number of scaffolds within a given percentge of the genome is plotted in purple originating at the centre of the plot
- White scale lines are drawn at successive orders of magnitude from 10 scaffolds onwards
- The fill colour of the circumferential axis indicates the percentage base composition of the assembly: AT = light blue; GC = dark blue; N = grey
- Contig length (if available) is indicated by darker grey segments overlaying the scaffold length plot
- Contig count (if available) may be toggled on to be shown in place of the scaffold count plot
- Complete, fragmented and duplicated BUSCO genes (if available) are shown in mid, light and dark green, respectively in the smaller plot in the upper right corner
- Partial and complete CEGMA values (if available) are shown in light and dark green, respectively in the smaller plot in the upper right corner

## basic usage

### input format

Data to be plotted must be supplied as a JSON format object.  As of version 1.1 data may be pre-binned to improve performance with assemblies containing potentially millions of contigs.  The simplest way to generate this is using the ``asm2stats.pl`` perl script in the ``pl`` folder:

```bash
perl asm2stats.pl genome_assembly.fa > output.json
```

This input format should be preferred as it improves performance and corrects for a bug in the javascript binning code by adjusting bin size to accommodate assembly spans that are not divisible by 1000, however the previous input format (with a full list of scaffold lengths is still supported).

### usage

The simplest plot requires a target div, an assembly span, a count of ACGT bases, the GC percentage and an array of scaffold lengths, however it is best to use the ``asm2stats.pl`` perl script described above to generate a richer, pre-processed input format.  See the ``Dplex.html`` file for a complete example using pre-binned data, basic usage is detailed below:

```html
<div id="assembly_stats">
<script>
  var stats = paste_json_object_from_asm2stats;
  var asm = new Assembly (stats);
  asm.drawPlot('assembly_stats');
</script>
```

The json object contains the following keys:
- ``assembly`` - the total assembly span
- ``ATGC`` - the assembly span without Ns (redundant if ``N`` is specified)
- ``GC`` - the GC percentage of the assembly
- ``N`` - the total number of Ns (redundant if ``ATGC`` is specified)
- ``scaffold_count`` - the total number of scaffolds in the assembly
- ``scaffolds`` - an array of scaffold lengths (only the longest scaffold is needed if ``binned_scaffold_lengths`` and ``binned_scaffold_counts`` are specified)
- ``binned_scaffold_lengths`` - an array of 1000 scaffold lengths representing the N0.1 to N100 scaffold lengths for the assembly
- ``binned_scaffold_counts`` - an array of 1000 scaffold counts representing the N0.1 to N100 scaffold numbers for the assembly
- ``contig_count`` - (optional) the total number of contigs in the assembly
- ``contigs`` - (optional) an array of contig lengths (only the longest contig is needed if ``binned_contig_lengths`` and ``binned_contig_counts`` are specified)
- ``binned_contig_lengths`` - (optional) an array of 1000 contig lengths representing the N0.1 to N100 contig lengths for the assembly
- ``binned_contig_counts`` - (optional) an array of 1000 contig counts representing the N0.1 to N100 contig numbers for the assembly
- ``binned_Ns`` - (optional) an array of 1000 values representing the N content of each bin based on size-sorted scaffold sequences
- ``binned_GCs`` - (optional) an array of 1000 values representing the GC content of each bin based on size-sorted scaffold sequences


Additional data will be plotted, if added to the stats object including:
- CEGMA scores

  ```json
  cegma_complete:    83.87,
  cegma_partial:    95.16
  ```

- BUSCO complete, duplicated, fragmented, missing and number of genes (will be plotted in place of CEGMA if both are present)

  ```json
  busco: { C:87.1,
           D:3.6,
           F:10.1,
           M:2.8,
           n:2675 }
  ```

While the plots were conceived as scale independent visualisations, there are occasions when it is useful to compare assemblies on the same radial (longest scaffold) or circumferential (assembly span) scales.  These scales may be modified on the plot by clicking the grey boxes under the scale heading.  Plots can also be drawn with an specific scale by supplying additional arguments to ``drawPlot()``.

For example to scale the radius to 10 Mb and the circumference to 400 Mb (values smaller than the default will be ignored):

```javascript
  asm.drawPlot('assembly_stats',10000000,400000000);
```

It is also possible to programmatically toggle the visibility of plot features by passing an array of classnames to ``toggleVisible()``:

```javascript
  asm.toggleVisible(['asm-longest_pie','asm-count']);
```
