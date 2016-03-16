# assembly_stats
An assembly metric visualisation to allow rapid visual assessment and comparison of assembly quality.  

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
- Partial and complete CEGMA values (if available) are shown in light and dark green, respectively in the smaller plot in the upper right corner

## basic usage
The simplest plot requires a target div, an assembly span, a count of ACGT bases, the GC percentage and an array of scaffold lengths:

```html
<div id="assembly_stats">
<script>
  var stats = { assembly:  298173436,
                ATGC:      289601338,
                GC:        34.4,
                scaffolds: [3082282,2858651,2756976,2480955,2126098,...]
              };
  var asm = new Assembly (stats);
  asm.drawPlot('assembly_stats');
</script>
```

Additional data will be plotted, if added to the stats object including:
- Contig lengths (as for scaffold lengths this must be an array and may be sorted or unsorted):

  ```json
  contigs: [141759,134655,133634,123664,121418,...]
  ```

- CEGMA (complete and partial) scores

  ```json
  cegma_complete:    83.87,
  cegma_partial:    95.16
  ```

- Binned GC and N content (must be presented as an array of 1000 bins for the size-sorted, concatenated scaffold sequences):

  ```json
  GCs: [36.4,35.3,35.5,...],
  Ns:  [4.9,3.5,3,2.5,2,...]
  ```

- Total N count (trivially N count may be given in place of or as well as AGCT count):
- N: 8572098

While the plots were concieved as scale independent visualisations, there are occassions when it is useful to compare assemblies on the same radial (longest scaffold) or circumferential (assembly span) scales.  These scales may be modified on the plot by clicking the grey boxes under the scale heading.  Plots can also be drawn with an specific scale by supplying additional arguments to ``drawPlot()``.

For example to scale the radius to 10 Mb and the circumference to 400 Mb (values smaller than the default will be ignored):

```javascript
  asm.drawPlot('assembly_stats',10000000,400000000);
```
