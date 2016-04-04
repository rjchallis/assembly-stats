#!/usr/bin/perl -w

use strict;
use List::Util qw(reduce);
use JSON;

my $break = 10; # number of consecutive Ns to break scaffolds into contigs
my $bins = 1000; # number of bins to place sequences into

# read scaffolds from file and split into contigs
my $i = -1;
my @scafs;
my @ctgs = ();
while (<>){
  chomp;
  if(/^>/){
    push @ctgs,split_scaf($scafs[$i],$break) if $scafs[$i];
    $i++;
    next;
  }
  $scafs[$i] .= $_;
}
push @ctgs,split_scaf($scafs[$i],$break);

my $output = bin_seqs(\@scafs,$bins,1);
my $extra = bin_seqs(\@ctgs,$bins);

$output->{contigs} = $extra->{scaffolds};
$output->{binned_contig_counts} = $extra->{binned_scaffold_counts};
$output->{contig_count} = $extra->{scaffold_count};
$output->{binned_contig_lengths} = $extra->{binned_scaffold_lengths};

my $json = JSON->new;
$json->pretty(1);
print $json->encode($output),"\n";

sub bin_seqs {
  my ($ref,$bins,$flag) = @_;
  my %return;
  # sort sequences (longest first)
  my @seqs = sort { length $b <=> length $a } @$ref;
  # count sequences
  my $count = scalar @seqs;
  # convert to array of lengths
  my @seq_lengths = map length, @seqs;
  # determine span
  my $span = reduce { $a + $b } @seq_lengths;
  # determine bin size
  my $fbin = $span / $bins;
  my $bin = int ($span / $bins);

  my $sum = 0;
  my $nsum = 0;
  my $gcsum = 0;
  my $catted;
  my $y = 0;
  my $z = 0;
  my @binned_lengths;
  my @binned_counts;
  my @binned_gcs;
  my @binned_ns;
  # place seq lengths and count in $bins bins
  for (my $x = 0; $x < $count; $x++){
    $sum += $seq_lengths[$x];
    $catted .= $seqs[$x];
    while ($sum >= ($y+1)*$fbin){
      $z += $bin;
      $binned_lengths[$y] = $seq_lengths[$x];
      $binned_counts[$y] = $x+1;
      if ($flag){
        # also bin gc and n content
        my $string = substr($catted,0,$bin,'');
        # apply a correction to accommodate non-integer bin-sizes
        my $correction = ($y+1)*$fbin - $z;
        my $extra = 0;
        if ($correction >= 0.5){
          $string .= substr($catted,0,int($correction+0.5),'');
          $extra = int($correction+0.5);
          $z += $correction;
        }
        $binned_ns[$y] = () = $string =~ /n/gi;
        $binned_gcs[$y] = () = $string =~ /[gc]/gi;
        $nsum += $binned_ns[$y];
        $gcsum += $binned_gcs[$y];
        $binned_gcs[$y] /= ($bin+$extra-$binned_ns[$y]) / 100;
        $binned_ns[$y] /= ($bin+$extra) / 100;
      }
      $y++;
    }
  }
  $return{assembly} = $span;
  $return{ATGC} = $span - $nsum;
  $return{GC} = int($gcsum / ($span - $nsum) * 10000) / 100;
  $return{N} = $nsum;
  $return{binned_GCs} = \@binned_gcs;
  $return{binned_Ns} = \@binned_ns;
  splice @seq_lengths, 1;
  $return{scaffolds} = \@seq_lengths;
  $return{scaffold_count} = $count;
  $return{binned_scaffold_counts} = \@binned_counts;
  $return{binned_scaffold_lengths} = \@binned_lengths;
  return \%return
}



sub split_scaf {
  # break scaffolds into contigs on runs of >= $break Ns
  my ($seq,$break) = @_;
  my @ctgs = split /N{$break,}/i,$seq;
  return @ctgs;
}
