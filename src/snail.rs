use std::cmp::Ordering;

use serde;
use serde::{Deserialize, Serialize};

use crate::assembly::AssemblyStats;
use crate::cli::Config;

mod compact_float {
    //! rounds a float to 3 decimal places, when serialized into a str, such as for JSON
    //! offsers space savings when such such precision is not needed.
    use serde::{Deserialize, Deserializer, Serializer};

    pub fn serialize<S>(float: &f32, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let s = format!("{:.3}", float);
        let parsed = s.parse::<f32>().unwrap();
        serializer.serialize_f32(parsed)
    }

    pub fn deserialize<'de, D>(deserializer: D) -> Result<f32, D::Error>
    where
        D: Deserializer<'de>,
    {
        f32::deserialize(deserializer)
    }
}

#[derive(Serialize, Deserialize)]
pub struct SummaryStats {
    #[serde(with = "compact_float")]
    min: f32,
    #[serde(with = "compact_float")]
    max: f32,
    #[serde(with = "compact_float")]
    mean: f32,
}

impl SummaryStats {
    pub fn min(&self) -> f32 {
        self.min
    }
    pub fn max(&self) -> f32 {
        self.max
    }
    pub fn mean(&self) -> f32 {
        self.mean
    }
}

#[derive(Serialize, Deserialize)]
pub struct SnailStats {
    #[serde(rename = "assembly")]
    span: u64,
    #[serde(rename = "ATGC")]
    atgc: u64,
    #[serde(skip)]
    gc: u64,
    #[serde(rename = "GC", with = "compact_float")]
    gc_proportion: f32,
    #[serde(rename = "N")]
    n: u64,
    #[serde(rename = "binned_GCs")]
    binned_gcs: Vec<SummaryStats>,
    #[serde(rename = "binned_Ns")]
    binned_ns: Vec<SummaryStats>,
    scaffolds: Vec<u64>,
    scaffold_count: u32,
    binned_scaffold_lengths: Vec<u64>,
    binned_scaffold_counts: Vec<u32>,
    contigs: Vec<u64>,
    contig_count: u32,
    binned_contig_lengths: Vec<u64>,
    binned_contig_counts: Vec<u32>,
}

impl SnailStats {
    pub fn span(&self) -> u64 {
        self.span
    }
    pub fn atgc(&self) -> u64 {
        self.atgc
    }
    pub fn gc(&self) -> u64 {
        self.gc
    }
    pub fn n(&self) -> u64 {
        self.n
    }
    pub fn binned_gcs(&self) -> &Vec<SummaryStats> {
        &self.binned_gcs
    }
    pub fn binned_ns(&self) -> &Vec<SummaryStats> {
        &self.binned_ns
    }
    pub fn scaffolds(&self) -> &Vec<u64> {
        &self.scaffolds
    }
    pub fn scaffold_count(&self) -> u32 {
        self.scaffold_count
    }
    pub fn binned_scaffold_lengths(&self) -> &Vec<u64> {
        &self.binned_scaffold_lengths
    }
    pub fn binned_scaffold_counts(&self) -> &Vec<u32> {
        &self.binned_scaffold_counts
    }
    pub fn contigs(&self) -> &Vec<u64> {
        &self.contigs
    }
    pub fn contig_count(&self) -> u32 {
        self.contig_count
    }
    pub fn binned_contig_lengths(&self) -> &Vec<u64> {
        &self.binned_contig_lengths
    }
    pub fn binned_contig_counts(&self) -> &Vec<u32> {
        &self.binned_contig_counts
    }
}

pub fn snail_stats(options: &Config, assembly_stats: AssemblyStats) -> SnailStats {
    let span = assembly_stats.span();
    let segment = span / options.segments as u64;
    let contig_span = assembly_stats.contigs().iter().sum::<u64>();
    let contig_segment = contig_span / options.segments as u64;
    // TODO: check span > segments
    let mut position: u64 = 0;
    let mut contig_position: u64 = 0;
    let mut binned_gcs: Vec<SummaryStats> = vec![];
    let mut binned_ns: Vec<SummaryStats> = vec![];
    let mut scaffold_index: usize = 0;
    let mut scaffold_sum: u64 = assembly_stats.seq_stats()[0].length();
    let mut binned_scaffold_lengths: Vec<u64> = vec![];
    let mut binned_scaffold_counts: Vec<u32> = vec![];
    let mut contig_index: usize = 0;
    let mut contig_sum: u64 = assembly_stats.contigs()[0];
    let mut binned_contig_lengths: Vec<u64> = vec![];
    let mut binned_contig_counts: Vec<u32> = vec![];
    for _ in 0..options.segments {
        position += segment;
        contig_position += contig_segment;
        let mut gcs: Vec<f32> =
            vec![assembly_stats.seq_stats()[scaffold_index].gc_proportion() * 100.0];
        let mut ns: Vec<f32> =
            vec![assembly_stats.seq_stats()[scaffold_index].n_proportion() * 100.0];
        while scaffold_sum < position {
            scaffold_index += 1;
            scaffold_sum += assembly_stats.seq_stats()[scaffold_index].length();
            gcs.push(assembly_stats.seq_stats()[scaffold_index].gc_proportion() * 100.0);
            ns.push(assembly_stats.seq_stats()[scaffold_index].n_proportion() * 100.0);
        }
        binned_scaffold_counts.push(scaffold_index as u32 + 1);
        binned_scaffold_lengths.push(assembly_stats.seq_stats()[scaffold_index].length());
        while contig_sum < contig_position {
            contig_index += 1;
            contig_sum += assembly_stats.contigs()[contig_index];
        }
        binned_contig_counts.push(contig_index as u32 + 1);
        binned_contig_lengths.push(assembly_stats.contigs()[contig_index]);
        gcs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        binned_gcs.push(SummaryStats {
            min: gcs[0],
            max: gcs[gcs.len() - 1],
            mean: gcs.iter().sum::<f32>() / gcs.len() as f32,
        });
        ns.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        binned_ns.push(SummaryStats {
            min: ns[0],
            max: ns[ns.len() - 1],
            mean: ns.iter().sum::<f32>() / ns.len() as f32,
        });
    }
    SnailStats {
        span,
        atgc: assembly_stats.atgc(),
        gc: assembly_stats.gc(),
        gc_proportion: assembly_stats.gc() as f32 / assembly_stats.atgc() as f32 * 100.0,
        n: assembly_stats.n(),
        binned_gcs,
        binned_ns,
        scaffolds: vec![assembly_stats.seq_stats()[0].length()],
        scaffold_count: assembly_stats.n_seqs(),
        binned_scaffold_lengths,
        binned_scaffold_counts,
        contigs: vec![assembly_stats.contigs()[0]],
        contig_count: assembly_stats.n_contigs(),
        binned_contig_lengths,
        binned_contig_counts,
    }
}
