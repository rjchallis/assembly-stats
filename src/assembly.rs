use crate::seq_stats::SeqStats;

pub struct AssemblyStats {
    seq_stats: Vec<SeqStats>,
    contigs: Vec<u64>,
    longest_contig: u64,
    longest_scaffold: u64,
    span: u64,
    n_seqs: u32,
    n_contigs: u32,
    at: u64,
    gc: u64,
    atgc: u64,
    n: u64,
}

impl AssemblyStats {
    pub fn seq_stats(&self) -> &Vec<SeqStats> {
        &self.seq_stats
    }
    pub fn span(&self) -> u64 {
        self.span
    }
    pub fn n_seqs(&self) -> u32 {
        self.n_seqs
    }
    pub fn n_contigs(&self) -> u32 {
        self.n_contigs
    }
    pub fn contigs(&self) -> &Vec<u64> {
        &self.contigs
    }
    pub fn longest_contig(&self) -> u64 {
        self.longest_contig
    }
    pub fn longest_scaffold(&self) -> u64 {
        self.longest_scaffold
    }
    pub fn at(&self) -> u64 {
        self.at
    }
    pub fn gc(&self) -> u64 {
        self.gc
    }
    pub fn atgc(&self) -> u64 {
        self.atgc
    }
    pub fn n(&self) -> u64 {
        self.n
    }
}

pub fn assembly_stats(mut entries: Vec<SeqStats>) -> AssemblyStats {
    // sort entries by sequence length
    entries.sort_by_key(|x| x.length);
    entries.reverse();
    let mut span: u64 = 0;
    let mut at: u64 = 0;
    let mut gc: u64 = 0;
    let mut atgc: u64 = 0;
    let mut n: u64 = 0;
    let mut n_contigs: u32 = 0;
    let longest_scaffold: u64 = entries[0].length();
    let mut contigs: Vec<u64> = vec![];
    for entry in &entries {
        span += entry.length();
        n_contigs += entry.contigs().len() as u32;
        at += entry.at();
        gc += entry.gc();
        atgc += entry.atgc();
        n += entry.n();
        contigs.extend(entry.contigs().iter().cloned());
    }
    contigs.sort();
    contigs.reverse();
    let longest_contig: u64 = contigs[0];

    let n_seqs = entries.len() as u32;
    AssemblyStats {
        seq_stats: entries,
        contigs,
        longest_contig,
        longest_scaffold,
        span,
        n_seqs,
        n_contigs,
        at,
        gc,
        atgc,
        n,
    }
}
