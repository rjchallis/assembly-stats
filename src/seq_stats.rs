use std::collections::HashMap;

use bio::io::fasta;

use crate::cli::Config;

pub struct SeqStats {
    at_count: u64,
    gc_count: u64,
    atgc_count: u64,
    n_count: u64,
    gc_proportion: f32,
    at_proportion: f32,
    n_proportion: f32,
    contigs: Vec<u64>,
    pub length: u64,
}

impl SeqStats {
    pub fn at(&self) -> u64 {
        self.at_count
    }
    pub fn gc(&self) -> u64 {
        self.gc_count
    }
    pub fn atgc(&self) -> u64 {
        self.atgc_count
    }
    pub fn n(&self) -> u64 {
        self.n_count
    }
    pub fn gc_proportion(&self) -> f32 {
        self.gc_proportion
    }
    pub fn at_proportion(&self) -> f32 {
        self.at_proportion
    }
    pub fn n_proportion(&self) -> f32 {
        self.n_proportion
    }
    pub fn contigs(&self) -> &Vec<u64> {
        &self.contigs
    }
    pub fn length(&self) -> u64 {
        self.length
    }
}

/// Count frequency of each base in a nucleotide sequence.
fn nucleotide_counts(dna: &[u8]) -> HashMap<&u8, u64> {
    let mut map = HashMap::new();
    for nucleotide in dna {
        let count = map.entry(nucleotide).or_insert(0);
        *count += 1;
    }
    map
}

pub fn multichar_count(counts: &HashMap<&u8, u64>, chars: &[u8]) -> u64 {
    let mut count: u64 = 0;
    for character in chars {
        count += counts.get(character).unwrap_or(&0)
    }
    count
}

pub fn extract_contigs(seq: &[u8], break_len: u16) -> Vec<u64> {
    let mut contigs: Vec<u64> = vec![];
    let mut gap = 0;
    let mut length = 0;
    for base in seq {
        match base {
            b'n' => gap += 1,
            b'N' => gap += 1,
            _ => {
                gap = 0;
                length += 1
            }
        }
        if gap == break_len {
            if length > 0 {
                contigs.push(length);
                length = 0;
            }
        }
    }
    if length > 0 {
        contigs.push(length);
    }
    contigs
}

pub fn seq_stats(options: &Config, fasta_record: fasta::Record) -> SeqStats {
    let seq = fasta_record.seq();
    let length: u64 = seq.len() as u64;
    let counts = nucleotide_counts(seq);
    let contigs = extract_contigs(&seq, options.split_contigs);
    // match upper or lower case bases
    let a_count = multichar_count(&counts, &[b'a', b'A']);
    let c_count = multichar_count(&counts, &[b'c', b'C']);
    let g_count = multichar_count(&counts, &[b'g', b'G']);
    let t_count = multichar_count(&counts, &[b't', b'T']);
    let n_count = multichar_count(&counts, &[b'n', b'N']);
    let s_count = multichar_count(&counts, &[b'S', b'S']);
    let w_count = multichar_count(&counts, &[b'w', b'W']);
    // count remaining ambiguous bases in a single bin
    let amb_count = multichar_count(
        &counts,
        &[
            b'b', b'B', b'd', b'D', b'h', b'H', b'k', b'K', b'm', b'M', b'r', b'R', b'v', b'V',
            b'y', b'Y',
        ],
    );

    // include S in GC count
    let gc_count = c_count + g_count + s_count;
    // include W in AT count
    let at_count = a_count + t_count + w_count;
    // include W and S in ACGT count
    let atgc_count = at_count + gc_count;
    // treat all other ambiguous bases as N
    let other_count = n_count + amb_count;
    SeqStats {
        at_count,
        gc_count,
        atgc_count,
        n_count: other_count,
        gc_proportion: ((gc_count) as f32 / atgc_count as f32),
        at_proportion: ((at_count) as f32 / atgc_count as f32),
        n_proportion: (other_count as f32 / length as f32),
        contigs,
        length,
    }
}

#[cfg(test)]
mod tests {

    use super::nucleotide_counts;

    const A: u8 = b'A';
    const C: u8 = b'C';
    const G: u8 = b'G';
    const T: u8 = b'T';

    #[test]
    fn test_nucleotide_counts() {
        let short_dna_string = "AAaCCcTTtGGgg".as_bytes();

        let nuc_counts = nucleotide_counts(short_dna_string);

        // Two of each upper case character
        assert_eq!(2, *nuc_counts.get(&A).unwrap());
        assert_eq!(2, *nuc_counts.get(&C).unwrap());
        assert_eq!(2, *nuc_counts.get(&G).unwrap());
        assert_eq!(2, *nuc_counts.get(&T).unwrap());
    }
}
