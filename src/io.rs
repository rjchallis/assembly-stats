use std::sync::mpsc;

use bio::io::fasta;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use serde_json;

use crate::cli::Config;
use crate::seq_stats::{seq_stats, SeqStats};
use crate::snail::SnailStats;

pub fn process_fasta(options: &Config) -> Vec<SeqStats> {
    // iterate over FASTA file to get number of sequences
    let mut n_seqs = 0;
    let mut reader = fasta::Reader::from_file(&options.fasta)
        .expect("[-]\tPath invalid.")
        .records();
    while let Some(Ok(_record)) = reader.next() {
        n_seqs += 1;
    }

    // set up progress bar
    let progress_bar = ProgressBar::new(n_seqs);
    let pb_style_result = ProgressStyle::with_template(
        "[+]\tProcessing records: {bar:40.cyan/blue} {pos:>7}/{len:12}",
    );
    let pb_style = match pb_style_result {
        Ok(style) => style,
        Err(error) => panic!("Problem with the progress bar: {:?}", error),
    };
    // pb_style.progress_chars(">>-");
    progress_bar.set_style(pb_style);

    // create channel for collecting output
    let (sender, receiver) = mpsc::channel();

    // iterate over FASTA sequences
    let reader = fasta::Reader::from_file(&options.fasta).expect("[-]\tPath invalid.");
    reader
        .records()
        .par_bridge()
        .for_each_with(sender, |s, record| {
            let fasta_record = record.expect("[-]\tError during fasta record parsing.");

            if let Err(e) = s.send(seq_stats(options, fasta_record)) {
                eprintln!("Writing error: {}", e.to_string());
            }
            progress_bar.inc(1);
        });
    progress_bar.finish();
    let entries: Vec<SeqStats> = receiver.iter().collect();

    entries
}

pub fn write_snail_json(options: &Config, snail_stats: SnailStats) -> serde_json::Result<()> {
    let snail_json: String = match options.pretty_print {
        true => serde_json::to_string_pretty(&snail_stats)?,
        _ => serde_json::to_string(&snail_stats)?,
    };
    println!("{}", snail_json);

    Ok(())
}
