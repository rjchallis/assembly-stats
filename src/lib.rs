use cli::Config;

use std::error::Error;

pub mod assembly;
pub mod cli;
pub mod io;
pub mod seq_stats;
pub mod snail;

pub fn run(options: Config) -> Result<(), Box<dyn Error>> {
    let entries = io::process_fasta(&options);
    let assembly_stats = assembly::assembly_stats(entries);
    let snail_stats = snail::snail_stats(&options, assembly_stats);
    io::write_snail_json(&options, snail_stats)?;
    Ok(())
}
