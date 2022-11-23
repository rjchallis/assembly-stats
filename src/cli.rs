use clap::Parser;
use clap_num::number_range;

fn up_to_100k(s: &str) -> Result<u32, String> {
    number_range(s, 1, 100000)
}

#[derive(Parser)]
#[command(author, version, about, long_about = None)] // Read from `Cargo.toml`
pub struct Config {
    /// FASTA file to generate snail plot for
    #[arg(long, short = 'f')]
    pub fasta: String,
    /// Number of segments to divide assembly into (1-100000)
    #[arg(long, short = 's', default_value_t = 1000, value_parser=up_to_100k)]
    pub segments: u32,
    /// Length of run of Ns to use to split scaffolds into contigs
    #[arg(long = "split-contigs", short = 'c', default_value_t = 10)]
    pub split_contigs: u16,
    /// JSON file to write output to
    #[arg(long, short = 'j', default_value_t = String::from("output.json"))]
    pub json: String,
    /// Flag to pretty print json output
    #[arg(long = "pretty-print", short = 'p', default_value_t = false)]
    pub pretty_print: bool,
}

pub fn parse() -> Config {
    Config::parse()
}
