use std::process;

use snail_plot;
use snail_plot::cli;

fn main() {
    let options = cli::parse();
    if let Err(e) = snail_plot::run(options) {
        println!("Application error: {e}");
        process::exit(1);
    }
}
