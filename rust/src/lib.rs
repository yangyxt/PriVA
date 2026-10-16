//! Source-mapped PriVA modules; production remains on the reference backend.
pub mod acmg_scoring;
mod cli;
pub mod spliceai2;

pub fn main_entry() {
    if let Err(error) = cli::run() {
        eprintln!("priva: {error}");
        std::process::exit(2);
    }
}
