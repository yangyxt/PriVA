//! Source-mapped PriVA modules; production remains on the reference backend.
pub mod acmg_consequence;
pub mod acmg_scoring;
pub mod acmg_variant_mechanism;
pub mod insilico_criteria;
mod cli;
pub mod pedigree_filtration;
pub mod pm4_bp3;
pub mod spliceai2;

pub fn main_entry() {
    if let Err(error) = cli::run() {
        eprintln!("priva: {error}");
        std::process::exit(2);
    }
}
