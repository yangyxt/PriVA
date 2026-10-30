//! Source-mapped PriVA modules; production remains on the reference backend.
pub mod acmg_consequence;
pub mod acmg_scoring;
pub mod acmg_variant_mechanism;
pub mod cadd_hub;
mod cli;
pub mod combine_annotations;
pub mod hpo_penetrance;
pub mod interpro_mapping_io;
mod pandas_table;
pub mod pedigree_filter;
mod pickle_reader;
pub mod protein_domain_mapping;
mod python_text;
pub mod spliceai2;
mod splicing_cutoffs;
pub mod splicing_var_analysis;
pub mod variant_ids;
mod vcf_text;

pub fn main_entry() {
    if let Err(error) = cli::run() {
        eprintln!("priva: {error}");
        std::process::exit(2);
    }
}
