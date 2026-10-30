//! Source-mapped PriVA modules; production remains on the reference backend.
pub mod acmg_consequence;
pub mod acmg_scoring;
pub mod acmg_variant_mechanism;
pub mod cadd_hub;
mod cli;
pub mod combine_annotations;
pub mod hpo_penetrance;
pub mod interpro_mapping_io;
mod json_text;
mod pandas_polars;
mod pandas_table;
pub mod pedigree_filter;
mod pickle_reader;
pub mod protein_domain_mapping;
mod python_args;
mod python_io;
mod python_text;
mod sha256;
pub mod spliceai2;
pub mod spliceai_cache_input;
pub mod spliceai_cutoffs;
mod splicing_cutoffs;
pub mod splicing_var_analysis;
mod tabix_index;
pub mod utr_anno_interpret;
pub mod variant_ids;
mod vcf_text;
mod yaml_config;

pub fn main_entry() {
    if let Err(error) = cli::run() {
        eprintln!("priva: {error}");
        std::process::exit(2);
    }
}
