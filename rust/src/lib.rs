//! Source-mapped PriVA modules; production remains on the reference backend.
pub mod acmg_consequence;
pub mod acmg_scoring;
mod acmg_step;
pub mod acmg_variant_mechanism;
mod alphamissense_vcf;
mod builder_error;
mod builder_table;
mod cache_convert;
mod cache_install;
pub mod cadd_hub;
mod cli;
mod clingen_map_builder;
mod clinvar_gene_af_builder;
pub mod combine_annotations;
mod config_text;
#[allow(dead_code)] // PM1 producer integration follows exact float64 comparison.
mod fdr;
#[allow(dead_code)] // Shared installer primitive; callers switch after their own parity tests.
mod file_lock;
mod gzip_writer;
mod hpo_assertions;
pub mod hpo_penetrance;
pub mod interpro_mapping_io;
mod json_text;
#[allow(
    dead_code,
    reason = "Publisher-checksum API awaiting installer-builder callers"
)]
mod md5;
mod mondo_scope;
mod pandas_polars;
mod pandas_table;
pub mod pedigree_filter;
mod pickle_reader;
mod pm1_regions_builder;
pub mod protein_domain_mapping;
mod python_args;
mod python_io;
mod python_text;
mod rust_cache;
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
        std::process::exit(if error.is::<builder_error::RuntimeError>() {
            1
        } else {
            2
        });
    }
}
