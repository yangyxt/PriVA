//! Source-mapped PriVA modules; production remains on the reference backend.
pub mod acmg_consequence;
pub mod acmg_scoring;
mod acmg_step;
pub mod acmg_variant_mechanism;
mod alpha_domain;
mod alphamissense_vcf;
mod builder_error;
mod builder_table;
mod cache_convert;
mod cache_install;
mod cache_resources;
pub mod cadd_hub;
mod cli;
mod clingen_map_builder;
mod clinvar_aa_splice_builder;
mod clinvar_exon_af_builder;
mod clinvar_gene_af_builder;
pub mod combine_annotations;
mod config_text;
mod conservation_resource;
mod das_input_continuity;
mod das_install_state;
mod exon_af_stats;
#[allow(dead_code)] // PM1 producer integration follows exact float64 comparison.
mod fdr;
// Qualified before the opt-in DAS procedures consume this exact random stream.
#[allow(dead_code)]
mod numpy_rng;
// Qualified independently before the opt-in DAS producer is connected.
#[allow(dead_code)]
mod ks_greater;
#[allow(dead_code)] // Shared installer primitive; callers switch after their own parity tests.
mod file_lock;
mod gzip_writer;
mod hpo_assertions;
pub mod hpo_penetrance;
mod interpro_mapping_builder;
pub mod interpro_mapping_io;
mod json_text;
#[allow(
    dead_code,
    reason = "Publisher-checksum API awaiting installer-builder callers"
)]
mod md5;
mod mechanism_source;
mod mondo_scope;
mod pandas_polars;
mod pandas_table;
pub mod pedigree_filter;
mod pickle_reader;
mod pm1_regions_builder;
pub mod protein_domain_mapping;
mod python_args;
mod python_io;
mod python_lines;
mod python_path;
mod python_set_order;
mod python_text;
mod rust_cache;
mod sha256;
pub mod spliceai2;
pub mod spliceai_cache_input;
pub mod spliceai_cutoffs;
mod splicevault_builder;
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

mod vcf_variant;
