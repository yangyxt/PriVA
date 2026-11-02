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
pub mod cadd_hub;
mod cli;
mod clingen_map_builder;
mod clinvar_aa_splice_builder;
mod clinvar_exon_af_builder;
mod clinvar_gene_af_builder;
pub mod combine_annotations;
mod config_text;
mod conservation_resource;
mod exon_af_stats;
#[allow(dead_code)] // PM1 producer integration follows exact float64 comparison.
mod fdr;
// Qualified before the opt-in DAS procedures consume this exact random stream.
#[allow(dead_code)]
mod numpy_rng;
// Qualified independently before the opt-in DAS producer is connected.
#[allow(dead_code)]
mod ks_greater;
// These numerical procedures are qualified before any DAS caller is enabled.
#[allow(dead_code)]
mod das_sampling;
// Candidate classification is validated before the full DAS builder is enabled.
mod das_benign;
#[allow(dead_code)] // Coldspot procedures await the complete opt-in DAS builder.
mod das_coldspot;
#[allow(dead_code)] // Common-AF procedures await the complete opt-in DAS builder.
mod das_common;
mod das_domain;
#[allow(dead_code)] // Retained exon ordering is qualified before DAS consumes it.
mod das_exon_pickle;
#[allow(dead_code)] // Reference/Family procedures await the complete opt-in DAS builder.
mod das_reference;
#[allow(dead_code)]
mod das_rescue;
#[allow(dead_code)] // Qualified task ordering precedes full DAS execution.
mod das_rescue_plan;
#[allow(dead_code)]
mod dip_stat;
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
#[allow(dead_code)] // Individual opt-in producers use distinct text/pair interfaces.
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
