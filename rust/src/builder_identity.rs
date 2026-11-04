//! Per-builder source identities for installer receipts. Each identity covers
//! only the named builder's own modules and Cargo.lock, embedded by rustc, so it
//! does not change with unrelated code, executable path or build directory.
use crate::cli::{Result, invalid};

type Sources = &'static [(&'static str, &'static [u8])];

macro_rules! sources {
    ($($path:literal),+ $(,)?) => {
        &[$(($path, include_bytes!($path).as_slice())),+]
    };
}

// The same files, in the same order, as alpha_domain::builder_identity.
const ALPHA_DOMAIN: Sources = sources![
    "alpha_domain/mod.rs",
    "alpha_domain/association.rs",
    "alpha_domain/key.rs",
    "alpha_domain/key_parquet.rs",
    "alpha_domain/model.rs",
    "alpha_domain/metadata.rs",
    "alpha_domain/parquet.rs",
    "alpha_domain/store.rs",
    "alpha_domain/vcf.rs",
    "vcf_variant.rs",
    "vcf_text.rs",
    "tabix_index.rs",
    "clinvar_aa_splice_builder.rs",
    "clingen_map_builder.rs",
    "python_text.rs",
    "json_text.rs",
    "builder_error.rs",
    "../Cargo.lock",
];
// The same files, in the same order, as alpha_domain::exon_history::identity.
const ALPHA_EXON_HISTORY: Sources = sources![
    "alpha_domain/exon_history/mod.rs",
    "alpha_domain/exon_history/model.rs",
    "alpha_domain/exon_history/store.rs",
    "alpha_domain/mod.rs",
    "alpha_domain/vcf.rs",
    "alpha_domain/model.rs",
    "alpha_domain/key.rs",
    "alpha_domain/key_parquet.rs",
    "alpha_domain/parquet.rs",
    "python_set_order.rs",
    "vcf_variant.rs",
    "vcf_text.rs",
    "tabix_index.rs",
    "clinvar_aa_splice_builder.rs",
    "clingen_map_builder.rs",
    "python_text.rs",
    "json_text.rs",
    "builder_error.rs",
    "../Cargo.lock",
];

static BUILDERS: [(&str, Sources); 10] = [
    ("alphamissense-domain", ALPHA_DOMAIN),
    ("domain-metadata", ALPHA_DOMAIN),
    ("alpha-exon-history", ALPHA_EXON_HISTORY),
    (
        "das-domains",
        sources![
            "das_domain.rs",
            "das_assembly.rs",
            "das_benign.rs",
            "das_benign_shards.rs",
            "das_coldspot_filter.rs",
            "das_coldspot_table.rs",
            "das_common.rs",
            "das_integer_ratio.rs",
            "das_reference.rs",
            "das_rescue.rs",
            "das_rescue_plan.rs",
            "das_sampling.rs",
            "dip_stat.rs",
            "ks_greater.rs",
            "exon_af_stats.rs",
            "fdr.rs",
            "numpy_rng.rs",
            "alpha_domain/das_projection.rs",
            "alpha_domain/das_transcripts.rs",
            "alpha_domain/mod.rs",
            "alpha_domain/model.rs",
            "alpha_domain/store.rs",
            "alpha_domain/parquet.rs",
            "alpha_domain/key.rs",
            "alpha_domain/key_parquet.rs",
            "alpha_domain/metadata.rs",
            "alpha_domain/vcf.rs",
            "protein_domain_mapping.rs",
            "interpro_mapping_io.rs",
            "cache_convert/interpro_mapping.rs",
            "cache_convert/intolerant_domains.rs",
            "cache_convert/transcript_domains.rs",
            "cache_convert/common.rs",
            "rust_cache.rs",
            "python_set_order.rs",
            "pandas_table.rs",
            "pandas_low_memory.rs",
            "vcf_variant.rs",
            "vcf_text.rs",
            "tabix_index.rs",
            "python_text.rs",
            "json_text.rs",
            "builder_error.rs",
            "../Cargo.lock",
        ],
    ),
    (
        "pm1-regions",
        sources![
            "pm1_regions_builder.rs",
            "cache_convert/pm1_regions.rs",
            "cache_convert/intolerant_domains.rs",
            "cache_convert/common.rs",
            "rust_cache.rs",
            "fdr.rs",
            "python_args.rs",
            "python_text.rs",
            "../Cargo.lock",
        ],
    ),
    (
        "clinvar-gene-stats",
        sources![
            "clinvar_gene_stat_builder.rs",
            "clinvar_gene_stat_builder/maps.rs",
            "clinvar_gene_stat_builder/statistics.rs",
            "cache_convert/clinvar_gene_stat.rs",
            "cache_convert/clinvar_gene_stat/typed.rs",
            "cache_convert/common.rs",
            "rust_cache.rs",
            "combine_annotations/gene_stats.rs",
            "clingen_map_builder.rs",
            "spliceai_cutoffs.rs",
            "pickle_reader.rs",
            "pandas_table.rs",
            "python_set_order.rs",
            "vcf_variant.rs",
            "vcf_text.rs",
            "python_text.rs",
            "json_text.rs",
            "builder_error.rs",
            "../Cargo.lock",
        ],
    ),
    (
        "clinvar-gene-af",
        sources![
            "clinvar_gene_af_builder.rs",
            "cache_convert/clinvar_gene_af.rs",
            "cache_convert/clinvar_gene_af/typed.rs",
            "cache_convert/common.rs",
            "rust_cache.rs",
            "tabix_index.rs",
            "vcf_variant.rs",
            "python_text.rs",
            "builder_error.rs",
            "../Cargo.lock",
        ],
    ),
    (
        "clinvar-exon-af",
        sources![
            "clinvar_exon_af_builder.rs",
            "exon_af_stats.rs",
            "cache_convert/clinvar_exon_af.rs",
            "cache_convert/common.rs",
            "rust_cache.rs",
            "clingen_map_builder.rs",
            "tabix_index.rs",
            "vcf_variant.rs",
            "vcf_text.rs",
            "python_text.rs",
            "builder_error.rs",
            "../Cargo.lock",
        ],
    ),
    (
        "clinvar-aa-splice",
        sources![
            "clinvar_aa_splice_builder.rs",
            "cache_convert.rs",
            "cache_convert/clinvar_splice_stat.rs",
            "cache_convert/clinvar_splice_stat/typed.rs",
            "cache_convert/common.rs",
            "rust_cache.rs",
            "tabix_index.rs",
            "vcf_variant.rs",
            "vcf_text.rs",
            "python_text.rs",
            "builder_error.rs",
            "../Cargo.lock",
        ],
    ),
    (
        "clingen-map",
        sources![
            "clingen_map_builder.rs",
            "cache_convert/clingen_map.rs",
            "cache_convert/clingen_map/typed.rs",
            "cache_convert/common.rs",
            "rust_cache.rs",
            "pandas_polars.rs",
            "pandas_table.rs",
            "python_args.rs",
            "python_text.rs",
            "vcf_text.rs",
            "../Cargo.lock",
        ],
    ),
];

const HELP: &str = "usage: priva builder-identity NAME\n\
Prints {\"builder\":NAME,\"builder_source_sha256\":SHA256}: SHA-256 over the\n\
builder's listed sources in order, each preceded by its byte length as a\n\
little-endian u64. Sources are embedded at compile time, so the value does not\n\
depend on executable path, debug symbols or build directory. Native\n\
intermediates keep the identity recorded when they were built and are rebuilt\n\
only on explicit request; their checks accept an earlier identity only when it\n\
is pinned. domain-metadata shares the alphamissense-domain builder and value.\n\
NAME and covered sources (relative to rust/):";

// Length-prefixed so that moving bytes between adjacent sources changes the value.
pub(crate) fn digest<'a>(sources: impl IntoIterator<Item = &'a [u8]>) -> String {
    let mut sha = crate::sha256::Sha256::new();
    for source in sources {
        sha.update(&(source.len() as u64).to_le_bytes());
        sha.update(source);
    }
    sha.finish_hex()
}

pub(crate) fn identity(name: &str) -> Option<String> {
    let (_, sources) = BUILDERS.iter().find(|(builder, _)| *builder == name)?;
    Some(digest(sources.iter().map(|(_, bytes)| *bytes)))
}

fn help() -> String {
    let mut text = String::from(HELP);
    for (name, sources) in &BUILDERS {
        text.push_str("\n  ");
        text.push_str(name);
        let mut line = String::from("   ");
        for (path, _) in *sources {
            let path = path
                .strip_prefix("../")
                .map_or_else(|| format!("src/{path}"), String::from);
            if line.len() + path.len() > 79 {
                text.push('\n');
                text.push_str(&line);
                line = String::from("   ");
            }
            line.push(' ');
            line.push_str(&path);
        }
        text.push('\n');
        text.push_str(&line);
    }
    text
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args.iter().any(|s| s == "--help" || s == "-h") {
        println!("{}", help());
        return Ok(());
    }
    let [name] = args else {
        return Err(invalid(help()));
    };
    let value = identity(name).ok_or_else(|| invalid(help()))?;
    println!("{{\"builder\":\"{name}\",\"builder_source_sha256\":\"{value}\"}}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{BUILDERS, digest, dispatch, identity};
    use std::collections::HashMap;

    #[test]
    fn identities_are_stable_and_differ_between_builders() {
        let mut seen = HashMap::new();
        for (name, _) in &BUILDERS {
            let value = identity(name).unwrap();
            assert_eq!(identity(name).unwrap(), value);
            assert_eq!(value.len(), 64);
            assert!(value.bytes().all(|b| b.is_ascii_hexdigit()));
            if let Some(other) = seen.insert(value, *name) {
                // The two dispatches of one builder record one footer identity.
                assert_eq!((other, *name), ("alphamissense-domain", "domain-metadata"));
            }
        }
        assert_eq!(seen.len(), BUILDERS.len() - 1);
    }

    #[test]
    fn alpha_identities_equal_the_values_recorded_in_footers() {
        let domain = crate::alpha_domain::builder_identity().unwrap();
        assert_eq!(identity("alphamissense-domain").unwrap(), domain);
        assert_eq!(identity("domain-metadata").unwrap(), domain);
        let history = crate::alpha_domain::exon_history::identity();
        assert_eq!(identity("alpha-exon-history").unwrap(), history);
    }

    #[test]
    fn digest_binds_source_boundaries_and_unknown_names_are_refused() {
        assert_ne!(
            digest([b"ab".as_slice(), b"c"]),
            digest([b"a".as_slice(), b"bc"])
        );
        assert_eq!(identity("cadd-hub"), None);
        assert!(dispatch(&["cadd-hub".to_owned()]).is_err());
        assert!(dispatch(&[]).is_err());
        assert!(dispatch(&["pm1-regions".to_owned()]).is_ok());
    }
}
