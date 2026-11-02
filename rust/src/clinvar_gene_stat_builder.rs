//! Source procedures for `scripts/clinvar_stat_variants.py`.
use crate::cli::{Result, invalid};
use crate::python_set_order::{StringHash, StringSet};
use std::collections::HashMap;
use std::fs;
use std::path::Path;

mod maps;
mod statistics;

const BUILD_USAGE: &str = "usage: priva clinvar-gene-stats build --vcf VCF --output RKYV --python-hash-seed UINT32 \\\n    (--clinvar-summary TSV --mane-summary TSV | --allele-transcript-map PICKLE) \\\n    [--threads N] [--config YAML | --donor-loss-cutoff FLOAT] \\\n    [--completion-order TRACE_JSON] [--completion-trace NEW_TRACE_JSON]";

const MERGE_USAGE: &str = "usage: priva clinvar-gene-stats merge --input SHARD_RKYV [--input SHARD_RKYV ...] --output RKYV";

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let usage = "usage: priva clinvar-gene-stats set-order --input LINES --output LINES --python-hash-seed UINT32 [--hashes OUTPUT]";
    if args.is_empty() || args == ["--help"] {
        println!("{BUILD_USAGE}\n{MERGE_USAGE}\n{usage}");
        return Ok(());
    }
    if args[0] == "build" {
        return build(&args[1..]);
    }
    if args[0] == "merge" {
        return merge(&args[1..]);
    }
    if args[0] != "set-order" {
        return Err(invalid(usage));
    }
    let mut options = HashMap::new();
    let mut at = 1;
    while at < args.len() {
        let name = args[at].as_str();
        if !["--input", "--output", "--python-hash-seed", "--hashes"].contains(&name) {
            return Err(invalid(format!("unknown option {name}\n{usage}")));
        }
        let value = args.get(at + 1).ok_or_else(|| invalid(usage))?;
        if options.insert(name, value.as_str()).is_some() {
            return Err(invalid(format!("duplicate option {name}")));
        }
        at += 2;
    }
    let required = |name: &str| options.get(name).copied().ok_or_else(|| invalid(usage));
    let seed = required("--python-hash-seed")?.parse::<u32>()?;
    let input = required("--input")?;
    let output = required("--output")?;
    crate::rust_cache::ensure_distinct_paths(Path::new(input), Path::new(output))?;
    if let Some(path) = options.get("--hashes") {
        crate::rust_cache::ensure_distinct_paths(Path::new(input), Path::new(path))?;
        distinct_destinations(Path::new(output), Path::new(path))?;
    }
    let text = fs::read_to_string(input)?;
    let mut set = StringSet::new(seed);
    let hash = StringHash::seeded(seed);
    let mut hashes = String::new();
    for line in text.lines() {
        set.insert(line.to_owned());
        hashes.push_str(&format!("{}\n", hash.hash(line)));
    }
    let order = set.into_order();
    let rendered = if order.is_empty() {
        String::new()
    } else {
        order.join("\n") + "\n"
    };
    fs::write(output, rendered)?;
    if let Some(path) = options.get("--hashes") {
        fs::write(path, hashes)?;
    }
    Ok(())
}

/// Collect only disjoint completed gene partitions. This does not combine
/// per-gene statistics from overlapping shards or choose a new symbol/order.
fn merge(args: &[String]) -> Result<()> {
    use crate::cache_convert::{Provenance, clinvar_gene_stat};
    let mut inputs = Vec::new();
    let mut output = None;
    if args.len() % 2 != 0 {
        return Err(invalid(MERGE_USAGE));
    }
    for pair in args.chunks_exact(2) {
        match pair[0].as_str() {
            "--input" => inputs.push(Path::new(&pair[1])),
            "--output" if output.is_none() => output = Some(Path::new(&pair[1])),
            _ => return Err(invalid(MERGE_USAGE)),
        }
    }
    let output = output.ok_or_else(|| invalid(MERGE_USAGE))?;
    if inputs.is_empty() {
        return Err(invalid(MERGE_USAGE));
    }
    let mut entries = std::collections::BTreeMap::new();
    let mut provenance = vec![Provenance {
        key: "producer".into(),
        value: "priva clinvar-gene-stats merge (disjoint genes)".into(),
    }];
    for (index, input) in inputs.iter().enumerate() {
        crate::rust_cache::ensure_distinct_paths(input, output)?;
        let path = input
            .to_str()
            .ok_or_else(|| invalid("cache path is not UTF-8"))?;
        let cache = clinvar_gene_stat::ClinvarGeneCache::open(path)?;
        for key in cache.keys() {
            if entries
                .insert(key.to_owned(), cache.plain(key).unwrap())
                .is_some()
            {
                return Err(invalid(
                    "gene partitions overlap; a disjoint merge cannot reproduce the source group",
                ));
            }
        }
        provenance.push(Provenance {
            key: format!("input_{index}_sha256"),
            value: crate::sha256::sha256_file(input)?,
        });
    }
    let count = entries.len();
    let version =
        clinvar_gene_stat::write_entries(output, entries.into_iter().collect(), provenance)?;
    println!(
        "{{\"schema\":5,\"schema_version\":{version},\"shards\":{},\"genes\":{count}}}",
        inputs.len()
    );
    Ok(())
}

fn build(args: &[String]) -> Result<()> {
    use crate::json_text::{self, Json};
    let mut options = HashMap::new();
    if args.len() % 2 != 0 {
        return Err(invalid(BUILD_USAGE));
    }
    for pair in args.chunks_exact(2) {
        let name = pair[0].as_str();
        if ![
            "--vcf",
            "--output",
            "--python-hash-seed",
            "--clinvar-summary",
            "--mane-summary",
            "--allele-transcript-map",
            "--threads",
            "--config",
            "--donor-loss-cutoff",
            "--completion-order",
            "--completion-trace",
        ]
        .contains(&name)
            || options.insert(name, pair[1].as_str()).is_some()
        {
            return Err(invalid(format!(
                "unknown or repeated option {name}\n{BUILD_USAGE}"
            )));
        }
    }
    let required = |name: &str| {
        options
            .get(name)
            .copied()
            .ok_or_else(|| invalid(BUILD_USAGE))
    };
    let vcf = Path::new(required("--vcf")?);
    let output = Path::new(required("--output")?);
    let seed = required("--python-hash-seed")?.parse::<u32>()?;
    let threads = options.get("--threads").unwrap_or(&"10").parse::<usize>()?;
    let prepared = options.get("--allele-transcript-map").copied();
    if prepared.is_some()
        && (options.contains_key("--clinvar-summary") || options.contains_key("--mane-summary"))
    {
        return Err(invalid(
            "supply the summary/MANE pair or a saved target map, not both",
        ));
    }
    if prepared.is_none() {
        required("--clinvar-summary")?;
        required("--mane-summary")?;
    }
    if options.contains_key("--config") && options.contains_key("--donor-loss-cutoff") {
        return Err(invalid(
            "supply a config or an explicit procedure cutoff, not both",
        ));
    }
    let cutoff = match options.get("--donor-loss-cutoff") {
        Some(value) => value.parse::<f64>()?,
        None => {
            crate::spliceai_cutoffs::cutoffs_for_config(options.get("--config").map(Path::new))
                .map_err(|failure| invalid(format!("SpliceAI cutoffs: {failure:?}")))?
                .donor_loss_putative_nmd
        }
    };
    let trace_output = options.get("--completion-trace").map(Path::new);
    let mut protected: Vec<_> = [
        "--vcf",
        "--clinvar-summary",
        "--mane-summary",
        "--allele-transcript-map",
        "--config",
        "--completion-order",
    ]
    .iter()
    .filter_map(|name| options.get(name).map(std::path::PathBuf::from))
    .collect();
    for suffix in [".tbi", ".csi"] {
        let index = std::path::PathBuf::from(format!("{}{suffix}", vcf.display()));
        if index.exists() {
            protected.push(index);
        }
    }
    for source in &protected {
        crate::rust_cache::ensure_distinct_paths(source, output)?;
        if let Some(trace) = trace_output {
            crate::rust_cache::ensure_distinct_paths(source, trace)?;
        }
    }
    if let Some(trace) = trace_output {
        distinct_destinations(output, trace)?;
        if trace.exists() {
            return Err(invalid("the completion-trace destination already exists"));
        }
    }
    let vcf_hash = if options.contains_key("--completion-order") || trace_output.is_some() {
        Some(crate::sha256::sha256_file(vcf)?)
    } else {
        None
    };
    let order = if let Some(path) = options.get("--completion-order") {
        let value = json_text::parse(&fs::read_to_string(path)?, true).map_err(invalid)?;
        if value.get("vcf_sha256").and_then(Json::as_str) != vcf_hash.as_deref() {
            return Err(invalid(
                "the completion trace does not name these VCF bytes",
            ));
        }
        let Some(Json::Array(order)) = value.get("order") else {
            return Err(invalid("the completion trace has no order array"));
        };
        Some(
            order
                .iter()
                .map(|value| match value {
                    Json::Int(value) => {
                        usize::try_from(*value).map_err(|_| invalid("invalid completion ordinal"))
                    }
                    _ => Err(invalid("completion ordinals must be integers")),
                })
                .collect::<Result<Vec<_>>>()?,
        )
    } else {
        None
    };
    let started = std::time::Instant::now();
    let (table, trace) =
        crate::combine_annotations::gene_stats::read(vcf, threads, order.as_deref())?;
    let rows = table.rows;
    let targets = match prepared {
        Some(path) => maps::from_pickle(Path::new(path))?,
        None => maps::from_tables(
            Path::new(required("--clinvar-summary")?),
            Path::new(required("--mane-summary")?),
        )?,
    };
    let alleles = maps::alleles(vcf)?;
    let entries = statistics::build(&table, &targets, &alleles, seed, cutoff)?;
    drop(table);
    drop(targets);
    drop(alleles);
    if let Some(path) = trace_output {
        use std::io::Write;
        let value = Json::Object(vec![
            ("vcf_sha256".into(), Json::Str(vcf_hash.unwrap())),
            ("python_hash_seed".into(), Json::Int(i128::from(seed))),
            (
                "order".into(),
                Json::Array(trace.iter().map(|at| Json::Int(*at as i128)).collect()),
            ),
        ]);
        let mut file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)?;
        file.write_all(json_text::dumps(&value, false).map_err(invalid)?.as_bytes())?;
        file.write_all(b"\n")?;
    }
    let count = entries.len();
    let provenance = vec![
        crate::cache_convert::Provenance {
            key: "producer".into(),
            value: "priva clinvar-gene-stats build".into(),
        },
        crate::cache_convert::Provenance {
            key: "python_hash_seed".into(),
            value: seed.to_string(),
        },
        crate::cache_convert::Provenance {
            key: "donor_loss_cutoff".into(),
            value: crate::python_text::float_repr(cutoff),
        },
        crate::cache_convert::Provenance {
            key: "collector_order".into(),
            value: if order.is_some() {
                "explicit completion trace"
            } else if threads == 1 {
                "source order"
            } else {
                "unordered worker completion"
            }
            .into(),
        },
    ];
    let version =
        crate::cache_convert::clinvar_gene_stat::write_entries(output, entries, provenance)?;
    println!(
        "{{\"schema\":5,\"schema_version\":{version},\"records\":{},\"annotation_rows\":{rows},\"genes\":{count},\"seconds\":{:.3}}}",
        trace.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

/// Both destinations may be new. The ordinary source/output guard requires its
/// source to exist, so compare canonical parent/name pairs in this one case.
fn distinct_destinations(first: &Path, second: &Path) -> Result<()> {
    if first.exists() {
        return crate::rust_cache::ensure_distinct_paths(first, second);
    }
    if second.exists() {
        return crate::rust_cache::ensure_distinct_paths(second, first);
    }
    let destination = |path: &Path| -> Result<std::path::PathBuf> {
        let parent = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let name = path
            .file_name()
            .ok_or_else(|| invalid("the output needs a file name"))?;
        Ok(fs::canonicalize(parent)?.join(name))
    };
    if destination(first)? == destination(second)? {
        return Err(invalid(
            "the archive and diagnostic destinations must differ",
        ));
    }
    Ok(())
}
