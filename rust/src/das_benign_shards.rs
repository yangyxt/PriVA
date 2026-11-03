//! Explicit indexed chromosome counts with complete, header-ordered merge.
//! The whole-file benign command remains unchanged. Shards retain occurrence
//! counts, never per-shard removal decisions; the threshold runs after summation.
use crate::alpha_domain;
use crate::alpha_domain::vcf::Context;
use crate::builder_error;
use crate::cli::{Result, invalid};
use crate::das_benign;
use crate::json_text::{self, Json};
use crate::sha256;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

const FORMAT: &str = "priva-das-benign-counts-v1";
const KEYS: [&str; 9] = [
    "format",
    "source_sha256",
    "input",
    "domains_sha256",
    "domain_rows",
    "contigs",
    "chromosome",
    "records",
    "counts",
];
const HELP: &str = "usage: priva das-domains benign-chromosome --domains JSONL --vcf INDEXED_VCF --chrom CONTIG --output NEW_JSON\n\
       priva das-domains benign-merge --domains JSONL --min-benign INTEGER --input CHR1_JSON --input CHR2_JSON ...\n\
Supply all header chromosomes in literal order, including empty chromosomes.\n\
The indexed procedure is explicit; benign retains the original sequential route.\n\
Merge checks the source, VCF/index and domain-file bindings, then emits complete\n\
counts/removal flags in requested domain order. Existing shard files are preserved.";

fn error(message: &str) -> Box<dyn std::error::Error> {
    builder_error::category("ValueError", invalid(message))
}
fn encode(value: &Json) -> Result<String> {
    json_text::dumps(value, false).map_err(invalid)
}
fn domains(path: &Path) -> Result<Vec<String>> {
    let mut result = Vec::new();
    for line in BufReader::new(std::fs::File::open(path)?).lines() {
        let Json::Str(value) = json_text::parse(&line?, true).map_err(invalid)? else {
            return Err(error("each domain-path line must be a JSON string"));
        };
        result.push(value);
    }
    Ok(result)
}
fn source_identity() -> String {
    let mut sha = sha256::Sha256::new();
    for source in [
        include_bytes!("das_benign.rs").as_slice(),
        include_bytes!("das_benign_shards.rs"),
        include_bytes!("alpha_domain/vcf.rs"),
        include_bytes!("vcf_variant.rs"),
        include_bytes!("vcf_text.rs"),
        include_bytes!("tabix_index.rs"),
        include_bytes!("clinvar_aa_splice_builder.rs"),
        include_bytes!("clingen_map_builder.rs"),
        include_bytes!("python_text.rs"),
        include_bytes!("../Cargo.lock"),
    ] {
        sha.update(&(source.len() as u64).to_le_bytes());
        sha.update(source);
    }
    sha.finish_hex()
}

fn chromosome(domain_path: &Path, vcf: &Path, chrom: &str, output: &Path) -> Result<()> {
    let requested = domains(domain_path)?;
    let domain_sha = sha256::sha256_file(domain_path)?;
    let wanted: HashSet<_> = requested.iter().cloned().collect();
    let input = alpha_domain::input_identity(vcf, true)?;
    let context = Context::with_fields(
        vcf,
        &["DOMAINS", "Gene", "Consequence", "Feature_type", "Feature"],
    )?;
    let (counts, records) = das_benign::counts_for(vcf, &wanted, Some(chrom))?;
    alpha_domain::verify_input(&input)?;
    if sha256::sha256_file(domain_path)? != domain_sha {
        return Err(error("domain input changed while counting"));
    }
    let mut seen = HashSet::new();
    let sparse = requested
        .iter()
        .filter_map(|domain| {
            let count = counts.get(domain)?;
            seen.insert(domain).then(|| {
                Json::Array(vec![
                    Json::Str(domain.clone()),
                    Json::Int(i128::from(*count)),
                ])
            })
        })
        .collect();
    let values = vec![
        Json::Str(FORMAT.into()),
        Json::Str(source_identity()),
        Json::Str(input),
        Json::Str(domain_sha),
        Json::Int(requested.len() as i128),
        Json::Array(context.header.contigs.into_iter().map(Json::Str).collect()),
        Json::Str(chrom.into()),
        Json::Int(i128::from(records)),
        Json::Array(sparse),
    ];
    let body = encode(&Json::Object(
        KEYS.iter()
            .zip(values)
            .map(|(key, value)| ((*key).into(), value))
            .collect(),
    ))?;
    // These are new immutable count artifacts, not an overwrite publication.
    // Complete computation precedes creation; write/flush/sync failures surface.
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(output)?;
    file.write_all(body.as_bytes())?;
    file.write_all(b"\n")?;
    file.flush()?;
    file.sync_all()?;
    Ok(())
}

fn merge(domain_path: &Path, inputs: &[PathBuf], minimum: num_bigint::BigInt) -> Result<()> {
    let requested = domains(domain_path)?;
    let domain_sha = sha256::sha256_file(domain_path)?;
    let wanted: HashSet<_> = requested.iter().cloned().collect();
    let mut counts = HashMap::<String, u64>::new();
    let mut common_header: Option<Vec<(String, Json)>> = None;
    for (ordinal, path) in inputs.iter().enumerate() {
        let Json::Object(fields) =
            json_text::parse(&std::fs::read_to_string(path)?, true).map_err(invalid)?
        else {
            return Err(error("invalid benign count artifact"));
        };
        if !fields.iter().map(|(key, _)| key.as_str()).eq(KEYS) {
            return Err(error("invalid benign count artifact fields"));
        }
        if fields[0].1.as_str() != Some(FORMAT)
            || fields[1].1.as_str() != Some(source_identity().as_str())
            || fields[3].1.as_str() != Some(domain_sha.as_str())
            || !matches!(fields[4].1, Json::Int(n) if n == requested.len() as i128)
        {
            return Err(error("benign shard source or domain binding differs"));
        }
        let Json::Array(contigs) = &fields[5].1 else {
            return Err(error("invalid benign shard chromosome catalogue"));
        };
        if contigs.len() != inputs.len()
            || contigs.iter().any(|value| value.as_str().is_none())
            || contigs.get(ordinal).and_then(Json::as_str) != fields[6].1.as_str()
        {
            return Err(error(
                "benign shards must cover every header chromosome in order",
            ));
        }
        if let Some(old) = &common_header {
            if encode(&Json::Object(old.clone()))? != encode(&Json::Object(fields[..6].to_vec()))? {
                return Err(error("benign shard headers differ"));
            }
        } else {
            let input = fields[2]
                .1
                .as_str()
                .ok_or_else(|| error("invalid source VCF binding"))?;
            alpha_domain::verify_input(input)?;
            common_header = Some(fields[..6].to_vec());
        }
        if !matches!(fields[7].1, Json::Int(n) if u64::try_from(n).is_ok()) {
            return Err(error("invalid benign source-record count"));
        }
        let Json::Array(rows) = &fields[8].1 else {
            return Err(error("invalid benign count rows"));
        };
        let mut seen = HashSet::new();
        for row in rows {
            let Json::Array(pair) = row else {
                return Err(error("invalid benign count row"));
            };
            if pair.len() != 2 {
                return Err(error("invalid benign count row width"));
            }
            let domain = pair[0]
                .as_str()
                .ok_or_else(|| error("invalid benign count domain"))?;
            if !wanted.contains(domain) || !seen.insert(domain) {
                return Err(error("duplicate or unrequested benign count domain"));
            }
            let Json::Int(value) = pair[1] else {
                return Err(error("benign occurrence count must be an integer"));
            };
            let value = u64::try_from(value)
                .ok()
                .filter(|n| *n > 0)
                .ok_or_else(|| error("invalid sparse benign occurrence count"))?;
            let total = counts.entry(domain.into()).or_default();
            *total = total.checked_add(value).ok_or_else(|| {
                builder_error::category(
                    "OverflowError",
                    invalid("benign annotation count exceeds u64"),
                )
            })?;
        }
    }
    let header = common_header.ok_or_else(|| error("no benign count shards supplied"))?;
    alpha_domain::verify_input(header[2].1.as_str().unwrap())?;
    if sha256::sha256_file(domain_path)? != domain_sha {
        return Err(error("domain input changed during merge"));
    }
    let mut output = std::io::BufWriter::new(std::io::stdout().lock());
    for domain in requested {
        let count = i128::from(counts.get(&domain).copied().unwrap_or(0));
        let row = Json::Object(vec![
            ("domain".into(), Json::Str(domain)),
            ("n_benign_missense".into(), Json::Int(count)),
            (
                "removed".into(),
                Json::Bool(num_bigint::BigInt::from(count) >= minimum),
            ),
        ]);
        writeln!(output, "{}", encode(&row)?)?;
    }
    output.flush()?;
    Ok(())
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args.is_empty() || args.iter().any(|s| s == "--help" || s == "-h") {
        println!("{HELP}");
        return Ok(());
    }
    if args[0] == "benign-chromosome" {
        if args.len() != 9 {
            return Err(invalid(HELP));
        }
        let mut options = HashMap::new();
        for pair in args[1..].chunks_exact(2) {
            if !["--domains", "--vcf", "--chrom", "--output"].contains(&pair[0].as_str())
                || options.insert(pair[0].as_str(), pair[1].as_str()).is_some()
            {
                return Err(invalid(HELP));
            }
        }
        return builder_error::finish(chromosome(
            Path::new(options["--domains"]),
            Path::new(options["--vcf"]),
            options["--chrom"],
            Path::new(options["--output"]),
        ));
    }
    if args[0] != "benign-merge" || args.len() < 7 || args.len() % 2 == 0 {
        return Err(invalid(HELP));
    }
    let mut domain_path = None;
    let mut minimum = None;
    let mut inputs = Vec::new();
    for pair in args[1..].chunks_exact(2) {
        match pair[0].as_str() {
            "--domains" if domain_path.is_none() => domain_path = Some(PathBuf::from(&pair[1])),
            "--min-benign" if minimum.is_none() => minimum = Some(das_benign::minimum(&pair[1])?),
            "--input" => inputs.push(PathBuf::from(&pair[1])),
            _ => return Err(invalid(HELP)),
        }
    }
    let (Some(domain_path), Some(minimum)) = (domain_path, minimum) else {
        return Err(invalid(HELP));
    };
    builder_error::finish(merge(&domain_path, &inputs, minimum))
}
