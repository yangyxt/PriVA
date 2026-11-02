//! Opt-in AlphaMissense statistics and domain metadata builders. Python/VEP/DAS
//! are not called. Parquet raw shards retain all keys needed by ordered merge;
//! final distributions retain exact Float64 bits and every residue observation.
mod association;
mod key;
mod key_parquet;
mod metadata;
mod model;
mod parquet;
mod store;
mod vcf;

use crate::builder_error;
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use model::{Model, Node};
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

const PRODUCER: &str = "cc7b2d60665e8052e6362107e359cc133e31c398f243589ee63ad7e8f90341c4";
const EXTRACTOR: &str = "37502dedf58f44eb67237c462a9a0c19a84bc696764ecc870edabcaa60b5e92c";
const HELP: &str = "usage: priva alphamissense-domain build --vcf VCF --output PARQUET [--threads N]\n\
       priva alphamissense-domain chromosome --vcf INDEXED_VCF --chrom CONTIG --output PARQUET\n\
       priva alphamissense-domain merge --input CHR1_PARQUET --input CHR2_PARQUET ... --output PARQUET\n\
       priva alphamissense-domain check --input PARQUET [--vcf VCF]\n\
       priva domain-metadata {build|chromosome|merge|check} ... [--assembly hg19|hg38|GRCh37|GRCh38]\n\
Metadata build/chromosome also accept repeated --target-domain PATH. Merge must\n\
receive every chromosome in literal header order, including empty chromosomes.\n\
Metadata build --assembly enables the source-associated cache wrapper and its\n\
adjacent .source.json receipt; without --assembly it always extracts fresh data.\n\
build defaults to the Python sequential route; --threads >1 selects its indexed\n\
chromosome route. For bounded parallel construction use separate chromosome jobs\n\
and merge. Output is a typed Parquet/zstd intermediate; no caller is activated.";

pub(super) fn exception(
    kind: &'static str,
    message: impl Into<String>,
) -> Box<dyn std::error::Error> {
    builder_error::category(kind, invalid(message.into()))
}
pub(super) fn unsupported(message: &str) -> Box<dyn std::error::Error> {
    exception("UnsupportedLayout", message)
}
pub(super) fn source_error(error: Box<dyn std::error::Error>) -> Box<dyn std::error::Error> {
    let message = error.to_string();
    if message == "failed to fill whole buffer" {
        return exception("OSError", message);
    }
    for kind in [
        "TypeError",
        "IndexError",
        "ValueError",
        "KeyError",
        "AttributeError",
        "NotImplementedError",
        "OSError",
        "UnicodeDecodeError",
        "UnsupportedLayout",
    ] {
        if let Some(rest) = message
            .strip_prefix(kind)
            .and_then(|s| s.strip_prefix(": "))
        {
            return exception(kind, rest);
        }
    }
    error
}
fn json(value: Json) -> Result<String> {
    json_text::dumps(&value, false).map_err(invalid)
}
fn texts(values: impl IntoIterator<Item = String>) -> Result<String> {
    json(Json::Array(values.into_iter().map(Json::Str).collect()))
}
fn parse_texts(value: &str) -> Result<Vec<String>> {
    let Json::Array(values) = json_text::parse(value, true).map_err(invalid)? else {
        return Err(invalid("invalid ordered text metadata"));
    };
    values
        .into_iter()
        .map(|v| match v {
            Json::Str(s) => Ok(s),
            _ => Err(invalid("invalid ordered text metadata")),
        })
        .collect()
}
fn file_identity(path: &Path) -> Result<Json> {
    let metadata = std::fs::metadata(path)?;
    let mut fields = vec![
        (
            "path".into(),
            Json::Str(std::path::absolute(path)?.to_string_lossy().into()),
        ),
        (
            "target".into(),
            Json::Str(std::fs::canonicalize(path)?.to_string_lossy().into()),
        ),
        ("bytes".into(), Json::Int(metadata.len() as i128)),
        (
            "mtime_ns".into(),
            Json::Str(
                metadata
                    .modified()?
                    .duration_since(std::time::UNIX_EPOCH)?
                    .as_nanos()
                    .to_string(),
            ),
        ),
    ];
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        fields.push(("inode".into(), Json::Str(metadata.ino().to_string())));
    }
    Ok(Json::Object(fields))
}
fn input_identity(path: &Path, indexes: bool) -> Result<String> {
    let mut fields = vec![
        ("vcf".into(), file_identity(path)?),
        ("indexes".into(), Json::Bool(indexes)),
    ];
    for suffix in if indexes {
        &[".tbi", ".csi"][..]
    } else {
        &[][..]
    } {
        let index = PathBuf::from(format!("{}{suffix}", path.display()));
        if index.exists() {
            fields.push((
                (*suffix).into(),
                Json::Object(vec![
                    ("identity".into(), file_identity(&index)?),
                    (
                        "sha256".into(),
                        Json::Str(crate::sha256::sha256_file(&index)?),
                    ),
                ]),
            ));
        }
    }
    json(Json::Object(fields))
}
fn verify_input(encoded: &str) -> Result<()> {
    let value = json_text::parse(encoded, true).map_err(invalid)?;
    let path = value
        .get("vcf")
        .and_then(|v| v.get("path"))
        .and_then(Json::as_str)
        .ok_or_else(|| invalid("input identity has no VCF path"))?;
    let indexes = matches!(value.get("indexes"), Some(Json::Bool(true)));
    if input_identity(Path::new(path), indexes)? != encoded {
        return Err(exception("ValueError", "source VCF/index identity changed"));
    }
    Ok(())
}
fn protect(input: &Path, output: &Path) -> Result<()> {
    if output.exists() {
        crate::rust_cache::ensure_distinct_paths(input, output)?;
    }
    Ok(())
}
fn builder_identity() -> Result<String> {
    // Content identity of this builder and its semantic primitives, embedded
    // by rustc. Cache reuse does not depend on executable path, inode, debug
    // symbols or build-directory names. The lockfile pins the Parquet codec.
    static DIGEST: std::sync::OnceLock<String> = std::sync::OnceLock::new();
    Ok(DIGEST
        .get_or_init(|| {
            let mut sha = crate::sha256::Sha256::new();
            for source in [
                include_bytes!("mod.rs").as_slice(),
                include_bytes!("association.rs"),
                include_bytes!("key.rs"),
                include_bytes!("key_parquet.rs"),
                include_bytes!("model.rs"),
                include_bytes!("metadata.rs"),
                include_bytes!("parquet.rs"),
                include_bytes!("store.rs"),
                include_bytes!("vcf.rs"),
                include_bytes!("../vcf_variant.rs"),
                include_bytes!("../vcf_text.rs"),
                include_bytes!("../tabix_index.rs"),
                include_bytes!("../clinvar_aa_splice_builder.rs"),
                include_bytes!("../clingen_map_builder.rs"),
                include_bytes!("../python_text.rs"),
                include_bytes!("../json_text.rs"),
                include_bytes!("../builder_error.rs"),
                include_bytes!("../../Cargo.lock"),
            ] {
                sha.update(&(source.len() as u64).to_le_bytes());
                sha.update(source);
            }
            sha.finish_hex()
        })
        .clone())
}

struct Options {
    command: String,
    metadata: bool,
    vcf: Option<PathBuf>,
    output: Option<PathBuf>,
    inputs: Vec<PathBuf>,
    chrom: Option<String>,
    threads: i64,
    assembly: Option<String>,
    targets: Option<HashSet<String>>,
}
fn common(options: &Options, context: &vcf::Context, input: &str) -> Result<Vec<(String, String)>> {
    let ordinal = options
        .chrom
        .as_ref()
        .map(|chrom| {
            context
                .header
                .contigs
                .iter()
                .position(|s| s == chrom)
                .ok_or_else(|| exception("ValueError", "invalid contig"))
        })
        .transpose()?;
    let mut values = vec![
        (
            "priva.kind".into(),
            if options.metadata {
                "metadata"
            } else {
                "statistics"
            }
            .into(),
        ),
        (
            "priva.stage".into(),
            if options.command == "chromosome" {
                "raw"
            } else {
                "final"
            }
            .into(),
        ),
        ("priva.input".into(), input.into()),
        ("priva.producer".into(), PRODUCER.into()),
        ("priva.extractor".into(), EXTRACTOR.into()),
        ("priva.builder_source_sha256".into(), builder_identity()?),
        (
            "priva.contigs".into(),
            texts(context.header.contigs.clone())?,
        ),
        (
            "priva.ordinal".into(),
            ordinal
                .map(|n| n.to_string())
                .unwrap_or_else(|| "all".into()),
        ),
    ];
    if let Some(assembly) = &options.assembly {
        values.push(("priva.assembly".into(), assembly.clone()));
    }
    if let Some(targets) = &options.targets {
        let mut targets: Vec<_> = targets.iter().cloned().collect();
        targets.sort();
        values.push(("priva.targets".into(), texts(targets)?));
    }
    Ok(values)
}
fn build(options: &Options) -> Result<()> {
    let path = options.vcf.as_ref().unwrap();
    let output = options.output.as_ref().unwrap();
    let associated = options.metadata && options.command == "build" && options.assembly.is_some();
    if associated {
        association::source(path)?;
    }
    protect(path, output)?;
    for suffix in [".tbi", ".csi"] {
        let index = PathBuf::from(format!("{}{suffix}", path.display()));
        if index.exists() {
            protect(&index, output)?;
        }
    }
    let input = input_identity(path, !options.metadata)?;
    // The Python metadata helper consults source association before invoking
    // its extractor. A matched but unreadable completed cache propagates its
    // read error; an unassociated file is replaced only after fresh success.
    if associated && association::matches(output, &input, options) {
        {
            let old = parquet::metadata(output)?;
            // The source wrapper associates VCF/assembly/implementation, not
            // extractor closure arguments. A warm cache retains its original
            // target selection; without --assembly extraction is always fresh.
            if !(old.get("priva.kind").map(String::as_str) == Some("metadata")
                && old.get("priva.stage").map(String::as_str) == Some("final")
                && old.get("priva.input") == Some(&input)
                && old.get("priva.producer").map(String::as_str) == Some(PRODUCER)
                && old.get("priva.extractor").map(String::as_str) == Some(EXTRACTOR)
                && old.get("priva.builder_source_sha256") == Some(&builder_identity()?)
                && old.get("priva.assembly") == options.assembly.as_ref())
            {
                return Err(exception(
                    "ValueError",
                    "associated metadata payload profile differs",
                ));
            }
            {
                let before = json(file_identity(output)?)?;
                let value = metadata::read(output)?;
                verify_input(&input)?;
                if json(file_identity(output)?)? != before
                    || !association::matches(output, &input, options)
                {
                    return Err(exception(
                        "ValueError",
                        "metadata output changed during loading",
                    ));
                }
                // An empty Python dict has no first prot_length and therefore
                // always regenerates. Zero is a present, reusable length.
                if !value.entries.is_empty() {
                    println!("{{\"keys\":{},\"reused\":true}}", value.entries.len());
                    return Ok(());
                }
            }
        }
    }
    let context = vcf::Context::new(path, options.metadata)?;
    let regions = if let Some(chrom) = &options.chrom {
        vec![Some(chrom.clone())]
    } else if options.threads > 1 {
        if context.header.contigs.is_empty() {
            return Err(exception(
                "ValueError",
                "Number of processes must be at least 1",
            ));
        }
        context.header.contigs.iter().cloned().map(Some).collect()
    } else {
        vec![None]
    };
    let mut footer = common(options, &context, &input)?;
    let rows;
    if options.metadata {
        let mut collected = metadata::Metadata::default();
        for region in regions {
            let mut part = metadata::Metadata::default();
            vcf::visit(path, &context, region.as_deref(), |record, definitions| {
                metadata::record(
                    &mut part,
                    record,
                    definitions,
                    &context,
                    options.targets.as_ref(),
                )
            })?;
            collected.merge(part)?;
        }
        verify_input(&input)?;
        footer.push((
            "priva.keys".into(),
            texts(collected.entries.iter().map(|(k, _)| k.clone()))?,
        ));
        if associated {
            if let Some(parent) = output.parent().filter(|p| !p.as_os_str().is_empty()) {
                std::fs::create_dir_all(parent)?;
            }
        }
        rows = parquet::publish(output, footer, |writer| metadata::write(writer, &collected))?;
        if associated {
            association::publish(output, &input, options)?;
        }
    } else {
        let mut collected = Model::default();
        for region in regions {
            let mut part = Model::default();
            vcf::visit(path, &context, region.as_deref(), |record, definitions| {
                vcf::record(&mut part, record, definitions, &context)
            })?;
            collected.merge(part)?;
        }
        verify_input(&input)?;
        footer.push((
            "priva.keys".into(),
            texts(collected.genes.members.iter().map(|(k, _)| k.clone()))?,
        ));
        rows = parquet::publish(output, footer, |writer| {
            store::write_model(writer, &collected, options.command == "chromosome")
        })?;
    }
    println!("{{\"rows\":{rows},\"format\":\"parquet-zstd\",\"schema\":2}}");
    Ok(())
}

fn merge(options: &Options) -> Result<()> {
    let output = options.output.as_ref().unwrap();
    let mut headers = Vec::new();
    for path in &options.inputs {
        protect(path, output)?;
        headers.push(parquet::metadata(path)?);
    }
    let first = &headers[0];
    let required = |key: &str| {
        first
            .get(key)
            .ok_or_else(|| exception("ValueError", format!("raw checkpoint missing {key}")))
    };
    let contigs = parse_texts(required("priva.contigs")?)?;
    if contigs.len() != headers.len() {
        return Err(exception(
            "ValueError",
            "merge requires every header chromosome",
        ));
    }
    let kind = if options.metadata {
        "metadata"
    } else {
        "statistics"
    };
    let input = required("priva.input")?.clone();
    if options
        .assembly
        .as_ref()
        .is_some_and(|value| first.get("priva.assembly") != Some(value))
    {
        return Err(exception(
            "ValueError",
            "selected assembly differs from checkpoints",
        ));
    }
    if required("priva.producer")? != PRODUCER
        || required("priva.extractor")? != EXTRACTOR
        || required("priva.builder_source_sha256")? != &builder_identity()?
    {
        return Err(exception(
            "ValueError",
            "checkpoint producer/implementation differs",
        ));
    }
    verify_input(&input)?;
    let mut global_keys = Vec::new();
    let mut first_at = HashMap::new();
    let mut last_at = HashMap::new();
    let mut inventories = Vec::new();
    for (ordinal, header) in headers.iter().enumerate() {
        if header.get("priva.stage").map(String::as_str) != Some("raw")
            || header.get("priva.kind").map(String::as_str) != Some(kind)
            || header.get("priva.ordinal") != Some(&ordinal.to_string())
        {
            return Err(exception(
                "ValueError",
                "wrong checkpoint stage/kind/chromosome order",
            ));
        }
        for key in [
            "priva.input",
            "priva.producer",
            "priva.extractor",
            "priva.builder_source_sha256",
            "priva.contigs",
            "priva.assembly",
            "priva.targets",
        ] {
            if header.get(key) != first.get(key) {
                return Err(exception("ValueError", "checkpoint source/profile differs"));
            }
        }
        let keys = parse_texts(
            header
                .get("priva.keys")
                .ok_or_else(|| invalid("checkpoint keys missing"))?,
        )?;
        let mut seen = HashSet::new();
        for key in &keys {
            if !seen.insert(key) {
                return Err(exception(
                    "ValueError",
                    "duplicate checkpoint inventory key",
                ));
            }
            if !first_at.contains_key(key) {
                first_at.insert(key.clone(), global_keys.len());
                global_keys.push(key.clone());
            }
            last_at.insert(key.clone(), ordinal);
        }
        inventories.push(keys);
    }
    let mut footer: Vec<_> = first
        .iter()
        .filter(|(key, _)| {
            !["priva.schema", "priva.ordinal", "priva.stage", "priva.keys"].contains(&key.as_str())
        })
        .map(|(key, value)| (key.clone(), value.clone()))
        .collect();
    footer.sort();
    footer.extend([
        ("priva.stage".into(), "final".into()),
        ("priva.ordinal".into(), "all".into()),
        ("priva.keys".into(), texts(global_keys.clone())?),
    ]);
    let rows = if options.metadata {
        let mut merged = metadata::Metadata::default();
        for (ordinal, path) in options.inputs.iter().enumerate() {
            let part = metadata::read(path)?;
            if part
                .entries
                .iter()
                .map(|(k, _)| k)
                .ne(inventories[ordinal].iter())
            {
                return Err(exception(
                    "ValueError",
                    "metadata checkpoint inventory differs",
                ));
            }
            merged.merge(part)?;
        }
        parquet::publish(output, footer, |writer| metadata::write(writer, &merged))?
    } else {
        let mut pending: HashMap<String, Node> = HashMap::new();
        parquet::publish(output, footer, |writer| {
            for (ordinal, path) in options.inputs.iter().enumerate() {
                let mut observed = Vec::new();
                store::scan_genes(path, |gene, gene_order, node| {
                    if gene_order != observed.len() as u64
                        || inventories[ordinal].get(observed.len()) != Some(&gene)
                    {
                        return Err(exception("ValueError", "raw gene inventory/order differs"));
                    }
                    observed.push(gene.clone());
                    if let Some(old) = pending.get_mut(&gene) {
                        old.merge_gene(node)?;
                    } else {
                        pending.insert(gene.clone(), node);
                    }
                    if last_at[&gene] == ordinal {
                        let node = pending.remove(&gene).unwrap();
                        store::write_gene(writer, &gene, first_at[&gene] as u64, &node, false)?;
                    }
                    Ok(())
                })?;
                if observed != inventories[ordinal] {
                    return Err(exception("ValueError", "missing raw checkpoint gene"));
                }
            }
            if !pending.is_empty() {
                return Err(exception("ValueError", "unmerged gene contribution"));
            }
            verify_input(&input)
        })?
    };
    if options.metadata && options.assembly.is_some() {
        verify_input(&input)?;
        association::publish(output, &input, options)?;
    }
    println!(
        "{{\"rows\":{rows},\"keys\":{},\"chromosomes\":{},\"format\":\"parquet-zstd\"}}",
        global_keys.len(),
        contigs.len()
    );
    Ok(())
}

fn check(options: &Options) -> Result<()> {
    let path = &options.inputs[0];
    let header = parquet::metadata(path)?;
    let expected = header
        .get("priva.input")
        .ok_or_else(|| invalid("input identity absent"))?;
    if let Some(vcf) = &options.vcf {
        if input_identity(vcf, !options.metadata)? != *expected {
            return Err(exception(
                "ValueError",
                "selected VCF/index association differs",
            ));
        }
    }
    if let Some(assembly) = &options.assembly {
        if header.get("priva.assembly") != Some(assembly) {
            return Err(exception(
                "ValueError",
                "selected assembly association differs",
            ));
        }
        if options.metadata {
            verify_input(expected)?;
            if !association::matches(path, expected, options) {
                return Err(exception(
                    "ValueError",
                    "metadata has no matching completed source association",
                ));
            }
        }
    }
    if header.get("priva.producer").map(String::as_str) != Some(PRODUCER)
        || header.get("priva.extractor").map(String::as_str) != Some(EXTRACTOR)
        || header.get("priva.builder_source_sha256") != Some(&builder_identity()?)
    {
        return Err(exception(
            "ValueError",
            "selected producer association differs",
        ));
    }
    let keys = parse_texts(
        header
            .get("priva.keys")
            .ok_or_else(|| invalid("key inventory absent"))?,
    )?;
    if options.metadata {
        let actual = metadata::read(path)?;
        if actual.entries.iter().map(|(key, _)| key).ne(keys.iter()) {
            return Err(exception("ValueError", "metadata key inventory differs"));
        }
    } else {
        if header.get("priva.kind").map(String::as_str) != Some("statistics") {
            return Err(exception("ValueError", "not a statistics intermediate"));
        }
        let mut observed = HashMap::new();
        store::scan_genes(path, |gene, ordinal, _| {
            if observed.insert(ordinal, gene).is_some() {
                return Err(exception("ValueError", "duplicate gene ordinal"));
            }
            Ok(())
        })?;
        if observed.len() != keys.len()
            || keys
                .iter()
                .enumerate()
                .any(|(i, key)| observed.get(&(i as u64)) != Some(key))
        {
            return Err(exception("ValueError", "statistics gene inventory differs"));
        }
    }
    println!("{{\"passed\":true,\"keys\":{}}}", keys.len());
    Ok(())
}

pub(crate) fn dispatch(args: &[String], metadata: bool) -> Result<()> {
    if args.iter().any(|s| s == "--help" || s == "-h") {
        println!("{HELP}");
        return Ok(());
    }
    let command = args
        .first()
        .filter(|s| ["build", "chromosome", "merge", "check"].contains(&s.as_str()))
        .ok_or_else(|| invalid(HELP))?
        .clone();
    let mut fields = HashMap::new();
    let mut inputs = Vec::new();
    let mut targets = None;
    let mut at = 1;
    while at < args.len() {
        let name = args[at].as_str();
        let value = args.get(at + 1).ok_or_else(|| invalid(HELP))?;
        let admitted = match command.as_str() {
            "build" => ["--vcf", "--output", "--threads"].contains(&name),
            "chromosome" => ["--vcf", "--output", "--chrom"].contains(&name),
            "merge" => ["--input", "--output"].contains(&name),
            _ => ["--input", "--vcf"].contains(&name),
        } || (metadata && name == "--assembly")
            || (metadata
                && ["build", "chromosome"].contains(&command.as_str())
                && name == "--target-domain");
        if !admitted {
            return Err(invalid(HELP));
        }
        match name {
            "--input" => inputs.push(PathBuf::from(value)),
            "--target-domain" => {
                targets
                    .get_or_insert_with(HashSet::new)
                    .insert(value.clone());
            }
            _ => {
                fields.insert(name, value.clone());
            }
        }
        at += 2;
    }
    let need = |name| fields.get(name).cloned().ok_or_else(|| invalid(HELP));
    let output = if command == "check" {
        None
    } else {
        Some(need("--output")?.into())
    };
    let vcf = if ["build", "chromosome"].contains(&command.as_str()) {
        Some(need("--vcf")?.into())
    } else {
        fields.get("--vcf").map(PathBuf::from)
    };
    let chrom = if command == "chromosome" {
        Some(need("--chrom")?)
    } else {
        None
    };
    if (command == "merge" && inputs.is_empty()) || (command == "check" && inputs.len() != 1) {
        return Err(invalid(HELP));
    }
    // The source uses Python int() and only distinguishes <=1 from >1 before
    // limiting workers to the number of contigs. Preserve Unicode digits,
    // whitespace, underscores and arbitrarily large admitted integer values.
    let threads = if let Some(text) = fields.get("--threads") {
        if !crate::python_text::python_int_is_valid(text) {
            return Err(invalid("--threads requires an integer"));
        }
        let digits: String = text
            .chars()
            .filter_map(crate::python_text::decimal_digit_value)
            .map(|d| char::from(b'0' + d as u8))
            .collect();
        let digits = digits.trim_start_matches('0');
        if text.trim_start().starts_with('-') || digits.is_empty() || digits == "1" {
            1
        } else {
            2
        }
    } else {
        1
    };
    let assembly = fields
        .get("--assembly")
        .map(|s| match s.as_str() {
            "GRCh37" | "hg19" => Ok("hg19".to_owned()),
            "GRCh38" | "hg38" => Ok("hg38".to_owned()),
            _ => Err(invalid("unsupported assembly")),
        })
        .transpose()?;
    let options = Options {
        command: command.clone(),
        metadata,
        vcf,
        output,
        inputs,
        chrom,
        threads,
        assembly,
        targets,
    };
    builder_error::finish(match command.as_str() {
        "merge" => merge(&options),
        "check" => check(&options),
        _ => build(&options),
    })
}
