//! Native construction/serialization history for fresh Alpha-to-DAS builds.
//! Raw chromosomes preserve first exon insertion; final files preserve the
//! explicit CPython-seed set iteration after the original transport/merge.
mod model;
mod store;

use super::{
    PRODUCER, Result, builder_accepted, exception, input_identity, json, parquet, parse_texts,
    protect, texts, vcf, verify_input,
};
use crate::builder_error;
use crate::cli;
use crate::json_text::Json;
use model::{Entry, Model, Node};
use std::collections::{HashMap, HashSet};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

pub(super) const SCHEMA: &str = "priva-alpha-exon-history-v1";
const HELP: &str = "usage: priva alpha-exon-history build --vcf VCF --hash-seed UINT32 --output PARQUET\n\
       priva alpha-exon-history chromosome --vcf INDEXED_VCF --chrom CONTIG --output PARQUET\n\
       priva alpha-exon-history merge --input CHR1_PARQUET --input CHR2_PARQUET ... --hash-seed UINT32 --output PARQUET\n\
       priva alpha-exon-history check --input PARQUET [--builder-sha256 SHA256]\n\
       priva alpha-exon-history identity --input PARQUET\n\
       priva alpha-exon-history export --input FINAL_PARQUET --output HISTORY\n\
Check prints the recorded and current builder_source_sha256 and by default\n\
requires the current builder; --builder-sha256 instead accepts history recorded\n\
by that pinned builder, with every other check unchanged. Identity prints the\n\
recorded value without qualifying the file. Merge and export require the\n\
current builder; rebuild only on explicit request.\n\
Chromosome jobs retain seed-independent first insertion order. Ordered merge\n\
requires all header chromosomes and reproduces CPython 3.12 worker transport\n\
and set.update under the explicit seed. Build uses original sequential order.\n\
Export emits native serialized-set history for das-domains transcripts; it\n\
does not call Python, write pickle, or activate an installer caller.";

fn invalid(message: &str) -> Box<dyn std::error::Error> {
    exception("ValueError", message)
}
pub(crate) fn identity() -> String {
    static DIGEST: std::sync::OnceLock<String> = std::sync::OnceLock::new();
    DIGEST
        .get_or_init(|| {
            let mut sha = crate::sha256::Sha256::new();
            for source in [
                include_bytes!("mod.rs").as_slice(),
                include_bytes!("model.rs"),
                include_bytes!("store.rs"),
                include_bytes!("../mod.rs"),
                include_bytes!("../vcf.rs"),
                include_bytes!("../model.rs"),
                include_bytes!("../key.rs"),
                include_bytes!("../key_parquet.rs"),
                include_bytes!("../parquet.rs"),
                include_bytes!("../../python_set_order.rs"),
                include_bytes!("../../vcf_variant.rs"),
                include_bytes!("../../vcf_text.rs"),
                include_bytes!("../../tabix_index.rs"),
                include_bytes!("../../clinvar_aa_splice_builder.rs"),
                include_bytes!("../../clingen_map_builder.rs"),
                include_bytes!("../../python_text.rs"),
                include_bytes!("../../json_text.rs"),
                include_bytes!("../../builder_error.rs"),
                include_bytes!("../../../Cargo.lock"),
            ] {
                sha.update(&(source.len() as u64).to_le_bytes());
                sha.update(source);
            }
            sha.finish_hex()
        })
        .clone()
}
fn checked(path: &Path, builder: Option<&str>) -> Result<HashMap<String, String>> {
    let metadata = parquet::metadata(path)?;
    for (key, wanted) in [
        ("priva.kind", "exon-history"),
        ("priva.history.schema", SCHEMA),
        ("priva.producer", PRODUCER),
    ] {
        if metadata.get(key).map(String::as_str) != Some(wanted) {
            return Err(invalid("exon history source/profile differs"));
        }
    }
    if !builder_accepted(
        metadata.get("priva.builder_source_sha256"),
        builder,
        &identity(),
    ) {
        return Err(invalid("exon history source/profile differs"));
    }
    let stage = metadata.get("priva.stage").map(String::as_str);
    let profile = metadata.get("priva.history.profile").map(String::as_str);
    match (stage, profile) {
        (Some("raw"), Some("first-insertion-v1"))
            if !metadata.contains_key("priva.history.seed") => {}
        (Some("final"), Some("cpython312-serial-v1" | "cpython312-pool-v1")) => {
            let _seed: u32 = metadata
                .get("priva.history.seed")
                .and_then(|s| s.parse().ok())
                .ok_or_else(|| invalid("final exon history seed absent/invalid"))?;
        }
        _ => return Err(invalid("invalid exon history stage/seed profile")),
    }
    Ok(metadata)
}
fn keys(model: &Model) -> Result<String> {
    texts(model.genes.members.iter().map(|(key, _)| key.clone()))
}
fn build(vcf_path: &Path, chrom: Option<&str>, seed: Option<u32>, output: &Path) -> Result<()> {
    protect(vcf_path, output)?;
    for suffix in [".tbi", ".csi"] {
        let index = PathBuf::from(format!("{}{suffix}", vcf_path.display()));
        if index.exists() {
            protect(&index, output)?;
        }
    }
    let input = input_identity(vcf_path, true)?;
    let context = vcf::Context::new(vcf_path, false)?;
    let ordinal = chrom
        .map(|name| {
            context
                .header
                .contigs
                .iter()
                .position(|c| c == name)
                .ok_or_else(|| invalid("invalid contig"))
        })
        .transpose()?;
    let mut model = Model::default();
    vcf::visit(vcf_path, &context, chrom, |record, definitions| {
        vcf::observations(record, definitions, &context, |item| model.observe(item))
    })?;
    if let Some(seed) = seed {
        model.prepare(seed, false);
    }
    let mut metadata = vec![
        ("priva.kind".into(), "exon-history".into()),
        ("priva.history.schema".into(), SCHEMA.into()),
        (
            "priva.stage".into(),
            if chrom.is_some() { "raw" } else { "final" }.into(),
        ),
        (
            "priva.history.profile".into(),
            if chrom.is_some() {
                "first-insertion-v1"
            } else {
                "cpython312-serial-v1"
            }
            .into(),
        ),
        ("priva.input".into(), input.clone()),
        ("priva.producer".into(), PRODUCER.into()),
        ("priva.builder_source_sha256".into(), identity()),
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
        ("priva.keys".into(), keys(&model)?),
    ];
    if let Some(seed) = seed {
        metadata.push(("priva.history.seed".into(), seed.to_string()));
    }
    let rows = parquet::publish(output, metadata, |writer| {
        store::write_model(writer, &model)?;
        verify_input(&input)
    })?;
    println!(
        "{{\"rows\":{rows},\"genes\":{},\"history_schema\":1}}",
        model.genes.members.len()
    );
    Ok(())
}
fn merge(inputs: &[PathBuf], seed: u32, output: &Path) -> Result<()> {
    let mut headers = vec![];
    for path in inputs {
        protect(path, output)?;
        headers.push(checked(path, None)?);
    }
    let first = &headers[0];
    let contigs = parse_texts(
        first
            .get("priva.contigs")
            .ok_or_else(|| invalid("history contigs absent"))?,
    )?;
    let input = first
        .get("priva.input")
        .ok_or_else(|| invalid("history input identity absent"))?;
    if contigs.len() != inputs.len() {
        return Err(invalid("history merge requires every header chromosome"));
    }
    verify_input(input)?;
    let mut global = vec![];
    let mut positions = HashMap::new();
    let mut last = HashMap::new();
    for (ordinal, header) in headers.iter().enumerate() {
        if header.get("priva.stage").map(String::as_str) != Some("raw")
            || header.get("priva.ordinal") != Some(&ordinal.to_string())
        {
            return Err(invalid("history merge chromosome order/stage differs"));
        }
        for key in [
            "priva.input",
            "priva.contigs",
            "priva.producer",
            "priva.builder_source_sha256",
            "priva.history.schema",
            "priva.history.profile",
        ] {
            if header.get(key) != first.get(key) {
                return Err(invalid("history checkpoint source/profile differs"));
            }
        }
        let inventory = parse_texts(
            header
                .get("priva.keys")
                .ok_or_else(|| invalid("history keys absent"))?,
        )?;
        let mut seen = HashSet::new();
        for key in inventory {
            if !seen.insert(key.clone()) {
                return Err(invalid("duplicate history inventory gene"));
            }
            if !positions.contains_key(&key) {
                positions.insert(key.clone(), global.len());
                global.push(key.clone());
            }
            last.insert(key, ordinal);
        }
    }
    let mut metadata: Vec<_> = first
        .iter()
        .filter(|(key, _)| {
            ![
                "priva.schema",
                "priva.stage",
                "priva.ordinal",
                "priva.keys",
                "priva.history.profile",
            ]
            .contains(&key.as_str())
        })
        .map(|(key, value)| (key.clone(), value.clone()))
        .collect();
    metadata.sort();
    metadata.extend([
        ("priva.stage".into(), "final".into()),
        ("priva.ordinal".into(), "all".into()),
        ("priva.keys".into(), texts(global.clone())?),
        ("priva.history.seed".into(), seed.to_string()),
        ("priva.history.profile".into(), "cpython312-pool-v1".into()),
    ]);
    let mut pending: HashMap<String, Node> = HashMap::new();
    let mut emitted = 0;
    let rows = parquet::publish(output, metadata, |writer| {
        for (ordinal, path) in inputs.iter().enumerate() {
            store::scan_genes(path, |gene, _, mut node| {
                node.prepare(seed, true);
                if let Some(old) = pending.get_mut(&gene) {
                    old.merge_gene(node)?;
                } else {
                    pending.insert(gene, node);
                }
                Ok(())
            })?;
            // Hold a later completed gene if an earlier gene still has a
            // contribution pending. Dictionary order is the first encounter,
            // never worker completion time or last-contribution order.
            while let Some(gene) = global.get(emitted).filter(|gene| last[*gene] <= ordinal) {
                let node = pending
                    .remove(gene)
                    .ok_or_else(|| invalid("history gene contribution absent"))?;
                store::write_gene(writer, gene, emitted as u64, &node)?;
                emitted += 1;
            }
        }
        if emitted != global.len() || !pending.is_empty() {
            return Err(invalid("history merge coverage differs"));
        }
        verify_input(input)
    })?;
    println!(
        "{{\"rows\":{rows},\"genes\":{},\"chromosomes\":{},\"history_schema\":1}}",
        global.len(),
        contigs.len()
    );
    Ok(())
}
fn visit_transcripts(
    node: &Node,
    path: &mut Vec<String>,
    visit: &mut impl FnMut(&[String], &model::Transcript) -> Result<()>,
) -> Result<()> {
    for (key, entry) in &node.members {
        path.push(key.clone());
        match entry {
            Entry::Branch(child) => visit_transcripts(child, path, visit)?,
            Entry::Transcript(value) => visit(path, value)?,
            Entry::Marker => {}
        }
        path.pop();
    }
    Ok(())
}
fn check(path: &Path, builder: Option<&str>) -> Result<()> {
    let metadata = checked(path, builder)?;
    let mut genes = 0;
    let mut sets = 0;
    let mut members = 0;
    store::scan_genes(path, |gene, _, node| {
        genes += 1;
        visit_transcripts(&node, &mut vec![gene], &mut |_, value| {
            sets += 1;
            members += value.exons.len();
            Ok(())
        })
    })?;
    println!(
        "{{\"passed\":true,\"genes\":{genes},\"sets\":{sets},\"members\":{members},{}}}",
        super::identity_fields(&metadata["priva.builder_source_sha256"], &identity())
    );
    Ok(())
}
fn hex(value: &str) -> String {
    if value.is_empty() {
        return "-".into();
    }
    value.bytes().map(|byte| format!("{byte:02x}")).collect()
}
fn export(input: &Path, output: &Path) -> Result<()> {
    protect(input, output)?;
    let metadata = checked(input, None)?;
    if metadata.get("priva.stage").map(String::as_str) != Some("final") {
        return Err(invalid("only final history can be exported"));
    }
    let before = super::file_identity(input)?;
    let mut sets = 0;
    let mut fields: Vec<_> = metadata
        .into_iter()
        .filter(|(key, _)| key != "priva.keys")
        .map(|(key, value)| (key, Json::Str(value)))
        .collect();
    fields.sort_by(|a, b| a.0.cmp(&b.0));
    let header = json(Json::Object(fields))?;
    parquet::atomic_file(output, |file, _| {
        let mut stream = BufWriter::new(file);
        writeln!(stream, "# {header}")?;
        store::scan_genes(input, |gene, _, node| {
            visit_transcripts(&node, &mut vec![gene], &mut |path, value| {
                write!(stream, "{}", path.len() + 1)?;
                for part in path {
                    write!(stream, " {}", hex(part))?;
                }
                write!(stream, " 65786f6e73 {}", value.exons.len())?;
                for exon in &value.exons {
                    write!(stream, " {}", hex(exon))?;
                }
                writeln!(stream)?;
                sets += 1;
                Ok(())
            })
        })?;
        if json(super::file_identity(input)?)? != json(before)? {
            return Err(invalid("history changed during export"));
        }
        stream.flush()?;
        stream.get_ref().sync_all()?;
        Ok(())
    })?;
    println!("{{\"sets\":{sets},\"history_schema\":1}}");
    Ok(())
}
pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args.iter().any(|s| s == "--help" || s == "-h") {
        println!("{HELP}");
        return Ok(());
    }
    let command = args
        .first()
        .map(String::as_str)
        .filter(|s| {
            [
                "build",
                "chromosome",
                "merge",
                "check",
                "identity",
                "export",
            ]
            .contains(s)
        })
        .ok_or_else(|| cli::invalid(HELP))?;
    let mut fields = HashMap::new();
    let mut inputs = vec![];
    let mut at = 1;
    while at < args.len() {
        let name = args[at].as_str();
        let value = args.get(at + 1).ok_or_else(|| cli::invalid(HELP))?;
        let admitted = match command {
            "build" => ["--vcf", "--hash-seed", "--output"].contains(&name),
            "chromosome" => ["--vcf", "--chrom", "--output"].contains(&name),
            "merge" => ["--input", "--hash-seed", "--output"].contains(&name),
            "export" => ["--input", "--output"].contains(&name),
            "check" => ["--input", "--builder-sha256"].contains(&name),
            _ => name == "--input",
        };
        if !admitted {
            return Err(cli::invalid(HELP));
        }
        if name == "--input" {
            inputs.push(PathBuf::from(value));
        } else if fields.insert(name, value.as_str()).is_some() {
            return Err(cli::invalid("duplicate history option"));
        }
        at += 2;
    }
    let need = |name| fields.get(name).copied().ok_or_else(|| cli::invalid(HELP));
    let output = if !["check", "identity"].contains(&command) {
        Some(Path::new(need("--output")?))
    } else {
        None
    };
    let seed = if ["build", "merge"].contains(&command) {
        Some(
            need("--hash-seed")?
                .parse::<u32>()
                .map_err(|_| cli::invalid("hash seed must be UINT32"))?,
        )
    } else {
        None
    };
    if (command == "merge" && inputs.is_empty())
        || (["check", "identity", "export"].contains(&command) && inputs.len() != 1)
    {
        return Err(cli::invalid(HELP));
    }
    let builder = fields
        .get("--builder-sha256")
        .copied()
        .map(super::builder_pin)
        .transpose()?;
    let result = match command {
        "build" => build(Path::new(need("--vcf")?), None, seed, output.unwrap()),
        "chromosome" => build(
            Path::new(need("--vcf")?),
            Some(need("--chrom")?),
            None,
            output.unwrap(),
        ),
        "merge" => merge(&inputs, seed.unwrap(), output.unwrap()),
        "check" => check(&inputs[0], builder.as_deref()),
        "identity" => super::recorded_identity(&inputs[0]),
        _ => export(&inputs[0], output.unwrap()),
    };
    builder_error::finish(result.map_err(super::source_error))
}

#[cfg(test)]
mod tests {
    use super::{PRODUCER, SCHEMA, check, checked, dispatch, identity, parquet};

    #[test]
    fn history_check_accepts_an_earlier_builder_only_when_pinned() {
        let path = std::env::temp_dir().join(format!(
            "priva_history_builder_pin_{}.parquet",
            std::process::id()
        ));
        let older = "b".repeat(64);
        let footer = [
            ("priva.kind", "exon-history"),
            ("priva.history.schema", SCHEMA),
            ("priva.stage", "final"),
            ("priva.history.profile", "cpython312-serial-v1"),
            ("priva.history.seed", "0"),
            ("priva.input", "{}"),
            ("priva.producer", PRODUCER),
            ("priva.builder_source_sha256", older.as_str()),
            ("priva.keys", "[]"),
        ];
        let footer = footer
            .into_iter()
            .map(|(key, value)| (key.to_owned(), value.to_owned()))
            .collect();
        parquet::publish(&path, footer, |_| Ok(())).unwrap();
        assert!(check(&path, None).is_err());
        assert!(check(&path, Some(older.as_str())).is_ok());
        assert!(check(&path, Some("c".repeat(64).as_str())).is_err());
        assert!(check(&path, Some(identity().as_str())).is_err());
        // Merge and export never take a pin and keep requiring the current builder.
        assert!(checked(&path, None).is_err());
        std::fs::remove_file(&path).unwrap();
    }

    #[test]
    fn builder_pin_is_a_sha256_check_option() {
        let message = |command: &str, pin: &str| {
            let args = [command, "--input", "absent", "--builder-sha256", pin];
            let args: Vec<String> = args.into_iter().map(String::from).collect();
            dispatch(&args).unwrap_err().to_string()
        };
        let refused = "--builder-sha256 requires SHA256";
        assert_eq!(message("check", "invalid"), refused);
        let pin = "0f".repeat(32);
        for command in ["identity", "merge", "export"] {
            assert!(message(command, &pin).starts_with("usage:"));
        }
    }
}
