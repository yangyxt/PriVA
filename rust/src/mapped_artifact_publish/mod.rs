//! Optional port of mapped_artifact_publish.py. Scientific commands remain
//! explicit external argument vectors. No installed caller is changed.
mod cli;
mod metadata;
mod transaction;
use crate::alpha_annotation_base::io::*;
use crate::biomart;
use crate::cli::{Result, invalid};
use crate::json_text::Json;
use crate::python_path::{absolute, is_file};
pub(crate) use cli::dispatch;
use metadata::{Value, require_equal};
use std::path::{Path, PathBuf};

const SCHEMA: &str = "priva-mapped-artifact-publication-v1";
/// The only native producer: `priva pm1-regions build`, run by this executable.
const NATIVE_PM1: &str = "pm1-regions";
/// Options of `priva pm1-regions build`; each must be given once as `--name value`.
const NATIVE_PM1_OPTIONS: [&str; 8] = [
    "assembly",
    "nm-to-enst",
    "intolerant-domains",
    "hcseeker",
    "rmc",
    "output",
    "rmc-oe-threshold",
    "rmc-p-bh-threshold",
];
fn error(category: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    crate::builder_error::category(category, invalid(message))
}
fn require(ok: bool, message: impl Into<String>) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(error("ValueError", message))
    }
}
fn object<const N: usize>(values: [(&str, Json); N]) -> Json {
    Json::Object(values.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}
/// Identity of the native producer recorded in the publication context and
/// receipt, as printed by `priva builder-identity pm1-regions`. Recovery
/// compares it, so a changed builder cannot finish another builder's pair.
fn native_producer_identity() -> Result<Json> {
    let identity = crate::builder_identity::identity(NATIVE_PM1)
        .ok_or_else(|| error("ValueError", "No pm1-regions builder identity"))?;
    let command = ["priva", NATIVE_PM1, "build"];
    Ok(object([
        ("builder", string(NATIVE_PM1)),
        ("builder_source_sha256", string(identity)),
        ("command", Json::Array(command.map(string).to_vec())),
    ]))
}
fn canonical_assembly(assembly: &str) -> &'static str {
    match assembly {
        "hg19" | "GRCh37" => "hg19",
        "hg38" | "GRCh38" => "hg38",
        _ => "",
    }
}
struct Context {
    artifact: PathBuf,
    mapping: PathBuf,
    inputs: Vec<PathBuf>,
    cache_root: PathBuf,
    output_dir: PathBuf,
    assembly: String,
    config: PathBuf,
    producer: PathBuf,
    /// This executable when `pm1-regions build` replaces the Python producer.
    native_producer: Option<PathBuf>,
    kind: String,
    annotation_base: Option<PathBuf>,
    native_lineage: Option<PathBuf>,
    helper: PathBuf,
    python: String,
    biomart: biomart::Context,
}
impl Context {
    fn receipt_check(&self, path: &Path) -> Result<()> {
        biomart::check_receipt(&self.biomart, path, &self.mapping, &self.inputs)
    }
    fn identity(&self) -> Result<Json> {
        biomart::mapping_identity(&self.biomart, &self.mapping, &self.inputs)
    }
    /// Bind the native builder's explicit inputs to this publication: the
    /// selected assembly, the selected mapping and its HCSeeker/RMC anchors.
    /// Only exact, distinct `--name value` pairs are admitted, so the values
    /// checked here are the values the builder reads.
    fn native_arguments(&self, argv: &[String]) -> Result<()> {
        use crate::python_path::path_value;
        let mut values = std::collections::HashMap::new();
        for pair in argv.chunks(2) {
            let name = pair[0].strip_prefix("--").unwrap_or_default();
            require(
                pair.len() == 2
                    && NATIVE_PM1_OPTIONS.contains(&name)
                    && values.insert(name, pair[1].as_str()).is_none(),
                "Native PM1 arguments must be distinct exact --name value pairs",
            )?;
        }
        let path = |key: &str| values.get(key).map(|v| path_value(v));
        let anchored = |key: &str| path(key).is_some_and(|p| self.inputs.contains(&p));
        require(
            values.get("assembly").copied() == Some(canonical_assembly(&self.assembly))
                && path("nm-to-enst").as_ref() == Some(&self.mapping)
                && values.contains_key("intolerant-domains")
                && anchored("hcseeker")
                && anchored("rmc"),
            "Native PM1 inputs must be the selected assembly, mapping and HCSeeker/RMC anchors",
        )
    }
    fn artifact_format(&self) -> &'static str {
        if self.native_producer.is_some() {
            "rkyv"
        } else {
            "pickle"
        }
    }
    fn expected(&self) -> Result<Json> {
        if self.native_producer.is_some() {
            require(
                self.kind == "pm1"
                    && self.annotation_base.is_none()
                    && self.native_lineage.is_none(),
                "The native producer builds only PM1 regions",
            )?;
            require(
                self.artifact.extension().is_some_and(|e| e == "rkyv"),
                "Native PM1 regions are published as an .rkyv artifact",
            )?;
        } else {
            let wanted = if self.kind == "motifs" {
                "am_pick_missense_constrained_segments.py"
            } else {
                "prepare_pm1_regions.py"
            };
            require(
                self.producer.is_absolute()
                    && self.producer.file_name().is_some_and(|n| n == wanted),
                "Use the selected motif/PM1 scientific producer",
            )?;
            require(
                !symlink(&self.producer)? && is_file(&self.producer)?,
                "Scientific producer is not a regular file",
            )?;
        }
        require(
            self.config.is_absolute() && !symlink(&self.config)? && is_file(&self.config)?,
            "Use the actual regular runtime config",
        )?;
        let mut result = object([
            ("kind", string(&self.kind)),
            ("assembly", string(&self.assembly)),
            ("uid", Json::Int(unsafe { libc::getuid() } as i128)),
            ("cache_root", pstr(&resolve(&self.cache_root)?)),
            ("output_dir", pstr(&resolve(&self.output_dir)?)),
        ]);
        if self.native_producer.is_some() {
            result.set("native_producer", native_producer_identity()?);
        } else {
            result.set("producer", pstr(&self.producer));
            result.set("producer_sha256", string(digest(&self.producer)?));
        }
        result.set("config", pstr(&self.config));
        result.set("config_sha256", string(digest(&self.config)?));
        result.set(
            "helper_sha256",
            string(digest(
                &self.helper.join("scripts/mapped_artifact_publish.py"),
            )?),
        );
        result.set(
            "mapping_helper_sha256",
            string(digest(&self.helper.join("scripts/biomart_nm_to_enst.py"))?),
        );
        result.set("mapping_identity", self.identity()?);
        if let Some(base) = &self.annotation_base {
            require(
                self.kind == "motifs",
                "Only motifs consume the complete AlphaMissense base",
            )?;
            let source = self
                .producer
                .parent()
                .and_then(Path::parent)
                .ok_or_else(|| error("IndexError", "1"))?;
            // The original imports these two source assets before reading the
            // selected config or base. Preserve missing/unreadable-asset order;
            // their implementation is native and no source code is executed.
            for name in ["alphamissense_annotation_base.py", "alphamissense_cache.py"] {
                let asset = self.producer.with_file_name(name);
                std::fs::read(&asset).map_err(|e| crate::python_path::os_error(e, &asset, None))?;
            }
            result.set(
                "scientific_input",
                crate::alpha_annotation_base::motif_binding(
                    &self.config,
                    source,
                    Some(base),
                    &self.python,
                )?,
            );
        }
        if let Some(lineage) = &self.native_lineage {
            result.set(
                "native_input_lineage",
                object([
                    ("path", string(absolute(lineage)?)),
                    ("sha256", string(digest(lineage)?)),
                ]),
            );
        }
        Ok(result)
    }
    fn target(&self) -> Result<PathBuf> {
        require(
            self.artifact.is_absolute()
                && !self
                    .artifact
                    .components()
                    .any(|c| c == std::path::Component::ParentDir),
            "Use an absolute selected artifact path",
        )?;
        biomart::managed_artifact_path(
            &self.cache_root,
            self.artifact.parent().unwrap_or(&self.artifact),
            Path::new(self.artifact.file_name().unwrap_or_default()),
        )
    }
    fn recovery_context(&self, target: &Path) -> Result<Json> {
        let mut expected = self.expected()?;
        let pending = transaction::load_pending(&suffix(target, ".publish.pending.json"))?;
        let pending = Value(pending.trim());
        let pinned = pending.get("context")?.get("native_input_lineage")?;
        if pinned.0 != "null" && expected.get("native_input_lineage").is_none() {
            require(
                self.kind == "motifs"
                    && expected.get("scientific_input").is_some()
                    && pinned.exact_object(&["path", "sha256"])?,
                "Unexpected native adoption context",
            )?;
            let selected = path(&pinned.item("path")?.json()?)?;
            require(
                selected.is_absolute()
                    && !symlink(&selected)?
                    && is_file(&selected)?
                    && pinned.item("sha256")?.equal(&string(digest(&selected)?))?,
                "Pinned native lineage changed",
            )?;
            let encoded = read_text(&selected)?;
            let lineage = Value::validate(&encoded)?;
            require(
                lineage
                    .get("schema")?
                    .equal(&string("priva-alphamissense-motif-input-v1"))?
                    && lineage.get("provenance")?.truth()
                    && lineage
                        .get("scientific_input")?
                        .equal(get(&expected, "scientific_input")?)?
                    && lineage
                        .get("mapping_identity")?
                        .equal(get(&expected, "mapping_identity")?)?
                    && lineage.get("artifact_sha256")?.equal(
                        &pending
                            .item("new")?
                            .item("artifact")?
                            .item("sha256")?
                            .json()?,
                    )?
                    && lineage.get("artifact_bytes")?.equal(
                        &pending
                            .item("new")?
                            .item("artifact")?
                            .item("size")?
                            .json()?,
                    )?,
                "Pinned native lineage no longer describes the pending motif",
            )?;
            let original = path(&lineage.item("artifact_path")?.json()?)?;
            self.receipt_check(&original)?;
            require(
                lineage
                    .item("artifact_sha256")?
                    .equal(get(&transaction::file_state(&original)?, "sha256")?)?
                    && lineage.item("receipt_sha256")?.equal(get(
                        &transaction::file_state(&suffix(&original, ".nm_to_enst.json"))?,
                        "sha256",
                    )?)?,
                "Preserved historical native pair changed",
            )?;
            expected.set("native_input_lineage", pinned.json()?);
        }
        Ok(expected)
    }
    fn recover(&self) -> Result<bool> {
        let mut target = self.artifact.clone();
        if !lexists(&suffix(&target, ".publish.pending.json"))? {
            if !(lexists(&target)? || lexists(&suffix(&target, ".nm_to_enst.json"))?) {
                return Ok(false);
            }
            if !self.cache_root.is_absolute() || !self.output_dir.is_absolute() {
                return Ok(false);
            }
            let identity = self.identity()?;
            let candidate = biomart::changed_target(
                &identity,
                &self.cache_root,
                &self.output_dir,
                &self.assembly,
                &self.kind,
                self.artifact_format(),
            )?;
            if !lexists(&suffix(&candidate, ".publish.pending.json"))? {
                return Ok(false);
            }
            require(
                !biomart::artifact_receipt(&target)?.matches(&identity)?,
                "Pending mapping namespace is not the selected artifact's changed mapping",
            )?;
            target = candidate;
        } else {
            target = self.target()?;
        }
        let _lock = lock(target.parent().unwrap())?;
        if lexists(&suffix(&target, ".publish.pending.json"))? {
            self.finish(&target, &self.recovery_context(&target)?)?;
            println!(
                "Recovered completed mapped artifact publication: {}",
                self.artifact.display()
            );
            return Ok(true);
        }
        Ok(false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::os::unix::fs::PermissionsExt;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static NEXT: AtomicUsize = AtomicUsize::new(0);
    /// Stands in for this executable: records its argv and writes `--output`.
    const PRODUCER: &str = "#!/bin/sh\n\
        printf '%s\\n' \"$@\" > \"${0%/*}/argv.txt\"\n\
        while [ $# -gt 1 ]; do\n\
        if [ \"$1\" = --output ]; then printf archive > \"$2\"; fi\n\
        shift\n\
        done\n\
        exit 0\n";

    struct Temp(PathBuf);
    impl Drop for Temp {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    fn fixture() -> (Temp, Context) {
        let root = std::env::temp_dir().join(format!(
            "priva-mapped-artifact-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&root).unwrap();
        let root = root.canonicalize().unwrap();
        let write = |name: &str, text: &str| {
            let path = root.join(name);
            std::fs::write(&path, text).unwrap();
            path
        };
        let producer = write("priva", PRODUCER);
        std::fs::set_permissions(&producer, std::fs::Permissions::from_mode(0o755)).unwrap();
        let helper = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
        let output = root.join("AlphaMissense");
        let context = Context {
            artifact: output.join("pm1_regions.hg38.rkyv"),
            mapping: write("selected.tsv", "Gene stable ID\tTranscript stable ID\n"),
            inputs: vec![
                write("HC_spots.hg38.tsv", "hotspots\n"),
                write("gnomad_v2.1.1_rmc_flat.hg38.tsv", "rmc\n"),
            ],
            cache_root: root.clone(),
            output_dir: output,
            assembly: "hg38".to_owned(),
            config: write("config.yaml", "assembly: hg38\n"),
            producer: PathBuf::new(),
            native_producer: Some(producer),
            kind: "pm1".to_owned(),
            annotation_base: None,
            native_lineage: None,
            helper: helper.to_path_buf(),
            python: "python".to_owned(),
            biomart: biomart::Context::new(helper).unwrap(),
        };
        (Temp(root), context)
    }

    fn native_argv(context: &Context) -> Vec<String> {
        let domains = context.cache_root.join("domains.rkyv");
        let pairs = [
            ("--assembly", Path::new("hg38")),
            ("--nm-to-enst", context.mapping.as_path()),
            ("--intolerant-domains", domains.as_path()),
            ("--hcseeker", context.inputs[0].as_path()),
            ("--rmc", context.inputs[1].as_path()),
            ("--output", context.artifact.as_path()),
        ];
        let mut argv = vec!["--".to_owned()];
        for (name, value) in pairs {
            argv.extend([name.to_owned(), value.display().to_string()]);
        }
        argv
    }

    #[test]
    fn native_context_replaces_the_python_producer_and_requires_rkyv() {
        let (_temp, mut context) = fixture();
        let expected = context.expected().unwrap();
        let identity = native_producer_identity().unwrap();
        let recorded = get(&expected, "native_producer").unwrap();
        assert!(equal(recorded, &identity));
        let builder = crate::builder_identity::identity("pm1-regions").unwrap();
        let source = get(&identity, "builder_source_sha256").unwrap();
        assert!(equal(source, &string(builder)));
        assert!(expected.get("producer").is_none());
        assert!(expected.get("producer_sha256").is_none());
        assert_eq!(context.artifact_format(), "rkyv");
        context.artifact.set_extension("pkl");
        assert!(context.expected().is_err());
        context.artifact.set_extension("rkyv");
        context.kind = "motifs".to_owned();
        assert!(context.expected().is_err());
    }

    #[test]
    fn native_arguments_bind_the_selected_mapping_assembly_and_anchors() {
        let (_temp, context) = fixture();
        let argv = native_argv(&context)[1..].to_vec();
        context.native_arguments(&argv).unwrap();
        let mut with_threshold = argv.clone();
        with_threshold.extend(["--rmc-oe-threshold".to_owned(), "0.4".to_owned()]);
        context.native_arguments(&with_threshold).unwrap();
        let mut repeated = with_threshold.clone();
        repeated.extend(["--rmc-oe-threshold".to_owned(), "0.5".to_owned()]);
        let replace = |at: usize, value: &str| {
            let mut changed = argv.clone();
            changed[at] = value.to_owned();
            changed
        };
        let other = context.cache_root.join("other.tsv").display().to_string();
        let rejected = [
            replace(1, "hg19"),
            replace(3, &other),
            replace(7, &other),
            replace(9, &other),
            replace(0, "--assem"),
            replace(2, "--assembly"),
            replace(4, "--threads"),
            argv[..argv.len() - 1].to_vec(),
            repeated,
        ];
        for changed in rejected {
            assert!(context.native_arguments(&changed).is_err());
        }
    }

    #[test]
    fn native_run_publishes_the_rkyv_pair_with_its_producer_identity() {
        let (temp, context) = fixture();
        let identity = context.identity().unwrap();
        let identity = crate::json_text::dumps(&identity, true).unwrap();
        context.run(&identity, &native_argv(&context)).unwrap();
        assert_eq!(std::fs::read(&context.artifact).unwrap(), b"archive");
        let argv = std::fs::read_to_string(temp.0.join("argv.txt")).unwrap();
        assert!(argv.starts_with("pm1-regions\nbuild\n--assembly\nhg38\n"));
        assert!(!argv.contains(&context.artifact.display().to_string()));
        let text = read_text(&suffix(&context.artifact, ".nm_to_enst.json")).unwrap();
        let receipt = Value::validate(&text).unwrap();
        let recorded = receipt.get("native_producer").unwrap();
        let identity = native_producer_identity().unwrap();
        assert!(recorded.equal(&identity).unwrap());
        let pending = suffix(&context.artifact, ".publish.pending.json");
        assert!(!lexists(&pending).unwrap());
        context.receipt_check(&context.artifact).unwrap();
        assert!(!context.recover().unwrap());
    }
}
