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
struct Context {
    artifact: PathBuf,
    mapping: PathBuf,
    inputs: Vec<PathBuf>,
    cache_root: PathBuf,
    output_dir: PathBuf,
    assembly: String,
    config: PathBuf,
    producer: PathBuf,
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
    fn expected(&self) -> Result<Json> {
        let wanted = if self.kind == "motifs" {
            "am_pick_missense_constrained_segments.py"
        } else {
            "prepare_pm1_regions.py"
        };
        require(
            self.producer.is_absolute() && self.producer.file_name().is_some_and(|n| n == wanted),
            "Use the selected motif/PM1 scientific producer",
        )?;
        require(
            !symlink(&self.producer)? && is_file(&self.producer)?,
            "Scientific producer is not a regular file",
        )?;
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
            ("producer", pstr(&self.producer)),
            ("producer_sha256", string(digest(&self.producer)?)),
            ("config", pstr(&self.config)),
            ("config_sha256", string(digest(&self.config)?)),
            (
                "helper_sha256",
                string(digest(
                    &self.helper.join("scripts/mapped_artifact_publish.py"),
                )?),
            ),
            (
                "mapping_helper_sha256",
                string(digest(&self.helper.join("scripts/biomart_nm_to_enst.py"))?),
            ),
            ("mapping_identity", self.identity()?),
        ]);
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
