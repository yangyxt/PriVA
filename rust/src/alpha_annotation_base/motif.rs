//! Explicit adoption of a completed motif pair with the native BioMart contract.
use super::io::*;
use super::publication::{Context, selected};
use super::{MOTIF_SCHEMA, error, object, profile, require};
use crate::cli::Result;
use crate::json_text::Json;
use crate::python_path::{absolute, is_file, metadata, os_error, tempfile_parent};
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::os::unix::fs::MetadataExt;
use std::path::{Path, PathBuf};

pub(super) fn binding(context: &Context<'_>) -> Result<Json> {
    let record = context.check_fields(false, &["counts"])?;
    Ok(object([
        ("schema", string(MOTIF_SCHEMA)),
        (
            "annotation_base_sha256",
            get(&record, "base_sha256")?.clone(),
        ),
        (
            "annotation_science_sha256",
            string(json_digest(get(&record, "science")?)?),
        ),
        ("counts", get(&record, "counts")?.clone()),
        (
            "producer_sha256",
            string(digest(
                &context
                    .source
                    .join("scripts/am_pick_missense_constrained_segments.py"),
            )?),
        ),
        ("runtime", profile::motif_runtime(context.python)?),
    ]))
}
pub(super) fn check(context: &Context<'_>, artifact: &Path) -> Result<()> {
    let encoded = read_text(&suffix(artifact, ".nm_to_enst.json"))?;
    // Validate the receipt before observing the scientific runtime, matching
    // the original json.loads -> motif_binding evaluation order.
    let record = compare_metadata(&encoded, &object([]))?;
    optional(&record, "scientific_input")?;
    let expected = binding(context)?;
    let record = compare_metadata(&encoded, &object([("scientific_input", expected.clone())]))?;
    require(
        equal(optional(&record, "scientific_input")?, &expected),
        "Motif receipt is not bound to this complete annotation and scientific producer",
    )
}
fn file_state(path: &Path) -> Result<Json> {
    if !lexists(path)? {
        return Ok(Json::Null);
    }
    let before = fs::symlink_metadata(path).map_err(|e| os_error(e, path, None))?;
    require(
        before.is_file() && before.nlink() == 1,
        format!(
            "Publication files must be regular, unlinked files: {}",
            path.display()
        ),
    )?;
    let hash = digest(path)?;
    let after = fs::symlink_metadata(path).map_err(|e| os_error(e, path, None))?;
    fn stable(m: &fs::Metadata) -> (u64, u64, i64, i64, i64, i64, u32, u64, u32) {
        (
            m.ino(),
            m.len(),
            m.mtime(),
            m.mtime_nsec(),
            m.ctime(),
            m.ctime_nsec(),
            m.mode(),
            m.nlink(),
            m.uid(),
        )
    }
    require(
        stable(&before) == stable(&after),
        format!("Publication file changed while reading: {}", path.display()),
    )?;
    Ok(object([
        ("inode", Json::Int(before.ino() as i128)),
        ("size", Json::Int(before.len() as i128)),
        ("mtime_ns", time(before.mtime(), before.mtime_nsec())),
        ("ctime_ns", time(before.ctime(), before.ctime_nsec())),
        ("mode", Json::Int(before.mode() as i128)),
        ("nlink", Json::Int(before.nlink() as i128)),
        ("uid", Json::Int(before.uid() as i128)),
        ("sha256", string(hash)),
    ]))
}
fn moved(actual: &Json, expected: &Json) -> Result<bool> {
    if matches!(actual, Json::Null) || matches!(expected, Json::Null) {
        return Ok(false);
    }
    for (k, v) in members(expected)? {
        if k != "ctime_ns" && !equal(get(actual, k)?, v) {
            return Ok(false);
        }
    }
    Ok(true)
}
fn exact_keys(value: &Json, keys: &[&str]) -> bool {
    value
        .as_object()
        .is_some_and(|v| v.len() == keys.len() && v.iter().all(|(k, _)| keys.contains(&k.as_str())))
}
fn size(value: &Json) -> Option<i128> {
    match value {
        Json::Int(v) => Some(*v),
        Json::Bool(v) => Some(i128::from(*v)),
        _ => None,
    }
}

struct Publisher<'a> {
    context: &'a Context<'a>,
    artifact: &'a Path,
    mapping: &'a Path,
    inputs: &'a [PathBuf],
    lineage_file: &'a Path,
    root: PathBuf,
    base: PathBuf,
    assembly: String,
}
impl Publisher<'_> {
    fn biomart(&self, artifact: Option<&Path>) -> Result<Json> {
        // The original helper location still controls bundle recognition and
        // provenance. The qualified native API preserves the identity/check
        // order without a Python implementation adapter.
        let context = crate::biomart::Context::new(self.context.helper)?;
        if let Some(artifact) = artifact {
            crate::biomart::check_receipt(&context, artifact, self.mapping, self.inputs)?;
            Ok(Json::Null)
        } else {
            crate::biomart::mapping_identity(&context, self.mapping, self.inputs)
        }
    }
    fn expected(&self) -> Result<Json> {
        let producer = self
            .context
            .source
            .join("scripts/am_pick_missense_constrained_segments.py");
        require(
            producer.is_absolute()
                && producer
                    .file_name()
                    .is_some_and(|n| n == "am_pick_missense_constrained_segments.py"),
            "Use the selected motif/PM1 scientific producer",
        )?;
        require(
            !symlink(&producer)? && is_file(&producer)?,
            "Scientific producer is not a regular file",
        )?;
        require(
            self.context.config.is_absolute()
                && !symlink(self.context.config)?
                && is_file(self.context.config)?,
            "Use the actual regular runtime config",
        )?;
        let mut result = object([
            ("kind", string("motifs")),
            ("assembly", string(&self.assembly)),
            ("uid", Json::Int(unsafe { libc::getuid() } as i128)),
            ("cache_root", pstr(&resolve(&self.root)?)),
            (
                "output_dir",
                pstr(&resolve(self.artifact.parent().unwrap_or(Path::new(".")))?),
            ),
            ("producer", pstr(&producer)),
            ("producer_sha256", string(digest(&producer)?)),
            ("config", pstr(self.context.config)),
            ("config_sha256", string(digest(self.context.config)?)),
            (
                "helper_sha256",
                string(digest(
                    &self
                        .context
                        .helper
                        .join("scripts/mapped_artifact_publish.py"),
                )?),
            ),
            (
                "mapping_helper_sha256",
                string(digest(
                    &self.context.helper.join("scripts/biomart_nm_to_enst.py"),
                )?),
            ),
            ("mapping_identity", self.biomart(None)?),
        ]);
        let bound = Context {
            config: self.context.config,
            source: self.context.source,
            base: Some(&self.base),
            python: self.context.python,
            helper: self.context.source,
        };
        result.set("scientific_input", binding(&bound)?);
        result.set(
            "native_input_lineage",
            object([
                ("path", string(absolute(self.lineage_file)?)),
                ("sha256", string(digest(self.lineage_file)?)),
            ]),
        );
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
        let output = self.artifact.parent().unwrap();
        require(
            self.root.is_absolute() && output.is_absolute(),
            "Changed mapping outputs require absolute cache_root and output_dir",
        )?;
        let logical_root = tempfile_parent(&self.root)?;
        let root = crate::python_io::resolve_strict(&logical_root)
            .map_err(|(e, p)| os_error(e, &p, None))?;
        require(
            root.is_dir() && root.parent() != Some(root.as_path()),
            "Changed mapping outputs require a configured cache_root directory",
        )?;
        let output = tempfile_parent(output)?;
        let name = self.artifact.file_name().unwrap_or_default();
        let target = if output.starts_with(&logical_root) {
            root.join(output.strip_prefix(&logical_root).unwrap())
                .join(name)
        } else if output.starts_with(&root) {
            output.join(name)
        } else {
            return Err(error(
                "ValueError",
                "Changed mapping output_dir escapes cache_root",
            ));
        };
        let mut current = root.clone();
        for part in target.strip_prefix(&root).unwrap().components() {
            current.push(part);
            require(
                !symlink(&current)?,
                format!(
                    "Linked managed mapping output is not allowed: {}",
                    current.display()
                ),
            )?;
            require(
                current == target || !exists(&current)? || current.is_dir(),
                format!(
                    "Managed mapping output parent is not a directory: {}",
                    current.display()
                ),
            )?;
        }
        Ok(target)
    }
    fn finish(&self, target: &Path, expected: &Json) -> Result<()> {
        let marker = suffix(target, ".publish.pending.json");
        let before = file_state(&marker)?;
        require(
            !matches!(before, Json::Null)
                && size(get(&before, "size")?).is_some_and(|v| v <= 65536),
            "Missing/oversized pending publication",
        )?;
        let value = read_json(&marker)?;
        require(
            equal(&file_state(&marker)?, &before),
            "Pending publication changed while reading",
        )?;
        require(
            exact_keys(
                &value,
                &[
                    "schema",
                    "context",
                    "target",
                    "stage",
                    "prior",
                    "new",
                    "producer_job_id",
                ],
            ) && equal(
                get(&value, "schema")?,
                &string("priva-mapped-artifact-publication-v1"),
            ),
            "Unknown pending publication format",
        )?;
        require(
            equal(get(&value, "context")?, expected)
                && equal(get(&value, "target")?, &pstr(target)),
            "Pending publication producer/config/mapping/target differs",
        )?;
        let stage = path(get(&value, "stage")?)?;
        require(
            stage.parent() == target.parent()
                && stage
                    .file_name()
                    .is_some_and(|s| s.to_string_lossy().starts_with(".mapped-artifact-motifs-"))
                && !symlink(&stage)?
                && stage.is_dir()
                && metadata(&stage)?.uid() == unsafe { libc::getuid() }
                && metadata(&stage)?.mode() & 0o777 == 0o700,
            "Pending publication stage is not owned beside target",
        )?;
        require(
            exact_keys(get(&value, "prior")?, &["artifact", "receipt"])
                && exact_keys(get(&value, "new")?, &["artifact", "receipt"]),
            "Pending pair fields differ",
        )?;
        let artifact = stage.join(target.file_name().unwrap());
        let paths = [
            ("artifact", artifact.clone(), target.to_owned()),
            (
                "receipt",
                suffix(&artifact, ".nm_to_enst.json"),
                suffix(target, ".nm_to_enst.json"),
            ),
        ];
        let mut locations = Vec::new();
        for (key, staged, final_path) in &paths {
            let staged_state = file_state(staged)?;
            let final_state = file_state(final_path)?;
            let new = get(get(&value, "new")?, key)?;
            require(
                exact_keys(
                    new,
                    &[
                        "inode", "size", "mtime_ns", "ctime_ns", "mode", "nlink", "uid", "sha256",
                    ],
                ) && size(optional(new, "size")?).is_some_and(|v| v > 0),
                "Pending staged file identity is incomplete",
            )?;
            if !matches!(staged_state, Json::Null) {
                require(
                    moved(&staged_state, new)?
                        && equal(&final_state, get(get(&value, "prior")?, key)?),
                    "Pending stage or previous selected artifact changed",
                )?;
                locations.push("stage");
            } else {
                require(
                    moved(&final_state, new)?,
                    "Pending published file changed or disappeared",
                )?;
                locations.push("final");
            }
        }
        require(
            locations[0] <= locations[1],
            "Pending publication is not an ordered two-file state",
        )?;
        let encoded = read_text(if locations[1] == "stage" {
            &paths[1].1
        } else {
            &paths[1].2
        })?;
        let wanted = artifact_record(expected, get(get(&value, "new")?, "artifact")?)?;
        let record = compare_metadata(&encoded, &wanted)?;
        require(
            equal(&record, &wanted),
            "Pending mapping receipt does not describe the staged artifact",
        )?;
        for (location, (_, staged, final_path)) in locations.iter().zip(&paths) {
            if *location == "stage" {
                rename(staged, final_path)?;
            }
        }
        self.biomart(Some(target))?;
        for (key, _, final_path) in &paths {
            require(
                moved(&file_state(final_path)?, get(get(&value, "new")?, key)?)?,
                "Published pair changed before commit",
            )?;
        }
        unlink(&marker)?;
        if fs::read_dir(&stage)
            .map_err(|e| os_error(e, &stage, None))?
            .next()
            .is_none()
        {
            fs::remove_dir(&stage).map_err(|e| os_error(e, &stage, None))?;
        }
        Ok(())
    }
}
fn artifact_record(expected: &Json, state: &Json) -> Result<Json> {
    let mut value = object([
        (
            "mapping_identity",
            get(expected, "mapping_identity")?.clone(),
        ),
        ("artifact_bytes", get(state, "size")?.clone()),
        ("artifact_sha256", get(state, "sha256")?.clone()),
    ]);
    for k in ["scientific_input", "native_input_lineage"] {
        if let Some(v) = expected.get(k) {
            value.set(k, v.clone());
        }
    }
    Ok(value)
}
// BioMart's existing atomic_json uses a .part file, not tempfile+fsync.
fn mapping_json(path: &Path, value: &Json) -> Result<()> {
    let temporary = suffix(path, ".part");
    require(
        !symlink(&temporary)?,
        format!(
            "Refusing to replace a symlink used as a staging file: {}",
            temporary.display()
        ),
    )?;
    let mut file = File::create(&temporary).map_err(|e| os_error(e, &temporary, None))?;
    file.write_all(formatted_json(value)?.as_bytes())
        .map_err(|e| os_error(e, &temporary, None))?;
    drop(file);
    rename(&temporary, path)
}
pub(super) fn adopt(
    context: &Context<'_>,
    artifact: &Path,
    lineage_file: &Path,
    mapping: &Path,
    inputs: &[PathBuf],
) -> Result<()> {
    let binding = binding(context)?;
    let lineage = read_json(lineage_file)?;
    require(
        equal(optional(&lineage, "schema")?, &string(MOTIF_SCHEMA))
            && equal(optional(&lineage, "scientific_input")?, &binding)
            && truth(optional(&lineage, "provenance")?),
        "Preserved native motif lineage differs from the complete base/producer",
    )?;
    let settings = context.settings()?;
    let publisher = Publisher {
        context,
        artifact,
        mapping,
        inputs,
        lineage_file,
        root: path(get(&settings.config, "cache_root")?)?,
        base: selected(&settings.paths, "base")?,
        assembly: profile::assembly(get(&settings.config, "assembly")?)?.into(),
    };
    let source = path(get(&lineage, "artifact_path")?)?;
    let target = publisher.target()?;
    require(
        resolve(&source)? != resolve(&target)?,
        "Native motif adoption must select a new managed target",
    )?;
    require(
        !symlink(&source)? && is_file(&source)?,
        "Preserved native motif must be an independent regular file",
    )?;
    mkdir(target.parent().unwrap())?;
    let _lock = lock(target.parent().unwrap())?;
    let expected = publisher.expected()?;
    require(
        equal(
            optional(&lineage, "scientific_input")?,
            optional(&expected, "scientific_input")?,
        ),
        "Preserved native motif base/producer binding differs",
    )?;
    let pending = suffix(&target, ".publish.pending.json");
    if lexists(&pending)? {
        return publisher.finish(&target, &expected);
    }
    if lexists(&target)? || lexists(&suffix(&target, ".nm_to_enst.json"))? {
        publisher.biomart(Some(&target))?;
        let completed = compare_metadata(
            &read_text(&suffix(&target, ".nm_to_enst.json"))?,
            &object([
                (
                    "scientific_input",
                    optional(&expected, "scientific_input")?.clone(),
                ),
                (
                    "native_input_lineage",
                    optional(&expected, "native_input_lineage")?.clone(),
                ),
                (
                    "artifact_sha256",
                    optional(&lineage, "artifact_sha256")?.clone(),
                ),
            ]),
        )?;
        require(
            equal(
                optional(&completed, "scientific_input")?,
                optional(&expected, "scientific_input")?,
            ) && equal(
                optional(&completed, "native_input_lineage")?,
                optional(&expected, "native_input_lineage")?,
            ) && equal(
                optional(&completed, "artifact_sha256")?,
                optional(&lineage, "artifact_sha256")?,
            ),
            "Native adoption must not replace an existing artifact or historical receipt",
        )?;
        return Ok(());
    }
    publisher.biomart(Some(&source))?;
    let source_before = file_state(&source)?;
    let source_receipt = suffix(&source, ".nm_to_enst.json");
    let receipt_before = file_state(&source_receipt)?;
    let original = compare_metadata(
        &read_text(&source_receipt)?,
        &object([
            (
                "mapping_identity",
                optional(&lineage, "mapping_identity")?.clone(),
            ),
            (
                "artifact_bytes",
                optional(&lineage, "artifact_bytes")?.clone(),
            ),
            (
                "artifact_sha256",
                optional(&lineage, "artifact_sha256")?.clone(),
            ),
        ]),
    )?;
    for k in ["mapping_identity", "artifact_bytes", "artifact_sha256"] {
        require(
            equal(optional(&lineage, k)?, optional(&original, k)?),
            "Preserved native motif identity differs",
        )?;
    }
    require(
        equal(
            optional(&lineage, "receipt_sha256")?,
            get(&receipt_before, "sha256")?,
        ),
        "Preserved native receipt SHA256 differs",
    )?;
    let stage = stage(target.parent().unwrap(), ".mapped-artifact-motifs-")?;
    let staged = stage.join(target.file_name().unwrap());
    // This original copy has no fsync, unlike VCF publication's copy_verified.
    let mut incoming = File::open(&source).map_err(|e| os_error(e, &source, None))?;
    let mut output = OpenOptions::new()
        .create_new(true)
        .write(true)
        .open(&staged)
        .map_err(|e| os_error(e, &staged, None))?;
    let mut buf = vec![0; 8 * 1024 * 1024];
    loop {
        let n = incoming
            .read(&mut buf)
            .map_err(|e| os_error(e, &source, None))?;
        if n == 0 {
            break;
        }
        output
            .write_all(&buf[..n])
            .map_err(|e| os_error(e, &staged, None))?;
    }
    drop(output);
    drop(incoming);
    let staged_state = file_state(&staged)?;
    require(
        equal(
            get(&staged_state, "sha256")?,
            get(&source_before, "sha256")?,
        ) && equal(&file_state(&source)?, &source_before)
            && equal(&file_state(&source_receipt)?, &receipt_before),
        "Preserved native motif changed while copying",
    )?;
    mapping_json(
        &suffix(&staged, ".nm_to_enst.json"),
        &artifact_record(&expected, &staged_state)?,
    )?;
    require(
        equal(&publisher.expected()?, &expected),
        "Native motif inputs changed before publication",
    )?;
    let value = object([
        ("schema", string("priva-mapped-artifact-publication-v1")),
        ("context", expected.clone()),
        ("target", pstr(&target)),
        ("stage", pstr(&stage)),
        (
            "prior",
            object([("artifact", Json::Null), ("receipt", Json::Null)]),
        ),
        (
            "new",
            object([
                ("artifact", file_state(&staged)?),
                ("receipt", file_state(&suffix(&staged, ".nm_to_enst.json"))?),
            ]),
        ),
        (
            "producer_job_id",
            std::env::var("SLURM_JOB_ID")
                .map(string)
                .unwrap_or(Json::Null),
        ),
    ]);
    mapping_json(&pending, &value)?;
    publisher.finish(&target, &expected)
}
