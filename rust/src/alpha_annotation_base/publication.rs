//! Ordered, restartable publication of converted/base/runtime VCF-index pairs.
//! A completed rename prefix is recoverable; no receipt authorizes regeneration.
use super::io::*;
use super::{SCHEMA, error, header, object, profile, require};
use crate::cli::Result;
use crate::json_text::Json;
use crate::python_path::{absolute, is_file, metadata};
use std::path::{Path, PathBuf};

pub(super) struct Context<'a> {
    pub config: &'a Path,
    pub source: &'a Path,
    pub base: Option<&'a Path>,
    pub python: &'a str,
    pub helper: &'a Path,
}
pub(super) fn selected(paths: &Json, key: &str) -> Result<PathBuf> {
    path(get(paths, key)?)
}
pub(super) fn receipt(base: &Path) -> PathBuf {
    suffix(base, ".annotation-base.json")
}
fn pending(base: &Path) -> PathBuf {
    suffix(base, ".annotation-base.pending.json")
}
fn parent(path: &Path) -> &Path {
    path.parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."))
}
fn output_state(paths: &Json) -> Result<Json> {
    let converted = selected(paths, "converted")?;
    let base = selected(paths, "base")?;
    Ok(object([
        ("converted", identity(&converted, true)?),
        (
            "converted_index",
            identity(&suffix(&converted, ".tbi"), true)?,
        ),
        ("base", identity(&base, true)?),
        ("base_index", identity(&suffix(&base, ".tbi"), true)?),
    ]))
}
pub(super) fn counts(path: &Path) -> Result<Json> {
    let output = run(
        &[
            "bcftools".into(),
            "index".into(),
            "--stats".into(),
            path.display().to_string(),
        ],
        false,
    )?;
    let rows = output
        .lines()
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .collect::<Vec<_>>();
    require(
        !rows.is_empty()
            && rows.iter().all(|r| {
                r.len() == 3 && !r[2].is_empty() && r[2].chars().all(|c| c.is_ascii_digit())
            }),
        format!(
            "Missing complete per-contig VCF index counts: {}",
            path.display()
        ),
    )?;
    let mut result = Json::Object(Vec::new());
    for row in rows {
        let n = row[2].parse::<i128>().map_err(|_| {
            super::error(
                "UnsupportedMetadata",
                "VCF index count exceeds supported integer range",
            )
        })?;
        if n != 0 {
            result.set(row[0], Json::Int(n));
        }
    }
    Ok(result)
}
fn moved(actual: &Json, expected: &Json) -> Result<bool> {
    if matches!(actual, Json::Null) || matches!(expected, Json::Null) {
        return Ok(false);
    }
    for key in ["inode", "size", "mtime_ns"] {
        if !equal(get(actual, key)?, get(expected, key)?) {
            return Ok(false);
        }
    }
    Ok(true)
}
fn operation(source: &Path, target: &Path) -> Result<Json> {
    Ok(object([
        ("source", pstr(source)),
        ("target", pstr(target)),
        ("new", identity(source, false)?),
        ("prior", identity(target, true)?),
    ]))
}
pub(super) fn validate_runtime(
    value: &Json,
    paths: &Json,
    science: &Json,
    base_sha: &Json,
) -> Result<()> {
    require(
        value.as_object().is_some()
            && equal(
                optional(value, "vcf")?,
                &identity(&selected(paths, "canonical")?, false)?,
            )
            && equal(
                optional(value, "index")?,
                &identity(&suffix(&selected(paths, "canonical")?, ".tbi"), false)?,
            )
            && equal(optional(value, "annotation_base_sha256")?, base_sha)
            && equal(
                optional(value, "science_sha256")?,
                &string(json_digest(science)?),
            ),
        "Closed native runtime is not bound to this exact complete annotation and scientific profile",
    )?;
    let empty = Json::Object(Vec::new());
    let proof = value.get("closed_native_proof").unwrap_or(&empty);
    require(
        truth(optional(proof, "path")?)
            && equal(
                optional(proof, "sha256")?,
                &string(digest(&path(get(proof, "path")?)?)?),
            ),
        "Closed native runtime proof is missing or changed",
    )
}
impl Context<'_> {
    pub(super) fn settings(&self) -> Result<profile::Settings> {
        profile::settings(self.config, self.source, self.base)
    }
    pub(super) fn capture(&self) -> Result<Json> {
        let settings = self.settings()?;
        let canonical = selected(&settings.paths, "canonical")?;
        Ok(object([
            ("science", profile::science(&settings)?),
            ("prior", output_state(&settings.paths)?),
            (
                "runtime_prior",
                Json::Array(vec![
                    identity(&canonical, true)?,
                    identity(&suffix(&canonical, ".tbi"), true)?,
                ]),
            ),
        ]))
    }
    pub(super) fn check(&self, allow_pending: bool) -> Result<Json> {
        self.check_fields(allow_pending, &[])
    }
    pub(super) fn check_fields(&self, allow_pending: bool, extra: &[&str]) -> Result<Json> {
        let settings = self.settings()?;
        let base = selected(&settings.paths, "base")?;
        require(
            allow_pending || !exists(&pending(&base))?,
            "Complete the pending base publication first",
        )?;
        let encoded = read_text(&receipt(&base))?;
        let schema = compare_metadata(&encoded, &object([("schema", string(SCHEMA))]))?;
        require(
            equal(optional(&schema, "schema")?, &string(SCHEMA)),
            "Complete annotation producer or scientific dependency changed",
        )?;
        let science = profile::science(&settings)?;
        let record = compare_metadata(
            &encoded,
            &object([("science", science.clone()), ("base_sha256", string(""))]),
        )?;
        require(
            profile::same_science(optional(&record, "science")?, &science),
            "Complete annotation producer or scientific dependency changed",
        )?;
        let outputs = output_state(&settings.paths)?;
        let recorded_outputs = compare_metadata(&encoded, &object([("outputs", outputs.clone())]))?;
        require(
            equal(optional(&recorded_outputs, "outputs")?, &outputs),
            "Complete annotation files changed or publication is incomplete",
        )?;
        let header = header::check(
            &base,
            get(&settings.config, "assembly")?,
            get(&settings.config, "vep_installed_version")?,
        )?;
        require(
            equal(
                get(&header, "vep_version")?,
                get(get(&record, "science")?, "vep_version")?,
            ),
            "Complete annotation VEP patch differs",
        )?;
        require(
            optional(&record, "base_sha256")?.as_str().is_some_and(|s| {
                s.len() == 64
                    && s.bytes()
                        .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
            }),
            "Complete annotation has no published content identity",
        )?;
        // The comparison-only projection is never published or exported.
        // After equality succeeds the selected scientific values are decoded
        // losslessly, preserving int/float spelling, negative zero and order.
        let mut fields = vec!["schema", "science", "base_sha256"];
        fields.extend_from_slice(extra);
        crate::json_text::parse_metadata_fields(&encoded, &fields)
            .map_err(|e| error("UnsupportedMetadata", e))
    }
    fn commit(&self) -> Result<()> {
        let settings = self.settings()?;
        let paths = &settings.paths;
        let base = selected(paths, "base")?;
        let marker = pending(&base);
        require(
            is_file(&marker)? && !symlink(&marker)? && metadata(&marker)?.len() < 4 * 1024 * 1024,
            "Missing or invalid base publication marker",
        )?;
        let value = read_json(&marker)?;
        require(
            equal(optional(&value, "schema")?, &string(SCHEMA))
                && profile::same_science(
                    optional(&value, "science")?,
                    &profile::science(&settings)?,
                ),
            "Pending annotation scientific inputs changed",
        )?;
        require(
            equal(optional(&value, "targets")?, paths),
            "Pending publication paths changed",
        )?;
        if !matches!(optional(&value, "native_runtime")?, Json::Null) {
            validate_runtime(
                get(&value, "native_runtime")?,
                paths,
                get(&value, "science")?,
                get(&value, "base_sha256")?,
            )?;
            let converted = selected(paths, "converted")?;
            require(
                equal(
                    get(&value, "preserved_converted")?,
                    &Json::Array(vec![
                        identity(&converted, false)?,
                        identity(&suffix(&converted, ".tbi"), false)?,
                    ]),
                ),
                "Preserved native converted files changed during import",
            )?;
        }
        let root = resolve(&path(get(&settings.config, "cache_root")?)?)?;
        let mut locations = Vec::new();
        for op in list(get(&value, "operations")?)? {
            let source = path(get(op, "source")?)?;
            let target = path(get(op, "target")?)?;
            require(
                resolve(&source)?.starts_with(&root)
                    && resolve(&target)?.starts_with(&root)
                    && !symlink(&source)?
                    && !symlink(&target)?,
                "Pending annotation path escapes managed storage",
            )?;
            if exists(&source)? {
                require(
                    equal(&identity(&source, false)?, get(op, "new")?)
                        && equal(&identity(&target, true)?, get(op, "prior")?),
                    "Staged annotation or previous target changed",
                )?;
                locations.push("stage");
            } else {
                require(
                    moved(&identity(&target, true)?, get(op, "new")?)?,
                    "Published annotation changed",
                )?;
                locations.push("final");
            }
        }
        require(
            locations.windows(2).all(|p| p[0] <= p[1]),
            "Pending publication is not an ordered prefix",
        )?;
        for (location, op) in locations.iter().zip(list(get(&value, "operations")?)?) {
            if *location == "stage" {
                rename(&path(get(op, "source")?)?, &path(get(op, "target")?)?)?;
            }
        }
        if equal(get(&value, "kind")?, &string("base")) {
            atomic_json(
                &receipt(&base),
                &object([
                    ("schema", string(SCHEMA)),
                    ("science", get(&value, "science")?.clone()),
                    ("outputs", output_state(paths)?),
                    ("base_sha256", get(&value, "base_sha256")?.clone()),
                    ("counts", get(&value, "counts")?.clone()),
                    ("lineage", optional(&value, "lineage")?.clone()),
                    (
                        "native_runtime",
                        optional(&value, "native_runtime")?.clone(),
                    ),
                    (
                        "producer_paths",
                        optional(&value, "producer_paths")?.clone(),
                    ),
                ]),
            )?;
        }
        unlink(&marker)
    }
    pub(super) fn recover(&self) -> Result<()> {
        let settings = self.settings()?;
        let base = selected(&settings.paths, "base")?;
        if exists(&pending(&base))? {
            let _lock = lock(parent(&base))?;
            self.commit()?;
        }
        Ok(())
    }
    pub(super) fn publish(
        &self,
        converted: &Path,
        annotation: &Path,
        expected: &Json,
        lineage: Json,
    ) -> Result<()> {
        let settings = self.settings()?;
        let paths = &settings.paths;
        require(
            equal(&self.capture()?, expected),
            "Annotation inputs or previous outputs changed during generation",
        )?;
        let root = resolve(&path(get(&settings.config, "cache_root")?)?)?;
        require(
            converted != annotation
                && resolve(converted)?.starts_with(&root)
                && resolve(annotation)?.starts_with(&root),
            "Scientific stages must remain inside cache_root",
        )?;
        let counts = counts(converted)?;
        require(
            equal(&counts, &self::counts(annotation)?),
            "VEP output lost or added input records, including unsupported AF chromosomes",
        )?;
        let header = header::check(
            annotation,
            get(&settings.config, "assembly")?,
            get(&settings.config, "vep_installed_version")?,
        )?;
        require(
            equal(
                get(&header, "vep_version")?,
                get(get(expected, "science")?, "vep_version")?,
            ),
            "Generated VEP patch differs from actual producer",
        )?;
        for (_, p) in members(paths)? {
            mkdir(parent(&path(p)?))?;
        }
        let base = selected(paths, "base")?;
        let canonical = selected(paths, "canonical")?;
        let _lock = lock(parent(&base))?;
        require(
            !exists(&pending(&base))?,
            "Recover the pending base publication before rebuilding",
        )?;
        require(
            equal(&self.capture()?, expected),
            "Annotation inputs changed before publication",
        )?;
        let stage = stage(parent(&canonical), ".alpha-base-runtime-")?;
        let runtime = stage.join(canonical.file_name().unwrap());
        let sha = copy_verified(annotation, &runtime, None)?;
        copy_verified(&suffix(annotation, ".tbi"), &suffix(&runtime, ".tbi"), None)?;
        let mut operations = Vec::new();
        for (source, target) in [
            (converted, selected(paths, "converted")?),
            (annotation, base.clone()),
            (runtime.as_path(), canonical),
        ] {
            for ext in ["", ".tbi"] {
                let source = suffix(source, ext);
                let target = suffix(&target, ext);
                require(
                    source != target,
                    "Publication needs an independent prepared source",
                )?;
                operations.push(operation(&source, &target)?);
            }
        }
        let producer_paths = object([
            ("source_root", string(absolute(self.source)?)),
            (
                "converter",
                pstr(&self.source.join("scripts/alphmis_tsv2vcf.py")),
            ),
            (
                "basic_vep_annotation",
                pstr(&self.source.join("scripts/install_utils.sh")),
            ),
            (
                "full_basic_procedure_sha256",
                string(crate::sha256::sha256_hex(
                    profile::procedure(
                        &self.source.join("scripts/install_utils.sh"),
                        "basic_vep_annotation",
                    )?
                    .as_bytes(),
                )),
            ),
        ]);
        atomic_json(
            &pending(&base),
            &object([
                ("schema", string(SCHEMA)),
                ("kind", string("base")),
                ("science", get(expected, "science")?.clone()),
                ("targets", paths.clone()),
                ("operations", Json::Array(operations)),
                ("base_sha256", string(sha)),
                ("counts", counts),
                ("lineage", lineage),
                ("producer_paths", producer_paths),
            ]),
        )?;
        self.commit()
    }
    pub(super) fn restore(&self) -> Result<()> {
        let record = self.check(false)?;
        let settings = self.settings()?;
        let paths = &settings.paths;
        let canonical = selected(paths, "canonical")?;
        let base = selected(paths, "base")?;
        mkdir(parent(&canonical))?;
        let _lock = lock(parent(&base))?;
        let stage = stage(parent(&canonical), ".alpha-base-runtime-")?;
        let runtime = stage.join(canonical.file_name().unwrap());
        copy_verified(&base, &runtime, Some(get(&record, "base_sha256")?))?;
        copy_verified(&suffix(&base, ".tbi"), &suffix(&runtime, ".tbi"), None)?;
        let operations = vec![
            operation(&runtime, &canonical)?,
            operation(&suffix(&runtime, ".tbi"), &suffix(&canonical, ".tbi"))?,
        ];
        atomic_json(
            &pending(&base),
            &object([
                ("schema", string(SCHEMA)),
                ("kind", string("runtime")),
                ("science", get(&record, "science")?.clone()),
                ("targets", paths.clone()),
                ("operations", Json::Array(operations)),
            ]),
        )?;
        self.commit()
    }
    pub(super) fn adopt(&self, lineage_file: &Path) -> Result<()> {
        let lineage = read_json(lineage_file)?;
        let expected = self.capture()?;
        require(
            equal(optional(&lineage, "schema")?, &string(SCHEMA))
                && profile::same_science(
                    optional(&lineage, "science")?,
                    get(&expected, "science")?,
                )
                && truth(optional(&lineage, "provenance")?),
            "Native complete-base lineage does not match the actual scientific producer",
        )?;
        let paths = self.settings()?.paths;
        let base = selected(&paths, "base")?;
        mkdir(parent(&base))?;
        let stage = stage(parent(&base), ".alpha-base-adopt-")?;
        let native_runtime = optional(&lineage, "native_runtime")?;
        if !matches!(native_runtime, Json::Null) {
            validate_runtime(
                native_runtime,
                &paths,
                get(&expected, "science")?,
                get(get(&lineage, "annotation")?, "sha256")?,
            )?;
            let converted = path(get(get(&lineage, "converted")?, "path")?)?;
            require(
                resolve(&converted)? == resolve(&selected(&paths, "converted")?)?,
                "Preserved native converted target differs",
            )?;
            let before = Json::Array(vec![
                identity(&converted, false)?,
                identity(&suffix(&converted, ".tbi"), false)?,
            ]);
            require(
                equal(
                    &string(digest(&converted)?),
                    get(get(&lineage, "converted")?, "sha256")?,
                ) && equal(
                    &string(digest(&suffix(&converted, ".tbi"))?),
                    get(get(&lineage, "converted")?, "index_sha256")?,
                ) && equal(&counts(&converted)?, get(&lineage, "counts")?),
                "Preserved native converted identity/counts differ",
            )?;
            let annotation = stage.join("annotation.vcf.gz");
            let source = get(&lineage, "annotation")?;
            let source_path = path(get(source, "path")?)?;
            let sha = copy_verified(&source_path, &annotation, Some(get(source, "sha256")?))?;
            copy_verified(
                &suffix(&source_path, ".tbi"),
                &suffix(&annotation, ".tbi"),
                Some(get(source, "index_sha256")?),
            )?;
            require(
                equal(&counts(&annotation)?, get(&lineage, "counts")?),
                "Native complete-base chromosome counts differ",
            )?;
            let science = get(&expected, "science")?;
            let header = header::check(
                &annotation,
                get(science, "assembly")?,
                &string(
                    text(get(science, "vep_version")?)
                        .split('.')
                        .next()
                        .unwrap(),
                ),
            )?;
            require(
                equal(get(&header, "vep_version")?, get(science, "vep_version")?),
                "Native complete-base VEP patch differs",
            )?;
            let _lock = lock(parent(&base))?;
            require(
                equal(&self.capture()?, &expected)
                    && equal(
                        &before,
                        &Json::Array(vec![
                            identity(&converted, false)?,
                            identity(&suffix(&converted, ".tbi"), false)?,
                        ]),
                    ),
                "Native import inputs changed",
            )?;
            validate_runtime(native_runtime, &paths, science, &string(&sha))?;
            require(
                !exists(&pending(&base))?,
                "Recover the pending base import first",
            )?;
            let operations = vec![
                operation(&annotation, &base)?,
                operation(&suffix(&annotation, ".tbi"), &suffix(&base, ".tbi"))?,
            ];
            atomic_json(
                &pending(&base),
                &object([
                    ("schema", string(SCHEMA)),
                    ("kind", string("base")),
                    ("science", science.clone()),
                    ("targets", paths),
                    ("operations", Json::Array(operations)),
                    ("base_sha256", string(sha)),
                    ("counts", get(&lineage, "counts")?.clone()),
                    ("native_runtime", native_runtime.clone()),
                    ("preserved_converted", before),
                    (
                        "producer_paths",
                        object([("source_root", string(absolute(self.source)?))]),
                    ),
                    (
                        "lineage",
                        object([
                            ("path", string(absolute(lineage_file)?)),
                            ("sha256", string(digest(lineage_file)?)),
                        ]),
                    ),
                ]),
            )?;
            return self.commit();
        }
        let mut selected = Vec::new();
        for key in ["converted", "annotation"] {
            let source = get(&lineage, key)?;
            let target = stage.join(format!("{key}.vcf.gz"));
            let source_path = path(get(source, "path")?)?;
            copy_verified(&source_path, &target, Some(get(source, "sha256")?))?;
            copy_verified(
                &suffix(&source_path, ".tbi"),
                &suffix(&target, ".tbi"),
                Some(get(source, "index_sha256")?),
            )?;
            require(
                equal(&counts(&target)?, get(&lineage, "counts")?),
                "Native complete-base chromosome counts differ",
            )?;
            selected.push(target);
        }
        self.publish(
            &selected[0],
            &selected[1],
            &expected,
            object([
                ("path", string(absolute(lineage_file)?)),
                ("sha256", string(digest(lineage_file)?)),
            ]),
        )
    }
}
