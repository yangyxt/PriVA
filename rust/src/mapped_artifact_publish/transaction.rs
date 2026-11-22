use super::*;
use crate::python_path::{metadata as stat_path, os_error, repr};
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};
use std::os::unix::fs::MetadataExt;
use std::os::unix::process::ExitStatusExt;
use std::process::Command;

pub(super) fn file_state(path: &Path) -> Result<Json> {
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
    fn state(m: &fs::Metadata) -> (u64, u64, i64, i64, i64, i64, u32, u64, u32) {
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
        state(&before) == state(&after),
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
pub(super) fn load_pending(marker: &Path) -> Result<String> {
    let before = file_state(marker)?;
    require(
        !matches!(before, Json::Null) && matches!(get(&before,"size")?,Json::Int(n) if *n<=65536),
        "Missing/oversized pending publication",
    )?;
    let text = read_text(marker)?;
    let value = Value::validate(&text)?;
    require(
        equal(&file_state(marker)?, &before),
        "Pending publication changed while reading",
    )?;
    require(
        value.exact_object(&[
            "schema",
            "context",
            "target",
            "stage",
            "prior",
            "new",
            "producer_job_id",
        ])? && value.item("schema")?.equal(&string(SCHEMA))?,
        "Unknown pending publication format",
    )?;
    Ok(text)
}
fn moved(actual: &Json, expected: Value<'_>) -> Result<bool> {
    if matches!(actual, Json::Null) || expected.0 == "null" {
        return Ok(false);
    }
    for (key, span) in expected.fields()?.ok_or_else(|| {
        error(
            "AttributeError",
            format!("'{}' object has no attribute 'items'", expected.kind()),
        )
    })? {
        if key != "ctime_ns" && !Value(&expected.0[span]).equal(get(actual, &key)?)? {
            return Ok(false);
        }
    }
    Ok(true)
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
    for name in [
        "scientific_input",
        "native_input_lineage",
        "native_producer",
    ] {
        if let Some(v) = expected.get(name) {
            value.set(name, v.clone());
        }
    }
    Ok(value)
}
fn raw_artifact_record(expected: &Json, state: Value<'_>) -> Result<Json> {
    artifact_record(
        expected,
        &object([
            ("size", state.item("size")?.json()?),
            ("sha256", state.item("sha256")?.json()?),
        ]),
    )
}
impl Context {
    pub(super) fn finish(&self, target: &Path, expected: &Json) -> Result<()> {
        let marker = suffix(target, ".publish.pending.json");
        let encoded = load_pending(&marker)?;
        let value = Value(encoded.trim());
        require(
            value.item("context")?.equal(expected)?
                && value.item("target")?.equal(&pstr(target))?,
            "Pending publication producer/config/mapping/target differs",
        )?;
        let stage_path = path(&value.item("stage")?.json()?)?;
        require(
            stage_path.parent() == target.parent()
                && stage_path.file_name().is_some_and(|n| {
                    n.to_string_lossy()
                        .starts_with(&format!(".mapped-artifact-{}-", self.kind))
                })
                && !symlink(&stage_path)?
                && stage_path.is_dir()
                && stat_path(&stage_path)?.uid() == unsafe { libc::getuid() }
                && stat_path(&stage_path)?.mode() & 0o777 == 0o700,
            "Pending publication stage is not owned beside target",
        )?;
        require(
            value.item("prior")?.keys(&["artifact", "receipt"])?
                && value.item("new")?.keys(&["artifact", "receipt"])?,
            "Pending pair fields differ",
        )?;
        let artifact = stage_path.join(target.file_name().unwrap_or_default());
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
            let new = value.item("new")?.item(key)?;
            require(
                new.exact_object(&[
                    "inode", "size", "mtime_ns", "ctime_ns", "mode", "nlink", "uid", "sha256",
                ])? && new.get("size")?.positive()?,
                "Pending staged file identity is incomplete",
            )?;
            if !matches!(staged_state, Json::Null) {
                require(
                    moved(&staged_state, new)?
                        && value.item("prior")?.item(key)?.equal(&final_state)?,
                    "Pending stage or previous selected artifact changed",
                )?;
                locations.push(false);
            } else {
                require(
                    moved(&final_state, new)?,
                    "Pending published file changed or disappeared",
                )?;
                locations.push(true);
            }
        }
        require(
            locations != [false, true],
            "Pending publication is not an ordered two-file state",
        )?;
        let raw = read_text(if locations[1] {
            &paths[1].2
        } else {
            &paths[1].1
        })?;
        require_equal(
            Value::validate(&raw)?,
            &raw_artifact_record(expected, value.item("new")?.item("artifact")?)?,
            "Pending mapping receipt does not describe the staged artifact",
        )?;
        for (published, (_, staged, final_path)) in locations.iter().zip(&paths) {
            if !published {
                rename(staged, final_path)?;
            }
        }
        self.receipt_check(target)?;
        for (key, _, final_path) in &paths {
            require(
                moved(&file_state(final_path)?, value.item("new")?.item(key)?)?,
                "Published pair changed before commit",
            )?;
        }
        unlink(&marker)?;
        if fs::read_dir(&stage_path)
            .map_err(|e| os_error(e, &stage_path, None))?
            .next()
            .is_none()
        {
            fs::remove_dir(&stage_path).map_err(|e| os_error(e, &stage_path, None))?;
        }
        Ok(())
    }
    pub(super) fn run(&self, expected_identity: &str, arguments: &[String]) -> Result<()> {
        let target = self.target()?;
        mkdir(target.parent().unwrap())?;
        let _lock = lock(target.parent().unwrap())?;
        let expected = self.expected()?;
        require_equal(
            Value::validate(expected_identity)?,
            get(&expected, "mapping_identity")?,
            "Selected mapping changed before generation",
        )?;
        if lexists(&suffix(&target, ".publish.pending.json"))? {
            return self.finish(&target, &self.recovery_context(&target)?);
        }
        if lexists(&target)? || lexists(&suffix(&target, ".nm_to_enst.json"))? {
            self.receipt_check(&target)?;
        }
        let prior = object([
            ("artifact", file_state(&target)?),
            ("receipt", file_state(&suffix(&target, ".nm_to_enst.json"))?),
        ]);
        let mut argv = arguments.to_vec();
        if argv.first().is_some_and(|v| v == "--") {
            argv.remove(0);
        }
        require(
            argv.iter().filter(|v| v.as_str() == "--output").count() == 1,
            "Scientific argv must contain exactly one --output",
        )?;
        let index = argv.iter().position(|v| v == "--output").unwrap() + 1;
        require(
            index < argv.len() && argv[index] == self.artifact.to_string_lossy(),
            "Scientific --output must be the selected configured artifact",
        )?;
        if self.native_producer.is_some() {
            self.native_arguments(&argv)?;
        }
        let stage_path = stage(
            target.parent().unwrap(),
            &format!(".mapped-artifact-{}-", self.kind),
        )?;
        let staged = stage_path.join(target.file_name().unwrap_or_default());
        argv[index] = staged.display().to_string();
        println!(
            "Generating mapped artifact in retained stage: {}",
            stage_path.display()
        );
        std::io::stdout().flush()?;
        let mut command = match &self.native_producer {
            Some(executable) => vec![
                executable.display().to_string(),
                NATIVE_PM1.to_owned(),
                "build".to_owned(),
            ],
            None => vec![self.python.clone(), self.producer.display().to_string()],
        };
        command.extend(argv);
        let status = Command::new(&command[0])
            .args(&command[1..])
            .status()
            .map_err(|e| os_error(e, Path::new(&command[0]), None))?;
        if !status.success() {
            let display = format!(
                "[{}]",
                command
                    .iter()
                    .map(|s| repr(s))
                    .collect::<Vec<_>>()
                    .join(", ")
            );
            return Err(error(
                "CalledProcessError",
                if let Some(code) = status.code() {
                    format!("Command '{display}' returned non-zero exit status {code}.")
                } else {
                    format!(
                        "Command '{display}' died with {}.",
                        crate::python_subprocess::signal_repr(status.signal().unwrap_or(1))
                    )
                },
            ));
        }
        let staged_state = file_state(&staged)?;
        require(
            !matches!(staged_state, Json::Null)
                && matches!(get(&staged_state,"size")?,Json::Int(n) if *n>0),
            if self.native_producer.is_some() {
                "Native producer created no nonempty archive"
            } else {
                "Scientific producer created no nonempty pickle"
            },
        )?;
        require(
            equal(&self.expected()?, &expected),
            "Producer/config/mapping changed during generation",
        )?;
        biomart::atomic_json(
            &suffix(&staged, ".nm_to_enst.json"),
            &artifact_record(&expected, &staged_state)?,
        )?;
        self.receipt_check(&staged)?;
        require(
            equal(
                &object([
                    ("artifact", file_state(&target)?),
                    ("receipt", file_state(&suffix(&target, ".nm_to_enst.json"))?),
                ]),
                &prior,
            ),
            "Selected artifact changed during generation",
        )?;
        let value = publication(&target, &stage_path, expected.clone(), prior, &staged)?;
        let temporary = stage_path.join("publication.json");
        let mut stream = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
            .map_err(|e| os_error(e, &temporary, None))?;
        stream
            .write_all(formatted_json(&value)?.as_bytes())
            .map_err(|e| os_error(e, &temporary, None))?;
        drop(stream);
        require(
            !lexists(&suffix(&target, ".publish.pending.json"))?,
            "Another pending publication exists",
        )?;
        rename(&temporary, &suffix(&target, ".publish.pending.json"))?;
        self.finish(&target, &expected)
    }
    pub(super) fn adopt(&self, source: &Path, lineage: Value<'_>) -> Result<()> {
        let target = self.target()?;
        require(
            resolve(source)? != resolve(&target)?,
            "Native motif adoption must select a new managed target",
        )?;
        require(
            !symlink(source)? && is_file(source)?,
            "Preserved native motif must be an independent regular file",
        )?;
        mkdir(target.parent().unwrap())?;
        let _lock = lock(target.parent().unwrap())?;
        let expected = self.expected()?;
        require_equal(
            lineage.get("scientific_input")?,
            optional(&expected, "scientific_input")?,
            "Preserved native motif base/producer binding differs",
        )?;
        if lexists(&suffix(&target, ".publish.pending.json"))? {
            return self.finish(&target, &expected);
        }
        if lexists(&target)? || lexists(&suffix(&target, ".nm_to_enst.json"))? {
            self.receipt_check(&target)?;
            let encoded = read_text(&suffix(&target, ".nm_to_enst.json"))?;
            let completed = Value::validate(&encoded)?;
            require(
                completed
                    .get("scientific_input")?
                    .equal(optional(&expected, "scientific_input")?)?
                    && completed
                        .get("native_input_lineage")?
                        .equal(optional(&expected, "native_input_lineage")?)?
                    && lineage
                        .get("artifact_sha256")?
                        .equal(&completed.get("artifact_sha256")?.json()?)?,
                "Native adoption must not replace an existing artifact or historical receipt",
            )?;
            return Ok(());
        }
        self.receipt_check(source)?;
        let source_before = file_state(source)?;
        let receipt_before = file_state(&suffix(source, ".nm_to_enst.json"))?;
        let encoded = read_text(&suffix(source, ".nm_to_enst.json"))?;
        let original = Value::validate(&encoded)?;
        for key in ["mapping_identity", "artifact_bytes", "artifact_sha256"] {
            require(
                lineage.get(key)?.equal(&original.get(key)?.json()?)?,
                "Preserved native motif identity differs",
            )?;
        }
        require_equal(
            lineage.get("receipt_sha256")?,
            get(&receipt_before, "sha256")?,
            "Preserved native receipt SHA256 differs",
        )?;
        let stage_path = stage(target.parent().unwrap(), ".mapped-artifact-motifs-")?;
        let staged = stage_path.join(target.file_name().unwrap_or_default());
        let mut input = File::open(source).map_err(|e| os_error(e, source, None))?;
        let mut output = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&staged)
            .map_err(|e| os_error(e, &staged, None))?;
        let mut buffer = vec![0; 8 * 1024 * 1024];
        loop {
            let count = input
                .read(&mut buffer)
                .map_err(|e| os_error(e, source, None))?;
            if count == 0 {
                break;
            }
            output
                .write_all(&buffer[..count])
                .map_err(|e| os_error(e, &staged, None))?;
        }
        drop(output);
        drop(input);
        let staged_state = file_state(&staged)?;
        require(
            equal(
                get(&staged_state, "sha256")?,
                get(&source_before, "sha256")?,
            ) && equal(&file_state(source)?, &source_before)
                && equal(
                    &file_state(&suffix(source, ".nm_to_enst.json"))?,
                    &receipt_before,
                ),
            "Preserved native motif changed while copying",
        )?;
        biomart::atomic_json(
            &suffix(&staged, ".nm_to_enst.json"),
            &artifact_record(&expected, &staged_state)?,
        )?;
        require(
            equal(&self.expected()?, &expected),
            "Native motif inputs changed before publication",
        )?;
        let value = publication(
            &target,
            &stage_path,
            expected.clone(),
            object([("artifact", Json::Null), ("receipt", Json::Null)]),
            &staged,
        )?;
        biomart::atomic_json(&suffix(&target, ".publish.pending.json"), &value)?;
        self.finish(&target, &expected)
    }
}
fn publication(
    target: &Path,
    stage: &Path,
    context: Json,
    prior: Json,
    staged: &Path,
) -> Result<Json> {
    Ok(object([
        ("schema", string(SCHEMA)),
        ("context", context),
        ("target", pstr(target)),
        ("stage", pstr(stage)),
        ("prior", prior),
        (
            "new",
            object([
                ("artifact", file_state(staged)?),
                ("receipt", file_state(&suffix(staged, ".nm_to_enst.json"))?),
            ]),
        ),
        (
            "producer_job_id",
            std::env::var("SLURM_JOB_ID")
                .map(string)
                .unwrap_or(Json::Null),
        ),
    ]))
}
