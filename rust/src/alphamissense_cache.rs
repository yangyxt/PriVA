//! Native `scripts/alphamissense_cache.py`: `install`, `validate` and
//! `check-annotation` for the selected original AlphaMissense prescore.
//!
//! The pinned publisher objects (size, MD5, generation), the bounded schema
//! probe, the install lock, adoption of a supplied or symlinked source, the
//! resumed download, retained replaced artifacts, the rebuilt tabix index and
//! its first-record probe, and the ready receipt follow the script. `bgzip -t`,
//! `tabix` and `curl` remain the external tools the script ran. Printed lines,
//! progress lines and exit statuses follow the script; a failure that the
//! script reported as `AlphaMissense installation failed: message` exits 1
//! here as `priva: CATEGORY: message`. `check-annotation` is the native check
//! that `alpha-annotation-base` already uses. Differences: the subcommand comes
//! first, and temporary indexing folders carry six random characters.
use crate::builder_error::{self, RuntimeError};
use crate::cli::{Result, invalid};
use crate::installer_config_paths::files::{
    exists, is_file, is_link, metadata, mkdir, native, os_error, rename, resolve, stamp,
};
use crate::installer_prerequisites::paths::{join, pathlib};
use crate::json_text::{self, DumpsOptions, Json, py_equal};
use crate::python_args::{parse_long_option_occurrences, require_options};
use crate::python_io::{GzipFailure, PythonGzipReader};
use crate::python_path::{absolute, repr};
use crate::python_subprocess::signal_repr;
use std::collections::HashMap;
use std::ffi::{CString, OsStr};
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufRead, BufReader, Write};
use std::os::fd::AsFd;
use std::os::unix::ffi::OsStrExt;
use std::os::unix::fs::MetadataExt;
use std::os::unix::process::ExitStatusExt;
use std::path::Path;
use std::process::{Command, Stdio};

const HELP: &str = "usage: priva alphamissense-cache {install,validate,check-annotation}\n\
       install --cache-dir CACHE_DIR --assembly ASSEMBLY [--source SOURCE]\n\
       validate --input INPUT --assembly ASSEMBLY\n\
       check-annotation --input INPUT --assembly ASSEMBLY --vep-major VEP_MAJOR\n\
Install the selected original AlphaMissense table into the managed cache; the output and\n\
exit statuses of scripts/alphamissense_cache.py. The subcommand comes first.";
const COMMANDS: [&str; 3] = ["install", "validate", "check-annotation"];
const HEADER: &str = "#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class";

/// An original publisher object, confirmed with the public Google Storage
/// metadata API on 2026-10-18.
struct Source {
    bytes: u64,
    md5: &'static str,
    generation: &'static str,
}

fn source(assembly: &str) -> Source {
    if assembly == "hg19" {
        Source {
            bytes: 622293310,
            md5: "3f14fba08c60b09a90a2ae14b7802e4c",
            generation: "1691073495414540",
        }
    } else {
        Source {
            bytes: 642961469,
            md5: "9fd167735f16a1b87da6eb3e4c25fcb5",
            generation: "1691073413649109",
        }
    }
}

fn publisher(assembly: &str) -> Json {
    let pinned = source(assembly);
    Json::Object(vec![
        ("bytes".into(), Json::Int(pinned.bytes as i128)),
        ("md5".into(), Json::Str(pinned.md5.into())),
        ("generation".into(), Json::Str(pinned.generation.into())),
    ])
}

fn error(category: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    builder_error::category(category, invalid(message))
}
fn require(condition: bool, message: impl Into<String>) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(error("ValueError", message))
    }
}
/// `str(error)` without the category prefix.
fn message(failure: &(dyn std::error::Error + 'static)) -> String {
    failure
        .downcast_ref::<RuntimeError>()
        .map_or_else(|| failure.to_string(), RuntimeError::message)
}
/// The failures the script's `except (ValueError, OSError, SubprocessError)`
/// clauses absorb; EOFError and zlib errors pass through.
fn absorbed(failure: &(dyn std::error::Error + 'static)) -> bool {
    !failure
        .downcast_ref::<RuntimeError>()
        .is_some_and(|runtime| matches!(runtime.category_name(), "EOFError" | "error" | "TypeError"))
}
fn progress(line: &str) {
    eprintln!("{line}");
}

fn assembly_name(value: &str) -> Result<&'static str> {
    match value {
        "hg19" | "GRCh37" => Ok("hg19"),
        "hg38" | "GRCh38" => Ok("hg38"),
        _ => Err(error(
            "ValueError",
            "AlphaMissense requires hg19/GRCh37 or hg38/GRCh38",
        )),
    }
}

/// `file_identity(path)`: local change detection without device numbers.
fn file_identity(path: &str) -> Result<Json> {
    let value = metadata(Path::new(path))?;
    Ok(Json::Object(vec![
        ("inode".into(), Json::Int(value.ino() as i128)),
        ("bytes".into(), Json::Int(value.len() as i128)),
        (
            "mtime_ns".into(),
            Json::Int(value.mtime() as i128 * 1_000_000_000 + value.mtime_nsec() as i128),
        ),
        (
            "ctime_ns".into(),
            Json::Int(value.ctime() as i128 * 1_000_000_000 + value.ctime_nsec() as i128),
        ),
    ]))
}
fn size(path: &str) -> Result<u64> {
    Ok(metadata(Path::new(path))?.len())
}
fn checksum(path: &str) -> Result<String> {
    crate::md5::md5_file(&native(Path::new(path))?)
        .map_err(|cause| os_error(cause, Path::new(path), None))
}

/// `shutil.which(name)`: an executable regular file on PATH.
fn which(name: &str) -> bool {
    let path = std::env::var_os("PATH").unwrap_or_else(|| "/bin:/usr/bin".into());
    std::env::split_paths(&path).any(|folder| {
        let candidate = folder.join(name);
        let executable = CString::new(candidate.as_os_str().as_bytes())
            .is_ok_and(|text| unsafe { libc::access(text.as_ptr(), libc::X_OK) } == 0);
        candidate.exists() && executable && !candidate.is_dir()
    })
}

/// `subprocess.run(command, check=True)` with the given standard output;
/// `capture` returns it as text, as `check_output(..., text=True)` does.
fn run(arguments: &[&str], output: Stdio, capture: bool) -> Result<String> {
    let mut command = Command::new(arguments[0]);
    command.args(&arguments[1..]).stdout(if capture {
        Stdio::piped()
    } else {
        output
    });
    let child = command
        .spawn()
        .map_err(|cause| os_error(cause, Path::new(arguments[0]), None))?;
    let finished = child
        .wait_with_output()
        .map_err(|cause| os_error(cause, Path::new(arguments[0]), None))?;
    if !finished.status.success() {
        let shown = arguments.iter().map(|value| repr(value)).collect::<Vec<_>>();
        let outcome = match (finished.status.code(), finished.status.signal()) {
            (Some(code), _) => format!("returned non-zero exit status {code}"),
            (None, Some(signal)) => format!("died with {}", signal_repr(signal)),
            (None, None) => format!("ended with {}", finished.status),
        };
        return Err(error(
            "CalledProcessError",
            format!("Command '[{}]' {outcome}.", shown.join(", ")),
        ));
    }
    let text = crate::python_io::decode_utf8(&finished.stdout)
        .map_err(|failure| error("UnicodeDecodeError", failure))?;
    Ok(text.replace("\r\n", "\n").replace('\r', "\n"))
}

fn gzip_error(cause: io::Error) -> Box<dyn std::error::Error> {
    let name = match GzipFailure::of(&cause) {
        Some(GzipFailure::EndOfFile) => "EOFError",
        Some(GzipFailure::Bad(_)) => "BadGzipFile",
        Some(GzipFailure::Zlib(_)) => "error",
        None => "OSError",
    };
    error(name, cause.to_string())
}

/// One line of `gzip.open(path, "rt")`: universal newlines, UTF-8.
fn text_line(reader: &mut impl BufRead) -> Result<Option<String>> {
    let mut line = Vec::new();
    loop {
        let available = reader.fill_buf().map_err(gzip_error)?;
        if available.is_empty() {
            break;
        }
        match available.iter().position(|b| matches!(b, b'\n' | b'\r')) {
            Some(at) => {
                let carriage = available[at] == b'\r';
                line.extend_from_slice(&available[..at]);
                reader.consume(at + 1);
                if carriage && reader.fill_buf().map_err(gzip_error)?.first() == Some(&b'\n') {
                    reader.consume(1);
                }
                line.push(b'\n');
                break;
            }
            None => {
                let length = available.len();
                line.extend_from_slice(available);
                reader.consume(length);
            }
        }
    }
    if line.is_empty() {
        return Ok(None);
    }
    crate::python_io::decode_utf8(&line)
        .map(|text| Some(text.to_owned()))
        .map_err(|failure| error("UnicodeDecodeError", failure))
}

/// `first_record(path, assembly)`: the bounded schema and assembly probe.
fn first_record(path: &str, assembly: &str) -> Result<(String, Vec<String>)> {
    let file = File::open(native(Path::new(path))?)
        .map_err(|cause| os_error(cause, Path::new(path), None))?;
    let mut reader = BufReader::new(PythonGzipReader::new(file));
    let mut header_seen = false;
    let mut number = 0;
    while let Some(line) = text_line(&mut reader)? {
        number += 1;
        require(
            number <= 100,
            "AlphaMissense header/first record was not found within 100 lines",
        )?;
        let line = line.trim_end_matches(['\r', '\n']).to_owned();
        if line.starts_with('#') {
            if line.starts_with("#CHROM") {
                require(line == HEADER, "Unexpected AlphaMissense raw TSV header")?;
                header_seen = true;
            }
            continue;
        }
        let fields: Vec<String> = line.split('\t').map(str::to_owned).collect();
        require(
            header_seen && fields.len() == 10,
            "AlphaMissense raw TSV has no expected ten-column record",
        )?;
        require(
            fields[4] == assembly,
            "AlphaMissense raw TSV assembly differs from requested assembly",
        )?;
        require(
            !fields[1].is_empty()
                && fields[1].bytes().all(|b| b.is_ascii_digit())
                && fields[1].bytes().any(|b| b != b'0'),
            "Invalid AlphaMissense position",
        )?;
        return Ok((line, fields));
    }
    Err(error("ValueError", "AlphaMissense raw TSV contains no records"))
}

/// `validate(path, assembly)`: the complete publisher bytes and BGZF integrity.
fn validate(path: &str, assembly: &str) -> Result<Json> {
    let assembly = assembly_name(assembly)?;
    let path = absolute(Path::new(path))?;
    let resolved = resolve(Path::new(&path))?.display().to_string();
    require(
        !resolved
            .split('/')
            .any(|part| part.to_lowercase().contains(".partial")),
        "An explicit AlphaMissense source must not resolve to a partial file",
    )?;
    require(
        is_file(Path::new(&path))? && size(&path)? > 0,
        format!("Missing/empty AlphaMissense source: {path}"),
    )?;
    let before = file_identity(&path)?;
    let expected = source(assembly);
    require(
        size(&path)? == expected.bytes,
        format!("AlphaMissense publisher size mismatch: {path}"),
    )?;
    require(
        checksum(&path)? == expected.md5,
        format!("AlphaMissense publisher MD5 mismatch: {path}"),
    )?;
    run(&["bgzip", "-t", &path], Stdio::null(), false)?;
    first_record(&path, assembly)?;
    require(
        py_equal(&file_identity(&path)?, &before),
        "AlphaMissense source changed during validation",
    )?;
    Ok(before)
}

/// `check_index(path, assembly)`: the index lists the first record's
/// chromosome and returns that record.
fn check_index(path: &str, assembly: &str) -> Result<Vec<String>> {
    let index = format!("{path}.tbi");
    require(
        is_file(Path::new(&index))? && size(&index)? > 0,
        "Missing/empty AlphaMissense tabix index",
    )?;
    let (record, fields) = first_record(path, assembly)?;
    let contigs: Vec<String> = run(&["tabix", "-l", path], Stdio::inherit(), true)?
        .lines()
        .map(str::to_owned)
        .collect();
    require(
        !contigs.is_empty() && contigs.contains(&fields[0]),
        "AlphaMissense index has no source chromosome",
    )?;
    let region = format!("{}:{}-{}", fields[0], fields[1], fields[1]);
    let query = run(&["tabix", path, &region], Stdio::inherit(), true)?;
    require(
        query.lines().any(|line| line == record),
        "AlphaMissense index does not recover its first source record",
    )?;
    Ok(contigs)
}

/// `preserve(path)`: retain a replaced managed artifact beside it.
fn preserve(path: &str) -> Result<()> {
    if exists(Path::new(path))? || is_link(Path::new(path))? {
        let retained = format!("{path}.invalid-{}", stamp()?);
        rename(Path::new(path), Path::new(&retained))?;
        progress(&format!("Retained previous managed artifact: {retained}"));
    }
    Ok(())
}

/// `shutil.copyfile(source, destination)`: content only.
fn copyfile(source: &str, destination: &str) -> Result<()> {
    let mut input = File::open(native(Path::new(source))?)
        .map_err(|cause| os_error(cause, Path::new(source), None))?;
    let mut output = File::create(native(Path::new(destination))?)
        .map_err(|cause| os_error(cause, Path::new(destination), None))?;
    io::copy(&mut input, &mut output).map_err(|cause| os_error(cause, Path::new(source), None))?;
    Ok(())
}

/// `download(path, assembly)`: resume the immutable object; one clean retry.
fn download(path: &str, assembly: &str) -> Result<String> {
    let expected = source(assembly);
    let url = format!(
        "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_{assembly}.tsv.gz?generation={}",
        expected.generation
    );
    let complete = |path: &str| -> Result<bool> {
        Ok(size(path)? == expected.bytes && checksum(path)? == expected.md5)
    };
    for _ in 0..2 {
        if exists(Path::new(path))? && size(path)? >= expected.bytes {
            if complete(path)? {
                return Ok(url);
            }
            preserve(path)?;
        }
        progress(&format!(
            "Downloading selected AlphaMissense {assembly} to {path}"
        ));
        let stderr = io::stderr()
            .as_fd()
            .try_clone_to_owned()
            .map_err(|cause| os_error(cause, Path::new("curl"), None))?;
        run(
            &[
                "curl",
                "--fail",
                "--location",
                "--retry",
                "3",
                "--connect-timeout",
                "60",
                "--continue-at",
                "-",
                "--output",
                path,
                &url,
            ],
            Stdio::from(stderr),
            false,
        )?;
        if is_file(Path::new(path))? && complete(path)? {
            return Ok(url);
        }
        preserve(path)?;
    }
    Err(error(
        "ValueError",
        "AlphaMissense download failed publisher byte verification after one clean retry",
    ))
}

/// `receipt_valid(receipt, path, assembly)`: any failure is false.
fn receipt_valid(receipt: &str, path: &str, assembly: &str) -> bool {
    let Ok(text) = crate::installer_config_paths::files::read_text(Path::new(receipt)) else {
        return false;
    };
    let Ok(data @ Json::Object(_)) = json_text::parse(&text, false) else {
        return false;
    };
    let field = |key: &str| data.get(key).cloned().unwrap_or(Json::Null);
    let (Ok(file), Ok(index)) = (file_identity(path), file_identity(&format!("{path}.tbi"))) else {
        return false;
    };
    py_equal(&field("schema_version"), &Json::Int(1))
        && py_equal(&field("assembly"), &Json::Str(assembly.into()))
        && py_equal(&field("publisher"), &publisher(assembly))
        && py_equal(&field("file_identity"), &file)
        && py_equal(&field("index_identity"), &index)
}

/// `tempfile.TemporaryDirectory(prefix=".indexing-", dir=folder)`
struct Indexing(String);

impl Indexing {
    fn new(folder: &str) -> Result<Self> {
        let template = format!("{folder}/.indexing-XXXXXX");
        let mut raw = CString::new(template.clone())
            .map_err(|_| error("ValueError", "embedded null byte"))?
            .into_bytes_with_nul();
        if unsafe { libc::mkdtemp(raw.as_mut_ptr().cast()) }.is_null() {
            return Err(os_error(
                io::Error::last_os_error(),
                Path::new(&template),
                None,
            ));
        }
        let made = OsStr::from_bytes(&raw[..raw.len() - 1]);
        Ok(Self(made.to_string_lossy().into_owned()))
    }
}

impl Drop for Indexing {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// `install(cache_dir=..., assembly=..., source=...)`
fn install(cache_dir: &str, assembly: &str, supplied: Option<&str>) -> Result<String> {
    let assembly = assembly_name(assembly)?;
    let cache_dir = absolute(Path::new(cache_dir))?;
    require(
        !cache_dir
            .split('/')
            .any(|part| part.to_lowercase().contains(".partial") || part == ".."),
        "AlphaMissense cache directory must be complete and canonical",
    )?;
    let folder = join(&cache_dir, "AlphaMissense");
    mkdir(Path::new(&folder))?;
    let destination = join(&folder, &format!("AlphaMissense_{assembly}.tsv.gz"));
    let index = format!("{destination}.tbi");
    let receipt = join(&folder, &format!("AlphaMissense_{assembly}.ready.json"));
    let supplied = supplied
        .map(|value| absolute(Path::new(value)))
        .transpose()?;
    let external = supplied.as_ref().filter(|value| **value != destination);
    for tool in ["bgzip", "tabix"] {
        require(
            which(tool),
            format!("Missing required AlphaMissense mature tool: {tool}"),
        )?;
    }
    let lock = join(&folder, &format!(".AlphaMissense_{assembly}.install.lock"));
    let lock_file = OpenOptions::new()
        .append(true)
        .create(true)
        .open(native(Path::new(&lock))?)
        .map_err(|cause| os_error(cause, Path::new(&lock), None))?;
    let _guard = crate::file_lock::ExclusiveLock::acquire(lock_file, crate::file_lock::Wait::Blocking)
        .map_err(|cause| os_error(cause, Path::new(&lock), None))?;
    let source_identity = external.map(|path| validate(path, assembly)).transpose()?;
    // A completed, unchanged receipt avoids rereading the large stream; the
    // current data and index still pass a real indexed record probe.
    if !is_link(Path::new(&destination))? && receipt_valid(&receipt, &destination, assembly) {
        match check_index(&destination, assembly) {
            Ok(_) => {
                progress(&format!("Reusing verified managed AlphaMissense {assembly}"));
                return Ok(destination);
            }
            Err(failure) if absorbed(failure.as_ref()) => {}
            Err(failure) => return Err(failure),
        }
    }
    let mut existing_valid = false;
    if is_file(Path::new(&destination))? && !is_link(Path::new(&destination))? {
        match validate(&destination, assembly) {
            Ok(_) => existing_valid = true,
            Err(failure) if absorbed(failure.as_ref()) => progress(&format!(
                "Recovering managed AlphaMissense cache: {}",
                message(failure.as_ref())
            )),
            Err(failure) => return Err(failure),
        }
    }
    let mut original_source = external.cloned();
    let source_url =
        format!("https://storage.googleapis.com/dm_alphamissense/AlphaMissense_{assembly}.tsv.gz");
    let contigs = {
        let temporary = Indexing::new(&folder)?;
        let staged = join(&temporary.0, &format!("AlphaMissense_{assembly}.tsv.gz"));
        let partial = format!("{destination}.partial");
        if existing_valid {
            let target = resolve(Path::new(&destination))?;
            std::os::unix::fs::symlink(&target, &staged)
                .map_err(|cause| os_error(cause, Path::new(&staged), None))?;
        } else {
            if let Some(external) = external {
                preserve(&partial)?;
                copyfile(external, &partial)?;
                let now = file_identity(external)?;
                require(
                    source_identity
                        .as_ref()
                        .is_some_and(|before| py_equal(&now, before)),
                    "Supplied AlphaMissense source changed while copying",
                )?;
            } else if is_link(Path::new(&destination))? && is_file(Path::new(&destination))? {
                // A prior compatibility view is an external source, not a
                // bundled resource. Adopt its verified bytes as a real file.
                match validate(&destination, assembly) {
                    Err(failure) if absorbed(failure.as_ref()) => {
                        // Recover the managed view without changing its old
                        // external target, even when that target is invalid.
                        download(&partial, assembly)?;
                    }
                    Err(failure) => return Err(failure),
                    Ok(_) => {
                        original_source =
                            Some(resolve(Path::new(&destination))?.display().to_string());
                        preserve(&partial)?;
                        copyfile(&destination, &partial)?;
                    }
                }
            } else {
                require(
                    which("curl"),
                    "Missing curl for automatic AlphaMissense acquisition",
                )?;
                download(&partial, assembly)?;
            }
            let expected = source(assembly);
            require(
                size(&partial)? == expected.bytes && checksum(&partial)? == expected.md5,
                "Staged AlphaMissense source failed publisher byte verification",
            )?;
            // Validation rejects partial input names; the staging alias is
            // created only after exact complete publisher bytes are present.
            copyfile(&partial, &staged)?;
            validate(&staged, assembly)?;
        }
        run(
            &["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", "-S", "1", &staged],
            Stdio::inherit(),
            false,
        )?;
        let contigs = check_index(&staged, assembly)?;
        if !existing_valid {
            preserve(&destination)?;
            rename(Path::new(&staged), Path::new(&destination))?;
            fs::remove_file(native(Path::new(&partial))?)
                .map_err(|cause| os_error(cause, Path::new(&partial), None))?;
        }
        preserve(&index)?;
        rename(Path::new(&format!("{staged}.tbi")), Path::new(&index))?;
        contigs
    };
    let pinned = source(assembly);
    let data = Json::Object(vec![
        ("schema_version".into(), Json::Int(1)),
        ("resource".into(), Json::Str("AlphaMissense".into())),
        ("assembly".into(), Json::Str(assembly.into())),
        ("publisher".into(), publisher(assembly)),
        ("source_url".into(), Json::Str(source_url)),
        ("source_generation".into(), Json::Str(pinned.generation.into())),
        (
            "supplied_source".into(),
            original_source.map_or(Json::Null, Json::Str),
        ),
        ("path".into(), Json::Str(destination.clone())),
        ("file_identity".into(), file_identity(&destination)?),
        ("index_identity".into(), file_identity(&index)?),
        (
            "contigs".into(),
            Json::Array(contigs.into_iter().map(Json::Str).collect()),
        ),
        ("scientific_acceptance".into(), Json::Bool(false)),
        (
            "scope".into(),
            Json::Str(
                "Exact original publisher prescore bytes and working tabix index; no VEP/statistics."
                    .into(),
            ),
        ),
    ]);
    let body = json_text::dumps_with_options(
        &data,
        DumpsOptions {
            indent: Some(2),
            ..DumpsOptions::default()
        },
    )
    .map_err(|failure| error("ValueError", failure))?;
    let (mut file, staged_receipt) =
        crate::installer_config_paths::files::temporary(Path::new(&folder), ".ready-")?;
    file.write_all(format!("{body}\n").as_bytes())
        .map_err(|cause| os_error(cause, &staged_receipt, None))?;
    drop(file);
    rename(&staged_receipt, Path::new(&receipt))?;
    Ok(destination)
}

fn options(command: &str) -> (&'static [&'static str], &'static [&'static str]) {
    match command {
        "install" => (&["cache-dir", "assembly", "source"], &["cache-dir", "assembly"]),
        "validate" => (&["input", "assembly"], &["input", "assembly"]),
        _ => (&["input", "assembly", "vep-major"], &["input", "assembly", "vep-major"]),
    }
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(command) = args.first() else {
        return Err(invalid(format!(
            "the following arguments are required: command\n{HELP}"
        )));
    };
    if matches!(command.as_str(), "-h" | "--help") {
        println!("{HELP}");
        return Ok(());
    }
    if !COMMANDS.contains(&command.as_str()) {
        return Err(invalid(format!(
            "argument command: invalid choice: {} (choose from {})",
            repr(command),
            COMMANDS.map(repr).join(", ")
        )));
    }
    let (names, required) = options(command);
    let Some(occurrences) = parse_long_option_occurrences(&args[1..], names, &[]).map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    let values: HashMap<String, String> = occurrences.into_iter().collect();
    require_options(&values, required).map_err(invalid)?;
    let path = |name: &str| values.get(name).map(|value| pathlib(value));
    let assembly = values.get("assembly").cloned().unwrap_or_default();
    let printed = builder_error::finish((|| -> Result<String> {
        match command.as_str() {
            "install" => install(
                &path("cache-dir").unwrap_or_default(),
                &assembly,
                path("source").as_deref(),
            ),
            "validate" => json_text::dumps(&validate(&path("input").unwrap_or_default(), &assembly)?, false)
                .map_err(|failure| error("ValueError", failure)),
            _ => {
                let result = crate::alpha_annotation_base::header::check(
                    Path::new(&path("input").unwrap_or_default()),
                    &Json::Str(assembly.clone()),
                    &Json::Str(values.get("vep-major").cloned().unwrap_or_default()),
                )?;
                json_text::dumps(&result, false).map_err(|failure| error("ValueError", failure))
            }
        }
    })())?;
    println!("{printed}");
    io::stdout().flush().ok();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(items: &[&str]) -> Vec<String> {
        items.iter().map(|item| (*item).to_owned()).collect()
    }

    #[test]
    fn admission_errors_exit_two_and_checks_exit_one() {
        for bad in [
            args(&[]),
            args(&["fetch"]),
            args(&["install", "--cache-dir", "/c"]),
            args(&["validate", "--input", "/x", "--assembly", "hg19", "--source", "/y"]),
        ] {
            let failure = dispatch(&bad).unwrap_err();
            assert!(!failure.is::<RuntimeError>(), "{failure}");
        }
        let failure = dispatch(&args(&["validate", "--input", "/x", "--assembly", "T2T"])).unwrap_err();
        assert_eq!(
            failure.to_string(),
            "ValueError: AlphaMissense requires hg19/GRCh37 or hg38/GRCh38"
        );
    }

    #[test]
    fn universal_newlines_split_lines() {
        let mut reader = io::Cursor::new(b"a\r\nb\rc\nd".to_vec());
        let mut lines = Vec::new();
        while let Some(line) = text_line(&mut reader).unwrap() {
            lines.push(line);
        }
        assert_eq!(lines, ["a\n", "b\n", "c\n", "d"]);
    }
}
