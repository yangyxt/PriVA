//! Verify a separately prepared DAS consumer-input equivalence certificate.
//! This procedure never computes equivalence, scans a VCF or modifies a cache.
//! Preserve the original check/publication order and no-overwrite hard link.
use crate::cli::{Result, invalid};
use crate::das_install_state::{self as state, Options};
use crate::json_text::{self, DumpsOptions, Json};
use crate::python_path::{absolute, os_error, path_value, repr, tempfile_parent};
use crate::sha256::Sha256;
use std::ffi::CString;
use std::fs::{self, File, Metadata};
use std::io::{Read, Write};
use std::os::fd::FromRawFd;
use std::os::unix::ffi::OsStrExt;
use std::os::unix::fs::{MetadataExt, OpenOptionsExt};
use std::path::{Path, PathBuf};

const KIND: &str = "das-clinvar-input-continuity-v1";
const PROOF_KIND: &str = "das-clinvar-consumer-equivalence-v1";
const MAX_BYTES: u64 = 16 * 1024 * 1024;
const CHECKER: &[u8] = include_bytes!("../../scripts/das_input_continuity.py");
const HELP: &str = "usage: priva das-input-continuity {record|check} --checker-source FILE --receipt JSON --continuity JSON --profile PROFILE --input FILE --output FILE [--source ROLE=FILE ...] [--bundled ROLE=FILE ...] [--parameter NAME=VALUE ...]\n\
record also requires --old-input FILE --new-input FILE --proof JSON --proof-sha256 HEX.\n\
The checker source is the pinned original Python source asset; it is hashed, never executed.\n\
Existing installer callers and scientific caches are unchanged.";

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
fn object(items: impl IntoIterator<Item = (&'static str, Json)>) -> Json {
    Json::Object(items.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}
fn string(value: impl Into<String>) -> Json {
    Json::Str(value.into())
}
fn kind(value: &Json) -> &'static str {
    match value {
        Json::Null => "NoneType",
        Json::Bool(_) => "bool",
        Json::Int(_) => "int",
        Json::Float(_) => "float",
        Json::Str(_) => "str",
        Json::Array(_) => "list",
        Json::Object(_) => "dict",
    }
}
fn get<'a>(value: &'a Json, name: &str) -> Result<&'a Json> {
    if !matches!(value, Json::Object(_)) {
        return Err(error(
            "AttributeError",
            format!("'{}' object has no attribute 'get'", kind(value)),
        ));
    }
    Ok(value.get(name).unwrap_or(&Json::Null))
}
fn at<'a>(value: &'a Json, name: &str) -> Result<&'a Json> {
    match value {
        Json::Object(_) => value.get(name).ok_or_else(|| error("KeyError", repr(name))),
        Json::Array(_) => Err(error(
            "TypeError",
            "list indices must be integers or slices, not str",
        )),
        Json::Str(_) => Err(error(
            "TypeError",
            "string indices must be integers, not 'str'",
        )),
        _ => Err(error(
            "TypeError",
            format!("'{}' object is not subscriptable", kind(value)),
        )),
    }
}
fn text(value: &Json) -> Result<&str> {
    value.as_str().ok_or_else(|| {
        error(
            "RuntimeError",
            "Unsupported non-string filesystem path in continuity metadata",
        )
    })
}
fn equal(a: &Json, b: &Json) -> bool {
    json_text::py_equal(a, b)
}
fn hash(raw: &[u8]) -> String {
    let mut h = Sha256::new();
    h.update(raw);
    h.finish_hex()
}
fn hex64(value: &Json) -> bool {
    value.as_str().is_some_and(|s| {
        s.len() == 64
            && s.bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    })
}
fn identity(info: &Metadata) -> (u64, u64, u64, i64, i64, i64, i64, u32) {
    (
        info.dev(),
        info.ino(),
        info.len(),
        info.mtime(),
        info.mtime_nsec(),
        info.ctime(),
        info.ctime_nsec(),
        info.mode(),
    )
}
fn metadata(path: &Path) -> Result<Vec<u8>> {
    let before = fs::symlink_metadata(path).map_err(|e| os_error(e, path, None))?;
    require(
        before.is_file() && before.len() > 0 && before.len() <= MAX_BYTES,
        format!(
            "Expected nonempty bounded regular metadata: {}",
            path.display()
        ),
    )?;
    let stream = fs::OpenOptions::new()
        .read(true)
        .custom_flags(libc::O_NOFOLLOW | libc::O_CLOEXEC)
        .open(path)
        .map_err(|e| os_error(e, path, None))?;
    require(
        identity(&stream.metadata().map_err(|e| os_error(e, path, None))?) == identity(&before),
        format!("Metadata changed before reading: {}", path.display()),
    )?;
    let mut bounded = stream.take(MAX_BYTES + 1);
    let mut raw = Vec::new();
    bounded
        .read_to_end(&mut raw)
        .map_err(|e| os_error(e, path, None))?;
    require(
        raw.len() as u64 == before.len()
            && identity(
                &bounded
                    .get_ref()
                    .metadata()
                    .map_err(|e| os_error(e, path, None))?,
            ) == identity(&before),
        format!("Metadata changed during reading: {}", path.display()),
    )?;
    drop(bounded);
    require(
        identity(&fs::symlink_metadata(path).map_err(|e| os_error(e, path, None))?)
            == identity(&before),
        format!("Metadata path changed during reading: {}", path.display()),
    )?;
    Ok(raw)
}
fn one_identity(path: &str) -> Result<Json> {
    let Json::Array(mut rows) = state::identities(&[path.to_owned()])? else {
        unreachable!()
    };
    Ok(rows.remove(0))
}
fn pin(path: &str) -> Result<Json> {
    let spelling = absolute(&path_value(path))?;
    let path = Path::new(&spelling);
    let before = one_identity(&spelling)?;
    let raw = metadata(path)?;
    require(
        equal(&one_identity(&spelling)?, &before),
        format!("Metadata identity changed: {}", path.display()),
    )?;
    Ok(object([
        ("path", string(spelling)),
        ("sha256", string(hash(&raw))),
        ("identity", before),
    ]))
}
fn pinned_bytes(value: &Json) -> Result<Vec<u8>> {
    let valid = value.as_object().is_some_and(|items| {
        items.len() == 3
            && ["path", "sha256", "identity"]
                .iter()
                .all(|name| value.get(name).is_some())
    }) && value
        .get("path")
        .and_then(Json::as_str)
        .is_some_and(|p| Path::new(p).is_absolute())
        && value.get("sha256").is_some_and(hex64);
    require(valid, "Invalid metadata proof pin")?;
    let path = text(at(value, "path")?)?;
    let before = one_identity(path)?;
    let raw = metadata(&path_value(path))?;
    require(
        equal(&before, at(value, "identity")?)
            && equal(at(value, "identity")?, &one_identity(path)?)
            && at(value, "sha256")?.as_str() == Some(hash(&raw).as_str()),
        format!("Pinned metadata changed: {path}"),
    )?;
    Ok(raw)
}
fn parse(raw: &[u8]) -> Result<Json> {
    // The original certificates are UTF-8 JSON. Other JSON byte encodings and
    // lone-surrogate values remain explicit admission limits, never replacement
    // text or silently accepted mismatches. A BOM in UTF-8 is handled by loads.
    if raw.starts_with(&[0xff, 0xfe])
        || raw.starts_with(&[0xfe, 0xff])
        || raw.starts_with(&[0, 0, 0xfe, 0xff])
        || ((raw.len() >= 4 || raw.len() == 2) && (raw[0] == 0 || raw[1] == 0))
    {
        return Err(error(
            "RuntimeError",
            "Unsupported continuity JSON byte encoding; UTF-8 is required",
        ));
    }
    let bytes = raw.strip_prefix(&[0xef, 0xbb, 0xbf]).unwrap_or(raw);
    let source =
        crate::python_io::decode_utf8(bytes).map_err(|e| error("UnicodeDecodeError", e))?;
    json_text::parse(source, false).map_err(|e| {
        if json_text::is_refusal(&e) {
            error("RuntimeError", format!("Unsupported JSON profile: {e}"))
        } else {
            error("JSONDecodeError", e)
        }
    })
}
fn read_pin(value: &Json) -> Result<Json> {
    parse(&pinned_bytes(value)?)
}
fn proof(
    path: &str,
    checksum: &str,
    profile: &str,
    old: &Json,
    new: &Json,
    sources: &Json,
) -> Result<(Json, Vec<Json>)> {
    let binding = pin(path)?;
    require(
        at(&binding, "sha256")?.as_str() == Some(checksum),
        "Explicit consumer-equivalence proof changed",
    )?;
    let value = read_pin(&binding)?;
    require(
        get(&value, "kind")?.as_str() == Some(PROOF_KIND)
            && matches!(get(&value, "schema_version")?, Json::Int(1))
            && get(&value, "status")?.as_str() == Some("completed")
            && get(&value, "profile")?.as_str() == Some(profile),
        "Consumer-equivalence proof is incomplete or belongs to another profile",
    )?;
    require(
        equal(get(&value, "old_input")?, old) && equal(get(&value, "new_input")?, new),
        "Consumer-equivalence proof input identities differ",
    )?;
    let gates = get(&value, "gates")?;
    require(
        matches!(gates, Json::Object(_))
            && [
                "full_record_census_equal",
                "full_global_order_equal",
                "protected_annotation_fields_equal",
                "headers_reviewed",
                "publication_identity_exact",
            ]
            .iter()
            .all(|name| matches!(gates.get(name), Some(Json::Bool(true)))),
        "Consumer-equivalence publication/header/row gates are incomplete",
    )?;
    let empty = Json::Object(Vec::new());
    let consumer = value.get("consumer").unwrap_or(&empty);
    let count = get(consumer, "ordered_rows")?;
    let digest = get(consumer, "old_sha256")?;
    let Json::Array(sources) = sources else {
        unreachable!()
    };
    let analysis: Vec<_> = sources
        .iter()
        .filter(|v| v.get("role").and_then(Json::as_str) == Some("analysis"))
        .map(|v| v.get("sha256").unwrap())
        .collect();
    require(
        get(consumer, "role")?.as_str() == Some("das_and_dependent_pm1")
            && matches!(get(consumer, "exact")?, Json::Bool(true))
            && matches!(get(consumer, "differing_rows")?, Json::Int(0))
            && matches!(count,Json::Int(n) if *n>0)
            && hex64(digest)
            && equal(get(consumer, "new_sha256")?, digest)
            && analysis.len() == 1
            && equal(get(consumer, "analysis_source_sha256")?, analysis[0]),
        "Exact ordered DAS consumer/source binding is missing",
    )?;
    let evidence = get(&value, "evidence")?;
    require(
        evidence.as_object().is_some_and(|_| {
            ["plan", "comparison", "publication"]
                .iter()
                .all(|name| evidence.get(name).is_some())
        }),
        "Underlying completed comparison/publication evidence is missing",
    )?;
    let evidence: Vec<_> = evidence
        .as_object()
        .unwrap()
        .iter()
        .map(|(_, v)| v.clone())
        .collect();
    for item in &evidence {
        pinned_bytes(item)?;
    }
    Ok((binding, evidence))
}
fn exists(path: &Path) -> Result<bool> {
    match fs::metadata(path) {
        Ok(_) => Ok(true),
        Err(e)
            if matches!(
                e.raw_os_error(),
                Some(libc::ENOENT | libc::ENOTDIR | libc::EBADF | libc::ELOOP)
            ) =>
        {
            Ok(false)
        }
        Err(e) => Err(os_error(e, path, None)),
    }
}
fn samefile(left: &Path, right: &Path) -> Result<bool> {
    let a = fs::metadata(left).map_err(|e| os_error(e, left, None))?;
    let b = fs::metadata(right).map_err(|e| os_error(e, right, None))?;
    Ok(a.dev() == b.dev() && a.ino() == b.ino())
}
fn candidate(
    options: &Options,
    old_path: &str,
    new_path: &str,
    proof_path: &str,
    checksum: &str,
) -> Result<Json> {
    let current = state::completion(state::input_state(options)?, &options.outputs)?;
    let original_pin = pin(&options.single["receipt"])?;
    let original = read_pin(&original_pin)?;
    let Json::Array(old_new) = state::identities(&[old_path.to_owned(), new_path.to_owned()])?
    else {
        unreachable!()
    };
    let (old, new) = (&old_new[0], &old_new[1]);
    require(
        !equal(old, new),
        "Continuity requires a distinct selected ClinVar input",
    )?;
    let Json::Array(inputs) = at(&current, "inputs")? else {
        unreachable!()
    };
    let positions: Vec<_> = inputs
        .iter()
        .enumerate()
        .filter(|(_, v)| equal(v, new))
        .map(|(i, _)| i)
        .collect();
    require(
        positions.len() == 1,
        "Exactly one selected new ClinVar input is required",
    )?;
    let ordinal = positions[0];
    let mut prior = current.clone();
    let mut old_inputs = inputs.clone();
    old_inputs[ordinal] = old.clone();
    prior.set("inputs", Json::Array(old_inputs));
    require(
        equal(&original, &prior),
        "Original DAS source, parameters, bundled inputs, other data or outputs differ",
    )?;
    let mut protected = vec![options.single["receipt"].as_str()];
    protected.extend(options.inputs.iter().map(String::as_str));
    protected.extend(options.outputs.iter().map(String::as_str));
    protected.push(proof_path);
    let destination = PathBuf::from(absolute(&path_value(&options.single["continuity"]))?);
    for selected in protected {
        let path = path_value(selected);
        require(
            destination.to_string_lossy() != absolute(&path)?
                && (!exists(&destination)? || !samefile(&destination, &path)?),
            "Continuity receipt aliases a protected file",
        )?;
    }
    let (proof_pin, evidence_pins) = proof(
        proof_path,
        checksum,
        &options.single["profile"],
        old,
        new,
        at(&current, "sources")?,
    )?;
    require(
        equal(
            &state::completion(state::input_state(options)?, &options.outputs)?,
            &current,
        ) && equal(
            &state::identities(&[old_path.to_owned(), new_path.to_owned()])?,
            &Json::Array(old_new.clone()),
        ) && equal(&pin(&options.single["receipt"])?, &original_pin),
        "DAS inputs, outputs or original receipt changed during proof verification",
    )?;
    for item in std::iter::once(&proof_pin).chain(evidence_pins.iter()) {
        pinned_bytes(item)?;
    }
    // Python pins __file__ here. An explicit matching source asset preserves
    // that original certificate field without executing a Python bridge.
    let checker = pin(&options.single["checker-source"])?;
    if at(&checker, "sha256")?.as_str() != Some(hash(CHECKER).as_str()) {
        return Err(error(
            "RuntimeError",
            "Unsupported continuity checker source; use the pinned unmodified source asset",
        ));
    }
    Ok(object([
        ("kind", string(KIND)),
        ("schema_version", Json::Int(1)),
        ("producer_execution", Json::Bool(false)),
        ("checker_sha256", at(&checker, "sha256")?.clone()),
        ("original_receipt", original_pin),
        ("old_input", old.clone()),
        ("new_input", new.clone()),
        ("input_ordinal", Json::Int(ordinal as i128)),
        ("current_state", current),
        ("equivalence_proof", proof_pin),
    ]))
}
fn dumps(value: &Json, indented: bool) -> Result<String> {
    json_text::dumps_with_options(
        value,
        DumpsOptions {
            sort_keys: true,
            indent: indented.then_some(2),
            allow_nan: false,
            ..DumpsOptions::default()
        },
    )
    .map_err(|e| error("ValueError", e))
}
fn record(options: &Options) -> Result<()> {
    let build = || {
        candidate(
            options,
            &options.single["old-input"],
            &options.single["new-input"],
            &options.single["proof"],
            &options.single["proof-sha256"],
        )
    };
    let value = build()?;
    let raw = (dumps(&value, true)? + "\n").into_bytes();
    let path = path_value(&options.single["continuity"]);
    // lexists includes a dangling link. Its metadata check must fail without
    // replacing the link or following its target.
    let lexists = match fs::symlink_metadata(&path) {
        Ok(_) => true,
        Err(e) if matches!(e.raw_os_error(), Some(libc::ENOENT | libc::ENOTDIR)) => false,
        Err(e) => return Err(os_error(e, &path, None)),
    };
    if lexists {
        return require(
            metadata(&path)? == raw,
            "Existing continuity differs; preserved",
        );
    }
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let name = path.file_name().unwrap_or_default().to_string_lossy();
    let template = tempfile_parent(parent)?.join(format!("{name}.XXXXXX"));
    let mut name = CString::new(template.as_os_str().as_bytes())
        .map_err(|_| error("ValueError", "embedded null byte"))?
        .into_bytes_with_nul();
    let fd = unsafe { libc::mkostemp(name.as_mut_ptr().cast(), libc::O_CLOEXEC) };
    let temporary = PathBuf::from(std::ffi::OsStr::from_bytes(&name[..name.len() - 1]));
    if fd < 0 {
        return Err(os_error(std::io::Error::last_os_error(), &temporary, None));
    }
    let result = (|| {
        let mut stream = unsafe { File::from_raw_fd(fd) };
        stream
            .write_all(&raw)
            .map_err(|e| os_error(e, &temporary, None))?;
        stream
            .sync_all()
            .map_err(|e| os_error(e, &temporary, None))?;
        drop(stream);
        require(
            equal(&build()?, &value),
            "DAS continuity changed before publication",
        )?;
        match fs::hard_link(&temporary, &path) {
            Ok(()) => Ok(()),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => require(
                metadata(&path)? == raw,
                "Concurrent continuity differs; preserved",
            ),
            Err(e) => Err(os_error(e, &temporary, Some(&path))),
        }
    })();
    fs::remove_file(&temporary).map_err(|e| os_error(e, &temporary, None))?;
    result
}
fn check(options: &Options) -> Result<()> {
    let value = parse(&metadata(&path_value(&options.single["continuity"]))?)?;
    require(
        get(&value, "kind")?.as_str() == Some(KIND),
        "Unknown DAS continuity receipt",
    )?;
    let expected = candidate(
        options,
        text(at(at(&value, "old_input")?, "path")?)?,
        text(at(at(&value, "new_input")?, "path")?)?,
        text(at(at(&value, "equivalence_proof")?, "path")?)?,
        text(at(at(&value, "equivalence_proof")?, "sha256")?)?,
    )?;
    require(
        dumps(&value, false)? == dumps(&expected, false)?,
        "DAS consumer-input continuity no longer matches",
    )?;
    println!("Reusing original DAS outputs through exact ClinVar consumer-input continuity");
    Ok(())
}
pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(command) = args.first().map(String::as_str) else {
        return Err(invalid(HELP));
    };
    if matches!(command, "-h" | "--help") {
        println!("{HELP}");
        return Ok(());
    }
    if !matches!(command, "record" | "check") {
        return Err(invalid(HELP));
    }
    let mut names = vec![
        "checker-source",
        "receipt",
        "continuity",
        "profile",
        "input",
        "output",
        "source",
        "bundled",
        "parameter",
    ];
    let mut required = vec![
        "checker-source",
        "receipt",
        "continuity",
        "profile",
        "input",
        "output",
    ];
    if command == "record" {
        names.extend(["old-input", "new-input", "proof", "proof-sha256"]);
        required.extend(["old-input", "new-input", "proof", "proof-sha256"]);
    }
    let Some(occurrences) =
        crate::python_args::parse_long_option_occurrences(&args[1..], &names, &[])
            .map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    let mut options = Options::default();
    for (name, value) in occurrences {
        match name.as_str() {
            "input" => options.inputs.push(value.clone()),
            "output" => options.outputs.push(value.clone()),
            "source" => options.sources.push(value.clone()),
            "bundled" => options.bundled.push(value.clone()),
            "parameter" => options.parameters.push(value.clone()),
            _ => (),
        }
        options.single.insert(name, value);
    }
    crate::python_args::require_options(&options.single, &required).map_err(invalid)?;
    crate::builder_error::finish(if command == "record" {
        record(&options)
    } else {
        check(&options)
    })
}
