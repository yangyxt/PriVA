//! Completion metadata for the multi-output DAS build, not scientific validity.
//! Keep the Python procedure's observation and publication order. In particular,
//! bundled resources are checked only for existence/nonzero size, never hashed.
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::python_args::parse_long_option_occurrences;
use crate::python_path::{absolute, is_file, metadata, os_error, path_value, tempfile_parent};
use crate::sha256::Sha256;
use std::collections::HashMap;
use std::ffi::CString;
use std::fs::{self, File};
use std::io::{Read, Write};
use std::os::fd::FromRawFd;
use std::os::unix::ffi::OsStrExt;
use std::os::unix::fs::MetadataExt;
use std::path::{Path, PathBuf};

const POLICY: &str = "data=path,target,inode,size,mtime_ns; node_local_st_dev_excluded; source=role,sha256; bundled=role,exists_nonempty; parameters=ordered";
const HELP: &str = "usage: priva das-install-state {capture|check|record} --profile PROFILE --input FILE [--input FILE ...] [--source ROLE=FILE ...] [--bundled ROLE=FILE ...] [--parameter NAME=VALUE ...]\n\
check and record also require --receipt JSON and --output FILE (repeatable).\n\
record requires --expected-inputs JSON from capture.\n\
This records completed publication only. Existing installer callers are unchanged.";

fn error(category: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    crate::builder_error::category(category, invalid(message))
}
fn object(items: impl IntoIterator<Item = (&'static str, Json)>) -> Json {
    Json::Object(items.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}
fn string(value: impl Into<String>) -> Json {
    Json::Str(value.into())
}
fn time_ns(seconds: i64, nanos: i64) -> i128 {
    seconds as i128 * 1_000_000_000 + nanos as i128
}
pub(crate) fn identities(paths: &[String]) -> Result<Json> {
    let mut result = Vec::new();
    for name in paths {
        let path = PathBuf::from(absolute(&path_value(name))?);
        let info = metadata(&path)?;
        if !info.is_file() || info.len() == 0 {
            return Err(error(
                "ValueError",
                format!("Missing/nonregular/empty DAS file: {}", path.display()),
            ));
        }
        // stat() above already established an existing regular target. The
        // resolved spelling is independent of its lexical (possibly linked) path.
        let target = fs::canonicalize(&path).map_err(|e| os_error(e, &path, None))?;
        result.push(object([
            ("path", string(path.to_string_lossy())),
            ("target", string(target.to_string_lossy())),
            ("inode", Json::Int(info.ino() as i128)),
            ("size", Json::Int(info.len() as i128)),
            (
                "mtime_ns",
                Json::Int(time_ns(info.mtime(), info.mtime_nsec())),
            ),
        ]));
    }
    Ok(Json::Array(result))
}
fn named(values: &[String]) -> Result<Vec<(String, String)>> {
    values
        .iter()
        .map(|value| {
            let pair = value.split_once('=');
            match pair {
                Some((name, selected)) if !name.is_empty() && !selected.is_empty() => {
                    Ok((name.to_owned(), selected.to_owned()))
                }
                _ => Err(error(
                    "ValueError",
                    format!("Use a nonempty logical NAME=VALUE: {value}"),
                )),
            }
        })
        .collect()
}
fn sources(values: &[String]) -> Result<Json> {
    let mut result = Vec::new();
    for (role, selected) in named(values)? {
        let path = path_value(&selected);
        let spelling = vec![path.to_string_lossy().into_owned()];
        let before = identities(&spelling)?;
        let info = metadata(&path)?;
        let before_ctime = time_ns(info.ctime(), info.ctime_nsec());
        let mut digest = Sha256::new();
        let mut stream = File::open(&path).map_err(|e| os_error(e, &path, None))?;
        let mut block = vec![0u8; 1024 * 1024];
        loop {
            let count = stream
                .read(&mut block)
                .map_err(|e| os_error(e, &path, None))?;
            if count == 0 {
                break;
            }
            digest.update(&block[..count]);
        }
        drop(stream);
        if !json_text::py_equal(&identities(&spelling)?, &before) || {
            let info = metadata(&path)?;
            time_ns(info.ctime(), info.ctime_nsec()) != before_ctime
        } {
            return Err(error(
                "ValueError",
                format!("DAS source changed while hashing: {}", path.display()),
            ));
        }
        result.push(object([
            ("role", string(role)),
            ("sha256", string(digest.finish_hex())),
        ]));
    }
    Ok(Json::Array(result))
}
fn bundled(values: &[String]) -> Result<Json> {
    let mut result = Vec::new();
    for (role, selected) in named(values)? {
        let path = path_value(&selected);
        if !is_file(&path)? || metadata(&path)?.len() == 0 {
            return Err(error(
                "ValueError",
                format!(
                    "Missing/nonregular/empty bundled DAS file: {}",
                    path.display()
                ),
            ));
        }
        result.push(object([
            ("role", string(role)),
            ("exists_nonempty", Json::Bool(true)),
        ]));
    }
    Ok(Json::Array(result))
}
#[derive(Default)]
pub(crate) struct Options {
    pub(crate) single: HashMap<String, String>,
    pub(crate) inputs: Vec<String>,
    pub(crate) sources: Vec<String>,
    pub(crate) bundled: Vec<String>,
    pub(crate) parameters: Vec<String>,
    pub(crate) outputs: Vec<String>,
}
pub(crate) fn input_state(options: &Options) -> Result<Json> {
    // Keep these evaluations sequential: a later malformed argument must not
    // hide the earlier input/source failure selected by the original procedure.
    Ok(object([
        ("profile", string(options.single["profile"].clone())),
        ("identity_policy", string(POLICY)),
        ("inputs", identities(&options.inputs)?),
        ("sources", sources(&options.sources)?),
        ("bundled", bundled(&options.bundled)?),
        (
            "parameters",
            Json::Array(
                named(&options.parameters)?
                    .into_iter()
                    .map(|(name, value)| Json::Array(vec![string(name), string(value)]))
                    .collect(),
            ),
        ),
    ]))
}
pub(crate) fn completion(current: Json, outputs: &[String]) -> Result<Json> {
    let Json::Object(items) = current else {
        unreachable!()
    };
    let mut values = vec![("version".to_owned(), Json::Int(3))];
    values.extend(items);
    values.push(("outputs".to_owned(), identities(outputs)?));
    Ok(Json::Object(values))
}
fn parse(text: &str) -> Result<Json> {
    json_text::parse(text, false).map_err(|message| {
        if json_text::is_refusal(&message) {
            error(
                "RuntimeError",
                format!("Unsupported JSON profile: {message}"),
            )
        } else {
            error("JSONDecodeError", message)
        }
    })
}
fn read_receipt(path: &Path) -> Result<Json> {
    let bytes = fs::read(path).map_err(|e| os_error(e, path, None))?;
    let text = crate::python_io::decode_utf8(&bytes)
        .map_err(|message| error("UnicodeDecodeError", message))?;
    parse(&text.replace("\r\n", "\n").replace('\r', "\n"))
}
fn publish(path: &Path, value: &Json) -> Result<()> {
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let name = path.file_name().unwrap_or_default().to_string_lossy();
    let template = tempfile_parent(parent)?.join(format!("{name}.XXXXXX"));
    let mut template = CString::new(template.as_os_str().as_bytes())
        .map_err(|_| error("ValueError", "embedded null byte"))?
        .into_bytes_with_nul();
    // mkostemp creates an independent 0600 file with O_EXCL, as mkstemp in
    // Python does. Only its unpredictable temporary basename differs.
    let fd = unsafe { libc::mkostemp(template.as_mut_ptr().cast(), libc::O_CLOEXEC) };
    let temporary = PathBuf::from(std::ffi::OsStr::from_bytes(&template[..template.len() - 1]));
    if fd < 0 {
        return Err(os_error(std::io::Error::last_os_error(), &temporary, None));
    }
    let result = (|| {
        let mut file = unsafe { File::from_raw_fd(fd) };
        let text = json_text::dumps_indented(value).map_err(invalid)? + "\n";
        file.write_all(text.as_bytes())
            .map_err(|e| os_error(e, &temporary, None))?;
        file.sync_all().map_err(|e| os_error(e, &temporary, None))?;
        drop(file);
        fs::rename(&temporary, path).map_err(|e| os_error(e, &temporary, Some(path)))?;
        Ok(())
    })();
    if temporary.exists() {
        fs::remove_file(&temporary).map_err(|e| os_error(e, &temporary, None))?;
    }
    result
}
fn execute(command: &str, options: &Options) -> Result<()> {
    let current = input_state(options)?;
    if command == "capture" {
        println!("{}", json_text::dumps(&current, false).map_err(invalid)?);
        return Ok(());
    }
    if command == "record"
        && !json_text::py_equal(&current, &parse(&options.single["expected-inputs"])?)
    {
        return Err(error(
            "ValueError",
            "DAS inputs changed during generation; no completion recorded",
        ));
    }
    let expected = completion(current, &options.outputs)?;
    let receipt = path_value(&options.single["receipt"]);
    if command == "record" {
        publish(&receipt, &expected)
    } else {
        if !json_text::py_equal(&read_receipt(&receipt)?, &expected) {
            return Err(error(
                "ValueError",
                "DAS completion metadata no longer matches selected files",
            ));
        }
        Ok(())
    }
}
pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(command) = args.first().map(String::as_str) else {
        return Err(invalid(HELP));
    };
    if matches!(command, "--help" | "-h") {
        println!("{HELP}");
        return Ok(());
    }
    if !matches!(command, "capture" | "record" | "check") {
        return Err(invalid(HELP));
    }
    let mut names = vec!["profile", "input", "source", "bundled", "parameter"];
    let mut required = vec!["profile", "input"];
    if command != "capture" {
        names.extend(["receipt", "output"]);
        required.extend(["receipt", "output"]);
    }
    if command == "record" {
        names.push("expected-inputs");
        required.push("expected-inputs");
    }
    let Some(occurrences) =
        parse_long_option_occurrences(&args[1..], &names, &[]).map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    let mut options = Options::default();
    for (name, value) in occurrences {
        let values = match name.as_str() {
            "input" => Some(&mut options.inputs),
            "output" => Some(&mut options.outputs),
            "source" => Some(&mut options.sources),
            "bundled" => Some(&mut options.bundled),
            "parameter" => Some(&mut options.parameters),
            _ => None,
        };
        if let Some(values) = values {
            values.push(value.clone());
        }
        options.single.insert(name, value);
    }
    crate::python_args::require_options(&options.single, &required).map_err(invalid)?;
    crate::builder_error::finish(execute(command, &options))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn named_roles_preserve_order_duplicates_spaces_and_embedded_equals() {
        let values = ["b=2", "a=1=2", "b=3", " = "].map(str::to_owned);
        assert_eq!(
            named(&values).unwrap(),
            vec![
                ("b".into(), "2".into()),
                ("a".into(), "1=2".into()),
                ("b".into(), "3".into()),
                (" ".into(), " ".into())
            ]
        );
        for value in ["", "missing", "=value", "role="] {
            assert!(named(&[value.into()]).is_err());
        }
    }
}
