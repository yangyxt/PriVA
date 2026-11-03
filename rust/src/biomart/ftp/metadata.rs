//! Original offline publisher admission. Keep checksum, decompression, row
//! selection and final XML-view publication in their source execution order.
use super::*;
use crate::python_path::is_file;
use std::io::{BufRead, BufReader};
use std::process::{Command, Stdio};

fn gzip(path: &Path) -> Result<crate::python_io::PythonGzipReader<File>> {
    Ok(crate::python_io::PythonGzipReader::with_zlib_diagnostics(
        File::open(path).map_err(|e| os_error(e, path, None))?,
    ))
}

fn pin_values(pin: &Json) -> Result<&[Json]> {
    let Json::Array(values) = pin else {
        return Err(error(
            "TypeError",
            "raw source pin must be a three-value sequence",
        ));
    };
    if values.len() < 3 {
        return Err(error(
            "ValueError",
            format!(
                "not enough values to unpack (expected 3, got {})",
                values.len()
            ),
        ));
    }
    if values.len() > 3 {
        return Err(error(
            "ValueError",
            "too many values to unpack (expected 3)",
        ));
    }
    Ok(values)
}

/// `int()` is used only to compare sum's first two tokens. Valid integers
/// outside i128 cannot equal the bounded publisher pins, but remain valid.
fn checksum_integer(token: &str) -> Result<Option<i128>> {
    use crate::python_text::{decimal_digit_value, is_python_space};
    let mapped: String = token
        .chars()
        .map(|c| {
            if c.is_ascii() {
                c
            } else if let Some(n) = decimal_digit_value(c) {
                char::from(b'0' + n as u8)
            } else if is_python_space(c) {
                ' '
            } else {
                c
            }
        })
        .collect();
    let body = mapped.trim_matches(|c| matches!(c, ' ' | '\t' | '\n' | '\u{b}' | '\u{c}' | '\r'));
    let digits = body
        .strip_prefix('-')
        .or_else(|| body.strip_prefix('+'))
        .unwrap_or(body);
    let valid = !digits.is_empty()
        && !digits.starts_with('_')
        && !digits.ends_with('_')
        && !digits.contains("__")
        && digits.bytes().all(|c| c.is_ascii_digit() || c == b'_');
    if !valid {
        return Err(error(
            "ValueError",
            format!("invalid literal for int() with base 10: {}", repr(token)),
        ));
    }
    let count = digits.bytes().filter(u8::is_ascii_digit).count();
    if count > 4300 {
        return Err(error(
            "ValueError",
            format!(
                "Exceeds the limit (4300 digits) for integer string conversion: value has {count} digits; use sys.set_int_max_str_digits() to increase the limit"
            ),
        ));
    }
    let negative = body.starts_with('-');
    let digits = digits.replace('_', "");
    let digits = digits.trim_start_matches('0');
    if digits.is_empty() {
        return Ok(Some(0));
    }
    Ok(format!("{}{digits}", if negative { "-" } else { "" })
        .parse()
        .ok())
}

fn check_bsd(path: &Path, expected: &Json) -> Result<()> {
    let args = ["-r", "--", path.to_str().expect("CLI paths are Unicode")];
    let output = Command::new("sum")
        .args(args)
        .stdin(Stdio::inherit())
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .output()
        .map_err(|e| os_error(e, Path::new("sum"), None))?;
    // check_output decodes its captured text before checking the exit status.
    let stdout = crate::python_io::decode_utf8(&output.stdout)
        .map_err(|e| error("UnicodeDecodeError", e))?;
    if !output.status.success() {
        use std::os::unix::process::ExitStatusExt;
        let command = format!("['sum', '-r', '--', {}]", repr(&path.display().to_string()));
        let message = if let Some(code) = output.status.code() {
            format!("Command '{command}' returned non-zero exit status {code}.")
        } else {
            format!(
                "Command '{command}' died with {}.",
                crate::python_subprocess::signal_repr(output.status.signal().unwrap_or(1))
            )
        };
        return Err(error("CalledProcessError", message));
    }
    let observed = stdout
        .split(crate::python_text::is_python_space)
        .filter(|s| !s.is_empty())
        .take(2)
        .map(checksum_integer)
        .collect::<Result<Vec<_>>>()?;
    let matches = match expected {
        Json::Array(values) => {
            values.len() == observed.len()
                && values.iter().zip(&observed).all(|(expected, value)| {
                    value.is_some_and(|value| equal(expected, &Json::Int(value)))
                })
        }
        _ => false,
    };
    require(
        matches,
        format!("Publisher BSD checksum mismatch: {}", path.display()),
    )
}

pub(in crate::biomart) fn verify_raw(path: &Path, pin: &Json) -> Result<Json> {
    let values = pin_values(pin)?;
    require(
        is_file(path)? && equal(&Json::Int(metadata(path)?.len().into()), &values[0]),
        format!("Raw size mismatch: {}", path.display()),
    )?;
    let digest = sha256(path)?;
    if truth(&values[2]) {
        require(
            equal(&string(&digest), &values[2]),
            format!("Pinned raw SHA256 mismatch: {}", path.display()),
        )?;
    }
    let mut result = object([
        ("path", string(path.display().to_string())),
        ("bytes", values[0].clone()),
        ("sha256", string(digest)),
    ]);
    if truth(&values[1]) {
        check_bsd(path, &values[1])?;
        result.set("publisher_bsd_sum_blocks", values[1].clone());
    }
    if name(path).ends_with(".gz") || name(path).ends_with(".gz.part") {
        let mut source = gzip(path)?;
        let mut buffer = vec![0; 1024 * 1024];
        let mut total = 0i128;
        loop {
            let count = source
                .read(&mut buffer)
                .map_err(|e| raw::read_error(e, path))?;
            if count == 0 {
                break;
            }
            total += count as i128;
        }
        result.set("full_gzip_integrity_verified", Json::Bool(true));
        result.set("uncompressed_bytes", Json::Int(total));
    }
    Ok(result)
}

pub(in crate::biomart) fn named_pin(name: &str) -> Result<Json> {
    let pin = array(get(get(facts(), "ftp_specification"), "raw_sources"))
        .iter()
        .find(|v| get(v, "name").as_str() == Some(name))
        .ok_or_else(|| error("KeyError", repr(name)))?;
    Ok(Json::Array(vec![
        get(pin, "bytes").clone(),
        get(pin, "publisher_bsd_sum_blocks").clone(),
        get(pin, "sha256").clone(),
    ]))
}

fn columns(schema: &Json, table: &str) -> Result<Json> {
    let fields = schema
        .get(table)
        .ok_or_else(|| error("KeyError", repr(table)))?;
    let Json::Array(fields) = fields else {
        return Err(error(
            "TypeError",
            "metadata schema columns must be a sequence",
        ));
    };
    let mut names = Vec::new();
    for field in fields {
        let Json::Array(parts) = field else {
            return Err(error("TypeError", "metadata schema entries must be pairs"));
        };
        if parts.len() < 2 {
            return Err(error(
                "ValueError",
                format!(
                    "not enough values to unpack (expected 2, got {})",
                    parts.len()
                ),
            ));
        }
        if parts.len() > 2 {
            return Err(error(
                "ValueError",
                "too many values to unpack (expected 2)",
            ));
        }
        names.push(parts[0].clone());
    }
    Ok(Json::Array(names))
}

fn each_line(path: &Path, mut visit: impl FnMut(&[u8]) -> Result<()>) -> Result<()> {
    let mut source = BufReader::with_capacity(8192, gzip(path)?);
    let mut line = Vec::new();
    loop {
        line.clear();
        if source
            .read_until(b'\n', &mut line)
            .map_err(|e| raw::read_error(e, path))?
            == 0
        {
            break;
        }
        visit(&line)?;
    }
    Ok(())
}
fn strip_lf(line: &[u8]) -> &[u8] {
    let end = line.iter().rposition(|c| *c != b'\n').map_or(0, |i| i + 1);
    &line[..end]
}
fn decompress(blob: &[u8]) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    crate::python_io::PythonGzipReader::with_zlib_diagnostics(blob).read_to_end(&mut bytes).map_err(|e| {
        // gzip.decompress has the same member/header checks but intentionally
        // omits the two CRC values from its own diagnostic.
        if matches!(crate::python_io::GzipFailure::of(&e), Some(crate::python_io::GzipFailure::Bad(text)) if text.starts_with("CRC check failed ")) {
            error("BadGzipFile", "CRC check failed")
        } else { raw::read_error(e, Path::new("")) }
    })?;
    Ok(bytes)
}

pub(in crate::biomart) fn publisher(raw: &Path, schemas: &Json, root: &Path) -> Result<Json> {
    let expected = [
        "dataset_id_key",
        "dataset",
        "display_name",
        "description",
        "type",
        "visible",
        "version",
        "modified",
    ];
    require(
        equal(
            &columns(schemas, "meta_conf__dataset__main")?,
            &Json::Array(expected.into_iter().map(string).collect()),
        ),
        "Dataset metadata column order changed",
    )?;
    let mut matches = Vec::new();
    each_line(&raw.join("meta_conf__dataset__main.txt.gz"), |line| {
        let cell = line
            .split(|c| *c == b'\t')
            .nth(1)
            .ok_or_else(|| error("IndexError", "list index out of range"))?;
        if cell == DATASET.as_bytes() {
            matches.push(
                strip_lf(line)
                    .split(|c| *c == b'\t')
                    .map(<[u8]>::to_vec)
                    .collect::<Vec<_>>(),
            );
        }
        Ok(())
    })?;
    require(
        matches.len() == 1,
        "Human dataset metadata row absent/nonunique",
    )?;
    let row = matches
        .remove(0)
        .iter()
        .map(|cell| {
            let value = raw::mysql_value(cell)?.ok_or_else(|| {
                error(
                    "AttributeError",
                    "'NoneType' object has no attribute 'decode'",
                )
            })?;
            raw::ascii(&value).map(string)
        })
        .collect::<Result<Vec<_>>>()?;
    require(
        equal(
            &columns(schemas, "meta_conf__xml__dm")?,
            &Json::Array(
                ["dataset_id_key", "xml", "compressed_xml", "message_digest"]
                    .into_iter()
                    .map(string)
                    .collect(),
            ),
        ),
        "XML metadata column order changed",
    )?;
    let mut found = Vec::new();
    each_line(&raw.join("meta_conf__xml__dm.txt.gz"), |line| {
        let fields = strip_lf(line).split(|c| *c == b'\t').collect::<Vec<_>>();
        require(fields.len() == 4, "Malformed publisher XML dump row")?;
        if fields[0] == text(&row[0]).as_bytes() {
            found.push(
                fields
                    .into_iter()
                    .map(raw::mysql_value)
                    .collect::<Result<Vec<_>>>()?,
            );
        }
        Ok(())
    })?;
    require(found.len() == 1, "Human XML row absent/nonunique")?;
    let fields = found.remove(0);
    let (Some(blob), Some(compressed)) = (&fields[1], &fields[2]) else {
        return Err(error(
            "ValueError",
            "Expected raw and compressed human XML copies",
        ));
    };
    require(
        decompress(compressed)? == *blob,
        "Raw/compressed human XML copies differ",
    )?;
    let result = super::super::verify::parse_xml(blob, &Json::Array(row))?;
    require(
        get(&result, "xml_sha256").as_str()
            == Some("a5afb502ad97bada17f9fd4cdc3886be18ba8b237f3b896568793eb0cdd74fa7"),
        "Pinned human XML changed",
    )?;
    let view = root.join("hsapiens_gene_ensembl.dataset.xml");
    if exists(&view)? {
        require(
            fs::read(&view).map_err(|e| os_error(e, &view, None))? == *blob,
            "Existing human XML view differs",
        )?;
    } else {
        fs::write(&view, blob).map_err(|e| os_error(e, &view, None))?;
    }
    Ok(result)
}
