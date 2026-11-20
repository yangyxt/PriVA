//! The publisher's five projections, 10,000-row inserts, original SQLite
//! relationship queries and ordered TSV export. No alternative join engine.
use super::*;
use sqlite::{Connection, Value};
use std::io::{BufRead, BufReader, BufWriter, Write};

fn bump(items: &mut Vec<(String, Json)>, name: &str) {
    if let Some((_, Json::Int(count))) = items.iter_mut().find(|(key, _)| key == name) {
        *count += 1;
    } else {
        items.push((name.to_owned(), Json::Int(1)));
    }
}
fn integer(bytes: &[u8]) -> Result<Value> {
    // Admission already required ASCII digits. int() has its own digit limit;
    // SQLite's signed-64-bit limit is deliberately deferred until batch bind.
    if bytes.len() > 4300 {
        return Err(error(
            "ValueError",
            format!(
                "Exceeds the limit (4300 digits) for integer string conversion: value has {} digits; use sys.set_int_max_str_digits() to increase the limit",
                bytes.len()
            ),
        ));
    }
    let digits = std::str::from_utf8(bytes).unwrap().trim_start_matches('0');
    if digits.is_empty() {
        return Ok(Value::Int(0));
    }
    Ok(digits
        .parse::<i64>()
        .map(Value::Int)
        .unwrap_or(Value::OversizeInteger))
}
fn create(db: &Connection, alias: &str, columns: &[String]) -> Result<()> {
    let declarations = columns
        .iter()
        .map(|name| {
            format!(
                "\"{name}\" {}",
                if name.ends_with("_key") {
                    "INTEGER"
                } else {
                    "TEXT"
                }
            )
        })
        .collect::<Vec<_>>()
        .join(", ");
    db.execute(&format!(
        "CREATE TABLE \"{alias}\" (__row INTEGER PRIMARY KEY, {declarations})"
    ))
}
pub(super) fn load(
    db: &Connection,
    raw: &Path,
    alias: &str,
    table: &str,
    columns: &[String],
    schema: &[(String, String)],
) -> Result<Json> {
    let names = schema.iter().map(|(name, _)| name).collect::<Vec<_>>();
    let indices = columns
        .iter()
        .map(|name| {
            names
                .iter()
                .position(|other| *other == name)
                .ok_or_else(|| error("ValueError", format!("{} is not in list", repr(name))))
        })
        .collect::<Result<Vec<_>>>()?;
    create(db, alias, columns)?;
    let mut batch = Vec::new();
    let mut row_number = 0i64;
    let mut nulls = Vec::new();
    let mut blanks = Vec::new();
    let sql = format!(
        "INSERT INTO \"{alias}\" VALUES ({})",
        vec!["?"; columns.len() + 1].join(",")
    );
    let path = raw.join(format!("{table}.txt.gz"));
    let mut source = BufReader::with_capacity(8192, raw::gzip(&path)?);
    let mut line = Vec::new();
    loop {
        line.clear();
        if source
            .read_until(b'\n', &mut line)
            .map_err(|e| raw::read_error(e, &path))?
            == 0
        {
            break;
        }
        row_number += 1;
        while line.last() == Some(&b'\n') {
            line.pop();
        }
        let cells = line.split(|c| *c == b'\t').collect::<Vec<_>>();
        require(
            cells.len() == names.len(),
            format!("Publisher dump column count differs at {table}:{row_number}"),
        )?;
        let mut projected = vec![Value::Int(row_number)];
        for (name, index) in columns.iter().zip(&indices) {
            let value = raw::mysql_value(cells[*index])?;
            let item = match value {
                None => {
                    bump(&mut nulls, name);
                    Value::Null
                }
                Some(bytes) if name.ends_with("_key") => {
                    require(
                        !bytes.is_empty() && bytes.iter().all(u8::is_ascii_digit),
                        format!("Noninteger key: {table}:{row_number}"),
                    )?;
                    integer(&bytes)?
                }
                Some(bytes) => {
                    let value = raw::ascii(&bytes)?;
                    require(
                        !value
                            .bytes()
                            .any(|c| matches!(c, b'\t' | b'\n' | b'\r' | 0)),
                        format!("Unsafe TSV field: {table}:{row_number}:{name}"),
                    )?;
                    if value.is_empty() {
                        bump(&mut blanks, name);
                    }
                    if matches!(name.as_str(), "stable_id_1023" | "stable_id_1066") {
                        let prefix = if name == "stable_id_1023" {
                            "ENSG"
                        } else {
                            "ENST"
                        };
                        require(
                            value.strip_prefix(prefix).is_some_and(|s| {
                                !s.is_empty() && s.bytes().all(|c| c.is_ascii_digit())
                            }),
                            format!("Malformed stable ID: {table}:{row_number}"),
                        )?;
                    }
                    Value::Text(value)
                }
            };
            projected.push(item);
        }
        batch.push(projected);
        if batch.len() >= 10000 {
            db.executemany(&sql, &batch)?;
            batch.clear();
        }
    }
    if !batch.is_empty() {
        db.executemany(&sql, &batch)?;
    }
    require(row_number > 0, format!("Empty human table: {table}"))?;
    for name in columns {
        if name.ends_with("_key") {
            db.execute(&format!(
                "CREATE INDEX \"{alias}_{name}\" ON \"{alias}\" (\"{name}\")"
            ))?;
        }
    }
    db.commit()?;
    println!("IMPORTED {alias} {row_number}");
    std::io::stdout().flush()?;
    Ok(object([
        ("publisher_table", string(table)),
        ("rows", Json::Int(row_number as i128)),
        (
            "selected_columns",
            Json::Array(columns.iter().map(string).collect()),
        ),
        (
            "source_columns",
            Json::Array(names.into_iter().map(string).collect()),
        ),
        ("sql_null_counts", Json::Object(nulls)),
        ("empty_string_counts", Json::Object(blanks)),
    ]))
}
fn queries() -> &'static Json {
    static QUERIES: OnceLock<Json> = OnceLock::new();
    QUERIES.get_or_init(|| {
        json_text::parse(include_str!("queries.json"), false).expect("unchanged source SQL")
    })
}
pub(super) fn audit(db: &Connection) -> Result<Json> {
    let mut report = Vec::new();
    let read = |key: &str, report: &mut Vec<(String, Json)>| -> Result<i64> {
        let value = db.scalar(text(get(queries(), key)))?;
        report.push((key.to_owned(), value.json()));
        match value {
            Value::Int(n) => Ok(n),
            _ => Err(error("TypeError", "SQL audit count is not an integer")),
        }
    };
    let transcripts = read("transcript_rows", &mut report)?;
    read("main_rows", &mut report)?;
    let distinct = read("distinct_transcript_keys", &mut report)?;
    require(
        transcripts == distinct && transcripts > 0,
        "Transcript main key is not unique",
    )?;
    require(
        db.scalar(text(get(queries(), "invalid_transcript_rows")))? == Value::Int(0),
        "Transcript main lacks stable IDs or positive keys",
    )?;
    require(
        read("transcripts_missing_from_translation_main", &mut report)? == 0,
        "Translation main omits transcripts; no candidate is promoted",
    )?;
    require(
        read(
            "translation_main_missing_or_changed_transcripts",
            &mut report,
        )? == 0,
        "Translation main references absent or inconsistent transcripts",
    )?;
    for (alias, key) in [
        ("hgnc", "gene_id_1020_key"),
        ("swissprot", "translation_id_1068_key"),
        ("refseq", "transcript_id_1064_key"),
    ] {
        let count=db.scalar(&format!("SELECT COUNT(*) FROM \"{alias}\" d WHERE d.\"{key}\">0 AND NOT EXISTS (SELECT 1 FROM main m WHERE m.\"{key}\"=d.\"{key}\")"))?;
        report.push((format!("{alias}_orphan_positive_key_rows"), count.json()));
        require(
            count == Value::Int(0),
            format!("Orphan publisher dimension keys: {alias}"),
        )?;
    }
    let untranslated = read("no_translation_main_rows", &mut report)?;
    let noncoding = read("noncoding_no_translation_transcripts", &mut report)?;
    require(
        untranslated > 0 && noncoding > 0,
        "No untranslated/noncoding transcript witness",
    )?;
    require(
        read("untranslated_with_nonblank_swissprot", &mut report)? == 0,
        "Untranslated transcript acquired SwissProt mapping",
    )?;
    read("expected_output_rows", &mut report)?;
    require(
        read("untranslated_output_rows", &mut report)? >= untranslated,
        "Left joins lost untranslated transcript rows",
    )?;
    Ok(Json::Object(report))
}
fn encode_ascii(text: &str) -> Result<&[u8]> {
    let chars = text.chars().collect::<Vec<_>>();
    if let Some(start) = chars.iter().position(|c| !c.is_ascii()) {
        let end = chars[start..]
            .iter()
            .position(char::is_ascii)
            .map_or(chars.len(), |n| start + n);
        let description = if end == start + 1 {
            let n = chars[start] as u32;
            let code = if n <= 255 {
                format!("\\x{n:02x}")
            } else if n <= 65535 {
                format!("\\u{n:04x}")
            } else {
                format!("\\U{n:08x}")
            };
            format!("character '{code}' in position {start}")
        } else {
            format!("characters in position {start}-{}", end - 1)
        };
        return Err(error(
            "UnicodeEncodeError",
            format!("'ascii' codec can't encode {description}: ordinal not in range(128)"),
        ));
    }
    Ok(text.as_bytes())
}
fn stream_error(cause: std::io::Error) -> Box<dyn std::error::Error> {
    let number = cause.raw_os_error().unwrap_or(libc::EIO);
    let category = if number == libc::EPIPE {
        "BrokenPipeError"
    } else {
        "OSError"
    };
    let description = cause
        .to_string()
        .split(" (os error ")
        .next()
        .unwrap_or("I/O error")
        .to_owned();
    error(category, format!("[Errno {number}] {description}"))
}
pub(super) fn export(db: &Connection, output: &Path, audit: &Json) -> Result<Json> {
    let name = name(output);
    require(
        !name.is_empty(),
        format!("{} has an empty name", repr(&output.display().to_string())),
    )?;
    let part = output.with_file_name(format!("{name}.part"));
    require(
        !exists(output)? && !is_link(output)?,
        "Preserve existing candidate output",
    )?;
    let mut rows = 0i128;
    let mut blanks = Vec::new();
    let file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&part)
        .map_err(|e| os_error(e, &part, None))?;
    let mut handle = BufWriter::with_capacity(8192, file);
    let headers = array(get(get(facts(), "ftp_specification"), "header"));
    let result = (|| -> Result<()> {
        let header = headers.iter().map(text).collect::<Vec<_>>().join("\t") + "\n";
        handle.write_all(header.as_bytes()).map_err(stream_error)?;
        let mut query = db.prepare(text(get(get(facts(), "ftp_specification"), "query_sql")))?;
        while let Some(values) = query.next()? {
            rows += 1;
            for (i, value) in values.iter().enumerate() {
                if matches!(value, Value::Null) || matches!(value,Value::Text(s) if s.is_empty()) {
                    bump(&mut blanks, text(&headers[i]));
                }
            }
            let fields = values
                .iter()
                .enumerate()
                .map(|(i, v)| match v {
                    Value::Null => Ok(""),
                    Value::Text(s) => Ok(s.as_str()),
                    other => Err(error(
                        "TypeError",
                        format!(
                            "sequence item {i}: expected str instance, {} found",
                            match other {
                                Value::Int(_) => "int",
                                Value::Float(_) => "float",
                                Value::Blob(_) => "bytes",
                                _ => "object",
                            }
                        ),
                    )),
                })
                .collect::<Result<Vec<_>>>()?;
            let line = fields.join("\t") + "\n";
            handle
                .write_all(encode_ascii(&line)?)
                .map_err(stream_error)?;
        }
        Ok(())
    })();
    // A with-block closes and flushes even when iteration/encoding fails. A
    // close failure replaces the active exception; never publish a partial.
    handle.flush().map_err(stream_error)?;
    drop(handle);
    result?;
    match audit {
        Json::Null => return Err(error("TypeError", "'NoneType' object is not subscriptable")),
        Json::Array(_) => {
            return Err(error(
                "TypeError",
                "list indices must be integers or slices, not str",
            ));
        }
        Json::Str(_) => {
            return Err(error(
                "TypeError",
                "string indices must be integers, not 'str'",
            ));
        }
        Json::Bool(_) => return Err(error("TypeError", "'bool' object is not subscriptable")),
        Json::Int(_) => return Err(error("TypeError", "'int' object is not subscriptable")),
        Json::Float(_) => return Err(error("TypeError", "'float' object is not subscriptable")),
        Json::Object(_) => {}
    }
    let expected = audit
        .get("expected_output_rows")
        .ok_or_else(|| error("KeyError", "'expected_output_rows'"))?;
    require(
        equal(&Json::Int(rows), expected),
        "Export row count differs from independent dimension product",
    )?;
    fs::rename(&part, output).map_err(|e| os_error(e, &part, Some(output)))?;
    Ok(object([
        ("rows", Json::Int(rows)),
        ("blank_field_counts", Json::Object(blanks)),
        ("bytes", Json::Int(metadata(output)?.len() as i128)),
        ("sha256", string(sha256(output)?)),
    ]))
}
