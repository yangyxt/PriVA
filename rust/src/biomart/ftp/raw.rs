//! Original MySQL dump decoding and bounded publisher SQL schema admission.
use super::*;
use crate::python_io::{GzipFailure, PythonGzipReader};
use std::collections::HashSet;
use std::io;
pub(super) type Schema = Vec<(String, Vec<(String, String)>)>;
pub(super) fn gzip(path: &Path) -> Result<PythonGzipReader<File>> {
    Ok(PythonGzipReader::new(
        File::open(path).map_err(|e| os_error(e, path, None))?,
    ))
}
pub(super) fn read_error(cause: io::Error, path: &Path) -> Box<dyn std::error::Error> {
    match GzipFailure::of(&cause) {
        Some(GzipFailure::EndOfFile) => error("EOFError", cause.to_string()),
        Some(GzipFailure::Bad(_)) => error("BadGzipFile", cause.to_string()),
        Some(GzipFailure::Zlib(_)) => error("error", cause.to_string()),
        None => os_error(cause, path, None),
    }
}
pub(super) fn line_error(
    cause: Box<dyn std::error::Error>,
    path: &Path,
) -> Box<dyn std::error::Error> {
    match cause.downcast::<io::Error>() {
        Ok(cause) => read_error(*cause, path),
        Err(cause) => cause,
    }
}
pub(super) fn ascii(bytes: &[u8]) -> Result<String> {
    if let Some(index) = bytes.iter().position(|b| !b.is_ascii()) {
        return Err(error(
            "UnicodeDecodeError",
            format!(
                "'ascii' codec can't decode byte 0x{:02x} in position {index}: ordinal not in range(128)",
                bytes[index]
            ),
        ));
    }
    Ok(std::str::from_utf8(bytes).unwrap().to_owned())
}
pub(super) fn mysql_value(value: &[u8]) -> Result<Option<Vec<u8>>> {
    if value == b"\\N" {
        return Ok(None);
    }
    let mut result = Vec::with_capacity(value.len());
    let mut index = 0;
    while index < value.len() {
        let mut current = value[index];
        if current == b'\\' {
            index += 1;
            require(index < value.len(), "Truncated MySQL field escape")?;
            current = match value[index] {
                b'0' => 0,
                b'b' => 8,
                b'n' => 10,
                b'r' => 13,
                b't' => 9,
                b'Z' => 26,
                b'\\' => b'\\',
                other => other,
            };
        }
        result.push(current);
        index += 1;
    }
    Ok(Some(result))
}
fn column(line: &str) -> Option<(String, String)> {
    let rest = line.strip_prefix("  `")?;
    let end = rest.find('`')?;
    if end == 0 {
        return None;
    }
    let name = &rest[..end];
    let rest = &rest[end + 1..];
    let mut whitespace = Vec::new();
    for (offset, c) in rest.char_indices() {
        if !crate::python_text::is_python_space(c) {
            break;
        }
        whitespace.push(offset + c.len_utf8());
    }
    // Greedy \s+ backtracks if needed so .+ can consume a final space.
    for start in whitespace.into_iter().rev() {
        let tail = &rest[start..];
        let stop = tail.find('\n').unwrap_or(tail.len());
        if stop != 0 {
            return Some((name.to_owned(), tail[..stop].to_owned()));
        }
    }
    None
}
pub(super) fn read_schema(path: &Path) -> Result<Schema> {
    let mut wanted: HashSet<String> = projections()
        .iter()
        .map(|(_, table, _)| table.clone())
        .collect();
    wanted.extend([
        "meta_conf__dataset__main".into(),
        "meta_conf__xml__dm".into(),
    ]);
    let mut schemas: Schema = Vec::new();
    let mut selected: Option<usize> = None;
    let mut total = 0usize;
    let mut lines = crate::python_lines::Lines::latin1(gzip(path)?);
    loop {
        let line = lines
            .readline(usize::MAX)
            .map_err(|e| line_error(e, path))?;
        if line.is_empty() {
            break;
        }
        total += line.chars().count();
        require(total <= 512_000_000, "Schema exceeds bounded size")?;
        if line.starts_with("CREATE TABLE ") {
            let name = line
                .strip_prefix("CREATE TABLE `")
                .and_then(|tail| tail.find('`').filter(|n| *n > 0).map(|end| &tail[..end]))
                .ok_or_else(|| {
                    error(
                        "AttributeError",
                        "'NoneType' object has no attribute 'group'",
                    )
                })?;
            selected = if wanted.contains(name) {
                require(
                    !schemas.iter().any(|(key, _)| key == name),
                    format!("Duplicate schema definition: {name}"),
                )?;
                schemas.push((name.to_owned(), Vec::new()));
                Some(schemas.len() - 1)
            } else {
                None
            };
        } else if let Some(index) = selected {
            if line.starts_with("  `") {
                schemas[index].1.push(
                    column(&line).ok_or_else(|| {
                        error("ValueError", "Cannot parse publisher column schema")
                    })?,
                );
            } else if line.starts_with(')') {
                require(
                    line.contains("CHARSET=latin1"),
                    format!("Unexpected table character set: {}", schemas[index].0),
                )?;
                selected = None;
            }
        }
    }
    require(
        schemas
            .iter()
            .map(|(key, _)| key.clone())
            .collect::<HashSet<_>>()
            == wanted,
        "Required publisher schema tables missing",
    )?;
    for (_, table, fields) in projections() {
        let columns = &schemas.iter().find(|(name, _)| name == table).unwrap().1;
        let names: HashSet<&str> = columns.iter().map(|(name, _)| name.as_str()).collect();
        require(
            names.len() == columns.len(),
            format!("Repeated source column: {table}"),
        )?;
        require(
            fields.iter().all(|field| names.contains(field.as_str())),
            format!("Required projected fields absent: {table}"),
        )?;
        for (name, declaration) in columns {
            if fields.contains(name) && name.ends_with("_key") {
                require(
                    declaration.starts_with("int("),
                    format!("Unexpected join key type: {table}.{name}"),
                )?;
            }
        }
    }
    Ok(schemas)
}
pub(super) fn schema_json(schema: &Schema) -> Json {
    Json::Object(
        schema
            .iter()
            .map(|(name, columns)| {
                (
                    name.clone(),
                    Json::Array(
                        columns
                            .iter()
                            .map(|(name, kind)| Json::Array(vec![string(name), string(kind)]))
                            .collect(),
                    ),
                )
            })
            .collect(),
    )
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mysql_null_empty_and_all_escaped_bytes_are_distinct() {
        assert_eq!(mysql_value(b"\\N").unwrap(), None);
        assert_eq!(mysql_value(b"").unwrap(), Some(vec![]));
        assert_eq!(mysql_value(b"\\\\N").unwrap(), Some(b"\\N".to_vec()));
        assert!(mysql_value(b"bad\\").is_err());
        for byte in 0..=255u8 {
            let expected = match byte {
                b'0' => 0,
                b'b' => 8,
                b'n' => 10,
                b'r' => 13,
                b't' => 9,
                b'Z' => 26,
                b'\\' => b'\\',
                other => other,
            };
            if byte != b'N' {
                assert_eq!(mysql_value(&[b'\\', byte]).unwrap(), Some(vec![expected]));
            }
        }
    }
    #[test]
    fn schema_column_whitespace_matches_greedy_regex_capture() {
        assert_eq!(
            column("  `x` int(10),\n"),
            Some(("x".into(), "int(10),".into()))
        );
        assert_eq!(column("  `x`  \n"), Some(("x".into(), " ".into())));
        assert!(column("  `x` \n").is_none());
        assert_eq!(
            column("  `x`\u{a0}varchar(2),\n"),
            Some(("x".into(), "varchar(2),".into()))
        );
    }
}
