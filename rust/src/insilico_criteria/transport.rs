//! One bounded, length-framed CSV session for complete source-position phases.
//! Protocol failures are fail-closed. Semantic unsupported responses are
//! explicit and let the host execute that actual mature operation once.
use super::{Cell, Value, calculate};
use std::io::{Read, Write};

const MAX_FRAME: usize = 256 * 1024 * 1024;

pub(super) fn decode_cell(tag: &str, raw: &str) -> Result<Cell, String> {
    match tag {
        "none" if raw.is_empty() => Ok(Cell::None),
        "na" if raw.is_empty() => Ok(Cell::Na),
        "bool" if raw == "0" || raw == "1" => Ok(Cell::Bool(raw == "1")),
        "int" => raw
            .parse()
            .map(Cell::Int)
            .map_err(|_| "invalid integer cell".into()),
        "float" if raw.len() == 16 => u64::from_str_radix(raw, 16)
            .map(|bits| Cell::Float(f64::from_bits(bits)))
            .map_err(|_| "invalid float bits".into()),
        "str" => Ok(Cell::Text(raw.into())),
        _ => Err("invalid typed cell".into()),
    }
}

pub(super) fn encode_cell(cell: &Cell) -> (&'static str, String) {
    match cell {
        Cell::None => ("none", String::new()),
        Cell::Na => ("na", String::new()),
        Cell::Bool(value) => ("bool", u8::from(*value).to_string()),
        Cell::Int(value) => ("int", value.to_string()),
        Cell::Float(value) => ("float", format!("{:016x}", value.to_bits())),
        Cell::Text(value) => ("str", value.clone()),
    }
}

pub(super) fn request(bytes: &[u8]) -> Result<(u64, String, Vec<Value>), String> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .flexible(false)
        .from_reader(bytes);
    let mut rows = reader.records();
    let row = rows
        .next()
        .ok_or("missing request")?
        .map_err(|error| error.to_string())?;
    if row.len() != 6 || &row[0] != "request" || &row[5] != "insilico-columns-v1" {
        return Err("request schema mismatch".into());
    }
    let sequence = row[1].parse::<u64>().map_err(|_| "invalid sequence")?;
    let operation = row[2].to_owned();
    let count = row[3]
        .parse::<usize>()
        .map_err(|_| "invalid operand count")?;
    if count > 32 || !row[4].is_empty() {
        return Err("invalid request count/reserved field".into());
    }
    let mut values: Vec<Value> = Vec::new();
    for row in rows {
        let row = row.map_err(|error| error.to_string())?;
        if row.len() != 6 {
            return Err("argument schema mismatch".into());
        }
        let argument = row[1]
            .parse::<usize>()
            .map_err(|_| "invalid operand ordinal")?;
        match &row[0] {
            "arg" => {
                if argument != values.len()
                    || argument >= count
                    || !row[4].is_empty()
                    || !row[5].is_empty()
                    || !["scalar", "series", "array", "list"].contains(&&row[2])
                {
                    return Err("nonconsecutive argument metadata".into());
                }
                values.push(Value {
                    kind: row[2].into(),
                    dtype: row[3].into(),
                    index: Vec::new(),
                    cells: Vec::new(),
                });
            }
            "index" | "cell" => {
                let value = values.get_mut(argument).ok_or("argument before metadata")?;
                let ordinal = row[2]
                    .parse::<usize>()
                    .map_err(|_| "invalid element ordinal")?;
                if !row[5].is_empty() {
                    return Err("nonempty reserved field".into());
                }
                if &row[0] == "index" {
                    if value.kind != "series" || ordinal != value.index.len() || &row[3] != "label"
                    {
                        return Err("invalid index record".into());
                    }
                    value.index.push(row[4].into());
                } else {
                    if ordinal != value.cells.len() {
                        return Err("nonconsecutive column ordinal".into());
                    }
                    value.cells.push(decode_cell(&row[3], &row[4])?);
                }
            }
            _ => return Err("unknown argument record".into()),
        }
    }
    if values.len() != count
        || values.iter().any(|value| {
            value.kind == "scalar" && value.cells.len() != 1
                || value.kind == "series" && value.index.len() != value.cells.len()
                || value.kind != "series" && !value.index.is_empty()
        })
    {
        return Err("incomplete argument denominator".into());
    }
    Ok((sequence, operation, values))
}

pub(super) fn response(sequence: u64, result: Result<Value, String>) -> Result<Vec<u8>, String> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(Vec::new());
    match result {
        Err(reason) => writer
            .write_record([
                "unsupported",
                &sequence.to_string(),
                &reason,
                "",
                "",
                "insilico-columns-v1",
            ])
            .map_err(|error| error.to_string())?,
        Ok(value) => {
            writer
                .write_record([
                    "result",
                    &sequence.to_string(),
                    &value.kind,
                    &value.dtype,
                    &value.cells.len().to_string(),
                    "insilico-columns-v1",
                ])
                .map_err(|error| error.to_string())?;
            for (i, label) in value.index.iter().enumerate() {
                writer
                    .write_record(["index", "0", &i.to_string(), "label", label, ""])
                    .map_err(|error| error.to_string())?;
            }
            for (i, cell) in value.cells.iter().enumerate() {
                let (tag, cell) = encode_cell(cell);
                writer
                    .write_record(["cell", "0", &i.to_string(), tag, &cell, ""])
                    .map_err(|error| error.to_string())?;
            }
        }
    }
    writer.into_inner().map_err(|error| error.to_string())
}

pub(super) fn read_frame(reader: &mut impl Read) -> Result<Option<Vec<u8>>, String> {
    let mut size = [0u8; 4];
    match reader.read(&mut size[..1]) {
        Ok(0) => return Ok(None),
        Ok(1) => (),
        _ => return Err("failed frame prefix".into()),
    }
    reader
        .read_exact(&mut size[1..])
        .map_err(|_| "incomplete frame prefix")?;
    let size = u32::from_be_bytes(size) as usize;
    if size == 0 || size > MAX_FRAME {
        return Err("invalid frame length".into());
    }
    let mut bytes = vec![0; size];
    reader
        .read_exact(&mut bytes)
        .map_err(|_| "incomplete frame body")?;
    Ok(Some(bytes))
}

fn session(reader: &mut impl Read, writer: &mut impl Write) -> Result<(), String> {
    let mut expected = 0;
    while let Some(bytes) = read_frame(reader)? {
        let (sequence, operation, values) = request(&bytes)?;
        if sequence != expected {
            return Err("request sequence differs".into());
        }
        expected += 1;
        let quit = operation == "quit";
        let result = if quit {
            if !values.is_empty() {
                return Err("quit operands forbidden".into());
            }
            Ok(Value::scalar(Cell::Text("closed".into())))
        } else {
            calculate(&operation, &values)
        };
        let bytes = response(sequence, result)?;
        if bytes.len() > MAX_FRAME {
            return Err("result exceeds frame limit".into());
        }
        writer
            .write_all(&(bytes.len() as u32).to_be_bytes())
            .map_err(|error| error.to_string())?;
        writer
            .write_all(&bytes)
            .map_err(|error| error.to_string())?;
        writer.flush().map_err(|error| error.to_string())?;
        if quit {
            return Ok(());
        }
    }
    Err("session ended without explicit quit".into())
}

pub(crate) fn dispatch(args: &[String]) -> crate::cli::Result<()> {
    if args == ["--help"] {
        println!(
            "priva insilico-criteria session|procedure-session\nTransitional column API: insilico-columns-v1. Rust-owned named PP3/BP4 and BP7 control: insilico-procedures-v2, source b8a093a9. Python remains the default; mature operations and callbacks remain explicit."
        );
        return Ok(());
    }
    if args == ["procedure-session"] {
        return super::procedure_transport::session(
            &mut std::io::stdin().lock(),
            &mut std::io::stdout().lock(),
        )
        .map_err(crate::cli::invalid);
    }
    if args != ["session"] {
        return Err(crate::cli::invalid(
            "expected session, procedure-session or --help",
        ));
    }
    session(&mut std::io::stdin().lock(), &mut std::io::stdout().lock())
        .map_err(crate::cli::invalid)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn malformed_ordinals_and_incomplete_frames_fail_closed() {
        let raw = b"request\t0\tsum\t1\t\tinsilico-columns-v1\narg\t0\tarray\tint64\t\t\ncell\t0\t1\tint\t4\t\n";
        assert!(request(raw).is_err());
        assert!(read_frame(&mut &[0, 0, 0, 4, 1][..]).is_err());
    }

    #[test]
    fn unicode_and_float_bits_roundtrip_without_json_coercion() {
        assert_eq!(
            decode_cell("str", "λ\t\n").unwrap(),
            Cell::Text("λ\t\n".into())
        );
        let original = Cell::Float(f64::from_bits(0x8000000000000000));
        let (tag, raw) = encode_cell(&original);
        let Cell::Float(decoded) = decode_cell(tag, &raw).unwrap() else {
            panic!()
        };
        assert_eq!(decoded.to_bits(), 0x8000000000000000);
    }
}
