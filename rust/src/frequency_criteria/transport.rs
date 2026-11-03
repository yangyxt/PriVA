//! One bounded, length-framed CSV session for complete source-position phases.
//! Protocol failures are fail-closed. Semantic unsupported responses are
//! explicit and let the host execute that actual mature operation once.
use super::{calculate, Cell, Value};
use std::io::{Read, Write};

const MAX_FRAME: usize = 256 * 1024 * 1024;

fn decode_cell(tag: &str, raw: &str) -> Result<Cell, String> {
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

fn encode_cell(cell: &Cell) -> (&'static str, String) {
    match cell {
        Cell::None => ("none", String::new()),
        Cell::Na => ("na", String::new()),
        Cell::Bool(value) => ("bool", u8::from(*value).to_string()),
        Cell::Int(value) => ("int", value.to_string()),
        Cell::Float(value) => ("float", format!("{:016x}", value.to_bits())),
        Cell::Text(value) => ("str", value.clone()),
    }
}

fn request(bytes: &[u8]) -> Result<(u64, String, Vec<Value>), String> {
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
    if row.len() != 6 || &row[0] != "request" || &row[5] != "frequency-columns-ef2e701-v1" {
        return Err("request schema mismatch".into());
    }
    let sequence = row[1].parse::<u64>().map_err(|_| "invalid sequence")?;
    let operation = row[2].to_owned();
    let count = row[3]
        .parse::<usize>()
        .map_err(|_| "invalid operand count")?;
    if count > 32 || !super::valid_site(&row[4], &row[2]) {
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
                    || !["scalar", "series", "array", "list", "pairs"].contains(&&row[2])
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
                if argument + 1 != values.len() {
                    return Err("interleaved argument records".into());
                }
                let value = values.get_mut(argument).ok_or("argument before metadata")?;
                let ordinal = row[2]
                    .parse::<usize>()
                    .map_err(|_| "invalid element ordinal")?;
                if !row[5].is_empty() {
                    return Err("nonempty reserved field".into());
                }
                if &row[0] == "index" {
                    if value.kind != "series" || !value.cells.is_empty() || ordinal != value.index.len() || &row[3] != "label"
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

fn response(sequence: u64, result: Result<Value, String>) -> Result<Vec<u8>, String> {
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
                "frequency-columns-ef2e701-v1",
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
                    "frequency-columns-ef2e701-v1",
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

fn read_frame(reader: &mut impl Read) -> Result<Option<Vec<u8>>, String> {
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

pub fn session(reader: &mut impl Read, writer: &mut impl Write) -> Result<(), String> {
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

pub fn dispatch(args: &[String]) -> crate::cli::Result<()> {
    if args == ["--help"] {
        println!("priva frequency-criteria session\nComplete frequency table/rule API with unchanged mature SciPy/codecs, protocol frequency-columns-ef2e701-v1. Python remains the default; unsupported mature operations remain explicit.");
        return Ok(());
    }
    if args != ["session"] {
        return Err(crate::cli::invalid("expected session or --help"));
    }
    session(&mut std::io::stdin().lock(), &mut std::io::stdout().lock())
        .map_err(crate::cli::invalid)
}

