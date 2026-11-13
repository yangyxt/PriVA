//! The frame file: a pandas DataFrame written exactly, value by value, so
//! that a Python function's input and output can be handed to its Rust
//! version and compared without any type inference in between. The Python
//! side is `rust/tests/acmg_step_oracle_common.py` (`write_frame`, `read_frame`).
//!
//! Layout (tab-separated text, `\n` line ends):
//!
//! ```text
//! priva-frame <TAB> 1 <TAB> <number of rows>
//! <column names, each as the hex of its UTF-8 bytes>
//! <pandas dtypes: int64, float64, bool, object, Int64 or boolean>
//! <one line per row, one cell per column>
//! ```
//!
//! A cell is `n` (the float NaN), `N` (`None`), `A` (`pd.NA`), `b0` or `b1`,
//! `i` and a decimal integer, `f` and the 16 hex digits of the IEEE bits of a
//! double, or `s` and the hex of the UTF-8 bytes of a text. The row index is
//! not stored: frames carry the default index `0..n`.
use crate::cli::{Result, invalid};
use crate::pandas_table::{Column, ObjectValue, Table};
use std::fs;
use std::path::Path;

const MAGIC: &str = "priva-frame";

fn hex(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        out.push_str(&format!("{byte:02x}"));
    }
    out
}

fn unhex(text: &str) -> Result<String> {
    if text.len() % 2 != 0 || !text.is_ascii() {
        return Err(invalid(format!("frame: {text:?} is not hex")));
    }
    let bytes = (0..text.len())
        .step_by(2)
        .map(|index| u8::from_str_radix(&text[index..index + 2], 16))
        .collect::<std::result::Result<Vec<u8>, _>>()
        .map_err(|_| invalid(format!("frame: {text:?} is not hex")))?;
    String::from_utf8(bytes).map_err(|_| invalid("frame: text is not UTF-8"))
}

fn object_cell(value: &ObjectValue) -> String {
    match value {
        ObjectValue::Nan => "n".to_owned(),
        ObjectValue::None => "N".to_owned(),
        ObjectValue::Text(text) => format!("s{}", hex(text.as_bytes())),
        ObjectValue::Int(number) => format!("i{number}"),
        ObjectValue::Float(number) if number.is_nan() => "n".to_owned(),
        ObjectValue::Float(number) => format!("f{:016x}", number.to_bits()),
        ObjectValue::Bool(flag) => format!("b{}", u8::from(*flag)),
    }
}

fn cell_text(column: &Column, row: usize) -> String {
    match column {
        Column::Int(values) => format!("i{}", values[row]),
        Column::Float(values) => object_cell(&ObjectValue::Float(values[row])),
        Column::Bool(values) => format!("b{}", u8::from(values[row])),
        Column::BoolOrMissing(values) => match values[row] {
            Some(flag) => format!("b{}", u8::from(flag)),
            None => "n".to_owned(),
        },
        Column::Text(values) => match &values[row] {
            Some(text) => format!("s{}", hex(text.as_bytes())),
            None => "n".to_owned(),
        },
        Column::Object(values) => object_cell(&values[row]),
        Column::NullableInt(values) => match values[row] {
            Some(number) => format!("i{number}"),
            None => "A".to_owned(),
        },
        Column::NullableBool(values) => match values[row] {
            Some(flag) => format!("b{}", u8::from(flag)),
            None => "A".to_owned(),
        },
    }
}

/// The table as frame-file text.
pub(crate) fn frame_text(table: &Table) -> String {
    let mut out = format!("{MAGIC}\t1\t{}\n", table.rows);
    let names: Vec<String> = table
        .names
        .iter()
        .map(|name| hex(name.as_bytes()))
        .collect();
    out.push_str(&names.join("\t"));
    out.push('\n');
    let dtypes: Vec<&str> = table.columns.iter().map(Column::dtype_name).collect();
    out.push_str(&dtypes.join("\t"));
    out.push('\n');
    for row in 0..table.rows {
        let cells: Vec<String> = table
            .columns
            .iter()
            .map(|column| cell_text(column, row))
            .collect();
        out.push_str(&cells.join("\t"));
        out.push('\n');
    }
    out
}

fn decode_cell(cell: &str) -> Result<ObjectValue> {
    if !cell.is_ascii() {
        return Err(invalid(format!("frame: bad cell {cell:?}")));
    }
    let (kind, payload) = cell.split_at(cell.len().min(1));
    Ok(match kind {
        "n" if payload.is_empty() => ObjectValue::Nan,
        "N" if payload.is_empty() => ObjectValue::None,
        "b" if payload == "0" => ObjectValue::Bool(false),
        "b" if payload == "1" => ObjectValue::Bool(true),
        "i" => ObjectValue::Int(
            payload
                .parse()
                .map_err(|_| invalid(format!("frame: bad integer cell {cell:?}")))?,
        ),
        "f" if payload.len() == 16 => ObjectValue::Float(f64::from_bits(
            u64::from_str_radix(payload, 16)
                .map_err(|_| invalid(format!("frame: bad float cell {cell:?}")))?,
        )),
        "s" => ObjectValue::Text(unhex(payload)?),
        _ => return Err(invalid(format!("frame: bad cell {cell:?}"))),
    })
}

fn wrong(dtype: &str, cell: &str) -> Box<dyn std::error::Error> {
    invalid(format!("frame: cell {cell:?} in a {dtype} column"))
}

fn build_column(dtype: &str, cells: &[&str]) -> Result<Column> {
    let decoded: Vec<ObjectValue> = cells
        .iter()
        .filter(|cell| **cell != "A")
        .map(|cell| decode_cell(cell))
        .collect::<Result<_>>()?;
    let nullable_values = |cells: &[&str]| -> Result<Vec<Option<ObjectValue>>> {
        cells
            .iter()
            .map(|cell| {
                if *cell == "A" {
                    Ok(None)
                } else {
                    decode_cell(cell).map(Some)
                }
            })
            .collect()
    };
    if dtype != "Int64" && dtype != "boolean" && decoded.len() != cells.len() {
        return Err(invalid(format!("frame: pd.NA in a {dtype} column")));
    }
    Ok(match dtype {
        "int64" => Column::Int(
            decoded
                .iter()
                .zip(cells)
                .map(|(value, cell)| match value {
                    ObjectValue::Int(number) => Ok(*number),
                    _ => Err(wrong(dtype, cell)),
                })
                .collect::<Result<_>>()?,
        ),
        "float64" => Column::Float(
            decoded
                .iter()
                .zip(cells)
                .map(|(value, cell)| match value {
                    ObjectValue::Float(number) => Ok(*number),
                    ObjectValue::Nan => Ok(f64::NAN),
                    _ => Err(wrong(dtype, cell)),
                })
                .collect::<Result<_>>()?,
        ),
        "bool" => Column::Bool(
            decoded
                .iter()
                .zip(cells)
                .map(|(value, cell)| match value {
                    ObjectValue::Bool(flag) => Ok(*flag),
                    _ => Err(wrong(dtype, cell)),
                })
                .collect::<Result<_>>()?,
        ),
        "object" => {
            let text_or_nan = decoded
                .iter()
                .all(|value| matches!(value, ObjectValue::Text(_) | ObjectValue::Nan));
            let bool_or_nan = decoded
                .iter()
                .all(|value| matches!(value, ObjectValue::Bool(_) | ObjectValue::Nan));
            let any_bool = decoded
                .iter()
                .any(|value| matches!(value, ObjectValue::Bool(_)));
            if text_or_nan {
                Column::Text(
                    decoded
                        .into_iter()
                        .map(|value| match value {
                            ObjectValue::Text(text) => Some(text),
                            _ => None,
                        })
                        .collect(),
                )
            } else if bool_or_nan && any_bool {
                Column::BoolOrMissing(
                    decoded
                        .into_iter()
                        .map(|value| match value {
                            ObjectValue::Bool(flag) => Some(flag),
                            _ => None,
                        })
                        .collect(),
                )
            } else {
                Column::Object(decoded)
            }
        }
        "Int64" => Column::NullableInt(
            nullable_values(cells)?
                .into_iter()
                .zip(cells)
                .map(|(value, cell)| match value {
                    None => Ok(None),
                    Some(ObjectValue::Int(number)) => Ok(Some(number)),
                    _ => Err(wrong(dtype, cell)),
                })
                .collect::<Result<_>>()?,
        ),
        "boolean" => Column::NullableBool(
            nullable_values(cells)?
                .into_iter()
                .zip(cells)
                .map(|(value, cell)| match value {
                    None => Ok(None),
                    Some(ObjectValue::Bool(flag)) => Ok(Some(flag)),
                    _ => Err(wrong(dtype, cell)),
                })
                .collect::<Result<_>>()?,
        ),
        other => {
            return Err(invalid(format!(
                "frame: dtype {other:?} is not one the frame file carries"
            )));
        }
    })
}

/// Reads frame-file text.
pub(crate) fn parse_frame(text: &str) -> Result<Table> {
    let mut lines = text.split('\n');
    let head: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    if head.len() != 3 || head[0] != MAGIC || head[1] != "1" {
        return Err(invalid("frame: not a version 1 frame file"));
    }
    let rows: usize = head[2]
        .parse()
        .map_err(|_| invalid("frame: bad row count"))?;
    let split = |line: &str| -> Vec<String> {
        if line.is_empty() {
            Vec::new()
        } else {
            line.split('\t').map(str::to_owned).collect()
        }
    };
    let names = split(
        lines
            .next()
            .ok_or_else(|| invalid("frame: no names line"))?,
    )
    .iter()
    .map(|name| unhex(name))
    .collect::<Result<Vec<String>>>()?;
    let dtypes = split(
        lines
            .next()
            .ok_or_else(|| invalid("frame: no dtypes line"))?,
    );
    if dtypes.len() != names.len() {
        return Err(invalid("frame: names and dtypes differ in number"));
    }
    let mut cells: Vec<Vec<&str>> = vec![Vec::with_capacity(rows); names.len()];
    for row in 0..rows {
        let line = lines
            .next()
            .ok_or_else(|| invalid(format!("frame: row {row} is missing")))?;
        let fields: Vec<&str> = if names.is_empty() {
            Vec::new()
        } else {
            line.split('\t').collect()
        };
        if fields.len() != names.len() {
            return Err(invalid(format!(
                "frame: row {row} has {} cells",
                fields.len()
            )));
        }
        for (column, field) in cells.iter_mut().zip(fields) {
            column.push(field);
        }
    }
    if lines.any(|line| !line.is_empty()) {
        return Err(invalid("frame: more lines than rows"));
    }
    let columns = dtypes
        .iter()
        .zip(&cells)
        .map(|(dtype, column)| build_column(dtype, column))
        .collect::<Result<Vec<Column>>>()?;
    Ok(Table {
        names,
        columns,
        rows,
    })
}

pub(crate) fn read_frame(path: &Path) -> Result<Table> {
    let text = fs::read_to_string(path)
        .map_err(|error| invalid(format!("{}: {error}", path.display())))?;
    parse_frame(&text)
}

pub(crate) fn write_frame(path: &Path, table: &Table) -> Result<()> {
    fs::write(path, frame_text(table))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_column_type_round_trips() {
        let table = Table {
            names: vec![
                "i".into(),
                "f".into(),
                "b".into(),
                "o_bool".into(),
                "o_text".into(),
                "o_mix".into(),
                "#Gene Symbol\t2".into(),
                "B".into(),
            ],
            columns: vec![
                Column::Int(vec![1, -9_223_372_036_854_775_808]),
                Column::Float(vec![f64::NAN, -0.0]),
                Column::Bool(vec![true, false]),
                Column::BoolOrMissing(vec![None, Some(false)]),
                Column::Text(vec![Some("é\t\"x\"".into()), None]),
                Column::Object(vec![ObjectValue::None, ObjectValue::Float(1e-5)]),
                Column::NullableInt(vec![None, Some(7)]),
                Column::NullableBool(vec![Some(true), None]),
            ],
            rows: 2,
        };
        let text = frame_text(&table);
        let back = parse_frame(&text).unwrap();
        assert_eq!(frame_text(&back), text);
        assert_eq!(back.names, table.names);
        assert_eq!(back.columns[2..], table.columns[2..]);
        assert!(text.contains("\tf8000000000000000\t"));
    }

    #[test]
    fn object_columns_take_the_narrowest_form() {
        let text =
            "priva-frame\t1\t2\n61\t62\t63\nobject\tobject\tobject\nb1\ts78\ti1\nn\tn\ts78\n";
        let table = parse_frame(text).unwrap();
        assert_eq!(
            table.columns[0],
            Column::BoolOrMissing(vec![Some(true), None])
        );
        assert_eq!(table.columns[1], Column::Text(vec![Some("x".into()), None]));
        assert_eq!(
            table.columns[2],
            Column::Object(vec![ObjectValue::Int(1), ObjectValue::Text("x".into())])
        );
    }

    #[test]
    fn malformed_frames_are_refused() {
        assert!(parse_frame("priva-frame\t2\t0\n\n\n").is_err());
        assert!(parse_frame("priva-frame\t1\t1\n61\nint64\nf0000000000000000\n").is_err());
        assert!(parse_frame("priva-frame\t1\t1\n61\nfloat64\nA\n").is_err());
        assert!(parse_frame("priva-frame\t1\t1\n61\nstring\nsab\n").is_err());
        assert!(parse_frame("priva-frame\t1\t0\n61\nint64\nextra\n").is_err());
        let empty = parse_frame("priva-frame\t1\t2\n\n\n\n\n").unwrap();
        assert_eq!((empty.names.len(), empty.rows), (0, 2));
    }
}
