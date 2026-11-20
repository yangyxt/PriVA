//! Pinned pandas C-parser block inference and concatenation.
//! Tokenization, implicit indexes and chunk boundaries remain caller-owned.
use crate::cli::Result;
use crate::pandas_table::{self, Cell, Column, ObjectValue};

#[derive(Clone)]
pub(crate) struct Values {
    pub(crate) column: Column,
    pub(crate) unsigned: Option<Vec<u64>>,
}

pub(crate) fn infer(tokens: Vec<String>) -> Result<Values> {
    match pandas_table::infer_column(tokens.clone(), true) {
        Ok(column) => Ok(Values {
            column,
            unsigned: None,
        }),
        Err(error)
            if error
                .to_string()
                .contains("outside the signed 64-bit integer range") =>
        {
            let mut numbers = Vec::new();
            let mut negative = false;
            let mut large = false;
            let mut valid = true;
            for token in tokens.iter().filter(|token| !pandas_table::is_na(token)) {
                let token = token.trim_matches([' ', '\t', '\n', '\r', '\x0b', '\x0c']);
                if token.starts_with('-') {
                    negative = true;
                    continue;
                }
                let digits = token.strip_prefix('+').unwrap_or(token);
                let value = (!digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()))
                    .then(|| digits.parse::<u64>().ok())
                    .flatten();
                if let Some(value) = value {
                    large |= value > i64::MAX as u64;
                    numbers.push(value);
                } else {
                    valid = false;
                    break;
                }
            }
            if valid && large && !negative && !tokens.iter().any(|token| pandas_table::is_na(token))
            {
                let column = Column::Text(numbers.iter().map(|n| Some(n.to_string())).collect());
                return Ok(Values {
                    column,
                    unsigned: Some(numbers),
                });
            }
            // pandas cancels uint64 conversion with NA/negative values and
            // then preserves their literal spelling; true overflow uses normal NA.
            let keep_na = valid && large;
            Ok(Values {
                column: Column::Text(
                    tokens
                        .into_iter()
                        .map(|s| {
                            if !keep_na && pandas_table::is_na(&s) {
                                None
                            } else {
                                Some(s)
                            }
                        })
                        .collect(),
                ),
                unsigned: None,
            })
        }
        Err(error) => Err(error),
    }
}

fn object(column: &Column, row: usize) -> ObjectValue {
    match column.cell(row) {
        Cell::Missing => ObjectValue::Nan,
        Cell::Text(s) => ObjectValue::Text(s.to_owned()),
        Cell::Int(n) => ObjectValue::Int(n),
        Cell::Float(n) => ObjectValue::Float(n),
        Cell::Bool(b) => ObjectValue::Bool(b),
    }
}

pub(crate) fn concatenate(blocks: Vec<Values>) -> Values {
    if blocks.len() == 1 {
        return blocks.into_iter().next().unwrap();
    }
    let all_integer = blocks
        .iter()
        .all(|b| matches!(b.column, Column::Int(_)) && b.unsigned.is_none());
    let all_numeric = blocks
        .iter()
        .all(|b| matches!(b.column, Column::Int(_) | Column::Float(_)) || b.unsigned.is_some());
    let all_unsigned = blocks.iter().all(|b| b.unsigned.is_some());
    let all_bool = blocks.iter().all(|b| matches!(b.column, Column::Bool(_)));
    let all_text = blocks
        .iter()
        .all(|b| matches!(b.column, Column::Text(_)) && b.unsigned.is_none());
    if all_unsigned {
        let values: Vec<u64> = blocks
            .into_iter()
            .flat_map(|b| b.unsigned.unwrap())
            .collect();
        return Values {
            column: Column::Text(values.iter().map(|n| Some(n.to_string())).collect()),
            unsigned: Some(values),
        };
    }
    let column = if all_integer {
        Column::Int(
            blocks
                .iter()
                .flat_map(|b| {
                    (0..b.column.len()).map(|r| match b.column.cell(r) {
                        Cell::Int(n) => n,
                        _ => unreachable!(),
                    })
                })
                .collect(),
        )
    } else if all_numeric {
        Column::Float(
            blocks
                .iter()
                .flat_map(|b| {
                    (0..b.column.len()).map(|r| {
                        if let Some(values) = &b.unsigned {
                            values[r] as f64
                        } else {
                            match b.column.cell(r) {
                                Cell::Int(n) => n as f64,
                                Cell::Float(n) => n,
                                Cell::Missing => f64::NAN,
                                _ => unreachable!(),
                            }
                        }
                    })
                })
                .collect(),
        )
    } else if all_bool {
        Column::Bool(
            blocks
                .iter()
                .flat_map(|b| (0..b.column.len()).map(|r| b.column.cell(r) == Cell::Bool(true)))
                .collect(),
        )
    } else if all_text {
        Column::Text(
            blocks
                .iter()
                .flat_map(|b| {
                    (0..b.column.len()).map(|r| match b.column.cell(r) {
                        Cell::Text(s) => Some(s.to_owned()),
                        _ => None,
                    })
                })
                .collect(),
        )
    } else {
        Column::Object(
            blocks
                .iter()
                .flat_map(|b| (0..b.column.len()).map(|r| object(&b.column, r)))
                .collect(),
        )
    };
    Values {
        column,
        unsigned: None,
    }
}
