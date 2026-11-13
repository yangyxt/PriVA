//! The coldspot producer's `read_csv(sep="\t")` and selected-row conversion.
//! Unlike the shared low_memory=False reader, pandas infers each C-parser
//! block before concatenation. Index labels are unused by this source loop.
use crate::builder_error;
use crate::cli::{Result, invalid};
use crate::json_text::Json;
use crate::pandas_low_memory::{Values, concatenate, infer};
use crate::pandas_table::{self, Cell, Column};
use crate::python_text::{decimal_digit_value, python_int_is_valid, python_strip};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

fn exception(name: &'static str, message: &str) -> Box<dyn std::error::Error> {
    builder_error::category(name, invalid(message))
}

#[derive(Clone)]
pub(crate) enum Gene {
    Text(String),
    Integer(String),
    Float(f64),
    Bool(bool),
}

impl Gene {
    pub(crate) fn json(&self) -> Json {
        let (kind, name, value) = match self {
            Self::Text(s) => ("str", "value", Json::Str(s.clone())),
            Self::Integer(s) => ("int", "value", Json::Str(s.clone())),
            Self::Float(f) => ("float", "bits", Json::Str(format!("{:016x}", f.to_bits()))),
            Self::Bool(b) => ("bool", "value", Json::Bool(*b)),
        };
        Json::Object(vec![
            ("kind".into(), Json::Str(kind.into())),
            (name.into(), value),
        ])
    }
}

pub(crate) struct Row {
    pub(crate) gene: Gene,
    // Python int has no machine-word limit. Keep its exact canonical decimal
    // spelling here; the downstream interval procedure owns its representation.
    pub(crate) start: String,
    pub(crate) end: String,
}

struct ParsedColumn {
    values: Values,
    // A uint64 block concatenated with object values becomes Python integers,
    // not strings. The shared CSV representation keeps their decimal spelling;
    // this caller also retains the original integer type for iterrows().
    object_unsigned: HashMap<usize, u64>,
}

impl ParsedColumn {
    fn cell(&self, row: usize) -> Gene {
        if let Some(values) = &self.values.unsigned {
            return Gene::Integer(values[row].to_string());
        }
        if let Some(value) = self.object_unsigned.get(&row) {
            return Gene::Integer(value.to_string());
        }
        match self.values.column.cell(row) {
            Cell::Missing => Gene::Float(f64::NAN),
            Cell::Text(s) => Gene::Text(s.to_owned()),
            Cell::Int(n) => Gene::Integer(n.to_string()),
            Cell::Float(f) => Gene::Float(f),
            Cell::Bool(b) => Gene::Bool(b),
        }
    }
}

fn integer(value: Gene) -> Result<String> {
    match value {
        Gene::Integer(n) => Ok(n),
        Gene::Bool(b) => Ok(u8::from(b).to_string()),
        Gene::Float(f) if f.is_nan() => Err(exception(
            "ValueError",
            "cannot convert float NaN to integer",
        )),
        Gene::Float(f) if !f.is_finite() => Err(exception(
            "OverflowError",
            "cannot convert float infinity to integer",
        )),
        Gene::Float(f) => Ok(if f == 0.0 {
            "0".into()
        } else {
            format!("{:.0}", f.trunc())
        }),
        Gene::Text(text) => {
            if !python_int_is_valid(&text) {
                return Err(exception("ValueError", "invalid literal for int()"));
            }
            let cleaned: String = python_strip(&text)
                .chars()
                .filter(|c| *c != '_')
                .map(|c| decimal_digit_value(c).map_or(c, |d| char::from(b'0' + d as u8)))
                .collect();
            let negative = cleaned.starts_with('-');
            let digits = cleaned
                .trim_start_matches(['+', '-'])
                .trim_start_matches('0');
            Ok(if digits.is_empty() {
                "0".into()
            } else if negative {
                format!("-{digits}")
            } else {
                digits.into()
            })
        }
    }
}

fn text(path: &Path) -> Result<String> {
    let file = File::open(path)?;
    let name = path.to_string_lossy().to_lowercase();
    let mut bytes = Vec::new();
    if name.ends_with(".gz") && !name.ends_with(".tar.gz") {
        let mut input = BufReader::new(file);
        loop {
            // gzip.py distinguishes an empty stream from a one-byte invalid
            // header. Its EOF path also skips zero padding between members.
            let mut magic = Vec::with_capacity(2);
            input.by_ref().take(2).read_to_end(&mut magic)?;
            if magic.is_empty() {
                break;
            }
            if magic != [0x1f, 0x8b] {
                return Err(exception("BadGzipFile", "Not a gzipped file"));
            }
            // The BufRead decoder consumes exactly one member. Preserve any
            // buffered following bytes for the next original-style header.
            flate2::bufread::GzDecoder::new(std::io::Cursor::new(magic).chain(&mut input))
                .read_to_end(&mut bytes)
                .map_err(builder_error::gzip)?;
            loop {
                let available = input.fill_buf()?;
                let zeros = available.iter().take_while(|b| **b == 0).count();
                if zeros == 0 {
                    break;
                }
                input.consume(zeros);
            }
        }
    } else if name.ends_with(".bz2") && !name.ends_with(".tar.bz2") {
        bzip2::read::MultiBzDecoder::new(file)
            .read_to_end(&mut bytes)
            .map_err(|e| {
                builder_error::category(
                    if e.kind() == std::io::ErrorKind::UnexpectedEof {
                        "EOFError"
                    } else {
                        "OSError"
                    },
                    e.into(),
                )
            })?;
    } else if name.ends_with(".xz") && !name.ends_with(".tar.xz") {
        xz2::read::XzDecoder::new_multi_decoder(file)
            .read_to_end(&mut bytes)
            .map_err(|e| {
                builder_error::category(
                    if e.kind() == std::io::ErrorKind::UnexpectedEof {
                        "EOFError"
                    } else {
                        "LZMAError"
                    },
                    e.into(),
                )
            })?;
    } else if name.ends_with(".zst") {
        let decoded = zstd::stream::read::Decoder::new(file)
            .map_err(|e| builder_error::category("ZstdError", e.into()))?
            .read_to_end(&mut bytes);
        // pandas uses python-zstandard's stream_reader, which returns the
        // decoded prefix at input EOF even if a frame remains incomplete.
        // Other frame/corruption failures retain their original category.
        // Pinned source: indygreg/python-zstandard 0.23.0,
        // c-ext/decompressionreader.c (read/readinto EOF paths).
        if let Err(error) = decoded {
            if error.kind() != std::io::ErrorKind::UnexpectedEof {
                return Err(builder_error::category("ZstdError", error.into()));
            }
        }
    } else if name.ends_with(".zip") {
        let mut archive = zip::ZipArchive::new(file)
            .map_err(|e| builder_error::category("BadZipFile", e.into()))?;
        if archive.len() != 1 {
            return Err(exception(
                "ValueError",
                "expected exactly one file in ZIP archive",
            ));
        }
        archive
            .by_index(0)
            .map_err(|e| builder_error::category("BadZipFile", e.into()))?
            .read_to_end(&mut bytes)
            .map_err(|e| builder_error::category("BadZipFile", e.into()))?;
    } else if [".tar", ".tar.gz", ".tar.bz2", ".tar.xz"]
        .iter()
        .any(|s| name.ends_with(s))
    {
        // tarfile.open(mode='r') detects its inner compression from the bytes,
        // independently of the suffix that selected the TAR container.
        let mut encoded = Vec::new();
        let mut file = file;
        file.read_to_end(&mut encoded)?;
        let reader: Box<dyn Read> = if encoded.starts_with(&[0x1f, 0x8b]) {
            Box::new(flate2::read::MultiGzDecoder::new(std::io::Cursor::new(
                encoded,
            )))
        } else if encoded.starts_with(b"BZh") {
            Box::new(bzip2::read::MultiBzDecoder::new(std::io::Cursor::new(
                encoded,
            )))
        } else if encoded.starts_with(b"\xfd7zXZ\0") {
            Box::new(xz2::read::XzDecoder::new_multi_decoder(
                std::io::Cursor::new(encoded),
            ))
        } else {
            Box::new(std::io::Cursor::new(encoded))
        };
        let mut archive = tar::Archive::new(reader);
        let mut entries = archive
            .entries()
            .map_err(|e| builder_error::category("ReadError", e.into()))?;
        let Some(entry) = entries.next() else {
            return Err(exception("ValueError", "zero files in TAR archive"));
        };
        let mut entry = entry.map_err(|e| builder_error::category("ReadError", e.into()))?;
        if !entry.header().entry_type().is_file() {
            return Err(exception("AssertionError", "TAR member is not a file"));
        }
        entry
            .read_to_end(&mut bytes)
            .map_err(|e| builder_error::category("ReadError", e.into()))?;
        if entries
            .next()
            .transpose()
            .map_err(|e| builder_error::category("ReadError", e.into()))?
            .is_some()
        {
            return Err(exception("ValueError", "multiple files in TAR archive"));
        }
    } else {
        let mut file = file;
        file.read_to_end(&mut bytes)?;
    }
    String::from_utf8(bytes).map_err(|e| builder_error::category("UnicodeDecodeError", e.into()))
}

pub(crate) fn read(path: &Path) -> Result<Vec<Row>> {
    let text = text(path)?;
    let mut names = None;
    let mut raw = Vec::new();
    let mut width = 0;
    let mut block_size = 1;
    // Reuse the existing builder traversal, including BOM and bare-CR rules.
    // Coldspot parsing owns row-width policy; the shared line index is unused.
    pandas_table::tokenize_builder_separated(&text, b'\t', true, |mut row, _| {
        for value in &mut row {
            if let Some(at) = value.find('\0') {
                value.truncate(at);
            }
        }
        let Some(header) = &names else {
            names = Some(pandas_table::column_names(row));
            return Ok(());
        };
        if raw.is_empty() {
            width = row.len().max(header.len());
            let heuristic = (1 << 20) / width.max(1);
            while block_size * 2 < heuristic {
                block_size *= 2;
            }
        }
        if row.len() > width && raw.len() % block_size != 0 {
            return Err(exception("ParserError", "too many fields in a table row"));
        }
        row.resize(width, String::new());
        raw.push(row);
        Ok(())
    })
    .map_err(|e| builder_error::category("ParserError", e))?;
    let names =
        names.ok_or_else(|| exception("EmptyDataError", "No columns to parse from file"))?;
    if raw.is_empty() {
        width = names.len();
    }
    let leading = width - names.len();
    let mut columns = Vec::new();
    for col in leading..width {
        let blocks = if raw.is_empty() {
            vec![infer(Vec::new())?]
        } else {
            raw.chunks(block_size)
                .map(|rows| infer(rows.iter().map(|row| row[col].clone()).collect()))
                .collect::<Result<Vec<_>>>()?
        };
        let mut uints = HashMap::new();
        let mut offset = 0;
        for block in &blocks {
            if let Some(values) = &block.unsigned {
                for (r, value) in values.iter().enumerate() {
                    uints.insert(offset + r, *value);
                }
            }
            offset += block.column.len();
        }
        let values = concatenate(blocks);
        if !matches!(values.column, Column::Object(_)) {
            uints.clear();
        }
        columns.push(ParsedColumn {
            values,
            object_unsigned: uints,
        });
    }
    let required = |name: &str| {
        names
            .iter()
            .position(|n| n == name)
            .ok_or_else(|| exception("KeyError", name))
    };
    let kind = required("type")?;
    let mut result = Vec::new();
    for row in 0..raw.len() {
        if !matches!(columns[kind].cell(row),Gene::Text(s) if s=="coldspot") {
            continue;
        }
        let gene = columns[required("gene")?].cell(row);
        let start = integer(columns[required("aa_start_pos")?].cell(row))?;
        let end = integer(columns[required("aa_end_pos")?].cell(row))?;
        result.push(Row { gene, start, end });
    }
    Ok(result)
}
