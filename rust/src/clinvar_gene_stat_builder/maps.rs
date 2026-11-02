//! Source-ordered summary/MANE and VCF AlleleID dictionaries.
use crate::cli::{Result, invalid};
use crate::pandas_table::{self, Cell, Column};
use crate::pickle_reader::{self, Value};
use crate::vcf_variant::{self, Scalar};
use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(super) enum Key {
    Missing,
    Int(i64),
    Float(u64),
    Text(String),
    Tuple(Vec<Key>),
}
impl Key {
    fn float(value: f64) -> Self {
        if value.is_nan() {
            Self::Missing
        } else if value.fract() == 0.0
            && (-((1_u64 << 63) as f64)..(1_u64 << 63) as f64).contains(&value)
        {
            Self::Int(value as i64)
        } else {
            Self::Float(value.to_bits())
        }
    }
    fn cell(value: Cell<'_>) -> Self {
        match value {
            Cell::Missing => Self::Missing,
            Cell::Int(value) => Self::Int(value),
            Cell::Float(value) => Self::float(value),
            Cell::Bool(value) => Self::Int(i64::from(value)),
            Cell::Text(value) => Self::Text(value.to_owned()),
        }
    }
    fn scalar(value: Scalar) -> Self {
        match value {
            Scalar::Missing => Self::Missing,
            Scalar::Int(value) => Self::Int(value),
            Scalar::Float(value) => Self::float(value),
            Scalar::Bool(value) => Self::Int(i64::from(value)),
            Scalar::Text(value) => Self::Text(value),
        }
    }
    fn pickle(value: &Value) -> Result<Self> {
        Ok(match value {
            Value::None => Self::Missing,
            Value::Bool(value) => Self::Int(i64::from(*value)),
            Value::Int(value) => Self::Int(*value),
            Value::Float(value) => Self::float(*value),
            Value::Str(value) => Self::Text(value.to_string()),
            _ => return Err(invalid("unsupported allele-to-transcript map key")),
        })
    }
    pub(super) fn promoted(&self) -> Self {
        if let Self::Int(value) = self {
            Self::float(*value as f64)
        } else {
            self.clone()
        }
    }
}

pub(super) type TargetMap = HashMap<Key, Option<String>>;

/// `Name.str.extract(r'(NM_\d+\.\d+)')`: the first match, with Python's
/// Unicode decimal-digit class. The matched text itself is not normalized.
fn refseq(text: &str) -> Option<String> {
    for (at, _) in text.match_indices("NM_") {
        let rest = &text[at + 3..];
        let digits = rest
            .char_indices()
            .take_while(|(_, c)| crate::python_text::decimal_digit_value(*c).is_some())
            .map(|(at, c)| at + c.len_utf8())
            .last()
            .unwrap_or(0);
        if digits == 0 || !rest[digits..].starts_with('.') {
            continue;
        }
        let version = &rest[digits + 1..];
        let count = version
            .char_indices()
            .take_while(|(_, c)| crate::python_text::decimal_digit_value(*c).is_some())
            .map(|(at, c)| at + c.len_utf8())
            .last()
            .unwrap_or(0);
        if count > 0 {
            return Some(text[at..at + 3 + digits + 1 + count].to_owned());
        }
    }
    None
}

pub(super) fn from_tables(summary: &Path, mane: &Path) -> Result<TargetMap> {
    let clinvar = pandas_table::read_table(summary)?;
    let assembly = clinvar.column("Assembly")?;
    let names = clinvar.column("Name")?;
    let ids = clinvar.column("#AlleleID")?;
    if !matches!(names, Column::Text(_) | Column::Object(_)) {
        return Err(invalid(
            "AttributeError: the ClinVar Name column has no .str accessor",
        ));
    }
    let mane = pandas_table::read_table(mane)?;
    let source = mane.column("RefSeq_nuc")?;
    let target = mane.column("Ensembl_nuc")?;
    // Iterating a float Series creates distinct NaN Python objects. A dict
    // therefore retains repeated NaN keys, and Series.map refuses their
    // non-unique pandas Index. Object-column missing values share np.nan.
    if matches!(source, Column::Float(values) if values.iter().filter(|v| v.is_nan()).count() > 1) {
        return Err(invalid("InvalidIndexError: repeated NaN MANE keys"));
    }
    let mut refseq_to_ensembl = HashMap::new();
    for row in 0..mane.rows {
        refseq_to_ensembl.insert(
            Key::cell(source.cell(row)),
            match target.cell(row) {
                Cell::Text(text) => Some(text.split('.').next().unwrap_or("").to_owned()),
                _ => None,
            },
        );
    }
    let mut result = HashMap::new();
    let mut has_text = false;
    let mut missing_numeric_ids = 0;
    for row in 0..clinvar.rows {
        if assembly.cell(row) != Cell::Text("GRCh38") {
            continue;
        }
        if matches!(ids, Column::Float(_)) && ids.cell(row) == Cell::Missing {
            missing_numeric_ids += 1;
        }
        let name = match names.cell(row) {
            Cell::Text(text) => refseq(text),
            _ => None,
        };
        let key = name.map_or(Key::Missing, Key::Text);
        let transcript = refseq_to_ensembl.get(&key).cloned().flatten();
        has_text |= transcript.is_some();
        result.insert(Key::cell(ids.cell(row)), transcript);
    }
    // The original .str.split(..., expand=True).iloc[:,0] also fails when no
    // mapped text survives; it does not silently create an empty target map.
    if !has_text {
        return Err(invalid("the mapped transcript column has no text to split"));
    }
    if missing_numeric_ids > 1 {
        return Err(invalid("InvalidIndexError: repeated NaN AlleleID keys"));
    }
    Ok(result)
}

/// Explicit reuse of a source-bound saved procedure result. The ordinary
/// restricted pickle reader accepts only plain data here, never reducers.
pub(super) fn from_pickle(path: &Path) -> Result<TargetMap> {
    let value = pickle_reader::load(&pickle_reader::read_file_by_magic(path)?)?;
    let Value::Dict(dict) = value else {
        return Err(invalid("the saved transcript map is not a dictionary"));
    };
    let dict = dict.borrow();
    let mut map = HashMap::with_capacity(dict.entries.len());
    for (key, value) in &dict.entries {
        let value = match value {
            Value::Str(text) => Some(text.to_string()),
            Value::Float(value) if value.is_nan() => None,
            _ => {
                return Err(invalid(
                    "the saved transcript map value is not text or float NaN",
                ));
            }
        };
        map.insert(Key::pickle(key)?, value);
    }
    Ok(map)
}

pub(super) struct Alleles {
    pub(super) values: HashMap<String, Key>,
    pub(super) numeric: bool,
    pub(super) nullable: bool,
}
pub(super) fn alleles(path: &Path) -> Result<Alleles> {
    let header = vcf_variant::header(path)?;
    let mut types = header.info.clone();
    let mut values = HashMap::new();
    let mut body = false;
    for line in crate::vcf_text::open_reader(path)?.lines() {
        let line = line?;
        if !body {
            if line.starts_with("#CHROM") {
                body = true;
            }
            continue;
        }
        let record = vcf_variant::decode_record(line.as_bytes(), &header, &mut types)?;
        let alt = record
            .fields
            .get(4)
            .filter(|value| **value != ".")
            .ok_or_else(|| invalid("TypeError: record.alts is None"))?
            .split(',')
            .next()
            .unwrap();
        let reference = record.fields.get(3).copied().unwrap_or(".");
        let reference = if reference.is_empty() { "." } else { reference };
        let alt = if alt.is_empty() { "." } else { alt };
        let variant = format!(
            "{}_{}_{}_{}",
            record.fields[0], record.position, reference, alt
        );
        if record.info.contains_key("ALLELEID") {
            let value = record.get("ALLELEID", &types)?.unwrap();
            let key = if value.scalar {
                Key::scalar(value.values.into_iter().next().unwrap())
            } else {
                Key::Tuple(value.values.into_iter().map(Key::scalar).collect())
            };
            values.insert(variant, key);
        }
    }
    let numeric = values
        .values()
        .all(|value| matches!(value, Key::Missing | Key::Int(_) | Key::Float(_)));
    let nullable = values
        .values()
        .any(|value| matches!(value, Key::Missing | Key::Float(_)));
    Ok(Alleles {
        values,
        numeric,
        nullable,
    })
}
