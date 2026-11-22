//! Source-ordered summary/MANE and VCF AlleleID dictionaries.
use crate::cli::{Result, invalid};
use crate::clingen_map_builder::{ProducerColumn, infer};
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
    Unsigned(u64),
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
        } else if value.fract() == 0.0 && ((1_u64 << 63) as f64..u64::MAX as f64).contains(&value) {
            Self::Unsigned(value as u64)
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
    fn column(column: &ProducerColumn, row: usize) -> Self {
        match column {
            ProducerColumn::Pandas(column) => Self::cell(column.cell(row)),
            ProducerColumn::Unsigned(values) => match i64::try_from(values[row]) {
                Ok(value) => Self::Int(value),
                Err(_) => Self::Unsigned(values[row]),
            },
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
        } else if let Self::Unsigned(value) = self {
            Self::float(*value as f64)
        } else {
            self.clone()
        }
    }
    fn large_numeric(&self) -> bool {
        match self {
            Self::Int(value) => value.unsigned_abs() >= 1_u64 << 53,
            Self::Unsigned(_) => true,
            Self::Float(bits) => f64::from_bits(*bits).abs() >= (1_u64 << 53) as f64,
            _ => false,
        }
    }
}

/// pandas coerces a numeric dict Index to float64 for floating query values.
/// Keep the original first-key order for collisions created by that coercion;
/// the last index slot wins even when an earlier dict key was updated later.
/// Ordinary small integer maps need no secondary values or ordering storage.
#[derive(Default)]
pub(super) struct TargetMap {
    values: HashMap<Key, Option<String>>,
    large_keys: Vec<Key>,
    rounded: HashMap<Key, Option<String>>,
}
impl TargetMap {
    fn insert(&mut self, key: Key, value: Option<String>) {
        let large = key.large_numeric();
        if self.values.insert(key.clone(), value).is_none() && large {
            self.large_keys.push(key);
        }
    }
    fn finish(mut self) -> Self {
        if self.values.keys().all(|key| {
            matches!(
                key,
                Key::Missing | Key::Int(_) | Key::Unsigned(_) | Key::Float(_)
            )
        }) {
            for key in &self.large_keys {
                self.rounded
                    .insert(key.promoted(), self.values[key].clone());
            }
        }
        self.large_keys.clear();
        self
    }
    pub(super) fn get(&self, key: &Key, floating: bool) -> Option<&Option<String>> {
        if floating && key.large_numeric() && !self.rounded.is_empty() {
            self.rounded.get(key)
        } else {
            self.values.get(key)
        }
    }
}

/// The source never uses either table's row index or unselected columns.
/// Retain pandas' implicit leading index and full-row syntax checks while
/// sharing the already supported uint64/object inference with ClinGen.
/// The ordinary table reader's restricted policy stays unchanged.
fn selected<const N: usize>(path: &Path, wanted: [&str; N]) -> Result<[ProducerColumn; N]> {
    let text = pandas_table::read_table_text(path)?;
    let text = text.strip_prefix('\u{feff}').unwrap_or(&text);
    let mut names = None;
    let mut positions = [0; N];
    let mut columns: [Vec<String>; N] = std::array::from_fn(|_| Vec::new());
    let mut width = None;
    pandas_table::tokenize(text, None, |record| {
        if names.is_none() {
            let header = pandas_table::column_names(record);
            for (at, name) in wanted.iter().enumerate() {
                positions[at] = header
                    .iter()
                    .position(|value| value == name)
                    .ok_or_else(|| invalid(format!("KeyError: '{name}'")))?;
            }
            names = Some(header.len());
            return Ok(());
        }
        let named = names.unwrap();
        let expected = *width.get_or_insert(record.len().max(named));
        if record.len() > expected {
            return Err(invalid("ParserError: table row has too many fields"));
        }
        let leading = expected - named;
        for (at, values) in columns.iter_mut().enumerate() {
            values.push(
                record
                    .get(leading + positions[at])
                    .cloned()
                    .unwrap_or_default(),
            );
        }
        Ok(())
    })?;
    if names.is_none() {
        return Err(invalid("EmptyDataError: no columns to parse"));
    }
    let columns = columns.into_iter().map(infer).collect::<Result<Vec<_>>>()?;
    Ok(columns.try_into().unwrap_or_else(|_| unreachable!()))
}

fn text_cell(column: &ProducerColumn, row: usize) -> Option<&str> {
    match column {
        ProducerColumn::Pandas(column) => match column.cell(row) {
            Cell::Text(value) => Some(value),
            _ => None,
        },
        ProducerColumn::Unsigned(_) => None,
    }
}

fn rows(column: &ProducerColumn) -> usize {
    match column {
        ProducerColumn::Pandas(column) => column.len(),
        ProducerColumn::Unsigned(values) => values.len(),
    }
}

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
    let [assembly, names, ids] = selected(summary, ["Assembly", "Name", "#AlleleID"])?;
    if !matches!(
        names,
        ProducerColumn::Pandas(Column::Text(_) | Column::Object(_))
    ) {
        return Err(invalid(
            "AttributeError: the ClinVar Name column has no .str accessor",
        ));
    }
    let [source, target] = selected(mane, ["RefSeq_nuc", "Ensembl_nuc"])?;
    // Iterating a float Series creates distinct NaN Python objects. A dict
    // therefore retains repeated NaN keys, and Series.map refuses their
    // non-unique pandas Index. Object-column missing values share np.nan.
    if matches!(&source, ProducerColumn::Pandas(Column::Float(values)) if values.iter().filter(|v| v.is_nan()).count() > 1)
    {
        return Err(invalid("InvalidIndexError: repeated NaN MANE keys"));
    }
    let mut refseq_to_ensembl = HashMap::new();
    for row in 0..rows(&source) {
        refseq_to_ensembl.insert(
            Key::column(&source, row),
            text_cell(&target, row).map(|text| text.split('.').next().unwrap_or("").to_owned()),
        );
    }
    let mut result = TargetMap::default();
    let mut has_text = false;
    let mut missing_numeric_ids = 0;
    for row in 0..rows(&ids) {
        if text_cell(&assembly, row) != Some("GRCh38") {
            continue;
        }
        if matches!(&ids, ProducerColumn::Pandas(Column::Float(values)) if values[row].is_nan()) {
            missing_numeric_ids += 1;
        }
        let name = text_cell(&names, row).and_then(refseq);
        let key = name.map_or(Key::Missing, Key::Text);
        let transcript = refseq_to_ensembl.get(&key).cloned().flatten();
        has_text |= transcript.is_some();
        result.insert(Key::column(&ids, row), transcript);
    }
    // The original .str.split(..., expand=True).iloc[:,0] also fails when no
    // mapped text survives; it does not silently create an empty target map.
    if !has_text {
        return Err(invalid(if result.values.is_empty() {
            "IndexError: the empty split frame has no first column"
        } else {
            "AttributeError: the mapped transcript column has no .str accessor"
        }));
    }
    if missing_numeric_ids > 1 {
        return Err(invalid("InvalidIndexError: repeated NaN AlleleID keys"));
    }
    Ok(result.finish())
}

/// Explicit reuse of a source-bound saved procedure result. The ordinary
/// restricted pickle reader accepts only plain data here, never reducers.
pub(super) fn from_pickle(path: &Path) -> Result<TargetMap> {
    let value = pickle_reader::load(&pickle_reader::read_file_by_magic(path)?)?;
    let Value::Dict(dict) = value else {
        return Err(invalid("the saved transcript map is not a dictionary"));
    };
    let dict = dict.borrow();
    let mut map = TargetMap {
        values: HashMap::with_capacity(dict.entries.len()),
        ..TargetMap::default()
    };
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
    Ok(map.finish())
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
    let mut floating = std::collections::HashSet::new();
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
            if value.scalar
                && value
                    .values
                    .iter()
                    .any(|value| matches!(value, Scalar::Float(_)))
            {
                floating.insert(variant.clone());
            } else {
                floating.remove(&variant);
            }
            let key = if value.scalar {
                Key::scalar(value.values.into_iter().next().unwrap())
            } else {
                Key::Tuple(value.values.into_iter().map(Key::scalar).collect())
            };
            values.insert(variant, key);
        }
    }
    let numeric = values.values().all(|value| {
        matches!(
            value,
            Key::Missing | Key::Int(_) | Key::Unsigned(_) | Key::Float(_)
        )
    });
    let nullable = !floating.is_empty()
        || values
            .values()
            .any(|value| matches!(value, Key::Missing | Key::Float(_)));
    Ok(Alleles {
        values,
        numeric,
        nullable,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn float_index_collisions_keep_first_dictionary_key_positions() {
        let base = 1_i64 << 53;
        let mut map = TargetMap::default();
        map.insert(Key::Int(base + 1), Some("first".into()));
        map.insert(Key::Int(base), Some("second".into()));
        map.insert(Key::Int(base + 1), Some("updated first".into()));
        let map = map.finish();
        assert_eq!(
            map.get(&Key::Int(base + 1), false),
            Some(&Some("updated first".into()))
        );
        assert_eq!(
            map.get(&Key::float(base as f64), true),
            Some(&Some("second".into()))
        );
    }

    #[test]
    fn object_index_keeps_exact_large_integer_keys() {
        let base = 1_i64 << 53;
        let mut map = TargetMap::default();
        map.insert(Key::Int(base + 1), Some("integer".into()));
        map.insert(Key::Text("other".into()), None);
        let map = map.finish();
        assert_eq!(map.get(&Key::float(base as f64), true), None);
    }

    /// Set all three paths to compare a prepared profile without rerunning its
    /// Python producer. The expected file must be that producer's saved map.
    #[test]
    #[ignore = "requires explicit ClinVar summary, MANE and source-produced map inputs"]
    fn prepared_summary_map_matches_every_source_entry() {
        let path = |name| {
            std::path::PathBuf::from(
                std::env::var_os(name).expect("missing explicit map-test input"),
            )
        };
        let actual = from_tables(&path("PRIVA_GENE_SUMMARY"), &path("PRIVA_GENE_MANE")).unwrap();
        let expected = from_pickle(&path("PRIVA_GENE_REFERENCE_MAP")).unwrap();
        let missing = expected
            .values
            .keys()
            .filter(|key| !actual.values.contains_key(*key))
            .count();
        let extra = actual
            .values
            .keys()
            .filter(|key| !expected.values.contains_key(*key))
            .count();
        let changed = expected
            .values
            .iter()
            .filter(|(key, value)| {
                actual
                    .values
                    .get(*key)
                    .is_some_and(|actual| actual != *value)
            })
            .count();
        println!(
            "{{\"map_entries\":{},\"missing\":{missing},\"extra\":{extra},\"value_mismatches\":{changed}}}",
            expected.values.len()
        );
        assert_eq!(
            (missing, extra, changed),
            (0, 0, 0),
            "complete map comparison differs; only counts are reported"
        );
    }
}
