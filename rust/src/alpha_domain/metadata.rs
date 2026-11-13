//! Only the metadata extraction boundary in `am_pick_intolerant_domains.py`.
//! No DAS, KS, reference selection, random sampling or intolerance policy lives here.
use super::parquet::{self, Row, Writer};
use super::vcf::{self, Context};
use super::{Result, exception, source_error};
use crate::vcf_variant::{Definition, Record, Scalar};
use std::collections::{HashMap, HashSet};
use std::path::Path;

#[derive(Clone)]
pub(crate) struct Entry {
    pub symbol: String,
    pub start: i64,
    pub end: i64,
    pub length: i64,
}
#[derive(Default)]
pub(crate) struct Metadata {
    pub entries: Vec<(String, Entry)>,
    index: HashMap<String, usize>,
}
impl Metadata {
    pub fn append_unique(&mut self, key: String, entry: Entry) -> Result<()> {
        if self.index.contains_key(&key) {
            return Err(exception("ValueError", "duplicate domain metadata key"));
        }
        self.index.insert(key.clone(), self.entries.len());
        self.entries.push((key, entry));
        Ok(())
    }
    pub fn observe(&mut self, key: &str, symbol: &str, position: i64, length: Option<i64>) {
        if let Some(&at) = self.index.get(key) {
            let entry = &mut self.entries[at].1;
            entry.start = entry.start.min(position);
            entry.end = entry.end.max(position);
            if let Some(length) = length {
                if length > entry.length {
                    entry.length = length;
                }
            }
        } else {
            self.index.insert(key.into(), self.entries.len());
            self.entries.push((
                key.into(),
                Entry {
                    symbol: symbol.into(),
                    start: position,
                    end: position,
                    length: length.unwrap_or(0),
                },
            ));
        }
    }
    pub fn merge(&mut self, other: Self) -> Result<()> {
        for (key, entry) in other.entries {
            if let Some(&at) = self.index.get(&key) {
                let old = &mut self.entries[at].1;
                old.start = old.start.min(entry.start);
                old.end = old.end.max(entry.end);
                old.length = old.length.max(entry.length);
            } else {
                self.append_unique(key, entry)?;
            }
        }
        Ok(())
    }
}

pub(super) fn record(
    out: &mut Metadata,
    record: &Record<'_>,
    definitions: &HashMap<String, Definition>,
    context: &Context,
    target: Option<&HashSet<String>>,
) -> Result<()> {
    if !definitions.contains_key("CSQ") {
        return Ok(());
    }
    let Some(value) = record.get("CSQ", definitions).map_err(source_error)? else {
        return Ok(());
    };
    if value.values.is_empty()
        || (value.scalar && matches!(value.values.first(), Some(Scalar::Missing)))
    {
        return Ok(());
    }
    if value.scalar
        && (matches!(
            value.values.first(),
            Some(Scalar::Bool(false)) | Some(Scalar::Int(0))
        ) || matches!(value.values.first(), Some(Scalar::Float(number)) if *number == 0.0))
    {
        return Ok(());
    }
    for annotation in vcf::csqs(value)? {
        let fields: Vec<_> = annotation.split('|').collect();
        if context.field(&fields, "Feature_type")? != "Transcript"
            || !context.field(&fields, "Feature")?.starts_with("ENST")
        {
            continue;
        }
        let domains = context.field(&fields, "DOMAINS")?;
        let gene = context.field(&fields, "Gene")?;
        let symbol = context.field(&fields, "SYMBOL")?;
        let protein = context.field(&fields, "Protein_position")?;
        if domains.is_empty() || gene.is_empty() || protein.is_empty() {
            continue;
        }
        let Some(position) = vcf::position(protein)? else {
            continue;
        };
        let length = if let Some(length) = protein.split('/').nth(1) {
            match vcf::integer(length)? {
                Some(value) => Some(value),
                None => continue,
            }
        } else {
            None
        };
        for entry in domains.split('&') {
            let parts: Vec<_> = entry.split(':').collect();
            if ["PDB-ENSP_mappings", "AFDB-ENSP_mappings", "MobiDB_lite"].contains(&parts[0]) {
                continue;
            }
            for end in 1..parts.len() {
                let key = format!("{gene}:{}", parts[..=end].join(":"));
                if target.is_some_and(|set| !set.contains(&key)) {
                    continue;
                }
                out.observe(&key, symbol, position, length);
            }
        }
    }
    Ok(())
}
pub(super) fn write(writer: &mut Writer, metadata: &Metadata) -> Result<()> {
    for (ordinal, (key, entry)) in metadata.entries.iter().enumerate() {
        writer.push(Row {
            gene: key.clone(),
            gene_order: ordinal as u64,
            key: entry.symbol.clone(),
            kind: "metadata".into(),
            minimum_int: Some(entry.start),
            maximum: Some(entry.end),
            auxiliary_int: Some(entry.length),
            ..Row::default()
        })?;
    }
    Ok(())
}
pub(crate) fn read(path: &Path) -> Result<Metadata> {
    let mut output = Metadata::default();
    let meta = parquet::scan(path, |row| {
        if row.kind != "metadata"
            || row.gene_order != output.entries.len() as u64
            || !row.path.is_empty()
            || !row.values.is_empty()
            || !row.keys.is_empty()
            || !row.typed_keys.is_empty()
            || !row.counts.is_empty()
            || row.minimum_float.is_some()
        {
            return Err(exception("ValueError", "invalid metadata-only Parquet row"));
        }
        output.append_unique(
            row.gene,
            Entry {
                symbol: row.key,
                start: row
                    .minimum_int
                    .ok_or_else(|| exception("ValueError", "metadata start absent"))?,
                end: row
                    .maximum
                    .ok_or_else(|| exception("ValueError", "metadata end absent"))?,
                length: row
                    .auxiliary_int
                    .ok_or_else(|| exception("ValueError", "metadata protein length absent"))?,
            },
        )
    })?;
    if meta.get("priva.kind").map(String::as_str) != Some("metadata") {
        return Err(exception(
            "ValueError",
            "not a domain metadata intermediate",
        ));
    }
    Ok(output)
}
