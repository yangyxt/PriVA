//! Typed, ordinal-preserving Parquet intermediate. Consecutive residue samples
//! use explicit positive counts beside exact Float64 values. All other lists
//! retain their item order; empty lists are distinct from nullable scalar fields.
use super::{Result, exception};
use super::{key::Key, key_parquet};
use crate::cli::invalid;
use crate::sha256::Sha256;
use crate::vcf_variant::Scalar;
use polars_core::prelude::*;
use polars_io::parquet::write::{BatchedWriter, KeyValueMetadata};
use polars_io::prelude::*;
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

#[cfg(test)]
#[path = "../../tests/support/alpha_readback_probe.rs"]
mod readback_probe;

pub(super) const SCHEMA: &str = "priva-alpha-domain-parquet-v2";
#[derive(Default, Clone, Debug)]
pub(super) struct Row {
    pub gene: String,
    pub gene_order: u64,
    pub path: Vec<String>,
    pub member_order: u64,
    pub item_order: u64,
    pub key: String,
    pub kind: String,
    pub keys: Vec<String>,
    pub typed_keys: Vec<Key>,
    pub values: Vec<f64>,
    pub counts: Vec<u64>,
    pub minimum_int: Option<i64>,
    pub minimum_float: Option<f64>,
    pub maximum: Option<i64>,
    pub auxiliary_int: Option<i64>,
}
fn strings(name: &str, values: impl Iterator<Item = Vec<String>>) -> Series {
    let values: Vec<_> = values
        .map(|values| Series::new("".into(), values))
        .collect();
    if values.is_empty() {
        Series::new_empty(name.into(), &DataType::List(Box::new(DataType::String)))
    } else {
        Series::new(name.into(), values)
    }
}
fn frame(rows: &[Row]) -> Result<DataFrame> {
    let numbers = |name: &str, items: Vec<Series>, dtype| {
        if items.is_empty() {
            Series::new_empty(name.into(), &DataType::List(Box::new(dtype)))
        } else {
            Series::new(name.into(), items)
        }
    };
    Ok(DataFrame::new(vec![
        Series::new(
            "gene".into(),
            rows.iter().map(|r| r.gene.as_str()).collect::<Vec<_>>(),
        )
        .into(),
        Series::new(
            "gene_order".into(),
            rows.iter().map(|r| r.gene_order).collect::<Vec<_>>(),
        )
        .into(),
        strings("path", rows.iter().map(|r| r.path.clone())).into(),
        Series::new(
            "member_order".into(),
            rows.iter().map(|r| r.member_order).collect::<Vec<_>>(),
        )
        .into(),
        Series::new(
            "item_order".into(),
            rows.iter().map(|r| r.item_order).collect::<Vec<_>>(),
        )
        .into(),
        Series::new(
            "key".into(),
            rows.iter().map(|r| r.key.as_str()).collect::<Vec<_>>(),
        )
        .into(),
        Series::new(
            "kind".into(),
            rows.iter().map(|r| r.kind.as_str()).collect::<Vec<_>>(),
        )
        .into(),
        strings("keys", rows.iter().map(|r| r.keys.clone())).into(),
        key_parquet::column(rows.iter().map(|r| r.typed_keys.as_slice()))?.into(),
        numbers(
            "values",
            rows.iter()
                .map(|r| Series::new("".into(), &r.values))
                .collect(),
            DataType::Float64,
        )
        .into(),
        numbers(
            "counts",
            rows.iter()
                .map(|r| Series::new("".into(), &r.counts))
                .collect(),
            DataType::UInt64,
        )
        .into(),
        Series::new(
            "minimum_int".into(),
            rows.iter().map(|r| r.minimum_int).collect::<Vec<_>>(),
        )
        .into(),
        Series::new(
            "minimum_float".into(),
            rows.iter().map(|r| r.minimum_float).collect::<Vec<_>>(),
        )
        .into(),
        Series::new(
            "maximum".into(),
            rows.iter().map(|r| r.maximum).collect::<Vec<_>>(),
        )
        .into(),
        Series::new(
            "auxiliary_int".into(),
            rows.iter().map(|r| r.auxiliary_int).collect::<Vec<_>>(),
        )
        .into(),
    ])?)
}

fn digest(row: &Row, sha: &mut Sha256) {
    fn text(value: &str, sha: &mut Sha256) {
        sha.update(&(value.len() as u64).to_le_bytes());
        sha.update(value.as_bytes());
    }
    text(&row.gene, sha);
    sha.update(&row.gene_order.to_le_bytes());
    sha.update(&(row.path.len() as u64).to_le_bytes());
    for value in &row.path {
        text(value, sha);
    }
    sha.update(&row.member_order.to_le_bytes());
    sha.update(&row.item_order.to_le_bytes());
    text(&row.key, sha);
    text(&row.kind, sha);
    sha.update(&(row.keys.len() as u64).to_le_bytes());
    for value in &row.keys {
        text(value, sha);
    }
    sha.update(&(row.typed_keys.len() as u64).to_le_bytes());
    for key in &row.typed_keys {
        match key {
            Key::Text(value) => {
                sha.update(&[0]);
                text(value, sha);
            }
            Key::Tuple(values) => {
                sha.update(&[1]);
                sha.update(&(values.len() as u64).to_le_bytes());
                for value in values.iter() {
                    match value {
                        Scalar::Missing => sha.update(&[0]),
                        Scalar::Text(v) => {
                            sha.update(&[1]);
                            text(v, sha);
                        }
                        Scalar::Int(v) => {
                            sha.update(&[2]);
                            sha.update(&v.to_le_bytes());
                        }
                        Scalar::Float(v) => {
                            sha.update(&[3]);
                            sha.update(&v.to_bits().to_le_bytes());
                        }
                        Scalar::Bool(v) => sha.update(&[4, u8::from(*v)]),
                    }
                }
            }
        }
    }
    sha.update(&(row.values.len() as u64).to_le_bytes());
    for value in &row.values {
        sha.update(&value.to_bits().to_le_bytes());
    }
    sha.update(&(row.counts.len() as u64).to_le_bytes());
    for value in &row.counts {
        sha.update(&value.to_le_bytes());
    }
    for value in [
        row.minimum_int.map(|n| n as u64),
        row.minimum_float.map(f64::to_bits),
        row.maximum.map(|n| n as u64),
        row.auxiliary_int.map(|n| n as u64),
    ] {
        sha.update(&[u8::from(value.is_some())]);
        if let Some(value) = value {
            sha.update(&value.to_le_bytes());
        }
    }
}

pub(super) struct Writer {
    writer: BatchedWriter<File>,
    rows: Vec<Row>,
    bytes: usize,
    sha: Sha256,
    count: u64,
}
impl Writer {
    fn new(file: File, mut metadata: Vec<(String, String)>) -> Result<Self> {
        metadata.insert(0, ("priva.schema".into(), SCHEMA.into()));
        let writer = ParquetWriter::new(file)
            .with_compression(ParquetCompression::Zstd(None))
            .set_parallel(false)
            .with_key_value_metadata(Some(KeyValueMetadata::from_static(metadata)))
            .batched(frame(&[])?.schema())?;
        Ok(Self {
            writer,
            rows: vec![],
            bytes: 0,
            sha: Sha256::new(),
            count: 0,
        })
    }
    pub fn push(&mut self, row: Row) -> Result<()> {
        digest(&row, &mut self.sha);
        self.count += 1;
        self.bytes += row.values.len() * 8
            + row.counts.len() * 8
            + row.keys.iter().map(String::len).sum::<usize>()
            + row.typed_keys.iter().map(Key::bytes).sum::<usize>()
            + row.path.iter().map(String::len).sum::<usize>()
            + row.gene.len()
            + row.key.len()
            + 128;
        self.rows.push(row);
        if self.rows.len() >= 16_384 || self.bytes >= 8 * 1024 * 1024 {
            self.flush()?;
        }
        Ok(())
    }
    fn flush(&mut self) -> Result<()> {
        if !self.rows.is_empty() {
            self.writer.write_batch(&frame(&self.rows)?)?;
            self.rows.clear();
            self.bytes = 0;
        }
        Ok(())
    }
    fn finish(mut self) -> Result<(u64, String)> {
        self.flush()?;
        self.writer.finish()?;
        Ok((self.count, self.sha.finish_hex()))
    }
}

/// Offsets and the flat child of one list column. A Series built for every
/// cell made whole-genome scans take hours.
pub(super) fn list_parts(
    series: &Series,
    null: impl FnOnce() -> Box<dyn std::error::Error>,
) -> Result<(Vec<usize>, Series)> {
    let series = series.rechunk();
    let lists = series.list()?;
    if lists.null_count() != 0 {
        return Err(null());
    }
    let mut arrays = lists.downcast_iter();
    let offsets: Vec<usize> = match (arrays.next(), arrays.next()) {
        (Some(array), None) => array.offsets().iter().map(|&at| at as usize).collect(),
        (None, None) => vec![0],
        _ => return Err(invalid("Parquet list column is not contiguous")),
    };
    let child = lists.get_inner();
    if offsets.len() != series.len() + 1 || offsets.last() > Some(&child.len()) {
        return Err(invalid("Parquet list offsets differ from their values"));
    }
    Ok((offsets, child))
}

fn span(offsets: &[usize], at: usize) -> std::ops::Range<usize> {
    offsets[at]..offsets[at + 1]
}

/// Columns are resolved once per row group; each row copies its own cells.
fn rows(frame: &DataFrame, mut visit: impl FnMut(Row) -> Result<()>) -> Result<()> {
    let column = |name: &str| -> Result<Series> {
        Ok(frame.column(name)?.as_materialized_series().rechunk())
    };
    let list = |name: &str| -> Result<(Vec<usize>, Series)> {
        list_parts(&column(name)?, || invalid("unexpected null Parquet list"))
    };
    let items = |child: &Series| -> Result<Vec<String>> {
        child
            .str()?
            .into_iter()
            .map(|value| {
                value
                    .map(str::to_owned)
                    .ok_or_else(|| invalid("unexpected null Parquet list item"))
            })
            .collect()
    };
    let gene = column("gene")?;
    let gene = gene.str()?;
    let key = column("key")?;
    let key = key.str()?;
    let kind = column("kind")?;
    let kind = kind.str()?;
    let gene_order = column("gene_order")?;
    let gene_order = gene_order.u64()?;
    let member_order = column("member_order")?;
    let member_order = member_order.u64()?;
    let item_order = column("item_order")?;
    let item_order = item_order.u64()?;
    let minimum_int = column("minimum_int")?;
    let minimum_int = minimum_int.i64()?;
    let minimum_float = column("minimum_float")?;
    let minimum_float = minimum_float.f64()?;
    let maximum = column("maximum")?;
    let maximum = maximum.i64()?;
    let auxiliary_int = column("auxiliary_int")?;
    let auxiliary_int = auxiliary_int.i64()?;
    let (path_at, path) = list("path")?;
    let path = items(&path)?;
    let (keys_at, keys) = list("keys")?;
    let keys = items(&keys)?;
    let (typed_at, typed) = list("typed_keys")?;
    let typed = key_parquet::decode(typed)?;
    let (values_at, values) = list("values")?;
    let values = values
        .f64()?
        .into_iter()
        .map(|v| v.ok_or_else(|| invalid("null sample")))
        .collect::<Result<Vec<_>>>()?;
    let (counts_at, counts) = list("counts")?;
    let counts = counts
        .u64()?
        .into_iter()
        .map(|v| v.ok_or_else(|| invalid("null repeat count")))
        .collect::<Result<Vec<_>>>()?;
    let text = |column: &StringChunked, at| -> Result<String> {
        Ok(column
            .get(at)
            .ok_or_else(|| invalid("unexpected null Parquet text"))?
            .to_owned())
    };
    let integer = |column: &UInt64Chunked, at| -> Result<u64> {
        column
            .get(at)
            .ok_or_else(|| invalid("unexpected null Parquet ordinal"))
    };
    for at in 0..frame.height() {
        visit(Row {
            gene: text(gene, at)?,
            gene_order: integer(gene_order, at)?,
            path: path[span(&path_at, at)].to_vec(),
            member_order: integer(member_order, at)?,
            item_order: integer(item_order, at)?,
            key: text(key, at)?,
            kind: text(kind, at)?,
            keys: keys[span(&keys_at, at)].to_vec(),
            typed_keys: typed[span(&typed_at, at)].to_vec(),
            values: values[span(&values_at, at)].to_vec(),
            counts: counts[span(&counts_at, at)].to_vec(),
            minimum_int: minimum_int.get(at),
            minimum_float: minimum_float.get(at),
            maximum: maximum.get(at),
            auxiliary_int: auxiliary_int.get(at),
        })?;
    }
    Ok(())
}

pub(super) fn metadata(path: &Path) -> Result<HashMap<String, String>> {
    let mut reader = ParquetReader::new(File::open(path)?);
    let metadata = reader.get_metadata()?.clone();
    let values: HashMap<_, _> = metadata
        .key_value_metadata
        .as_deref()
        .unwrap_or_default()
        .iter()
        .filter_map(|v| v.value.as_ref().map(|value| (v.key.clone(), value.clone())))
        .collect();
    if values.get("priva.schema").map(String::as_str) != Some(SCHEMA) {
        return Err(exception(
            "ValueError",
            "unsupported AlphaMissense Parquet schema",
        ));
    }
    Ok(values)
}

pub(super) fn scan(
    path: &Path,
    mut visit: impl FnMut(Row) -> Result<()>,
) -> Result<HashMap<String, String>> {
    let values = metadata(path)?;
    let mut reader = ParquetReader::new(File::open(path)?);
    let metadata = reader.get_metadata()?.clone();
    if metadata.row_groups.is_empty() {
        let empty = reader.finish()?;
        if empty.schema().as_ref() != frame_schema()?.as_ref() || empty.height() != 0 {
            return Err(exception(
                "ValueError",
                "invalid empty AlphaMissense Parquet schema",
            ));
        }
    }
    let mut offset = 0;
    for group in &metadata.row_groups {
        let count = group.num_rows();
        // polars-io 0.51's None/Columns path also deserializes preceding row
        // groups for an empty slice, retaining their empty frames. RowGroups
        // skips nonoverlapping groups. This exact slice intersects one group,
        // so its rows keep source order without prefix work on every read.
        let mut reader = ParquetReader::new(File::open(path)?)
            .with_slice(Some((offset, count)))
            .read_parallel(ParallelStrategy::RowGroups);
        reader.set_metadata(metadata.clone());
        let frame = reader.finish()?;
        if frame.schema().as_ref() != frame_schema()?.as_ref() || frame.height() != count {
            return Err(exception(
                "ValueError",
                "invalid typed AlphaMissense Parquet row group",
            ));
        }
        rows(&frame, &mut visit)?;
        offset += count;
    }
    Ok(values)
}
fn frame_schema() -> Result<SchemaRef> {
    Ok(frame(&[])?.schema().clone())
}

/// One staged file, explicit Parquet finish, fsync, exact typed readback, then
/// rename. The previous destination survives every prepublication failure.
pub(super) fn publish(
    path: &Path,
    metadata: Vec<(String, String)>,
    emit: impl FnOnce(&mut Writer) -> Result<()>,
) -> Result<u64> {
    atomic_file(path, |file, temporary| {
        let sync = file.try_clone()?;
        let mut writer = Writer::new(file, metadata)?;
        emit(&mut writer)?;
        let (count, expected) = writer.finish()?;
        sync.sync_all()?;
        let mut sha = Sha256::new();
        let mut rows = 0;
        scan(temporary, |row| {
            digest(&row, &mut sha);
            rows += 1;
            Ok(())
        })?;
        if count != rows || sha.finish_hex() != expected {
            return Err(invalid("Parquet typed readback differs from source values"));
        }
        Ok(count)
    })
}

/// The metadata association uses the same staged-file ownership and explicit
/// finish boundary as its Parquet payload. Existing destination modes survive.
pub(super) fn atomic_file<T>(
    path: &Path,
    write: impl FnOnce(File, &Path) -> Result<T>,
) -> Result<T> {
    let permissions = match std::fs::metadata(path) {
        Ok(value) => Some(value.permissions()),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => None,
        Err(error) => return Err(error.into()),
    };
    static SERIAL: AtomicU64 = AtomicU64::new(0);
    let name = path
        .file_name()
        .ok_or_else(|| invalid("output needs a filename"))?
        .to_string_lossy();
    let mut selected = None;
    for _ in 0..128 {
        let candidate = path.with_file_name(format!(
            ".{name}.{}.{}.tmp",
            std::process::id(),
            SERIAL.fetch_add(1, Ordering::Relaxed)
        ));
        match OpenOptions::new()
            .read(true)
            .write(true)
            .create_new(true)
            .open(&candidate)
        {
            Ok(file) => {
                selected = Some((candidate, file));
                break;
            }
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(e.into()),
        }
    }
    let (temporary, file) = selected.ok_or_else(|| invalid("cannot reserve Parquet output"))?;
    let result = (|| {
        let value = write(file, &temporary)?;
        if let Some(permissions) = permissions {
            std::fs::set_permissions(&temporary, permissions)?;
        }
        std::fs::rename(&temporary, path)?;
        Ok(value)
    })();
    if result.is_err() {
        let _ = std::fs::remove_file(&temporary);
    }
    result
}
