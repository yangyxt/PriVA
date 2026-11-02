//! Gene-builder entry to the original whole-result collector. Ordinary sample
//! conversion keeps its existing serial entry point and parsing profile.
use super::*;
use crate::cache_convert::clinvar_gene_stat::Symbol;
use crate::vcf_variant;
use std::collections::HashSet;
use std::io::{Read, Seek, SeekFrom};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, mpsc};

pub(crate) const GENE_FIELDS: [&str; 11] = [
    "SYMBOL",
    "Gene",
    "BIOTYPE",
    "Feature",
    "DOMAINS",
    "Consequence",
    "LoF",
    "NMD",
    "SpliceVault_top_events",
    "SpliceVault_SpliceAI_delta",
    "SpliceAI_pred_DS_DL",
];

enum FinalColumn {
    Position(Vec<i64>),
    Numeric { values: Vec<f32>, whole: bool },
    Objects(Vec<Object>),
}

/// Values after pandas' final nullable Float32/Int32 conversion. Object floats
/// retain their binary64 bits, including earlier chunk-specific rounding.
#[derive(Clone, Copy, Debug)]
pub(crate) enum Scalar<'a> {
    PandasNA,
    None,
    Int64(i64),
    Int32(i32),
    Float32(f32),
    Float64(f64),
    Text(&'a str),
}
impl Scalar<'_> {
    pub(crate) fn missing(self) -> bool {
        matches!(self, Self::PandasNA | Self::None)
            || matches!(self, Self::Float64(value) if value.is_nan())
    }
    pub(crate) fn text(self) -> String {
        match self {
            Self::PandasNA => "<NA>".into(),
            Self::None => "None".into(),
            Self::Int64(value) => value.to_string(),
            Self::Int32(value) => value.to_string(),
            Self::Float32(value) => numpy_f32_repr(value),
            Self::Float64(value) => python_repr(value),
            Self::Text(value) => value.to_owned(),
        }
    }
    pub(crate) fn text_or_nan(self) -> String {
        if self.missing() {
            "nan".into()
        } else {
            self.text()
        }
    }
    pub(crate) fn number_or_nan(self) -> Result<f64> {
        if self.missing() {
            return Ok(f64::NAN);
        }
        match self {
            Self::Int64(value) => Ok(value as f64),
            Self::Int32(value) => Ok(value as f64),
            Self::Float32(value) => Ok(value as f64),
            Self::Float64(value) => Ok(value),
            Self::Text(value) => {
                py_float(value).ok_or_else(|| invalid("ValueError: a gene score is not numeric"))
            }
            _ => unreachable!(),
        }
    }
    pub(crate) fn symbol(self) -> Result<Symbol> {
        Ok(match self {
            Self::PandasNA => Symbol::PandasNA,
            Self::Int32(value) => Symbol::NumpyI32(value),
            Self::Float32(value) => Symbol::NumpyF32(value.to_bits()),
            Self::Float64(value) => Symbol::Float(value.to_bits()),
            Self::Text(value) => Symbol::Text(value.to_owned()),
            _ => {
                return Err(invalid(
                    "the source gene symbol has no supported schema-5 representation",
                ));
            }
        })
    }
}

pub(crate) struct Table {
    index: HashMap<String, usize>,
    columns: Vec<FinalColumn>,
    pub(crate) rows: usize,
}
impl Table {
    pub(crate) fn column(&self, name: &str) -> Result<usize> {
        self.index
            .get(name)
            .copied()
            .ok_or_else(|| invalid(format!("KeyError: gene annotation column {name}")))
    }
    pub(crate) fn cell(&self, column: usize, row: usize) -> Scalar<'_> {
        match &self.columns[column] {
            FinalColumn::Position(values) => Scalar::Int64(values[row]),
            FinalColumn::Numeric { values, whole } => {
                let value = values[row];
                if value.is_nan() {
                    Scalar::PandasNA
                } else if *whole {
                    Scalar::Int32(wrap_int32(value) as i32)
                } else {
                    Scalar::Float32(value)
                }
            }
            FinalColumn::Objects(values) => match &values[row] {
                Object::None => Scalar::None,
                Object::Float(value) => Scalar::Float64(*value),
                Object::Str(value) => Scalar::Text(value),
            },
        }
    }
    pub(crate) fn get(&self, name: &str, row: usize) -> Result<Scalar<'_>> {
        Ok(self.cell(self.column(name)?, row))
    }
    pub(crate) fn numeric_column(&self, name: &str) -> Result<bool> {
        Ok(!matches!(
            self.columns[self.column(name)?],
            FinalColumn::Objects(_)
        ))
    }
}

struct Collector {
    table: Option<Annotations>,
    chunk: Chunk,
}
impl Collector {
    fn new() -> Self {
        Self {
            table: None,
            chunk: Chunk {
                keys: vec![],
                values: vec![],
                numeric: vec![],
                rows: 0,
            },
        }
    }
    fn push(&mut self, rows: Vec<Row>) -> Result<()> {
        if rows.is_empty() {
            return Ok(());
        }
        let table = self.table.get_or_insert_with(|| {
            let names: Vec<_> = rows[0].items.iter().map(|(key, _)| key.clone()).collect();
            Annotations {
                index: names
                    .iter()
                    .enumerate()
                    .map(|(at, name)| (name.clone(), at))
                    .collect(),
                columns: names
                    .iter()
                    .map(|name| {
                        if name == "pos" {
                            Column::Position(vec![])
                        } else {
                            Column::Objects(vec![])
                        }
                    })
                    .collect(),
                names,
                count: 0,
            }
        });
        for row in rows {
            if self.chunk.rows == 0 {
                self.chunk.keys = row.items.iter().map(|(key, _)| key.clone()).collect();
                self.chunk.values = vec![vec![]; self.chunk.keys.len()];
                self.chunk.numeric = vec![false; self.chunk.keys.len()];
            }
            for (at, key) in self.chunk.keys.iter().enumerate() {
                let (value, numeric) = collect(key == "pos", row.get(key), key)?;
                self.chunk.values[at].push(value);
                self.chunk.numeric[at] |= numeric;
            }
            self.chunk.rows += 1;
        }
        // The original checks the threshold after each complete worker result.
        if self.chunk.rows >= CHUNK_SIZE {
            table.store(&mut self.chunk)?;
        }
        Ok(())
    }
    fn finish(self) -> Result<Table> {
        let mut table = self
            .table
            .ok_or_else(|| invalid("no record has an Ensembl feature in CSQ"))?;
        let mut chunk = self.chunk;
        if chunk.rows > 0 {
            table.store(&mut chunk)?;
        }
        let mut columns = Vec::with_capacity(table.columns.len());
        for (name, mut column) in table.names.iter().zip(table.columns) {
            column.fill_to(table.count);
            let values = match column {
                Column::Position(values) if name == "pos" => {
                    columns.push(FinalColumn::Position(values));
                    continue;
                }
                Column::Position(values) => values.into_iter().map(|value| value as f32).collect(),
                Column::Single(values) => values,
                Column::Objects(values) => {
                    if values.iter().any(|value| matches!(value, Object::Str(_))) {
                        columns.push(FinalColumn::Objects(values));
                        continue;
                    }
                    values
                        .into_iter()
                        .map(|value| match value {
                            Object::Float(value) => value as f32,
                            _ => f32::NAN,
                        })
                        .collect()
                }
            };
            let whole = values
                .iter()
                .filter(|value| !value.is_nan())
                .all(|value| value.is_finite() && value.fract() == 0.0);
            columns.push(FinalColumn::Numeric { values, whole });
        }
        Ok(Table {
            index: table.index,
            columns,
            rows: table.count,
        })
    }
}

struct Task {
    ordinal: usize,
    line: String,
    header: Arc<VcfHeader>,
}

fn checked_header(path: &Path) -> Result<(vcf_variant::Header, Arc<VcfHeader>, Vec<String>)> {
    let checked = vcf_variant::header(path)?;
    let mut file = File::open(path)?;
    let mut magic = [0; 2];
    if file.read_exact(&mut magic).is_ok() && magic == [0x1f, 0x8b] {
        let mut tail = [0; 28];
        if file.seek(SeekFrom::End(-28)).is_err()
            || file.read_exact(&mut tail).is_err()
            || tail != crate::vcf_text::EOF_BLOCK
        {
            return Err(invalid("OSError: BGZF EOF marker is absent"));
        }
    }
    let mut infos = HashMap::new();
    let mut samples = Vec::new();
    for line in crate::vcf_text::open_reader(path)?.lines() {
        let line = line?;
        if line.starts_with("#CHROM") {
            samples = line.split('\t').skip(9).map(str::to_owned).collect();
            break;
        }
        if let Some(rest) = line.strip_prefix("##INFO=") {
            let attrs = structured(rest);
            let get = |name: &str| attrs.get(name).cloned().unwrap_or_default();
            infos.entry(get("ID")).or_insert(Definition {
                number: get("Number"),
                kind: get("Type"),
                description: get("Description"),
            });
        }
    }
    let header = Arc::new(VcfHeader { infos, samples });
    let declaration = header
        .infos
        .get("CSQ")
        .and_then(|value| value.description.split_once("Format: "))
        .ok_or_else(|| invalid("ClinVar INFO/CSQ lacks its Format declaration"))?
        .1;
    let fields: Vec<_> = declaration.split('|').collect();
    if fields.len() != fields.iter().collect::<HashSet<_>>().len()
        || GENE_FIELDS.iter().any(|name| !fields.contains(name))
    {
        return Err(invalid("ClinVar gene CSQ fields are missing or duplicated"));
    }
    let mut required: Vec<String> = [
        "CSQ",
        "CLNSIG",
        "CLNREVSTAT",
        "ALLELEID",
        "AC_joint",
        "AN_joint",
        "AF_joint",
        "nhomalt_joint",
        "AC_grpmax_joint",
        "AF_grpmax_joint",
        "AN_grpmax_joint",
        "nhomalt_grpmax_joint",
    ]
    .iter()
    .map(|s| (*s).into())
    .collect();
    for pop in POPULATIONS {
        for suffix in [
            format!("AF_joint_{pop}_XX"),
            format!("AF_joint_{pop}_XY"),
            format!("nhomalt_joint_{pop}_XX"),
        ] {
            required.push(suffix);
        }
    }
    if required.iter().any(|name| !header.infos.contains_key(name)) {
        return Err(invalid(
            "ClinVar gene annotation lacks required INFO fields",
        ));
    }
    let csq = header.csq_fields("CSQ")?;
    if header.infos.contains_key("CLNCSQ") {
        header.csq_fields("CLNCSQ")?;
    }
    Ok((checked, header, csq))
}

fn tasks(
    path: &Path,
    checked: &vcf_variant::Header,
    mut header: Arc<VcfHeader>,
    cancelled: &AtomicBool,
    mut emit: impl FnMut(Task) -> Result<()>,
) -> Result<usize> {
    let mut info = checked.info.clone();
    let mut body = false;
    let mut count = 0;
    for line in crate::vcf_text::open_reader(path)?.lines() {
        if cancelled.load(Ordering::Relaxed) {
            break;
        }
        let line = line?;
        if !body {
            if line.starts_with("#CHROM") {
                body = true;
            }
            continue;
        }
        let record = vcf_variant::decode_record(line.as_bytes(), checked, &mut info)?;
        for key in record.info.keys() {
            if !header.infos.contains_key(*key) {
                Arc::make_mut(&mut header).infos.insert(
                    (*key).to_owned(),
                    Definition {
                        number: "1".into(),
                        kind: "String".into(),
                        description: "Dummy".into(),
                    },
                );
            }
        }
        let mut fields: Vec<String> = record
            .fields
            .iter()
            .map(|value| (*value).to_owned())
            .collect();
        if fields.len() < 5 {
            return Err(invalid("invalid VCF record"));
        }
        fields.resize(fields.len().max(8), ".".into());
        fields[1] = record.position.to_string();
        emit(Task {
            ordinal: count,
            line: fields.join("\t"),
            header: Arc::clone(&header),
        })?;
        count += 1;
    }
    Ok(count)
}

pub(crate) fn read(
    path: &Path,
    threads: usize,
    order: Option<&[usize]>,
) -> Result<(Table, Vec<usize>)> {
    if threads == 0 {
        return Err(invalid("the worker count must be positive"));
    }
    if let Some(order) = order {
        if order.iter().copied().collect::<HashSet<_>>().len() != order.len()
            || order.iter().any(|at| *at >= order.len())
        {
            return Err(invalid(
                "completion order must be a permutation of zero-based record ordinals",
            ));
        }
    }
    let (checked, header, csq) = checked_header(path)?;
    let feature = csq.iter().position(|name| name == "Feature");
    let source = header.infos.contains_key("VARIANT_SOURCE");
    let control = header.infos.contains_key("control_AC");
    let cadd = CaddScores::new();
    let hpo = HpoAnnotations::new();
    let worker = |task: &Task| {
        record_rows_profile(
            &task.line,
            &task.header,
            &csq,
            feature,
            source,
            control,
            &cadd,
            &hpo,
            true,
        )
    };
    let mut collector = Collector::new();
    let mut trace = Vec::new();
    let mut pending = HashMap::new();
    let mut accept = |ordinal: usize, rows: Vec<Row>| -> Result<()> {
        if let Some(order) = order {
            if ordinal >= order.len() {
                return Err(invalid("completion trace has fewer records than the VCF"));
            }
            pending.insert(ordinal, rows);
            while let Some(&next) = order.get(trace.len()) {
                let Some(rows) = pending.remove(&next) else {
                    break;
                };
                collector.push(rows)?;
                trace.push(next);
            }
        } else {
            collector.push(rows)?;
            trace.push(ordinal);
        }
        Ok(())
    };
    let cancelled = AtomicBool::new(false);
    let count = if threads == 1 {
        tasks(path, &checked, header, &cancelled, |task| {
            accept(task.ordinal, worker(&task)?)
        })?
    } else {
        std::thread::scope(|scope| -> Result<usize> {
            let (input, receive) = mpsc::sync_channel::<Task>(threads * 2);
            let receive = Arc::new(Mutex::new(receive));
            let (output, results) = mpsc::sync_channel(threads * 2);
            let cancelled = &cancelled;
            let worker = &worker;
            for _ in 0..threads {
                let receive = Arc::clone(&receive);
                let output = output.clone();
                scope.spawn(move || {
                    while !cancelled.load(Ordering::Relaxed) {
                        let Ok(task) = receive.lock().unwrap().recv() else {
                            break;
                        };
                        let result = worker(&task).map_err(|error| error.to_string());
                        if output.send((task.ordinal, result)).is_err() {
                            break;
                        }
                    }
                });
            }
            drop(receive);
            drop(output);
            let reader = scope.spawn(move || {
                tasks(path, &checked, header, cancelled, |task| {
                    input
                        .send(task)
                        .map_err(|_| invalid("annotation workers stopped"))
                })
                .map_err(|error| error.to_string())
            });
            let mut error = None;
            for (ordinal, result) in results {
                if error.is_some() {
                    continue;
                }
                if let Err(failure) = result
                    .map_err(invalid)
                    .and_then(|rows| accept(ordinal, rows))
                {
                    cancelled.store(true, Ordering::Relaxed);
                    error = Some(failure);
                }
            }
            let produced = reader
                .join()
                .map_err(|_| invalid("VCF reader thread panicked"))?;
            if let Some(error) = error {
                return Err(error);
            }
            produced.map_err(invalid)
        })?
    };
    if order.is_some_and(|order| order.len() != count)
        || trace.len() != count
        || !pending.is_empty()
    {
        return Err(invalid("completion trace and VCF record counts differ"));
    }
    Ok((collector.finish()?, trace))
}
