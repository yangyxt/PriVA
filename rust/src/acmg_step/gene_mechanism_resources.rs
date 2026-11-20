//! `scripts/gene_mechanism_resources.py` (1,028 lines): `GeneMechanismHub` and the local HGNC
//! resolver.
//!
//! Thread T2 (gene-condition hub). Converted, what the ACMG step uses (it builds the hub with
//! `use_hgnc_package=False`): `_LocalHgncResolver` (`HgncResolver`), the hub's `resolve_symbol`
//! and `_try_resolve_symbol` (the `types::SymbolResolver` the exact-allele matcher calls),
//! `_load_condition_cache`, `_load_clingen` and `_load_loeuf`. `_load_mechanisms` is the
//! matcher's (`gene_mechanism_variants::mechanisms_by_symbol`). The DDG2P, HPO and diagnostic
//! summaries (`mechanism_history`, `known_inheritance_mode`, `gene_summary`, ...) serve group 4.
//!
//! The HGNC table is read as `csv.DictReader` reads it (tab-separated). The shipped table has no
//! quote characters, carriage returns or short rows; a table with any of them is refused rather
//! than parsed with the csv module's quoting rules. The ClinGen and LOEUF tables are read as
//! `pd.read_csv(dtype=str)` reads them (`pandas_table::read_table_with`, `all_text`); the
//! `.fillna("")` that follows is `unwrap_or_default`.
use super::clinvar_vcv::open_text;
use super::gene_mechanism_common::{
    DEFAULT_CLINGEN_DOSAGE, DEFAULT_HGNC_TABLE, DEFAULT_HPO_CONDITION_MECHANISM_CACHE,
    DEFAULT_LOEUF_TABLE, HPO_CONDITION_MECHANISM_SCHEMA_VERSION, LOOKUP_FIELD_PRIORITY, clean,
    clean_json, clean_text, default_data_path, python_error, python_path_text, python_text_repr,
    refused, safe_float, safe_int, show_outcome, show_text, split_multi,
};
use super::types::{HarnessCall, SymbolResolver, not_converted};
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::pandas_table::{Cell, Column, ReadOptions, Table, UseColumns, read_table_with};
use std::collections::HashMap;
use std::path::Path;
use std::rc::Rc;

/// The csv module's default field size limit.
const CSV_FIELD_LIMIT: usize = 131_072;

/// A dict in Python's insertion order: a key keeps its first place.
#[derive(Clone, Debug)]
pub(crate) struct OrderedMap<V> {
    entries: Vec<(String, V)>,
    places: HashMap<String, usize>,
}

impl<V> Default for OrderedMap<V> {
    fn default() -> Self {
        OrderedMap {
            entries: Vec::new(),
            places: HashMap::new(),
        }
    }
}

impl<V> OrderedMap<V> {
    /// `d[key] = value`.
    pub(crate) fn insert(&mut self, key: String, value: V) {
        match self.places.get(&key) {
            Some(place) => self.entries[*place].1 = value,
            None => {
                self.places.insert(key.clone(), self.entries.len());
                self.entries.push((key, value));
            }
        }
    }

    /// `d.setdefault(key, value)`.
    pub(crate) fn set_default(&mut self, key: String, value: V) {
        if !self.places.contains_key(&key) {
            self.places.insert(key.clone(), self.entries.len());
            self.entries.push((key, value));
        }
    }

    pub(crate) fn get(&self, key: &str) -> Option<&V> {
        self.places.get(key).map(|place| &self.entries[*place].1)
    }

    pub(crate) fn contains_key(&self, key: &str) -> bool {
        self.places.contains_key(key)
    }

    pub(crate) fn entries(&self) -> &[(String, V)] {
        &self.entries
    }
}

/// `open(path)`'s error for a missing file, or a refusal for other failures.
fn open_error(path: &str, error: std::io::Error) -> Box<dyn std::error::Error> {
    if error.kind() == std::io::ErrorKind::NotFound {
        match python_text_repr(&python_path_text(path)) {
            Ok(shown) => python_error(
                "FileNotFoundError",
                format!("[Errno 2] No such file or directory: {shown}"),
            ),
            Err(refusal) => refusal,
        }
    } else {
        refused(format!("reading {path}: {error}"))
    }
}

// ------------------------------------------------------------------------- the resolver

/// `_LocalHgncResolver`: the local HGNC table, looked up by symbol, HGNC id,
/// Ensembl and Entrez ids, previous and alias symbols, RefSeq, UniProt and MANE
/// ids, best field first.
#[derive(Clone, Debug, Default)]
pub(crate) struct HgncResolver {
    /// The normalised value (`_norm`) to the (field priority, row) pairs that hold it.
    lookup: HashMap<String, Vec<(u8, u32)>>,
    /// Each row's `symbol` (empty when the table has no `symbol` column).
    symbols: Vec<String>,
}

impl HgncResolver {
    /// `_LocalHgncResolver(hgnc_table)`: reads the table.
    pub(crate) fn load(path: &str) -> Result<HgncResolver> {
        let bytes = std::fs::read(path).map_err(|error| open_error(path, error))?;
        let text = String::from_utf8(bytes).map_err(|_| refused(format!("{path} is not UTF-8")))?;
        Self::parse(&text, path)
    }

    fn parse(text: &str, path: &str) -> Result<HgncResolver> {
        if text.contains(['"', '\r', '\0']) {
            return Err(refused(format!(
                "{path} holds quote characters, carriage returns or NUL bytes, which the csv module reads by its own rules"
            )));
        }
        let mut resolver = HgncResolver::default();
        let mut lines = text.split('\n');
        let Some(header) = lines.next().filter(|line| !line.is_empty()) else {
            if text.is_empty() {
                return Ok(resolver);
            }
            return Err(refused(format!("{path} starts with a blank line")));
        };
        let names: Vec<&str> = header.split('\t').collect();
        // `dict(zip(fieldnames, row))`: a repeated name takes its last column.
        let column = |name: &str| names.iter().rposition(|candidate| *candidate == name);
        let fields: Vec<(u8, usize)> = LOOKUP_FIELD_PRIORITY
            .iter()
            .enumerate()
            .filter_map(|(priority, name)| column(name).map(|place| (priority as u8, place)))
            .collect();
        let symbol_column = column("symbol");
        for line in lines {
            // csv.DictReader skips blank rows; the split after a final line break is one.
            if line.is_empty() {
                continue;
            }
            let values: Vec<&str> = line.split('\t').collect();
            if values.len() < names.len() {
                return Err(refused(format!(
                    "{path} has a row shorter than its header (csv.DictReader fills it with None)"
                )));
            }
            if values
                .iter()
                .any(|value| value.chars().count() > CSV_FIELD_LIMIT)
            {
                return Err(refused(format!(
                    "{path} has a field past the csv module's size limit"
                )));
            }
            let row =
                u32::try_from(resolver.symbols.len()).map_err(|_| invalid("too many rows"))?;
            for (priority, place) in &fields {
                for value in split_multi(&clean_text(values[*place])) {
                    let key = clean_text(&value).to_uppercase();
                    if !key.is_empty() {
                        resolver
                            .lookup
                            .entry(key)
                            .or_default()
                            .push((*priority, row));
                    }
                }
            }
            resolver
                .symbols
                .push(symbol_column.map_or(String::new(), |place| values[place].to_owned()));
        }
        Ok(resolver)
    }

    /// `resolve(query, warn_ambiguous=...)`: the one symbol of the rows that hold
    /// the query in their best field; the cleaned query itself when nothing
    /// matches or those rows name more than one symbol (with a warning, as
    /// Python logs one).
    pub(crate) fn resolve(&self, query: Option<Cell<'_>>, warn_ambiguous: bool) -> String {
        let cleaned = clean(query);
        if cleaned.is_empty() {
            return cleaned;
        }
        let Some(matches) = self.lookup.get(&cleaned.to_uppercase()) else {
            return cleaned;
        };
        let best = matches
            .iter()
            .map(|(priority, _)| *priority)
            .min()
            .expect("a lookup entry is never empty");
        let mut symbols: Vec<&str> = matches
            .iter()
            .filter(|(priority, _)| *priority == best)
            .map(|(_, row)| self.symbols[*row as usize].as_str())
            .filter(|symbol| !symbol.is_empty())
            .collect();
        symbols.sort_unstable();
        symbols.dedup();
        match symbols.as_slice() {
            [] => cleaned,
            [symbol] => (*symbol).to_owned(),
            _ => {
                if warn_ambiguous {
                    eprintln!(
                        "WARNING:resolve:Ambiguous HGNC query {cleaned} resolves to {}; keeping input symbol",
                        symbols.join(",")
                    );
                }
                cleaned
            }
        }
    }
}

impl SymbolResolver for HgncResolver {
    /// `hub.resolve_symbol(gene_symbol)` (the resolver warns about ambiguous queries).
    fn resolve_symbol(&self, gene_symbol: Cell<'_>) -> Result<String> {
        Ok(self.resolve(Some(gene_symbol), true))
    }

    /// `hub._try_resolve_symbol(gene_symbol)` (no warning).
    fn try_resolve_symbol(&self, gene_symbol: Cell<'_>) -> String {
        self.resolve(Some(gene_symbol), false)
    }
}

/// The resolver a harness call uses: option `hgnc_table`, else the repository's table.
pub(crate) fn harness_resolver(options: &Json) -> Result<HgncResolver> {
    let path = match options.get("hgnc_table") {
        Some(Json::Str(path)) => path.clone(),
        Some(_) => return Err(invalid("the option hgnc_table must be a path")),
        None => default_data_path(DEFAULT_HGNC_TABLE),
    };
    HgncResolver::load(&path)
}

// ----------------------------------------------------------------------------- the hub

/// The integrated condition cache as `_load_condition_cache` keeps it: the
/// genes (cache key, record) and each symbol's gene, under its cache spelling
/// and its resolved symbol.
#[derive(Clone, Debug, Default)]
pub(crate) struct ConditionCache {
    pub(crate) genes: Vec<(String, Json)>,
    pub(crate) by_symbol: OrderedMap<usize>,
}

impl ConditionCache {
    /// `condition_records.get(symbol, {})`: `None` stands for the empty dict.
    pub(crate) fn gene(&self, symbol: &str) -> Option<&Json> {
        self.by_symbol
            .get(symbol)
            .map(|index| &self.genes[*index].1)
    }
}

/// One ClinGen dosage record (`_load_clingen`).
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ClinGenRecord {
    pub(crate) genomic_location: String,
    pub(crate) haploinsufficiency_score: Option<i128>,
    pub(crate) haploinsufficiency_description: String,
    pub(crate) triplosensitivity_score: Option<i128>,
    pub(crate) triplosensitivity_description: String,
}

/// The files of a `GeneMechanismHub(...)`; the defaults are the repository's.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct HubPaths {
    pub(crate) condition_cache: String,
    pub(crate) clingen_dosage: String,
    pub(crate) loeuf_table: String,
    pub(crate) hgnc_table: String,
}

impl Default for HubPaths {
    fn default() -> Self {
        HubPaths {
            condition_cache: default_data_path(DEFAULT_HPO_CONDITION_MECHANISM_CACHE),
            clingen_dosage: default_data_path(DEFAULT_CLINGEN_DOSAGE),
            loeuf_table: default_data_path(DEFAULT_LOEUF_TABLE),
            hgnc_table: default_data_path(DEFAULT_HGNC_TABLE),
        }
    }
}

/// `GeneMechanismHub(..., use_hgnc_package=False)`: the resolver is read at
/// construction, every other resource on first use, as in Python.
pub(crate) struct GeneMechanismHub {
    pub(crate) paths: HubPaths,
    pub(crate) resolver: HgncResolver,
    condition_cache: Option<Rc<ConditionCache>>,
    clingen: Option<Rc<OrderedMap<ClinGenRecord>>>,
    loeuf: Option<Rc<OrderedMap<f64>>>,
}

impl SymbolResolver for GeneMechanismHub {
    fn resolve_symbol(&self, gene_symbol: Cell<'_>) -> Result<String> {
        self.resolver.resolve_symbol(gene_symbol)
    }

    fn try_resolve_symbol(&self, gene_symbol: Cell<'_>) -> String {
        self.resolver.try_resolve_symbol(gene_symbol)
    }
}

fn value_error(message: String) -> Box<dyn std::error::Error> {
    python_error("ValueError", message)
}

/// `pd.read_csv(path, sep="\t", dtype=str, low_memory=False, ...)`; a missing
/// file is the FileNotFoundError pandas raises.
fn read_text_table(path: &str, use_columns: Option<UseColumns>) -> Result<Table> {
    if let Err(error) = std::fs::metadata(path) {
        return Err(open_error(path, error));
    }
    read_table_with(
        Path::new(path),
        &ReadOptions {
            all_text: true,
            use_columns,
            ..ReadOptions::default()
        },
    )
}

/// `row.get(name)` of a `dtype=str` table after `.fillna("")`: the text, `""`
/// for a missing value, `None` for an absent column.
fn text_cell<'a>(table: &'a Table, name: &str, row: usize) -> Option<Cell<'a>> {
    table
        .column(name)
        .ok()
        .map(|column| match column.cell(row) {
            Cell::Missing => Cell::Text(""),
            cell => cell,
        })
}

impl GeneMechanismHub {
    pub(crate) fn new(paths: HubPaths) -> Result<GeneMechanismHub> {
        let resolver = HgncResolver::load(&paths.hgnc_table)?;
        Ok(GeneMechanismHub {
            paths,
            resolver,
            condition_cache: None,
            clingen: None,
            loeuf: None,
        })
    }

    /// `_try_resolve_symbol` of a cleaned symbol.
    fn try_resolve(&self, symbol: &str) -> String {
        self.resolver.resolve(Some(Cell::Text(symbol)), false)
    }

    /// `_load_condition_cache()`: the versioned integrated cache, its genes
    /// indexed by cache spelling and by resolved symbol.
    pub(crate) fn load_condition_cache(&mut self) -> Result<Rc<ConditionCache>> {
        if self.condition_cache.is_none() {
            let cache = self.read_condition_cache()?;
            self.condition_cache = Some(Rc::new(cache));
        }
        Ok(Rc::clone(
            self.condition_cache.as_ref().expect("loaded above"),
        ))
    }

    fn read_condition_cache(&self) -> Result<ConditionCache> {
        let path = &self.paths.condition_cache;
        let shown = python_path_text(path);
        if !Path::new(path).is_file() {
            return Err(python_error(
                "FileNotFoundError",
                format!(
                    "Integrated HPO condition-mechanism cache not found: {shown}. Run `bash scripts/install_utils.sh hpo_condition_mechanism_cache_install config.yaml`."
                ),
            ));
        }
        let text = open_text(path)?;
        let payload = json_text::parse(&text, false).map_err(|message| {
            if json_text::is_refusal(&message) {
                refused(format!("{path}: {message}"))
            } else {
                python_error("JSONDecodeError", message)
            }
        })?;
        let Json::Object(mut items) = payload else {
            return Err(value_error(format!("{shown} must contain a JSON object")));
        };
        let meta = items
            .iter()
            .find(|(key, _)| key == "_meta")
            .map(|(_, value)| value);
        let Some(meta @ Json::Object(_)) = meta else {
            return Err(value_error(format!("{shown} is missing object _meta")));
        };
        let schema_version = clean_json(meta.get("schema_version"))?;
        if schema_version != HPO_CONDITION_MECHANISM_SCHEMA_VERSION {
            return Err(value_error(format!(
                "Unsupported condition cache schema {}; expected '{HPO_CONDITION_MECHANISM_SCHEMA_VERSION}'",
                python_text_repr(&schema_version)?
            )));
        }
        let Some(position) = items.iter().position(|(key, _)| key == "genes") else {
            return Err(value_error(format!("{shown} is missing object genes")));
        };
        let Json::Object(genes) = items.swap_remove(position).1 else {
            return Err(value_error(format!("{shown} is missing object genes")));
        };
        let mut cache = ConditionCache::default();
        for (cache_symbol, gene) in genes {
            if !matches!(gene, Json::Object(_)) {
                return Err(value_error(format!(
                    "Condition cache gene {} must be an object",
                    python_text_repr(&cache_symbol)?
                )));
            }
            let symbol = clean_text(&cache_symbol);
            let index = cache.genes.len();
            cache.genes.push((cache_symbol, gene));
            if symbol.is_empty() {
                continue;
            }
            let resolved = self.try_resolve(&symbol);
            cache.by_symbol.insert(symbol, index);
            if !resolved.is_empty() {
                cache.by_symbol.set_default(resolved, index);
            }
        }
        Ok(cache)
    }

    /// `_load_clingen()`: the ClinGen dosage records by symbol (the first record
    /// of a symbol wins), each also under its resolved symbol.
    pub(crate) fn load_clingen(&mut self) -> Result<Rc<OrderedMap<ClinGenRecord>>> {
        if self.clingen.is_none() {
            let table = read_text_table(&self.paths.clingen_dosage, None)?;
            let mut by_symbol = OrderedMap::default();
            for row in 0..table.rows {
                let symbol = clean(text_cell(&table, "#Gene Symbol", row));
                if symbol.is_empty() {
                    continue;
                }
                let record = ClinGenRecord {
                    genomic_location: clean(text_cell(&table, "Genomic Location", row)),
                    haploinsufficiency_score: safe_int(text_cell(
                        &table,
                        "Haploinsufficiency Score",
                        row,
                    ))?,
                    haploinsufficiency_description: clean(text_cell(
                        &table,
                        "Haploinsufficiency Description",
                        row,
                    )),
                    triplosensitivity_score: safe_int(text_cell(
                        &table,
                        "Triplosensitivity Score",
                        row,
                    ))?,
                    triplosensitivity_description: clean(text_cell(
                        &table,
                        "Triplosensitivity Description",
                        row,
                    )),
                };
                let resolved = self.try_resolve(&symbol);
                by_symbol.set_default(symbol, record.clone());
                if !resolved.is_empty() {
                    by_symbol.set_default(resolved, record);
                }
            }
            self.clingen = Some(Rc::new(by_symbol));
        }
        Ok(Rc::clone(self.clingen.as_ref().expect("loaded above")))
    }

    /// `_load_loeuf()`: each gene's lowest LOEUF (`oe_lof_upper`) over its
    /// canonical transcripts, else over all of them; genes in sorted order,
    /// canonical ones first; then each value also under its resolved symbol.
    pub(crate) fn load_loeuf(&mut self) -> Result<Rc<OrderedMap<f64>>> {
        if self.loeuf.is_none() {
            let wanted = ["#gene", "canonical", "oe_lof_upper"].map(str::to_owned);
            let table = read_text_table(
                &self.paths.loeuf_table,
                Some(UseColumns::Names(wanted.to_vec())),
            )?;
            let text = |name: &str, row: usize| -> Result<String> {
                Ok(match table.column(name)?.cell(row) {
                    Cell::Text(text) => text.to_owned(),
                    Cell::Missing => String::new(),
                    _ => return Err(invalid("a dtype=str column holds text")),
                })
            };
            // groupby("#gene") (sorted keys) over all rows and over the canonical ones.
            let mut all: std::collections::BTreeMap<String, Vec<f64>> = Default::default();
            let mut canonical: std::collections::BTreeMap<String, Vec<f64>> = Default::default();
            for row in 0..table.rows {
                let gene = text("#gene", row)?;
                let value = safe_float(Some(Cell::Text(&text("oe_lof_upper", row)?)))?;
                let is_canonical = text("canonical", row)?.to_lowercase() == "true";
                all.entry(gene.clone()).or_default().push(value);
                if is_canonical {
                    canonical.entry(gene).or_default().push(value);
                }
            }
            // `min(values)` keeps the first of equal values.
            let lowest = |values: &[f64]| -> Option<f64> {
                values
                    .iter()
                    .copied()
                    .filter(|value| !value.is_nan())
                    .fold(None, |best, value| match best {
                        Some(best) if best <= value => Some(best),
                        _ => Some(value),
                    })
            };
            let mut out: OrderedMap<f64> = OrderedMap::default();
            for (gene, values) in &canonical {
                if let Some(value) = lowest(values) {
                    out.set_default(gene.clone(), value);
                }
            }
            for (gene, values) in &all {
                if out.contains_key(gene) {
                    continue;
                }
                if let Some(value) = lowest(values) {
                    out.insert(gene.clone(), value);
                }
            }
            let mut resolved_out = out.clone();
            for (gene, value) in out.entries() {
                let resolved = self.try_resolve(gene);
                if !resolved.is_empty() {
                    resolved_out.set_default(resolved, *value);
                }
            }
            self.loeuf = Some(Rc::new(resolved_out));
        }
        Ok(Rc::clone(self.loeuf.as_ref().expect("loaded above")))
    }
}

// ------------------------------------------------------------------------------- harness

/// The hub of a harness call: options `condition_cache`, `clingen_dosage`,
/// `loeuf_table`, `hgnc_table` (absent ones are the repository's files).
pub(crate) fn harness_hub(options: &Json) -> Result<GeneMechanismHub> {
    GeneMechanismHub::new(harness_paths(options)?)
}

/// The hub files a harness call names.
pub(crate) fn harness_paths(options: &Json) -> Result<HubPaths> {
    let mut paths = HubPaths::default();
    for (name, slot) in [
        ("condition_cache", &mut paths.condition_cache),
        ("clingen_dosage", &mut paths.clingen_dosage),
        ("loeuf_table", &mut paths.loeuf_table),
        ("hgnc_table", &mut paths.hgnc_table),
    ] {
        match options.get(name) {
            Some(Json::Str(path)) => *slot = path.clone(),
            Some(_) => return Err(invalid(format!("the option {name} must be a path"))),
            None => {}
        }
    }
    Ok(paths)
}

fn texts_column(values: Vec<String>) -> Column {
    Column::Text(values.into_iter().map(Some).collect())
}

fn show_score(score: Option<i128>) -> String {
    score.map_or("None".to_owned(), |score| format!("i:{score}"))
}

/// `priva acmg-step gene_mechanism_resources <function>`:
/// `resolve` (column `query`, optional bool column `warn`; one result per row),
/// `_load_condition_cache` (columns `symbol`, `cache_symbol`), `_load_clingen`
/// (`symbol` and the record's fields), `_load_loeuf` (`symbol`, `loeuf`).
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    match function {
        "resolve" => {
            let resolver = harness_resolver(&call.options)?;
            let table = &call.table;
            let queries = table.column("query")?;
            let warn = table.column("warn").ok();
            let results = (0..table.rows)
                .map(|row| {
                    let warn_ambiguous =
                        !matches!(warn.map(|column| column.cell(row)), Some(Cell::Bool(false)));
                    show_outcome(Ok(show_text(
                        &resolver.resolve(Some(queries.cell(row)), warn_ambiguous),
                    )))
                })
                .collect::<Result<Vec<String>>>()?;
            Ok(super::gene_mechanism_common::result_table(results))
        }
        "_load_condition_cache" => {
            let mut hub = harness_hub(&call.options)?;
            let cache = hub.load_condition_cache()?;
            let (symbols, keys): (Vec<String>, Vec<String>) = cache
                .by_symbol
                .entries()
                .iter()
                .map(|(symbol, index)| (symbol.clone(), cache.genes[*index].0.clone()))
                .unzip();
            Ok(Table {
                names: vec!["symbol".to_owned(), "cache_symbol".to_owned()],
                rows: symbols.len(),
                columns: vec![texts_column(symbols), texts_column(keys)],
            })
        }
        "_load_clingen" => {
            let mut hub = harness_hub(&call.options)?;
            let clingen = hub.load_clingen()?;
            let records = clingen.entries();
            let pick = |field: fn(&ClinGenRecord) -> String| {
                texts_column(records.iter().map(|(_, record)| field(record)).collect())
            };
            Ok(Table {
                names: [
                    "symbol",
                    "genomic_location",
                    "haploinsufficiency_score",
                    "haploinsufficiency_description",
                    "triplosensitivity_score",
                    "triplosensitivity_description",
                ]
                .map(str::to_owned)
                .to_vec(),
                rows: records.len(),
                columns: vec![
                    texts_column(records.iter().map(|(symbol, _)| symbol.clone()).collect()),
                    pick(|record| record.genomic_location.clone()),
                    pick(|record| show_score(record.haploinsufficiency_score)),
                    pick(|record| record.haploinsufficiency_description.clone()),
                    pick(|record| show_score(record.triplosensitivity_score)),
                    pick(|record| record.triplosensitivity_description.clone()),
                ],
            })
        }
        "_load_loeuf" => {
            let mut hub = harness_hub(&call.options)?;
            let loeuf = hub.load_loeuf()?;
            let values = loeuf.entries();
            Ok(Table {
                names: vec!["symbol".to_owned(), "loeuf".to_owned()],
                rows: values.len(),
                columns: vec![
                    texts_column(values.iter().map(|(symbol, _)| symbol.clone()).collect()),
                    Column::Float(values.iter().map(|(_, value)| *value).collect()),
                ],
            })
        }
        _ => Err(not_converted("gene_mechanism_resources", function)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TABLE: &str = "hgnc_id\tsymbol\talias_symbol\tprev_symbol\tentrez_id\n\
        HGNC:1\tAAA\tX1|shared\tOLD\t11\n\
        HGNC:2\tBBB\tshared| nan |X2\tOLD\t22\n\
        \n\
        HGNC:3\tCCC\t\tAAA\t33\n";

    #[test]
    fn resolver_follows_priority_and_ambiguity() {
        let resolver = HgncResolver::parse(TABLE, "t.tsv").unwrap();
        let resolve = |query: &str| resolver.resolve(Some(Cell::Text(query)), false);
        assert_eq!(resolve(" aaa "), "AAA");
        assert_eq!(resolve("x2"), "BBB");
        assert_eq!(resolve("hgnc:3"), "CCC");
        // Both rows hold "shared" as an alias: ambiguous, the input stays.
        assert_eq!(resolve("Shared"), "Shared");
        assert_eq!(resolve("old"), "old");
        assert_eq!(resolve("nan"), "");
        assert_eq!(resolve("ZZZ"), "ZZZ");
        assert_eq!(resolver.resolve(Some(Cell::Int(22)), false), "BBB");
        assert_eq!(resolver.resolve(Some(Cell::Float(22.0)), false), "22.0");
        assert!(HgncResolver::parse("a\tb\n1\n", "t.tsv").is_err());
        assert!(HgncResolver::parse("a\tb\n\"1\"\t2\n", "t.tsv").is_err());
        assert!(HgncResolver::parse("", "t.tsv").is_ok());
    }

    #[test]
    fn ordered_map_keeps_first_places() {
        let mut map = OrderedMap::default();
        map.insert("a".to_owned(), 1);
        map.set_default("b".to_owned(), 2);
        map.set_default("a".to_owned(), 3);
        map.insert("b".to_owned(), 4);
        assert_eq!(map.entries(), &[("a".to_owned(), 1), ("b".to_owned(), 4)]);
    }
}
