//! Whole DAS step in the order of `analyze_domain_data`: domain score arrays,
//! curated functional references, fixed composites, per-query KS with BH FDR,
//! the four post-KS filters, both rescues, then the published sets and tables.
//! Each step calls its separately qualified native procedure. Deliberate
//! differences: KS rows follow query order rather than worker completion
//! order; rescues use the approved domain-sorted-rescue-v1 order; outputs are
//! published only after the whole procedure succeeds.
use crate::alpha_domain::{das_projection, metadata};
use crate::builder_error;
use crate::cache_convert::Provenance;
use crate::cache_convert::interpro_mapping::InterproCache;
use crate::cache_convert::intolerant_domains;
use crate::cli::{Result, invalid};
use crate::das_coldspot_filter;
use crate::das_domain::DomainLookup;
use crate::das_reference::{self, Reference};
use crate::das_rescue::Entry;
use crate::das_rescue_plan::{self, Inputs, Mode};
use crate::das_sampling;
use crate::fdr;
use crate::json_text::{self, Json};
use crate::numpy_rng::Pcg64;
use crate::protein_domain_mapping::DomainNormalizer;
use crate::python_text::float_repr;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

pub(crate) const HELP: &str = "usage: priva das-domains analyze --statistics PARQUET --source-builder-sha256 SHA256 \
--interpro RKYV --normalization-json JSON --python-hash-seed UINT32 --threads INTEGER --assembly NAME \
--output-dir DIRECTORY [--metadata PARQUET] [--spots TSV] [--clinvar-vcf VCF] [--fdr-threshold FLOAT] \
[--coldspot-overlap-threshold FLOAT] [--common-af-threshold FLOAT]\n\
Runs the whole DAS step of am_pick_intolerant_domains.py. --metadata is the native\n\
domain metadata (the source's --am_vcf); the coldspot filter also needs --spots.\n\
Writes domain_tolerance_analysis, domain_filter_results, the native\n\
all_intolerant_domains set, and the curated reference and excluded sets as JSON.";

/// `analyze_domain_data` literals.
const SEED: u64 = 42;
const T_CAP: usize = 200_000;
const T_MIN: usize = 50;
const MIN_REF_FOR_TYPE: usize = 100;
const MIN_BENIGN: u64 = 3;
/// `identify_functional_domain` and `filter_whole_protein_family` defaults.
const REFERENCE_COVERAGE: f64 = 0.9;
const REFERENCE_MIN_LENGTH: i64 = 300;
const FAMILY_COVERAGE: f64 = 0.9;
const FAMILY_MIN_LENGTH: i64 = 250;
/// Rescue defaults: permutations, minimum median AM score and seed.
const PERMUTATIONS: i64 = 500;
const GO_MF: (f64, u32) = (0.564, 42);
const NON_GO_MF: (f64, u32) = (0.75, 99);
const RESCUE_PROFILE: &str = "domain-sorted-rescue-v1";

struct Options {
    statistics: PathBuf,
    source_builder: String,
    interpro: String,
    normalization: PathBuf,
    hash_seed: u32,
    threads: usize,
    assembly: String,
    output: PathBuf,
    metadata: Option<PathBuf>,
    spots: Option<PathBuf>,
    clinvar: Option<PathBuf>,
    fdr: f64,
    overlap: f64,
    common: f64,
}

fn options(args: &[String]) -> Result<Options> {
    const REQUIRED: [&str; 8] = [
        "--statistics",
        "--source-builder-sha256",
        "--interpro",
        "--normalization-json",
        "--python-hash-seed",
        "--threads",
        "--assembly",
        "--output-dir",
    ];
    const OPTIONAL: [&str; 6] = [
        "--metadata",
        "--spots",
        "--clinvar-vcf",
        "--fdr-threshold",
        "--coldspot-overlap-threshold",
        "--common-af-threshold",
    ];
    if args.len() % 2 != 0 {
        return Err(invalid(HELP));
    }
    let mut values = HashMap::new();
    for pair in args.chunks_exact(2) {
        let name = pair[0].as_str();
        if !(REQUIRED.contains(&name) || OPTIONAL.contains(&name))
            || values.insert(name, pair[1].as_str()).is_some()
        {
            return Err(invalid(HELP));
        }
    }
    if !REQUIRED.iter().all(|name| values.contains_key(name)) {
        return Err(invalid(HELP));
    }
    let source = values["--source-builder-sha256"];
    if source.len() != 64 || !source.bytes().all(|b| b.is_ascii_hexdigit()) {
        return Err(invalid("--source-builder-sha256 requires SHA256"));
    }
    let threads = values["--threads"]
        .parse::<usize>()
        .ok()
        .filter(|&threads| threads > 0)
        .ok_or_else(|| invalid("--threads requires a positive integer"))?;
    let float = |name: &str, default: f64| -> Result<f64> {
        match values.get(name) {
            None => Ok(default),
            Some(text) => text
                .parse::<f64>()
                .ok()
                .filter(|value| value.is_finite())
                .ok_or_else(|| invalid(format!("{name} requires a finite FLOAT"))),
        }
    };
    let assembly = values["--assembly"];
    if assembly.is_empty() || !assembly.bytes().all(|b| b.is_ascii_alphanumeric()) {
        return Err(invalid("--assembly requires letters and digits"));
    }
    Ok(Options {
        statistics: values["--statistics"].into(),
        source_builder: source.to_owned(),
        interpro: values["--interpro"].to_owned(),
        normalization: values["--normalization-json"].into(),
        hash_seed: values["--python-hash-seed"]
            .parse()
            .map_err(|_| invalid("--python-hash-seed requires UINT32"))?,
        threads,
        assembly: assembly.to_owned(),
        output: values["--output-dir"].into(),
        metadata: values.get("--metadata").map(PathBuf::from),
        spots: values.get("--spots").map(PathBuf::from),
        clinvar: values.get("--clinvar-vcf").map(PathBuf::from),
        fdr: float("--fdr-threshold", 0.01)?,
        overlap: float("--coldspot-overlap-threshold", 0.6)?,
        common: float("--common-af-threshold", das_projection::AF_COMMON)?,
    })
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args.first().is_some_and(|arg| arg == "--help") {
        println!("{HELP}");
        return Ok(());
    }
    let options = options(args)?;
    let counts = builder_error::finish(analyze(&options))?;
    println!("{counts}");
    Ok(())
}

fn sampling(error: das_sampling::Error) -> Box<dyn std::error::Error> {
    builder_error::category(error.category, invalid(error.message))
}

fn reference_failure(failure: das_reference::Failure) -> Box<dyn std::error::Error> {
    builder_error::category(failure.category(), invalid("domain coverage failed"))
}

/// Results in item order; workers take the next unclaimed index.
fn parallel<T: Sync, R: Send>(
    items: &[T],
    threads: usize,
    work: impl Fn(&T) -> std::result::Result<R, das_sampling::Error> + Sync,
) -> Result<Vec<R>> {
    let next = AtomicUsize::new(0);
    let done = Mutex::new(Vec::with_capacity(items.len()));
    std::thread::scope(|scope| {
        for _ in 0..threads.min(items.len()) {
            scope.spawn(|| {
                let mut local = Vec::new();
                loop {
                    let at = next.fetch_add(1, Ordering::Relaxed);
                    let Some(item) = items.get(at) else { break };
                    local.push((at, work(item)));
                }
                done.lock().unwrap().extend(local);
            });
        }
    });
    let mut done = done.into_inner().unwrap();
    done.sort_unstable_by_key(|(at, _)| *at);
    // The first failing item in input order is reported, as a sequential loop would.
    done.into_iter()
        .map(|(_, result)| result.map_err(sampling))
        .collect()
}

/// `max(50, min(value, 200000))`, the source's composite size rule.
fn composite_size(value: usize) -> i64 {
    value.clamp(T_MIN, T_CAP) as i64
}

fn analyze(options: &Options) -> Result<String> {
    // 1. Domain score arrays in `visualize_domain_distribution` order.
    let projection =
        das_projection::project(&options.statistics, &options.source_builder, options.common)?;
    let order: Vec<String> = projection
        .scores
        .iter()
        .map(|(domain, _)| domain.clone())
        .collect();
    let scores: HashMap<String, Vec<f64>> = projection.scores.iter().cloned().collect();
    let metadata: Option<HashMap<String, metadata::Entry>> = match &options.metadata {
        Some(path) => Some(metadata::read(path)?.entries.into_iter().collect()),
        None => None,
    };
    let lookup = DomainLookup::new(
        DomainNormalizer::from_parts(
            &std::fs::read_to_string(&options.normalization)?,
            HashMap::new(),
        )?,
        InterproCache::open(&options.interpro)?,
        options.hash_seed,
    );
    // Every later lookup is of a scored domain, and the first is made in this order.
    let mut entries: HashMap<&str, Option<Entry>> = HashMap::with_capacity(order.len());
    for domain in &order {
        let entry = lookup
            .first(domain)
            .map_err(invalid)?
            .map(|found| found.entry);
        entries.insert(domain, entry);
    }
    let kind = |domain: &str| entries[domain].as_ref().map(|entry| entry.kind.clone());
    let located = |domain: &str| -> das_reference::Metadata {
        let meta = &metadata.as_ref().unwrap()[domain];
        das_reference::Metadata {
            start: meta.start,
            end: meta.end,
            protein_length: Some(meta.length),
        }
    };
    let has_metadata = |domain: &str| metadata.as_ref().is_some_and(|m| m.contains_key(domain));

    // 2. Curated functional references, in domain-score order.
    let mut references = Vec::new();
    for domain in &order {
        let meta = has_metadata(domain.as_str()).then(|| located(domain.as_str()));
        let decision = das_reference::reference(
            entries[domain.as_str()].as_ref(),
            meta,
            REFERENCE_COVERAGE,
            REFERENCE_MIN_LENGTH,
        )
        .map_err(reference_failure)?;
        if decision == Reference::Functional {
            references.push(domain.as_str());
        }
    }
    let reference_set: HashSet<&str> = references.iter().copied().collect();
    let queries: Vec<&str> = order
        .iter()
        .map(String::as_str)
        .filter(|domain| !reference_set.contains(domain))
        .collect();
    if references.is_empty() || queries.is_empty() {
        return Err(invalid(
            "No references or no query domains found; DAS analysis did not complete",
        ));
    }

    // 3. One generator: the global composite, then type composites in the
    // order each type first appears among the references.
    let t = composite_size(queries.iter().map(|d| scores[*d].len()).max().unwrap());
    let mut rng = Pcg64::new(SEED);
    let arrays = |domains: &[&str]| -> Vec<Vec<f64>> {
        domains
            .iter()
            .map(|domain| scores[*domain].clone())
            .collect()
    };
    let global =
        das_sampling::fixed_composite(&arrays(&references), t, &mut rng).map_err(sampling)?;
    let mut by_type_order: Vec<(String, Vec<&str>)> = Vec::new();
    for &domain in &references {
        let Some(kind) = kind(domain) else { continue };
        match by_type_order.iter_mut().find(|(name, _)| *name == kind) {
            Some((_, members)) => members.push(domain),
            None => by_type_order.push((kind, vec![domain])),
        }
    }
    let mut by_type: HashMap<String, Vec<f64>> = HashMap::new();
    for (kind, members) in &by_type_order {
        if members.len() < MIN_REF_FOR_TYPE {
            continue;
        }
        let total = members.iter().map(|d| scores[*d].len()).sum::<usize>();
        let composite =
            das_sampling::fixed_composite(&arrays(members), composite_size(total), &mut rng)
                .map_err(sampling)?;
        by_type.insert(kind.clone(), composite);
    }

    // 4. KS of each query against its type composite, else the global one.
    let tasks: Vec<(&str, &[f64], u64)> = queries
        .iter()
        .enumerate()
        .map(|(at, domain)| {
            let composite = kind(*domain)
                .and_then(|kind| by_type.get(&kind))
                .unwrap_or(&global);
            (
                *domain,
                composite.as_slice(),
                (SEED + at as u64) & 0xFFFF_FFFF,
            )
        })
        .collect();
    let ks = parallel(&tasks, options.threads, |(domain, composite, seed)| {
        let cap = (t as usize).min(composite.len()) as i64;
        das_sampling::ks_sample(
            &scores[*domain],
            composite,
            &mut Pcg64::new(*seed),
            Some(cap),
        )
    })?;

    // 5. Benjamini-Hochberg over the raw KS p-values.
    let raw: Vec<f64> = ks.iter().map(|result| result.pvalue).collect();
    let (tolerant_flags, adjusted) = fdr::statsmodels_bh(&raw, options.fdr).map_err(invalid)?;
    let mut tolerance = Table::new(&[
        "domain",
        "ks_stat",
        "p_value",
        "fdr_corrected_p_value",
        "is_more_tolerant",
        "n_query",
        "n_comp_used",
        "composite_T",
        "median_query",
        "median_composite",
    ]);
    for (at, result) in ks.iter().enumerate() {
        tolerance.push(vec![
            Cell::Text(queries[at].to_owned()),
            Cell::Float(result.statistic),
            Cell::Float(result.pvalue),
            Cell::Float(adjusted[at]),
            Cell::Bool(tolerant_flags[at]),
            Cell::Int(result.n_query as i64),
            Cell::Int(result.n_composite as i64),
            Cell::Int(global.len() as i64),
            Cell::Float(result.median_query),
            Cell::Float(result.median_composite),
        ]);
    }
    let ks_intolerant: BTreeSet<&str> = queries
        .iter()
        .zip(&tolerant_flags)
        .filter_map(|(domain, tolerant)| (!tolerant).then_some(*domain))
        .collect();
    let tolerant: Vec<&str> = queries
        .iter()
        .zip(&tolerant_flags)
        .filter_map(|(domain, tolerant)| tolerant.then_some(*domain))
        .collect();

    // 6. Filters, each over the set the previous one left.
    let mut kept = ks_intolerant.clone();
    let mut excluded: Vec<(&str, BTreeSet<&str>)> = Vec::new();
    let mut coldspot = None;
    if let (Some(spots), Some(meta)) = (&options.spots, &metadata) {
        let coldspots = das_coldspot_filter::read_coldspots(spots)?;
        let mut rows = HashMap::new();
        for domain in &kept {
            let summary =
                das_coldspot_filter::summary(domain, meta, &coldspots, &lookup, options.overlap)?;
            rows.insert(*domain, summary);
        }
        coldspot = Some(remove(&mut kept, &mut excluded, "coldspot", rows, |s| {
            s.removed
        }));
    }
    let mut benign = None;
    if let Some(vcf) = &options.clinvar {
        let domains: HashSet<String> = kept.iter().map(|d| d.to_string()).collect();
        let counts = crate::das_benign::counts(vcf, &domains)?;
        let rows = kept
            .iter()
            .map(|d| (*d, counts.get(*d).copied().unwrap_or(0)))
            .collect();
        benign = Some(remove(&mut kept, &mut excluded, "benign", rows, |n| {
            *n >= MIN_BENIGN
        }));
    }
    let rows = kept.iter().map(|d| (*d, projection.common(d))).collect();
    let common = remove(&mut kept, &mut excluded, "common", rows, |s| s.removed);
    let mut family = None;
    if metadata.is_some() {
        let mut rows = HashMap::new();
        for domain in &kept {
            let meta = has_metadata(*domain).then(|| located(*domain));
            let summary = das_reference::family(
                entries[*domain].as_ref(),
                meta,
                FAMILY_COVERAGE,
                FAMILY_MIN_LENGTH,
            )
            .map_err(reference_failure)?;
            rows.insert(*domain, summary);
        }
        family = Some(remove(&mut kept, &mut excluded, "family", rows, |s| {
            s.removed
        }));
    }

    // 7. Both rescues of the KS-tolerant domains, then the references.
    let mut rescued = Vec::new();
    for (mode, (minimum_median, seed)) in [(Mode::GoMf, GO_MF), (Mode::NonGoMf, NON_GO_MF)] {
        let inputs = Inputs {
            scores: &scores,
            global: &global,
            by_type: &by_type,
            minimum_median,
            seed,
        };
        let plan = das_rescue_plan::plan(mode, tolerant.iter().copied(), inputs, |domain| {
            Ok::<_, Box<dyn std::error::Error>>(entries[domain].clone())
        })?;
        let workers = plan.workers(options.threads as i64).map_err(invalid)?;
        let dips = parallel(&plan.tasks, workers.max(1), |task| {
            das_sampling::dip_rescue(
                task.query,
                task.composite,
                PERMUTATIONS,
                &mut Pcg64::new(task.seed as u64),
            )
        })?;
        let probabilities: Vec<f64> = dips.iter().map(|dip| dip.pvalue).collect();
        let decisions = plan.finish(&probabilities, options.fdr).map_err(invalid)?;
        let mut found: BTreeSet<&str> = decisions
            .tested
            .iter()
            .filter(|decision| decision.rescued)
            .map(|decision| plan.tasks[decision.task_index].domain)
            .collect();
        if decisions.include_direct {
            found.extend(plan.direct.iter().copied());
        }
        rescued.push(found);
    }
    let mut intolerant: BTreeSet<&str> = kept.clone();
    for found in &rescued {
        intolerant.extend(found.iter().copied());
    }
    intolerant.extend(references.iter().copied());

    // 8. Filter results for every KS-intolerant domain, in sorted order.
    let mut columns = vec![
        "domain_path",
        "gene_symbol",
        "aa_start",
        "aa_end",
        "in_original_intolerant",
    ];
    if coldspot.is_some() {
        columns.extend([
            "coldspot_max_frac_of_coldspot",
            "coldspot_max_frac_of_domain",
            "coldspot_removed",
        ]);
    }
    if benign.is_some() {
        columns.extend(["n_benign_missense", "benign_removed"]);
    }
    columns.extend(["n_common_missense", "max_af", "common_removed"]);
    if family.is_some() {
        columns.extend([
            "interpro_type",
            "domain_coverage",
            "prot_length",
            "family_removed",
        ]);
    }
    columns.push("final_intolerant");
    let mut filters = Table::new(&columns);
    for domain in &ks_intolerant {
        let meta = metadata.as_ref().and_then(|m| m.get(*domain));
        let mut row = vec![
            Cell::Text(domain.to_string()),
            Cell::Text(meta.map_or(String::new(), |m| m.symbol.clone())),
            meta.map_or(Cell::Missing, |m| Cell::Int(m.start)),
            meta.map_or(Cell::Missing, |m| Cell::Int(m.end)),
            Cell::Bool(true),
        ];
        if let Some(rows) = &coldspot {
            row.extend(match rows.get(domain) {
                Some(s) => [
                    Cell::Float(s.of_coldspot),
                    Cell::Float(s.of_domain),
                    Cell::Bool(s.removed),
                ],
                None => [Cell::Missing, Cell::Missing, Cell::Bool(false)],
            });
        }
        if let Some(rows) = &benign {
            row.extend(match rows.get(domain) {
                Some(n) => [Cell::Int(*n as i64), Cell::Bool(*n >= MIN_BENIGN)],
                None => [Cell::Missing, Cell::Bool(false)],
            });
        }
        row.extend(match common.get(domain) {
            Some(s) => [
                Cell::Int(s.fully_common as i64),
                match s.fraction {
                    crate::das_common::Fraction::IntegerZero => Cell::Int(0),
                    crate::das_common::Fraction::Float(value) => Cell::Float(value),
                },
                Cell::Bool(s.removed),
            ],
            None => [Cell::Missing, Cell::Missing, Cell::Bool(false)],
        });
        if let Some(rows) = &family {
            row.extend(match rows.get(domain) {
                Some(s) => [
                    s.kind.clone().map_or(Cell::Missing, Cell::Text),
                    s.coverage.map_or(Cell::Missing, Cell::Float),
                    s.protein_length.map_or(Cell::Missing, Cell::Int),
                    Cell::Bool(s.removed),
                ],
                None => [
                    Cell::Text(String::new()),
                    Cell::Missing,
                    Cell::Missing,
                    Cell::Bool(false),
                ],
            });
        }
        row.push(Cell::Bool(intolerant.contains(domain)));
        filters.push(row);
    }

    // 9. Publish only after every step succeeded.
    std::fs::create_dir_all(&options.output)?;
    let assembly = &options.assembly;
    let path = |name: &str, extension: &str| {
        options
            .output
            .join(format!("{name}.{assembly}.{extension}"))
    };
    let mut tolerance_bytes = Vec::new();
    tolerance.write(&mut tolerance_bytes)?;
    let mut filter_bytes = Vec::new();
    filters.write(&mut filter_bytes)?;
    publish(&path("domain_tolerance_analysis", "tsv"), &tolerance_bytes)?;
    publish(&path("domain_filter_results", "tsv"), &filter_bytes)?;
    let set = |domains: &mut dyn Iterator<Item = &str>| {
        Json::Array(domains.map(|domain| Json::Str(domain.to_owned())).collect())
    };
    let curated: BTreeSet<&str> = references.iter().copied().collect();
    let curated_json =
        json_text::dumps(&set(&mut curated.iter().copied()), false).map_err(invalid)?;
    publish(
        &path("curated_functional_references", "json"),
        format!("{curated_json}\n").as_bytes(),
    )?;
    let excluded_json = Json::Object(
        excluded
            .iter()
            .map(|(name, domains)| (name.to_string(), set(&mut domains.iter().copied())))
            .collect(),
    );
    let excluded_json = json_text::dumps(&excluded_json, false).map_err(invalid)?;
    publish(
        &path("excluded_domains", "json"),
        format!("{excluded_json}\n").as_bytes(),
    )?;
    let provenance = [
        (
            "producer",
            format!("priva das-domains analyze {}", env!("CARGO_PKG_VERSION")),
        ),
        ("statistics_builder_sha256", options.source_builder.clone()),
        ("python_hash_seed", options.hash_seed.to_string()),
        ("fdr_threshold", float_repr(options.fdr)),
        ("rescue_profile", RESCUE_PROFILE.to_owned()),
    ]
    .into_iter()
    .map(|(key, value)| Provenance {
        key: key.to_owned(),
        value,
    })
    .collect();
    let written = intolerant_domains::write_domains(
        &path("all_intolerant_domains", "rkyv"),
        intolerant.iter().map(|domain| domain.to_string()).collect(),
        provenance,
    )?;
    Ok(format!(
        "{{\"domains\":{},\"references\":{},\"queries\":{},\"ks_intolerant\":{},\"coldspot_removed\":{},\"benign_removed\":{},\"common_removed\":{},\"family_removed\":{},\"rescued_go_mf\":{},\"rescued_non_go_mf\":{},\"final_intolerant\":{written}}}",
        order.len(),
        references.len(),
        queries.len(),
        ks_intolerant.len(),
        removed_count(&excluded, "coldspot"),
        removed_count(&excluded, "benign"),
        removed_count(&excluded, "common"),
        removed_count(&excluded, "family"),
        rescued[0].len(),
        rescued[1].len(),
    ))
}

/// Records one filter's summaries and removes its excluded domains from `kept`.
fn remove<'a, S>(
    kept: &mut BTreeSet<&'a str>,
    excluded: &mut Vec<(&'static str, BTreeSet<&'a str>)>,
    name: &'static str,
    rows: HashMap<&'a str, S>,
    removed: impl Fn(&S) -> bool,
) -> HashMap<&'a str, S> {
    let gone: BTreeSet<&str> = rows
        .iter()
        .filter_map(|(domain, summary)| removed(summary).then_some(*domain))
        .collect();
    kept.retain(|domain| !gone.contains(domain));
    excluded.push((name, gone));
    rows
}

fn removed_count(excluded: &[(&str, BTreeSet<&str>)], name: &str) -> String {
    excluded
        .iter()
        .find(|(filter, _)| *filter == name)
        .map_or("null".to_owned(), |(_, domains)| domains.len().to_string())
}

/// A staged file renamed into place.
fn publish(path: &Path, bytes: &[u8]) -> Result<()> {
    let mut staged = path.as_os_str().to_owned();
    staged.push(".partial");
    let mut file = std::fs::File::create(&staged)?;
    file.write_all(bytes)?;
    file.sync_all()?;
    std::fs::rename(&staged, path)?;
    Ok(())
}

enum Cell {
    Text(String),
    Bool(bool),
    Int(i64),
    Float(f64),
    Missing,
}

/// `DataFrame(rows).to_csv(sep="\t", index=False)`: a column holding a float or
/// a missing value beside integers is float64, so its integers print as
/// floats; missing values and NaN print empty; Booleans print as True/False.
struct Table {
    columns: Vec<String>,
    rows: Vec<Vec<Cell>>,
}

impl Table {
    fn new(columns: &[&str]) -> Self {
        Self {
            columns: columns.iter().map(|name| name.to_string()).collect(),
            rows: Vec::new(),
        }
    }

    fn push(&mut self, row: Vec<Cell>) {
        assert_eq!(row.len(), self.columns.len());
        self.rows.push(row);
    }

    fn write(&self, out: &mut impl Write) -> Result<()> {
        // A frame built from no records has no columns either.
        if self.rows.is_empty() {
            writeln!(out)?;
            return Ok(());
        }
        let floats: Vec<bool> = (0..self.columns.len())
            .map(|at| {
                let column = || self.rows.iter().map(move |row| &row[at]);
                column().all(|cell| matches!(cell, Cell::Int(_) | Cell::Float(_) | Cell::Missing))
                    && column().any(|cell| matches!(cell, Cell::Float(_) | Cell::Missing))
            })
            .collect();
        let header: Vec<String> = self.columns.iter().map(|name| field(name)).collect();
        writeln!(out, "{}", header.join("\t"))?;
        for row in &self.rows {
            let cells: Vec<String> = row
                .iter()
                .zip(&floats)
                .map(|(cell, &float)| match cell {
                    Cell::Text(text) => field(text),
                    Cell::Bool(value) => String::from(if *value { "True" } else { "False" }),
                    Cell::Int(value) if float => float_repr(*value as f64),
                    Cell::Int(value) => value.to_string(),
                    Cell::Float(value) if value.is_nan() => String::new(),
                    Cell::Float(value) => float_repr(*value),
                    Cell::Missing => String::new(),
                })
                .collect();
            writeln!(out, "{}", cells.join("\t"))?;
        }
        Ok(())
    }
}

/// csv.QUOTE_MINIMAL with a tab separator.
fn field(text: &str) -> String {
    if text.contains(['\t', '"', '\n', '\r']) {
        format!("\"{}\"", text.replace('"', "\"\""))
    } else {
        text.to_owned()
    }
}
