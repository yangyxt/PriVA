//! JSON/binary64 transport for candidate and aggregation comparisons.
#[allow(dead_code)]
#[path = "../../src/das_rescue.rs"]
mod das_rescue;
#[path = "../../src/das_rescue_plan.rs"]
mod das_rescue_plan;
#[allow(dead_code)]
#[path = "../../src/exon_af_stats.rs"]
mod exon_af_stats;
#[allow(dead_code)]
#[path = "../../src/fdr.rs"]
mod fdr;
#[allow(dead_code)]
#[path = "../../src/json_text.rs"]
mod json_text;
#[allow(dead_code)]
#[path = "../../src/sha256.rs"]
mod sha256;
use das_rescue_plan::{Inputs, Mode};
use json_text::Json;
use std::collections::{BTreeSet, HashMap};
use std::io::{self, BufRead};
type Error = Box<dyn std::error::Error>;

fn get<'a>(value: &'a Json, key: &str) -> Result<&'a Json, Error> {
    value
        .get(key)
        .ok_or_else(|| format!("missing {key}").into())
}
fn array(value: &Json) -> Result<&[Json], Error> {
    match value {
        Json::Array(values) => Ok(values),
        _ => Err("expected array".into()),
    }
}
fn text(value: &Json) -> Result<&str, Error> {
    value.as_str().ok_or_else(|| "expected text".into())
}
fn integer(value: &Json) -> Result<i128, Error> {
    match value {
        Json::Int(value) => Ok(*value),
        _ => Err("expected integer".into()),
    }
}
fn float(value: &Json) -> Result<f64, Error> {
    Ok(f64::from_bits(u64::from_str_radix(text(value)?, 16)?))
}
fn floats(value: &Json) -> Result<Vec<f64>, Error> {
    array(value)?.iter().map(float).collect()
}
fn digest(values: &[f64]) -> String {
    let mut hash = sha256::Sha256::new();
    for value in values {
        hash.update(&value.to_bits().to_be_bytes());
    }
    hash.finish_hex()
}
fn object(values: impl IntoIterator<Item = (&'static str, Json)>) -> Json {
    Json::Object(
        values
            .into_iter()
            .map(|(name, value)| (name.into(), value))
            .collect(),
    )
}
fn strings<'a>(values: impl IntoIterator<Item = &'a str>) -> Json {
    Json::Array(values.into_iter().map(|s| Json::Str(s.into())).collect())
}
fn error(category: &str, calls: &[String]) -> Json {
    object([
        ("error", Json::Str(category.into())),
        ("calls", strings(calls.iter().map(String::as_str))),
    ])
}

fn evaluate(settings: &Json, rows: &[Json]) -> Result<Json, Error> {
    let mode = match text(get(settings, "mode")?)? {
        "go" => Mode::GoMf,
        "non" => Mode::NonGoMf,
        _ => return Err("unknown mode".into()),
    };
    let mut scores = HashMap::new();
    let mut entries = HashMap::new();
    let mut pvalues = HashMap::new();
    let mut paths = Vec::new();
    let mut errors = BTreeSet::new();
    for row in rows {
        let domain = text(get(row, "domain")?)?.to_owned();
        paths.push(domain.clone());
        let raw = get(row, "scores")?;
        if !matches!(raw, Json::Null) {
            scores.insert(domain.clone(), floats(raw)?);
        }
        if matches!(row.get("lookup_error"), Some(Json::Bool(true))) {
            errors.insert(domain.clone());
        }
        let raw = get(row, "entry")?;
        if !matches!(raw, Json::Null) {
            let fields = array(raw)?;
            if fields.len() != 2 {
                return Err("entry needs type and GO flag".into());
            }
            entries.insert(
                domain.clone(),
                das_rescue::Entry {
                    kind: text(&fields[0])?.into(),
                    has_go_molecular_function: matches!(fields[1], Json::Bool(true)),
                },
            );
        }
        pvalues.insert(domain, float(get(row, "pvalue")?)?);
    }
    let global = floats(get(settings, "global")?)?;
    let mut by_type = HashMap::new();
    for row in array(get(settings, "composites")?)? {
        let row = array(row)?;
        if row.len() != 2 {
            return Err("composite needs type and values".into());
        }
        by_type.insert(text(&row[0])?.to_owned(), floats(&row[1])?);
    }
    let mut calls = Vec::new();
    let plan = das_rescue_plan::plan(
        mode,
        paths.iter().map(String::as_str),
        Inputs {
            scores: &scores,
            global: &global,
            by_type: &by_type,
            minimum_median: float(get(settings, "minimum")?)?,
            seed: integer(get(settings, "seed")?)? as u32,
        },
        |domain| {
            calls.push(
                domain
                    .split_once(':')
                    .map_or(domain, |(_, rest)| rest)
                    .to_owned(),
            );
            if errors.contains(domain) {
                Err("KeyError")
            } else {
                Ok(entries.get(domain).cloned())
            }
        },
    );
    let plan = match plan {
        Ok(value) => value,
        Err(category) => return Ok(error(category, &calls)),
    };
    let threads = i64::try_from(integer(get(settings, "threads")?)?)?;
    let workers = match plan.workers(threads) {
        Ok(value) => value,
        Err(_) => return Ok(error("ValueError", &calls)),
    };
    let probabilities: Vec<_> = plan.tasks.iter().map(|task| pvalues[task.domain]).collect();
    let decisions = plan.finish(&probabilities, float(get(settings, "alpha")?)?)?;
    let mut rescued = BTreeSet::new();
    if decisions.include_direct {
        rescued.extend(plan.direct.iter().copied());
    }
    let mut summaries = Vec::new();
    for decision in &decisions.tested {
        let task = &plan.tasks[decision.task_index];
        if decision.rescued {
            rescued.insert(task.domain);
        }
        summaries.push(Json::Array(vec![
            Json::Str(task.domain.into()),
            Json::Str(format!(
                "{:016x}",
                probabilities[decision.task_index].to_bits()
            )),
            Json::Str(format!("{:016x}", decision.adjusted.to_bits())),
            Json::Bool(decision.rescued),
        ]));
    }
    let global_hash = digest(&global);
    let type_hashes: HashMap<_, _> = by_type
        .iter()
        .map(|(k, v)| (k.as_str(), digest(v)))
        .collect();
    let tasks = plan
        .tasks
        .iter()
        .map(|task| {
            Json::Array(vec![
                Json::Str(task.domain.into()),
                Json::Int(task.query.len() as i128),
                Json::Str(digest(task.query)),
                task.composite_type
                    .map_or(Json::Null, |s| Json::Str(s.into())),
                Json::Int(task.composite.len() as i128),
                Json::Str(
                    task.composite_type
                        .map_or(&global_hash, |kind| &type_hashes[kind])
                        .clone(),
                ),
                Json::Int(integer(get(settings, "permutations").unwrap()).unwrap()),
                Json::Int(task.seed as i128),
            ])
        })
        .collect();
    Ok(object([
        ("error", Json::Null),
        ("calls", strings(calls.iter().map(String::as_str))),
        ("eligible", Json::Int(plan.eligible as i128)),
        ("workers", Json::Int(workers as i128)),
        ("tasks", Json::Array(tasks)),
        ("direct", strings(plan.direct)),
        ("rescued", strings(rescued)),
        ("summary", Json::Array(summaries)),
    ]))
}

fn main() -> Result<(), Error> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args != ["plan"] && args != ["stream"] {
        return Err("usage: probe plan|stream".into());
    }
    let mut lines = io::stdin().lock().lines();
    if args == ["stream"] {
        let settings = json_text::parse(&lines.next().ok_or("missing settings")??, true)?;
        let rows: Vec<_> = lines
            .map(|line| Ok(json_text::parse(&line?, true)?))
            .collect::<Result<_, Error>>()?;
        println!("{}", json_text::dumps(&evaluate(&settings, &rows)?, false)?);
    } else {
        for line in lines {
            let value = json_text::parse(&line?, true)?;
            println!(
                "{}",
                json_text::dumps(&evaluate(&value, array(get(&value, "domains")?)?)?, false)?
            );
        }
    }
    Ok(())
}
