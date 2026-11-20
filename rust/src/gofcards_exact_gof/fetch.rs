//! The workbook download and the per-allele summary requests. `requests`
//! becomes one `curl` process per request, with the same headers, redirects,
//! status errors and retries; transport failures carry curl's message.
use super::*;
use crate::gene_nonlof_mechanism::{SheetCell, get, read_first_sheet, truthy};
use std::collections::{HashMap, HashSet};
use std::process::{Command, Stdio};

const PUBLIC_COLUMN_MAP: [(&str, &str); 15] = [
    ("genesymbol", "Gene_Symbol"),
    ("transcript", "Transcript"),
    ("chr", "Chr"),
    ("hg19start", "Start"),
    ("hg19end", "End"),
    ("ref", "Ref"),
    ("alt", "Alt"),
    ("function", "Function"),
    ("pathways proteins involved", "Pathways_proteins_involved"),
    ("disorder involved", "Phenotype"),
    ("pmid", "PMID"),
    ("animal model", "Animal_model"),
    ("cell model", "Cell_model"),
    ("pscore", "Pscore"),
    ("order numbe", "source_order"),
];

const SUMMARY_FIELDS: [&str; 2] = ["hg38_start", "Accession"];
const API_HEADERS: [(&str, &str); 6] = [
    (
        "User-Agent",
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) \
Chrome/125 Safari/537.36",
    ),
    ("Accept", "application/json, text/plain, */*"),
    ("Referer", "https://www.genemed.tech/gofcards/"),
    ("Origin", "https://www.genemed.tech"),
    ("Admin-Work-Version", "x"),
    ("Content-Type", "application/json; charset=UTF-8"),
];

struct Reply {
    status: u32,
    reason: String,
    url: String,
    headers: Vec<(String, String)>,
    body: Vec<u8>,
}

impl Reply {
    /// `raise_for_status`.
    fn check(&self) -> std::result::Result<(), Failure> {
        let kind = match self.status {
            400..=499 => "Client",
            500..=599 => "Server",
            _ => return Ok(()),
        };
        Err(Failure {
            category: "HTTPError",
            message: format!(
                "{} {kind} Error: {} for url: {}",
                self.status, self.reason, self.url
            ),
        })
    }

    /// `headers.get(name, "")`, repeated headers joined as requests does.
    fn header(&self, name: &str) -> String {
        self.headers
            .iter()
            .filter(|(key, _)| key.eq_ignore_ascii_case(name))
            .map(|(_, value)| value.as_str())
            .collect::<Vec<_>>()
            .join(", ")
    }
}

struct Failure {
    category: &'static str,
    message: String,
}

impl Failure {
    fn raise(self) -> Box<dyn std::error::Error> {
        error(self.category, self.message)
    }
}

/// The status line and headers of the last response curl dumped.
fn parse_headers(text: &str) -> (String, Vec<(String, String)>) {
    let mut reason = String::new();
    let mut headers = Vec::new();
    for line in text.lines() {
        let line = line.trim_end_matches('\r');
        if line.starts_with("HTTP/") {
            reason = line.splitn(3, ' ').nth(2).unwrap_or("").to_owned();
            headers.clear();
        } else if let Some((name, value)) = line.split_once(':') {
            headers.push((name.to_owned(), value.trim().to_owned()));
        }
    }
    (reason, headers)
}

/// `session.get(url, headers=headers, timeout=timeout)`.
fn request(
    url: &str,
    headers: &[(&str, &str)],
    timeout: i128,
    scratch: &Path,
) -> Result<std::result::Result<Reply, Failure>> {
    let stem = format!(".gofcards_request.{}", std::process::id());
    let header_file = scratch.join(format!("{stem}.headers"));
    let body_file = scratch.join(format!("{stem}.body"));
    let mut command = Command::new("curl");
    let seconds = timeout.to_string();
    command.args([
        "--silent",
        "--show-error",
        "--location",
        "--compressed",
        "--http1.1",
        "--connect-timeout",
        &seconds,
        "--speed-limit",
        "1",
        "--speed-time",
        &seconds,
        "--write-out",
        "%{http_code} %{url_effective}",
    ]);
    command.arg("--dump-header").arg(&header_file);
    command.arg("--output").arg(&body_file);
    for (name, value) in headers {
        command.args(["--header", &format!("{name}: {value}")]);
    }
    command.arg(url);
    let output = command
        .stdin(Stdio::null())
        .output()
        .map_err(|e| os_error(e, Path::new("curl"), None))?;
    let dumped = fs::read(&header_file).unwrap_or_default();
    let body = fs::read(&body_file).unwrap_or_default();
    let _ = fs::remove_file(&header_file);
    let _ = fs::remove_file(&body_file);
    if !output.status.success() {
        let category = if output.status.code() == Some(28) {
            "Timeout"
        } else {
            "ConnectionError"
        };
        return Ok(Err(Failure {
            category,
            message: crate::python_text::python_strip(&String::from_utf8_lossy(&output.stderr))
                .to_owned(),
        }));
    }
    let written = String::from_utf8_lossy(&output.stdout);
    let (status, final_url) = written.split_once(' ').unwrap_or((&written, url));
    let (reason, headers) = parse_headers(&String::from_utf8_lossy(&dumped));
    Ok(Ok(Reply {
        status: status.parse().unwrap_or(0),
        reason,
        url: final_url.to_owned(),
        headers,
        body,
    }))
}

/// `urllib.parse.quote_plus`.
fn quote_plus(text: &str) -> String {
    let mut out = String::new();
    for byte in text.bytes() {
        match byte {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'_' | b'.' | b'-' | b'~' => {
                out.push(char::from(byte));
            }
            b' ' => out.push('+'),
            other => out.push_str(&format!("%{other:02X}")),
        }
    }
    out
}

/// One allele as the DataFrame columns the builder reads describe it.
pub(super) struct Source {
    pub(super) allele_key: String,
    pub(super) locus_key: String,
    pub(super) variant_type: &'static str,
    pub(super) chrom: String,
    pub(super) start: String,
    pub(super) end: String,
    pub(super) reference: String,
    pub(super) alternate: String,
    pub(super) symbol: String,
    pub(super) evidence: Json,
}

fn cell_text(cell: Option<&SheetCell>) -> Option<String> {
    match cell? {
        SheetCell::Nan => None,
        SheetCell::NaT => Some("NaT".into()),
        SheetCell::Text(text) | SheetCell::Date(text) => Some(text.clone()),
        SheetCell::Int(number) => Some(number.to_string()),
        SheetCell::Float(number) if number.is_nan() => None,
        SheetCell::Float(number) => Some(crate::python_text::float_repr(*number)),
        SheetCell::Bool(flag) => Some((if *flag { "True" } else { "False" }).into()),
    }
}

fn clean_cell(cell: Option<&SheetCell>) -> String {
    cell_text(cell).map_or_else(String::new, |text| clean_text(&text))
}

fn int_cell(cell: Option<&SheetCell>) -> String {
    match cell {
        Some(SheetCell::Float(number)) if number.is_finite() && number.fract() == 0.0 => {
            whole_float(*number)
        }
        other => cell_text(other).map_or_else(String::new, |text| strip_point_zero(&text)),
    }
}

/// `fetch_sources`: the workbook, its rows, and the provenance it starts.
pub(super) fn fetch_sources(
    work: &Path,
    url: &str,
    timeout: i128,
) -> Result<(Vec<Source>, Vec<(String, Json)>)> {
    fs::create_dir_all(work).map_err(|e| os_error(e, work, None))?;
    let public_path = child(work, "gofcards_public.xlsx");
    log(&format!(
        "Downloading public workbook -> {}",
        public_path.display()
    ));
    let reply =
        request(url, &[("User-Agent", "Mozilla/5.0")], timeout, work)?.map_err(Failure::raise)?;
    reply.check().map_err(Failure::raise)?;
    fs::write(&public_path, &reply.body).map_err(|e| os_error(e, &public_path, None))?;

    let (names, rows) = read_first_sheet(&public_path)?;
    let renamed: Vec<String> = names
        .iter()
        .map(|name| {
            let key = crate::python_text::python_strip(name).to_lowercase();
            PUBLIC_COLUMN_MAP
                .iter()
                .find(|(from, _)| *from == key)
                .map_or_else(|| name.clone(), |(_, to)| (*to).to_owned())
        })
        .collect();
    let column = |name: &str| -> Result<Option<usize>> {
        let mut found = renamed.iter().enumerate().filter(|(_, have)| *have == name);
        let first = found.next().map(|(index, _)| index);
        if found.next().is_some() {
            return Err(error(
                "ValueError",
                format!("The workbook has more than one {name} column"),
            ));
        }
        Ok(first)
    };
    let [
        gene,
        transcript,
        chrom,
        start,
        end,
        reference,
        alternate,
        function,
        pathway,
        phenotype,
        pmid,
        animal,
        cell_model,
        pscore,
        order,
    ] = [
        column("Gene_Symbol")?,
        column("Transcript")?,
        column("Chr")?,
        column("Start")?,
        column("End")?,
        column("Ref")?,
        column("Alt")?,
        column("Function")?,
        column("Pathways_proteins_involved")?,
        column("Phenotype")?,
        column("PMID")?,
        column("Animal_model")?,
        column("Cell_model")?,
        column("Pscore")?,
        column("source_order")?,
    ];
    let mut sources = Vec::with_capacity(rows.len());
    for row in &rows {
        let at = |index: Option<usize>| index.and_then(|index| row.get(index));
        let reference = normalize_allele(&clean_cell(at(reference)));
        let alternate = normalize_allele(&clean_cell(at(alternate)));
        let variant_type = if reference.chars().count() == 1 && alternate.chars().count() == 1 {
            "SNV"
        } else {
            "Indel"
        };
        let chrom = normalize_chrom(&clean_cell(at(chrom)));
        let start = int_cell(at(start));
        let end = int_cell(at(end));
        let locus_key = [
            chrom.as_str(),
            start.as_str(),
            end.as_str(),
            reference.as_str(),
            alternate.as_str(),
        ]
        .join("|");
        let source_transcript = clean_cell(at(transcript));
        sources.push(Source {
            allele_key: format!("{variant_type}|{locus_key}"),
            locus_key,
            variant_type,
            chrom,
            start,
            end,
            reference,
            alternate,
            symbol: clean_cell(at(gene)),
            evidence: Json::Object(
                [
                    ("source_order", clean_cell(at(order))),
                    ("pmid", clean_cell(at(pmid))),
                    ("pscore", clean_cell(at(pscore))),
                    ("disease", clean_cell(at(phenotype))),
                    ("function", clean_cell(at(function))),
                    ("pathway", clean_cell(at(pathway))),
                    ("animal_model", clean_cell(at(animal))),
                    ("cell_model", clean_cell(at(cell_model))),
                    ("source_transcript", source_transcript),
                ]
                .into_iter()
                .map(|(key, value)| (key.to_owned(), Json::Str(value)))
                .collect(),
            ),
        });
    }
    let unique: HashSet<&str> = sources.iter().map(|s| s.locus_key.as_str()).collect();
    let provenance = vec![
        ("public_excel_url".to_owned(), Json::Str(url.to_owned())),
        (
            "public_excel_last_modified".to_owned(),
            Json::Str(reply.header("Last-Modified")),
        ),
        (
            "public_workbook_rows".to_owned(),
            Json::Int(sources.len() as i128),
        ),
        (
            "public_unique_alleles".to_owned(),
            Json::Int(unique.len() as i128),
        ),
    ];
    log(&format!(
        "public workbook: {} rows, {} unique alleles",
        sources.len(),
        unique.len()
    ));
    Ok((sources, provenance))
}

/// One attempt at an allele: the reply's `data` block.
fn summary_attempt(
    url: &str,
    timeout: i128,
    work: &Path,
) -> Result<std::result::Result<Json, String>> {
    let reply = match request(url, &API_HEADERS, timeout, work)? {
        Ok(reply) => reply,
        Err(failure) => return Ok(Err(failure.message)),
    };
    if let Err(failure) = reply.check() {
        return Ok(Err(failure.message));
    }
    let parsed = match json_text::parse(&String::from_utf8_lossy(&reply.body), false) {
        Ok(parsed) => parsed,
        Err(message) => return Ok(Err(message)),
    };
    if !truthy(&parsed) {
        return Ok(Ok(Json::Object(vec![])));
    }
    if !matches!(parsed, Json::Object(_)) {
        return Ok(Err(format!(
            "'{}' object has no attribute 'get'",
            crate::gene_nonlof_mechanism::type_name(&parsed)
        )));
    }
    Ok(Ok(parsed
        .get("data")
        .filter(|data| truthy(data))
        .cloned()
        .unwrap_or(Json::Object(vec![]))))
}

fn read_summary_cache(path: &Path) -> Result<HashMap<String, Json>> {
    let mut cache = HashMap::new();
    for line in read_lines(path)? {
        let line = crate::python_text::python_strip(&line);
        if line.is_empty() {
            continue;
        }
        let row = json_text::parse(line, false).map_err(|e| error("JSONDecodeError", e))?;
        let Json::Object(_) = &row else {
            return Err(error(
                "TypeError",
                format!(
                    "{} indices must be integers",
                    crate::gene_nonlof_mechanism::type_name(&row)
                ),
            ));
        };
        let key = row
            .get("allele_key")
            .ok_or_else(|| error("KeyError", "'allele_key'"))?;
        let data = row
            .get("data")
            .filter(|data| truthy(data))
            .cloned()
            .unwrap_or(Json::Object(vec![]));
        cache.insert(crate::gene_nonlof_mechanism::py_str(key), data);
    }
    log(&format!(
        "summary cache already holds {} alleles",
        cache.len()
    ));
    Ok(cache)
}

/// `fetch_summary_annotations`.
pub(super) fn fetch_summary_annotations(
    records: &[Source],
    cache_jsonl: &Path,
    endpoint: &str,
    timeout: i128,
    work: &Path,
) -> Result<HashMap<String, Json>> {
    let mut cache = if cache_jsonl.exists() {
        read_summary_cache(cache_jsonl)?
    } else {
        HashMap::new()
    };
    if let Some(parent) = cache_jsonl.parent() {
        fs::create_dir_all(parent).map_err(|e| os_error(e, parent, None))?;
    }
    let mut seen = HashSet::new();
    let pending: Vec<&Source> = records
        .iter()
        .filter(|record| seen.insert(record.locus_key.as_str()))
        .filter(|record| !cache.contains_key(&record.allele_key))
        .collect();
    if !pending.is_empty() {
        log(&format!(
            "querying the summary endpoint for {} uncached alleles",
            pending.len()
        ));
    }
    let mut handle = OpenOptions::new()
        .create(true)
        .append(true)
        .open(cache_jsonl)
        .map_err(|e| os_error(e, cache_jsonl, None))?;
    for (index, record) in pending.iter().enumerate() {
        let index = index + 1;
        let level = if record.variant_type == "Indel" {
            "Indel"
        } else {
            "SNV"
        };
        let query = [
            ("projectCode", "GoFCards"),
            ("variantLevelType", level),
            ("chr", record.chrom.as_str()),
            ("start", record.start.as_str()),
            ("end", record.end.as_str()),
            ("ref", record.reference.as_str()),
            ("alt", record.alternate.as_str()),
        ]
        .iter()
        .map(|(key, value)| format!("{key}={}", quote_plus(value)))
        .collect::<Vec<_>>()
        .join("&");
        let url = format!("{endpoint}?{query}");
        let mut data = Json::Object(vec![]);
        for attempt in 1..=4u64 {
            match summary_attempt(&url, timeout, work)? {
                Ok(found) => {
                    data = found;
                    break;
                }
                Err(message) if attempt == 4 => {
                    return Err(error(
                        "RuntimeError",
                        format!(
                            "summary endpoint failed for {} after 4 attempts: {message}",
                            record.allele_key
                        ),
                    ));
                }
                Err(_) => std::thread::sleep(std::time::Duration::from_secs(5 * attempt)),
            }
        }
        let mut kept = Vec::new();
        for field in SUMMARY_FIELDS {
            if let Some(value) = get(&data, field)?.filter(|value| truthy(value)) {
                kept.push((field.to_owned(), value.clone()));
            }
        }
        let kept = Json::Object(kept);
        let line = json_text::dumps_with_options(
            &Json::Object(vec![
                ("allele_key".into(), Json::Str(record.allele_key.clone())),
                ("data".into(), kept.clone()),
            ]),
            DumpsOptions {
                indent: None,
                sort_keys: true,
                ensure_ascii: false,
                compact: false,
                allow_nan: true,
            },
        )
        .map_err(|e| error("ValueError", e))?;
        cache.insert(record.allele_key.clone(), kept);
        handle
            .write_all(format!("{line}\n").as_bytes())
            .map_err(|e| os_error(e, cache_jsonl, None))?;
        if index == 1 || index % 200 == 0 {
            log(&format!("  summary {index}/{}", pending.len()));
        }
    }
    let mut with_hg38 = 0;
    let mut with_clinvar = 0;
    for value in cache.values() {
        with_hg38 += usize::from(get(value, "hg38_start")?.is_some_and(truthy));
    }
    for value in cache.values() {
        with_clinvar += usize::from(get(value, "Accession")?.is_some_and(truthy));
    }
    log(&format!(
        "summary endpoint: {} alleles; {with_hg38} carry a GRCh38 position; \
{with_clinvar} carry a ClinVar accession",
        cache.len()
    ));
    Ok(cache)
}
