//! `plan` and `release`.
use super::digest::Installer;
use super::{
    Options, Result, absolute, delete_key, fail, identity, object, output_of, pretty, read_json,
    step, steps, string, strings, today_utc, truthy, write_file, yaml_keys, yaml_value,
};
use crate::json_text::Json;
use crate::sha256::Sha256;
use std::collections::HashSet;
use std::path::Path;
use std::process::{Command, Stdio};

/// Bash's `[[ TEXT == PATTERN ]]` for the patterns the registry's settings
/// use: `*`, `?` and `[...]`.
pub(crate) fn glob_match(pattern: &str, text: &str) -> bool {
    fn class(pattern: &[char], at: usize, c: char) -> Option<(bool, usize)> {
        let mut i = at + 1;
        let negate = matches!(pattern.get(i), Some('!') | Some('^'));
        if negate {
            i += 1;
        }
        let mut matched = false;
        let mut first = true;
        while i < pattern.len() && (first || pattern[i] != ']') {
            first = false;
            if pattern.get(i + 1) == Some(&'-') && pattern.get(i + 2).is_some_and(|e| *e != ']') {
                if pattern[i] <= c && c <= pattern[i + 2] {
                    matched = true;
                }
                i += 3;
            } else {
                if pattern[i] == c {
                    matched = true;
                }
                i += 1;
            }
        }
        (i < pattern.len()).then_some((matched != negate, i + 1))
    }
    fn go(p: &[char], pi: usize, t: &[char], ti: usize) -> bool {
        if pi == p.len() {
            return ti == t.len();
        }
        match p[pi] {
            '*' => (ti..=t.len()).any(|k| go(p, pi + 1, t, k)),
            '?' => ti < t.len() && go(p, pi + 1, t, ti + 1),
            '[' if ti < t.len() => match class(p, pi, t[ti]) {
                Some((true, next)) => go(p, next, t, ti + 1),
                Some((false, _)) => false,
                None => p[pi] == t[ti] && go(p, pi + 1, t, ti + 1),
            },
            '\\' if pi + 1 < p.len() => {
                ti < t.len() && p[pi + 1] == t[ti] && go(p, pi + 2, t, ti + 1)
            }
            c => ti < t.len() && c == t[ti] && go(p, pi + 1, t, ti + 1),
        }
    }
    let p: Vec<char> = pattern.chars().collect();
    let t: Vec<char> = text.chars().collect();
    go(&p, 0, &t, 0)
}

fn listed(list: &str, id: &str) -> bool {
    format!(",{list},").contains(&format!(",{id},"))
}

fn days_in_month(year: i64, month: i64) -> i64 {
    match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        _ if (year % 4 == 0 && year % 100 != 0) || year % 400 == 0 => 29,
        _ => 28,
    }
}

/// YYYY-MM-DD and N: N calendar months later, the day clamped to the month.
pub(crate) fn add_months(date: &str, months: i64) -> Option<String> {
    let year: i64 = date.get(0..4)?.parse().ok()?;
    let month: i64 = date.get(5..7)?.parse().ok()?;
    let day: i64 = date.get(8..10)?.parse().ok()?;
    let total = year * 12 + month - 1 + months;
    let (year, month) = (total.div_euclid(12), total.rem_euclid(12) + 1);
    let day = day.min(days_in_month(year, month));
    Some(format!("{year:04}-{month:02}-{day:02}"))
}

/// A pin taken fewer than MONTHS calendar months ago.
fn pin_is_fresh(pin: &Json, months: i64) -> bool {
    let recorded = pin.get("recorded").and_then(Json::as_str).unwrap_or("");
    !recorded.is_empty() && add_months(recorded, months).is_some_and(|limit| today_utc() < limit)
}

/// A pins.json entry: `{"release": ..., "recorded": "YYYY-MM-DD"}`, or an old
/// entry that is the release text alone.
pub(crate) fn normalize_pin(pin: Option<&Json>) -> Json {
    match pin {
        None | Some(Json::Null) => Json::Null,
        Some(Json::Str(text)) => object(vec![("release".into(), string(text.clone()))]),
        Some(other) => other.clone(),
    }
}

fn int(value: Option<&Json>) -> i64 {
    match value {
        Some(Json::Int(n)) => *n as i64,
        Some(Json::Float(f)) => *f as i64,
        _ => 0,
    }
}

/// `{KEY}` replaced by the value of KEY, first match first, every occurrence.
fn expand_keys(mut text: String, value: impl Fn(&str) -> Result<String>) -> Result<String> {
    loop {
        let Some((start, end)) = first_key(&text) else {
            return Ok(text);
        };
        let key = text[start + 1..end - 1].to_owned();
        let replacement = value(&key)?;
        text = text.replace(&format!("{{{key}}}"), &replacement);
    }
}

fn first_key(text: &str) -> Option<(usize, usize)> {
    let bytes = text.as_bytes();
    let mut at = 0;
    while let Some(open) = bytes[at..].iter().position(|b| *b == b'{').map(|o| at + o) {
        let mut end = open + 1;
        while end < bytes.len() && (bytes[end].is_ascii_alphanumeric() || bytes[end] == b'_') {
            end += 1;
        }
        if end > open + 1 && bytes.get(end) == Some(&b'}') {
            return Some((open, end + 1));
        }
        at = open + 1;
    }
    None
}

pub(crate) fn plan(options: &Options) -> Result<()> {
    let registry_path = options.path("registry")?;
    let config = options.path("config")?;
    let out = options.path("out")?;
    if !config.is_file() {
        return Err(fail("plan needs --registry, --config and --out"));
    }
    let registry = read_json(&registry_path)?;
    let root = options.source_root()?;
    std::fs::create_dir_all(&out)
        .map_err(|e| fail(format!("cannot create {}: {e}", out.display())))?;
    let out = absolute(&out);
    let installer = Installer::load(&root)?;

    // Keys a step writes, and folder keys, are removed so each step resolves
    // its own destinations under its own cache_root; files you supply stay.
    let base = out.join("base.yaml");
    std::fs::copy(&config, &base)
        .map_err(|e| fail(format!("cannot copy {}: {e}", config.display())))?;
    let keep: HashSet<String> = steps(&registry)?
        .iter()
        .flat_map(|s| strings(s, "user_inputs"))
        .collect();
    let mut written: Vec<String> = strings(&registry, "folder_keys");
    for entry in steps(&registry)? {
        written.extend(strings(entry, "outputs"));
        written.extend(strings(entry, "optional_outputs"));
    }
    written.sort();
    written.dedup();
    for key in written.iter().filter(|k| !keep.contains(*k)) {
        delete_key(&base, key)?;
    }

    let publish_digest = publish_digest(&root)?;
    let all_keys = yaml_keys(&config)?;
    let adopt_config = options.get("adopt-config").unwrap_or("");
    let adopt = options.get("adopt").unwrap_or("");
    let run = options.get("run").unwrap_or("");
    let refresh = options.get("refresh").unwrap_or("");
    let pins = match options.get("pins") {
        Some(path) if !path.is_empty() && Path::new(path).is_file() => {
            Some(read_json(Path::new(path))?)
        }
        _ => None,
    };
    let release_files: Vec<Json> = options
        .all("releases")
        .iter()
        .map(|file| read_json(Path::new(file)))
        .collect::<Result<_>>()?;

    let mut planned: Vec<(String, Json)> = Vec::new();
    for entry in steps(&registry)? {
        let id = super::text(entry, "id");
        let digest = installer
            .digest(entry)
            .map_err(|e| fail(format!("digest of {id} failed: {e}")))?;
        let mut settings: Vec<(String, Json)> = Vec::new();
        for pattern in strings(entry, "settings") {
            for key in all_keys.iter().filter(|k| glob_match(&pattern, k)) {
                let value = string(yaml_value(&config, key)?);
                match settings.iter_mut().find(|(name, _)| name == key) {
                    Some(slot) => slot.1 = value,
                    None => settings.push((key.clone(), value)),
                }
            }
        }
        let mut files = Vec::new();
        let mut paths = Vec::new();
        for key in strings(entry, "user_inputs") {
            let value = yaml_value(&config, &key)?;
            if !value.is_empty() {
                files.push(Json::Array(vec![
                    string(key.clone()),
                    string(identity(&value)),
                ]));
                paths.push((key, string(value)));
            }
        }

        // The step's release: checked now, or for a pinned step its pin. A
        // pin older than max_age_months calendar months gives way to the check.
        let mut release = Json::Null;
        let has_release =
            matches!(entry.get("release"), Some(Json::Array(items)) if !items.is_empty());
        let pinned = truthy(entry.get("pinned"));
        if has_release {
            let pin = normalize_pin(pins.as_ref().and_then(|p| p.get(&id)));
            let months = int(entry.get("max_age_months"));
            let fresh = months == 0 || (!matches!(pin, Json::Null) && pin_is_fresh(&pin, months));
            if pinned && !listed(refresh, &id) && fresh {
                release = match pin.get("release") {
                    None | Some(Json::Null) | Some(Json::Bool(false)) => string("initial"),
                    Some(value) => value.clone(),
                };
            } else {
                for file in &release_files {
                    if file.get("step").and_then(Json::as_str) == Some(id.as_str()) {
                        release = file.get("release").cloned().unwrap_or(Json::Null);
                    }
                }
                if matches!(release, Json::Null) {
                    return Err(fail(format!("no release check result for step {id}")));
                }
            }
        }

        let mut adopted = Json::Null;
        if !adopt_config.is_empty() && ((adopt == "all" && !listed(run, &id)) || listed(adopt, &id))
        {
            let old = Path::new(adopt_config);
            let mut keys = Vec::new();
            let mut identities = Vec::new();
            let mut outputs = strings(entry, "outputs");
            outputs.extend(strings(entry, "optional_outputs"));
            for key in outputs {
                let mut value = yaml_value(old, &key)?;
                // An older installation may not name this output (registry
                // adopt_paths gives where its installer put the file).
                if value.is_empty() {
                    let template = entry
                        .get("adopt_paths")
                        .and_then(|paths| paths.get(&key))
                        .and_then(Json::as_str)
                        .unwrap_or("")
                        .to_owned();
                    value = expand_keys(template, |name| yaml_value(old, name))?;
                    if value.is_empty() || !Path::new(&value).exists() {
                        value.clear();
                    }
                }
                if value.is_empty() {
                    continue;
                }
                identities.push(string(identity(&value)));
                keys.push((key, string(value)));
            }
            adopted = object(vec![
                ("keys".into(), object(keys)),
                ("identity".into(), Json::Array(identities)),
            ]);
        }
        planned.push((
            id,
            object(vec![
                ("digest".into(), string(digest)),
                ("settings".into(), object(settings)),
                ("user_files".into(), Json::Array(files)),
                // The same files by path, for a step that runs as its own
                // process: the workflow stages them into its task. The
                // lineage counts user_files only.
                ("user_paths".into(), object(paths)),
                ("adopt".into(), adopted),
                ("release".into(), release),
                ("pinned".into(), Json::Bool(pinned)),
            ]),
        ));
    }
    let plan = object(vec![
        ("schema".into(), Json::Int(1)),
        ("base_config".into(), string(base.to_string_lossy())),
        ("publish_digest".into(), string(publish_digest)),
        ("steps".into(), object(planned)),
    ]);
    write_file(&out.join("plan.json"), &pretty(&plan)?)
}

/// What the publish task runs besides the receipts it reads: the Rust files
/// behind `priva cache-install` (prepare and publish), by content. Without
/// the Rust sources (an installed program), the SHA-256 of the program.
fn publish_digest(root: &Path) -> Result<String> {
    let table = crate::nextflow_task::rust_reach::commands(root);
    if let Some(entry) = table.get("cache-install") {
        let files = crate::nextflow_task::rust_reach::command_files(root, entry)?;
        let text: String = files
            .iter()
            .map(|(file, sum)| format!("rust {file} {sum}\n"))
            .collect();
        return Ok(crate::sha256::sha256_hex(text.as_bytes()));
    }
    let program = std::env::var("PRIVA_BINARY")
        .ok()
        .filter(|p| !p.is_empty())
        .ok_or_else(|| fail("PRIVA_BINARY must name the priva program"))?;
    crate::sha256::sha256_file(Path::new(&program)).map_err(|e| fail(format!("{program}: {e}")))
}

/// A release URL with `{grch}`, `{loftee_branch}` and `{KEY}` filled in.
fn expand_url(template: &str, config: &Path) -> Result<String> {
    let assembly = yaml_value(config, "assembly")?;
    let (grch, branch) = match assembly.as_str() {
        "hg19" | "GRCh37" => ("GRCh37", "master"),
        _ => ("GRCh38", "grch38"),
    };
    let text = template
        .replace("{grch}", grch)
        .replace("{loftee_branch}", branch);
    expand_keys(text, |key| {
        let value = yaml_value(config, key)?;
        if value.is_empty() {
            return Err(fail(format!("{key} is not set in {}", config.display())));
        }
        Ok(value)
    })
}

/// The ETag and Last-Modified header lines of the final response after
/// redirects, lower-cased: `etag: ...|last-modified: ...|`.
pub(crate) fn header_release(headers: &[u8]) -> String {
    let text: Vec<u8> = headers.iter().copied().filter(|b| *b != b'\r').collect();
    let (mut etag, mut modified) = (String::new(), String::new());
    for line in String::from_utf8_lossy(&text).split('\n') {
        let lower = line.to_ascii_lowercase();
        if lower.starts_with("http/") {
            etag.clear();
            modified.clear();
        }
        if lower.starts_with("etag:") {
            etag = lower.clone();
        }
        if lower.starts_with("last-modified:") {
            modified = lower;
        }
    }
    format!("{etag}|{modified}|")
}

/// The SHA-256 of the lines sorted bytewise (`LC_ALL=C sort | sha256sum`).
pub(crate) fn sorted_lines_sha256(content: &[u8]) -> String {
    let mut lines: Vec<&[u8]> = content.split(|b| *b == b'\n').collect();
    if content.ends_with(b"\n") || content.is_empty() {
        lines.pop();
    }
    lines.sort();
    let mut hash = Sha256::new();
    for line in lines {
        hash.update(line);
        hash.update(b"\n");
    }
    hash.finish_hex()
}

pub(crate) fn release(options: &Options) -> Result<()> {
    let registry = read_json(&options.path("registry")?)?;
    let config = options.path("config")?;
    let id = options.one("step")?;
    let out = options.path("out")?;
    let entry = step(&registry, id)?;
    let mut parts = Vec::new();
    if let Some(Json::Array(specs)) = entry.get("release") {
        for spec in specs {
            let kind = super::text(spec, "kind");
            let url = expand_url(&super::text(spec, "url"), &config)?;
            let answer = if options.stub {
                "stub".to_owned()
            } else if kind == "git" {
                let reference = spec.get("ref").and_then(Json::as_str).unwrap_or("HEAD");
                let reference = expand_url(reference, &config)?;
                let listing = output_of("git", ["ls-remote", &url, &reference]).unwrap_or_default();
                let first = String::from_utf8_lossy(&listing)
                    .lines()
                    .next()
                    .and_then(|line| line.split('\t').next())
                    .unwrap_or("")
                    .to_owned();
                if first.is_empty() {
                    return Err(fail(format!("{id}: no commit for {reference} at {url}")));
                }
                first
            } else {
                // A source that sends neither ETag nor Last-Modified (an API
                // answer) is named by the SHA-256 of its sorted lines, up to
                // 500 MB: the ClinGen evidence repository returns the same
                // rows in varying order. The size is left out: Orphadata sends
                // Content-Length with some answers and not with others, which
                // gave F4 another release, and so another folder, from one
                // launch to the next (job 1546377).
                let headers = Command::new("curl")
                    .args(["-sIL", "--max-time", "120", &url])
                    .stderr(Stdio::null())
                    .output()
                    .map(|o| o.stdout)
                    .unwrap_or_default();
                let answer = header_release(&headers);
                if answer.starts_with("||") {
                    let body = output_of(
                        "curl",
                        [
                            "-sfL",
                            "--max-time",
                            "600",
                            "--max-filesize",
                            "500000000",
                            &url,
                        ],
                    )
                    .ok_or_else(|| fail(format!("{id}: cannot read {url}")))?;
                    format!("sorted-lines-sha256:{}", sorted_lines_sha256(&body))
                } else {
                    answer
                }
            };
            parts.push(object(vec![
                ("url".into(), string(url)),
                ("release".into(), string(answer)),
            ]));
        }
    }
    let parts = Json::Array(parts);
    let record = object(vec![
        ("step".into(), string(id)),
        ("release".into(), string(super::compact(&parts)?)),
        ("sources".into(), parts),
    ]);
    write_file(&out, &pretty(&record)?)
}
