//! Installation boundary between legacy cache builders and native runtime readers.
//!
//! Builders still produce/read their original pickles. `prepare` gives them a
//! private config; `publish` converts completed files and changes the user's
//! existing runtime keys only after every selected cache succeeds. Per-archive
//! receipts retain the exact builder inputs for the next installation. A native
//! builder may instead write an archive directly; `publish` checks and adopts it
//! under the same managed naming, with a receipt whose origin is `direct`.
//! Nothing here changes a scientific producer, its order, or its cache refresh
//! policy.
use crate::cache_convert as caches;
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::sha256::{sha256_file, sha256_hex};
use crate::yaml_config::{self, ConfigDocument, Parsed, Scalar};
use std::collections::BTreeMap;
use std::fs::{self, OpenOptions};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

const KEYS: [&str; 10] = [
    "clinvar_aa_stat",
    "clinvar_splice_stat",
    "alphamissense_tranx_domain_map",
    "clinvar_patho_exon_af_stat",
    "clinvar_gene_stat",
    "clinvar_patho_af_stat",
    "pm1_regions_pkl",
    "all_intolerant_domains",
    "interpro_mapping_pickle",
    "clingen_map",
];
const HELP: &str = "usage: priva cache-install {prepare|publish} --config CONFIG --build-config PRIVATE_CONFIG --source-root REPO\n\
       priva cache-install resolve --key KEY --build-config PRIVATE_CONFIG --source-root REPO\n\
prepare creates a private builder config; publish converts/reuses completed caches, adopts checked archives written\n\
directly by native builders, and atomically updates CONFIG.\n\
resolve prints the checked native archive for KEY; a builder pickle is converted once and publish reuses that conversion.\n\
Keep the private config and its adjacent receipt for repeat installations. No annotation is rerun by this command.";
static NEXT: AtomicU64 = AtomicU64::new(0);

fn text(path: &Path) -> Result<&str> {
    path.to_str()
        .ok_or_else(|| invalid("cache-install requires UTF-8 paths"))
}

fn nonempty(path: &Path) -> Result<()> {
    let info = fs::metadata(path)?;
    if !info.is_file() || info.len() == 0 {
        return Err(invalid(format!(
            "{}: required nonempty file is missing",
            path.display()
        )));
    }
    Ok(())
}

fn document(contents: &str) -> Result<ConfigDocument> {
    match yaml_config::parse_config(contents).map_err(|e| invalid(e.to_string()))? {
        Parsed::Mapping(doc) if !doc.has_duplicate_keys() => Ok(doc),
        _ => Err(invalid(
            "cache-install needs one YAML mapping without duplicate keys",
        )),
    }
}

fn setting<'a>(doc: &'a ConfigDocument, key: &str) -> Result<Option<&'a str>> {
    match doc.scalar(key) {
        None if doc.get(key).is_none() => Ok(None),
        Some(Scalar::Null) => Ok(None),
        Some(Scalar::Str(value)) if value.is_empty() => Ok(None),
        Some(Scalar::Str(value)) => Ok(Some(value)),
        _ => Err(invalid(format!("{key} must be a path string or null"))),
    }
}

fn required<'a>(doc: &'a ConfigDocument, key: &str) -> Result<&'a str> {
    setting(doc, key)?.ok_or_else(|| invalid(format!("configure {key} first")))
}

fn profile(doc: &ConfigDocument) -> Result<String> {
    let assembly = match required(doc, "assembly")? {
        "hg19" | "GRCh37" => "GRCh37",
        "hg38" | "GRCh38" => "GRCh38",
        _ => {
            return Err(invalid(
                "unsupported assembly; expected hg19/GRCh37 or hg38/GRCh38",
            ));
        }
    };
    let release = match doc.scalar("vep_installed_version") {
        Some(Scalar::Str(value) | Scalar::Int(value))
            if !value.is_empty() && value.bytes().all(|b| b.is_ascii_digit()) =>
        {
            value
        }
        _ => return Err(invalid("vep_installed_version must be a numeric release")),
    };
    Ok(format!("{assembly}/VEP{release}"))
}

fn resolve(value: &str, repo: &Path) -> PathBuf {
    let path = Path::new(value);
    if path.is_absolute() {
        path.to_owned()
    } else {
        repo.join(path)
    }
}

fn rewrite(contents: &str, updates: &BTreeMap<String, String>) -> Result<String> {
    let doc = document(contents)?;
    let mut spans = Vec::new();
    for (key, value) in updates {
        let entry = doc
            .entry(key)
            .ok_or_else(|| invalid(format!("missing original key {key}")))?;
        spans.push((
            entry.value_start,
            entry.value_end,
            yaml_config::quoted_for_yaml(value),
        ));
    }
    spans.sort_by_key(|x| std::cmp::Reverse(x.0));
    let mut result = contents.to_owned();
    for (start, end, value) in spans {
        result.replace_range(start..end, &value);
    }
    Ok(result)
}

fn sidecar(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
}

fn write_new(path: &Path, data: &[u8]) -> Result<()> {
    let mut file = OpenOptions::new().write(true).create_new(true).open(path)?;
    if let Err(error) = file.write_all(data).and_then(|_| file.sync_all()) {
        let _ = fs::remove_file(path);
        return Err(error.into());
    }
    Ok(())
}

/// Replace by rename, never by truncating an archive/config a running reader uses.
fn atomic_write(path: &Path, data: &[u8]) -> Result<()> {
    atomic_write_checked(path, data, || Ok(()))
}

fn parent(path: &Path) -> &Path {
    path.parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."))
}

fn temporary(path: &Path, suffix: &str) -> PathBuf {
    sidecar(
        path,
        &format!(
            ".{}.{}.{}.{suffix}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap_or_default()
                .as_nanos(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ),
    )
}

fn atomic_write_checked(
    path: &Path,
    data: &[u8],
    verify: impl FnOnce() -> Result<()>,
) -> Result<()> {
    let temporary = temporary(path, "tmp");
    write_new(&temporary, data)?;
    let result = (|| -> Result<()> {
        if path.exists() {
            fs::set_permissions(&temporary, fs::metadata(path)?.permissions())?;
        }
        verify()?;
        fs::rename(&temporary, path)?;
        fs::File::open(parent(path))?.sync_all()?;
        Ok(())
    })();
    if result.is_err() {
        let _ = fs::remove_file(&temporary);
    }
    result
}

fn json(path: &Path) -> Result<Json> {
    json_text::parse(&fs::read_to_string(path)?, true).map_err(invalid)
}

fn field<'a>(value: &'a Json, key: &str) -> Result<&'a str> {
    value
        .get(key)
        .and_then(Json::as_str)
        .ok_or_else(|| invalid(format!("receipt lacks string {key}")))
}

fn receipt(fields: &[(&str, &str)]) -> Result<Vec<u8>> {
    let value = Json::Object(
        fields
            .iter()
            .map(|(k, v)| ((*k).to_owned(), Json::Str((*v).to_owned())))
            .collect(),
    );
    Ok((json_text::dumps_indented(&value).map_err(invalid)? + "\n").into_bytes())
}

fn bundle(key: &str, repo: &Path, selected_profile: &str) -> Option<(PathBuf, PathBuf)> {
    let name = match key {
        "interpro_mapping_pickle" => "data/InterPro/Interpro_entry_mapping".to_owned(),
        "clingen_map" => format!(
            "data/clingen/clingen_map.{}",
            if selected_profile.starts_with("GRCh37/") {
                "hg19"
            } else {
                "hg38"
            }
        ),
        // The full-value ClinVar gene archive exceeds the distribution limit.
        // Keep its small pickle in Git and convert under cache_root at install.
        _ => return None,
    };
    Some((
        repo.join(format!("{name}.pkl.gz")),
        repo.join(format!("{name}.v1.rkyv")),
    ))
}

fn same_file(left: &Path, right: &Path) -> bool {
    matches!((fs::canonicalize(left), fs::canonicalize(right)), (Ok(a), Ok(b)) if a == b)
}

fn prepare(config: &Path, private: &Path, repo: &Path) -> Result<()> {
    let contents = fs::read_to_string(config)?;
    let doc = document(&contents)?;
    let selected_profile = profile(&doc)?;
    let original = fs::canonicalize(config)?;
    let original_hash = sha256_hex(contents.as_bytes());
    let context_path = sidecar(private, ".cache-install.json");
    let mut replace_completed = false;
    let mut previous = None;
    if private.exists() || context_path.exists() {
        nonempty(private)?;
        let context = json(&context_path)?;
        if field(&context, "config")? != text(&original)?
            || field(&context, "source_root")? != text(repo)?
        {
            return Err(invalid(
                "private config belongs to another source root/config; retained unchanged",
            ));
        }
        replace_completed = context.get("published_sha256").and_then(Json::as_str)
            == Some(&original_hash)
            || context.get("state").and_then(Json::as_str) == Some("complete");
        if !replace_completed && field(&context, "config_sha256")? == original_hash {
            if profile(&document(&fs::read_to_string(private)?)?)? != selected_profile {
                return Err(invalid(
                    "pending builder profile changed; private config retained unchanged",
                ));
            }
            println!(
                "Resuming the existing private builder config; its completed updates are retained."
            );
            return Ok(());
        }
        // The planned runtime digest is written before config publication, so
        // a crash immediately after the rename is recoverable. Once a prior
        // publication is complete, explicit subsequent user edits are inputs
        // for a fresh builder view; pending edits are never guessed away.
        if !replace_completed {
            return Err(invalid(
                "original config changed during a pending installation; preserve/reconcile its private config first",
            ));
        }
        previous = Some(context);
    }
    let mut updates = BTreeMap::new();
    for key in KEYS {
        let Some(value) = setting(&doc, key)? else {
            continue;
        };
        let native = resolve(value, repo);
        if native.extension().and_then(|s| s.to_str()) != Some("rkyv") {
            continue;
        }
        nonempty(&native)?;
        let source = if let Some((source, archive)) =
            bundle(key, repo, &selected_profile).filter(|(_, archive)| same_file(archive, &native))
        {
            // Repository bundles need existence/nonzero only. Actual mapped
            // reads still enforce the format's memory-safety contract.
            nonempty(&archive)?;
            source
        } else {
            let record = json(&sidecar(&native, ".source.json"))?;
            if field(&record, "key")? != key
                || field(&record, "profile")? != selected_profile
                || !same_file(Path::new(field(&record, "archive")?), &native)
            {
                return Err(invalid(format!(
                    "{key}: native cache receipt does not match this path/profile"
                )));
            }
            match record.get("origin").and_then(Json::as_str) {
                // A native builder's own archive has no builder input to
                // restore; the private view keeps selecting the archive.
                Some("direct") => continue,
                None => PathBuf::from(field(&record, "source")?),
                Some(_) => {
                    return Err(invalid(format!(
                        "{key}: native cache receipt has an unknown origin"
                    )));
                }
            }
        };
        // Both managed and bundled archives can have explicitly selected
        // source aliases. Preserve each configuration's choice without changing
        // shared archive provenance or adding content checks to bundle reuse.
        let source = match previous.as_ref().filter(|old| {
            old.get(&format!("native_{key}"))
                .and_then(Json::as_str)
                .is_some_and(|path| same_file(Path::new(path), &native))
        }) {
            Some(old) => PathBuf::from(field(old, &format!("source_{key}"))?),
            None => source,
        };
        nonempty(&source)?;
        updates.insert(key.to_owned(), text(&source)?.to_owned());
    }
    let context = receipt(&[
        ("config", text(&original)?),
        ("config_sha256", &original_hash),
        ("source_root", text(repo)?),
        ("profile", &selected_profile),
    ])?;
    let builder = rewrite(&contents, &updates)?;
    if replace_completed {
        atomic_write(private, builder.as_bytes())?;
        atomic_write(&context_path, &context)?;
    } else {
        write_new(private, builder.as_bytes())?;
        if let Err(error) = write_new(&context_path, &context) {
            let _ = fs::remove_file(private);
            return Err(error);
        }
    }
    println!("Builder config prepared; runtime config remains unchanged.");
    Ok(())
}

fn check_context(config: &Path, private: &Path, repo: &Path) -> Result<Json> {
    let context = json(&sidecar(private, ".cache-install.json"))?;
    if field(&context, "config")? != text(&fs::canonicalize(config)?)?
        || field(&context, "source_root")? != text(repo)?
        || field(&context, "config_sha256")? != sha256_file(config)?
    {
        return Err(invalid(
            "original config/source root changed after prepare; refusing to overwrite it",
        ));
    }
    Ok(context)
}

fn directory(root: &Path, requested: &Path) -> Result<PathBuf> {
    let relative = requested
        .strip_prefix(root)
        .map_err(|_| invalid("native cache destination must be inside cache_root"))?;
    let mut path = root.to_owned();
    for part in relative.components() {
        let std::path::Component::Normal(name) = part else {
            return Err(invalid("cache destination contains a non-normal component"));
        };
        path.push(name);
        match fs::symlink_metadata(&path) {
            Ok(info) if info.is_dir() && !info.file_type().is_symlink() => (),
            Ok(_) => {
                return Err(invalid(format!(
                    "{}: managed cache directory is not a real directory",
                    path.display()
                )));
            }
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => fs::create_dir(&path)?,
            Err(e) => return Err(e.into()),
        }
    }
    Ok(path)
}

fn target(
    key: &str,
    source_hash: &str,
    doc: &ConfigDocument,
    repo: &Path,
    schema_version: u32,
) -> Result<PathBuf> {
    let root_input = Path::new(required(doc, "cache_root")?);
    if !root_input.is_absolute() {
        return Err(invalid("cache_root must be an existing absolute directory"));
    }
    let root = fs::canonicalize(root_input)?;
    if !root.is_dir() || root.parent().is_none() {
        return Err(invalid("cache_root must be a dedicated existing directory"));
    }
    let family = if key.starts_with("clinvar_") {
        setting(doc, "clinvar_vcf_dir")?
            .map(|v| resolve(v, repo))
            .unwrap_or_else(|| root_input.join("ClinVar"))
    } else if [
        "alphamissense_tranx_domain_map",
        "pm1_regions_pkl",
        "all_intolerant_domains",
    ]
    .contains(&key)
    {
        setting(doc, "vep_plugins_cachedir")?
            .map(|v| resolve(v, repo))
            .unwrap_or_else(|| root_input.join("VEP_plugin_cache"))
            .join("AlphaMissense")
    } else {
        root_input.join(if key == "clingen_map" {
            "ClinGen"
        } else {
            "InterPro"
        })
    };
    // A user may select a symlink for cache_root itself, but not redirect a
    // managed child. Translate the one approved root alias before checking.
    let family = if let Ok(relative) = family.strip_prefix(root_input) {
        root.join(relative)
    } else {
        family
    };
    let family = directory(&root, &family)?;
    // Source-specific names preserve the previous runtime profile during a
    // refresh. The full digest is verified in the adjacent receipt.
    let selected_profile = profile(doc)?.replace('/', ".");
    Ok(family.join(format!(
        "{key}.{selected_profile}.{}.v{schema_version}.rkyv",
        &source_hash[..16]
    )))
}

fn validate(key: &str, path: &Path, schema_version: u32) -> Result<()> {
    let path = text(path)?;
    match key {
        "clinvar_aa_stat" => {
            caches::ClinvarAaCache::open(path)?;
        }
        "clinvar_splice_stat" => {
            let cache = caches::clinvar_splice_stat::typed::Cache::open(path)?;
            if cache.version() != schema_version {
                return Err(invalid(
                    "splice archive schema version differs from its managed generation",
                ));
            }
        }
        "alphamissense_tranx_domain_map" => {
            caches::transcript_domains::TranscriptDomainCache::open(path)?;
        }
        "clinvar_patho_exon_af_stat" => {
            caches::clinvar_exon_af::ClinvarExonAfCache::open(path)?;
        }
        "clinvar_gene_stat" => {
            let cache = caches::clinvar_gene_stat::ClinvarGeneCache::open(path)?;
            if cache.schema_version() != schema_version {
                return Err(invalid(
                    "gene-statistics archive schema version differs from its managed generation",
                ));
            }
        }
        "clinvar_patho_af_stat" => {
            let cache = caches::clinvar_gene_af::ClinvarGeneAfCache::open(path)?;
            if cache.schema_version() != schema_version {
                return Err(invalid(
                    "gene-AF archive schema version differs from its managed generation",
                ));
            }
        }
        "pm1_regions_pkl" => {
            caches::pm1_regions::Pm1RegionCache::open(path)?;
        }
        "all_intolerant_domains" => {
            caches::intolerant_domains::IntolerantDomainCache::open(path)?;
        }
        "interpro_mapping_pickle" => {
            caches::interpro_mapping::InterproCache::open(path)?;
        }
        "clingen_map" => {
            let cache = caches::clingen_map::ClingenCache::open(path)?;
            if cache.schema_version() != schema_version {
                return Err(invalid(
                    "ClinGen archive schema version differs from its managed generation",
                ));
            }
        }
        _ => return Err(invalid("unknown native cache")),
    }
    Ok(())
}

/// Check an archive a native builder wrote and return its own schema version.
/// Schema ids follow `KEYS` order (`CACHE_FORMATS.md`, section 4), so an
/// archive of another cache is refused before any reader maps it.
fn check_direct(key: &str, archive: &Path) -> Result<u32> {
    if !archive.is_file() {
        return Err(invalid(format!(
            "{key}: native archive {} named by the builder config does not exist",
            archive.display()
        )));
    }
    nonempty(archive)?;
    let Some(position) = KEYS.iter().position(|name| *name == key) else {
        return Err(invalid("unknown native cache"));
    };
    let expected = position as u32 + 1;
    let mut header = [0_u8; crate::rust_cache::HEADER_LEN];
    let mut file = fs::File::open(archive)?;
    if file.read_exact(&mut header).is_err() {
        return Err(invalid(format!(
            "{key}: {} is too short to be a native archive",
            archive.display()
        )));
    }
    if &header[0..8] != crate::rust_cache::MAGIC {
        return Err(invalid(format!(
            "{key}: {} is not a PRIVARKV archive",
            archive.display()
        )));
    }
    let word = |at: usize| u32::from_le_bytes(header[at..at + 4].try_into().unwrap());
    if word(12) != expected {
        return Err(invalid(format!(
            "{key}: {} holds cache schema {}, not schema {expected}",
            archive.display(),
            word(12)
        )));
    }
    validate(key, archive, word(16))?;
    Ok(word(16))
}

fn install_one(
    key: &str,
    source: &Path,
    doc: &ConfigDocument,
    repo: &Path,
    selected_profile: &str,
) -> Result<PathBuf> {
    nonempty(source)?;
    if let Some((pickle, native)) =
        bundle(key, repo, selected_profile).filter(|(pickle, _)| same_file(pickle, source))
    {
        nonempty(&pickle)?;
        nonempty(&native)?;
        // Progress goes to stderr: `resolve` prints only the archive path.
        eprintln!("{key}: reuse repository-bundled native cache");
        return Ok(native);
    }
    // This fingerprint identifies the bytes converted; it does not check a
    // repository-bundled pickle against an expected publisher checksum.
    let source_hash = sha256_file(source)?;
    let schema_version = if key == "clingen_map" {
        caches::clingen_map::pickle_schema_version(text(source)?)?
    } else if key == "clinvar_patho_af_stat" {
        caches::clinvar_gene_af::pickle_schema_version(text(source)?)?
    } else if key == "clinvar_splice_stat" {
        caches::clinvar_splice_stat::pickle_schema_version(text(source)?)?
    } else if key == "clinvar_gene_stat" {
        caches::clinvar_gene_stat::pickle_schema_version(text(source)?)?
    } else {
        1
    };
    let schema_text = schema_version.to_string();
    let output = target(key, &source_hash, doc, repo, schema_version)?;
    let record_path = sidecar(&output, ".source.json");
    if output.exists() || record_path.exists() {
        // Do not overwrite a supplied/unregistered file or silently adopt a
        // partial archive. Every reusable managed archive has a closed receipt.
        let record = json(&record_path)?;
        if field(&record, "key")? != key
            || field(&record, "format")? != "PRIVARKV-v1"
            || field(&record, "source_sha256")? != source_hash
            || field(&record, "profile")? != selected_profile
            || field(&record, "archive")? != text(&output)?
            || (schema_version != 1
                && record.get("schema_version").and_then(Json::as_str)
                    != Some(schema_text.as_str()))
            || record
                .get("schema_version")
                .is_some_and(|value| value.as_str() != Some(schema_text.as_str()))
        {
            return Err(invalid(format!(
                "{}: conflicting cache receipt; original retained",
                output.display()
            )));
        }
        finish_archive(key, &output, &record, schema_version)?;
        nonempty(&output)?;
        if sha256_file(&output)? != field(&record, "archive_sha256")? {
            return Err(invalid(format!(
                "{}: archive checksum differs from completed conversion; original retained",
                output.display()
            )));
        }
        validate(key, &output, schema_version)?;
        eprintln!("{key}: reuse completed native conversion");
        return Ok(output);
    }
    let pending = temporary(&output, "pending");
    caches::dispatch(&[
        key.to_owned(),
        "--pickle".to_owned(),
        text(source)?.to_owned(),
        "--output".to_owned(),
        text(&pending)?.to_owned(),
    ])?;
    if sha256_file(source)? != source_hash {
        return Err(invalid(format!(
            "{} changed while converting; runtime config was not published",
            source.display()
        )));
    }
    let archive_hash = sha256_file(&pending)?;
    // Preserve the selected logical path, including DataLad's annex link.
    // Restoring an annex object's physical path would bypass its managed path.
    let selected_source = std::path::absolute(source)?;
    let source_bytes = fs::metadata(source)?.len().to_string();
    let archive_bytes = fs::metadata(&pending)?.len().to_string();
    let mut fields = vec![
        ("format", "PRIVARKV-v1"),
        ("key", key),
        ("profile", selected_profile),
        ("source", text(&selected_source)?),
        ("source_sha256", &source_hash),
        ("source_bytes", &source_bytes),
        ("archive", text(&output)?),
        ("archive_sha256", &archive_hash),
        ("archive_bytes", &archive_bytes),
        ("pending_archive", text(&pending)?),
    ];
    // Keep existing version-1 receipt bytes/fields unchanged. Container v1 is
    // independent of the optional ClinGen schema v2 payload.
    if schema_version != 1 {
        fields.push(("schema_version", &schema_text));
    }
    publish_generation(key, &output, &pending, &fields, schema_version)?;
    Ok(output)
}

fn publish_generation(
    key: &str,
    output: &Path,
    pending: &Path,
    fields: &[(&str, &str)],
    schema_version: u32,
) -> Result<()> {
    let record_bytes = receipt(fields)?;
    let pending_receipt = sidecar(pending, ".source.json");
    write_new(&pending_receipt, &record_bytes)?;
    // The checked staging file is never a selectable final archive. Publish
    // its recovery receipt first, without replacing anything supplied. If the
    // process stops here, the next invocation can finish this exact generation.
    fs::hard_link(&pending_receipt, sidecar(output, ".source.json"))?;
    fs::File::open(parent(output))?.sync_all()?;
    let record = json_text::parse(std::str::from_utf8(&record_bytes)?, true).map_err(invalid)?;
    finish_archive(key, output, &record, schema_version)?;
    fs::remove_file(&pending_receipt)?;
    Ok(())
}

/// Adopt an archive a native builder wrote into the private config. Its own
/// bytes take the place of a source pickle: the managed name carries the
/// archive digest, and the receipt records `origin: direct` and no source.
fn install_direct(
    key: &str,
    archive: &Path,
    doc: &ConfigDocument,
    repo: &Path,
    selected_profile: &str,
) -> Result<PathBuf> {
    let schema_version = check_direct(key, archive)?;
    let schema_text = schema_version.to_string();
    let archive_hash = sha256_file(archive)?;
    let own_record_path = sidecar(archive, ".source.json");
    if own_record_path.exists() {
        // A managed archive selected again (a direct one passed through by
        // `prepare`, or a completed conversion) is reused at its own path.
        let record = json(&own_record_path)?;
        if field(&record, "key")? != key
            || field(&record, "format")? != "PRIVARKV-v1"
            || field(&record, "profile")? != selected_profile
            || !same_file(Path::new(field(&record, "archive")?), archive)
            || field(&record, "archive_sha256")? != archive_hash
        {
            return Err(invalid(format!(
                "{}: native archive differs from its managed receipt; original retained",
                archive.display()
            )));
        }
        eprintln!("{key}: reuse managed native archive");
        return Ok(archive.to_owned());
    }
    let output = target(key, &archive_hash, doc, repo, schema_version)?;
    let record_path = sidecar(&output, ".source.json");
    if output.exists() || record_path.exists() {
        let record = json(&record_path)?;
        if record.get("origin").and_then(Json::as_str) != Some("direct")
            || field(&record, "key")? != key
            || field(&record, "format")? != "PRIVARKV-v1"
            || field(&record, "profile")? != selected_profile
            || field(&record, "archive")? != text(&output)?
            || field(&record, "archive_sha256")? != archive_hash
            || field(&record, "schema_version")? != schema_text
        {
            return Err(invalid(format!(
                "{}: conflicting cache receipt; original retained",
                output.display()
            )));
        }
        finish_archive(key, &output, &record, schema_version)?;
        nonempty(&output)?;
        if sha256_file(&output)? != archive_hash {
            return Err(invalid(format!(
                "{}: archive checksum differs from its direct publication; original retained",
                output.display()
            )));
        }
        validate(key, &output, schema_version)?;
        eprintln!("{key}: reuse published direct native archive");
        return Ok(output);
    }
    // Copy rather than link: the builder's file stays its own, and the managed
    // generation cannot change if a builder later rewrites that path.
    let pending = temporary(&output, "pending");
    if let Err(error) = fs::copy(archive, &pending) {
        let _ = fs::remove_file(&pending);
        return Err(error.into());
    }
    if sha256_file(&pending)? != archive_hash {
        let _ = fs::remove_file(&pending);
        return Err(invalid(format!(
            "{} changed while it was copied; runtime config was not published",
            archive.display()
        )));
    }
    let archive_bytes = fs::metadata(&pending)?.len().to_string();
    let fields = [
        ("format", "PRIVARKV-v1"),
        ("origin", "direct"),
        ("key", key),
        ("profile", selected_profile),
        ("archive", text(&output)?),
        ("archive_sha256", &archive_hash),
        ("archive_bytes", &archive_bytes),
        ("schema_version", &schema_text),
        ("pending_archive", text(&pending)?),
    ];
    publish_generation(key, &output, &pending, &fields, schema_version)?;
    Ok(output)
}

fn finish_archive(key: &str, output: &Path, record: &Json, schema_version: u32) -> Result<()> {
    if output.exists() {
        return Ok(());
    }
    let pending = Path::new(field(record, "pending_archive")?);
    if parent(pending) != parent(output)
        || !text(Path::new(
            pending
                .file_name()
                .ok_or_else(|| invalid("pending archive lacks filename"))?,
        ))?
        .starts_with(&format!(
            "{}.",
            text(Path::new(
                output
                    .file_name()
                    .ok_or_else(|| invalid("archive lacks filename"))?
            ))?
        ))
        || pending.extension().and_then(|e| e.to_str()) != Some("pending")
    {
        return Err(invalid("pending archive is outside the owned generation"));
    }
    nonempty(pending)?;
    if sha256_file(pending)? != field(record, "archive_sha256")? {
        return Err(invalid(
            "pending archive checksum mismatch; retained without publication",
        ));
    }
    validate(key, pending, schema_version)?;
    fs::hard_link(pending, output)?;
    fs::File::open(parent(output))?.sync_all()?;
    fs::remove_file(pending)?;
    Ok(())
}

fn publish(config: &Path, private: &Path, repo: &Path) -> Result<()> {
    publish_impl(config, private, repo, || Ok(()))
}

/// One publisher at a time creates managed generations under `cache_root`.
/// The lock is held until the returned file is dropped.
fn lock_cache_root(doc: &ConfigDocument) -> Result<fs::File> {
    let root = fs::canonicalize(required(doc, "cache_root")?)?;
    if !root.is_dir() || root.parent().is_none() {
        return Err(invalid("cache_root must be a dedicated directory"));
    }
    let lock_path = root.join(".priva-native-caches.lock");
    if fs::symlink_metadata(&lock_path).is_ok_and(|m| !m.is_file() || m.file_type().is_symlink()) {
        return Err(invalid("native cache lock is not a regular file"));
    }
    let lock = OpenOptions::new()
        .create(true)
        .truncate(false)
        .read(true)
        .write(true)
        .open(&lock_path)?;
    lock.try_lock().map_err(|e| {
        invalid(format!(
            "another native-cache publisher holds cache_root: {e}"
        ))
    })?;
    Ok(lock)
}

/// The native archive a builder can read for `key` before publication. A
/// direct archive is returned after its checks; a pickle goes through the
/// managed conversion of `install_one`, which `publish` later reuses.
fn resolve_native(key: &str, private: &Path, repo: &Path) -> Result<PathBuf> {
    if !KEYS.contains(&key) {
        return Err(invalid(format!("{key}: not a native cache key")));
    }
    let doc = document(&fs::read_to_string(private)?)?;
    let selected_profile = profile(&doc)?;
    let Some(value) = setting(&doc, key)? else {
        return Err(invalid(format!(
            "{key}: not configured in the builder config"
        )));
    };
    let source = resolve(value, repo);
    if source.extension().and_then(|s| s.to_str()) == Some("rkyv") {
        check_direct(key, &source)?;
        return Ok(source);
    }
    let _lock = lock_cache_root(&doc)?;
    install_one(key, &source, &doc, repo, &selected_profile)
}

fn unchanged_metadata(a: &fs::Metadata, b: &fs::Metadata) -> bool {
    let same = a.len() == b.len() && a.modified().ok() == b.modified().ok();
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        same && (a.dev(), a.ino(), a.ctime(), a.ctime_nsec())
            == (b.dev(), b.ino(), b.ctime(), b.ctime_nsec())
    }
    #[cfg(not(unix))]
    {
        same
    }
}

fn publish_impl(
    config: &Path,
    private: &Path,
    repo: &Path,
    after_convert: impl FnOnce() -> Result<()>,
) -> Result<()> {
    let canonical_config = fs::canonicalize(config)?;
    if canonical_config
        .components()
        .collect::<Vec<_>>()
        .windows(3)
        .any(|parts| {
            parts[0].as_os_str() == ".git"
                && parts[1].as_os_str() == "annex"
                && parts[2].as_os_str() == "objects"
        })
    {
        return Err(invalid(
            "configuration is a locked DataLad object; use an unlocked or copied configuration",
        ));
    }
    let config = canonical_config.as_path();
    let mut context = check_context(config, private, repo)?;
    let contents = fs::read_to_string(private)?;
    let doc = document(&contents)?;
    let selected_profile = profile(&doc)?;
    if selected_profile != field(&context, "profile")? {
        return Err(invalid("assembly/VEP profile changed during installation"));
    }
    let _lock = lock_cache_root(&doc)?;
    let mut updates = BTreeMap::new();
    let mut sources = Vec::new();
    for key in KEYS {
        let Some(value) = setting(&doc, key)? else {
            continue;
        };
        let source = resolve(value, repo);
        let direct = source.extension().and_then(|s| s.to_str()) == Some("rkyv");
        let output = if direct {
            install_direct(key, &source, &doc, repo, &selected_profile)?
        } else {
            install_one(key, &source, &doc, repo, &selected_profile)?
        };
        if direct || bundle(key, repo, &selected_profile).is_none_or(|(_, native)| native != output)
        {
            // A direct archive is its own builder input; recheck its bytes.
            let record = json(&sidecar(&output, ".source.json"))?;
            let digest = if direct {
                "archive_sha256"
            } else {
                "source_sha256"
            };
            sources.push((source.clone(), field(&record, digest)?.to_owned()));
        }
        context.set(
            &format!("native_{key}"),
            Json::Str(text(&output)?.to_owned()),
        );
        context.set(
            &format!("source_{key}"),
            Json::Str(text(&std::path::absolute(&source)?)?.to_owned()),
        );
        updates.insert(key.to_owned(), text(&output)?.to_owned());
    }
    after_convert()?;
    // A later conversion can take minutes. Recheck the exact input bytes for
    // every managed generation before publishing the combined configuration.
    // Small native repository bundles retain existence/nonzero-only adoption.
    let mut source_states = Vec::new();
    for (source, expected) in sources {
        let before = fs::metadata(&source)?;
        let digest = sha256_file(&source)?;
        let after = fs::metadata(&source)?;
        if digest != expected || !unchanged_metadata(&before, &after) {
            return Err(invalid(format!(
                "{} changed before runtime publication",
                source.display()
            )));
        }
        source_states.push((source, after));
    }
    let runtime = rewrite(&contents, &updates)?;
    // Builders may take hours. A concurrent user edit must not be discarded.
    check_context(config, private, repo)?;
    if fs::read_to_string(private)? != contents {
        return Err(invalid("builder config changed during conversion"));
    }
    context.set(
        "published_sha256",
        Json::Str(sha256_hex(runtime.as_bytes())),
    );
    let context_path = sidecar(private, ".cache-install.json");
    atomic_write(
        &context_path,
        (json_text::dumps_indented(&context).map_err(invalid)? + "\n").as_bytes(),
    )?;
    atomic_write_checked(config, runtime.as_bytes(), || {
        for (source, expected) in &source_states {
            if !unchanged_metadata(expected, &fs::metadata(source)?) {
                return Err(invalid(format!(
                    "{} changed during final publication checks",
                    source.display()
                )));
            }
        }
        if sha256_file(config)? != field(&context, "config_sha256")? {
            return Err(invalid(
                "original configuration changed before publication; user edit retained",
            ));
        }
        Ok(())
    })?;
    context.set("state", Json::Str("complete".to_owned()));
    atomic_write(
        &context_path,
        (json_text::dumps_indented(&context).map_err(invalid)? + "\n").as_bytes(),
    )?;
    println!(
        "Runtime config published with {} native cache paths; original builder sources retained.",
        updates.len()
    );
    Ok(())
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args == ["--help"] || args == ["-h"] {
        println!("{HELP}");
        return Ok(());
    }
    if args.len() != 7 || !["prepare", "publish", "resolve"].contains(&args[0].as_str()) {
        return Err(invalid(HELP));
    }
    let first = if args[0] == "resolve" {
        "--key"
    } else {
        "--config"
    };
    let mut options = BTreeMap::new();
    for pair in args[1..].chunks_exact(2) {
        if ![first, "--build-config", "--source-root"].contains(&pair[0].as_str())
            || options
                .insert(pair[0].as_str(), PathBuf::from(&pair[1]))
                .is_some()
        {
            return Err(invalid(HELP));
        }
    }
    let private = options.get("--build-config").ok_or_else(|| invalid(HELP))?;
    let repo = fs::canonicalize(options.get("--source-root").ok_or_else(|| invalid(HELP))?)?;
    if args[0] == "resolve" {
        let key = options.get("--key").ok_or_else(|| invalid(HELP))?;
        // The one stdout line is the archive path; progress goes to stderr.
        println!("{}", text(&resolve_native(text(key)?, private, &repo)?)?);
        return Ok(());
    }
    let config = options.get("--config").ok_or_else(|| invalid(HELP))?;
    crate::rust_cache::ensure_distinct_paths(config, private)?;
    match args[0].as_str() {
        "prepare" => prepare(config, private, &repo),
        _ => publish(config, private, &repo),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    struct Temp(PathBuf);
    impl Temp {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "priva-cache-install-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir(&path).unwrap();
            Self(path)
        }
        fn config(&self, extra: &str) -> (PathBuf, PathBuf) {
            let config = self.0.join("config.yaml");
            fs::write(&config, format!("# retain comments\nassembly: hg38\nvep_installed_version: 116\ncache_root: {}\n{extra}", self.0.display())).unwrap();
            (config, self.0.join("private.yaml"))
        }
    }
    impl Drop for Temp {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    #[test]
    fn install_roundtrip_restores_builder_inputs_and_reuses_immutable_archive() {
        let tmp = Temp::new();
        let pickle = tmp.0.join("gene.pkl");
        fs::write(&pickle, b"\x80\x04}\x94.").unwrap();
        let (config, private) = tmp.config(&format!(
            "clinvar_patho_af_stat: '{}' # keep\nother: [a, b]\n",
            pickle.display()
        ));
        prepare(&config, &private, &tmp.0).unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        let runtime = fs::read_to_string(&config).unwrap();
        assert!(runtime.contains("# keep\nother: [a, b]"));
        let doc = document(&runtime).unwrap();
        let archive = PathBuf::from(required(&doc, "clinvar_patho_af_stat").unwrap());
        assert!(archive.starts_with(tmp.0.join("ClinVar")));
        let before = fs::metadata(&archive).unwrap().modified().unwrap();
        let private2 = private.clone();
        prepare(&config, &private2, &tmp.0).unwrap();
        assert_eq!(
            setting(
                &document(&fs::read_to_string(&private2).unwrap()).unwrap(),
                "clinvar_patho_af_stat"
            )
            .unwrap(),
            pickle.to_str()
        );
        publish(&config, &private2, &tmp.0).unwrap();
        assert_eq!(before, fs::metadata(&archive).unwrap().modified().unwrap());
        // A new completed builder result gets a separate generation. The old
        // archive remains usable until the final runtime config rename.
        prepare(&config, &private2, &tmp.0).unwrap();
        fs::write(&pickle, b"\x80\x04}\x94\x8c\x04GENE\x94Ns.").unwrap();
        publish(&config, &private2, &tmp.0).unwrap();
        let updated = fs::read_to_string(&config).unwrap();
        assert_ne!(runtime, updated);
        assert_eq!(before, fs::metadata(&archive).unwrap().modified().unwrap());
    }

    fn splice_pickle(nullable: bool) -> Vec<u8> {
        fn string(bytes: &mut Vec<u8>, value: &str) {
            bytes.extend_from_slice(&[0x8c, value.len() as u8]);
            bytes.extend_from_slice(value.as_bytes());
        }
        let mut bytes = b"\x80\x04}\x8c\x01T]}(".to_vec();
        for (key, value) in [
            ("chrom", "1"),
            ("pos", "42"),
            ("ref", "A"),
            ("alt", "T"),
            ("exon", "1/2"),
            ("intron", ""),
            ("hgvsc", "T:c.1+1A>T"),
            ("consequence", "splice_region_variant"),
            ("clinvar_sig", "Pathogenic"),
            ("clinvar_review", "reviewed_by_expert_panel"),
            ("splice_ai", ""),
        ] {
            string(&mut bytes, key);
            match key {
                "pos" => bytes.extend_from_slice(b"K*"),
                "splice_ai" => bytes.push(b'}'),
                "alt" if nullable => bytes.push(b'N'),
                _ => string(&mut bytes, value),
            }
        }
        bytes.extend_from_slice(b"uas.");
        bytes
    }

    #[test]
    fn splice_versioned_publication_preserves_v1_and_reuses_v2() {
        let tmp = Temp::new();
        let source = tmp.0.join("splice.pkl");
        fs::write(&source, splice_pickle(false)).unwrap();
        let (config, private) =
            tmp.config(&format!("clinvar_splice_stat: '{}'\n", source.display()));
        let selected = || {
            PathBuf::from(
                required(
                    &document(&fs::read_to_string(&config).unwrap()).unwrap(),
                    "clinvar_splice_stat",
                )
                .unwrap(),
            )
        };
        prepare(&config, &private, &tmp.0).unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        let old = selected();
        let old_bytes = fs::read(&old).unwrap();
        let old_receipt = fs::read(sidecar(&old, ".source.json")).unwrap();
        assert!(text(&old).unwrap().ends_with(".v1.rkyv"));
        assert!(
            json(&sidecar(&old, ".source.json"))
                .unwrap()
                .get("schema_version")
                .is_none()
        );
        prepare(&config, &private, &tmp.0).unwrap();
        fs::write(&source, splice_pickle(true)).unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        let typed = selected();
        assert!(text(&typed).unwrap().ends_with(".v2.rkyv"));
        let cache = caches::clinvar_splice_stat::typed::Cache::open(text(&typed).unwrap()).unwrap();
        assert_eq!(cache.version(), 2);
        assert!(cache.plain("T").unwrap()[0].alternative.is_none());
        let record = json(&sidecar(&typed, ".source.json")).unwrap();
        assert_eq!(field(&record, "schema_version").unwrap(), "2");
        assert_eq!(field(&record, "format").unwrap(), "PRIVARKV-v1");
        let modified = fs::metadata(&typed).unwrap().modified().unwrap();
        let typed_bytes = fs::read(&typed).unwrap();
        prepare(&config, &private, &tmp.0).unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        assert_eq!(selected(), typed);
        assert_eq!(fs::metadata(&typed).unwrap().modified().unwrap(), modified);
        assert_eq!(fs::read(&typed).unwrap(), typed_bytes);
        assert_eq!(fs::read(&old).unwrap(), old_bytes);
        assert_eq!(
            fs::read(sidecar(&old, ".source.json")).unwrap(),
            old_receipt
        );
    }

    #[test]
    fn splice_recovery_checks_schema_and_preserves_failed_publications() {
        let tmp = Temp::new();
        let source = tmp.0.join("splice.pkl");
        fs::write(&source, splice_pickle(true)).unwrap();
        let (config, _) = tmp.config("");
        let doc = document(&fs::read_to_string(&config).unwrap()).unwrap();
        let install = || {
            install_one(
                "clinvar_splice_stat",
                &source,
                &doc,
                &tmp.0,
                "GRCh38/VEP116",
            )
        };
        let output = install().unwrap();
        let bytes = fs::read(&output).unwrap();
        let pending = temporary(&output, "pending");
        fs::rename(&output, &pending).unwrap();
        let receipt = sidecar(&output, ".source.json");
        let mut record = json(&receipt).unwrap();
        record.set(
            "pending_archive",
            Json::Str(text(&pending).unwrap().to_owned()),
        );
        fs::write(&receipt, json_text::dumps_indented(&record).unwrap()).unwrap();
        assert_eq!(install().unwrap(), output);
        assert_eq!(fs::read(&output).unwrap(), bytes);
        record.set("schema_version", Json::Str("1".to_owned()));
        fs::write(&receipt, json_text::dumps_indented(&record).unwrap()).unwrap();
        assert!(install().is_err());
        assert_eq!(fs::read(&output).unwrap(), bytes);
        record.set("schema_version", Json::Str("2".to_owned()));
        let wrong = tmp.0.join("text.pkl");
        fs::write(&wrong, splice_pickle(false)).unwrap();
        caches::clinvar_splice_stat::convert(text(&wrong).unwrap(), &pending).unwrap();
        fs::remove_file(&output).unwrap();
        record.set("archive_sha256", Json::Str(sha256_file(&pending).unwrap()));
        fs::write(&receipt, json_text::dumps_indented(&record).unwrap()).unwrap();
        let preserved = fs::read(&pending).unwrap();
        assert!(install().is_err());
        assert!(!output.exists());
        assert_eq!(fs::read(&pending).unwrap(), preserved);
    }

    #[test]
    fn clingen_versioned_publication_preserves_v1_and_reuses_v2() {
        let tmp = Temp::new();
        let pickle = tmp.0.join("curated.pkl");
        fs::write(&pickle, b"\x80\x04}\x8c\x01a\x8c\x03PS3s.").unwrap();
        let (config, private) = tmp.config(&format!("clingen_map: '{}'\n", pickle.display()));
        let selected = || {
            PathBuf::from(
                required(
                    &document(&fs::read_to_string(&config).unwrap()).unwrap(),
                    "clingen_map",
                )
                .unwrap(),
            )
        };
        prepare(&config, &private, &tmp.0).unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        let old = selected();
        assert!(text(&old).unwrap().ends_with(".v1.rkyv"));
        let old_bytes = fs::read(&old).unwrap();
        let old_receipt = fs::read(sidecar(&old, ".source.json")).unwrap();
        assert!(
            json(&sidecar(&old, ".source.json"))
                .unwrap()
                .get("schema_version")
                .is_none()
        );
        let independent = tmp.0.join("independent.rkyv");
        caches::clingen_map::convert(text(&pickle).unwrap(), &independent).unwrap();
        assert_eq!(fs::read(&independent).unwrap(), old_bytes);

        prepare(&config, &private, &tmp.0).unwrap();
        fs::write(&pickle, b"\x80\x04}\x8c\x01aK\x11s.").unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        let typed = selected();
        assert!(text(&typed).unwrap().ends_with(".v2.rkyv"));
        let cache = caches::clingen_map::ClingenCache::open(text(&typed).unwrap()).unwrap();
        assert_eq!(cache.schema_version(), 2);
        assert_eq!(
            cache.code("a"),
            Some(caches::clingen_map::Code::Integer(17))
        );
        assert_eq!(cache.get("a"), Some(None));
        let typed_record = json(&sidecar(&typed, ".source.json")).unwrap();
        assert_eq!(field(&typed_record, "format").unwrap(), "PRIVARKV-v1");
        assert_eq!(field(&typed_record, "schema_version").unwrap(), "2");
        let typed_bytes = fs::read(&typed).unwrap();
        let typed_receipt = fs::read(sidecar(&typed, ".source.json")).unwrap();
        let modified = fs::metadata(&typed).unwrap().modified().unwrap();
        prepare(&config, &private, &tmp.0).unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        assert_eq!(selected(), typed);
        assert_eq!(fs::metadata(&typed).unwrap().modified().unwrap(), modified);
        assert_eq!(fs::read(&typed).unwrap(), typed_bytes);
        assert_eq!(
            fs::read(sidecar(&typed, ".source.json")).unwrap(),
            typed_receipt
        );
        assert_eq!(fs::read(&old).unwrap(), old_bytes);
        assert_eq!(
            fs::read(sidecar(&old, ".source.json")).unwrap(),
            old_receipt
        );
        // An unsupported new source cannot change the published config or either archive.
        prepare(&config, &private, &tmp.0).unwrap();
        let before = fs::read(&config).unwrap();
        fs::write(&pickle, b"\x80\x04}\x8c\x01aNs.").unwrap();
        assert!(publish(&config, &private, &tmp.0).is_err());
        assert_eq!(fs::read(&config).unwrap(), before);
        assert_eq!(fs::read(&typed).unwrap(), typed_bytes);
        assert_eq!(fs::read(&old).unwrap(), old_bytes);
    }

    #[test]
    fn clingen_recovery_checks_schema_and_preserves_mislabeled_archives() {
        let tmp = Temp::new();
        let source = tmp.0.join("curated.pkl");
        fs::write(&source, b"\x80\x04}\x8c\x01aK\x11s.").unwrap();
        let (config, _) = tmp.config("");
        let doc = document(&fs::read_to_string(&config).unwrap()).unwrap();
        let install = || install_one("clingen_map", &source, &doc, &tmp.0, "GRCh38/VEP116");
        let output = install().unwrap();
        let expected = fs::read(&output).unwrap();
        let pending = temporary(&output, "pending");
        fs::rename(&output, &pending).unwrap();
        let record_path = sidecar(&output, ".source.json");
        let mut record = json(&record_path).unwrap();
        record.set(
            "pending_archive",
            Json::Str(text(&pending).unwrap().to_owned()),
        );
        fs::write(&record_path, json_text::dumps_indented(&record).unwrap()).unwrap();
        assert_eq!(install().unwrap(), output);
        assert_eq!(fs::read(&output).unwrap(), expected);
        assert!(!pending.exists());

        let valid_receipt = fs::read(&record_path).unwrap();
        record.set("schema_version", Json::Str("1".to_owned()));
        fs::write(&record_path, json_text::dumps_indented(&record).unwrap()).unwrap();
        assert!(install().is_err());
        assert_eq!(fs::read(&output).unwrap(), expected);
        fs::write(&record_path, valid_receipt).unwrap();
        record.set("schema_version", Json::Str("2".to_owned()));

        // Even a matching checksum must not let a v1 payload be published as v2.
        let wrong_source = tmp.0.join("text.pkl");
        fs::write(&wrong_source, b"\x80\x04}\x8c\x01a\x8c\x03PS3s.").unwrap();
        caches::clingen_map::convert(text(&wrong_source).unwrap(), &pending).unwrap();
        fs::remove_file(&output).unwrap();
        record.set("archive_sha256", Json::Str(sha256_file(&pending).unwrap()));
        fs::write(&record_path, json_text::dumps_indented(&record).unwrap()).unwrap();
        let retained = fs::read(&pending).unwrap();
        assert!(install().is_err());
        assert!(!output.exists());
        assert_eq!(fs::read(&pending).unwrap(), retained);

        // Likewise a valid v2 payload under an owned v1 name cannot be reused.
        let old = install_one("clingen_map", &wrong_source, &doc, &tmp.0, "GRCh38/VEP116").unwrap();
        fs::write(&old, &expected).unwrap();
        let old_record_path = sidecar(&old, ".source.json");
        let mut old_record = json(&old_record_path).unwrap();
        old_record.set("archive_sha256", Json::Str(sha256_file(&old).unwrap()));
        fs::write(
            &old_record_path,
            json_text::dumps_indented(&old_record).unwrap(),
        )
        .unwrap();
        assert!(install_one("clingen_map", &wrong_source, &doc, &tmp.0, "GRCh38/VEP116").is_err());
        assert_eq!(fs::read(&old).unwrap(), expected);
    }

    fn gene_stat_pickle(symbol: &[u8]) -> Vec<u8> {
        let mut out = b"\x80\x05}\x8c\x01G}(".to_vec();
        let text = |out: &mut Vec<u8>, name: &str| {
            out.extend_from_slice(&[0x8c, name.len() as u8]);
            out.extend_from_slice(name.as_bytes());
        };
        text(&mut out, "gene_symbol");
        out.extend_from_slice(symbol);
        text(&mut out, "ensembl_id");
        text(&mut out, "G");
        for name in caches::clinvar_gene_stat::SET_NAMES {
            text(&mut out, name);
            out.push(0x8f);
        }
        for name in ["domains", "transcripts"] {
            text(&mut out, name);
            out.push(b'}');
        }
        out.extend_from_slice(b"us.");
        out
    }

    #[test]
    fn gene_statistics_versioned_publication_preserves_v1_and_reuses_v2() {
        let tmp = Temp::new();
        let source = tmp.0.join("curated.pkl");
        fs::write(&source, gene_stat_pickle(b"\x8c\x04TEXT")).unwrap();
        let (config, private) = tmp.config(&format!("clinvar_gene_stat: '{}'\n", source.display()));
        let selected = || {
            PathBuf::from(
                required(
                    &document(&fs::read_to_string(&config).unwrap()).unwrap(),
                    "clinvar_gene_stat",
                )
                .unwrap(),
            )
        };
        prepare(&config, &private, &tmp.0).unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        let old = selected();
        assert!(text(&old).unwrap().ends_with(".v1.rkyv"));
        let old_bytes = fs::read(&old).unwrap();
        let old_receipt = fs::read(sidecar(&old, ".source.json")).unwrap();
        let independent = tmp.0.join("independent.rkyv");
        caches::clinvar_gene_stat::convert(text(&source).unwrap(), &independent).unwrap();
        assert_eq!(fs::read(&independent).unwrap(), old_bytes);
        assert!(
            json(&sidecar(&old, ".source.json"))
                .unwrap()
                .get("schema_version")
                .is_none()
        );

        prepare(&config, &private, &tmp.0).unwrap();
        fs::write(&source, gene_stat_pickle(b"G\x80\0\0\0\0\0\0\0")).unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        let typed = selected();
        assert!(text(&typed).unwrap().ends_with(".v2.rkyv"));
        let cache =
            caches::clinvar_gene_stat::ClinvarGeneCache::open(text(&typed).unwrap()).unwrap();
        assert_eq!(cache.schema_version(), 2);
        assert_eq!(
            cache.plain("G").unwrap().gene_symbol,
            caches::clinvar_gene_stat::Symbol::Float(1 << 63)
        );
        let receipt = json(&sidecar(&typed, ".source.json")).unwrap();
        assert_eq!(field(&receipt, "format").unwrap(), "PRIVARKV-v1");
        assert_eq!(field(&receipt, "schema_version").unwrap(), "2");
        let typed_bytes = fs::read(&typed).unwrap();
        let typed_receipt = fs::read(sidecar(&typed, ".source.json")).unwrap();
        let modified = fs::metadata(&typed).unwrap().modified().unwrap();
        prepare(&config, &private, &tmp.0).unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        assert_eq!(selected(), typed);
        assert_eq!(fs::metadata(&typed).unwrap().modified().unwrap(), modified);
        assert_eq!(fs::read(&typed).unwrap(), typed_bytes);
        assert_eq!(
            fs::read(sidecar(&typed, ".source.json")).unwrap(),
            typed_receipt
        );
        assert_eq!(fs::read(&old).unwrap(), old_bytes);
        assert_eq!(
            fs::read(sidecar(&old, ".source.json")).unwrap(),
            old_receipt
        );

        prepare(&config, &private, &tmp.0).unwrap();
        let before = fs::read(&config).unwrap();
        fs::write(
            &source,
            gene_stat_pickle(b"cposix\nsystem\n\x8c\x04echo\x85R"),
        )
        .unwrap();
        assert!(publish(&config, &private, &tmp.0).is_err());
        assert_eq!(fs::read(&config).unwrap(), before);
        assert_eq!(fs::read(&old).unwrap(), old_bytes);
        assert_eq!(fs::read(&typed).unwrap(), typed_bytes);
    }

    #[test]
    fn gene_statistics_recovery_validates_schema_and_keeps_failed_payloads() {
        let tmp = Temp::new();
        let source = tmp.0.join("curated.pkl");
        fs::write(&source, gene_stat_pickle(b"G\x80\0\0\0\0\0\0\0")).unwrap();
        let (config, _) = tmp.config("");
        let doc = document(&fs::read_to_string(&config).unwrap()).unwrap();
        let install = || install_one("clinvar_gene_stat", &source, &doc, &tmp.0, "GRCh38/VEP116");
        let output = install().unwrap();
        let expected = fs::read(&output).unwrap();
        let pending = temporary(&output, "pending");
        fs::rename(&output, &pending).unwrap();
        let receipt = sidecar(&output, ".source.json");
        let mut record = json(&receipt).unwrap();
        record.set(
            "pending_archive",
            Json::Str(text(&pending).unwrap().to_owned()),
        );
        fs::write(&receipt, json_text::dumps_indented(&record).unwrap()).unwrap();
        assert_eq!(install().unwrap(), output);
        assert_eq!(fs::read(&output).unwrap(), expected);
        assert!(!pending.exists());
        let good_receipt = fs::read(&receipt).unwrap();
        record.set("schema_version", Json::Str("1".into()));
        fs::write(&receipt, json_text::dumps_indented(&record).unwrap()).unwrap();
        assert!(install().is_err());
        assert_eq!(fs::read(&output).unwrap(), expected);
        fs::write(&receipt, good_receipt).unwrap();
        record.set("schema_version", Json::Str("2".into()));

        let text_source = tmp.0.join("text.pkl");
        fs::write(&text_source, gene_stat_pickle(b"\x8c\x04TEXT")).unwrap();
        caches::clinvar_gene_stat::convert(text(&text_source).unwrap(), &pending).unwrap();
        fs::remove_file(&output).unwrap();
        record.set("archive_sha256", Json::Str(sha256_file(&pending).unwrap()));
        fs::write(&receipt, json_text::dumps_indented(&record).unwrap()).unwrap();
        let retained = fs::read(&pending).unwrap();
        assert!(install().is_err());
        assert!(!output.exists());
        assert_eq!(fs::read(&pending).unwrap(), retained);

        let old = install_one(
            "clinvar_gene_stat",
            &text_source,
            &doc,
            &tmp.0,
            "GRCh38/VEP116",
        )
        .unwrap();
        fs::write(&old, &expected).unwrap();
        let old_receipt = sidecar(&old, ".source.json");
        let mut record = json(&old_receipt).unwrap();
        record.set("archive_sha256", Json::Str(sha256_file(&old).unwrap()));
        fs::write(&old_receipt, json_text::dumps_indented(&record).unwrap()).unwrap();
        assert!(
            install_one(
                "clinvar_gene_stat",
                &text_source,
                &doc,
                &tmp.0,
                "GRCh38/VEP116"
            )
            .is_err()
        );
        assert_eq!(fs::read(&old).unwrap(), expected);
    }

    #[test]
    fn gene_af_versioned_publication_preserves_v1_and_reuses_v2() {
        let tmp = Temp::new();
        let pickle = tmp.0.join("curated.pkl");
        fs::write(&pickle, b"\x80\x04}\x8c\x01a}\x8c\x02afK\x01ss.").unwrap();
        let (config, private) =
            tmp.config(&format!("clinvar_patho_af_stat: '{}'\n", pickle.display()));
        let selected = || {
            PathBuf::from(
                required(
                    &document(&fs::read_to_string(&config).unwrap()).unwrap(),
                    "clinvar_patho_af_stat",
                )
                .unwrap(),
            )
        };
        prepare(&config, &private, &tmp.0).unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        let old = selected();
        assert!(text(&old).unwrap().ends_with(".v1.rkyv"));
        let old_bytes = fs::read(&old).unwrap();
        let old_receipt = fs::read(sidecar(&old, ".source.json")).unwrap();
        assert!(
            json(&sidecar(&old, ".source.json"))
                .unwrap()
                .get("schema_version")
                .is_none()
        );
        let independent = tmp.0.join("independent.rkyv");
        caches::clinvar_gene_af::convert(text(&pickle).unwrap(), &independent).unwrap();
        assert_eq!(fs::read(&independent).unwrap(), old_bytes);

        prepare(&config, &private, &tmp.0).unwrap();
        fs::write(&pickle, b"\x80\x04}\x8c\x01a}\x8c\x02af\x88ss.").unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        let typed = selected();
        assert!(text(&typed).unwrap().ends_with(".v2.rkyv"));
        let cache =
            caches::clinvar_gene_af::ClinvarGeneAfCache::open(text(&typed).unwrap()).unwrap();
        assert_eq!(cache.schema_version(), 2);
        assert!(matches!(
            cache
                .get("a")
                .unwrap()
                .unwrap()
                .iter()
                .next()
                .unwrap()
                .value,
            caches::clinvar_gene_af::ValueRef::Bool(true)
        ));
        let typed_record = json(&sidecar(&typed, ".source.json")).unwrap();
        assert_eq!(field(&typed_record, "format").unwrap(), "PRIVARKV-v1");
        assert_eq!(field(&typed_record, "schema_version").unwrap(), "2");
        let typed_bytes = fs::read(&typed).unwrap();
        let typed_receipt = fs::read(sidecar(&typed, ".source.json")).unwrap();
        let modified = fs::metadata(&typed).unwrap().modified().unwrap();
        prepare(&config, &private, &tmp.0).unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        assert_eq!(selected(), typed);
        assert_eq!(fs::metadata(&typed).unwrap().modified().unwrap(), modified);
        assert_eq!(fs::read(&typed).unwrap(), typed_bytes);
        assert_eq!(
            fs::read(sidecar(&typed, ".source.json")).unwrap(),
            typed_receipt
        );
        assert_eq!(fs::read(&old).unwrap(), old_bytes);
        assert_eq!(
            fs::read(sidecar(&old, ".source.json")).unwrap(),
            old_receipt
        );
        // An unsupported new source cannot change the published config or either archive.
        prepare(&config, &private, &tmp.0).unwrap();
        let before = fs::read(&config).unwrap();
        fs::write(&pickle, b"\x80\x04}\x8c\x01a}\x8c\x02afNss.").unwrap();
        assert!(publish(&config, &private, &tmp.0).is_err());
        assert_eq!(fs::read(&config).unwrap(), before);
        assert_eq!(fs::read(&typed).unwrap(), typed_bytes);
        assert_eq!(fs::read(&old).unwrap(), old_bytes);
    }

    #[test]
    fn gene_af_recovery_checks_schema_and_preserves_mislabeled_archives() {
        let tmp = Temp::new();
        let source = tmp.0.join("curated.pkl");
        fs::write(&source, b"\x80\x04}\x8c\x01a}\x8c\x02af\x88ss.").unwrap();
        let (config, _) = tmp.config("");
        let doc = document(&fs::read_to_string(&config).unwrap()).unwrap();
        let install = || {
            install_one(
                "clinvar_patho_af_stat",
                &source,
                &doc,
                &tmp.0,
                "GRCh38/VEP116",
            )
        };
        let output = install().unwrap();
        let expected = fs::read(&output).unwrap();
        let pending = temporary(&output, "pending");
        fs::rename(&output, &pending).unwrap();
        let record_path = sidecar(&output, ".source.json");
        let mut record = json(&record_path).unwrap();
        record.set(
            "pending_archive",
            Json::Str(text(&pending).unwrap().to_owned()),
        );
        fs::write(&record_path, json_text::dumps_indented(&record).unwrap()).unwrap();
        assert_eq!(install().unwrap(), output);
        assert_eq!(fs::read(&output).unwrap(), expected);
        assert!(!pending.exists());

        let valid_receipt = fs::read(&record_path).unwrap();
        record.set("schema_version", Json::Str("1".to_owned()));
        fs::write(&record_path, json_text::dumps_indented(&record).unwrap()).unwrap();
        assert!(install().is_err());
        assert_eq!(fs::read(&output).unwrap(), expected);
        fs::write(&record_path, valid_receipt).unwrap();
        record.set("schema_version", Json::Str("2".to_owned()));

        // Even a matching checksum must not let a v1 payload be published as v2.
        let wrong_source = tmp.0.join("text.pkl");
        fs::write(&wrong_source, b"\x80\x04}\x8c\x01a}\x8c\x02afK\x01ss.").unwrap();
        caches::clinvar_gene_af::convert(text(&wrong_source).unwrap(), &pending).unwrap();
        fs::remove_file(&output).unwrap();
        record.set("archive_sha256", Json::Str(sha256_file(&pending).unwrap()));
        fs::write(&record_path, json_text::dumps_indented(&record).unwrap()).unwrap();
        let retained = fs::read(&pending).unwrap();
        assert!(install().is_err());
        assert!(!output.exists());
        assert_eq!(fs::read(&pending).unwrap(), retained);

        // Likewise a valid v2 payload under an owned v1 name cannot be reused.
        let old = install_one(
            "clinvar_patho_af_stat",
            &wrong_source,
            &doc,
            &tmp.0,
            "GRCh38/VEP116",
        )
        .unwrap();
        fs::write(&old, &expected).unwrap();
        let old_record_path = sidecar(&old, ".source.json");
        let mut old_record = json(&old_record_path).unwrap();
        old_record.set("archive_sha256", Json::Str(sha256_file(&old).unwrap()));
        fs::write(
            &old_record_path,
            json_text::dumps_indented(&old_record).unwrap(),
        )
        .unwrap();
        assert!(
            install_one(
                "clinvar_patho_af_stat",
                &wrong_source,
                &doc,
                &tmp.0,
                "GRCh38/VEP116"
            )
            .is_err()
        );
        assert_eq!(fs::read(&old).unwrap(), expected);
    }

    #[test]
    fn clingen_bundle_remains_v1_without_inspection_or_receipt_rewrite() {
        let tmp = Temp::new();
        fs::create_dir_all(tmp.0.join("data/clingen")).unwrap();
        let pickle = tmp.0.join("data/clingen/clingen_map.hg38.pkl.gz");
        let archive = tmp.0.join("data/clingen/clingen_map.hg38.v1.rkyv");
        fs::write(&pickle, b"nonempty bundled fixture").unwrap();
        fs::write(&archive, b"existing bundled bytes").unwrap();
        let (config, private) = tmp.config("clingen_map: data/clingen/clingen_map.hg38.pkl.gz\n");
        prepare(&config, &private, &tmp.0).unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        assert_eq!(fs::read(&archive).unwrap(), b"existing bundled bytes");
        assert!(!sidecar(&archive, ".source.json").exists());
        assert_eq!(
            required(
                &document(&fs::read_to_string(config).unwrap()).unwrap(),
                "clingen_map"
            )
            .unwrap(),
            text(&archive).unwrap()
        );
    }

    #[test]
    fn refused_publish_keeps_user_edits_and_previous_config() {
        let tmp = Temp::new();
        let (config, private) = tmp.config("");
        prepare(&config, &private, &tmp.0).unwrap();
        fs::write(&config, "user edit\n").unwrap();
        assert!(publish(&config, &private, &tmp.0).is_err());
        assert_eq!(fs::read_to_string(&config).unwrap(), "user edit\n");
        let (config, private2) = tmp.config("clinvar_gene_stat: missing.pkl\n");
        let private2 = private2.with_file_name("missing.yaml");
        let before = fs::read(&config).unwrap();
        prepare(&config, &private2, &tmp.0).unwrap();
        assert!(publish(&config, &private2, &tmp.0).is_err());
        assert_eq!(fs::read(&config).unwrap(), before);
    }

    #[test]
    fn pending_build_config_is_preserved_but_user_profile_changes_are_refused() {
        let tmp = Temp::new();
        let (config, private) = tmp.config("");
        prepare(&config, &private, &tmp.0).unwrap();
        let changed = fs::read_to_string(&private).unwrap() + "new_builder_output: new.vcf.gz\n";
        fs::write(&private, &changed).unwrap();
        prepare(&config, &private, &tmp.0).unwrap();
        assert_eq!(fs::read_to_string(&private).unwrap(), changed);
        let original = fs::read_to_string(&config).unwrap().replace("hg38", "hg19");
        fs::write(&config, original).unwrap();
        assert!(prepare(&config, &private, &tmp.0).is_err());
        assert_eq!(fs::read_to_string(&private).unwrap(), changed);
    }

    #[test]
    fn damaged_managed_archive_is_not_reused_or_overwritten() {
        let tmp = Temp::new();
        let pickle = tmp.0.join("gene.pkl");
        fs::write(&pickle, b"\x80\x04}\x94.").unwrap();
        let (config, private) =
            tmp.config(&format!("clinvar_patho_af_stat: '{}'\n", pickle.display()));
        prepare(&config, &private, &tmp.0).unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        let runtime = fs::read_to_string(&config).unwrap();
        let doc = document(&runtime).unwrap();
        let archive = PathBuf::from(required(&doc, "clinvar_patho_af_stat").unwrap());
        fs::write(&archive, "damaged").unwrap();
        prepare(&config, &private, &tmp.0).unwrap();
        assert!(publish(&config, &private, &tmp.0).is_err());
        assert_eq!(fs::read_to_string(&config).unwrap(), runtime);
        assert_eq!(fs::read_to_string(&archive).unwrap(), "damaged");
    }

    #[test]
    fn small_repository_bundle_uses_only_existence_and_nonzero_checks() {
        let tmp = Temp::new();
        fs::create_dir_all(tmp.0.join("data/InterPro")).unwrap();
        let pickle = tmp.0.join("data/InterPro/Interpro_entry_mapping.pkl.gz");
        let archive = tmp.0.join("data/InterPro/Interpro_entry_mapping.v1.rkyv");
        fs::write(&pickle, "fixture").unwrap();
        fs::write(&archive, "fixture").unwrap();
        let (config, private) =
            tmp.config("interpro_mapping_pickle: data/InterPro/Interpro_entry_mapping.pkl.gz\n");
        prepare(&config, &private, &tmp.0).unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        assert_eq!(
            setting(
                &document(&fs::read_to_string(&config).unwrap()).unwrap(),
                "interpro_mapping_pickle"
            )
            .unwrap(),
            archive.to_str()
        );
        let second = tmp.0.join("second.yaml");
        prepare(&config, &second, &tmp.0).unwrap();
        assert_eq!(
            setting(
                &document(&fs::read_to_string(&second).unwrap()).unwrap(),
                "interpro_mapping_pickle"
            )
            .unwrap(),
            pickle.to_str()
        );
        #[cfg(unix)]
        {
            let alias = tmp.0.join("selected_snapshot.pkl.gz");
            std::os::unix::fs::symlink(&pickle, &alias).unwrap();
            let alias_config = tmp.0.join("alias.yaml");
            let alias_private = tmp.0.join("alias.private.yaml");
            fs::write(&alias_config, format!("assembly: hg38\nvep_installed_version: 116\ncache_root: {}\ninterpro_mapping_pickle: {}\n",tmp.0.display(),alias.display())).unwrap();
            prepare(&alias_config, &alias_private, &tmp.0).unwrap();
            publish(&alias_config, &alias_private, &tmp.0).unwrap();
            prepare(&alias_config, &alias_private, &tmp.0).unwrap();
            assert_eq!(
                setting(
                    &document(&fs::read_to_string(&alias_private).unwrap()).unwrap(),
                    "interpro_mapping_pickle"
                )
                .unwrap(),
                alias.to_str()
            );
        }
    }

    #[test]
    fn reject_unknown_native_inputs_duplicate_keys_and_redirected_cache_dirs() {
        let tmp = Temp::new();
        fs::write(tmp.0.join("unregistered.rkyv"), "x").unwrap();
        let (config, private) = tmp.config("clinvar_gene_stat: unregistered.rkyv\n");
        assert!(prepare(&config, &private, &tmp.0).is_err());
        assert!(!private.exists());
        assert!(document("assembly: hg38\nassembly: hg19\n").is_err());
        #[cfg(unix)]
        {
            std::os::unix::fs::symlink(&tmp.0, tmp.0.join("redirect")).unwrap();
            assert!(directory(&tmp.0, &tmp.0.join("redirect/out")).is_err());
        }
    }

    #[test]
    fn basename_paths_publish_in_an_isolated_working_directory() {
        const CASE: &str = "PRIVA_INSTALL_BASENAME_TEST";
        if std::env::var_os(CASE).is_some() {
            let repo = std::env::current_dir().unwrap();
            fs::write(
                "config.yaml",
                format!(
                    "assembly: hg38\nvep_installed_version: 116\ncache_root: {}\n",
                    repo.display()
                ),
            )
            .unwrap();
            prepare(Path::new("config.yaml"), Path::new("private.yaml"), &repo).unwrap();
            publish(Path::new("config.yaml"), Path::new("private.yaml"), &repo).unwrap();
            return;
        }
        let tmp = Temp::new();
        let status = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "cache_install::tests::basename_paths_publish_in_an_isolated_working_directory",
            ])
            .current_dir(&tmp.0)
            .env(CASE, "1")
            .status()
            .unwrap();
        assert!(status.success());
    }

    #[test]
    fn each_identical_warm_cycle_starts_a_new_pending_transaction() {
        let tmp = Temp::new();
        let (config, private) = tmp.config("");
        for _ in 0..3 {
            prepare(&config, &private, &tmp.0).unwrap();
            assert!(
                json(&sidecar(&private, ".cache-install.json"))
                    .unwrap()
                    .get("state")
                    .is_none()
            );
            publish(&config, &private, &tmp.0).unwrap();
        }
        prepare(&config, &private, &tmp.0).unwrap();
        let builder = fs::read_to_string(&private).unwrap() + "new_output: finished.vcf.gz\n";
        fs::write(&private, &builder).unwrap();
        let original = fs::read_to_string(&config).unwrap() + "user_edit: keep\n";
        fs::write(&config, &original).unwrap();
        assert!(prepare(&config, &private, &tmp.0).is_err());
        assert_eq!(fs::read_to_string(&private).unwrap(), builder);
        assert_eq!(fs::read_to_string(&config).unwrap(), original);
    }

    #[test]
    fn shared_archive_keeps_each_configs_original_source_selection() {
        let tmp = Temp::new();
        let a = tmp.0.join("source_a.pkl");
        let b = tmp.0.join("source_b.pkl");
        fs::write(&a, b"\x80\x04}\x94.").unwrap();
        fs::copy(&a, &b).unwrap();
        let (config_a, private_a) =
            tmp.config(&format!("clinvar_patho_af_stat: '{}'\n", a.display()));
        let config_b = tmp.0.join("b.yaml");
        let private_b = tmp.0.join("b.private.yaml");
        fs::write(
            &config_b,
            fs::read_to_string(&config_a)
                .unwrap()
                .replace(text(&a).unwrap(), text(&b).unwrap()),
        )
        .unwrap();
        prepare(&config_a, &private_a, &tmp.0).unwrap();
        publish(&config_a, &private_a, &tmp.0).unwrap();
        prepare(&config_b, &private_b, &tmp.0).unwrap();
        publish(&config_b, &private_b, &tmp.0).unwrap();
        let doc_a = document(&fs::read_to_string(&config_a).unwrap()).unwrap();
        let doc_b = document(&fs::read_to_string(&config_b).unwrap()).unwrap();
        assert_eq!(
            setting(&doc_a, "clinvar_patho_af_stat").unwrap(),
            setting(&doc_b, "clinvar_patho_af_stat").unwrap()
        );
        prepare(&config_b, &private_b, &tmp.0).unwrap();
        assert_eq!(
            setting(
                &document(&fs::read_to_string(&private_b).unwrap()).unwrap(),
                "clinvar_patho_af_stat"
            )
            .unwrap(),
            b.to_str()
        );
        fs::remove_file(&b).unwrap();
        prepare(&config_a, &private_a, &tmp.0).unwrap();
        assert_eq!(
            setting(
                &document(&fs::read_to_string(&private_a).unwrap()).unwrap(),
                "clinvar_patho_af_stat"
            )
            .unwrap(),
            a.to_str()
        );
    }

    #[test]
    fn copied_archive_with_stale_location_receipt_is_refused() {
        let tmp = Temp::new();
        let source = tmp.0.join("gene.pkl");
        fs::write(&source, b"\x80\x04}\x94.").unwrap();
        let (config, private) =
            tmp.config(&format!("clinvar_patho_af_stat: '{}'\n", source.display()));
        prepare(&config, &private, &tmp.0).unwrap();
        publish(&config, &private, &tmp.0).unwrap();
        let original = fs::read_to_string(&config).unwrap();
        let old = PathBuf::from(
            required(&document(&original).unwrap(), "clinvar_patho_af_stat").unwrap(),
        );
        let new_root = tmp.0.join("moved");
        fs::create_dir(&new_root).unwrap();
        let changed = original.replace(
            &format!("cache_root: {}", tmp.0.display()),
            &format!("cache_root: {}", new_root.display()),
        );
        fs::write(&config, &changed).unwrap();
        prepare(&config, &private, &tmp.0).unwrap();
        let doc = document(&fs::read_to_string(&private).unwrap()).unwrap();
        let new = target(
            "clinvar_patho_af_stat",
            &sha256_file(&source).unwrap(),
            &doc,
            &tmp.0,
            1,
        )
        .unwrap();
        fs::copy(&old, &new).unwrap();
        fs::copy(sidecar(&old, ".source.json"), sidecar(&new, ".source.json")).unwrap();
        assert!(publish(&config, &private, &tmp.0).is_err());
        assert_eq!(fs::read_to_string(&config).unwrap(), changed);
    }

    #[test]
    fn interrupted_owned_generation_recovers_from_its_verified_staging_receipt() {
        let tmp = Temp::new();
        let source = tmp.0.join("gene.pkl");
        fs::write(&source, b"\x80\x04}\x94.").unwrap();
        let (config, _) = tmp.config("");
        let doc = document(&fs::read_to_string(config).unwrap()).unwrap();
        let output = install_one(
            "clinvar_patho_af_stat",
            &source,
            &doc,
            &tmp.0,
            "GRCh38/VEP116",
        )
        .unwrap();
        let expected = fs::read(&output).unwrap();
        let pending = temporary(&output, "pending");
        fs::rename(&output, &pending).unwrap();
        let record_path = sidecar(&output, ".source.json");
        let mut record = json(&record_path).unwrap();
        record.set(
            "pending_archive",
            Json::Str(text(&pending).unwrap().to_owned()),
        );
        fs::write(&record_path, json_text::dumps_indented(&record).unwrap()).unwrap();
        let restored = install_one(
            "clinvar_patho_af_stat",
            &source,
            &doc,
            &tmp.0,
            "GRCh38/VEP116",
        )
        .unwrap();
        assert_eq!(restored, output);
        assert_eq!(fs::read(&output).unwrap(), expected);
        assert!(!pending.exists());
    }

    #[cfg(unix)]
    #[test]
    fn preserve_configuration_alias_and_logical_source_link() {
        let tmp = Temp::new();
        let source = tmp.0.join("source.pkl");
        let alias = tmp.0.join("source_alias.pkl");
        fs::write(&source, b"\x80\x04}\x94.").unwrap();
        std::os::unix::fs::symlink(&source, &alias).unwrap();
        let (config, private) =
            tmp.config(&format!("clinvar_patho_af_stat: '{}'\n", alias.display()));
        let config_alias = tmp.0.join("config_alias.yaml");
        std::os::unix::fs::symlink(&config, &config_alias).unwrap();
        prepare(&config_alias, &private, &tmp.0).unwrap();
        publish(&config_alias, &private, &tmp.0).unwrap();
        assert!(config_alias.is_symlink());
        prepare(&config_alias, &private, &tmp.0).unwrap();
        assert_eq!(
            setting(
                &document(&fs::read_to_string(&private).unwrap()).unwrap(),
                "clinvar_patho_af_stat"
            )
            .unwrap(),
            alias.to_str()
        );
    }

    #[test]
    fn config_edit_during_temporary_write_is_preserved() {
        let tmp = Temp::new();
        let file = tmp.0.join("config.yaml");
        fs::write(&file, "original").unwrap();
        let before = sha256_file(&file).unwrap();
        let result = atomic_write_checked(&file, b"new runtime", || {
            fs::write(&file, "user edit")?;
            if sha256_file(&file)? != before {
                return Err(invalid("changed"));
            }
            Ok(())
        });
        assert!(result.is_err());
        assert_eq!(fs::read_to_string(&file).unwrap(), "user edit");
    }

    #[test]
    fn earlier_source_change_during_later_cache_work_blocks_publication() {
        let tmp = Temp::new();
        let source = tmp.0.join("gene.pkl");
        fs::write(&source, b"\x80\x04}\x94.").unwrap();
        let (config, private) =
            tmp.config(&format!("clinvar_patho_af_stat: '{}'\n", source.display()));
        let before = fs::read(&config).unwrap();
        prepare(&config, &private, &tmp.0).unwrap();
        let result = publish_impl(&config, &private, &tmp.0, || {
            fs::write(&source, b"\x80\x04}\x94\x8c\x04GENE\x94Ns.")?;
            Ok(())
        });
        assert!(result.is_err());
        assert_eq!(fs::read(&config).unwrap(), before);
        // The old generation remains evidence; a quiescent retry converts the
        // actual new bytes and can publish without changing source files.
        publish(&config, &private, &tmp.0).unwrap();
    }

    fn select_in_builder(private: &Path, from: &Path, to: &Path) {
        let builder = fs::read_to_string(private).unwrap();
        assert!(builder.contains(text(from).unwrap()));
        fs::write(
            private,
            builder.replace(text(from).unwrap(), text(to).unwrap()),
        )
        .unwrap();
    }

    #[test]
    fn direct_archive_publishes_once_and_passes_through_prepare() {
        let tmp = Temp::new();
        let pickle = tmp.0.join("curated.pkl");
        fs::write(&pickle, b"\x80\x04}\x8c\x01a}\x8c\x02afK\x01ss.").unwrap();
        let built = tmp.0.join("built.rkyv");
        caches::clinvar_gene_af::convert(text(&pickle).unwrap(), &built).unwrap();
        let (config, private) =
            tmp.config(&format!("clinvar_patho_af_stat: '{}'\n", pickle.display()));
        let selected = || {
            PathBuf::from(
                required(
                    &document(&fs::read_to_string(&config).unwrap()).unwrap(),
                    "clinvar_patho_af_stat",
                )
                .unwrap(),
            )
        };
        prepare(&config, &private, &tmp.0).unwrap();
        select_in_builder(&private, &pickle, &built);
        publish(&config, &private, &tmp.0).unwrap();
        let archive = selected();
        assert_ne!(archive, built);
        assert!(archive.starts_with(tmp.0.join("ClinVar")));
        assert!(text(&archive).unwrap().ends_with(".v1.rkyv"));
        assert_eq!(fs::read(&archive).unwrap(), fs::read(&built).unwrap());
        let record = json(&sidecar(&archive, ".source.json")).unwrap();
        assert_eq!(field(&record, "origin").unwrap(), "direct");
        assert_eq!(field(&record, "key").unwrap(), "clinvar_patho_af_stat");
        assert_eq!(field(&record, "profile").unwrap(), "GRCh38/VEP116");
        assert_eq!(
            field(&record, "archive_sha256").unwrap(),
            sha256_file(&built).unwrap()
        );
        assert!(record.get("source").is_none());
        let receipt = fs::read(sidecar(&archive, ".source.json")).unwrap();
        let modified = fs::metadata(&archive).unwrap().modified().unwrap();

        // The next builder view keeps selecting the archive; no pickle is restored.
        prepare(&config, &private, &tmp.0).unwrap();
        assert_eq!(
            setting(
                &document(&fs::read_to_string(&private).unwrap()).unwrap(),
                "clinvar_patho_af_stat"
            )
            .unwrap(),
            archive.to_str()
        );
        publish(&config, &private, &tmp.0).unwrap();
        assert_eq!(selected(), archive);
        // The same bytes written again by a builder reuse the published generation.
        prepare(&config, &private, &tmp.0).unwrap();
        select_in_builder(&private, &archive, &built);
        publish(&config, &private, &tmp.0).unwrap();
        assert_eq!(selected(), archive);
        assert_eq!(
            fs::metadata(&archive).unwrap().modified().unwrap(),
            modified
        );
        assert_eq!(
            fs::read(sidecar(&archive, ".source.json")).unwrap(),
            receipt
        );
    }

    #[test]
    fn damaged_wrong_kind_or_missing_direct_archives_are_refused() {
        let tmp = Temp::new();
        let pickle = tmp.0.join("gene.pkl");
        fs::write(&pickle, b"\x80\x04}\x94.").unwrap();
        let splice = tmp.0.join("splice.pkl");
        fs::write(&splice, splice_pickle(false)).unwrap();
        let wrong = tmp.0.join("splice.rkyv");
        caches::clinvar_splice_stat::convert(text(&splice).unwrap(), &wrong).unwrap();
        let damaged = tmp.0.join("damaged.rkyv");
        caches::clinvar_gene_af::convert(text(&pickle).unwrap(), &damaged).unwrap();
        let bytes = fs::read(&damaged).unwrap();
        fs::write(&damaged, &bytes[..bytes.len() - 1]).unwrap();
        let short = tmp.0.join("short.rkyv");
        fs::write(&short, "damaged").unwrap();
        let missing = tmp.0.join("missing.rkyv");
        let (config, private) =
            tmp.config(&format!("clinvar_patho_af_stat: '{}'\n", pickle.display()));
        prepare(&config, &private, &tmp.0).unwrap();
        let builder = fs::read_to_string(&private).unwrap();
        let before = fs::read(&config).unwrap();
        for (archive, message) in [
            (&wrong, "schema"),
            (&damaged, "length"),
            (&short, "too short"),
            (&missing, "does not exist"),
        ] {
            fs::write(&private, &builder).unwrap();
            select_in_builder(&private, &pickle, archive);
            let error = publish(&config, &private, &tmp.0).unwrap_err().to_string();
            assert!(error.contains(message), "{error}");
            assert_eq!(fs::read(&config).unwrap(), before);
            assert!(resolve_native("clinvar_patho_af_stat", &private, &tmp.0).is_err());
        }
        assert!(!tmp.0.join("ClinVar").exists());
    }

    #[test]
    fn resolved_conversion_is_reused_by_publish() {
        let tmp = Temp::new();
        let pickle = tmp.0.join("gene.pkl");
        fs::write(&pickle, b"\x80\x04}\x94.").unwrap();
        let (config, private) =
            tmp.config(&format!("clinvar_patho_af_stat: '{}'\n", pickle.display()));
        prepare(&config, &private, &tmp.0).unwrap();
        assert!(resolve_native("clinvar_gene_stat", &private, &tmp.0).is_err());
        assert!(resolve_native("unknown", &private, &tmp.0).is_err());
        let resolved = resolve_native("clinvar_patho_af_stat", &private, &tmp.0).unwrap();
        let bytes = fs::read(&resolved).unwrap();
        let receipt = fs::read(sidecar(&resolved, ".source.json")).unwrap();
        let modified = fs::metadata(&resolved).unwrap().modified().unwrap();
        assert_eq!(
            resolve_native("clinvar_patho_af_stat", &private, &tmp.0).unwrap(),
            resolved
        );
        publish(&config, &private, &tmp.0).unwrap();
        assert_eq!(
            setting(
                &document(&fs::read_to_string(&config).unwrap()).unwrap(),
                "clinvar_patho_af_stat"
            )
            .unwrap(),
            resolved.to_str()
        );
        assert_eq!(fs::read(&resolved).unwrap(), bytes);
        assert_eq!(
            fs::metadata(&resolved).unwrap().modified().unwrap(),
            modified
        );
        assert_eq!(
            fs::read(sidecar(&resolved, ".source.json")).unwrap(),
            receipt
        );
        let generations = fs::read_dir(resolved.parent().unwrap())
            .unwrap()
            .map(|entry| entry.unwrap().file_name())
            .filter(|name| name.to_string_lossy().ends_with(".rkyv"))
            .count();
        assert_eq!(generations, 1);
        // A converted archive named directly resolves to itself after its checks.
        prepare(&config, &private, &tmp.0).unwrap();
        select_in_builder(&private, &pickle, &resolved);
        assert_eq!(
            resolve_native("clinvar_patho_af_stat", &private, &tmp.0).unwrap(),
            resolved
        );
    }
}
