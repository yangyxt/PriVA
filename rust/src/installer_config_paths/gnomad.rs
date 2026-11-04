//! Only the original read-only dependency used during source adoption: exact
//! family names, header/assembly checks and one indexed record per chromosome.
//! No download, index creation, full stream scan or scientific regeneration.
use super::*;
use crate::vcf_variant;
use std::collections::HashMap;

fn assembly(value: &str) -> Result<&'static str> {
    match value {
        "hg19" | "GRCh37" => Ok("hg19"),
        "hg38" | "GRCh38" => Ok("hg38"),
        _ => Err(error("ValueError", "Select hg19/GRCh37 or hg38/GRCh38")),
    }
}
fn chromosomes() -> Vec<String> {
    (1..=22)
        .map(|n| format!("chr{n}"))
        .chain(["chrX".into(), "chrY".into()])
        .collect()
}
pub(crate) fn default_path(cache_dir: &Path, selected: &str) -> Result<PathBuf> {
    let a = assembly(selected)?;
    Ok(join(
        Path::new(&absolute(cache_dir)?),
        if a == "hg19" {
            "gnomad.joint.v4.1.sites.hg19.chrX.vcf.gz"
        } else {
            "gnomad.joint.v4.1.sites.chrX.vcf.bgz"
        },
    ))
}
pub(crate) fn family_paths(path: &Path) -> Result<Vec<(String, PathBuf)>> {
    let path = path_value(&absolute(path)?);
    let filename = name(&path);
    let (prefix, suffix) = if let Some(v) = filename.strip_suffix(".chrX.vcf.gz") {
        (v, ".gz")
    } else if let Some(v) = filename.strip_suffix(".chrX.vcf.bgz") {
        (v, ".bgz")
    } else {
        ("", "")
    };
    require(
        !prefix.is_empty() && !prefix.contains('\n'),
        "Selected gnomAD filename must retain the caller's .chrX.vcf.gz or .chrX.vcf.bgz template; refusing to substitute another dataset",
    )?;
    require(
        !path.components().any(|p| {
            let s = p.as_os_str().to_string_lossy();
            s.to_lowercase().contains(".partial") || s.starts_with(".building")
        }),
        "Selected gnomAD family must not be a partial/building cache",
    )?;
    Ok(chromosomes()
        .into_iter()
        .map(|c| {
            let p = path.with_file_name(format!("{prefix}.{c}.vcf{suffix}"));
            (c, p)
        })
        .collect())
}
fn reference_lengths(reference: Value<'_>) -> Result<HashMap<String, i128>> {
    let mut result = HashMap::new();
    if reference.null() {
        return Ok(result);
    }
    let index = PathBuf::from(format!("{}.fai", reference.text()));
    if !is_file(&index)? {
        return Ok(result);
    }
    let text = files::read_text(&index)?;
    for line in text.split_inclusive('\n') {
        let fields = line
            .trim_end_matches(['\r', '\n'])
            .split('\t')
            .collect::<Vec<_>>();
        require(
            fields.len() >= 2
                && !fields[1].is_empty()
                && fields[1].chars().all(|c| c.is_ascii_digit()),
            format!("Malformed selected reference index: {}", index.display()),
        )?;
        let key = format!("chr{}", fields[0].strip_prefix("chr").unwrap_or(fields[0]));
        if !chromosomes().contains(&key) {
            continue;
        }
        let length = fields[1].parse().unwrap_or(i128::MAX);
        require(
            result.get(&key).is_none_or(|v| *v == length),
            "Conflicting selected reference contig lengths",
        )?;
        result.insert(key, length);
    }
    Ok(result)
}
struct Header {
    typed: vcf_variant::Header,
    lengths: HashMap<String, i128>,
}
fn header(path: &Path) -> Result<Header> {
    let raw = crate::clinvar_contract::header::header_text(&files::native(path)?)?;
    let text = crate::python_io::decode_utf8(&raw).map_err(|e| error("UnicodeDecodeError", e))?;
    let mut lengths = HashMap::new();
    for line in text.lines() {
        if let Some(attrs) = vcf_variant::attributes(line, "##contig=<") {
            if let Some(id) = attrs.get("ID") {
                if let Some(length) = attrs
                    .get("length")
                    .and_then(|v| v.parse::<i128>().ok())
                    .filter(|v| *v != 0)
                {
                    lengths.entry(id.clone()).or_insert(length);
                }
            }
        }
    }
    let typed = vcf_variant::read_header(&mut std::io::Cursor::new(raw))?;
    Ok(Header { typed, lengths })
}
/// `validate_family(...)` with the reference as `gnomad_cache.py validate`
/// passes it: a path, or none.
pub(crate) fn validate_selected(
    path: &Path,
    selected: &str,
    reference: Option<&str>,
) -> Result<Json> {
    validate_family(
        path,
        selected,
        reference.map_or(Value::Missing, Value::Text),
    )
}
pub(super) fn validate_family(path: &Path, selected: &str, reference: Value<'_>) -> Result<Json> {
    let a = assembly(selected)?;
    let family = family_paths(path)?;
    let mut tracked = vec![];
    for (_, path) in &family {
        require(
            is_file(path)? && metadata(path)?.len() > 0,
            format!(
                "Selected gnomAD family member is missing/empty: {}; restore that selected family",
                path.display()
            ),
        )?;
        tracked.push((path.clone(), files::cache_identity(path)?));
        let mut indexes = vec![];
        for suffix in [".tbi", ".csi"] {
            let p = files::suffix(path, suffix);
            if is_file(&p)? {
                indexes.push(p);
            }
        }
        let mut valid = !indexes.is_empty();
        for p in &indexes {
            valid &= metadata(p)?.len() > 0;
        }
        require(
            valid,
            format!(
                "Selected gnomAD family requires its nonempty matching index: {}",
                path.display()
            ),
        )?;
        for p in indexes {
            let id = files::cache_identity(&p)?;
            tracked.push((p, id));
        }
    }
    let lengths = reference_lengths(reference)?;
    let chr1 = if a == "hg19" { 249250621 } else { 248956422 };
    require(
        lengths.is_empty() || lengths.get("chr1") == Some(&chr1),
        "Selected reference index does not identify the requested gnomAD assembly",
    )?;
    let result = (|| {
        let anchor = header(&family[0].1)?;
        let mut expected = anchor.lengths;
        require(
            expected.get("chr1") == Some(&chr1),
            "Selected gnomAD chr1 header does not identify the requested assembly",
        )?;
        expected.extend(lengths);
        let mut required = [
            "AC_joint",
            "AN_joint",
            "AF_joint",
            "nhomalt_joint",
            "AC_grpmax_joint",
            "AF_grpmax_joint",
            "AN_grpmax_joint",
            "nhomalt_grpmax_joint",
            "AC_joint_XX",
            "AN_joint_XX",
            "AC_joint_XY",
            "AN_joint_XY",
            "AF_joint_XX",
            "AF_joint_XY",
            "nhomalt_joint_XX",
            "nhomalt_joint_XY",
        ]
        .into_iter()
        .map(str::to_owned)
        .collect::<Vec<_>>();
        for pop in [
            "nfe",
            "eas",
            "afr",
            "amr",
            "asj",
            "fin",
            "sas",
            "mid",
            "remaining",
        ] {
            for field in [
                format!("AF_joint_{pop}_XX"),
                format!("AF_joint_{pop}_XY"),
                format!("nhomalt_joint_{pop}_XX"),
            ] {
                required.push(field);
            }
        }
        let mut records = vec![];
        for (chromosome, path) in &family {
            let h = header(path)?;
            let missing = required
                .iter()
                .filter(|k| !h.typed.info.contains_key(*k))
                .cloned()
                .collect::<Vec<_>>();
            require(
                missing.is_empty(),
                format!(
                    "Selected gnomAD header lacks caller-required INFO fields in {}: {}",
                    path.display(),
                    missing.join(",")
                ),
            )?;
            for key in &required {
                let k = h.typed.info[key].kind;
                require(
                    if key.starts_with("AF_") {
                        matches!(k, vcf_variant::InfoType::Float)
                    } else {
                        matches!(k, vcf_variant::InfoType::Integer)
                    },
                    format!(
                        "Selected gnomAD INFO/{key} has incompatible type in {}",
                        path.display()
                    ),
                )?;
            }
            require(
                h.typed.contigs.contains(chromosome) && expected.contains_key(chromosome),
                format!(
                    "Selected gnomAD header/reference lacks {chromosome} assembly length: {}",
                    path.display()
                ),
            )?;
            require(
                h.lengths.get(chromosome) == expected.get(chromosome),
                format!(
                    "Selected gnomAD {chromosome} header differs from the selected assembly/reference: {}",
                    path.display()
                ),
            )?;
            let native_path = files::native(path)?;
            let index = vcf_variant::index(&native_path)
                .map_err(|_| error("ValueError", "fetch requires an index"))?;
            let mut first = None;
            let fetched = crate::tabix_index::fetch_contig_lines(
                &native_path,
                &index,
                chromosome,
                &mut |line| {
                    first = Some(line.to_vec());
                    Err("first-record-ready".into())
                },
            );
            if let Err(cause) = fetched {
                if cause != "first-record-ready" {
                    return Err(error("OSError", "unable to fetch next record"));
                }
            }
            if let Some(line) = first {
                let record = vcf_variant::decode_record(&line, &h.typed, &mut h.typed.info.clone())
                    .map_err(|_| error("OSError", "unable to parse next record"))?;
                require(
                    record.fields.first().copied() == Some(chromosome)
                        && record.position > 0
                        && record.position as i128 <= expected[chromosome],
                    format!(
                        "Selected gnomAD indexed record is incompatible with its chromosome: {}",
                        path.display()
                    ),
                )?;
            }
            records.push(Json::Object(vec![
                ("chromosome".into(), Json::Str(chromosome.clone())),
                ("path".into(), Json::Str(path.display().to_string())),
                ("bytes".into(), Json::Int(metadata(path)?.len() as i128)),
            ]));
        }
        Ok(records)
    })();
    let records = result.map_err(|e: Box<dyn std::error::Error>| {
        let text = e.to_string();
        let message = text
            .split_once(": ")
            .filter(|(k, _)| {
                matches!(
                    *k,
                    "ValueError" | "OSError" | "KeyError" | "FileNotFoundError" | "PermissionError"
                )
            })
            .map(|(_, m)| m)
            .unwrap_or(&text);
        error(
            "ValueError",
            format!("Cannot reuse selected gnomAD family: {message}"),
        )
    })?;
    for (path, before) in tracked {
        require(
            json_text::py_equal(&files::cache_identity(&path)?, &before),
            "Selected gnomAD family changed during inspection",
        )?;
    }
    Ok(Json::Object(vec![("assembly".into(),Json::Str(a.into())),("chr_x".into(),Json::Str(family[22].1.display().to_string())),("members".into(),Json::Array(records)),("scientific_acceptance".into(),Json::Bool(false)),("scope".into(),Json::Str("Original caller headers, coherent assembly and indexed access; no data/release substitution.".into()))]))
}
