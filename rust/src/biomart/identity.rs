use super::*;
use crate::clinvar_contract::resolve_non_strict;
use crate::python_path::{is_file, tempfile_parent};
use std::io::Write;
use std::os::unix::fs::MetadataExt;

pub(super) fn is_bundled(context: &Context, mapping: &Path) -> Result<bool> {
    let root = context.source_root()?;
    let logical = tempfile_parent(mapping)?;
    let basename = name(&logical);
    Ok(source_for_basename(&basename).is_some()
        && logical == tempfile_parent(&root.join("data/hgnc").join(&basename))?
        && !is_link(&logical)?
        && resolve_non_strict(&logical)?
            == resolve_non_strict(&root)?.join("data/hgnc").join(basename))
}
pub(crate) fn mapping_identity(
    context: &Context,
    mapping: &Path,
    inputs: &[PathBuf],
) -> Result<Json> {
    let sidecar = provenance_path(mapping);
    if is_bundled(context, mapping)? {
        require(
            is_file(mapping)? && metadata(mapping)?.len() > 0,
            format!("Missing or empty bundled mapping: {}", mapping.display()),
        )?;
    } else if exists(&sidecar)?
        || is_link(&sidecar)?
        || source_for_basename(&name(mapping)).is_some()
    {
        verify::pinned(mapping, None)?;
    }
    let mapping_hash = sha256(mapping)?;
    let mapping_bytes = metadata(mapping)?.len();
    let provenance_hash = if exists(&sidecar)? {
        string(sha256(&sidecar)?)
    } else {
        Json::Null
    };
    let mut anchors = Vec::with_capacity(inputs.len());
    for path in inputs {
        anchors.push(object([
            ("name", string(name(path))),
            ("bytes", Json::Int(metadata(path)?.len() as i128)),
            ("sha256", string(sha256(path)?)),
        ]));
    }
    Ok(object([
        ("mapping_sha256", string(mapping_hash)),
        ("mapping_bytes", Json::Int(mapping_bytes as i128)),
        ("provenance_sha256", provenance_hash),
        ("anchor_inputs", Json::Array(anchors)),
    ]))
}
/// Exactly the source's `(OSError, ValueError)` catch family. The mature XML
/// ParseError is deliberately outside this wrapper, as are unexpected errors.
fn caught_message(cause: &dyn std::error::Error) -> Option<String> {
    let rendered = cause.to_string();
    [
        "OSError",
        "FileNotFoundError",
        "PermissionError",
        "FileExistsError",
        "IsADirectoryError",
        "NotADirectoryError",
        "BlockingIOError",
        "BrokenPipeError",
        "BadGzipFile",
        "ValueError",
        "UnicodeDecodeError",
        "UnicodeEncodeError",
        "JSONDecodeError",
    ]
    .into_iter()
    .find_map(|name| {
        rendered
            .strip_prefix(&format!("{name}: "))
            .map(str::to_owned)
    })
}
pub(crate) fn check_receipt(
    context: &Context,
    artifact: &Path,
    mapping: &Path,
    inputs: &[PathBuf],
) -> Result<()> {
    let result = (|| {
        let recorded = receipt::artifact_receipt(artifact)?;
        require(
            recorded.matches(&mapping_identity(context, mapping, inputs)?)?,
            "mapping content/provenance or HCSeeker/RMC anchors changed",
        )
    })();
    result.map_err(|e|match caught_message(e.as_ref()) {
        Some(message)=>error("ValueError",format!("Cannot reuse or overwrite existing {}: missing/stale NM-to-ENST receipt ({message}). Select a new output path or explicitly move the existing artifact and its receipt before rebuilding.",artifact.display())),
        None=>e,
    })
}
pub(super) fn atomic_json(path: &Path, value: &Json) -> Result<()> {
    let temporary = PathBuf::from(format!("{}.part", path.display()));
    require(
        !is_link(&temporary)?,
        format!(
            "Refusing to replace a symlink used as a staging file: {}",
            temporary.display()
        ),
    )?;
    let mut stream = File::create(&temporary).map_err(|e| os_error(e, &temporary, None))?;
    // Python opens/truncates the staging file before attempting serialization.
    let data = json_text::dumps_with_options(
        value,
        json_text::DumpsOptions {
            indent: Some(2),
            sort_keys: true,
            ensure_ascii: true,
            compact: false,
            allow_nan: true,
        },
    )
    .map_err(|e| error("ValueError", e))?;
    stream
        .write_all(format!("{data}\n").as_bytes())
        .map_err(|e| os_error(e, &temporary, None))?;
    drop(stream);
    fs::rename(&temporary, path).map_err(|e| os_error(e, &temporary, Some(path)))?;
    Ok(())
}
pub(super) fn write_receipt(
    context: &Context,
    artifact: &Path,
    mapping: &Path,
    inputs: &[PathBuf],
    expected: &str,
) -> Result<()> {
    require(
        is_file(artifact)? && metadata(artifact)?.len() > 0,
        format!(
            "Scientific command did not create a nonempty output; no mapping receipt written: {}",
            artifact.display()
        ),
    )?;
    let expected = parse_json(expected)?;
    require(
        equal(&expected, &mapping_identity(context, mapping, inputs)?),
        "Mapping content/provenance changed during scientific generation; output has no valid mapping receipt",
    )?;
    atomic_json(
        &PathBuf::from(format!("{}.nm_to_enst.json", artifact.display())),
        &object([
            ("mapping_identity", expected),
            (
                "artifact_bytes",
                Json::Int(metadata(artifact)?.len() as i128),
            ),
            ("artifact_sha256", string(sha256(artifact)?)),
        ]),
    )
}
pub(super) fn changed_target(
    identity: &Json,
    cache_root: &Path,
    output_dir: &Path,
    assembly: &str,
    kind: &str,
) -> Result<PathBuf> {
    require(
        matches!(assembly, "hg19" | "GRCh37" | "hg38" | "GRCh38")
            && matches!(kind, "motifs" | "pm1"),
        "Select an explicit supported assembly and motif/PM1 artifact kind",
    )?;
    let assembly = if matches!(assembly, "hg19" | "GRCh37") {
        "hg19"
    } else {
        "hg38"
    };
    let filename = if kind == "motifs" {
        format!("AlphaMissense_{assembly}.kde.pkl")
    } else {
        format!("pm1_regions.{assembly}.pkl")
    };
    let serialized = json_text::dumps_with_options(
        identity,
        json_text::DumpsOptions {
            indent: None,
            sort_keys: true,
            ensure_ascii: true,
            compact: true,
            allow_nan: true,
        },
    )
    .map_err(|e| error("ValueError", e))?;
    let relative = Path::new("nm-to-enst")
        .join(assembly)
        .join(kind)
        .join(hash(serialized.as_bytes()))
        .join(filename);
    require(
        cache_root.is_absolute() && output_dir.is_absolute(),
        "Changed mapping outputs require absolute cache_root and output_dir",
    )?;
    let logical_root = tempfile_parent(cache_root)?;
    let root =
        crate::python_io::resolve_strict(&logical_root).map_err(|(e, p)| os_error(e, &p, None))?;
    require(
        root.is_dir() && root.parent().is_some(),
        "Changed mapping outputs require a configured cache_root directory",
    )?;
    let logical_output = tempfile_parent(output_dir)?;
    let target = if let Ok(tail) = logical_output.strip_prefix(&logical_root) {
        root.join(tail).join(relative)
    } else if logical_output.starts_with(&root) {
        logical_output.join(relative)
    } else {
        return Err(error(
            "ValueError",
            "Changed mapping output_dir escapes cache_root",
        ));
    };
    let mut current = root.clone();
    for part in target
        .strip_prefix(&root)
        .expect("managed target")
        .components()
    {
        current.push(part);
        require(
            !is_link(&current)?,
            format!(
                "Linked managed mapping output is not allowed: {}",
                current.display()
            ),
        )?;
        require(
            current == target || !exists(&current)? || current.is_dir(),
            format!(
                "Managed mapping output parent is not a directory: {}",
                current.display()
            ),
        )?;
    }
    Ok(target)
}
pub(crate) fn artifact_route(
    context: &Context,
    artifact: &Path,
    mapping: &Path,
    cache_root: &Path,
    output_dir: &Path,
    assembly: &str,
    kind: &str,
    inputs: &[PathBuf],
) -> Result<PathBuf> {
    let recorded=receipt::artifact_receipt(artifact).map_err(|e|match caught_message(e.as_ref()) {
        Some(message)=>error("ValueError",format!("Cannot reuse or overwrite existing {}: missing/stale NM-to-ENST receipt ({message}). Retain the artifact and recover its receipt explicitly.",artifact.display())),None=>e,
    })?;
    require(
        is_file(mapping)? && metadata(mapping)?.len() > 0,
        format!("Missing or empty selected mapping: {}", mapping.display()),
    )?;
    let identity = mapping_identity(context, mapping, inputs)?;
    if recorded.matches(&identity)? {
        return Ok(artifact.to_path_buf());
    }
    let target = changed_target(&identity, cache_root, output_dir, assembly, kind)?;
    let receipt = PathBuf::from(format!("{}.nm_to_enst.json", target.display()));
    require(
        !is_link(&receipt)?,
        format!(
            "Linked managed mapping receipt is not allowed: {}",
            receipt.display()
        ),
    )?;
    if exists(&target)? || exists(&receipt)? {
        require(
            !is_file(&target)? || metadata(&target)?.nlink() == 1,
            format!(
                "Shared managed mapping artifact is not allowed: {}",
                target.display()
            ),
        )?;
        check_receipt(context, &target, mapping, inputs)?;
    }
    Ok(target)
}
