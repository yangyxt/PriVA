//! `ensure_utr_reporting(directory, release)`: reproduce only the approved
//! VEP 116 UTRAnnotator reporting profile (all already-calculated UTR effects,
//! reported deterministically; the five biological methods are unchanged).
//! Pinned upstream and installed bytes, the original and the patch kept as
//! backups, and a pending receipt that makes an interrupted run resumable.
use super::{
    RECEIPT, Result, digest, difflib, is, load_json, load_receipt, local_patch_matches, object,
    parent, read_bytes, setdefault, sha256_hex, source, text, value_error, write_receipt,
};
use crate::file_lock::{ExclusiveLock, Wait};
use crate::installer_config_paths::files::{
    exists, is_link, metadata, mkdir, os_error, rename, temporary, unlink_missing,
};
use crate::installer_prerequisites::paths::join;
use crate::json_text::{Json, py_equal};
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::os::unix::fs::PermissionsExt;
use std::path::Path;

const PROFILE: &str = "utrannotator-116-all-effects-20261024-v1";
const UPSTREAM_SHA: &str = "a002fe5a654845627ce2d1f08b4f475db41a051e24002f862f4fc4f4c18d60b0";
const INSTALLED_SHA: &str = "0e4814fb71c9cc1af48a70830ae33434ba620e517f94ed8436979510484e9138";
const MODULE: &str = "UTRAnnotator.pm";

fn outcome(status: &str) -> Json {
    object(vec![("status", text(status)), ("profile", text(PROFILE))])
}

/// The recovery receipt holds exactly `receipt_before` and `receipt_after`.
fn transaction_shape(value: &Json) -> bool {
    match value {
        Json::Object(items) => {
            let mut keys: Vec<&str> = items.iter().map(|(name, _)| name.as_str()).collect();
            keys.sort_unstable();
            keys == ["receipt_after", "receipt_before"]
        }
        _ => false,
    }
}

fn remove(path: &str) -> Result<()> {
    fs::remove_file(path).map_err(|cause| os_error(cause, Path::new(path), None))
}

pub(super) fn ensure(directory: &str, release: &str) -> Result<Json> {
    if release != "116" {
        return Ok(object(vec![
            ("status", text("not_applicable")),
            ("release", text(release)),
        ]));
    }
    let target = join(directory, MODULE);
    let receipt_path = join(directory, RECEIPT);
    let pending = join(directory, ".priva-utr-reporting.pending.json");
    let lock = join(directory, ".priva-utr-reporting.lock");
    for path in [&target, &receipt_path, &pending, &lock] {
        if is_link(Path::new(path))? {
            return Err(value_error(
                "Refusing symlinked UTR reporting installation files",
            ));
        }
    }
    // A valid completed installation is a read-only no-op, including an
    // existing recorded patch. Do not create a new lock or backup there.
    let prior = load_receipt(directory);
    let recorded = prior
        .get("local_patches")
        .filter(|patches| matches!(patches, Json::Object(_)))
        .and_then(|patches| patches.get(MODULE))
        .filter(|patch| matches!(patch, Json::Object(_)));
    if is(prior.get("release"), release)
        && is(prior.get("source"), &source(release))
        && matches!(prior.get("plugins"), Some(Json::Object(_)))
        && is(
            prior.get("plugins").and_then(|plugins| plugins.get(MODULE)),
            UPSTREAM_SHA,
        )
        && recorded.is_some_and(|patch| {
            is(patch.get("profile"), PROFILE) && is(patch.get("installed_sha256"), INSTALLED_SHA)
        })
        && !exists(Path::new(&pending))?
        && local_patch_matches(directory, &prior, MODULE)
    {
        return Ok(outcome("reused"));
    }
    let guard = OpenOptions::new()
        .append(true)
        .create(true)
        .open(&lock)
        .map_err(|cause| os_error(cause, Path::new(&lock), None))?;
    let guard = ExclusiveLock::acquire(guard, Wait::Blocking)
        .map_err(|cause| os_error(cause, Path::new(&lock), None))?;
    let result = locked(directory, release, &target, &receipt_path, &pending);
    drop(guard);
    result
}

fn locked(
    directory: &str,
    release: &str,
    target: &str,
    receipt_path: &str,
    pending: &str,
) -> Result<Json> {
    let receipt = load_receipt(directory);
    if !is(receipt.get("release"), release)
        || !is(receipt.get("source"), &source(release))
        || !matches!(receipt.get("plugins"), Some(Json::Object(_)))
        || !is(
            receipt.get("plugins").and_then(|plugins| plugins.get(MODULE)),
            UPSTREAM_SHA,
        )
    {
        return Err(value_error(
            "UTR reporting fix requires the pinned official VEP116 module; files were preserved",
        ));
    }
    let empty = Json::Object(Vec::new());
    let patches = receipt.get("local_patches").unwrap_or(&empty);
    let recorded = patches.get(MODULE);
    if !matches!(patches, Json::Object(_))
        || recorded.is_some_and(|patch| !matches!(patch, Json::Object(_)))
    {
        return Err(value_error(
            "Malformed UTR reporting patch receipt; files were preserved",
        ));
    }
    if let Some(recorded @ Json::Object(fields)) = recorded
        && !fields.is_empty()
    {
        if !is(recorded.get("profile"), PROFILE)
            || !is(recorded.get("installed_sha256"), INSTALLED_SHA)
            || !local_patch_matches(directory, &receipt, MODULE)
        {
            return Err(value_error(
                "Existing UTR reporting patch differs from the approved profile; files were preserved",
            ));
        }
        if exists(Path::new(pending))? {
            let recovery = load_json(pending)?;
            let before = recovery.get("receipt_before");
            if !transaction_shape(&recovery)
                || !recovery
                    .get("receipt_after")
                    .is_some_and(|after| py_equal(after, &receipt))
                || !matches!(before, Some(Json::Object(_)))
            {
                return Err(value_error(
                    "Conflicting UTR reporting recovery receipt; files were preserved",
                ));
            }
            let before = before.unwrap_or(&empty);
            let old_patches = before.get("local_patches").unwrap_or(&empty);
            if !matches!(old_patches, Json::Object(_)) || old_patches.get(MODULE).is_some() {
                return Err(value_error(
                    "Conflicting UTR reporting recovery origin; files were preserved",
                ));
            }
            let mut expected = before.clone();
            setdefault(&mut expected, "local_patches").set(MODULE, recorded.clone());
            if !py_equal(&expected, &receipt) {
                return Err(value_error(
                    "Conflicting UTR reporting recovery origin; files were preserved",
                ));
            }
            remove(pending)?;
        }
        return Ok(outcome("reused"));
    }
    let before_sha = digest(receipt_path)?;
    let observed_sha = digest(target)?;
    if observed_sha != UPSTREAM_SHA && observed_sha != INSTALLED_SHA {
        return Err(value_error(
            "UTRAnnotator differs from the pinned original; files were preserved",
        ));
    }
    let recovery = if exists(Path::new(pending))? {
        Some(load_json(pending)?)
    } else {
        None
    };
    if recovery.as_ref().is_some_and(|value| !transaction_shape(value)) {
        return Err(value_error(
            "Malformed UTR reporting recovery receipt; files were preserved",
        ));
    }
    if observed_sha == INSTALLED_SHA && recovery.is_none() {
        return Err(value_error(
            "Unrecorded UTR reporting customization; files were preserved",
        ));
    }
    // Bundled executable source, not a downloaded annotation table. Its hash
    // proves identity to the already approved plugin profile.
    let module = crate::protein_domain_mapping::data_directory()
        .join("vep")
        .join("UTRAnnotator116_all_effects.pm");
    let candidate = read_bytes(&module.display().to_string())?;
    if sha256_hex(&candidate) != INSTALLED_SHA {
        return Err(value_error(
            "Bundled UTR reporting implementation does not reproduce the approved profile",
        ));
    }
    let folder = join(directory, ".priva-local-patches/UTRAnnotator116_all_effects");
    if is_link(Path::new(&folder))? || is_link(Path::new(&parent(&folder)))? {
        return Err(value_error(
            "Refusing a symlinked UTR patch backup directory",
        ));
    }
    let original = join(&folder, "UTRAnnotator.original.pm");
    let patch_file = join(&folder, "UTRAnnotator_all_effects.patch");
    let upstream = if observed_sha == UPSTREAM_SHA {
        read_bytes(target)?
    } else {
        read_bytes(&original)?
    };
    if sha256_hex(&upstream) != UPSTREAM_SHA {
        return Err(value_error(
            "UTR reporting recovery original differs; files were preserved",
        ));
    }
    let decode = |bytes: &[u8]| {
        crate::python_io::decode_utf8(bytes)
            .map(str::to_owned)
            .map_err(|message| super::error("UnicodeDecodeError", message))
    };
    let patch = difflib::unified_diff(
        &decode(&upstream)?,
        &decode(&candidate)?,
        "UTRAnnotator.original.pm",
        MODULE,
    )
    .into_bytes();
    let metadata_record = object(vec![
        ("profile", text(PROFILE)),
        ("upstream_sha256", text(UPSTREAM_SHA)),
        ("installed_sha256", text(INSTALLED_SHA)),
        ("original_file", text(&original)),
        ("patch_file", text(&patch_file)),
        ("patch_sha256", text(&sha256_hex(&patch))),
    ]);
    let mut updated = receipt.clone();
    setdefault(&mut updated, "local_patches").set(MODULE, metadata_record);
    let transaction = object(vec![
        ("receipt_before", receipt.clone()),
        ("receipt_after", updated.clone()),
    ]);
    if recovery
        .as_ref()
        .is_some_and(|value| !py_equal(value, &transaction))
    {
        return Err(value_error(
            "UTR reporting recovery context changed; files were preserved",
        ));
    }
    mkdir(Path::new(&folder))?;
    for (path, contents) in [(&original, &upstream), (&patch_file, &patch)] {
        let at = Path::new(path);
        if is_link(at)? || (exists(at)? && read_bytes(path)? != *contents) {
            return Err(value_error(format!(
                "Conflicting UTR reporting backup; file was preserved: {path}"
            )));
        }
        if !exists(at)? {
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(path)
                .and_then(|mut file| file.write_all(contents))
                .map_err(|cause| os_error(cause, at, None))?;
        }
    }
    // Record the exact recoverable transaction before either live write.
    write_receipt(pending, &transaction)?;
    if digest(receipt_path)? != before_sha || digest(target)? != observed_sha {
        return Err(value_error(
            "UTR reporting files changed during preparation; live files were preserved",
        ));
    }
    if observed_sha == UPSTREAM_SHA {
        let (mut file, staged) = temporary(Path::new(directory), ".priva-utr-")?;
        let written = file
            .write_all(&candidate)
            .map_err(|cause| os_error(cause, &staged, None));
        drop(file);
        written?;
        let published = metadata(Path::new(target))
            .and_then(|current| {
                fs::set_permissions(
                    &staged,
                    fs::Permissions::from_mode(current.permissions().mode() & 0o777),
                )
                .map_err(|cause| os_error(cause, &staged, None))
            })
            .and_then(|()| rename(&staged, Path::new(target)));
        let removed = unlink_missing(&staged);
        published.and(removed)?;
    }
    write_receipt(receipt_path, &updated)?;
    if !local_patch_matches(directory, &load_receipt(directory), MODULE) {
        return Err(value_error(
            "UTR reporting publication is incomplete; recovery receipt was retained",
        ));
    }
    remove(pending)?;
    Ok(outcome("installed"))
}
