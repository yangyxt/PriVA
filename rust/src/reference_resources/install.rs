//! `install`: each resource under its own lock, missing managed defaults
//! downloaded and published by link (never replacing a concurrent writer's
//! file), indexes repaired with `samtools faidx`, receipts written for managed
//! defaults, and finally the configuration rewritten when every resource is
//! ready.
use super::files::{
    Algorithm, digest, file_size, identity, identity_json, make_directories, publish_new,
    read_text, replace, unlink_missing_ok, write_text,
};
use super::plan::{Kind, Resource, checked_path};
use super::verify::{receipt_path, source_digest, valid_index, valid_index_at, validate_chain};
use super::{Result, Source, dumps_sorted, error, require, source};
use crate::file_lock::{ExclusiveLock, Wait};
use crate::installer_prerequisites::paths::{exists, is_file, is_symlink, parent};
use crate::json_text::Json;
use crate::python_path::{metadata, os_error, repr};
use crate::python_subprocess::signal_repr;
use std::ffi::CString;
use std::fs::{self, File, OpenOptions};
use std::os::unix::fs::PermissionsExt;
use std::os::unix::process::ExitStatusExt;
use std::path::Path;
use std::process::Command;

/// `shutil.which("samtools") is not None`
pub(super) fn samtools_on_path() -> bool {
    let path = std::env::var_os("PATH").unwrap_or_else(|| "/bin:/usr/bin".into());
    if path.is_empty() {
        return false;
    }
    std::env::split_paths(&path).any(|folder| {
        let candidate = folder.join("samtools");
        let executable = CString::new(candidate.as_os_str().as_encoded_bytes())
            .is_ok_and(|text| unsafe { libc::access(text.as_ptr(), libc::X_OK) } == 0);
        candidate.exists() && executable && !candidate.is_dir()
    })
}

/// `staging(path)`: a staging, lock or receipt path that is a regular file
/// when it exists, and never a symbolic link.
fn staging(path: &str) -> Result<String> {
    let path = checked_path(path)?;
    require(
        !is_symlink(&path)? && (!exists(&path)? || is_file(&path)?),
        format!("Reference staging path is not a regular file: {path}"),
    )?;
    Ok(path)
}

/// `download(source, destination)`: the verified pinned archive, reusing a
/// complete one and resuming an exact partial.
fn download(source: &Source, destination: &str) -> Result<String> {
    let destination = staging(destination)?;
    if exists(&destination)? {
        require(
            file_size(&destination)? == source.bytes
                && digest(&destination, Algorithm::Md5)? == source.md5,
            format!(
                "Existing UCSC source archive checksum mismatch: {destination}; existing bytes \
                 retained"
            ),
        )?;
        return Ok(destination);
    }
    let partial = staging(&format!("{destination}.partial"))?;
    let offset = if exists(&partial)? {
        file_size(&partial)?
    } else {
        0
    };
    require(
        offset <= source.bytes,
        format!(
            "UCSC partial is larger than the pinned source: {partial}; move it aside explicitly"
        ),
    )?;
    if offset < source.bytes {
        super::http::fetch(source.url, source.bytes, offset, &partial)?;
    }
    require(
        file_size(&partial)? == source.bytes && digest(&partial, Algorithm::Md5)? == source.md5,
        format!(
            "UCSC source archive checksum/size mismatch: {partial}; partial retained; move it \
             aside explicitly before restarting"
        ),
    )?;
    publish_new(&partial, &destination)?;
    Ok(destination)
}

/// `subprocess.run(["samtools", "faidx", "--fai-idx", index, path], check=True)`
fn samtools_faidx(index: &str, path: &str) -> Result<()> {
    let arguments = ["faidx", "--fai-idx", index, path];
    let status = Command::new("samtools")
        .args(arguments)
        .status()
        .map_err(|failure| os_error(failure, Path::new("samtools"), None))?;
    if status.success() {
        return Ok(());
    }
    let shown = std::iter::once("samtools")
        .chain(arguments)
        .map(repr)
        .collect::<Vec<_>>()
        .join(", ");
    let outcome = match (status.code(), status.signal()) {
        (Some(code), _) => format!("returned non-zero exit status {code}"),
        (None, Some(signal)) => format!("died with {}", signal_repr(signal)),
        (None, None) => format!("ended with {status}"),
    };
    Err(error(
        "CalledProcessError",
        format!("Command '[{shown}]' {outcome}."),
    ))
}

/// `index_fasta(path, assembly)`: repairs a missing or wrong index only; the
/// FASTA bytes are unchanged.
fn index_fasta(path: &str, assembly: &str) -> Result<()> {
    if valid_index(path, assembly)? {
        return Ok(());
    }
    let index = checked_path(&format!("{path}.fai"))?;
    let partial = staging(&format!("{index}.partial"))?;
    unlink_missing_ok(&partial)?;
    let before = identity(path)?;
    samtools_faidx(&partial, path)?;
    require(
        identity(path)? == before,
        format!("Reference changed while indexing: {path}; staged index retained"),
    )?;
    require(
        valid_index_at(&partial, assembly)?,
        format!("samtools index does not identify required {assembly} chromosome 1: {path}"),
    )?;
    replace(&partial, &index)
}

/// `install_resource(resource)`
pub(super) fn install_resource(resource: &Resource) -> Result<()> {
    let path = resource.path.as_str();
    let source = source(resource.assembly);
    make_directories(&parent(path))?;
    let lock_path = staging(&format!("{path}.install.lock"))?;
    let mut options = OpenOptions::new();
    let opened = options.append(true).create(true).open(&lock_path);
    let file = opened.map_err(|failure| os_error(failure, Path::new(&lock_path), None))?;
    let acquired = ExclusiveLock::acquire(file, Wait::Blocking);
    let lock = acquired.map_err(|failure| os_error(failure, Path::new(&lock_path), None))?;
    // Recheck after waiting for a concurrent writer.
    super::verify::preflight(std::slice::from_ref(resource))?;
    let output_sha256 = if exists(path)? {
        if resource.managed {
            Some(super::verify::verify_managed(resource)?)
        } else {
            None
        }
    } else {
        let archive = download(source, &format!("{path}.gz"))?;
        let partial = staging(&format!("{path}.partial"))?;
        let created = File::create(&partial);
        let mut output = created.map_err(|failure| os_error(failure, Path::new(&partial), None))?;
        let sha256 = source_digest(&archive, Some((&mut output, partial.as_str())))?;
        drop(output);
        if resource.kind == Kind::Fasta {
            unlink_missing_ok(&staging(&format!("{partial}.fai"))?)?;
            index_fasta(&partial, resource.assembly)?;
        } else {
            validate_chain(&partial)?;
        }
        publish_new(&partial, path)?;
        if resource.kind == Kind::Fasta {
            publish_new(&format!("{partial}.fai"), &format!("{path}.fai"))?;
        }
        Some(sha256)
    };
    if resource.kind == Kind::Fasta {
        index_fasta(path, resource.assembly)?;
    }
    if resource.managed {
        let receipt = receipt_path(path)?;
        let partial_receipt = staging(&format!("{receipt}.partial"))?;
        let metadata = Json::Object(vec![
            ("source".to_owned(), source.json()),
            ("output".to_owned(), Json::Str(path.to_owned())),
            (
                "sha256".to_owned(),
                output_sha256.map_or(Json::Null, Json::Str),
            ),
            ("identity".to_owned(), identity_json(&identity(path)?)),
        ]);
        let serialized = dumps_sorted(&metadata)? + "\n";
        if !exists(&receipt)? || read_text(&receipt)? != serialized {
            write_text(&partial_receipt, &serialized)?;
            replace(&partial_receipt, &receipt)?;
        }
    }
    eprintln!(
        "Validated {} {}: {path}",
        resource.assembly,
        resource.kind.name()
    );
    drop(lock);
    Ok(())
}

/// `write_config_values(path, original, updated)`: saved only when the
/// configuration is still the text that was planned, keeping its mode.
pub(super) fn write_config_values(path: &str, original: &str, updated: &str) -> Result<()> {
    require(
        read_text(path)? == original,
        format!(
            "Configuration changed during reference preparation: {path}; existing configuration \
             retained"
        ),
    )?;
    if updated != original {
        let partial = staging(&format!("{path}.reference.partial"))?;
        write_text(&partial, updated)?;
        let mode = metadata(Path::new(path))?.permissions().mode() & 0o7777;
        fs::set_permissions(&partial, fs::Permissions::from_mode(mode))
            .map_err(|failure| os_error(failure, Path::new(&partial), None))?;
        replace(&partial, path)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder_error::RuntimeError;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT: AtomicU64 = AtomicU64::new(0);

    fn folder() -> String {
        let path = std::env::temp_dir().join(format!(
            "priva-reference-install-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&path).unwrap();
        crate::installer_prerequisites::paths::abspath(&path.display().to_string()).unwrap()
    }

    fn message(failure: Box<dyn std::error::Error>) -> String {
        failure.downcast_ref::<RuntimeError>().unwrap().message()
    }

    #[test]
    fn staging_paths_must_be_regular_files() {
        let root = folder();
        let link = format!("{root}/link");
        std::os::unix::fs::symlink(format!("{root}/absent"), &link).unwrap();
        let failure = staging(&link).unwrap_err();
        let expected = format!("Reference staging path is not a regular file: {link}");
        assert_eq!(message(failure), expected);
        assert!(staging(&root).is_err());
        assert_eq!(
            staging(&format!("{root}/new")).unwrap(),
            format!("{root}/new")
        );
    }

    #[test]
    fn a_wrong_existing_archive_or_oversized_partial_is_never_fetched_over() {
        let root = folder();
        let chain = source("hg19_to_hg38");
        let archive = format!("{root}/chain.gz");
        fs::write(&archive, "wrong").unwrap();
        let failure = download(chain, &archive).unwrap_err();
        assert!(message(failure).starts_with("Existing UCSC source archive checksum mismatch"));
        fs::remove_file(&archive).unwrap();
        let partial = File::create(format!("{archive}.partial")).unwrap();
        partial.set_len(chain.bytes + 1).unwrap();
        let failure = download(chain, &archive).unwrap_err();
        assert!(message(failure).starts_with("UCSC partial is larger than the pinned source"));
        partial.set_len(chain.bytes).unwrap();
        let failure = download(chain, &archive).unwrap_err();
        assert!(message(failure).starts_with("UCSC source archive checksum/size mismatch"));
        assert!(!Path::new(&archive).exists());
    }

    #[test]
    fn configuration_rewrites_keep_the_mode_and_refuse_concurrent_changes() {
        let root = folder();
        let path = format!("{root}/config.yaml");
        fs::write(&path, "a: 1\n").unwrap();
        fs::set_permissions(&path, fs::Permissions::from_mode(0o640)).unwrap();
        write_config_values(&path, "a: 1\n", "a: 2\n").unwrap();
        assert_eq!(fs::read_to_string(&path).unwrap(), "a: 2\n");
        let mode = fs::metadata(&path).unwrap().permissions().mode() & 0o7777;
        assert_eq!(mode, 0o640);
        let failure = write_config_values(&path, "a: 1\n", "a: 3\n").unwrap_err();
        assert!(message(failure).starts_with("Configuration changed during reference preparation"));
        assert_eq!(fs::read_to_string(&path).unwrap(), "a: 2\n");
    }
}
