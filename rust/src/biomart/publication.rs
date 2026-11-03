//! Original BioMart provenance and receipt publication, without acquisition.
use super::publication_json::Record;
use super::*;
use crate::clinvar_contract::resolve_non_strict;
use crate::python_path::tempfile_parent;
use std::io::{self, BufWriter, Write};
use std::os::unix::fs::OpenOptionsExt;
use std::time::{SystemTime, UNIX_EPOCH};

pub(super) fn external_path(context: &Context, path: &Path, ftp: bool) -> Result<()> {
    let root = if ftp {
        resolve_non_strict(
            &context
                .helper_root
                .join("scripts/biomart_ftp_nm_to_enst.py"),
        )?
        .parent()
        .and_then(Path::parent)
        .ok_or_else(|| error("IndexError", "1"))?
        .to_path_buf()
    } else {
        context.source_root()?
    };
    require(
        !resolve_non_strict(path)?.starts_with(resolve_non_strict(&root)?),
        "BioMart acquisition workspace must be outside the checkout",
    )
}

fn utc_now() -> Result<String> {
    let duration = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|e| error("OSError", e.to_string()))?;
    let seconds = duration.as_secs() as libc::time_t;
    let mut tm = std::mem::MaybeUninit::<libc::tm>::uninit();
    // The mature platform calendar routine performs the UTC conversion; no
    // locale-dependent formatting or executable-location inference is used.
    if unsafe { libc::gmtime_r(&seconds, tm.as_mut_ptr()) }.is_null() {
        return Err(io_error(io::Error::last_os_error()));
    }
    let tm = unsafe { tm.assume_init() };
    let fraction = if duration.subsec_micros() == 0 {
        String::new()
    } else {
        format!(".{:06}", duration.subsec_micros())
    };
    Ok(format!(
        "{:04}-{:02}-{:02}T{:02}:{:02}:{:02}{fraction}+00:00",
        tm.tm_year + 1900,
        tm.tm_mon + 1,
        tm.tm_mday,
        tm.tm_hour,
        tm.tm_min,
        tm.tm_sec
    ))
}

pub(super) fn make_provenance(context: &Context, mapping: &Path, operation: &str) -> Result<Json> {
    let observed = verify::ftp_table_identity(mapping)?;
    let mut value = get(facts(), "ftp_specification").clone();
    let selected = source(&string(FTP))?;
    value.set("bytes", get(selected, "expected_bytes").clone());
    value.set("sha256", get(selected, "expected_sha256").clone());
    value.set("statistics", observed);
    value.set("operation", string(operation));
    value.set("registered_utc", string(utc_now()?));
    let producer = tempfile_parent(
        &context
            .helper_root
            .join("scripts/biomart_ftp_nm_to_enst.py"),
    )?;
    value.set(
        "registration_producer_path",
        string(producer.display().to_string()),
    );
    value.set("registration_producer_sha256", string(sha256(&producer)?));
    value.set(
        "slurm_job_id",
        std::env::var("SLURM_JOB_ID").map_or(Json::Null, string),
    );
    value.set("provenance_basis",string("Exact reviewed table snapshot and fixed FTP source/method pins, first extracted by Slurm1515184; registration does not imply a new download or scientific acceptance."));
    Ok(value)
}

fn io_error(cause: io::Error) -> Box<dyn std::error::Error> {
    // TextIOWrapper.write/flush and os.fsync exceptions have no filename.
    let category = match cause.kind() {
        io::ErrorKind::PermissionDenied => "PermissionError",
        io::ErrorKind::BrokenPipe => "BrokenPipeError",
        _ => "OSError",
    };
    let message = if let Some(n) = cause.raw_os_error() {
        let raw = cause.to_string();
        format!(
            "[Errno {n}] {}",
            raw.strip_suffix(&format!(" (os error {n})"))
                .unwrap_or(&raw)
        )
    } else {
        cause.to_string()
    };
    error(category, message)
}
fn write_record(file: File, record: &Record, sync: bool) -> Result<()> {
    let mut out = BufWriter::with_capacity(8192, file);
    let result = (|| {
        record.write(true, true, |s| {
            out.write_all(s.as_bytes()).map_err(io_error)
        })?;
        out.write_all(b"\n").map_err(io_error)?;
        out.flush().map_err(io_error)?;
        if sync {
            out.get_ref().sync_all().map_err(io_error)?;
        }
        Ok(())
    })();
    // The source context manager closes/flushed the stream on failures too;
    // retain every already-written byte and any close-time error.
    out.flush().map_err(io_error)?;
    result
}

pub(super) fn atomic(path: &Path, record: &Record, ftp: bool) -> Result<()> {
    let temporary = if ftp {
        let basename = name(path);
        if basename.is_empty() {
            return Err(error(
                "ValueError",
                format!(
                    "PosixPath({}) has an empty name",
                    repr(&path.display().to_string())
                ),
            ));
        }
        path.with_file_name(format!("{basename}.part"))
    } else {
        PathBuf::from(format!("{}.part", path.display()))
    };
    if ftp {
        require(
            !exists(path)? && !is_link(path)?,
            format!("Preserve existing receipt: {}", path.display()),
        )?;
    } else {
        require(
            !is_link(&temporary)?,
            format!(
                "Refusing to replace a symlink used as a staging file: {}",
                temporary.display()
            ),
        )?;
    }
    let file = if ftp {
        fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
    } else {
        File::create(&temporary)
    }
    .map_err(|e| os_error(e, &temporary, None))?;
    write_record(file, record, false)?;
    fs::rename(&temporary, path).map_err(|e| os_error(e, &temporary, Some(path)))?;
    Ok(())
}

fn random_name() -> Result<String> {
    const ALPHABET: &[u8] = b"abcdefghijklmnopqrstuvwxyz0123456789_";
    let mut result = String::with_capacity(8);
    while result.len() < 8 {
        let mut bytes = [0_u8; 16];
        let count = unsafe { libc::getrandom(bytes.as_mut_ptr().cast(), bytes.len(), 0) };
        if count < 0 {
            let cause = io::Error::last_os_error();
            if cause.kind() == io::ErrorKind::Interrupted {
                continue;
            }
            return Err(io_error(cause));
        }
        for &b in &bytes[..count as usize] {
            if b < 222 {
                result.push(ALPHABET[(b % 37) as usize] as char);
                if result.len() == 8 {
                    break;
                }
            }
        }
    }
    Ok(result)
}
fn temporary(workspace: &Path, basename: &str) -> Result<(File, PathBuf)> {
    let parent = tempfile_parent(workspace)?;
    // Match tempfile's eight-character name alphabet, exclusive create,
    // permission mode, path joining and retry scope. Random tokens are bound
    // explicitly in comparisons; their literal bytes are not deterministic.
    for _ in 0..238_328 {
        let path = parent.join(format!(
            "{basename}.acquisition.{}.json.part",
            random_name()?
        ));
        match fs::OpenOptions::new()
            .read(true)
            .write(true)
            .create_new(true)
            .mode(0o600)
            .custom_flags(libc::O_NOFOLLOW)
            .open(&path)
        {
            Ok(file) => return Ok((file, path)),
            Err(e) if e.kind() == io::ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(os_error(e, &path, None)),
        }
    }
    Err(error(
        "FileExistsError",
        "[Errno 17] No usable temporary file name found",
    ))
}
pub(super) fn runtime(
    context: &Context,
    record: &Record,
    workspace: &Path,
    basename: &str,
    ftp: bool,
) -> Result<PathBuf> {
    if !ftp {
        external_path(context, workspace, false)?;
    }
    external_path(context, workspace, true)?;
    let (file, partial) = temporary(workspace, basename)?;
    write_record(file, record, true)?;
    let destination = partial.with_extension("");
    fs::hard_link(&partial, &destination).map_err(|e| os_error(e, &partial, Some(&destination)))?;
    fs::remove_file(&partial).map_err(|e| os_error(e, &partial, None))?;
    Ok(destination)
}
