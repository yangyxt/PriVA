//! Read-only committed annex provenance; never unlock a recorded export.
use super::*;
use std::ffi::OsStr;
use std::fs::{File, OpenOptions};
use std::io::Read;
use std::os::unix::ffi::OsStrExt;
use std::os::unix::fs::{MetadataExt, OpenOptionsExt};
const REFUSAL: &str = "Refusing a linked installed-package export destination without stable committed git-annex content";
#[derive(PartialEq, Debug)]
struct Identity {
    dev: u64,
    ino: u64,
    mode: u32,
    size: u64,
    mtime: (i64, i64),
    ctime: (i64, i64),
    nlink: u64,
}
fn identity(m: fs::Metadata) -> Identity {
    Identity {
        dev: m.dev(),
        ino: m.ino(),
        mode: m.mode(),
        size: m.size(),
        mtime: (m.mtime(), m.mtime_nsec()),
        ctime: (m.ctime(), m.ctime_nsec()),
        nlink: m.nlink(),
    }
}
fn stat(path: &Path) -> Result<Identity> {
    Ok(identity(
        fs::symlink_metadata(path).map_err(|e| os_error(e, path, None))?,
    ))
}
fn fstat(file: &File) -> Result<Identity> {
    Ok(identity(file.metadata()?))
}
fn git(directory: &Path, args: &[&str]) -> Result<Vec<u8>> {
    process::git(directory, &args.iter().map(OsStr::new).collect::<Vec<_>>())
}
fn decoded(bytes: Vec<u8>) -> Result<String> {
    Ok(crate::python_io::decode_utf8(&bytes)
        .map_err(|m| error("UnicodeDecodeError", m))?
        .trim_matches(crate::python_text::is_python_space)
        .into())
}
fn owner_valid(owner: &str) -> bool {
    let value = owner.replace("urn:", "").replace("uuid:", "");
    let value = value.trim_matches(['{', '}']).replace('-', "");
    value.len() == 32 && value.bytes().all(|b| b.is_ascii_hexdigit())
}
#[derive(PartialEq, Debug)]
pub(super) struct Snapshot {
    source: PathBuf,
    before: Identity,
    literal: PathBuf,
    target: PathBuf,
    target_before: Identity,
    repo: PathBuf,
    common: PathBuf,
    entry: Vec<u8>,
    owner: String,
    pub payload: Vec<u8>,
}
fn checked(output: &Path) -> Result<Snapshot> {
    let source = files::realpath(output.parent().unwrap(), true)?.join(output.file_name().unwrap());
    let before = stat(&source)?;
    require(
        before.mode & libc::S_IFMT == libc::S_IFLNK,
        "Export link changed",
    )?;
    let literal = fs::read_link(&source).map_err(|e| os_error(e, &source, None))?;
    let repo = files::realpath(
        &path_value(&decoded(git(
            source.parent().unwrap(),
            &["rev-parse", "--show-toplevel"],
        )?)?),
        true,
    )?;
    let common = path_value(&decoded(git(&repo, &["rev-parse", "--git-common-dir"])?)?);
    let common = files::realpath(
        &if common.is_absolute() {
            common
        } else {
            repo.join(common)
        },
        true,
    )?;
    let owner = decoded(git(&repo, &["config", "--local", "--get", "annex.uuid"])?)?;
    require(owner_valid(&owner), "Invalid annex UUID")?;
    let objects = common.join("annex/objects");
    let target = crate::python_path::tempfile_parent(&source.parent().unwrap().join(&literal))?;
    require(
        !files::is_link(&objects)? && files::is_dir(&objects)? && target.starts_with(&objects),
        "Foreign export target",
    )?;
    let strict = crate::python_io::resolve_strict(&target).map_err(|(e, p)| {
        if e.raw_os_error() == Some(libc::ELOOP) {
            files::loop_error(&p)
        } else {
            os_error(e, &p, None)
        }
    })?;
    require(strict == target, "Foreign export target")?;
    let key = target.file_name().unwrap_or_default().to_string_lossy();
    let (algorithm, body) = key
        .split_once("-s")
        .ok_or_else(|| invalid("Unsupported annex content key"))?;
    let (size, digest) = body
        .split_once("--")
        .ok_or_else(|| invalid("Unsupported annex content key"))?;
    let (digest, _suffix) = digest
        .split_once('.')
        .map_or((digest, None), |(a, b)| (a, Some(b)));
    require(
        matches!(algorithm, "MD5E" | "SHA256E")
            && !size.is_empty()
            && size.bytes().all(|b| b.is_ascii_digit())
            && !digest.is_empty()
            && digest
                .bytes()
                .all(|b| b.is_ascii_digit() || matches!(b, b'a'..=b'f'))
            && target.parent().and_then(Path::file_name) == target.file_name(),
        "Unsupported annex content key",
    )?;
    let relative = source.strip_prefix(&repo)?;
    let tree = || {
        process::git(
            &repo,
            &[
                OsStr::new("ls-tree"),
                OsStr::new("-z"),
                OsStr::new("HEAD"),
                OsStr::new("--"),
                relative.as_os_str(),
            ],
        )
    };
    let entry = tree()?;
    let parts = entry
        .strip_suffix(&[0])
        .unwrap_or(&entry)
        .split(|b| *b == b'\t')
        .collect::<Vec<_>>();
    require(parts.len() == 2, "Invalid committed export pointer")?;
    let meta = parts[0]
        .split(|b| b.is_ascii_whitespace())
        .filter(|v| !v.is_empty())
        .collect::<Vec<_>>();
    require(
        meta.len() == 3
            && meta[0] == b"120000"
            && meta[1] == b"blob"
            && parts[1] == relative.as_os_str().as_bytes(),
        "Export pointer is not the committed literal link",
    )?;
    require(meta[2].is_ascii(), "Non-ASCII committed object ID")?;
    let blob = process::git(
        &repo,
        &[
            OsStr::new("cat-file"),
            OsStr::new("blob"),
            OsStr::from_bytes(meta[2]),
        ],
    )?;
    require(
        blob == literal.as_os_str().as_bytes(),
        "Export pointer is not the committed literal link",
    )?;
    let target_before = stat(&target)?;
    require(
        target_before.mode & libc::S_IFMT == libc::S_IFREG && target_before.size <= 1024 * 1024,
        "Invalid installed-package object",
    )?;
    let mut file = OpenOptions::new()
        .read(true)
        .custom_flags(libc::O_NOFOLLOW)
        .open(&target)
        .map_err(|e| os_error(e, &target, None))?;
    require(
        fstat(&file)? == target_before,
        "Export object changed before read",
    )?;
    let mut payload = vec![];
    (&mut file)
        .take(1024 * 1024 + 1)
        .read_to_end(&mut payload)?;
    require(
        fstat(&file)? == target_before,
        "Export object changed during read",
    )?;
    drop(file);
    let hash = if algorithm == "MD5E" {
        crate::md5::md5_hex(&payload)
    } else {
        crate::sha256::sha256_hex(&payload)
    };
    let current_entry = tree()?;
    let current_owner = decoded(git(&repo, &["config", "--local", "--get", "annex.uuid"])?)?;
    require(
        size.parse::<u64>().ok() == Some(payload.len() as u64)
            && hash == digest
            && current_entry == entry
            && current_owner == owner
            && stat(&source)? == before
            && fs::read_link(&source)? == literal
            && stat(&target)? == target_before,
        "Export link, key, object or committed pointer changed",
    )?;
    Ok(Snapshot {
        source,
        before,
        literal,
        target,
        target_before,
        repo,
        common,
        entry,
        owner,
        payload,
    })
}
pub(super) fn snapshot(output: &Path) -> Result<Snapshot> {
    checked(output).map_err(|cause| {
        if cause.is::<files::PathLoop>() {
            cause
        } else {
            error("ValueError", REFUSAL)
        }
    })
}
