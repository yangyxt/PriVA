//! The container of the branch's Rust-native keyed caches (`rust/CACHE_FORMATS.md`): a 32-byte
//! header followed by an rkyv 0.8 archive, read through a memory map.
//!
//! Header, all integers little-endian:
//!
//! | bytes | field |
//! | --- | --- |
//! | 0..8 | magic `PRIVARKV` |
//! | 8..12 | container version, 1 |
//! | 12..16 | schema id (one per cache, listed in `CACHE_FORMATS.md`) |
//! | 16..20 | schema version |
//! | 20..24 | archive layout, 1 = rkyv 0.8, little-endian, aligned, 32-bit relative pointers |
//! | 24..32 | archive length in bytes (the file is exactly header plus archive) |
//!
//! A reader checks every header field against what it was built for, then validates the whole
//! archive once with rkyv's `bytecheck` (`rkyv::access`), and only then hands out references.
//! Files are written to a temporary name in the target folder and renamed into place, so a
//! reader never maps a half-written cache.
use crate::cli::{Result, invalid};
use rkyv::api::high::HighValidator;
use rkyv::bytecheck::CheckBytes;
use rkyv::{Archive, Portable, rancor};
use std::io::Write;
use std::marker::PhantomData;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

pub(crate) const MAGIC: &[u8; 8] = b"PRIVARKV";
pub(crate) const CONTAINER_VERSION: u32 = 1;
pub(crate) const LAYOUT_RKYV08_LE_ALIGNED_P32: u32 = 1;
pub(crate) const HEADER_LEN: usize = 32;

/// What a reader expects to find in a cache file.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Schema {
    pub(crate) name: &'static str,
    pub(crate) id: u32,
    pub(crate) version: u32,
}

pub(crate) fn header(schema: Schema, archive_len: usize) -> [u8; HEADER_LEN] {
    let mut out = [0_u8; HEADER_LEN];
    out[0..8].copy_from_slice(MAGIC);
    out[8..12].copy_from_slice(&CONTAINER_VERSION.to_le_bytes());
    out[12..16].copy_from_slice(&schema.id.to_le_bytes());
    out[16..20].copy_from_slice(&schema.version.to_le_bytes());
    out[20..24].copy_from_slice(&LAYOUT_RKYV08_LE_ALIGNED_P32.to_le_bytes());
    out[24..32].copy_from_slice(&(archive_len as u64).to_le_bytes());
    out
}

/// Checks the header of `bytes` (a whole cache file) against `schema`; returns the archive.
pub(crate) fn archive_bytes<'a>(bytes: &'a [u8], schema: Schema, path: &str) -> Result<&'a [u8]> {
    let refuse = |what: String| invalid(format!("{path}: {what}; not a {} cache", schema.name));
    if bytes.len() < HEADER_LEN || &bytes[0..8] != MAGIC {
        return Err(refuse(
            "no PRIVARKV header (a Python pickle must first be converted with priva cache-convert)"
                .to_owned(),
        ));
    }
    let word = |at: usize| u32::from_le_bytes(bytes[at..at + 4].try_into().unwrap());
    if word(8) != CONTAINER_VERSION {
        return Err(refuse(format!("container version {}", word(8))));
    }
    if word(12) != schema.id || word(16) != schema.version {
        return Err(refuse(format!(
            "schema {} version {} (expected {} version {})",
            word(12),
            word(16),
            schema.id,
            schema.version
        )));
    }
    if word(20) != LAYOUT_RKYV08_LE_ALIGNED_P32 {
        return Err(refuse(format!("archive layout {}", word(20))));
    }
    let length = u64::from_le_bytes(bytes[24..32].try_into().unwrap());
    if length != (bytes.len() - HEADER_LEN) as u64 {
        return Err(refuse(format!(
            "archive length {length} but {} bytes follow the header",
            bytes.len() - HEADER_LEN
        )));
    }
    Ok(&bytes[HEADER_LEN..])
}

/// Writes header and archive to `path` through a temporary file in the same folder.
#[cfg_attr(not(test), allow(dead_code))]
pub(crate) fn write_cache(path: &Path, schema: Schema, archive: &[u8]) -> Result<()> {
    write_cache_checked(path, schema, archive, |_| Ok(()))
}

/// Stages a complete cache, checks it, and publishes only after the check succeeds.
/// Converters use `validate` for a fresh mmap and complete comparison with their source.
/// A failed check preserves an existing destination byte for byte.
pub(crate) fn write_cache_checked(
    path: &Path,
    schema: Schema,
    archive: &[u8],
    validate: impl FnOnce(&Path) -> Result<()>,
) -> Result<()> {
    if archive.len() > i32::MAX as usize {
        return Err(invalid(
            "the archive exceeds the PRIVARKV 32-bit layout limit",
        ));
    }
    let name = path
        .file_name()
        .ok_or_else(|| invalid(format!("{}: not a file path", path.display())))?
        .to_string_lossy()
        .into_owned();
    // Exclusive creation prevents another writer's temporary file (or a stale symlink)
    // from being opened and truncated. Each successful caller owns its own name.
    static NEXT_TEMPORARY: AtomicU64 = AtomicU64::new(0);
    let mut selected = None;
    for _ in 0..128 {
        let serial = NEXT_TEMPORARY.fetch_add(1, Ordering::Relaxed);
        let temporary = path.with_file_name(format!(".{name}.{}.{serial}.tmp", std::process::id()));
        match std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)
        {
            Ok(file) => {
                selected = Some((temporary, file));
                break;
            }
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(error.into()),
        }
    }
    let (temporary, mut file) =
        selected.ok_or_else(|| invalid("cannot reserve a unique cache temporary file"))?;
    let outcome = (|| -> Result<()> {
        file.write_all(&header(schema, archive.len()))?;
        file.write_all(archive)?;
        file.sync_all()?;
        drop(file);
        validate(&temporary)?;
        std::fs::rename(&temporary, path)?;
        Ok(())
    })();
    if outcome.is_err() {
        let _ = std::fs::remove_file(&temporary);
    }
    outcome
}

/// Refuse an output which names the source, including symlink and hard-link aliases.
/// This is deliberately called before reading/converting the source cache.
pub(crate) fn ensure_distinct_paths(source: &Path, output: &Path) -> Result<()> {
    let source_metadata = std::fs::metadata(source)?;
    let source_path = std::fs::canonicalize(source)?;
    match std::fs::metadata(output) {
        Ok(output_metadata) => {
            #[cfg(unix)]
            {
                use std::os::unix::fs::MetadataExt;
                if source_metadata.dev() == output_metadata.dev()
                    && source_metadata.ino() == output_metadata.ino()
                {
                    return Err(invalid(
                        "the output must not be the pickle itself or its alias",
                    ));
                }
            }
            #[cfg(not(unix))]
            let _ = (source_metadata, output_metadata);
            if source_path == std::fs::canonicalize(output)? {
                return Err(invalid(
                    "the output must not be the pickle itself or its alias",
                ));
            }
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            let parent = output
                .parent()
                .filter(|p| !p.as_os_str().is_empty())
                .unwrap_or(Path::new("."));
            let name = output
                .file_name()
                .ok_or_else(|| invalid("the output needs a file name"))?;
            if source_path == std::fs::canonicalize(parent)?.join(name) {
                return Err(invalid(
                    "the output must not be the pickle itself or its alias",
                ));
            }
        }
        Err(error) => return Err(error.into()),
    }
    Ok(())
}

/// A validated cache file, memory-mapped; `root()` is the archived value.
pub(crate) struct MappedCache<T: Archive> {
    map: memmap2::Mmap,
    _value: PhantomData<T>,
}

impl<T> MappedCache<T>
where
    T: Archive,
    T::Archived: Portable + for<'a> CheckBytes<HighValidator<'a, rancor::Error>>,
{
    pub(crate) fn open(path: &str, schema: Schema) -> Result<MappedCache<T>> {
        let file = std::fs::File::open(path)
            .map_err(|error| invalid(format!("FileNotFoundError: {path}: {error}")))?;
        Self::open_file(&file, path, schema)
    }

    /// Versioned readers inspect and map the same open file, so an atomic path
    /// replacement cannot separate header selection from archive validation.
    pub(crate) fn open_file(
        file: &std::fs::File,
        path: &str,
        schema: Schema,
    ) -> Result<MappedCache<T>> {
        // SAFETY: cache files are immutable once renamed into place (see `write_cache`); the
        // whole archive is validated below before any reference is handed out.
        let map = unsafe { memmap2::Mmap::map(file)? };
        let archive = archive_bytes(&map, schema, path)?;
        rkyv::access::<T::Archived, rancor::Error>(archive).map_err(|error| {
            invalid(format!(
                "{path}: the {} archive is damaged: {error}",
                schema.name
            ))
        })?;
        Ok(MappedCache {
            map,
            _value: PhantomData,
        })
    }

    pub(crate) fn root(&self) -> &T::Archived {
        // SAFETY: `open` validated these exact bytes with `rkyv::access`, and the map is
        // read-only and lives as long as `self`.
        unsafe { rkyv::access_unchecked::<T::Archived>(&self.map[HEADER_LEN..]) }
    }
}

/// The position of `key` in a vector of archived texts sorted by their bytes.
pub(crate) fn find_sorted(keys: &[rkyv::string::ArchivedString], key: &str) -> Option<usize> {
    keys.binary_search_by(|probe| probe.as_str().cmp(key)).ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    const SCHEMA: Schema = Schema {
        name: "test",
        id: 99,
        version: 2,
    };

    #[test]
    fn failed_validation_preserves_existing_output() {
        let folder = std::env::temp_dir().join(format!("priva_publish_{}", std::process::id()));
        std::fs::create_dir_all(&folder).unwrap();
        let output = folder.join("checked.rkyv");
        std::fs::write(&output, b"previous cache").unwrap();
        let failed = write_cache_checked(&output, SCHEMA, b"candidate", |temporary| {
            assert_ne!(temporary, output);
            assert_eq!(
                archive_bytes(&std::fs::read(temporary)?, SCHEMA, "staged")?,
                b"candidate"
            );
            Err(invalid("a value differs"))
        });
        assert!(failed.unwrap_err().to_string().contains("a value differs"));
        assert_eq!(std::fs::read(&output).unwrap(), b"previous cache");
        assert_eq!(std::fs::read_dir(&folder).unwrap().count(), 1);
        write_cache_checked(&output, SCHEMA, b"candidate", |_| Ok(())).unwrap();
        assert_eq!(
            archive_bytes(&std::fs::read(&output).unwrap(), SCHEMA, "published").unwrap(),
            b"candidate"
        );
        std::fs::remove_file(output).unwrap();
        std::fs::remove_dir(folder).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn source_aliases_are_refused() {
        let folder = std::env::temp_dir().join(format!("priva_alias_{}", std::process::id()));
        std::fs::create_dir_all(&folder).unwrap();
        let source = folder.join("source.pkl");
        let hard = folder.join("hard.pkl");
        let link = folder.join("link.pkl");
        std::fs::write(&source, b"pickle").unwrap();
        std::fs::hard_link(&source, &hard).unwrap();
        std::os::unix::fs::symlink(&source, &link).unwrap();
        for path in [&source, &hard, &link] {
            assert!(ensure_distinct_paths(&source, path).is_err());
        }
        assert!(ensure_distinct_paths(&source, &folder.join("new.rkyv")).is_ok());
        for path in [&source, &hard, &link] {
            std::fs::remove_file(path).unwrap();
        }
        std::fs::remove_dir(folder).unwrap();
    }

    #[test]
    fn header_fields_are_checked() {
        let mut bytes = header(SCHEMA, 4).to_vec();
        bytes.extend_from_slice(&[1, 2, 3, 4]);
        assert_eq!(archive_bytes(&bytes, SCHEMA, "f").unwrap(), &[1, 2, 3, 4]);
        let other = Schema { id: 98, ..SCHEMA };
        assert!(
            archive_bytes(&bytes, other, "f")
                .unwrap_err()
                .to_string()
                .contains("schema 99")
        );
        let newer = Schema {
            version: 3,
            ..SCHEMA
        };
        assert!(archive_bytes(&bytes, newer, "f").is_err());
        assert!(
            archive_bytes(&bytes[..35], SCHEMA, "f")
                .unwrap_err()
                .to_string()
                .contains("archive length 4")
        );
        let mut pickle = vec![0x80, 0x04];
        pickle.extend_from_slice(&[0; 40]);
        assert!(
            archive_bytes(&pickle, SCHEMA, "f")
                .unwrap_err()
                .to_string()
                .contains("cache-convert")
        );
        let mut layout = bytes.clone();
        layout[20] = 2;
        assert!(archive_bytes(&layout, SCHEMA, "f").is_err());
    }
}
