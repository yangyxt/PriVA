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
pub(crate) fn write_cache(path: &Path, schema: Schema, archive: &[u8]) -> Result<()> {
    let name = path
        .file_name()
        .ok_or_else(|| invalid(format!("{}: not a file path", path.display())))?
        .to_string_lossy()
        .into_owned();
    let temporary = path.with_file_name(format!(".{name}.{}.tmp", std::process::id()));
    let outcome = (|| -> Result<()> {
        let mut file = std::fs::File::create(&temporary)?;
        file.write_all(&header(schema, archive.len()))?;
        file.write_all(archive)?;
        file.sync_all()?;
        std::fs::rename(&temporary, path)?;
        Ok(())
    })();
    if outcome.is_err() {
        let _ = std::fs::remove_file(&temporary);
    }
    outcome
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
        // SAFETY: cache files are immutable once renamed into place (see `write_cache`); the
        // whole archive is validated below before any reference is handed out.
        let map = unsafe { memmap2::Mmap::map(&file)? };
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
