//! pandas' output-compression adapter. Codecs and archive parsing are provided
//! by mature libraries; this module owns pandas' filename and chunk semantics.
//! In particular, ZIP appends duplicate members, while compressed TAR fails
//! when pandas reopens it in *uncompressed* append mode for the second chunk.
use super::table;
use crate::cli::{Result, invalid};
use crate::pandas_table::Table;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum Format {
    Bzip2,
    Xz,
    Zstd,
    Zip,
    Tar,
}

impl Format {
    // pandas infers compression case-insensitively, with TAR compound suffixes
    // before their individual codec suffixes. Plain/gzip keep their old writer.
    pub(super) fn infer(lowercase: &str) -> Option<Self> {
        if [".tar", ".tar.gz", ".tar.bz2", ".tar.xz"]
            .iter()
            .any(|suffix| lowercase.ends_with(suffix))
        {
            Some(Self::Tar)
        } else if lowercase.ends_with(".bz2") {
            Some(Self::Bzip2)
        } else if lowercase.ends_with(".xz") {
            Some(Self::Xz)
        } else if lowercase.ends_with(".zip") {
            Some(Self::Zip)
        } else if lowercase.ends_with(".zst") {
            Some(Self::Zstd)
        } else {
            None
        }
    }
}

pub(super) fn publish(path: &Path, data: &Table, first: bool, format: Format) -> Result<()> {
    match format {
        Format::Zip => return zip(path, data, first),
        Format::Tar => return tar(path, data, first),
        _ => (),
    }
    let file = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(first)
        .append(!first)
        .open(path)?;
    let out = BufWriter::new(file);
    // Each to_csv call creates a complete stream, even for an empty later
    // chunk. Finish on serialization errors too, as Python's context manager
    // does; any close error takes precedence over the earlier write error.
    match format {
        Format::Bzip2 => {
            let mut writer = bzip2::write::BzEncoder::new(out, bzip2::Compression::best());
            let result = table::write(&mut writer, data, first);
            writer.finish()?.flush()?;
            result
        }
        Format::Xz => {
            let mut writer = xz2::write::XzEncoder::new(out, 6);
            let result = table::write(&mut writer, data, first);
            writer.finish()?.flush()?;
            result
        }
        Format::Zstd => {
            let mut writer = zstd::stream::write::Encoder::new(out, 3)?;
            let result = table::write(&mut writer, data, first);
            writer.finish()?.flush()?;
            result
        }
        _ => unreachable!(),
    }
}

fn basename(path: &Path) -> String {
    path.file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .into_owned()
}

// pathlib treats a dot-only basename (e.g. ".zip") as having no suffix.
fn without_suffix(name: &str, suffix: &str) -> String {
    name.strip_suffix(suffix)
        .filter(|stem| !stem.is_empty())
        .unwrap_or(name)
        .to_owned()
}

fn zip_time() -> Result<zip::DateTime> {
    let seconds: libc::time_t = SystemTime::now()
        .duration_since(UNIX_EPOCH)?
        .as_secs()
        .try_into()?;
    let mut tm = std::mem::MaybeUninit::<libc::tm>::uninit();
    // localtime_r writes the result into caller-owned storage. Python's
    // ZipInfo uses local wall time, with ZIP's two-second resolution.
    if unsafe { libc::localtime_r(&seconds, tm.as_mut_ptr()) }.is_null() {
        return Err(invalid("ValueError: ZIP local time is unavailable"));
    }
    let tm = unsafe { tm.assume_init() };
    zip::DateTime::from_date_and_time(
        (tm.tm_year + 1900).try_into()?,
        (tm.tm_mon + 1).try_into()?,
        tm.tm_mday.try_into()?,
        tm.tm_hour.try_into()?,
        tm.tm_min.try_into()?,
        tm.tm_sec.try_into()?,
    )
    .map_err(|_| invalid("ValueError: ZIP timestamp outside the supported range").into())
}

fn zip(path: &Path, data: &Table, first: bool) -> Result<()> {
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(first)
        .open(path)?;
    if !first {
        // Validate the existing directory on open, before serializing CSV.
        // Defer the writer itself: its Drop rewrites that directory even when
        // no member was added, while pandas leaves an empty append untouched.
        zip::ZipArchive::new(&mut file)?;
    }
    // pandas opens the archive first, then buffers one CSV chunk and adds a
    // member only when that buffer is nonempty (including the first header).
    let mut bytes = Vec::new();
    let result = table::write(&mut bytes, data, first);
    if !first && bytes.is_empty() {
        return result;
    }
    let mut writer = if first {
        zip::ZipWriter::new(file)
    } else {
        zip::ZipWriter::new_append(file)?
    };
    if !bytes.is_empty() {
        let name = without_suffix(&basename(path), ".zip");
        let options = zip::write::FileOptions::default()
            .compression_method(zip::CompressionMethod::Deflated)
            .unix_permissions(0o600)
            .last_modified_time(zip_time()?)
            // Python reserves ZIP64 before writing when size * 1.05 exceeds
            // its signed 32-bit ZIP64_LIMIT; no narrowing of large members.
            .large_file(bytes.len() as f64 * 1.05 > i32::MAX as f64);
        writer.start_file(name, options)?;
        writer.write_all(&bytes)?;
    }
    let mut file = writer.finish()?;
    // zip 0.6's writer retains duplicate names, as required by pandas. Its
    // portable defaults add S_IFREG and creator-version 46. Python writestr
    // uses permissions-only 0600 and the version needed by the member. Adapt
    // only these two central-header metadata fields; the crate owns all
    // compression, CRC, offsets, ZIP64 and append bookkeeping.
    let headers = {
        let mut archive = zip::ZipArchive::new(&mut file)?;
        (0..archive.len())
            .map(|i| archive.by_index_raw(i).map(|f| f.central_header_start()))
            .collect::<std::result::Result<Vec<_>, _>>()?
    };
    for header in headers {
        file.seek(SeekFrom::Start(header + 6))?;
        let mut needed = [0; 2];
        file.read_exact(&mut needed)?;
        file.seek(SeekFrom::Start(header + 4))?;
        file.write_all(&[needed[0], 3])?;
        file.seek(SeekFrom::Start(header + 38))?;
        file.write_all(&(0o600u32 << 16).to_le_bytes())?;
    }
    file.flush()?;
    result
}

fn tar_header(name: &str, size: u64, kind: tar::EntryType, mode: u32) -> Result<tar::Header> {
    let mut header = tar::Header::new_ustar();
    // TarInfo's defaults are fixed metadata, not the source file's stat.
    header.set_mode(mode);
    header.set_uid(0);
    header.set_gid(0);
    header.set_mtime(0);
    header.set_size(size);
    header.set_entry_type(kind);
    for (dest, character) in header.as_mut_bytes()[..100].iter_mut().zip(name.chars()) {
        *dest = if character.is_ascii() {
            character as u8
        } else {
            b'?'
        };
    }
    header.set_cksum();
    // Python spells the checksum as six octal digits, NUL, space.
    let checksum = format!("{:06o}\0 ", header.cksum()?);
    header.as_mut_bytes()[148..156].copy_from_slice(checksum.as_bytes());
    Ok(header)
}

fn pax_record(key: &str, value: &str) -> Vec<u8> {
    let body = format!(" {key}={value}\n");
    let mut size = body.len() + 1;
    loop {
        let next = body.len() + size.to_string().len();
        if next == size {
            return format!("{size}{body}").into_bytes();
        }
        size = next;
    }
}

fn tar_chunk(out: &mut impl Write, name: &str, bytes: &[u8], offset: u64) -> Result<()> {
    let mut builder = tar::Builder::new(&mut *out);
    let mut end = offset;
    if !bytes.is_empty() {
        let mut pax = Vec::new();
        if !name.is_ascii() || name.len() > 100 {
            pax.extend(pax_record("path", name));
        }
        let size = bytes.len() as u64;
        let size_in_header = if size >= 8u64.pow(11) {
            pax.extend(pax_record("size", &size.to_string()));
            0
        } else {
            size
        };
        if !pax.is_empty() {
            let header = tar_header(
                "././@PaxHeader",
                pax.len() as u64,
                tar::EntryType::XHeader,
                0,
            )?;
            builder.append(&header, pax.as_slice())?;
            end += 512 + (pax.len() as u64).div_ceil(512) * 512;
        }
        let header = tar_header(name, size_in_header, tar::EntryType::Regular, 0o644)?;
        builder.append(&header, bytes)?;
        end += 512 + size.div_ceil(512) * 512;
    }
    builder.finish()?;
    drop(builder);
    end += 1024;
    let padding = (10240 - end % 10240) % 10240;
    out.write_all(&[0; 10240][..padding as usize])?;
    out.flush()?;
    Ok(())
}

fn tar_append_offset(file: &mut File) -> Result<u64> {
    let mut archive = tar::Archive::new(file);
    let mut offset = 0;
    for entry in archive.entries_with_seek()? {
        let entry = entry.map_err(|e| invalid(format!("ReadError: {e}")))?;
        offset = entry.raw_file_position() + entry.size().div_ceil(512) * 512;
    }
    Ok(offset)
}

fn tar(path: &Path, data: &Table, first: bool) -> Result<()> {
    let mut file = OpenOptions::new()
        .read(!first)
        .write(true)
        .create(true)
        .truncate(first)
        .open(path)?;
    let offset = if first {
        0
    } else {
        tar_append_offset(&mut file)?
    };
    file.seek(SeekFrom::Start(offset))?;
    let name = without_suffix(&basename(path), ".tar");
    let mut bytes = Vec::new();
    let result = table::write(&mut bytes, data, first);
    let out = BufWriter::new(file);
    // pandas' TAR-specific extension check is case-sensitive, unlike its
    // initial format inference. Compressed member names retain the full
    // basename. Append never selects a codec, even for an empty next chunk.
    match if first {
        path.extension().and_then(|v| v.to_str())
    } else {
        None
    } {
        Some("gz") => {
            let mut gzip = crate::gzip_writer::GzipWriter::new(
                out,
                crate::gzip_writer::GzipOptions {
                    level: 9,
                    mtime: SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs() as u32,
                    filename: basename(path)
                        .strip_suffix(".gz")
                        .unwrap()
                        .chars()
                        .map(|c| u8::try_from(u32::from(c)).ok())
                        .collect(),
                },
            )?;
            tar_chunk(&mut gzip, &name, &bytes, offset)?;
            gzip.finish()?;
        }
        Some("bz2") => {
            let mut writer = bzip2::write::BzEncoder::new(out, bzip2::Compression::best());
            tar_chunk(&mut writer, &name, &bytes, offset)?;
            writer.finish()?.flush()?;
        }
        Some("xz") => {
            let mut writer = xz2::write::XzEncoder::new(out, 6);
            tar_chunk(&mut writer, &name, &bytes, offset)?;
            writer.finish()?.flush()?;
        }
        _ => tar_chunk(&mut { out }, &name, &bytes, offset)?,
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compound_suffix_precedes_the_individual_codec() {
        for name in ["out.tar.gz", "out.tar.bz2", "out.tar.xz", "out.tar"] {
            assert_eq!(Format::infer(name), Some(Format::Tar));
        }
        assert_eq!(Format::infer("out.bz2"), Some(Format::Bzip2));
        assert_eq!(Format::infer("out.xz"), Some(Format::Xz));
        assert_eq!(Format::infer("out.zst"), Some(Format::Zstd));
        assert_eq!(Format::infer("out.zip"), Some(Format::Zip));
        assert_eq!(Format::infer("out.gz"), None);
        assert_eq!(Format::infer("out.tsv"), None);
    }

    #[test]
    fn archive_name_keeps_pandas_case_and_dotfile_rules() {
        assert_eq!(without_suffix("a.zip", ".zip"), "a");
        assert_eq!(without_suffix(".zip", ".zip"), ".zip");
        assert_eq!(without_suffix("a.ZIP", ".zip"), "a.ZIP");
        assert_eq!(without_suffix("a.tar.gz", ".tar"), "a.tar.gz");
    }

    #[test]
    fn pax_record_length_counts_utf8_bytes_and_decimal_length() {
        for value in ["é", "path", &"x".repeat(90), &"x".repeat(990)] {
            let record = pax_record("path", value);
            let text = std::str::from_utf8(&record).unwrap();
            let length: usize = text.split_once(' ').unwrap().0.parse().unwrap();
            assert_eq!(length, record.len());
            assert!(text.ends_with(&format!(" path={value}\n")));
        }
    }
}
