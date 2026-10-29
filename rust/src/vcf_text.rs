//! Plain-text VCF input and output for converted PriVA scripts.
//!
//! Records are handled as tab-separated text, so columns a script does not
//! change are copied byte for byte. Input may be plain or gzip/BGZF
//! compressed; output ending in `.gz` or `.bgz` is written as BGZF so that
//! tabix and bcftools can index and read it.
use crate::cli::Result;
use flate2::read::MultiGzDecoder;
use flate2::write::DeflateEncoder;
use flate2::{Compression, Crc};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

/// Uncompressed bytes per BGZF block, as written by htslib.
const BLOCK_SIZE: usize = 0xff00;
const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

pub(crate) fn open_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 2];
    let read = file.read(&mut magic)?;
    let file = File::open(path)?;
    if read == 2 && magic == [0x1f, 0x8b] {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

pub(crate) fn create_writer(path: &Path) -> Result<Box<dyn TextSink>> {
    let file = BufWriter::new(File::create(path)?);
    let name = path.to_string_lossy();
    if name.ends_with(".gz") || name.ends_with(".bgz") {
        Ok(Box::new(BgzfWriter {
            inner: file,
            buffer: Vec::with_capacity(BLOCK_SIZE),
        }))
    } else {
        Ok(Box::new(PlainWriter { inner: file }))
    }
}

/// A line sink that must be finished explicitly, so a BGZF end-of-file block
/// is written only after every record succeeded.
pub(crate) trait TextSink {
    fn write_line(&mut self, line: &[u8]) -> io::Result<()>;
    fn finish(self: Box<Self>) -> io::Result<()>;
}

struct PlainWriter {
    inner: BufWriter<File>,
}

impl TextSink for PlainWriter {
    fn write_line(&mut self, line: &[u8]) -> io::Result<()> {
        self.inner.write_all(line)
    }

    fn finish(mut self: Box<Self>) -> io::Result<()> {
        self.inner.flush()
    }
}

struct BgzfWriter {
    inner: BufWriter<File>,
    buffer: Vec<u8>,
}

impl BgzfWriter {
    fn write_block(&mut self, data: &[u8]) -> io::Result<()> {
        let mut encoder = DeflateEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(data)?;
        let compressed = encoder.finish()?;
        let mut crc = Crc::new();
        crc.update(data);
        let block_size = 18 + compressed.len() + 8;
        if block_size > 0x10000 {
            return Err(io::Error::other("BGZF block exceeds 64 KiB"));
        }
        let total = (block_size - 1) as u16;
        self.inner.write_all(&[
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00,
        ])?;
        self.inner.write_all(&total.to_le_bytes())?;
        self.inner.write_all(&compressed)?;
        self.inner.write_all(&crc.sum().to_le_bytes())?;
        self.inner.write_all(&(data.len() as u32).to_le_bytes())
    }
}

impl TextSink for BgzfWriter {
    fn write_line(&mut self, mut line: &[u8]) -> io::Result<()> {
        while !line.is_empty() {
            let room = BLOCK_SIZE - self.buffer.len();
            let take = room.min(line.len());
            self.buffer.extend_from_slice(&line[..take]);
            line = &line[take..];
            if self.buffer.len() == BLOCK_SIZE {
                let block = std::mem::take(&mut self.buffer);
                self.write_block(&block)?;
                self.buffer = block;
                self.buffer.clear();
            }
        }
        Ok(())
    }

    fn finish(mut self: Box<Self>) -> io::Result<()> {
        if !self.buffer.is_empty() {
            let block = std::mem::take(&mut self.buffer);
            self.write_block(&block)?;
        }
        self.inner.write_all(&EOF_BLOCK)?;
        self.inner.flush()
    }
}

/// Keys present in a VCF INFO column ("." means none).
pub(crate) fn info_value<'a>(info: &'a str, key: &str) -> Option<Option<&'a str>> {
    if info == "." {
        return None;
    }
    info.split(';').find_map(|item| match item.split_once('=') {
        Some((name, value)) if name == key => Some(Some(value)),
        None if item == key => Some(None),
        _ => None,
    })
}

/// The canonical PASS filter line htslib adds to every header it writes.
const PASS_LINE: &str = "##FILTER=<ID=PASS,Description=\"All filters passed\">";

/// Header lines as pysam (htslib) writes them: the PASS filter line comes
/// second, right after `##fileformat`, and any PASS line elsewhere is dropped.
pub(crate) fn pysam_header(lines: &[String]) -> Vec<String> {
    let mut out = Vec::with_capacity(lines.len() + 1);
    for (index, line) in lines.iter().enumerate() {
        if line.starts_with("##FILTER=<ID=PASS,") {
            continue;
        }
        out.push(line.clone());
        if index == 0 {
            out.push(PASS_LINE.to_owned());
        }
    }
    out
}

/// INFO column after pysam's END synchronisation on write: htslib takes the
/// record length from INFO END when END is not before POS, else from REF.
/// END is removed when that length equals REF's and no allele is symbolic;
/// otherwise END is set to the last covered position (updated in place, or
/// appended at the end of INFO).
pub(crate) fn pysam_sync_end(pos: i64, reference: &str, alts: &[&str], info: &str) -> String {
    let symbolic = alts
        .iter()
        .any(|alt| alt.starts_with('<') && alt.ends_with('>'));
    let ref_len = reference.len() as i64;
    let end = match info_value(info, "END") {
        Some(Some(value)) => value.parse::<i64>().ok().filter(|end| *end >= pos),
        _ => None,
    };
    let length = end.map_or(ref_len, |end| end - pos + 1);
    let items: Vec<&str> = if info == "." {
        Vec::new()
    } else {
        info.split(';').collect()
    };
    let is_end = |item: &&str| *item == "END" || item.starts_with("END=");
    let mut out: Vec<String> = Vec::with_capacity(items.len() + 1);
    if length == ref_len && !symbolic {
        out.extend(
            items
                .iter()
                .filter(|item| !is_end(item))
                .map(|item| item.to_string()),
        );
    } else {
        let value = format!("END={}", pos + length - 1);
        let mut placed = false;
        for item in &items {
            if is_end(item) {
                out.push(value.clone());
                placed = true;
            } else {
                out.push(item.to_string());
            }
        }
        if !placed {
            out.push(value);
        }
    }
    if out.is_empty() {
        ".".to_owned()
    } else {
        out.join(";")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sync_end_matches_pysam_on_write() {
        // Expected values: pysam (priva_acmg environment on HGSC) writing the same records.
        assert_eq!(
            pysam_sync_end(10, "N", &["<DEL>"], "END=500;SVTYPE=DEL"),
            "END=500;SVTYPE=DEL"
        );
        assert_eq!(
            pysam_sync_end(20, "ACGT", &["A"], "END=23;SVTYPE=DEL"),
            "SVTYPE=DEL"
        );
        assert_eq!(
            pysam_sync_end(30, "A", &["C"], "SVTYPE=X;END=40"),
            "SVTYPE=X;END=40"
        );
        assert_eq!(
            pysam_sync_end(50, "ACG", &["<INV>"], "SVTYPE=INV"),
            "SVTYPE=INV;END=52"
        );
        assert_eq!(pysam_sync_end(60, "A", &["<*>"], "."), "END=60");
        assert_eq!(
            pysam_sync_end(70, "A", &["G", "<DUP>"], "SVTYPE=DUP"),
            "SVTYPE=DUP;END=70"
        );
        assert_eq!(pysam_sync_end(80, "A", &["G"], "."), ".");
        assert_eq!(pysam_sync_end(90, "AT", &["A"], "END=85"), ".");
    }

    #[test]
    fn header_puts_pass_second() {
        let lines: Vec<String> = [
            "##fileformat=VCFv4.2",
            "##FILTER=<ID=q10,Description=\"low\">",
            "##FILTER=<ID=PASS,Description=\"All filters passed\">",
            "#CHROM",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        assert_eq!(
            pysam_header(&lines),
            [
                "##fileformat=VCFv4.2",
                PASS_LINE,
                "##FILTER=<ID=q10,Description=\"low\">",
                "#CHROM"
            ]
        );
    }
}
