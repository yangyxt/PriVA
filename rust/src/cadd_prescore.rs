//! The CADD prescored lookup on a compact copy of a prescored table, and the
//! builder of that copy.
//!
//! CADD's prescore step (`Snakefile`, checkpoint `prescore`) runs
//! `python src/scripts/extract_scored.py --header -p TABLE --found_out=FILE`
//! once per prescored table, with the variants on standard input. PriVA's
//! copy of that script (`scripts/cadd_scripts/extract_scored.py`) hands the
//! work to `priva cadd-prescore` when a compact copy of TABLE exists and
//! `PRIVA_CADD_PRESCORE_BIN` names this program; otherwise it runs its own
//! Python lookup on the full table.
//!
//! The compact copy keeps only the six columns PriVA reads from CADD's output
//! (`#Chrom`, `Pos`, `Ref`, `Alt`, `FeatureID`, `PHRED`), in the same row
//! order, in BGZF blocks of about 4 KB of text instead of 64 KB. Its index,
//! `TABLE.priva_compact.idx`, lists the first position and the file offset
//! of every block, so a lookup decompresses only the blocks that hold a
//! wanted position. A block never holds two chromosomes, and a row never
//! spans two blocks.
//!
//! `priva cadd-prescore` writes what the Python lookup writes, with one
//! planned difference: a found row has all of the table's columns, but only
//! the six kept columns carry values; the others are empty. CADD's `join`
//! rule sorts and concatenates these rows with the rows CADD scores itself,
//! so the joined table keeps its column layout. Everything else matches the
//! Python lookup: the stderr line `Opening TABLE...`, the table's header lines
//! with `--header`, a row is found for a variant only when its `Pos`, `Ref`
//! and `Alt` strings equal the variant's (its second, second-to-last and last
//! fields), found rows in table order for each variant in input order, every
//! other variant line unchanged on standard output, `#` lines passed through,
//! and `Encountered uncovered chromosome` for a variant whose chromosome is not
//! in the table or whose position is not a number.
//!
//! Lookups run on several threads (`PRIVA_CADD_PRESCORE_THREADS`, default the
//! available CPUs up to 8); the output does not depend on the thread count.
//!
//! Index layout (little-endian): MAGIC, full table size (u64), compact file
//! size (u64), target block text size (u32), number of full-table columns
//! (u32), number of kept columns (u32) and their full-table indices (u32
//! each), header text length (u32) and the header lines (each ending in a
//! newline), number of chromosomes (u32); then per chromosome in file order:
//! name length (u32), name, number of blocks (u64), first positions (u32
//! each), block file offsets (u64 each).
use crate::cli::{Result, invalid};
use flate2::{Compress, Compression, Crc, Decompress, FlushCompress, FlushDecompress, Status};
use std::collections::{BTreeMap, HashMap};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Condvar, Mutex};
use std::time::Instant;

const MAGIC: &[u8; 16] = b"PRIVA-CADD6IDX\x01\x00";
const KEPT: [&str; 6] = ["#Chrom", "Pos", "Ref", "Alt", "FeatureID", "PHRED"];
const COMPACT_SUFFIX: &str = ".priva_compact";
const INDEX_SUFFIX: &str = ".priva_compact.idx";
const DEFAULT_BLOCK: usize = 4096;
/// Text per BGZF block as htslib writes it; a compact block never exceeds it.
const MAX_BLOCK_TEXT: usize = 0xff00;
const BATCH_LINES: usize = 200_000;
const ITEM_POSITIONS: usize = 2048;
/// Positions the Python lookup sends through its merged reader (tabix range).
const MAX_POSITION: i128 = 1 << 29;
/// pysam refuses a region end above this; the Python lookup then reports an
/// uncovered chromosome.
const PYSAM_MAX_POS: i128 = (1 << 31) - 1;
const READ_AHEAD: u64 = 256 << 10;
const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];
const LOOKUP_USAGE: &str =
    "usage: priva cadd-prescore -p TABLE [-i INPUT] [--found_out FILE] [--header]";
const BUILD_USAGE: &str = "usage: priva cadd-prescore-compact {build,check} TABLE [--threads N] [--block-size BYTES] [--level 1-9]";

fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
}

fn find(text: &[u8], byte: u8) -> Option<usize> {
    text.iter().position(|&value| value == byte)
}

// ---------------------------------------------------------------- BGZF blocks

/// Length of the BGZF block whose first 18 bytes are `header`.
fn block_length(header: &[u8]) -> Option<usize> {
    let valid = header.len() >= 18
        && header[..4] == [0x1f, 0x8b, 0x08, 0x04]
        && header[10..16] == [0x06, 0x00, b'B', b'C', 0x02, 0x00];
    valid.then(|| u16::from_le_bytes([header[16], header[17]]) as usize + 1)
}

/// Read the next whole block from a sequential reader; None at end of file.
fn read_block(reader: &mut impl Read, block: &mut Vec<u8>) -> Result<bool> {
    block.resize(18, 0);
    let mut filled = 0;
    while filled < 18 {
        let read = reader.read(&mut block[filled..])?;
        if read == 0 {
            if filled == 0 {
                return Ok(false);
            }
            return Err(invalid("the table ends inside a BGZF block header"));
        }
        filled += read;
    }
    let length = block_length(block).ok_or_else(|| invalid("not a BGZF block"))?;
    if length < 26 {
        return Err(invalid("BGZF block shorter than its header"));
    }
    block.resize(length, 0);
    reader.read_exact(&mut block[18..])?;
    Ok(true)
}

/// Decompress one whole block into `text`, checking its length and CRC.
fn inflate(block: &[u8], text: &mut Vec<u8>, inflater: &mut Decompress) -> Result<()> {
    let length = block.len();
    if length < 26 {
        return Err(invalid("truncated BGZF block"));
    }
    let size = u32::from_le_bytes(block[length - 4..].try_into().unwrap()) as usize;
    let expected_crc = u32::from_le_bytes(block[length - 8..length - 4].try_into().unwrap());
    text.clear();
    text.reserve(size.max(1));
    inflater.reset(false);
    let status = inflater
        .decompress_vec(&block[18..length - 8], text, FlushDecompress::Finish)
        .map_err(|error| invalid(format!("a BGZF block does not inflate: {error}")))?;
    if status != Status::StreamEnd || text.len() != size {
        return Err(invalid("a BGZF block inflates to a different length"));
    }
    let mut crc = Crc::new();
    crc.update(text);
    if crc.sum() != expected_crc {
        return Err(invalid("a BGZF block fails its CRC check"));
    }
    Ok(())
}

/// Append `text` to `out` as one BGZF block.
fn deflate(text: &[u8], compressor: &mut Compress, out: &mut Vec<u8>) -> Result<()> {
    compressor.reset();
    let mut data = Vec::with_capacity(text.len() + text.len() / 8 + 64);
    let status = compressor
        .compress_vec(text, &mut data, FlushCompress::Finish)
        .map_err(|error| invalid(format!("deflate failed: {error}")))?;
    let total = 18 + data.len() + 8;
    if status != Status::StreamEnd || total > 0x10000 {
        return Err(invalid("a compact block does not fit in one BGZF block"));
    }
    out.extend_from_slice(&EOF_BLOCK[..16]);
    out.extend_from_slice(&((total - 1) as u16).to_le_bytes());
    out.extend_from_slice(&data);
    let mut crc = Crc::new();
    crc.update(text);
    out.extend_from_slice(&crc.sum().to_le_bytes());
    out.extend_from_slice(&(text.len() as u32).to_le_bytes());
    Ok(())
}

/// Offset of the first BGZF block that starts in [start, end), or end.
/// A candidate counts only when another block or the end of file follows it.
fn first_block(file: &File, start: u64, end: u64, file_size: u64) -> Result<u64> {
    if start == 0 {
        return Ok(0);
    }
    let mut window_start = start;
    while window_start < end {
        let length = (1u64 << 20).min(file_size - window_start) as usize;
        let mut window = vec![0u8; length];
        file.read_exact_at(&mut window, window_start)?;
        let mut from = 0;
        while from + 16 <= window.len() {
            let Some(found) = window[from..]
                .windows(16)
                .position(|candidate| candidate == &EOF_BLOCK[..16])
            else {
                break;
            };
            let offset = window_start + (from + found) as u64;
            if offset >= end {
                return Ok(end);
            }
            let mut header = [0u8; 18];
            if offset + 18 <= file_size {
                file.read_exact_at(&mut header, offset)?;
                if let Some(length) = block_length(&header) {
                    let next = offset + length as u64;
                    if next == file_size {
                        return Ok(offset);
                    }
                    if next + 16 <= file_size {
                        let mut following = [0u8; 16];
                        file.read_exact_at(&mut following, next)?;
                        if following == EOF_BLOCK[..16] {
                            return Ok(offset);
                        }
                    }
                }
            }
            from += found + 1;
        }
        // Overlap windows so a signature that straddles two windows is seen.
        window_start += (length as u64).saturating_sub(15).max(1);
    }
    Ok(end)
}

/// Header lines (`#` lines at the start, without newlines) of a BGZF table.
fn table_header(table: &Path) -> Result<Vec<Vec<u8>>> {
    let mut reader = BufReader::with_capacity(1 << 20, File::open(table)?);
    let (mut block, mut text) = (Vec::new(), Vec::new());
    let mut inflater = Decompress::new(false);
    let (mut pending, mut lines) = (Vec::new(), Vec::new());
    while read_block(&mut reader, &mut block)? {
        inflate(&block, &mut text, &mut inflater)?;
        pending.extend_from_slice(&text);
        while let Some(end) = find(&pending, b'\n') {
            if pending.first() != Some(&b'#') {
                return Ok(lines);
            }
            lines.push(pending[..end].to_vec());
            pending.drain(..=end);
        }
        if !pending.is_empty() && pending[0] != b'#' {
            return Ok(lines);
        }
        if pending.len() > 1 << 24 {
            return Err(invalid("the table's header is not followed by rows"));
        }
    }
    if pending.first() == Some(&b'#') {
        lines.push(pending);
    }
    Ok(lines)
}

/// The Ref and Alt columns the Python lookup compares, from the header lines:
/// each line naming Ref sets the Ref column, and also the Alt column when it
/// names Alt (defaults 2 and 3).
fn ref_alt_columns(header: &[Vec<u8>]) -> (usize, usize) {
    let (mut reference, mut alternative) = (2, 3);
    for line in header {
        let fields: Vec<&[u8]> = line.split(|&byte| byte == b'\t').collect();
        if let Some(found) = fields.iter().position(|field| *field == b"Ref") {
            reference = found;
            if let Some(found) = fields.iter().position(|field| *field == b"Alt") {
                alternative = found;
            }
        }
    }
    (reference, alternative)
}

// ------------------------------------------------------------- compact index

struct CompactIndex {
    compact_size: u64,
    columns: usize,
    kept: Vec<usize>,
    header: Vec<Vec<u8>>,
    chromosomes: Vec<Vec<u8>>,
    by_name: HashMap<Vec<u8>, usize>,
    positions: Vec<Vec<u32>>,
    offsets: Vec<Vec<u64>>,
}

struct Bytes<'a> {
    data: &'a [u8],
    at: usize,
}

impl Bytes<'_> {
    fn take(&mut self, count: usize) -> Result<&[u8]> {
        let end = self
            .at
            .checked_add(count)
            .filter(|end| *end <= self.data.len())
            .ok_or_else(|| invalid("the compact index is truncated"))?;
        let slice = &self.data[self.at..end];
        self.at = end;
        Ok(slice)
    }
    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
}

/// True when TABLE's compact copy and index exist and belong to TABLE.
fn compact_matches(table: &Path) -> bool {
    let read = || -> Result<bool> {
        let mut head = [0u8; 32];
        File::open(with_suffix(table, INDEX_SUFFIX))?.read_exact(&mut head)?;
        let table_size = fs::metadata(table)?.len();
        let compact_size = fs::metadata(with_suffix(table, COMPACT_SUFFIX))?.len();
        Ok(&head[..16] == MAGIC
            && u64::from_le_bytes(head[16..24].try_into().unwrap()) == table_size
            && u64::from_le_bytes(head[24..32].try_into().unwrap()) == compact_size)
    };
    read().unwrap_or(false)
}

fn load_index(table: &Path) -> Result<CompactIndex> {
    let path = with_suffix(table, INDEX_SUFFIX);
    if !compact_matches(table) {
        return Err(invalid(format!(
            "{} is missing or belongs to a different table; build it with priva cadd-prescore-compact build",
            path.display()
        )));
    }
    let data = fs::read(&path)?;
    let mut bytes = Bytes { data: &data, at: 0 };
    bytes.take(16)?;
    bytes.u64()?;
    let compact_size = bytes.u64()?;
    bytes.u32()?;
    let columns = bytes.u32()? as usize;
    let kept_count = bytes.u32()? as usize;
    let kept = (0..kept_count)
        .map(|_| bytes.u32().map(|value| value as usize))
        .collect::<Result<Vec<_>>>()?;
    let header_length = bytes.u32()? as usize;
    let header_text = bytes.take(header_length)?;
    let header = header_text
        .split(|&byte| byte == b'\n')
        .filter(|line| !line.is_empty())
        .map(<[u8]>::to_vec)
        .collect();
    let count = bytes.u32()? as usize;
    let (mut chromosomes, mut by_name) = (Vec::new(), HashMap::new());
    let (mut positions, mut offsets) = (Vec::new(), Vec::new());
    for index in 0..count {
        let name_length = bytes.u32()? as usize;
        let name = bytes.take(name_length)?.to_vec();
        let blocks = bytes.u64()? as usize;
        let mut first = Vec::with_capacity(blocks);
        for _ in 0..blocks {
            first.push(bytes.u32()?);
        }
        let mut offset = Vec::with_capacity(blocks);
        for _ in 0..blocks {
            offset.push(bytes.u64()?);
        }
        by_name.insert(name.clone(), index);
        chromosomes.push(name);
        positions.push(first);
        offsets.push(offset);
    }
    if bytes.at != data.len() || kept.len() != KEPT.len() || kept.iter().any(|&k| k >= columns) {
        return Err(invalid("the compact index is malformed"));
    }
    Ok(CompactIndex {
        compact_size,
        columns,
        kept,
        header,
        chromosomes,
        by_name,
        positions,
        offsets,
    })
}

// ------------------------------------------------------------------- builder

struct Layout {
    columns: usize,
    kept: Vec<usize>,
}

fn layout(header: &[Vec<u8>]) -> Result<Layout> {
    let names = header
        .iter()
        .rev()
        .find(|line| !line.starts_with(b"##"))
        .ok_or_else(|| invalid("the table has no column header line"))?;
    let fields: Vec<&[u8]> = names.split(|&byte| byte == b'\t').collect();
    let kept = KEPT
        .iter()
        .map(|name| {
            fields
                .iter()
                .position(|field| *field == name.as_bytes())
                .ok_or_else(|| invalid(format!("the table has no {name} column")))
        })
        .collect::<Result<Vec<_>>>()?;
    if kept[0] != 0 || kept[1] != 1 || kept.windows(2).any(|pair| pair[0] >= pair[1]) {
        return Err(invalid(
            "the table's #Chrom and Pos must be its first two columns and the kept columns in order",
        ));
    }
    Ok(Layout {
        columns: fields.len(),
        kept,
    })
}

#[derive(Default)]
struct Part {
    data: Vec<u8>,
    /// (chromosome within this part, first position, offset within data)
    blocks: Vec<(usize, u32, u64)>,
    chromosomes: Vec<Vec<u8>>,
    last: Option<u32>,
    rows: u64,
}

struct Packer<'a> {
    layout: &'a Layout,
    block_size: usize,
    compressor: Compress,
    part: Part,
    block: Vec<u8>,
    block_first: Option<(usize, u32)>,
    row: Vec<u8>,
    tabs: Vec<usize>,
}

fn parse_position(text: &[u8]) -> Option<u32> {
    if text.is_empty() || text.len() > 10 || !text.iter().all(u8::is_ascii_digit) {
        return None;
    }
    std::str::from_utf8(text).ok()?.parse().ok()
}

impl Packer<'_> {
    fn flush_block(&mut self) -> Result<()> {
        if let Some((chromosome, position)) = self.block_first.take() {
            let offset = self.part.data.len() as u64;
            deflate(&self.block, &mut self.compressor, &mut self.part.data)?;
            self.part.blocks.push((chromosome, position, offset));
            self.block.clear();
        }
        Ok(())
    }

    fn row(&mut self, line: &[u8]) -> Result<()> {
        if line.first() == Some(&b'#') {
            return Ok(());
        }
        self.tabs.clear();
        self.tabs.extend(
            line.iter()
                .enumerate()
                .filter(|(_, b)| **b == b'\t')
                .map(|(i, _)| i),
        );
        if self.tabs.len() + 1 != self.layout.columns {
            return Err(invalid(format!(
                "a table row has {} columns, not {}",
                self.tabs.len() + 1,
                self.layout.columns
            )));
        }
        let bounds = |index: usize, tabs: &[usize]| -> (usize, usize) {
            let start = if index == 0 { 0 } else { tabs[index - 1] + 1 };
            (start, tabs.get(index).copied().unwrap_or(line.len()))
        };
        let mut kept = [(0usize, 0usize); KEPT.len()];
        for (k, &column) in self.layout.kept.iter().enumerate() {
            kept[k] = bounds(column, &self.tabs);
        }
        let chromosome = &line[kept[0].0..kept[0].1];
        let position_text = &line[kept[1].0..kept[1].1];
        let position = parse_position(position_text).ok_or_else(|| {
            invalid(format!(
                "position {:?} is not a whole number",
                String::from_utf8_lossy(position_text)
            ))
        })?;
        let same = self
            .part
            .chromosomes
            .last()
            .is_some_and(|last| last == chromosome);
        if same {
            if self.part.last.is_some_and(|last| position < last) {
                return Err(invalid(format!(
                    "the table is not sorted by position at {}:{position}",
                    String::from_utf8_lossy(chromosome)
                )));
            }
        } else {
            if self.part.chromosomes.iter().any(|seen| seen == chromosome) {
                return Err(invalid(format!(
                    "the table is not grouped by chromosome at {}",
                    String::from_utf8_lossy(chromosome)
                )));
            }
            self.flush_block()?;
            self.part.chromosomes.push(chromosome.to_vec());
        }
        let local = self.part.chromosomes.len() - 1;
        self.row.clear();
        for (k, &(start, end)) in kept.iter().enumerate() {
            if k > 0 {
                self.row.push(b'\t');
            }
            self.row.extend_from_slice(&line[start..end]);
        }
        self.row.push(b'\n');
        if self.row.len() > MAX_BLOCK_TEXT {
            return Err(invalid("a compact row is longer than one BGZF block"));
        }
        if !self.block.is_empty()
            && (self.block.len() + self.row.len() > self.block_size
                || self.block.len() + self.row.len() > MAX_BLOCK_TEXT)
        {
            self.flush_block()?;
        }
        if self.block_first.is_none() {
            self.block_first = Some((local, position));
        }
        let row = std::mem::take(&mut self.row);
        self.block.extend_from_slice(&row);
        self.row = row;
        self.part.last = Some(position);
        self.part.rows += 1;
        Ok(())
    }
}

/// Compact rows that start in the blocks of [start, end). A part drops the
/// text before its first newline; it owns the rows after it, and the row that
/// ends at the first newline after its range.
fn build_part(
    table: &Path,
    start: u64,
    end: u64,
    file_size: u64,
    layout: &Layout,
    block_size: usize,
    level: u32,
) -> Result<Part> {
    let file = File::open(table)?;
    let first = first_block(&file, start, end, file_size)?;
    let mut packer = Packer {
        layout,
        block_size,
        compressor: Compress::new(Compression::new(level), false),
        part: Part::default(),
        block: Vec::with_capacity(block_size + 1024),
        block_first: None,
        row: Vec::new(),
        tabs: Vec::new(),
    };
    if first >= end {
        return Ok(packer.part);
    }
    let mut reader = BufReader::with_capacity(4 << 20, file);
    reader.seek(SeekFrom::Start(first))?;
    let (mut block, mut text) = (Vec::new(), Vec::new());
    let mut inflater = Decompress::new(false);
    let mut offset = first;
    let mut drop_head = start != 0;
    let mut pending: Vec<u8> = Vec::new();
    loop {
        if !read_block(&mut reader, &mut block)? {
            if !drop_head && !pending.is_empty() {
                packer.row(&pending)?;
            }
            break;
        }
        let in_range = offset < end;
        offset += block.len() as u64;
        inflate(&block, &mut text, &mut inflater)?;
        let mut rest = &text[..];
        if drop_head {
            if !in_range {
                break;
            }
            match find(rest, b'\n') {
                Some(newline) => {
                    rest = &rest[newline + 1..];
                    drop_head = false;
                }
                None => continue,
            }
        }
        if !in_range {
            match find(rest, b'\n') {
                Some(newline) => {
                    pending.extend_from_slice(&rest[..newline]);
                    if !pending.is_empty() {
                        packer.row(&pending)?;
                    }
                    break;
                }
                None => {
                    pending.extend_from_slice(rest);
                    continue;
                }
            }
        }
        if !pending.is_empty() {
            match find(rest, b'\n') {
                Some(newline) => {
                    pending.extend_from_slice(&rest[..newline]);
                    let line = std::mem::take(&mut pending);
                    packer.row(&line)?;
                    pending = line;
                    pending.clear();
                    rest = &rest[newline + 1..];
                }
                None => {
                    pending.extend_from_slice(rest);
                    continue;
                }
            }
        }
        while let Some(newline) = find(rest, b'\n') {
            packer.row(&rest[..newline])?;
            rest = &rest[newline + 1..];
        }
        pending.extend_from_slice(rest);
    }
    packer.flush_block()?;
    Ok(packer.part)
}

struct Schedule {
    next: usize,
    written: usize,
    done: BTreeMap<usize, Part>,
    error: Option<String>,
}

#[derive(Default)]
struct Assembled {
    order: Vec<Vec<u8>>,
    positions: Vec<Vec<u32>>,
    offsets: Vec<Vec<u64>>,
    last: Option<u32>,
    rows: u64,
}

impl Assembled {
    fn add(&mut self, part: Part, base: u64) -> Result<()> {
        for (local, position, offset) in &part.blocks {
            let name = &part.chromosomes[*local];
            if self.order.last() != Some(name) {
                if self.order.contains(name) {
                    return Err(invalid(format!(
                        "the table is not grouped by chromosome at {}",
                        String::from_utf8_lossy(name)
                    )));
                }
                self.order.push(name.clone());
                self.positions.push(Vec::new());
                self.offsets.push(Vec::new());
                self.last = None;
            }
            if self.last.is_some_and(|last| *position < last) {
                return Err(invalid(format!(
                    "the table is not sorted by position at {}:{position}",
                    String::from_utf8_lossy(name)
                )));
            }
            self.last = Some(*position);
            self.positions.last_mut().unwrap().push(*position);
            self.offsets.last_mut().unwrap().push(base + offset);
        }
        if let (Some(last), Some(_)) = (part.last, part.blocks.last()) {
            self.last = Some(last);
        }
        self.rows += part.rows;
        Ok(())
    }
}

/// Write TABLE.priva_compact and its index; returns (rows, blocks, bytes).
fn build(
    table: &Path,
    threads: usize,
    block_size: usize,
    level: u32,
    parts: usize,
) -> Result<(u64, u64, u64)> {
    let file_size = fs::metadata(table)?.len();
    let header = table_header(table)?;
    let layout = layout(&header)?;
    let parts = parts.max(1);
    let threads = threads.max(1);
    let bounds: Vec<u64> = (0..=parts)
        .map(|part| (file_size as u128 * part as u128 / parts as u128) as u64)
        .collect();
    let directory = table.parent().unwrap_or(Path::new("."));
    let stem = table.file_name().unwrap().to_string_lossy();
    let pid = std::process::id();
    let compact_tmp = directory.join(format!(".{stem}.priva_compact.{pid}.tmp"));
    let index_tmp = directory.join(format!(".{stem}.priva_compact.idx.{pid}.tmp"));
    let outcome = (|| -> Result<(u64, u64, u64)> {
        let mut out = BufWriter::with_capacity(8 << 20, File::create(&compact_tmp)?);
        let mut header_text = Vec::new();
        for line in &header {
            header_text.extend_from_slice(line);
            header_text.push(b'\n');
        }
        let mut compressor = Compress::new(Compression::new(level), false);
        let mut written_bytes = 0u64;
        for chunk in header_text.chunks(MAX_BLOCK_TEXT) {
            let mut block = Vec::new();
            deflate(chunk, &mut compressor, &mut block)?;
            out.write_all(&block)?;
            written_bytes += block.len() as u64;
        }
        let window = threads * 3;
        let schedule = Mutex::new(Schedule {
            next: 0,
            written: 0,
            done: BTreeMap::new(),
            error: None,
        });
        let wake = Condvar::new();
        let mut assembled = Assembled::default();
        let started = Instant::now();
        std::thread::scope(|scope| {
            for _ in 0..threads {
                scope.spawn(|| {
                    loop {
                        let part_index = {
                            let mut state = schedule.lock().unwrap();
                            while state.error.is_none()
                                && state.next < parts
                                && state.next >= state.written + window
                            {
                                state = wake.wait(state).unwrap();
                            }
                            if state.error.is_some() || state.next >= parts {
                                return;
                            }
                            state.next += 1;
                            state.next - 1
                        };
                        let outcome = build_part(
                            table,
                            bounds[part_index],
                            bounds[part_index + 1],
                            file_size,
                            &layout,
                            block_size,
                            level,
                        );
                        let mut state = schedule.lock().unwrap();
                        match outcome {
                            Ok(part) => {
                                state.done.insert(part_index, part);
                            }
                            Err(error) => {
                                state.error.get_or_insert(error.to_string());
                            }
                        }
                        wake.notify_all();
                    }
                });
            }
            let mut reported = 0;
            loop {
                let next_part = {
                    let mut state = schedule.lock().unwrap();
                    loop {
                        if state.error.is_some() || state.written >= parts {
                            break None;
                        }
                        let want = state.written;
                        if let Some(part) = state.done.remove(&want) {
                            break Some(part);
                        }
                        state = wake.wait(state).unwrap();
                    }
                };
                let Some(part) = next_part else { break };
                let length = part.data.len() as u64;
                let appended = out
                    .write_all(&part.data)
                    .map_err(|error| error.to_string())
                    .and_then(|_| {
                        assembled
                            .add(part, written_bytes)
                            .map_err(|e| e.to_string())
                    });
                let mut state = schedule.lock().unwrap();
                match appended {
                    Ok(()) => {
                        written_bytes += length;
                        state.written += 1;
                    }
                    Err(error) => {
                        state.error.get_or_insert(error);
                    }
                }
                let written = state.written;
                wake.notify_all();
                drop(state);
                if written * 50 / parts > reported {
                    reported = written * 50 / parts;
                    eprintln!(
                        "{written}/{parts} parts, {} rows, {:.0} s",
                        assembled.rows,
                        started.elapsed().as_secs_f64()
                    );
                }
            }
        });
        if let Some(error) = schedule.into_inner().unwrap().error {
            return Err(invalid(error));
        }
        out.write_all(&EOF_BLOCK)?;
        let file = out
            .into_inner()
            .map_err(|error| invalid(error.to_string()))?;
        file.sync_all()?;
        drop(file);
        let compact_size = fs::metadata(&compact_tmp)?.len();
        if fs::metadata(table)?.len() != file_size {
            return Err(invalid("the table changed while it was read"));
        }
        let mut index = BufWriter::new(File::create(&index_tmp)?);
        index.write_all(MAGIC)?;
        index.write_all(&file_size.to_le_bytes())?;
        index.write_all(&compact_size.to_le_bytes())?;
        index.write_all(&(block_size as u32).to_le_bytes())?;
        index.write_all(&(layout.columns as u32).to_le_bytes())?;
        index.write_all(&(layout.kept.len() as u32).to_le_bytes())?;
        for &column in &layout.kept {
            index.write_all(&(column as u32).to_le_bytes())?;
        }
        index.write_all(&(header_text.len() as u32).to_le_bytes())?;
        index.write_all(&header_text)?;
        index.write_all(&(assembled.order.len() as u32).to_le_bytes())?;
        let mut blocks = 0u64;
        for (chromosome, name) in assembled.order.iter().enumerate() {
            index.write_all(&(name.len() as u32).to_le_bytes())?;
            index.write_all(name)?;
            let positions = &assembled.positions[chromosome];
            index.write_all(&(positions.len() as u64).to_le_bytes())?;
            for position in positions {
                index.write_all(&position.to_le_bytes())?;
            }
            for offset in &assembled.offsets[chromosome] {
                index.write_all(&offset.to_le_bytes())?;
            }
            blocks += positions.len() as u64;
        }
        index
            .into_inner()
            .map_err(|error| invalid(error.to_string()))?
            .sync_all()?;
        fs::rename(&compact_tmp, with_suffix(table, COMPACT_SUFFIX))?;
        fs::rename(&index_tmp, with_suffix(table, INDEX_SUFFIX))?;
        Ok((assembled.rows, blocks, compact_size))
    })();
    if outcome.is_err() {
        let _ = fs::remove_file(&compact_tmp);
        let _ = fs::remove_file(&index_tmp);
    }
    outcome
}

// -------------------------------------------------------------------- lookup

#[derive(Clone, Copy)]
struct Span {
    start: usize,
    end: usize,
}

struct Merged {
    chromosome: usize,
    position: u32,
    line: Vec<u8>,
    pos_field: Span,
    reference: Span,
    alternative: Span,
}

enum Entry {
    Header(Vec<u8>),
    /// The Python lookup prints "Encountered uncovered chromosome" for these.
    Uncovered(Vec<u8>),
    NotFound(Vec<u8>),
    Merged(Merged),
}

/// Python's int() on ASCII text: None where it raises ValueError.
fn python_int(text: &[u8]) -> Option<i128> {
    let is_space = |byte: &u8| {
        matches!(
            byte,
            b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c | 0x1c..=0x1f
        )
    };
    let start = text.iter().position(|b| !is_space(b))?;
    let end = text.iter().rposition(|b| !is_space(b))? + 1;
    let mut digits = &text[start..end];
    let negative = match digits.first() {
        Some(b'-') => {
            digits = &digits[1..];
            true
        }
        Some(b'+') => {
            digits = &digits[1..];
            false
        }
        _ => false,
    };
    if digits.is_empty()
        || !digits[0].is_ascii_digit()
        || !digits[digits.len() - 1].is_ascii_digit()
    {
        return None;
    }
    let mut value: i128 = 0;
    let mut previous_underscore = false;
    for &byte in digits {
        if byte == b'_' {
            if previous_underscore {
                return None;
            }
            previous_underscore = true;
            continue;
        }
        if !byte.is_ascii_digit() {
            return None;
        }
        previous_underscore = false;
        value = value
            .saturating_mul(10)
            .saturating_add(i128::from(byte - b'0'));
    }
    Some(if negative { -value } else { value })
}

fn field_spans(line: &[u8]) -> Vec<Span> {
    let mut spans = Vec::new();
    let mut start = 0;
    for (index, &byte) in line.iter().enumerate() {
        if byte == b'\t' {
            spans.push(Span { start, end: index });
            start = index + 1;
        }
    }
    spans.push(Span {
        start,
        end: line.len(),
    });
    spans
}

/// Classify one input line (without its line ending). Err is the Python
/// lookup's IndexError for a line without a second field.
fn classify(line: Vec<u8>, index: &CompactIndex) -> std::result::Result<Entry, Vec<u8>> {
    if line.first() == Some(&b'#') {
        return Ok(Entry::Header(line));
    }
    let spans = field_spans(&line);
    if spans.len() < 2 {
        return Err(line);
    }
    let field = |span: Span| &line[span.start..span.end];
    let chromosome = index.by_name.get(field(spans[0])).copied();
    let Some(value) = python_int(field(spans[1])) else {
        return Ok(Entry::Uncovered(line));
    };
    match chromosome {
        Some(chromosome) if (1..MAX_POSITION).contains(&value) => Ok(Entry::Merged(Merged {
            chromosome,
            position: value as u32,
            pos_field: spans[1],
            reference: spans[spans.len() - 2],
            alternative: spans[spans.len() - 1],
            line,
        })),
        None => Ok(Entry::Uncovered(line)),
        Some(_) if value < 1 || value > PYSAM_MAX_POS => Ok(Entry::Uncovered(line)),
        Some(_) => Ok(Entry::NotFound(line)),
    }
}

/// Reads compact rows block by block, forwards only.
struct Cursor<'a> {
    file: &'a File,
    file_size: u64,
    positions: &'a [u32],
    offsets: &'a [u64],
    cache: Vec<u8>,
    cache_start: u64,
    text: Vec<u8>,
    inflater: Decompress,
    loaded: Option<usize>,
    at: usize,
}

impl<'a> Cursor<'a> {
    fn new(file: &'a File, file_size: u64, positions: &'a [u32], offsets: &'a [u64]) -> Self {
        Cursor {
            file,
            file_size,
            positions,
            offsets,
            cache: Vec::new(),
            cache_start: 0,
            text: Vec::new(),
            inflater: Decompress::new(false),
            loaded: None,
            at: 0,
        }
    }

    fn bytes(&mut self, offset: u64, length: usize) -> Result<&[u8]> {
        let end = offset + length as u64;
        if offset < self.cache_start || end > self.cache_start + self.cache.len() as u64 {
            let want = (length as u64).max(READ_AHEAD).min(self.file_size - offset);
            if want < length as u64 {
                return Err(invalid("the compact table is truncated"));
            }
            self.cache.resize(want as usize, 0);
            self.file.read_exact_at(&mut self.cache, offset)?;
            self.cache_start = offset;
        }
        let from = (offset - self.cache_start) as usize;
        Ok(&self.cache[from..from + length])
    }

    fn load(&mut self, block: usize) -> Result<()> {
        let offset = self.offsets[block];
        let length = block_length(self.bytes(offset, 18)?)
            .ok_or_else(|| invalid("the compact index points outside a BGZF block"))?;
        self.bytes(offset, length)?;
        let from = (offset - self.cache_start) as usize;
        let mut text = std::mem::take(&mut self.text);
        inflate(
            &self.cache[from..from + length],
            &mut text,
            &mut self.inflater,
        )?;
        self.text = text;
        self.loaded = Some(block);
        self.at = 0;
        Ok(())
    }

    /// Call `emit` with every row at `position`, in table order. Positions
    /// must come in increasing order.
    fn rows_at(&mut self, position: u32, mut emit: impl FnMut(&[u8])) -> Result<()> {
        if self.positions.is_empty() {
            return Ok(());
        }
        let start = self
            .positions
            .partition_point(|&first| first < position)
            .saturating_sub(1);
        if self.loaded.is_none_or(|loaded| start > loaded) {
            self.load(start)?;
        }
        loop {
            if self.at >= self.text.len() {
                let next = self.loaded.unwrap() + 1;
                if next >= self.offsets.len() {
                    return Ok(());
                }
                self.load(next)?;
                continue;
            }
            let rest = &self.text[self.at..];
            let length = find(rest, b'\n').unwrap_or(rest.len());
            let row = &rest[..length];
            let tab = find(row, b'\t').ok_or_else(|| invalid("a compact row has no tab"))?;
            let after = &row[tab + 1..];
            let row_position = parse_position(&after[..find(after, b'\t').unwrap_or(after.len())])
                .ok_or_else(|| invalid("a compact row has no position"))?;
            if row_position > position {
                return Ok(());
            }
            if row_position == position {
                emit(row);
            }
            self.at += length + 1;
        }
    }
}

struct Matcher {
    /// Compact field indices of Pos, Ref and Alt.
    pos: usize,
    reference: usize,
    alternative: usize,
}

/// Look up one work item: sorted variant entries of one chromosome.
fn lookup_item(
    index: &CompactIndex,
    file: &File,
    matcher: &Matcher,
    batch: &[Entry],
    chromosome: usize,
    variants: &[usize],
) -> Result<Vec<(usize, Vec<u8>)>> {
    let mut cursor = Cursor::new(
        file,
        index.compact_size,
        &index.positions[chromosome],
        &index.offsets[chromosome],
    );
    let merged = |entry: usize| match &batch[entry] {
        Entry::Merged(merged) => merged,
        _ => unreachable!(),
    };
    let mut found = Vec::new();
    let mut at = 0;
    let mut fields: Vec<Span> = Vec::new();
    while at < variants.len() {
        let position = merged(variants[at]).position;
        let mut end = at;
        while end < variants.len() && merged(variants[end]).position == position {
            end += 1;
        }
        let group = &variants[at..end];
        cursor.rows_at(position, |row| {
            fields.clear();
            fields.extend(field_spans(row));
            if fields.len() <= matcher.pos.max(matcher.reference).max(matcher.alternative) {
                return;
            }
            let value = |span: Span| &row[span.start..span.end];
            for &entry in group {
                let variant = merged(entry);
                let text = |span: Span| &variant.line[span.start..span.end];
                if value(fields[matcher.reference]) == text(variant.reference)
                    && value(fields[matcher.alternative]) == text(variant.alternative)
                    && value(fields[matcher.pos]) == text(variant.pos_field)
                {
                    found.push((entry, row.to_vec()));
                }
            }
        })?;
        at = end;
    }
    Ok(found)
}

struct Expander {
    /// Tabs to write before each kept field, then after the last one.
    gaps: Vec<usize>,
}

impl Expander {
    fn new(columns: usize, kept: &[usize]) -> Self {
        let mut gaps = Vec::new();
        let mut previous = 0;
        for (k, &column) in kept.iter().enumerate() {
            gaps.push(if k == 0 { column } else { column - previous });
            previous = column;
        }
        gaps.push(columns - 1 - previous);
        Expander { gaps }
    }

    fn write(&self, out: &mut dyn Write, row: &[u8]) -> io::Result<()> {
        const TABS: [u8; 64] = [b'\t'; 64];
        let tabs = |out: &mut dyn Write, mut count: usize| -> io::Result<()> {
            while count > 0 {
                let take = count.min(TABS.len());
                out.write_all(&TABS[..take])?;
                count -= take;
            }
            Ok(())
        };
        for (k, field) in row.split(|&byte| byte == b'\t').enumerate() {
            tabs(&mut *out, self.gaps[k])?;
            out.write_all(field)?;
        }
        tabs(&mut *out, self.gaps[self.gaps.len() - 1])?;
        out.write_all(b"\n")
    }
}

struct Lookup<'a> {
    index: &'a CompactIndex,
    file: &'a File,
    matcher: Matcher,
    expander: Expander,
    threads: usize,
}

impl Lookup<'_> {
    fn flush(
        &self,
        batch: &[Entry],
        out: &mut dyn Write,
        found_out: &mut Option<&mut dyn Write>,
        errors: &mut dyn Write,
    ) -> Result<()> {
        let mut by_chromosome: Vec<Vec<usize>> = vec![Vec::new(); self.index.chromosomes.len()];
        for (entry, item) in batch.iter().enumerate() {
            if let Entry::Merged(merged) = item {
                by_chromosome[merged.chromosome].push(entry);
            }
        }
        let mut items: Vec<(usize, Vec<usize>)> = Vec::new();
        for (chromosome, mut entries) in by_chromosome.into_iter().enumerate() {
            if entries.is_empty() {
                continue;
            }
            entries.sort_by_key(|&entry| match &batch[entry] {
                Entry::Merged(merged) => merged.position,
                _ => 0,
            });
            let position = |entry: usize| match &batch[entry] {
                Entry::Merged(merged) => merged.position,
                _ => 0,
            };
            let (mut start, mut distinct) = (0, 0);
            for at in 0..entries.len() {
                if at > start && position(entries[at]) != position(entries[at - 1]) {
                    distinct += 1;
                    if distinct >= ITEM_POSITIONS {
                        items.push((chromosome, entries[start..at].to_vec()));
                        start = at;
                        distinct = 0;
                    }
                }
            }
            items.push((chromosome, entries[start..].to_vec()));
        }
        type Found = std::result::Result<Vec<(usize, Vec<u8>)>, String>;
        let results: Vec<Mutex<Option<Found>>> = items.iter().map(|_| Mutex::new(None)).collect();
        let next = AtomicUsize::new(0);
        std::thread::scope(|scope| {
            for _ in 0..self.threads.min(items.len()).max(1) {
                scope.spawn(|| {
                    loop {
                        let item = next.fetch_add(1, Ordering::Relaxed);
                        if item >= items.len() {
                            return;
                        }
                        let (chromosome, entries) = &items[item];
                        let result = lookup_item(
                            self.index,
                            self.file,
                            &self.matcher,
                            batch,
                            *chromosome,
                            entries,
                        );
                        *results[item].lock().unwrap() = Some(result.map_err(|e| e.to_string()));
                    }
                });
            }
        });
        let mut found: Vec<Vec<Vec<u8>>> = (0..batch.len()).map(|_| Vec::new()).collect();
        for result in results {
            for (entry, row) in result.into_inner().unwrap().unwrap().map_err(invalid)? {
                found[entry].push(row);
            }
        }
        for (entry, item) in batch.iter().enumerate() {
            match item {
                Entry::Header(line) | Entry::NotFound(line) => {
                    out.write_all(line)?;
                    out.write_all(b"\n")?;
                }
                Entry::Uncovered(line) => {
                    errors.write_all(b"Encountered uncovered chromosome\n")?;
                    out.write_all(line)?;
                    out.write_all(b"\n")?;
                }
                Entry::Merged(merged) => {
                    if found[entry].is_empty() {
                        out.write_all(&merged.line)?;
                        out.write_all(b"\n")?;
                    } else {
                        for row in &found[entry] {
                            match found_out {
                                Some(writer) => self.expander.write(&mut **writer, row)?,
                                None => self.expander.write(&mut *out, row)?,
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

fn default_threads() -> usize {
    std::env::var("PRIVA_CADD_PRESCORE_THREADS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .filter(|&value| value > 0)
        .unwrap_or_else(|| {
            std::thread::available_parallelism()
                .map(|count| count.get().min(8))
                .unwrap_or(1)
        })
}

/// The lookup itself, on any input and outputs (tests use memory buffers).
fn run_lookup(
    table: &Path,
    input: &mut dyn BufRead,
    out: &mut dyn Write,
    mut found_out: Option<&mut dyn Write>,
    header: bool,
    threads: usize,
    errors: &mut dyn Write,
) -> Result<()> {
    if !table.exists() || !with_suffix(table, ".tbi").exists() {
        return Err(invalid("No valid file with pre-scored variants."));
    }
    writeln!(errors, "Opening {}...", table.display())?;
    let index = load_index(table)?;
    if header {
        for line in &index.header {
            let writer: &mut dyn Write = match &mut found_out {
                Some(writer) => &mut **writer,
                None => &mut *out,
            };
            writer.write_all(line)?;
            writer.write_all(b"\n")?;
        }
    }
    let (reference, alternative) = ref_alt_columns(&index.header);
    let compact_field = |column: usize, name: &str| {
        index
            .kept
            .iter()
            .position(|&kept| kept == column)
            .ok_or_else(|| {
                invalid(format!(
                    "the compared {name} column is not in the compact table"
                ))
            })
    };
    let matcher = Matcher {
        pos: compact_field(1, "Pos")?,
        reference: compact_field(reference, "Ref")?,
        alternative: compact_field(alternative, "Alt")?,
    };
    let file = File::open(with_suffix(table, COMPACT_SUFFIX))?;
    let lookup = Lookup {
        index: &index,
        file: &file,
        matcher,
        expander: Expander::new(index.columns, &index.kept),
        threads: threads.max(1),
    };
    let mut batch = Vec::with_capacity(BATCH_LINES);
    let mut raw = Vec::new();
    loop {
        raw.clear();
        if input.read_until(b'\n', &mut raw)? == 0 {
            break;
        }
        while matches!(raw.last(), Some(b'\n' | b'\r')) {
            raw.pop();
        }
        match classify(raw.clone(), &index) {
            Ok(entry) => batch.push(entry),
            Err(line) => {
                // The Python lookup writes every earlier result, then fails.
                lookup.flush(&batch, out, &mut found_out, errors)?;
                return Err(invalid(format!(
                    "IndexError: list index out of range (input line {:?})",
                    String::from_utf8_lossy(&line)
                )));
            }
        }
        if batch.len() >= BATCH_LINES {
            lookup.flush(&batch, out, &mut found_out, errors)?;
            batch.clear();
        }
    }
    lookup.flush(&batch, out, &mut found_out, errors)?;
    out.flush()?;
    if let Some(writer) = found_out {
        writer.flush()?;
    }
    Ok(())
}

// ---------------------------------------------------------------- commands

/// Arguments in optparse's forms: `-p TABLE`, `-pTABLE`, `--path TABLE`,
/// `--path=TABLE`; positional arguments are ignored, as optparse keeps them.
fn lookup_arguments(
    args: &[String],
) -> Result<(Option<String>, Option<String>, Option<String>, bool)> {
    let (mut path, mut input, mut found, mut header) = (None, None, None, false);
    let mut at = 0;
    while at < args.len() {
        let argument = &args[at];
        let (name, inline) = if let Some(rest) = argument.strip_prefix("--") {
            match rest.split_once('=') {
                Some((name, value)) => (format!("--{name}"), Some(value.to_owned())),
                None => (argument.clone(), None),
            }
        } else if argument.len() > 2 && argument.starts_with('-') && !argument.starts_with("--") {
            (argument[..2].to_owned(), Some(argument[2..].to_owned()))
        } else {
            (argument.clone(), None)
        };
        let mut value = || -> Result<String> {
            if let Some(value) = inline.clone() {
                return Ok(value);
            }
            at += 1;
            args.get(at).cloned().ok_or_else(|| {
                invalid(format!(
                    "{name} option requires an argument\n{LOOKUP_USAGE}"
                ))
            })
        };
        match name.as_str() {
            "-p" | "--path" => path = Some(value()?),
            "-i" | "--input" => input = Some(value()?),
            "--found_out" => found = Some(value()?),
            "--header" => {
                if inline.is_some() {
                    return Err(invalid(format!(
                        "--header option does not take a value\n{LOOKUP_USAGE}"
                    )));
                }
                header = true;
            }
            other if other.starts_with('-') && other != "-" => {
                return Err(invalid(format!("no such option: {other}\n{LOOKUP_USAGE}")));
            }
            _ => {}
        }
        at += 1;
    }
    Ok((path, input, found, header))
}

pub(crate) fn dispatch_lookup(args: &[String]) -> Result<()> {
    if args == ["--help"] || args == ["-h"] {
        println!("{LOOKUP_USAGE}");
        return Ok(());
    }
    let (path, input, found, header) = lookup_arguments(args)?;
    let stdout = io::stdout();
    let mut out = BufWriter::with_capacity(1 << 20, stdout.lock());
    let mut errors = io::stderr();
    let outcome = (|| -> Result<()> {
        let mut reader: Box<dyn BufRead> = match &input {
            Some(path) => Box::new(BufReader::with_capacity(1 << 20, File::open(path)?)),
            None => Box::new(BufReader::with_capacity(1 << 20, io::stdin().lock())),
        };
        let mut found_file = match &found {
            Some(path) => Some(BufWriter::with_capacity(1 << 20, File::create(path)?)),
            None => None,
        };
        let path = path.ok_or_else(|| invalid("TypeError: no prescored table given (-p)"))?;
        run_lookup(
            Path::new(&path),
            reader.as_mut(),
            &mut out,
            found_file.as_mut().map(|writer| writer as &mut dyn Write),
            header,
            default_threads(),
            &mut errors,
        )
    })();
    let _ = out.flush();
    if let Err(error) = outcome {
        eprintln!("priva cadd-prescore: {error}");
        std::process::exit(1);
    }
    Ok(())
}

pub(crate) fn dispatch_build(args: &[String]) -> Result<()> {
    let usage = || invalid(BUILD_USAGE);
    let procedure = args.first().ok_or_else(usage)?;
    if procedure == "--help" || procedure == "-h" {
        println!("{BUILD_USAGE}");
        return Ok(());
    }
    let table = PathBuf::from(args.get(1).ok_or_else(usage)?);
    let (mut threads, mut block_size, mut level) = (default_threads(), DEFAULT_BLOCK, 6u32);
    let mut at = 2;
    while at < args.len() {
        let value = args.get(at + 1).ok_or_else(usage)?;
        let number: usize = value.parse().map_err(|_| usage())?;
        match args[at].as_str() {
            "--threads" if number > 0 => threads = number,
            "--block-size" if (256..=MAX_BLOCK_TEXT).contains(&number) => block_size = number,
            "--level" if (1..=9).contains(&number) => level = number as u32,
            _ => return Err(usage()),
        }
        at += 2;
    }
    match procedure.as_str() {
        "check" => std::process::exit(if compact_matches(&table) { 0 } else { 1 }),
        "build" => {
            let size = fs::metadata(&table)?.len();
            let parts = (threads * 16).max(size.div_ceil(256 << 20) as usize);
            let started = Instant::now();
            match build(&table, threads, block_size, level, parts) {
                Ok((rows, blocks, bytes)) => {
                    println!(
                        "Wrote {} ({rows} rows, {blocks} blocks, {bytes} bytes) and its index in {:.0} s",
                        with_suffix(&table, COMPACT_SUFFIX).display(),
                        started.elapsed().as_secs_f64()
                    );
                    Ok(())
                }
                Err(error) => {
                    eprintln!("priva cadd-prescore-compact: {error}");
                    std::process::exit(1);
                }
            }
        }
        _ => Err(usage()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const COLUMNS: [&str; 9] = [
        "#Chrom",
        "Pos",
        "Ref",
        "Alt",
        "Type",
        "FeatureID",
        "Extra",
        "RawScore",
        "PHRED",
    ];

    struct Folder(PathBuf);

    impl Drop for Folder {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    fn folder(name: &str) -> Folder {
        let path =
            std::env::temp_dir().join(format!("priva_cadd_prescore_{name}_{}", std::process::id()));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Folder(path)
    }

    /// BGZF with `block` bytes of text per block, so rows straddle blocks.
    fn bgzf(text: &[u8], block: usize) -> Vec<u8> {
        let mut out = Vec::new();
        let mut compressor = Compress::new(Compression::new(6), false);
        for chunk in text.chunks(block) {
            deflate(chunk, &mut compressor, &mut out).unwrap();
        }
        out.extend_from_slice(&EOF_BLOCK);
        out
    }

    fn rows() -> Vec<Vec<String>> {
        let mut rows = Vec::new();
        for (chrom, last) in [("1", 120u32), ("2", 40), ("X", 15)] {
            for position in (3..=last).step_by(3) {
                let reference = ["A", "C", "G", "T"][(position % 4) as usize];
                for alt in ["A", "C", "G", "T"].into_iter().filter(|a| *a != reference) {
                    // Two features at every fifth position.
                    let features = if position % 5 == 0 { 2 } else { 1 };
                    for feature in 0..features {
                        rows.push(
                            [
                                chrom.to_owned(),
                                position.to_string(),
                                reference.to_owned(),
                                alt.to_owned(),
                                "SNV".to_owned(),
                                format!("ENST{chrom}{position}{feature}"),
                                format!("extra{}", position * 7),
                                format!("0.{position}"),
                                format!("{}.{feature}", position % 40),
                            ]
                            .to_vec(),
                        );
                    }
                }
            }
        }
        rows
    }

    fn table(directory: &Path) -> PathBuf {
        let mut text = b"##CADD test release\n".to_vec();
        text.extend_from_slice(COLUMNS.join("\t").as_bytes());
        text.push(b'\n');
        for row in rows() {
            text.extend_from_slice(row.join("\t").as_bytes());
            text.push(b'\n');
        }
        let path = directory.join("table.tsv.gz");
        fs::write(&path, bgzf(&text, 333)).unwrap();
        fs::write(with_suffix(&path, ".tbi"), b"").unwrap();
        path
    }

    /// The Python lookup's result, computed row by row.
    fn expected(input: &str) -> (String, String, String) {
        let rows = rows();
        let (mut out, mut found, mut errors) = (String::new(), String::new(), String::new());
        found.push_str("##CADD test release\n");
        found.push_str(&COLUMNS.join("\t"));
        found.push('\n');
        for line in input.lines() {
            if line.starts_with('#') {
                out.push_str(line);
                out.push('\n');
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            let known = ["1", "2", "X"].contains(&fields[0]);
            let number = python_int(fields[1].as_bytes());
            if !known || number.is_none() || number.unwrap() < 1 {
                errors.push_str("Encountered uncovered chromosome\n");
                out.push_str(line);
                out.push('\n');
                continue;
            }
            let (reference, alternative) = (fields[fields.len() - 2], fields[fields.len() - 1]);
            let hits: Vec<&Vec<String>> = rows
                .iter()
                .filter(|row| {
                    row[0] == fields[0]
                        && row[1] == fields[1]
                        && row[2] == reference
                        && row[3] == alternative
                })
                .collect();
            if hits.is_empty() {
                out.push_str(line);
                out.push('\n');
            }
            for row in hits {
                let expanded: Vec<&str> = (0..COLUMNS.len())
                    .map(|column| match column {
                        0 | 1 | 2 | 3 | 5 | 8 => row[column].as_str(),
                        _ => "",
                    })
                    .collect();
                found.push_str(&expanded.join("\t"));
                found.push('\n');
            }
        }
        (out, found, errors)
    }

    fn input() -> String {
        let mut lines = vec![
            "##fileformat=VCFv4.1".to_owned(),
            "#CHROM\tPOS\tID\tREF\tALT".to_owned(),
        ];
        for (chrom, position, reference, alt) in [
            ("1", "3", "T", "A"),
            ("1", "3", "T", "A"),
            ("1", "4", "A", "C"),
            ("1", "15", "T", "G"),
            ("1", "15", "T", "T"),
            ("1", "30", "C", "A"),
            ("1", "120", "A", "G"),
            ("1", "121", "C", "G"),
            ("1", " 12", "A", "C"),
            ("1", "012", "A", "C"),
            ("1", "0", "A", "C"),
            ("1", "abc", "A", "C"),
            ("7", "3", "T", "A"),
            ("2", "1", "A", "C"),
            ("2", "39", "T", "A"),
            ("2", "40", "A", "G"),
            ("X", "15", "T", "C"),
            ("X", "18", "C", "A"),
        ] {
            lines.push(format!("{chrom}\t{position}\t.\t{reference}\t{alt}"));
        }
        // Every position on chromosome 1, so many work items and blocks.
        for position in (1..=125).rev() {
            lines.push(format!("1\t{position}\t.\tG\tT"));
        }
        lines.join("\n") + "\n"
    }

    fn run(table: &Path, input: &str, threads: usize) -> (String, String, String) {
        let (mut out, mut found, mut errors) = (Vec::new(), Vec::new(), Vec::new());
        let mut reader = io::Cursor::new(input.as_bytes().to_vec());
        run_lookup(
            table,
            &mut reader,
            &mut out,
            Some(&mut found),
            true,
            threads,
            &mut errors,
        )
        .unwrap();
        let errors = String::from_utf8(errors).unwrap();
        let errors = errors.split_once('\n').unwrap().1.to_owned();
        (
            String::from_utf8(out).unwrap(),
            String::from_utf8(found).unwrap(),
            errors,
        )
    }

    #[test]
    fn blocks_round_trip() {
        let text = b"chr1\t1\tA\n".repeat(500);
        let data = bgzf(&text, 1000);
        let mut reader = io::Cursor::new(data);
        let (mut block, mut part, mut all) = (Vec::new(), Vec::new(), Vec::new());
        let mut inflater = Decompress::new(false);
        while read_block(&mut reader, &mut block).unwrap() {
            inflate(&block, &mut part, &mut inflater).unwrap();
            all.extend_from_slice(&part);
        }
        assert_eq!(all, text);
    }

    #[test]
    fn python_int_follows_python() {
        assert_eq!(python_int(b" 12 "), Some(12));
        assert_eq!(python_int(b"+7"), Some(7));
        assert_eq!(python_int(b"-3"), Some(-3));
        assert_eq!(python_int(b"1_000"), Some(1000));
        for bad in [&b""[..], b"1__0", b"_1", b"1_", b"1.0", b"abc", b"--1"] {
            assert_eq!(python_int(bad), None);
        }
    }

    #[test]
    fn compact_rows_are_the_kept_columns_in_order() {
        let directory = folder("rows");
        let table = table(&directory.0);
        for parts in [1, 2, 7, 40] {
            build(&table, 3, 300, 6, parts).unwrap();
            let data = fs::read(with_suffix(&table, COMPACT_SUFFIX)).unwrap();
            let mut reader = io::Cursor::new(data);
            let (mut block, mut part, mut text) = (Vec::new(), Vec::new(), Vec::new());
            let mut inflater = Decompress::new(false);
            while read_block(&mut reader, &mut block).unwrap() {
                inflate(&block, &mut part, &mut inflater).unwrap();
                text.extend_from_slice(&part);
            }
            let mut wanted = "##CADD test release\n".to_owned() + &COLUMNS.join("\t") + "\n";
            for row in rows() {
                wanted.push_str(&[0, 1, 2, 3, 5, 8].map(|c| row[c].as_str()).join("\t"));
                wanted.push('\n');
            }
            assert_eq!(String::from_utf8(text).unwrap(), wanted, "parts {parts}");
            assert!(compact_matches(&table));
            let index = load_index(&table).unwrap();
            assert_eq!(
                index.chromosomes,
                [b"1".to_vec(), b"2".to_vec(), b"X".to_vec()]
            );
            assert!(index.positions[0].len() > 10);
        }
    }

    #[test]
    fn lookup_matches_the_python_rules_with_any_thread_count() {
        let directory = folder("lookup");
        let table = table(&directory.0);
        build(&table, 2, 300, 6, 9).unwrap();
        let input = input();
        let wanted = expected(&input);
        for threads in [1, 2, 5] {
            assert_eq!(run(&table, &input, threads), wanted, "threads {threads}");
        }
    }

    #[test]
    fn a_line_without_a_second_field_fails_after_earlier_results() {
        let directory = folder("short");
        let table = table(&directory.0);
        build(&table, 1, 300, 6, 3).unwrap();
        let (mut out, mut errors) = (Vec::new(), Vec::new());
        let mut reader = io::Cursor::new(b"1\t4\t.\tA\tC\nlonely\n1\t3\t.\tT\tA\n".to_vec());
        let outcome = run_lookup(&table, &mut reader, &mut out, None, false, 1, &mut errors);
        assert!(outcome.is_err());
        assert_eq!(out, b"1\t4\t.\tA\tC\n");
    }

    #[test]
    fn a_stale_index_is_refused() {
        let directory = folder("stale");
        let table = table(&directory.0);
        build(&table, 1, 300, 6, 2).unwrap();
        let mut data = fs::read(&table).unwrap();
        data.extend_from_slice(&EOF_BLOCK);
        fs::write(&table, data).unwrap();
        assert!(!compact_matches(&table));
        let (mut out, mut errors) = (Vec::new(), Vec::new());
        let mut reader = io::Cursor::new(Vec::new());
        assert!(run_lookup(&table, &mut reader, &mut out, None, false, 1, &mut errors).is_err());
    }
}
