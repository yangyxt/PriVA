//! Reading a tabix index (`.tbi`) and fetching the lines of a bgzip-compressed
//! VCF that overlap a region, as htslib does for `pysam.TabixFile.fetch`.
//!
//! The index is bgzip-compressed itself; it lists, per contig, bins of 16 kb
//! upwards with the virtual file offsets (compressed block offset shifted left
//! by 16 bits, plus the offset inside the decompressed block) of the lines in
//! each bin, and a linear index of the first offset per 16 kb window. A query
//! gathers the chunks of the bins that overlap the region, drops chunks that
//! end before the window's first offset, and reads lines from each chunk until
//! it reaches a line that starts at or after the end of the region.
use flate2::read::{DeflateDecoder, MultiGzDecoder};
use std::collections::HashMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

pub(crate) struct TabixIndex {
    pub(crate) names: Vec<String>,
    format: i32,
    column_sequence: i32,
    column_begin: i32,
    column_end: i32,
    meta_character: u8,
    references: Vec<ReferenceIndex>,
}

struct ReferenceIndex {
    bins: HashMap<u32, Vec<(u64, u64)>>,
    linear: Vec<u64>,
}

struct Cursor<'a> {
    data: &'a [u8],
    position: usize,
}

impl<'a> Cursor<'a> {
    fn take(&mut self, count: usize) -> Result<&'a [u8], String> {
        let end = self
            .position
            .checked_add(count)
            .filter(|end| *end <= self.data.len())
            .ok_or_else(|| "the tabix index is truncated".to_owned())?;
        let slice = &self.data[self.position..end];
        self.position = end;
        Ok(slice)
    }

    fn int32(&mut self) -> Result<i32, String> {
        let bytes = self.take(4)?;
        Ok(i32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
    }

    fn uint32(&mut self) -> Result<u32, String> {
        let bytes = self.take(4)?;
        Ok(u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
    }

    fn uint64(&mut self) -> Result<u64, String> {
        let bytes = self.take(8)?;
        let mut array = [0u8; 8];
        array.copy_from_slice(bytes);
        Ok(u64::from_le_bytes(array))
    }

    fn count(&mut self) -> Result<usize, String> {
        let value = self.int32()?;
        usize::try_from(value).map_err(|_| "the tabix index has a negative count".to_owned())
    }
}

pub(crate) fn read_index(path: &Path) -> Result<TabixIndex, String> {
    let file = File::open(path).map_err(|error| format!("{}: {error}", path.display()))?;
    let mut data = Vec::new();
    MultiGzDecoder::new(file)
        .read_to_end(&mut data)
        .map_err(|error| format!("could not read the tabix index {}: {error}", path.display()))?;
    let mut cursor = Cursor {
        data: &data,
        position: 0,
    };
    if cursor.take(4)? != b"TBI\x01" {
        return Err(format!("{} is not a tabix index", path.display()));
    }
    let reference_count = cursor.count()?;
    let format = cursor.int32()?;
    let column_sequence = cursor.int32()?;
    let column_begin = cursor.int32()?;
    let column_end = cursor.int32()?;
    let meta = cursor.int32()?;
    let _skip_lines = cursor.int32()?;
    let names_length = cursor.count()?;
    let names_block = cursor.take(names_length)?;
    let mut names: Vec<String> = Vec::new();
    for part in names_block.split(|byte| *byte == 0).take(reference_count) {
        names.push(String::from_utf8_lossy(part).into_owned());
    }
    if names.len() != reference_count {
        return Err("the tabix index lists fewer contig names than contigs".to_owned());
    }
    let mut references = Vec::new();
    for _ in 0..reference_count {
        let bin_count = cursor.count()?;
        let mut bins = HashMap::new();
        for _ in 0..bin_count {
            let bin = cursor.uint32()?;
            let chunk_count = cursor.count()?;
            let mut chunks = Vec::new();
            for _ in 0..chunk_count {
                let start = cursor.uint64()?;
                let stop = cursor.uint64()?;
                chunks.push((start, stop));
            }
            bins.insert(bin, chunks);
        }
        let window_count = cursor.count()?;
        let mut linear = Vec::new();
        for _ in 0..window_count {
            linear.push(cursor.uint64()?);
        }
        references.push(ReferenceIndex { bins, linear });
    }
    Ok(TabixIndex {
        names,
        format,
        column_sequence,
        column_begin,
        column_end,
        meta_character: u8::try_from(meta).unwrap_or(b'#'),
        references,
    })
}

/// The bins a region `[begin, end)` can fall in (16 kb leaves, five levels).
fn region_bins(begin: i64, end: i64) -> Vec<u32> {
    let last = end - 1;
    let mut bins = vec![0u32];
    for (offset, shift) in [(1i64, 26), (9, 23), (73, 20), (585, 17), (4681, 14)] {
        let mut bin = offset + (begin >> shift);
        while bin <= offset + (last >> shift) {
            bins.push(bin as u32);
            bin += 1;
        }
    }
    bins
}

/// A reader of one bgzip file that can seek to a virtual offset.
struct BgzfReader {
    file: File,
    block: Vec<u8>,
    block_start: u64,
    next_start: u64,
    position: usize,
}

impl BgzfReader {
    fn open(path: &Path) -> Result<BgzfReader, String> {
        let file = File::open(path).map_err(|error| format!("{}: {error}", path.display()))?;
        Ok(BgzfReader {
            file,
            block: Vec::new(),
            block_start: 0,
            next_start: 0,
            position: 0,
        })
    }

    /// Loads the block that starts at `offset`. An empty block (the end-of-file
    /// marker, or the end of the file) leaves `block` empty.
    fn load_block(&mut self, offset: u64) -> Result<(), String> {
        self.block.clear();
        self.position = 0;
        self.block_start = offset;
        self.next_start = offset;
        self.file
            .seek(SeekFrom::Start(offset))
            .map_err(|error| error.to_string())?;
        let mut fixed = [0u8; 12];
        let mut filled = 0;
        while filled < fixed.len() {
            let read = self
                .file
                .read(&mut fixed[filled..])
                .map_err(|error| error.to_string())?;
            if read == 0 {
                break;
            }
            filled += read;
        }
        if filled == 0 {
            return Ok(());
        }
        if filled < fixed.len()
            || fixed[0] != 0x1f
            || fixed[1] != 0x8b
            || fixed[2] != 8
            || fixed[3] & 4 == 0
        {
            return Err("the VCF is not a bgzip file".to_owned());
        }
        let extra_length = u16::from_le_bytes([fixed[10], fixed[11]]) as usize;
        let mut extra = vec![0u8; extra_length];
        self.file
            .read_exact(&mut extra)
            .map_err(|error| error.to_string())?;
        let mut total_size = None;
        let mut at = 0;
        while at + 4 <= extra.len() {
            let sub_length = u16::from_le_bytes([extra[at + 2], extra[at + 3]]) as usize;
            if extra[at] == b'B'
                && extra[at + 1] == b'C'
                && sub_length == 2
                && at + 6 <= extra.len()
            {
                total_size = Some(u16::from_le_bytes([extra[at + 4], extra[at + 5]]) as usize + 1);
            }
            at += 4 + sub_length;
        }
        let total_size = total_size.ok_or_else(|| "the VCF is not a bgzip file".to_owned())?;
        if total_size < 20 + extra_length {
            return Err("a bgzip block is too small".to_owned());
        }
        let compressed_length = total_size - 20 - extra_length;
        let mut compressed = vec![0u8; compressed_length];
        self.file
            .read_exact(&mut compressed)
            .map_err(|error| error.to_string())?;
        let mut footer = [0u8; 8];
        self.file
            .read_exact(&mut footer)
            .map_err(|error| error.to_string())?;
        let expected = u32::from_le_bytes([footer[4], footer[5], footer[6], footer[7]]) as usize;
        let mut decoded = Vec::with_capacity(expected);
        DeflateDecoder::new(&compressed[..])
            .read_to_end(&mut decoded)
            .map_err(|error| format!("a bgzip block is corrupt: {error}"))?;
        if decoded.len() != expected {
            return Err("a bgzip block has the wrong length".to_owned());
        }
        self.block = decoded;
        self.next_start = offset + total_size as u64;
        Ok(())
    }

    fn seek_virtual(&mut self, virtual_offset: u64) -> Result<(), String> {
        self.load_block(virtual_offset >> 16)?;
        let inside = (virtual_offset & 0xffff) as usize;
        if inside > self.block.len() {
            return Err("a virtual offset points past the end of its block".to_owned());
        }
        self.position = inside;
        Ok(())
    }

    /// The virtual offset of the next byte; a fully read block counts as the start of the next.
    fn tell(&self) -> u64 {
        if self.position >= self.block.len() {
            self.next_start << 16
        } else {
            (self.block_start << 16) | self.position as u64
        }
    }

    /// The next line without its line ending, or `None` at the end of the file.
    fn read_line(&mut self) -> Result<Option<Vec<u8>>, String> {
        let mut line: Vec<u8> = Vec::new();
        loop {
            if self.position >= self.block.len() {
                let next = self.next_start;
                self.load_block(next)?;
                if self.block.is_empty() {
                    break;
                }
            }
            let rest = &self.block[self.position..];
            match rest.iter().position(|byte| *byte == b'\n') {
                Some(found) => {
                    line.extend_from_slice(&rest[..found]);
                    self.position += found + 1;
                    if line.last() == Some(&b'\r') {
                        line.pop();
                    }
                    return Ok(Some(line));
                }
                None => {
                    line.extend_from_slice(rest);
                    self.position = self.block.len();
                }
            }
        }
        if line.is_empty() {
            Ok(None)
        } else {
            Ok(Some(line))
        }
    }
}

fn column<'a>(fields: &[&'a str], wanted: i32) -> Result<&'a str, String> {
    let slot =
        usize::try_from(wanted - 1).map_err(|_| "the tabix index names no column".to_owned())?;
    fields
        .get(slot)
        .copied()
        .ok_or_else(|| "a VCF line has too few columns".to_owned())
}

/// Where a line lies: its contig, and its start and end as htslib computes them.
fn line_interval(index: &TabixIndex, line: &[u8]) -> Result<(String, i64, i64), String> {
    let text = std::str::from_utf8(line).map_err(|_| "a VCF line is not UTF-8 text".to_owned())?;
    let fields: Vec<&str> = text.split('\t').collect();
    let contig = column(&fields, index.column_sequence)?.to_owned();
    let position: i64 = column(&fields, index.column_begin)?
        .parse()
        .map_err(|_| "a VCF line has a position that is not an integer".to_owned())?;
    let ucsc = index.format & 0x10000 != 0;
    let begin = if ucsc { position } else { position - 1 };
    let preset = index.format & 0xffff;
    let end = if preset == 2 {
        let mut end = begin + column(&fields, 4)?.len() as i64;
        if let Ok(info) = column(&fields, 8) {
            for item in info.split(';') {
                if let Some(value) = item.strip_prefix("END=") {
                    if let Ok(parsed) = value.parse::<i64>() {
                        if parsed > begin {
                            end = parsed;
                        }
                    }
                }
            }
        }
        end
    } else if preset == 0 {
        if index.column_end > 0 {
            column(&fields, index.column_end)?
                .parse()
                .map_err(|_| "a line has an end that is not an integer".to_owned())?
        } else {
            begin + 1
        }
    } else {
        return Err("only generic and VCF tabix indexes are supported".to_owned());
    };
    Ok((contig, begin, end))
}

/// Calls `visit(number, line)` for each line that overlaps `[begin, end)` on
/// `contig`, in file order, until it returns `Ok(true)` or the lines run out.
pub(crate) fn fetch_lines(
    data_path: &Path,
    index: &TabixIndex,
    contig: &str,
    begin: i64,
    end: i64,
    visit: &mut dyn FnMut(usize, &[u8]) -> Result<bool, String>,
) -> Result<(), String> {
    let id = index
        .names
        .iter()
        .position(|name| name == contig)
        .ok_or_else(|| format!("the index has no contig {contig}"))?;
    let reference = &index.references[id];
    let minimum = if reference.linear.is_empty() {
        0
    } else {
        let window = (begin >> 14) as usize;
        reference.linear[window.min(reference.linear.len() - 1)]
    };
    let mut chunks: Vec<(u64, u64)> = Vec::new();
    for bin in region_bins(begin, end) {
        if let Some(list) = reference.bins.get(&bin) {
            for (start, stop) in list {
                if *stop > minimum {
                    chunks.push((*start, *stop));
                }
            }
        }
    }
    chunks.sort();
    let mut merged: Vec<(u64, u64)> = Vec::new();
    for (start, stop) in chunks {
        if let Some(last) = merged.last_mut() {
            if start <= last.1 {
                if stop > last.1 {
                    last.1 = stop;
                }
                continue;
            }
        }
        merged.push((start, stop));
    }
    let mut reader = BgzfReader::open(data_path)?;
    let mut number = 0usize;
    for (start, stop) in merged {
        reader.seek_virtual(start)?;
        while reader.tell() < stop {
            let line = match reader.read_line()? {
                Some(line) => line,
                None => break,
            };
            if line.is_empty() || line.first() == Some(&index.meta_character) {
                continue;
            }
            let (line_contig, line_begin, line_end) = line_interval(index, &line)?;
            if line_contig != contig || line_begin >= end {
                return Ok(());
            }
            if line_end > begin {
                if visit(number, &line)? {
                    return Ok(());
                }
                number += 1;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bins_cover_the_region_at_every_level() {
        // Position 10 (0-based 9..10) lies in leaf bin 4681 and its ancestors.
        assert_eq!(region_bins(9, 10), vec![0, 1, 9, 73, 585, 4681]);
        // A region spanning a 16 kb boundary names two leaf bins.
        let bins = region_bins(16380, 16390);
        assert!(bins.contains(&4681) && bins.contains(&4682));
    }
}
