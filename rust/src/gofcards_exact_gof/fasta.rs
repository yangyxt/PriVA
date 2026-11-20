//! `pyfaidx.Fasta(path, sequence_always_upper=True, rebuild=False)` for the
//! two calls the builder makes: contig length and a sequence slice. A slice
//! past the end of a contig returns what is there, as with `strict_bounds`
//! off. A missing `.fai` is written by `samtools faidx`.
use super::*;
use std::collections::HashMap;
use std::os::unix::fs::FileExt;
use std::process::{Command, Stdio};

struct Contig {
    length: i128,
    offset: i128,
    line_bases: i128,
    line_bytes: i128,
    end_byte: i128,
}

pub(super) struct Fasta {
    file: File,
    contigs: HashMap<String, Contig>,
}

fn index_field(text: &str) -> Result<i128> {
    crate::python_text::python_int(text).ok_or_else(|| {
        error(
            "ValueError",
            format!(
                "invalid literal for int() with base 10: {}",
                crate::python_path::repr(text)
            ),
        )
    })
}

impl Fasta {
    pub(super) fn open(path: &Path) -> Result<Self> {
        let name = path.to_string_lossy();
        let file = File::open(path).map_err(|_| {
            error(
                "FastaNotFoundError",
                format!("Cannot read FASTA from file {name}"),
            )
        })?;
        let lower = name.to_lowercase();
        if lower.ends_with(".gz") || lower.ends_with(".bgz") {
            return Err(error(
                "UnsupportedCompressionFormat",
                format!("Compressed FASTA is not supported: {name}"),
            ));
        }
        let index = PathBuf::from(format!("{name}.fai"));
        if !index.exists() {
            let status = Command::new("samtools")
                .arg("faidx")
                .arg(path)
                .stdin(Stdio::null())
                .status()
                .map_err(|e| os_error(e, Path::new("samtools"), None))?;
            if !status.success() {
                return Err(error(
                    "FastaIndexingError",
                    format!("samtools faidx failed for {name}"),
                ));
            }
        }
        let text = fs::read_to_string(&index).map_err(|e| os_error(e, &index, None))?;
        let mut contigs = HashMap::new();
        for line in text.lines() {
            let fields: Vec<&str> = line.trim_end().split('\t').collect();
            let [contig, length, offset, line_bases, line_bytes] = fields[..] else {
                return Err(error(
                    "ValueError",
                    format!(
                        "not enough values to unpack (expected 5, got {})",
                        fields.len()
                    ),
                ));
            };
            let length = index_field(length)?;
            let offset = index_field(offset)?;
            let line_bases = index_field(line_bases)?;
            let line_bytes = index_field(line_bytes)?;
            let newlines = if line_bases != 0 {
                (length + line_bases - 1) / line_bases * (line_bytes - line_bases)
            } else {
                0
            };
            let record = Contig {
                length,
                offset,
                line_bases,
                line_bytes,
                end_byte: offset + newlines + length,
            };
            if contigs.insert(contig.to_owned(), record).is_some() {
                return Err(error("ValueError", format!("Duplicate key \"{contig}\"")));
            }
        }
        Ok(Self { file, contigs })
    }

    /// The script's `fasta_contig_name`.
    pub(super) fn contig_name(&self, chrom: &str) -> Option<String> {
        [
            chrom.to_owned(),
            format!("chr{chrom}"),
            chrom.strip_prefix("chr").unwrap_or(chrom).to_owned(),
        ]
        .into_iter()
        .find(|candidate| self.contigs.contains_key(candidate))
    }

    pub(super) fn length(&self, contig: &str) -> i128 {
        self.contigs[contig].length
    }

    /// `Faidx.from_file(contig, start, end)` with 1-based coordinates.
    fn read(&self, contig: &str, start: i128, end: i128) -> Result<String> {
        let record = &self.contigs[contig];
        let start0 = start - 1;
        let bases = record.line_bases;
        let newline_bytes = record.line_bytes - bases;
        let before = if start0 > 0 && bases != 0 {
            (start0 - 1) / bases
        } else {
            0
        };
        let to_end = if bases != 0 { end / bases } else { 0 };
        let mut wanted = (to_end - before) * newline_bytes + end - start0;
        let first = record.offset + before * newline_bytes + start0;
        if first + wanted > record.end_byte {
            wanted = record.end_byte - first;
        }
        if wanted <= 0 {
            return Ok(String::new());
        }
        let mut bytes = vec![0u8; wanted as usize];
        let mut filled = 0;
        while filled < bytes.len() {
            let read = self
                .file
                .read_at(&mut bytes[filled..], first as u64 + filled as u64)?;
            if read == 0 {
                break;
            }
            filled += read;
        }
        bytes.truncate(filled);
        Ok(String::from_utf8_lossy(&bytes)
            .replace(['\n', '\r'], "")
            .to_uppercase())
    }

    /// The script's `_fasta_slice`; the error is `str(exception)`.
    pub(super) fn slice(
        &self,
        chrom: &str,
        start: i128,
        length: i128,
    ) -> std::result::Result<String, String> {
        let Some(name) = self.contig_name(chrom) else {
            return Err(crate::python_path::repr(&format!(
                "Chromosome {} not present in FASTA",
                crate::python_path::repr(chrom)
            )));
        };
        if start < 1 {
            return Err(format!("Invalid coordinate {start}"));
        }
        self.read(&name, start, start - 1 + length)
            .map_err(|e| e.to_string())
    }
}
