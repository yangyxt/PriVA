//! Literal streaming projection of `scripts/alphmis_tsv2vcf.py`.
//!
//! Index contig order, first-line skipping, Python whitespace and every selected
//! field lexeme are retained. The original ignores bgzip/tabix exit failures;
//! this opt-in entry point reports them but retains that outcome. Its caller
//! must validate the resulting VCF, as the original installer does. Child
//! programs receive literal argv paths: shell expansion is never performed.
use crate::cli::{Result, invalid};
use crate::python_args::{Parsed, parse_long_options};
use crate::python_io::PythonGzipReader;
use crate::python_text::is_python_space;
use std::fs::{self, File};
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const HELP: &str = "usage: priva alphamissense-vcf build --input INDEXED_TSV_GZ --output VCF_GZ\n\
Project the original AlphaMissense columns without interpreting scores.\n\
Requires the input's .tbi path and existing bgzip/tabix executables. Output ends in .gz.\n\
Paths are literal; shell expansion is not supported. As in the original Python\n\
procedure, child command failures are reported but do not change its exit status.\n\
Validate the output and index before using or publishing them.";

const HEADER: &str = "##fileformat=VCFv4.2\n\
##INFO=<ID=UNIPROT,Number=1,Type=String,Description=\"UniProt ID\">\n\
##INFO=<ID=TRANSCRIPT,Number=1,Type=String,Description=\"Ensembl Transcript ID\">\n\
##INFO=<ID=PVAR,Number=1,Type=String,Description=\"Protein Variant\">\n\
##INFO=<ID=AM_PATHOGENICITY,Number=1,Type=Float,Description=\"Pathogenicity Score\">\n\
##INFO=<ID=AM_CLASS,Number=1,Type=String,Description=\"Classification\">\n";
const COLUMNS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

fn appended(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(suffix);
    name.into()
}

/// Read just the contig catalogue used by pysam.TabixFile.contigs. The source
/// explicitly chooses INPUT.tbi; that file can contain TBI or CSI index bytes.
fn contigs(path: &Path) -> Result<Vec<String>> {
    struct Cursor<'a>(&'a [u8]);
    impl<'a> Cursor<'a> {
        fn take(&mut self, count: usize) -> Result<&'a [u8]> {
            if self.0.len() < count {
                return Err(invalid("the tabix index is truncated"));
            }
            let (part, tail) = self.0.split_at(count);
            self.0 = tail;
            Ok(part)
        }
        fn count(&mut self) -> Result<usize> {
            let value = i32::from_le_bytes(self.take(4)?.try_into().unwrap());
            usize::try_from(value).map_err(|_| invalid("negative tabix index count"))
        }
        fn array(&mut self, width: usize) -> Result<()> {
            let bytes = self
                .count()?
                .checked_mul(width)
                .ok_or_else(|| invalid("index count overflow"))?;
            self.take(bytes)?;
            Ok(())
        }
    }
    let mut data = Vec::new();
    PythonGzipReader::new(File::open(path)?).read_to_end(&mut data)?;
    let mut cursor = Cursor(&data);
    let magic = cursor.take(4)?;
    let (csi, references, names) = match magic {
        b"TBI\x01" => {
            let references = cursor.count()?;
            cursor.take(24)?;
            let length = cursor.count()?;
            (false, references, cursor.take(length)?)
        }
        b"CSI\x01" => {
            cursor.count()?; // min_shift
            cursor.count()?; // depth
            let length = cursor.count()?;
            let mut aux = Cursor(cursor.take(length)?);
            aux.take(24)?;
            let length = aux.count()?;
            let names = aux.take(length)?;
            (true, cursor.count()?, names)
        }
        _ => return Err(invalid("the input .tbi file is not a TBI or CSI index")),
    };
    let names = names
        .split(|byte| *byte == 0)
        .filter(|name| !name.is_empty())
        .map(|name| {
            std::str::from_utf8(name)
                .map(str::to_owned)
                .map_err(|_| invalid("index contig is not UTF-8"))
        })
        .collect::<Result<Vec<_>>>()?;
    if names.len() != references {
        return Err(invalid("tabix index contig catalogue has the wrong length"));
    }
    for _ in 0..references {
        let bins = cursor.count()?;
        for _ in 0..bins {
            cursor.take(if csi { 12 } else { 4 })?;
            cursor.array(16)?;
        }
        if !csi {
            cursor.array(8)?;
        }
    }
    Ok(names)
}

/// Pysam opens the data before loading its index. For generic TSV, htslib
/// detects text from the first 1024 decoded bytes; ordinary gzip lacks the
/// virtual tell position required by TabixFile. Plain text opens here and is
/// subsequently refused by the source's separate gzip reader.
fn check_tabix_data(path: &Path) -> Result<()> {
    let mut prefix = [0_u8; 18];
    let length = File::open(path)?.read(&mut prefix)?;
    let compressed = length >= 2 && prefix[..2] == [31, 139];
    if compressed && (length < 18 || prefix[3] & 4 == 0 || &prefix[12..16] != b"BC\x02\x00") {
        return Err(invalid("TabixFile requires BGZF for compressed input"));
    }
    let mut reader: Box<dyn Read> = if compressed {
        Box::new(PythonGzipReader::new(File::open(path)?))
    } else {
        Box::new(File::open(path)?)
    };
    let mut peek = Vec::new();
    reader.by_ref().take(1024).read_to_end(&mut peek)?;
    if peek
        .iter()
        .any(|byte| *byte < b' ' && !matches!(byte, b'\t' | b'\r' | b'\n'))
    {
        return Err(invalid("the indexed input is not generic text"));
    }
    Ok(())
}

/// TextIOWrapper's UTF-8 decoding and universal CR/LF/CRLF line boundaries.
/// Decode each input block before returning its lines, so invalid UTF-8 does
/// not silently publish earlier lines from the same decoded block.
struct Lines<R: Read> {
    source: R,
    text: String,
    at: usize,
    pending: Vec<u8>,
    eof: bool,
    skip_lf: bool,
}
impl<R: Read> Lines<R> {
    fn new(source: R) -> Self {
        Self {
            source,
            text: String::new(),
            at: 0,
            pending: Vec::new(),
            eof: false,
            skip_lf: false,
        }
    }
    fn fill(&mut self) -> Result<bool> {
        if self.eof {
            return Ok(false);
        }
        let mut block = [0_u8; 8192];
        let count = self.source.read(&mut block)?;
        self.eof = count == 0;
        self.pending.extend_from_slice(&block[..count]);
        let complete = match std::str::from_utf8(&self.pending) {
            Ok(_) => self.pending.len(),
            Err(error) if error.error_len().is_none() && !self.eof => error.valid_up_to(),
            Err(_) => return Err(invalid("UnicodeDecodeError: input is not UTF-8")),
        };
        self.text = std::str::from_utf8(&self.pending[..complete])
            .unwrap()
            .to_owned();
        self.pending.drain(..complete);
        self.at = 0;
        Ok(!self.eof || !self.text.is_empty())
    }
    fn next(&mut self) -> Result<Option<String>> {
        let mut line = String::new();
        loop {
            if self.at == self.text.len() && !self.fill()? {
                return Ok((!line.is_empty()).then_some(line));
            }
            if self.at == self.text.len() {
                continue;
            }
            if self.skip_lf {
                self.skip_lf = false;
                if self.text.as_bytes()[self.at] == b'\n' {
                    self.at += 1;
                    continue;
                }
            }
            let rest = &self.text[self.at..];
            if let Some(end) = rest.bytes().position(|byte| matches!(byte, b'\r' | b'\n')) {
                line.push_str(&rest[..end]);
                self.skip_lf = rest.as_bytes()[end] == b'\r';
                self.at += end + 1;
                return Ok(Some(line));
            }
            line.push_str(rest);
            self.at = self.text.len();
        }
    }
}

fn project(input: impl Read, output: &mut impl Write, names: &[String]) -> Result<usize> {
    output.write_all(HEADER.as_bytes())?;
    for name in names {
        writeln!(output, "##contig=<ID={name}>")?;
    }
    if names.is_empty() {
        writeln!(output)?;
    }
    output.write_all(COLUMNS.as_bytes())?;
    let mut lines = Lines::new(input);
    let _ = lines.next()?; // The source always discards the first physical text line.
    let mut count = 0;
    while let Some(line) = lines.next()? {
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line
            .split(is_python_space)
            .filter(|field| !field.is_empty())
            .collect();
        if fields.len() < 10 {
            return Err(invalid("IndexError: list index out of range"));
        }
        writeln!(
            output,
            "{}\t{}\t.\t{}\t{}\t.\tPASS\tUNIPROT={};TRANSCRIPT={};PVAR={};AM_PATHOGENICITY={};AM_CLASS={}",
            fields[0],
            fields[1],
            fields[2],
            fields[3],
            fields[5],
            fields[6],
            fields[7],
            fields[8],
            fields[9]
        )?;
        count += 1;
    }
    Ok(count)
}

fn child_status(command: &mut Command) -> i32 {
    match command.status() {
        Ok(status) => status.code().unwrap_or(128),
        Err(error) => {
            eprintln!("priva alphamissense-vcf: child command failed: {error}");
            if error.kind() == std::io::ErrorKind::NotFound {
                127
            } else {
                126
            }
        }
    }
}

fn build(input: &Path, output: &Path) -> Result<String> {
    if !output.to_string_lossy().ends_with(".gz") {
        return Err(invalid("ValueError: Output VCF file must end with .gz"));
    }
    let index = appended(input, ".tbi");
    fs::metadata(input)?;
    fs::metadata(&index)?;
    check_tabix_data(input)?;
    let names = contigs(&index)?;
    let temporary = appended(output, ".tmp.vcf");
    let source = PythonGzipReader::new(File::open(input)?);
    let mut writer = BufWriter::new(File::create(&temporary)?);
    let rows = project(source, &mut writer, &names)?;
    writer.flush()?;
    drop(writer);
    // Absolute argv paths prevent filenames beginning with '-' from becoming
    // child options. There is no shell, quoting, expansion, or command string.
    let bgzip = match File::create(output) {
        Ok(file) => child_status(
            Command::new("bgzip")
                .args(["-f", "-c"])
                .arg(std::path::absolute(&temporary)?)
                .stdout(Stdio::from(file)),
        ),
        Err(error) => {
            eprintln!("priva alphamissense-vcf: output redirection failed: {error}");
            2
        }
    };
    let tabix = if bgzip == 0 {
        Some(child_status(
            Command::new("tabix")
                .args(["-p", "vcf", "-f"])
                .arg(std::path::absolute(output)?),
        ))
    } else {
        None
    };
    // Deliberately source-matched: subprocess.run(..., shell=True) has no
    // check=True, and the source removes the temporary even after child failure.
    fs::remove_file(&temporary)?;
    Ok(format!(
        "{{\"rows\":{rows},\"contigs\":{},\"bgzip_status\":{bgzip},\"tabix_status\":{}}}",
        names.len(),
        tabix.map_or("null".to_owned(), |status| status.to_string())
    ))
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args.is_empty() || args == ["--help"] || args == ["-h"] {
        println!("{HELP}");
        return Ok(());
    }
    if args[0] != "build" {
        return Err(invalid(HELP));
    }
    match parse_long_options(&args[1..], &["input", "output"]).map_err(invalid)? {
        Parsed::Help => println!("{HELP}"),
        Parsed::Values(values) => {
            let required = |name: &str| {
                values
                    .get(name)
                    .map(PathBuf::from)
                    .ok_or_else(|| invalid(format!("missing --{name}\n{HELP}")))
            };
            match build(&required("input")?, &required("output")?) {
                Ok(summary) => println!("{summary}"),
                Err(error) => {
                    eprintln!("priva alphamissense-vcf: {error}");
                    std::process::exit(1);
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn projection_keeps_literal_fields_duplicates_and_index_order() {
        let row = "chr2\u{1c}020\u{a0}G C - P ENST p.G2C -0.000 class extra\r";
        let source = format!("discarded\r#ignored\r{row}{row}");
        let mut output = Vec::new();
        assert_eq!(
            project(
                source.as_bytes(),
                &mut output,
                &["chr2".into(), "chr1".into()]
            )
            .unwrap(),
            2
        );
        let text = String::from_utf8(output).unwrap();
        assert!(text.contains("##contig=<ID=chr2>\n##contig=<ID=chr1>"));
        assert_eq!(text.matches("chr2\t020\t.\tG\tC\t.\tPASS\tUNIPROT=P;TRANSCRIPT=ENST;PVAR=p.G2C;AM_PATHOGENICITY=-0.000;AM_CLASS=class\n").count(), 2);
    }
    #[test]
    fn universal_newlines_cross_decode_blocks_and_utf8_failure_is_not_lossy() {
        let data = format!("{}\r\nβ\rz\n", "a".repeat(8191));
        let mut lines = Lines::new(data.as_bytes());
        assert_eq!(lines.next().unwrap().unwrap().len(), 8191);
        assert_eq!(lines.next().unwrap().as_deref(), Some("β"));
        assert_eq!(lines.next().unwrap().as_deref(), Some("z"));
        assert_eq!(lines.next().unwrap(), None);
        assert!(Lines::new(&b"first\n\xff"[..]).next().is_err());
    }
}
