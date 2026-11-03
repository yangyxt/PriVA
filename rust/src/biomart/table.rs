use super::*;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, BufWriter, Write};

fn digits(value: &str) -> bool {
    !value.is_empty()
        && value
            .chars()
            .all(|c| crate::python_text::decimal_digit_value(c).is_some())
}
fn stable(value: &str, prefix: &str) -> bool {
    value.strip_prefix(prefix).is_some_and(digits)
}
fn refseq(value: &str) -> bool {
    let bytes = value.as_bytes();
    if bytes.len() < 4
        || !bytes[0].is_ascii_uppercase()
        || !bytes[1].is_ascii_uppercase()
        || bytes[2] != b'_'
    {
        return false;
    }
    match value[3..].split_once('.') {
        Some((a, b)) => digits(a) && digits(b),
        None => digits(&value[3..]),
    }
}
/// Binary lines and per-line UTF-8 decode intentionally differ from a universal
/// newline TSV reader. Publish every original accepted row, including repeats.
pub(super) fn validate(path: &Path, raw: bool, output: Option<&Path>) -> Result<Json> {
    // The original truncates its destination before it even opens the input.
    let mut destination = output
        .map(|p| {
            File::create(p)
                .map(BufWriter::new)
                .map_err(|e| os_error(e, p, None))
        })
        .transpose()?;
    let result: Result<(Json, bool, usize)> = (|| {
        let input = File::open(path).map_err(|e| os_error(e, path, None))?;
        let mut input = BufReader::with_capacity(8192, input);
        let mut line = Vec::new();
        input
            .read_until(b'\n', &mut line)
            .map_err(|e| os_error(e, path, None))?;
        let header =
            crate::python_io::decode_utf8(&line).map_err(|e| error("UnicodeDecodeError", e))?;
        let observed = Json::Array(
            header
                .trim_end_matches(['\r', '\n'])
                .split('\t')
                .map(string)
                .collect(),
        );
        require(
            equal(&observed, get(facts(), "headers")),
            format!(
                "BioMart table must have exactly the pinned five-column schema: {}",
                path.display()
            ),
        )?;
        if let Some(out) = &mut destination {
            out.write_all(&line)
                .map_err(|e| os_error(e, output.unwrap(), None))?;
        }
        let mut counts = [0usize; 5];
        let (mut genes, mut transcripts, mut accessions, mut pairs, mut full_rows, mut ignored) = (
            HashSet::new(),
            HashSet::new(),
            HashSet::new(),
            HashSet::new(),
            HashSet::new(),
            HashSet::new(),
        );
        let mut mappings: HashMap<String, HashSet<String>> = HashMap::new();
        let mut completed = false;
        let mut line_number = 1;
        loop {
            line.clear();
            if input
                .read_until(b'\n', &mut line)
                .map_err(|e| os_error(e, path, None))?
                == 0
            {
                break;
            }
            line_number += 1;
            if raw && line.trim_ascii_end_matches_cr_lf() == b"[success]" {
                let mut rest = Vec::new();
                input
                    .read_to_end(&mut rest)
                    .map_err(|e| os_error(e, path, None))?;
                require(
                    !completed
                        && rest
                            .iter()
                            .all(|b| matches!(b, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)),
                    "BioMart completion stamp is not the single terminal record",
                )?;
                completed = true;
                break;
            }
            let decoded =
                crate::python_io::decode_utf8(&line).map_err(|e| error("UnicodeDecodeError", e))?;
            let cells = decoded
                .trim_end_matches(['\r', '\n'])
                .split('\t')
                .map(str::to_owned)
                .collect::<Vec<_>>();
            require(
                cells.len() == 5 && stable(&cells[0], "ENSG") && stable(&cells[1], "ENST"),
                format!(
                    "Malformed/truncated BioMart record at line {line_number}: {}",
                    path.display()
                ),
            )?;
            let nm = &cells[4];
            require(
                !nm.starts_with("NM_") || stable(nm, "NM_"),
                format!(
                    "BioMart NM accessions must be bare IDs for the unchanged motif reader (line {line_number})"
                ),
            )?;
            require(
                nm.is_empty() || refseq(nm),
                format!(
                    "Malformed BioMart RefSeq accession at line {line_number}: {}",
                    path.display()
                ),
            )?;
            counts[0] += 1;
            genes.insert(cells[0].clone());
            transcripts.insert(cells[1].clone());
            full_rows.insert(cells.clone());
            if nm.starts_with("NM_") {
                counts[1] += 1;
                accessions.insert(nm.clone());
                pairs.insert((nm.clone(), cells[1].clone()));
                mappings
                    .entry(nm.clone())
                    .or_default()
                    .insert(cells[1].clone());
            } else {
                counts[2] += 1;
                if nm.is_empty() {
                    counts[4] += 1;
                } else {
                    counts[3] += 1;
                    ignored.insert(nm.clone());
                }
            }
            if let Some(out) = &mut destination {
                out.write_all(&line)
                    .map_err(|e| os_error(e, output.unwrap(), None))?;
            }
        }
        // These checks occur after both handles close in the Python source;
        // finish them below, after destination flush, using the returned state.
        let mut ignored = ignored.into_iter().collect::<Vec<String>>();
        ignored.sort();
        let mut result = Json::Object(
            [
                "rows",
                "nm_rows",
                "rows_without_nm",
                "non_nm_refseq_rows",
                "rows_without_refseq",
            ]
            .into_iter()
            .zip(counts)
            .map(|(k, v)| (k.into(), Json::Int(v as i128)))
            .collect(),
        );
        result.set(
            "non_nm_refseq_accessions",
            Json::Array(ignored.into_iter().map(string).collect()),
        );
        for (key, value) in [
            ("genes", genes.len()),
            ("transcripts", transcripts.len()),
            ("nm_accessions", accessions.len()),
            ("unique_nm_enst_pairs", pairs.len()),
            ("duplicate_full_rows", counts[0] - full_rows.len()),
            ("duplicate_nm_enst_rows", counts[1] - pairs.len()),
            (
                "nm_with_multiple_enst",
                mappings.values().filter(|v| v.len() > 1).count(),
            ),
            (
                "max_enst_per_nm",
                mappings.values().map(HashSet::len).max().unwrap_or(0),
            ),
        ] {
            result.set(key, Json::Int(value as i128));
        }
        result.set(
            "completion_stamp_verified",
            if raw {
                Json::Bool(completed)
            } else {
                Json::Null
            },
        );
        Ok((result, completed, counts[1]))
    })();
    if let Some(mut out) = destination {
        out.flush()
            .map_err(|e| os_error(e, output.unwrap(), None))?;
    }
    let (result, completed, nm_rows) = result?;
    require(
        !raw || completed,
        "BioMart export is incomplete: missing terminal [success] completion stamp",
    )?;
    require(
        nm_rows > 0,
        "BioMart export contains no usable NM-to-ENST mappings",
    )?;
    Ok(result)
}

trait TrimCrLf {
    fn trim_ascii_end_matches_cr_lf(&self) -> &[u8];
}
impl TrimCrLf for [u8] {
    fn trim_ascii_end_matches_cr_lf(&self) -> &[u8] {
        let n = self
            .iter()
            .rposition(|b| !matches!(b, b'\r' | b'\n'))
            .map_or(0, |i| i + 1);
        &self[..n]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn python_decimal_digits_are_distinct_from_numeric_characters() {
        assert!(stable("ENSG١２0", "ENSG"));
        assert!(!stable("ENSG²", "ENSG"));
        assert!(refseq("XM_١.２"));
        assert!(!refseq("XM_1.2.3"));
        assert!(!refseq("ΑΜ_1"));
    }
}
