//! `scripts/clinvar_vcv.py` (1,061 lines): the ClinVar VCV helpers; the step uses only `open_text`.
//!
//! Thread T2 (gene-condition hub). Converted: `open_text` in its reading mode, which the hub's
//! caches (the condition cache and the non-LOF mechanism JSON) are read through. The XML parser
//! and the writing mode stay with the group 4 builders.
use super::gene_mechanism_common::{python_error, python_path_text, python_text_repr, refused};
use super::types::{HarnessCall, not_converted};
use crate::cli::Result;
use crate::pandas_table::Table;
use std::io::Read;

/// `with open_text(path) as handle: handle.read()`: the file's text, gunzipped
/// when the name ends in `.gz`, decoded as UTF-8, with Python's universal
/// newlines (`\r\n` and `\r` read as `\n`). A missing file is Python's
/// FileNotFoundError (naming the path as `str(Path(path))` prints it); text
/// that is not UTF-8 is refused.
pub(crate) fn open_text(path: &str) -> Result<String> {
    let mut raw = Vec::new();
    match std::fs::File::open(path) {
        Ok(mut file) => {
            file.read_to_end(&mut raw)
                .map_err(|error| refused(format!("reading {path}: {error}")))?;
        }
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            return Err(python_error(
                "FileNotFoundError",
                format!(
                    "[Errno 2] No such file or directory: {}",
                    python_text_repr(&python_path_text(path))?
                ),
            ));
        }
        Err(error) => return Err(refused(format!("opening {path}: {error}"))),
    }
    if path.ends_with(".gz") {
        let mut plain = Vec::new();
        flate2::read::MultiGzDecoder::new(&raw[..])
            .read_to_end(&mut plain)
            .map_err(|error| refused(format!("decompressing {path}: {error}")))?;
        raw = plain;
    }
    let text = String::from_utf8(raw).map_err(|_| refused(format!("{path} is not UTF-8")))?;
    Ok(if text.contains('\r') {
        text.replace("\r\n", "\n").replace('\r', "\n")
    } else {
        text
    })
}

/// `priva acmg-step clinvar_vcv <function>`: runs one converted function on a frame.
pub(crate) fn run_harness(function: &str, call: &HarnessCall) -> Result<Table> {
    let _ = call;
    Err(not_converted("clinvar_vcv", function))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plain_and_gzip_text_with_universal_newlines() {
        let folder = std::env::temp_dir().join(format!("t2_open_text_{}", std::process::id()));
        std::fs::create_dir_all(&folder).unwrap();
        let plain = folder.join("a.json");
        std::fs::write(&plain, "{\"a\":\r\n1}\r").unwrap();
        assert_eq!(open_text(plain.to_str().unwrap()).unwrap(), "{\"a\":\n1}\n");
        let packed = folder.join("a.json.gz");
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        std::io::Write::write_all(&mut encoder, b"[1]").unwrap();
        std::fs::write(&packed, encoder.finish().unwrap()).unwrap();
        assert_eq!(open_text(packed.to_str().unwrap()).unwrap(), "[1]");
        let missing = folder.join("missing.json");
        assert!(
            open_text(missing.to_str().unwrap())
                .unwrap_err()
                .to_string()
                .starts_with("FileNotFoundError: [Errno 2] No such file or directory: '")
        );
        std::fs::remove_dir_all(&folder).unwrap();
    }
}
