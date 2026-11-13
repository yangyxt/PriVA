use super::*;
use crate::python_lines::Lines;
use std::cell::RefCell;
use std::collections::{BTreeMap, HashSet};
use std::fs::File;
use std::io::Read;
use std::rc::Rc;

/// Python 3.12 / Unicode 15.0 `str.isdigit`, including digits that int()
/// subsequently refuses. The FAI offset field is checked but never converted.
pub(crate) fn isdigit(text: &str) -> bool {
    !text.is_empty() && text.chars().all(|c| crate::python_text::decimal_digit_value(c).is_some()
        || matches!(c as u32,0xb2..=0xb3 | 0xb9 | 0x1369..=0x1371 | 0x19da | 0x2070 | 0x2074..=0x2079
            | 0x2080..=0x2089 | 0x2460..=0x2468 | 0x2474..=0x247c | 0x2488..=0x2490 | 0x24ea
            | 0x24f5..=0x24fd | 0x24ff | 0x2776..=0x277e | 0x2780..=0x2788 | 0x278a..=0x2792
            | 0x10a40..=0x10a43 | 0x10e60..=0x10e68 | 0x11052..=0x1105a | 0x1f100..=0x1f10a))
}
pub(super) fn decimal(text: &str) -> Result<String> {
    let mut digits = String::new();
    for c in text.chars() {
        let value = crate::python_text::decimal_digit_value(c).ok_or_else(|| {
            error(
                "ValueError",
                format!("invalid literal for int() with base 10: {}", repr(text)),
            )
        })?;
        digits.push(char::from(b'0' + value as u8));
    }
    if digits.len() > 4300 {
        return Err(error(
            "ValueError",
            format!(
                "Exceeds the limit (4300 digits) for integer string conversion: value has {} digits; use sys.set_int_max_str_digits() to increase the limit",
                digits.len()
            ),
        ));
    }
    let trimmed = digits.trim_start_matches('0');
    Ok(if trimmed.is_empty() {
        "0".into()
    } else {
        trimmed.into()
    })
}
fn at_least(a: &str, b: &str) -> bool {
    a.len() > b.len() || (a.len() == b.len() && a >= b)
}

/// Observe only the bytes actually read. This supplies the source-specific
/// UnicodeDecodeError diagnostic while the shared line reader keeps Alpha's
/// existing method bodies/error handling unchanged. Only a trailing incomplete
/// UTF-8 character is retained between reads; no second file read is made.
struct DiagnosticRead<R> {
    source: R,
    pending: Vec<u8>,
    diagnostic: Rc<RefCell<Option<String>>>,
}
impl<R: Read> Read for DiagnosticRead<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        let count = self.source.read(out)?;
        self.pending.extend_from_slice(&out[..count]);
        match std::str::from_utf8(&self.pending) {
            Ok(_) => self.pending.clear(),
            Err(e) if e.error_len().is_none() && count > 0 => {
                self.pending.drain(..e.valid_up_to());
            }
            Err(_) => {
                *self.diagnostic.borrow_mut() = crate::python_io::decode_utf8(&self.pending).err();
            }
        }
        Ok(count)
    }
}
pub(super) fn lengths(
    reference: &Path,
    assembly: &str,
) -> Result<(BTreeMap<String, String>, Json)> {
    require(
        is_file(reference)? && metadata(reference)?.len() > 0,
        format!(
            "Missing or empty configured reference: {}",
            reference.display()
        ),
    )?;
    let index = PathBuf::from(format!("{}.fai", reference.display()));
    require(
        is_file(&index)? && metadata(&index)?.len() > 0,
        format!(
            "Missing configured reference FASTA index: {}; index the selected reference first",
            index.display()
        ),
    )?;
    let before = identity(&index)?;
    let source = File::open(&index).map_err(|e| os_error(e, &index, None))?;
    let diagnostic = Rc::new(RefCell::new(None));
    let mut lines = Lines::new(DiagnosticRead {
        source,
        pending: Vec::new(),
        diagnostic: diagnostic.clone(),
    });
    let mut seen = HashSet::new();
    let mut chromosomes = Vec::new();
    loop {
        let line = lines.next().map_err(|e| {
            if let Some(message) = diagnostic.borrow().as_ref() {
                error("UnicodeDecodeError", message.clone())
            } else if let Some(cause) = e.downcast_ref::<io::Error>() {
                if let Some(code) = cause.raw_os_error() {
                    os_error(io::Error::from_raw_os_error(code), &index, None)
                } else {
                    e
                }
            } else {
                e
            }
        })?;
        let Some(line) = line else {
            break;
        };
        let columns: Vec<_> = line.split('\t').collect();
        let malformed = format!(
            "Malformed configured reference FASTA index: {}",
            index.display()
        );
        require(
            columns.len() >= 5
                && !columns[0].is_empty()
                && columns[1..5].iter().all(|v| isdigit(v)),
            malformed.clone(),
        )?;
        let length = decimal(columns[1])?;
        require(length != "0", malformed.clone())?;
        let bases = decimal(columns[3])?;
        require(bases != "0", malformed.clone())?;
        let bytes = decimal(columns[4])?;
        require(at_least(&bytes, &bases), malformed)?;
        require(
            seen.insert(columns[0].to_owned()),
            format!("Duplicate reference chromosome: {}", columns[0]),
        )?;
        chromosomes.push((columns[0].to_owned(), length));
    }
    let mut normalized = BTreeMap::new();
    for (name, length) in chromosomes {
        let key = normalize(&name);
        require(
            normalized.get(&key).is_none_or(|old| *old == length),
            format!("Conflicting normalized chromosome lengths in reference: {name}"),
        )?;
        normalized.insert(key, length);
    }
    require(
        normalized.get("1").map(String::as_str) == Some(chr1(assembly)),
        format!(
            "Configured reference chr1 length does not identify {assembly}: {}",
            reference.display()
        ),
    )?;
    let record = object([
        ("path", string(absolute(reference)?)),
        ("identity", identity(reference)?),
        ("fai", string(absolute(&index)?)),
        ("fai_sha256", string(sha256(&index)?)),
    ]);
    require(
        json_text::py_equal(&identity(&index)?, &before),
        "Configured reference index changed during validation",
    )?;
    Ok((normalized, record))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn decimal_and_digit_checks_are_separate_and_not_machine_integer_limited() {
        assert!(isdigit("²"));
        assert!(decimal("²").is_err());
        assert_eq!(decimal("００٢４").unwrap(), "24");
        assert_eq!(decimal(&"9".repeat(4300)).unwrap().len(), 4300);
        assert!(decimal(&"9".repeat(4301)).is_err());
        assert!(!isdigit("1_0"));
        assert!(!isdigit("Ⅻ"));
        assert!(at_least(&"1".repeat(200), &"9".repeat(199)));
    }
}
