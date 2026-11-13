//! Lossless internal transport for Python strings containing surrogate code
//! points. This opt-in encoding never appears in published JSON, paths or text.
//! Ordinary reserved Unicode characters are escaped too, so the mapping is
//! injective; no chosen sentinel can collide with a user's literal character.
use std::os::unix::ffi::{OsStrExt, OsStringExt};
use std::path::{Path, PathBuf};
const START: u32 = 0xf0000;
const ESCAPE: u32 = START + 0x800;

pub(crate) fn encode_points(points: impl IntoIterator<Item = u32>) -> String {
    let mut result = String::new();
    for point in points {
        if (0xd800..=0xdfff).contains(&point) {
            result.push(char::from_u32(START + point - 0xd800).unwrap());
        } else {
            if (START..=ESCAPE).contains(&point) {
                result.push(char::from_u32(ESCAPE).unwrap());
            }
            result.push(char::from_u32(point).unwrap());
        }
    }
    result
}
pub(crate) fn encode(text: &str) -> String {
    encode_points(text.chars().map(u32::from))
}
pub(crate) fn points(text: &str) -> Vec<u32> {
    let mut result = Vec::new();
    let mut chars = text.chars().map(u32::from);
    while let Some(point) = chars.next() {
        result.push(if point == ESCAPE {
            chars.next().unwrap_or(ESCAPE)
        } else if (START..ESCAPE).contains(&point) {
            0xd800 + point - START
        } else {
            point
        });
    }
    result
}
pub(crate) fn display(text: &str) -> String {
    points(text)
        .into_iter()
        .map(|p| char::from_u32(p).map_or_else(|| format!("\\u{p:04x}"), |c| c.to_string()))
        .collect()
}
pub(crate) fn encoding_error(points: &[u32], at: usize) -> String {
    let count = points[at..]
        .iter()
        .take_while(|p| (0xd800..=0xdfff).contains(*p))
        .count();
    if count == 1 {
        format!(
            "'utf-8' codec can't encode character '\\u{:04x}' in position {at}: surrogates not allowed",
            points[at]
        )
    } else {
        format!(
            "'utf-8' codec can't encode characters in position {at}-{}: surrogates not allowed",
            at + count - 1
        )
    }
}
/// Match POSIX os.fsencode: surrogateescape is admitted only for DC80–DCFF.
pub(crate) fn native(path: &Path) -> Result<PathBuf, String> {
    let Some(text) = path.to_str() else {
        return Ok(path.to_path_buf());
    };
    let values = points(text);
    let mut bytes = Vec::new();
    for (at, &point) in values.iter().enumerate() {
        if (0xdc80..=0xdcff).contains(&point) {
            bytes.push((point - 0xdc00) as u8);
        } else if let Some(c) = char::from_u32(point) {
            bytes.extend_from_slice(c.encode_utf8(&mut [0; 4]).as_bytes());
        } else {
            return Err(encoding_error(&values, at));
        }
    }
    Ok(PathBuf::from(std::ffi::OsString::from_vec(bytes)))
}
/// Re-enter the transport after a native filesystem operation (including
/// directory entries whose original POSIX names are not UTF-8).
pub(crate) fn from_native(path: &Path) -> PathBuf {
    let mut raw = path.as_os_str().as_bytes();
    let mut values = Vec::new();
    while !raw.is_empty() {
        match std::str::from_utf8(raw) {
            Ok(text) => {
                values.extend(text.chars().map(u32::from));
                break;
            }
            Err(cause) => {
                let n = cause.valid_up_to();
                values.extend(
                    std::str::from_utf8(&raw[..n])
                        .unwrap()
                        .chars()
                        .map(u32::from),
                );
                let end = n + cause.error_len().unwrap_or(raw.len() - n);
                values.extend(raw[n..end].iter().map(|b| 0xdc00 + u32::from(*b)));
                raw = &raw[end..];
            }
        }
    }
    PathBuf::from(encode_points(values))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn all_codepoints_and_reserved_unicode_roundtrip_without_collision() {
        let values = [
            0, 65, 0xd800, 0xdbff, 0xdc80, 0xdfff, START, ESCAPE, 0x10ffff,
        ];
        assert_eq!(points(&encode_points(values)), values);
        assert_ne!(encode_points([0xd800]), encode_points([START]));
        assert_ne!(encode_points([ESCAPE, START]), encode_points([START]));
    }
    #[test]
    fn native_names_keep_surrogateescape_bytes_and_reject_other_surrogates() {
        let path = PathBuf::from(encode_points([47, 97, 0xdc80, 0xdcff]));
        let native = native(&path).unwrap();
        assert_eq!(native.as_os_str().as_bytes(), b"/a\x80\xff");
        assert_eq!(from_native(&native), path);
        assert!(
            self::native(Path::new(&encode_points([0xd800])))
                .unwrap_err()
                .contains("position 0")
        );
    }
}
