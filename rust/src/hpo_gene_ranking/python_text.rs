//! Python 3.12 rules that the converted scripts' outputs depend on: how a float
//! is printed (`repr`, `format(x, ".6f")`, `format(x, "g")`), `round(x, n)`,
//! `sum()` of floats, `str.strip()`/`str.split()` white space and how pathlib
//! names files. The float printing and `sum()` follow PriVA's
//! `rust/src/python_text.rs`, where they were checked against Python.

/// `str.isspace()` for one character: Unicode white space plus U+001C to U+001F.
pub fn is_python_space(c: char) -> bool {
    matches!(
        c,
        '\t' | '\n'
            | '\u{b}'
            | '\u{c}'
            | '\r'
            | '\u{1c}'
            | '\u{1d}'
            | '\u{1e}'
            | '\u{1f}'
            | ' '
            | '\u{85}'
            | '\u{a0}'
            | '\u{1680}'
            | '\u{2000}'
            ..='\u{200a}' | '\u{2028}' | '\u{2029}' | '\u{202f}' | '\u{205f}' | '\u{3000}'
    )
}

/// `text.strip()`.
pub fn strip(text: &str) -> &str {
    text.trim_matches(is_python_space)
}

/// `text.split(maxsplit=1)` of a text that has no white space at either end:
/// the first word and the rest after the first run of white space.
pub fn split_first_word(text: &str) -> (&str, &str) {
    match text.find(is_python_space) {
        None => (text, ""),
        Some(position) => {
            let rest = text[position..].trim_start_matches(is_python_space);
            (&text[..position], rest)
        }
    }
}

/// The shortest digits that read back as `value` (positive and finite) and
/// the decimal exponent of the first digit.
fn shortest_digits(value: f64) -> (String, i32) {
    let scientific = format!("{value:e}");
    let (mantissa, exponent) = scientific
        .split_once('e')
        .expect("Rust prints floats in scientific form with an e");
    let exponent: i32 = exponent.parse().expect("the exponent is an integer");
    (mantissa.chars().filter(|c| *c != '.').collect(), exponent)
}

/// Python's `repr` takes the even digit when the value lies exactly halfway
/// between two shortest spellings; Rust's shortest form rounds such a tie up.
/// Only spellings of 15 or more digits can be affected.
fn tie_to_even(value: f64, digits: String, exponent: i32) -> (String, i32) {
    let count = digits.len();
    if count < 15 {
        return (digits, exponent);
    }
    let is_tie = |expansion: &str| {
        let (all, _) = expansion.split_once('e').unwrap_or((expansion, "0"));
        let all: String = all.chars().filter(|c| *c != '.').collect();
        all.len() > count
            && all.as_bytes()[count] == b'5'
            && all[count + 1..].bytes().all(|byte| byte == b'0')
    };
    if !is_tie(&format!("{:.*e}", count + 9, value)) || !is_tie(&format!("{value:.780e}")) {
        return (digits, exponent);
    }
    let expansion = format!("{value:.780e}");
    let (mantissa, exact_exponent) = expansion.split_once('e').expect("scientific form");
    let mut kept: Vec<u8> = mantissa
        .bytes()
        .filter(|byte| *byte != b'.')
        .take(count)
        .collect();
    let mut exponent: i32 = exact_exponent.parse().expect("the exponent is an integer");
    if (kept[count - 1] - b'0') % 2 == 1 {
        let mut index = count;
        loop {
            if index == 0 {
                kept.insert(0, b'1');
                kept.pop();
                exponent += 1;
                break;
            }
            index -= 1;
            if kept[index] == b'9' {
                kept[index] = b'0';
            } else {
                kept[index] += 1;
                break;
            }
        }
    }
    while kept.len() > 1 && kept.last() == Some(&b'0') {
        kept.pop();
    }
    let chosen = String::from_utf8(kept).expect("digits are ASCII");
    let spelled = format!("{}.{}e{exponent}", &chosen[..1], &chosen[1..]);
    if spelled.parse::<f64>() == Ok(value) {
        (chosen, exponent)
    } else {
        (digits, exponent)
    }
}

/// `repr(value)` (and `str(value)`) of a Python float: the shortest digits
/// that read back as the same number, in positional form for exponents from
/// -4 up to 15 and in scientific form (`1e-05`, `1.5e+16`) otherwise.
pub fn float_repr(value: f64) -> String {
    if value.is_nan() {
        return "nan".to_owned();
    }
    if value.is_infinite() {
        return (if value > 0.0 { "inf" } else { "-inf" }).to_owned();
    }
    if value == 0.0 {
        return (if value.is_sign_negative() {
            "-0.0"
        } else {
            "0.0"
        })
        .to_owned();
    }
    let (digits, exponent) = shortest_digits(value.abs());
    let (digits, exponent) = tie_to_even(value.abs(), digits, exponent);
    // The value is 0.DIGITS x 10^decimal_point.
    let decimal_point = exponent + 1;
    let count = digits.len() as i32;
    let mut text = String::new();
    if value < 0.0 {
        text.push('-');
    }
    if decimal_point <= -4 || decimal_point > 16 {
        text.push_str(&digits[..1]);
        if count > 1 {
            text.push('.');
            text.push_str(&digits[1..]);
        }
        let shown = decimal_point - 1;
        text.push('e');
        text.push(if shown < 0 { '-' } else { '+' });
        text.push_str(&format!("{:02}", shown.abs()));
    } else if decimal_point <= 0 {
        text.push_str("0.");
        for _ in 0..(-decimal_point) {
            text.push('0');
        }
        text.push_str(&digits);
    } else if decimal_point >= count {
        text.push_str(&digits);
        for _ in 0..(decimal_point - count) {
            text.push('0');
        }
        text.push_str(".0");
    } else {
        text.push_str(&digits[..decimal_point as usize]);
        text.push('.');
        text.push_str(&digits[decimal_point as usize..]);
    }
    text
}

/// `sum(values)` of Python floats in Python 3.12, which adds back the rounding
/// error of each step (Neumaier summation); an empty list sums to 0.
pub fn sum_floats(values: &[f64]) -> f64 {
    let mut iterator = values.iter();
    let Some(first) = iterator.next() else {
        return 0.0;
    };
    let mut total = *first;
    let mut compensation = 0.0_f64;
    for &item in iterator {
        let step = total + item;
        if total.abs() >= item.abs() {
            compensation += (total - step) + item;
        } else {
            compensation += (item - step) + total;
        }
        total = step;
    }
    if compensation != 0.0 && compensation.is_finite() {
        total += compensation;
    }
    total
}

/// `format(value, f".{digits}f")`: the decimal with `digits` places nearest to
/// the exact binary value, an exact tie going to the even last digit.
///
/// Rust's `format!` rounds every value that is not an exact tie correctly. An
/// exact tie needs `value * 10^digits` to end in `.5`; for a binary fraction
/// that happens only when `value * 2^(digits + 1)` is an odd whole number,
/// and those values are settled here.
pub fn fixed(value: f64, digits: usize) -> String {
    if value.is_nan() {
        return "nan".to_owned();
    }
    if value.is_infinite() {
        return (if value > 0.0 { "inf" } else { "-inf" }).to_owned();
    }
    let magnitude = value.abs();
    let scaled = magnitude * 2f64.powi(digits as i32 + 1);
    let body = if scaled.fract() == 0.0 && scaled < 9.0e15 && (scaled as u64) % 2 == 1 {
        // magnitude * 10^digits = odd * 5^digits / 2: take the even neighbour.
        let product = scaled as u128 * 5u128.pow(digits as u32);
        let lower = (product - 1) / 2;
        let even = if lower % 2 == 0 { lower } else { lower + 1 };
        let mut text = even.to_string();
        if text.len() <= digits {
            text = format!("{}{text}", "0".repeat(digits + 1 - text.len()));
        }
        let (whole, fraction) = text.split_at(text.len() - digits);
        if digits == 0 {
            whole.to_owned()
        } else {
            format!("{whole}.{fraction}")
        }
    } else {
        format!("{magnitude:.digits$}")
    };
    if value.is_sign_negative() {
        format!("-{body}")
    } else {
        body
    }
}

/// `round(value, digits)` of a Python float.
pub fn round(value: f64, digits: usize) -> f64 {
    if !value.is_finite() {
        return value;
    }
    fixed(value, digits)
        .parse()
        .expect("a formatted float reads back")
}

/// `count` significant digits of a positive finite value, rounded half to even
/// on the exact binary value, and the decimal exponent of the first digit.
fn significant_digits(magnitude: f64, count: usize) -> (Vec<u8>, i32) {
    // A 64-bit float has at most 767 significant digits, so this is exact.
    let exact = format!("{magnitude:.800e}");
    let (mantissa, exponent) = exact.split_once('e').expect("scientific form");
    let mut exponent: i32 = exponent.parse().expect("the exponent is an integer");
    let all: Vec<u8> = mantissa.bytes().filter(|byte| *byte != b'.').collect();
    let mut kept = all[..count].to_vec();
    let next = all[count];
    let rest_is_zero = all[count + 1..].iter().all(|byte| *byte == b'0');
    let round_up =
        next > b'5' || (next == b'5' && (!rest_is_zero || (kept[count - 1] - b'0') % 2 == 1));
    if round_up {
        let mut index = count;
        loop {
            if index == 0 {
                kept.insert(0, b'1');
                kept.pop();
                exponent += 1;
                break;
            }
            index -= 1;
            if kept[index] == b'9' {
                kept[index] = b'0';
            } else {
                kept[index] += 1;
                break;
            }
        }
    }
    (kept, exponent)
}

/// `format(value, "g")`: six significant digits, positional form for decimal
/// exponents from -4 to 5, trailing zeros (and a bare point) removed.
pub fn general(value: f64) -> String {
    if value.is_nan() {
        return "nan".to_owned();
    }
    if value.is_infinite() {
        return (if value > 0.0 { "inf" } else { "-inf" }).to_owned();
    }
    let sign = if value.is_sign_negative() { "-" } else { "" };
    if value == 0.0 {
        return format!("{sign}0");
    }
    let precision = 6usize;
    let (digits, exponent) = significant_digits(value.abs(), precision);
    let digits = String::from_utf8(digits).expect("digits are ASCII");
    let body = if (-4..precision as i32).contains(&exponent) {
        let text = if exponent >= 0 {
            let whole = exponent as usize + 1;
            format!("{}.{}", &digits[..whole], &digits[whole..])
        } else {
            format!("0.{}{}", "0".repeat((-exponent - 1) as usize), digits)
        };
        let text = text.trim_end_matches('0');
        text.trim_end_matches('.').to_owned()
    } else {
        let rest = digits[1..].trim_end_matches('0');
        let mut text = digits[..1].to_owned();
        if !rest.is_empty() {
            text.push('.');
            text.push_str(rest);
        }
        format!(
            "{text}e{}{:02}",
            if exponent < 0 { '-' } else { '+' },
            exponent.abs()
        )
    };
    format!("{sign}{body}")
}

/// `str(PurePosixPath(text))`: repeated slashes and `.` parts dropped, no
/// trailing slash, `.` for an empty path, two leading slashes kept.
pub fn path_text(text: &str) -> String {
    let (root, rest) = if text.starts_with("//") && !text.starts_with("///") {
        ("//", &text[2..])
    } else if text.starts_with('/') {
        ("/", text.trim_start_matches('/'))
    } else {
        ("", text)
    };
    let parts: Vec<&str> = rest
        .split('/')
        .filter(|part| !part.is_empty() && *part != ".")
        .collect();
    let joined = parts.join("/");
    if root.is_empty() && joined.is_empty() {
        ".".to_owned()
    } else {
        format!("{root}{joined}")
    }
}

/// `str(PurePosixPath(directory) / name)`.
pub fn path_join(directory: &str, name: &str) -> String {
    let base = path_text(directory);
    if base == "." {
        path_text(name)
    } else if base.ends_with('/') {
        path_text(&format!("{base}{name}"))
    } else {
        path_text(&format!("{base}/{name}"))
    }
}

/// `PurePosixPath(path).with_suffix(suffix)`: the last `.part` of the final
/// name (when the dot is neither its first nor its last character) is
/// replaced, otherwise the suffix is added.
pub fn with_suffix(path: &str, suffix: &str) -> String {
    let normalized = path_text(path);
    let (directory, name) = match normalized.rfind('/') {
        Some(position) => (&normalized[..=position], &normalized[position + 1..]),
        None => ("", normalized.as_str()),
    };
    let stem = match name.rfind('.') {
        Some(position) if position > 0 && position < name.len() - 1 => &name[..position],
        _ => name,
    };
    format!("{directory}{stem}{suffix}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn float_repr_matches_python() {
        assert_eq!(float_repr(0.1), "0.1");
        assert_eq!(float_repr(1.0), "1.0");
        assert_eq!(float_repr(-3.0), "-3.0");
        assert_eq!(float_repr(0.0001), "0.0001");
        assert_eq!(float_repr(0.00001), "1e-05");
        assert_eq!(float_repr(1.5e-7), "1.5e-07");
        assert_eq!(float_repr(1e15), "1000000000000000.0");
        assert_eq!(float_repr(1e16), "1e+16");
        assert_eq!(float_repr(0.30000000000000004), "0.30000000000000004");
        assert_eq!(float_repr(129.142852783203125), "129.14285278320312");
        assert_eq!(float_repr(-0.0), "-0.0");
    }

    #[test]
    fn fixed_and_round_take_the_even_digit_on_exact_ties() {
        // 1/1024 = 0.0009765625 is an exact tie at nine places.
        assert_eq!(fixed(0.0009765625, 9), "0.000976562");
        assert_eq!(round(0.0009765625, 9), 0.000976562);
        // 3/1024 = 0.0029296875: the kept digit 7 is odd, so it goes up.
        assert_eq!(fixed(0.0029296875, 9), "0.002929688");
        assert_eq!(fixed(0.125, 2), "0.12");
        assert_eq!(fixed(0.375, 2), "0.38");
        assert_eq!(fixed(2.5, 0), "2");
        assert_eq!(fixed(1.0 / 3.0, 6), "0.333333");
        assert_eq!(fixed(2.675, 2), "2.67");
        assert_eq!(fixed(-0.0000001, 6), "-0.000000");
        assert_eq!(round(1.0 / 3.0, 6), 0.333333);
        assert_eq!(float_repr(round(0.1 + 0.2, 9)), "0.3");
    }

    #[test]
    fn general_matches_python_g() {
        assert_eq!(general(3.0), "3");
        assert_eq!(general(12.0), "12");
        assert_eq!(general(2.5), "2.5");
        assert_eq!(general(0.0), "0");
        assert_eq!(general(1e-5), "1e-05");
        assert_eq!(general(0.0001), "0.0001");
        assert_eq!(general(123456.0), "123456");
        assert_eq!(general(1234567.0), "1.23457e+06");
        assert_eq!(general(999999.5), "1e+06");
        assert_eq!(general(0.1015625), "0.101562");
        assert_eq!(general(0.1171875), "0.117188");
    }

    #[test]
    fn sum_adds_back_rounding_error_like_python_3_12() {
        assert_eq!(sum_floats(&[]), 0.0);
        assert_eq!(sum_floats(&[0.1, 0.2, 0.3]), 0.6);
        assert_eq!(sum_floats(&[0.1; 10]), 1.0);
    }

    #[test]
    fn paths_follow_pathlib() {
        assert_eq!(path_text("a//b/./c/"), "a/b/c");
        assert_eq!(path_text("./a"), "a");
        assert_eq!(path_text(""), ".");
        assert_eq!(path_join("data/raw/", "x.txt"), "data/raw/x.txt");
        assert_eq!(path_join(".", "x.txt"), "x.txt");
        assert_eq!(
            with_suffix("out/PID18-023", ".gene_ranked.tsv"),
            "out/PID18-023.gene_ranked.tsv"
        );
        assert_eq!(
            with_suffix("out/S.4", ".gene_ranked.tsv"),
            "out/S.gene_ranked.tsv"
        );
        assert_eq!(with_suffix("out/.hidden", ".tsv"), "out/.hidden.tsv");
        assert_eq!(with_suffix("out/name.", ".tsv"), "out/name..tsv");
    }

    #[test]
    fn white_space_follows_python() {
        assert_eq!(strip("\u{a0} GO:1 name \n"), "GO:1 name");
        assert_eq!(split_first_word("GO:1  two words"), ("GO:1", "two words"));
        assert_eq!(split_first_word("GO:1"), ("GO:1", ""));
    }
}
