//! Python text rules that converted scripts depend on: `float()` and `int()` of
//! a string, `repr()` of a float and `sum()` of floats.
//!
//! Python 3.12 (the version of PriVA's conda environment) adds the rounding
//! error of each step back when `sum()` adds floats: `sum([0.1] * 10)` is
//! `1.0` there but `0.9999999999999999` for a plain loop, so a comparison such
//! as `total > 0.3` can differ between the two. `sum_floats` reproduces the
//! compensated sum.
//!
//! Known gap: Python also reads non-ASCII decimal digits (for example
//! full-width digits) in `int()` and `float()`; these readers do not.
use crate::acmg_scoring::python_strip;

/// Removes the underscores Python allows between two digits (`1_000`); any
/// other underscore makes the text invalid.
fn without_underscores(text: &str) -> Option<String> {
    let chars: Vec<char> = text.chars().collect();
    let mut cleaned = String::with_capacity(text.len());
    for (index, &character) in chars.iter().enumerate() {
        if character == '_' {
            let before = index.checked_sub(1).map(|position| chars[position]);
            let after = chars.get(index + 1).copied();
            if !before.is_some_and(|c| c.is_ascii_digit())
                || !after.is_some_and(|c| c.is_ascii_digit())
            {
                return None;
            }
        } else {
            cleaned.push(character);
        }
    }
    Some(cleaned)
}

/// `int(text)` for a decimal integer; `None` when Python raises ValueError
/// (or when the value does not fit in 64 bits).
pub(crate) fn parse_int(text: &str) -> Option<i64> {
    let text = python_strip(text);
    let (negative, rest) = match text.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, text.strip_prefix('+').unwrap_or(text)),
    };
    let digits = without_underscores(rest)?;
    if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
        return None;
    }
    let magnitude: i128 = digits.parse().ok()?;
    i64::try_from(if negative { -magnitude } else { magnitude }).ok()
}

/// `float(text)`; `None` when Python raises ValueError. Accepts `inf`,
/// `infinity` and `nan` in any case with an optional sign.
pub(crate) fn parse_float(text: &str) -> Option<f64> {
    let text = python_strip(text);
    if text.contains('_') {
        without_underscores(text)?.parse::<f64>().ok()
    } else {
        text.parse::<f64>().ok()
    }
}

/// `repr(value)` of a Python float: the shortest digits that read back as the
/// same number, in positional form for exponents from -4 up to 15 and in
/// scientific form (`1e-05`, `1.5e+16`) otherwise.
pub(crate) fn float_repr(value: f64) -> String {
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
    let scientific = format!("{:e}", value.abs());
    let (mantissa, exponent) = scientific
        .split_once('e')
        .expect("Rust prints floats in scientific form with an e");
    let exponent: i32 = exponent.parse().expect("the exponent is an integer");
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
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

/// `sum(values)` of Python floats in Python 3.12 (Neumaier compensated
/// summation); an empty list sums to 0.
pub(crate) fn sum_floats(values: &[f64]) -> f64 {
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
    // Do not let the compensation turn an infinite or overflowed sum into NaN.
    if compensation != 0.0 && compensation.is_finite() {
        total += compensation;
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integers_follow_python_int() {
        assert_eq!(parse_int(" 12 "), Some(12));
        assert_eq!(parse_int("+7"), Some(7));
        assert_eq!(parse_int("-0"), Some(0));
        assert_eq!(parse_int("007"), Some(7));
        assert_eq!(parse_int("1_000"), Some(1000));
        assert_eq!(parse_int("1__0"), None);
        assert_eq!(parse_int("_1"), None);
        assert_eq!(parse_int("1_"), None);
        assert_eq!(parse_int(""), None);
        assert_eq!(parse_int("-"), None);
        assert_eq!(parse_int("1.0"), None);
        assert_eq!(parse_int("3-4"), None);
        assert_eq!(parse_int("\u{1c}5\u{1f}"), Some(5));
    }

    #[test]
    fn floats_follow_python_float() {
        assert_eq!(parse_float("1e3"), Some(1000.0));
        assert_eq!(parse_float(" .5 "), Some(0.5));
        assert_eq!(parse_float("1_0.5"), Some(10.5));
        assert_eq!(parse_float("inf"), Some(f64::INFINITY));
        assert_eq!(parse_float("-Infinity"), Some(f64::NEG_INFINITY));
        assert!(parse_float("NaN").unwrap().is_nan());
        assert_eq!(parse_float("1,5"), None);
        assert_eq!(parse_float("1e"), None);
        assert_eq!(parse_float("_1"), None);
        assert_eq!(parse_float(""), None);
        assert_eq!(parse_float("0x10"), None);
    }

    #[test]
    fn float_repr_matches_python() {
        assert_eq!(float_repr(0.1), "0.1");
        assert_eq!(float_repr(5.0), "5.0");
        assert_eq!(float_repr(-3.0), "-3.0");
        assert_eq!(float_repr(100.0), "100.0");
        assert_eq!(float_repr(12345.678), "12345.678");
        assert_eq!(float_repr(0.0001), "0.0001");
        assert_eq!(float_repr(0.00001), "1e-05");
        assert_eq!(float_repr(1.5e-7), "1.5e-07");
        assert_eq!(float_repr(1e15), "1000000000000000.0");
        assert_eq!(float_repr(1e16), "1e+16");
        assert_eq!(float_repr(1.2345678901234568e17), "1.2345678901234568e+17");
        assert_eq!(float_repr(1e100), "1e+100");
        assert_eq!(float_repr(0.30000000000000004), "0.30000000000000004");
        assert_eq!(float_repr(f64::NAN), "nan");
        assert_eq!(float_repr(f64::NEG_INFINITY), "-inf");
        assert_eq!(float_repr(-0.0), "-0.0");
    }

    #[test]
    fn sum_adds_back_rounding_error_like_python_3_12() {
        // In Python 3.12 sum([0.1, 0.2, 0.3]) is 0.6, where a plain loop gives
        // 0.6000000000000001.
        assert_eq!(sum_floats(&[]), 0.0);
        assert_eq!(sum_floats(&[0.1, 0.2, 0.3]), 0.6);
        assert_eq!(sum_floats(&[0.1; 10]), 1.0);
        assert_eq!(sum_floats(&[0.02, 0.006]), 0.02 + 0.006);
    }
}
