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
//! full-width digits) in `int()` and `float()`; these readers do not, except
//! `python_int_is_valid`, which only tells whether `int()` accepts a text.
pub(crate) use crate::acmg_scoring::python_strip;

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
/// between two shortest spellings (`129.142852783203125`, which a 32-bit float
/// can hold, prints as `129.14285278320312`); Rust's shortest form rounds
/// such a tie up. Only spellings of 15 or more digits can be affected.
fn tie_to_even(value: f64, digits: String, exponent: i32) -> (String, i32) {
    let count = digits.len();
    if count < 15 {
        return (digits, exponent);
    }
    // The digits after the shortest ones: a tie is a 5 followed by nothing but zeros.
    let is_tie = |expansion: &str| {
        let (all, _) = expansion.split_once('e').unwrap_or((expansion, "0"));
        let all: String = all.chars().filter(|c| *c != '.').collect();
        all.len() > count
            && all.as_bytes()[count] == b'5'
            && all[count + 1..].bytes().all(|byte| byte == b'0')
    };
    // A cheap look at a few more digits first; the full expansion (a 64-bit
    // float has at most 767 significant digits) only when that looks like a tie.
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
    // Both neighbours read back as the value for a true tie; keep Rust's digits otherwise.
    if spelled.parse::<f64>() == Ok(value) {
        (chosen, exponent)
    } else {
        (digits, exponent)
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

/// `str.isspace()` for one character: Unicode white space plus U+001C to U+001F.
pub(crate) fn is_python_space(c: char) -> bool {
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

/// `int(text)` for a base-10 string: surrounding white space, an optional
/// sign, ASCII digits with single underscores between them. Digits of other
/// scripts, which Python also accepts, are not supported and give `None`.
pub(crate) fn python_int(text: &str) -> Option<i128> {
    let body = python_strip(text);
    let (negative, digits) = match body.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, body.strip_prefix('+').unwrap_or(body)),
    };
    if digits.is_empty()
        || digits.starts_with('_')
        || digits.ends_with('_')
        || digits.contains("__")
    {
        return None;
    }
    if !digits
        .bytes()
        .all(|byte| byte.is_ascii_digit() || byte == b'_')
    {
        return None;
    }
    let cleaned: String = digits.chars().filter(|c| *c != '_').collect();
    let value: i128 = cleaned.parse().ok()?;
    Some(if negative { -value } else { value })
}

/// The first character (the digit zero) of each run of ten characters that
/// Unicode lists as decimal digits (category Nd), in the Unicode 15.0 data of
/// Python 3.12; read from `unicodedata.decimal` of the interpreter PriVA runs on.
const DECIMAL_DIGIT_ZEROS: [u32; 68] = [
    0x30, 0x660, 0x6f0, 0x7c0, 0x966, 0x9e6, 0xa66, 0xae6, 0xb66, 0xbe6, 0xc66, 0xce6, 0xd66,
    0xde6, 0xe50, 0xed0, 0xf20, 0x1040, 0x1090, 0x17e0, 0x1810, 0x1946, 0x19d0, 0x1a80, 0x1a90,
    0x1b50, 0x1bb0, 0x1c40, 0x1c50, 0xa620, 0xa8d0, 0xa900, 0xa9d0, 0xa9f0, 0xaa50, 0xabf0, 0xff10,
    0x104a0, 0x10d30, 0x11066, 0x110f0, 0x11136, 0x111d0, 0x112f0, 0x11450, 0x114d0, 0x11650,
    0x116c0, 0x11730, 0x118e0, 0x11950, 0x11c50, 0x11d50, 0x11da0, 0x11f50, 0x16a60, 0x16ac0,
    0x16b50, 0x1d7ce, 0x1d7d8, 0x1d7e2, 0x1d7ec, 0x1d7f6, 0x1e140, 0x1e2f0, 0x1e4f0, 0x1e950,
    0x1fbf0,
];

/// The value of a decimal digit of any script, as Python's `int()` reads it.
pub(crate) fn decimal_digit_value(c: char) -> Option<u32> {
    let code = c as u32;
    let after = DECIMAL_DIGIT_ZEROS.partition_point(|zero| *zero <= code);
    let offset = code - *DECIMAL_DIGIT_ZEROS.get(after.checked_sub(1)?)?;
    (offset < 10).then_some(offset)
}

/// Python 3.12 refuses to convert a string of more than this many digits.
const MAX_INT_STRING_DIGITS: usize = 4300;

/// Whether `int(text)` succeeds for a base-10 string of any length: white
/// space around it, an optional sign, and digits of any script with single
/// underscores between them. Python first turns every non-ASCII white space
/// into a blank and every non-ASCII decimal digit into its ASCII digit; the
/// ASCII text is then read with C white space (so U+001C to U+001F, which
/// `str.strip` removes, are not blanks here).
pub(crate) fn python_int_is_valid(text: &str) -> bool {
    let converted: String = text
        .chars()
        .map(|c| {
            if c.is_ascii() {
                c
            } else if let Some(value) = decimal_digit_value(c) {
                char::from(b'0' + value as u8)
            } else if is_python_space(c) {
                ' '
            } else {
                c
            }
        })
        .collect();
    let body =
        converted.trim_matches(|c| matches!(c, ' ' | '\t' | '\n' | '\u{b}' | '\u{c}' | '\r'));
    let digits = body
        .strip_prefix('-')
        .or_else(|| body.strip_prefix('+'))
        .unwrap_or(body);
    if digits.is_empty()
        || digits.starts_with('_')
        || digits.ends_with('_')
        || digits.contains("__")
        || !digits
            .bytes()
            .all(|byte| byte.is_ascii_digit() || byte == b'_')
    {
        return false;
    }
    digits.bytes().filter(u8::is_ascii_digit).count() <= MAX_INT_STRING_DIGITS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn digits_of_other_scripts_are_read_by_int() {
        assert_eq!(decimal_digit_value('7'), Some(7));
        assert_eq!(decimal_digit_value('\u{661}'), Some(1));
        assert_eq!(decimal_digit_value('\u{669}'), Some(9));
        assert_eq!(decimal_digit_value('\u{66a}'), None);
        assert_eq!(decimal_digit_value('\u{ff13}'), Some(3));
        assert_eq!(decimal_digit_value('\u{1d7ff}'), Some(9));
        assert_eq!(decimal_digit_value('\u{1fbf9}'), Some(9));
        assert_eq!(
            decimal_digit_value('\u{b2}'),
            None,
            "superscript two is a digit, not a decimal digit"
        );
        assert_eq!(decimal_digit_value('\u{bd}'), None);
        assert_eq!(decimal_digit_value('a'), None);
        assert_eq!(decimal_digit_value('\u{10ffff}'), None);
        // Every run holds ten consecutive digits and the runs are sorted.
        assert!(
            DECIMAL_DIGIT_ZEROS
                .windows(2)
                .all(|pair| pair[1] >= pair[0] + 10)
        );
        assert!(python_int_is_valid(" \u{661}\u{662} "));
        assert!(python_int_is_valid("-\u{e53}_\u{ff11}"));
        assert!(python_int_is_valid("\u{a0}+12\u{2003}"));
        assert!(!python_int_is_valid("\u{661}x"));
        assert!(!python_int_is_valid("\u{b2}"));
        assert!(!python_int_is_valid("1\u{66b}5"));
    }

    #[test]
    fn int_text_of_any_size_up_to_the_digit_limit_is_valid() {
        for (text, valid) in [
            ("12", true),
            (" -7 ", true),
            ("1_0", true),
            ("1__0", false),
            ("_1", false),
            ("1_", false),
            ("", false),
            ("-", false),
            ("+ 1", false),
            ("1.0", false),
            ("1e3", false),
            ("0x10", false),
            ("1 2", false),
            ("\u{1c}5\u{1f}", false),
            ("\t5\u{b}\u{c}\r\n", true),
        ] {
            assert_eq!(python_int_is_valid(text), valid, "{text:?}");
        }
        assert!(python_int_is_valid(&"9".repeat(4300)));
        assert!(!python_int_is_valid(&"9".repeat(4301)));
        // Underscores are not digits.
        assert!(python_int_is_valid(&format!("{}1", "1_".repeat(2150))));
    }

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
    fn float_repr_takes_the_even_digit_on_an_exact_tie() {
        // Values that a 32-bit float holds exactly, so that the 18th digit is a 5.
        assert_eq!(float_repr(129.142852783203125), "129.14285278320312");
        assert_eq!(float_repr(-129.142852783203125), "-129.14285278320312");
        assert_eq!(float_repr(100.047149658203125), "100.04714965820312");
        assert_eq!(float_repr(100.106781005859375), "100.10678100585938");
        assert_eq!(float_repr(100.166412353515625), "100.16641235351562");
        // Values that are no tie keep their shortest digits.
        assert_eq!(float_repr(0.1f32 as f64), "0.10000000149011612");
        assert_eq!(float_repr(0.30000000000000004), "0.30000000000000004");
        assert_eq!(float_repr(2.675), "2.675");
        assert_eq!(float_repr(1.0 / 3.0), "0.3333333333333333");
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

    #[test]
    fn strip_uses_python_white_space() {
        assert_eq!(python_strip("\u{1c} a b\u{a0}\n"), "a b");
        assert_eq!(python_strip("   "), "");
    }

    #[test]
    fn int_follows_python_rules() {
        assert_eq!(python_int(" 12 "), Some(12));
        assert_eq!(python_int("-0_7"), Some(-7));
        assert_eq!(python_int("+5"), Some(5));
        assert_eq!(python_int("1__0"), None);
        assert_eq!(python_int("_1"), None);
        assert_eq!(python_int("1_"), None);
        assert_eq!(python_int("1.0"), None);
        assert_eq!(python_int(""), None);
        assert_eq!(python_int("0x10"), None);
        assert_eq!(python_int("007"), Some(7));
    }
}
