//! Long-option parsing that follows Python's `argparse` for the converted scripts.
//!
//! Rules reproduced from argparse:
//! - `--name value` and `--name=value` both work, and a unique prefix of a
//!   long option stands for it (`--var` for `--variants_table`); an exact name
//!   wins over longer names that start with it;
//! - a value may not look like an option: it starts with `-`, is not just `-`,
//!   is not a negative number and holds no space;
//! - a repeated option keeps its last value;
//! - arguments that no option takes are reported together at the end;
//! - `-h` or `--help` stops the parsing and asks for the usage text.
//! The caller checks which options are required and exits with code 2 on an error.
use std::collections::HashMap;

/// What the command line asked for.
pub(crate) enum Parsed {
    Help,
    Values(HashMap<String, String>),
}

/// `-5`, `-5.5` and `-.5` are values, not options (argparse's negative number rule).
fn is_negative_number(token: &str) -> bool {
    let body = match token.strip_prefix('-') {
        Some(body) => body,
        None => return false,
    };
    if body.is_empty() {
        return false;
    }
    if body.chars().all(|c| c.is_ascii_digit()) {
        return true;
    }
    match body.split_once('.') {
        Some((whole, fraction)) => {
            whole.chars().all(|c| c.is_ascii_digit())
                && !fraction.is_empty()
                && fraction.chars().all(|c| c.is_ascii_digit())
        }
        None => false,
    }
}

fn looks_like_option(token: &str) -> bool {
    token.starts_with('-') && token != "-" && !is_negative_number(token) && !token.contains(' ')
}

/// The option an argument names: the exact name, else the only name it prefixes.
fn match_option<'a>(given: &str, names: &[&'a str]) -> Result<Option<&'a str>, String> {
    let bare = given.trim_start_matches('-');
    if let Some(name) = names.iter().find(|name| **name == bare) {
        return Ok(Some(*name));
    }
    let candidates: Vec<&'a str> = names
        .iter()
        .filter(|name| name.starts_with(bare))
        .copied()
        .collect();
    match candidates.len() {
        0 => Ok(None),
        1 => Ok(Some(candidates[0])),
        _ => Err(format!(
            "ambiguous option: {given} could match {}",
            candidates
                .iter()
                .map(|name| format!("--{name}"))
                .collect::<Vec<_>>()
                .join(", ")
        )),
    }
}

/// Parses `args` against the long option `names` (written without the dashes).
pub(crate) fn parse_long_options(args: &[String], names: &[&str]) -> Result<Parsed, String> {
    let mut values: HashMap<String, String> = HashMap::new();
    let mut unrecognized: Vec<String> = Vec::new();
    let mut index = 0;
    while index < args.len() {
        let token = &args[index];
        index += 1;
        if token == "-h" || token == "--help" {
            return Ok(Parsed::Help);
        }
        if !token.starts_with("--") || token == "--" {
            unrecognized.push(token.clone());
            continue;
        }
        let (given, inline) = match token.split_once('=') {
            Some((name, value)) => (name.to_owned(), Some(value.to_owned())),
            None => (token.clone(), None),
        };
        let name = match match_option(&given, names)? {
            Some(name) => name,
            None => {
                unrecognized.push(token.clone());
                continue;
            }
        };
        let value = match inline {
            Some(value) => value,
            None => match args.get(index) {
                Some(next) if !looks_like_option(next) => {
                    index += 1;
                    next.clone()
                }
                _ => return Err(format!("argument --{name}: expected one argument")),
            },
        };
        values.insert(name.to_owned(), value);
    }
    if !unrecognized.is_empty() {
        return Err(format!("unrecognized arguments: {}", unrecognized.join(" ")));
    }
    Ok(Parsed::Values(values))
}

/// The error for required options that are absent, in the order given.
pub(crate) fn require_options(
    values: &HashMap<String, String>,
    required: &[&str],
) -> Result<(), String> {
    let missing: Vec<String> = required
        .iter()
        .filter(|name| !values.contains_key(**name))
        .map(|name| format!("--{name}"))
        .collect();
    if missing.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "the following arguments are required: {}",
            missing.join(", ")
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(items: &[&str]) -> Vec<String> {
        items.iter().map(|item| item.to_string()).collect()
    }

    fn values(items: &[&str], names: &[&str]) -> Result<HashMap<String, String>, String> {
        match parse_long_options(&strings(items), names)? {
            Parsed::Values(values) => Ok(values),
            Parsed::Help => Err("help".to_owned()),
        }
    }

    #[test]
    fn options_take_a_value_in_both_spellings_and_by_prefix() {
        let names = ["variants_table", "domain_map"];
        let parsed = values(&["--var", "a.tsv", "--domain_map=m.pkl"], &names).unwrap();
        assert_eq!(parsed["variants_table"], "a.tsv");
        assert_eq!(parsed["domain_map"], "m.pkl");
        let repeated = values(&["--domain_map", "x", "--domain_map", "y"], &names).unwrap();
        assert_eq!(repeated["domain_map"], "y");
    }

    #[test]
    fn an_exact_name_beats_a_longer_name_that_starts_with_it() {
        let names = ["snv", "snv_extra"];
        assert_eq!(values(&["--snv", "a"], &names).unwrap()["snv"], "a");
        assert!(values(&["--sn", "a"], &names).is_err());
    }

    #[test]
    fn a_missing_value_or_an_option_as_value_is_an_error() {
        let names = ["a", "b"];
        assert!(values(&["--a"], &names).is_err());
        assert!(values(&["--a", "--b", "x"], &names).is_err());
        assert_eq!(values(&["--a", "-5"], &names).unwrap()["a"], "-5");
        assert_eq!(values(&["--a", "x y"], &names).unwrap()["a"], "x y");
        assert_eq!(values(&["--a", ""], &names).unwrap()["a"], "");
    }

    #[test]
    fn unknown_arguments_are_reported_together() {
        let error = values(&["--a", "1", "stray", "--zzz"], &["a"]).unwrap_err();
        assert_eq!(error, "unrecognized arguments: stray --zzz");
    }

    #[test]
    fn help_stops_parsing_and_required_options_are_listed() {
        assert!(matches!(
            parse_long_options(&strings(&["--a", "1", "-h"]), &["a"]),
            Ok(Parsed::Help)
        ));
        let mut present = HashMap::new();
        present.insert("a".to_owned(), "1".to_owned());
        assert_eq!(
            require_options(&present, &["a", "b", "c"]).unwrap_err(),
            "the following arguments are required: --b, --c"
        );
    }
}
