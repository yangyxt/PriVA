//! Read-only resource lifecycle inventory; no cache payload is opened.
//! Preserve the Python registry's object order and validation/traversal order.
use crate::cli::{Result, invalid};
use crate::json_text::{self, Json};
use crate::python_args::{Parsed, parse_long_options_with_flags, require_options};
use crate::python_set_order::StringSet;
use crate::python_text::is_python_space;
use std::collections::HashSet;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

const HELP: &str = "usage: priva cache-resources {list|show|refresh-plan|validate} --registry JSON [--resource NAME] [--refreshable] [--check-routes]\n\
Read only the selected resource registry and, for --check-routes, its adjacent scripts/install_utils.sh.\n\
Use an explicit PYTHONHASHSEED (0..4294967295) to preserve Python's dependency/route set order.\n\
The registry path is explicit; the binary's install location never selects a checkout.\n\
Existing Python callers remain unchanged. A refresh plan is an inventory, never execution authorization.";
const PLAN_NOTE: &str = "Review affected resources for the selected assembly/profile. This list is not an execution order or authorization to refresh.";

fn error(category: &'static str, message: impl Into<String>) -> Box<dyn std::error::Error> {
    crate::builder_error::category(category, invalid(message))
}
fn require(condition: bool, message: impl Into<String>) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(error("ValueError", message))
    }
}
fn kind(value: &Json) -> &'static str {
    match value {
        Json::Null => "NoneType",
        Json::Bool(_) => "bool",
        Json::Int(_) => "int",
        Json::Float(_) => "float",
        Json::Str(_) => "str",
        Json::Array(_) => "list",
        Json::Object(_) => "dict",
    }
}
fn get<'a>(value: &'a Json, key: &str, default: &'a Json) -> Result<&'a Json> {
    if value.as_object().is_none() {
        return Err(error(
            "AttributeError",
            format!("'{}' object has no attribute 'get'", kind(value)),
        ));
    }
    Ok(value.get(key).unwrap_or(default))
}
fn strings(value: &Json) -> Option<Vec<String>> {
    match value {
        Json::Array(values) => values
            .iter()
            .map(|v| v.as_str().map(str::to_owned))
            .collect(),
        _ => None,
    }
}
fn object(items: impl IntoIterator<Item = (&'static str, Json)>) -> Json {
    Json::Object(items.into_iter().map(|(k, v)| (k.to_owned(), v)).collect())
}
fn id_valid(text: &str) -> bool {
    let mut bytes = text.bytes();
    bytes.next().is_some_and(|c| c.is_ascii_lowercase())
        && bytes.all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == b'_')
}

/// Only membership is needed from visiting/visited/affected sets. Iterated
/// dependencies are constructed exactly as set(list), then update(each list).
fn dependencies(entry: &Json, seed: u32) -> Vec<String> {
    let mut result = StringSet::new(seed);
    for value in strings(entry.get("depends_on").unwrap()).unwrap() {
        result.insert(value);
    }
    if let Some(assemblies) = entry.get("assembly_dependencies") {
        for (_, values) in assemblies.as_object().unwrap() {
            for value in strings(values).unwrap() {
                result.insert(value);
            }
        }
    }
    result.into_order()
}

fn validate_registry(registry: &Json, seed: u32) -> Result<()> {
    require(
        json_text::py_equal(get(registry, "schema_version", &Json::Null)?, &Json::Int(1)),
        "Unsupported cache resource schema_version",
    )?;
    let resources = get(registry, "resources", &Json::Null)?;
    let entries = resources.as_object();
    require(
        entries.is_some_and(|v| !v.is_empty()),
        "The resource map must be a nonempty object",
    )?;
    let entries = entries.unwrap();
    let empty = Json::Object(Vec::new());
    for (name, entry) in entries {
        require(id_valid(name), format!("Invalid resource ID: {name}"))?;
        let lifecycle = get(entry, "lifecycle", &Json::Null)?;
        if matches!(lifecycle, Json::Array(_) | Json::Object(_)) {
            return Err(error(
                "TypeError",
                format!("unhashable type: '{}'", kind(lifecycle)),
            ));
        }
        require(
            lifecycle
                .as_str()
                .is_some_and(|s| matches!(s, "maintained" | "snapshot" | "derived")),
            format!("Invalid lifecycle for {name}"),
        )?;
        for key in ["config_keys", "bundled_paths", "depends_on"] {
            require(
                strings(get(entry, key, &Json::Null)?).is_some(),
                format!("{name}.{key} must be a list of strings"),
            )?;
        }
        let assemblies = get(entry, "assembly_dependencies", &empty)?;
        require(
            assemblies
                .as_object()
                .is_some_and(|pairs| pairs.iter().all(|(_, values)| strings(values).is_some())),
            format!("Invalid assembly dependencies for {name}"),
        )?;
        let refresh = get(entry, "refresh", &empty)?;
        for key in ["automatic", "keep_function"] {
            require(
                matches!(get(refresh, key, &Json::Null)?, Json::Bool(_)),
                format!("{name}.refresh.{key} must be boolean"),
            )?;
        }
        let functions = strings(get(refresh, "entrypoints", &Json::Null)?);
        require(
            functions.is_some(),
            format!("Invalid refresh entrypoints for {name}"),
        )?;
        require(
            json_text::py_equal(
                refresh.get("keep_function").unwrap(),
                &Json::Bool(!functions.unwrap().is_empty()),
            ),
            format!("Refresh-function inventory is inconsistent for {name}"),
        )?;
        for dependency in dependencies(entry, seed) {
            require(
                resources.get(&dependency).is_some(),
                format!("Unknown dependency {dependency} of {name}"),
            )?;
        }
    }
    let mut visiting = HashSet::new();
    let mut visited = HashSet::new();
    // An explicit stack retains recursive DFS's enter/parent/leave sequence
    // without risking a native stack overflow on a large supplied registry.
    for (name, _) in entries {
        let mut pending = vec![(name.clone(), false)];
        while let Some((name, leave)) = pending.pop() {
            if leave {
                visiting.remove(&name);
                visited.insert(name);
                continue;
            }
            require(
                !visiting.contains(&name),
                format!("Resource dependency cycle at {name}"),
            )?;
            if visited.contains(&name) {
                continue;
            }
            visiting.insert(name.clone());
            pending.push((name.clone(), true));
            for parent in dependencies(resources.get(&name).unwrap(), seed)
                .into_iter()
                .rev()
            {
                pending.push((parent, false));
            }
        }
    }
    Ok(())
}

/// The source regex is literal except for Python whitespace and its ASCII
/// identifier. Whitespace may span lines; only the initial anchor is multiline.
fn shell_functions(shell: &str) -> HashSet<String> {
    let mut names = HashSet::new();
    for start in std::iter::once(0).chain(shell.match_indices('\n').map(|(i, _)| i + 1)) {
        let Some(after) = shell[start..].strip_prefix("function") else {
            continue;
        };
        let tail = after.trim_start_matches(is_python_space);
        if tail.len() == after.len() {
            continue;
        }
        let mut chars = tail.char_indices();
        if !chars
            .next()
            .is_some_and(|(_, c)| c.is_ascii_alphabetic() || c == '_')
        {
            continue;
        }
        let end = chars
            .find(|(_, c)| !c.is_ascii_alphanumeric() && *c != '_')
            .map_or(tail.len(), |(i, _)| i);
        if tail[end..]
            .trim_start_matches(is_python_space)
            .starts_with('(')
        {
            names.insert(tail[..end].to_owned());
        }
    }
    names
}

fn check_routes(registry: &Json, source_root: &Path, seed: u32) -> Result<()> {
    let shell = read_text(&source_root.join("scripts/install_utils.sh"))?;
    let functions = shell_functions(&shell);
    for (name, entry) in registry.get("resources").unwrap().as_object().unwrap() {
        let mut names = strings(entry.get("refresh").unwrap().get("entrypoints").unwrap()).unwrap();
        let mut nonstrings = Vec::new();
        if let Some(install) = entry.get("install_entrypoints") {
            match install {
                Json::Array(values) => {
                    // set(list) consumes the complete iterable before any
                    // missing-route check. Its first unhashable member wins
                    // even when an earlier string names a missing function.
                    for value in values {
                        if matches!(value, Json::Array(_) | Json::Object(_)) {
                            return Err(error(
                                "TypeError",
                                format!("unhashable type: '{}'", kind(value)),
                            ));
                        }
                        if let Some(value) = value.as_str() {
                            names.push(value.to_owned());
                        } else {
                            nonstrings.push(value);
                        }
                    }
                }
                value => {
                    return Err(error(
                        "TypeError",
                        format!(
                            "can only concatenate list (not \"{}\") to list",
                            kind(value)
                        ),
                    ));
                }
            }
        }
        if !nonstrings.is_empty() {
            // A lone non-string is necessarily absent from a string function
            // catalogue. If all string members exist, the deciding error is
            // independent of mixed-object set hashing/iteration. Multiple
            // candidates remain outside the proven string-set profile.
            if nonstrings.len() != 1 || names.iter().any(|name| !functions.contains(name)) {
                return Err(error(
                    "RuntimeError",
                    "Unsupported profile: mixed route sets with multiple missing candidates require Python mixed-object set order",
                ));
            }
            let missing = match nonstrings[0] {
                Json::Null => "None".to_owned(),
                Json::Bool(true) => "True".to_owned(),
                Json::Bool(false) => "False".to_owned(),
                Json::Int(value) => value.to_string(),
                Json::Float(value) => json_text::python_float_repr(*value),
                _ => unreachable!("strings and unhashable objects were handled above"),
            };
            return Err(error(
                "ValueError",
                format!("Missing installer entrypoint {missing} for {name}"),
            ));
        }
        let mut ordered = StringSet::new(seed);
        for name in names {
            ordered.insert(name);
        }
        for function in ordered.into_order() {
            require(
                functions.contains(&function),
                format!("Missing installer entrypoint {function} for {name}"),
            )?;
        }
    }
    Ok(())
}

fn refresh_plan(resources: &Json, start: &str, seed: u32) -> Json {
    let mut affected = HashSet::from([start.to_owned()]);
    let mut ordered = vec![start.to_owned()];
    loop {
        let following: Vec<_> = resources
            .as_object()
            .unwrap()
            .iter()
            .filter(|(name, entry)| {
                !affected.contains(name)
                    && dependencies(entry, seed)
                        .iter()
                        .any(|name| affected.contains(name))
            })
            .map(|(name, _)| name.clone())
            .collect();
        if following.is_empty() {
            break;
        }
        affected.extend(following.iter().cloned());
        ordered.extend(following);
    }
    object([
        ("mode", Json::Str("read_only_inventory".into())),
        ("requested_resource", Json::Str(start.into())),
        ("note", Json::Str(PLAN_NOTE.into())),
        (
            "affected",
            Json::Array(
                ordered
                    .into_iter()
                    .map(|name| {
                        let entry = resources.get(&name).unwrap();
                        object([
                            ("resource", Json::Str(name)),
                            ("lifecycle", entry.get("lifecycle").unwrap().clone()),
                            ("refresh", entry.get("refresh").unwrap().clone()),
                        ])
                    })
                    .collect(),
            ),
        ),
    ])
}

fn path_value(text: &str) -> PathBuf {
    let prefix = if text.starts_with("//") && !text.starts_with("///") {
        "//"
    } else if text.starts_with('/') {
        "/"
    } else {
        ""
    };
    let parts: Vec<_> = text
        .split('/')
        .filter(|s| !s.is_empty() && *s != ".")
        .collect();
    let path = format!("{prefix}{}", parts.join("/"));
    PathBuf::from(if path.is_empty() { "." } else { &path })
}
fn path_repr(path: &Path) -> String {
    let text = path.to_string_lossy();
    let quote = if text.contains('\'') && !text.contains('"') {
        '"'
    } else {
        '\''
    };
    let mut out = String::from(quote);
    for c in text.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c == quote => {
                out.push('\\');
                out.push(c);
            }
            c if c.is_control() && (c as u32) <= 255 => {
                out.push_str(&format!("\\x{:02x}", c as u32))
            }
            c if c.is_control() || matches!(c, '\u{2028}' | '\u{2029}') => {
                out.push_str(&format!("\\u{:04x}", c as u32))
            }
            c => out.push(c),
        }
    }
    out.push(quote);
    out
}
fn file_error(cause: io::Error, path: &Path) -> Box<dyn std::error::Error> {
    let category = match cause.kind() {
        io::ErrorKind::NotFound => "FileNotFoundError",
        io::ErrorKind::PermissionDenied => "PermissionError",
        io::ErrorKind::IsADirectory => "IsADirectoryError",
        io::ErrorKind::NotADirectory => "NotADirectoryError",
        _ => "OSError",
    };
    let message = if let Some(errno) = cause.raw_os_error() {
        let message = cause.to_string();
        let plain = message
            .strip_suffix(&format!(" (os error {errno})"))
            .unwrap_or(&message);
        format!("[Errno {errno}] {plain}: {}", path_repr(path))
    } else {
        cause.to_string()
    };
    error(category, message)
}
fn read_text(path: &Path) -> Result<String> {
    let bytes = fs::read(path).map_err(|e| file_error(e, path))?;
    crate::python_io::decode_utf8(&bytes)
        .map(|s| s.replace("\r\n", "\n").replace('\r', "\n"))
        .map_err(|message| error("UnicodeDecodeError", message))
}

/// Python's object_pairs_hook runs when each object closes, rather than at
/// the second occurrence of a key. Keep the original successful parse path
/// unchanged. Only a duplicate error needs this local precedence adapter.
///
/// The shared parser establishes the first syntax-error offset. Before that
/// offset the token stream is valid; a small quote/brace walk merely locates
/// completed object spans. Parsing those spans in completion order preserves
/// child-hook errors before parent hooks or later syntax errors. Each byte is
/// reparsed at most once per enclosing object (the shared depth bound applies),
/// rather than repeatedly parsing a growing document prefix.
fn parse_registry(text: &str) -> std::result::Result<Json, String> {
    let first = match json_text::parse(text, true) {
        Ok(value) => return Ok(value),
        Err(message) if message.starts_with("Duplicate manifest key: ") => message,
        Err(message) => return Err(message),
    };
    let syntax = json_text::parse(text, false).err();
    let limit = if let Some(message) = &syntax {
        if json_text::is_refusal(message) {
            return Err(message.clone());
        }
        let character = message
            .rsplit_once("(char ")
            .and_then(|(_, s)| s.strip_suffix(')'))
            .and_then(|s| s.parse::<usize>().ok());
        let Some(character) = character else {
            return Err(message.clone());
        };
        text.char_indices()
            .nth(character)
            .map_or(text.len(), |(i, _)| i)
    } else {
        text.len()
    };
    let mut starts = Vec::new();
    let mut quoted = false;
    let mut escaped = false;
    for (at, byte) in text.as_bytes()[..limit].iter().copied().enumerate() {
        if quoted {
            if escaped {
                escaped = false;
            } else if byte == b'\\' {
                escaped = true;
            } else if byte == b'"' {
                quoted = false;
            }
        } else {
            match byte {
                b'"' => quoted = true,
                b'{' => starts.push(at),
                b'}' => {
                    if let Some(start) = starts.pop() {
                        json_text::parse(&text[start..=at], true)?;
                    }
                }
                _ => {}
            }
        }
    }
    Err(syntax.unwrap_or(first))
}
fn load_registry(path: &Path, seed: u32) -> Result<Json> {
    let text = read_text(path)?;
    let registry = parse_registry(&text).map_err(|message| {
        if let Some(key) = message.strip_prefix("Duplicate manifest key: ") {
            error("ValueError", format!("Duplicate registry key: {key}"))
        } else if json_text::is_refusal(&message) {
            error(
                "RuntimeError",
                format!("Unsupported JSON profile: {message}"),
            )
        } else {
            error("JSONDecodeError", message)
        }
    })?;
    validate_registry(&registry, seed)?;
    Ok(registry)
}
fn seed() -> Result<u32> {
    std::env::var("PYTHONHASHSEED").ok().filter(|s| !s.is_empty() && s.bytes().all(|c| c.is_ascii_digit()))
        .and_then(|s| s.parse().ok()).ok_or_else(|| error("RuntimeError", "Unsupported profile: set an explicit PYTHONHASHSEED in 0..4294967295; another Python process's random secret is unknown"))
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(command) = args.first().map(String::as_str) else {
        return Err(invalid(HELP));
    };
    if matches!(command, "--help" | "-h") {
        println!("{HELP}");
        return Ok(());
    }
    let (names, flags): (&[&str], &[&str]) = match command {
        "list" => (&["registry"], &["refreshable"]),
        "show" | "refresh-plan" => (&["registry", "resource"], &[]),
        "validate" => (&["registry"], &["check-routes"]),
        _ => return Err(invalid(HELP)),
    };
    let Parsed::Values(options) =
        parse_long_options_with_flags(&args[1..], names, flags).map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    require_options(&options, names).map_err(invalid)?;
    let path = path_value(&options["registry"]);
    crate::builder_error::finish((|| {
        let seed = seed()?;
        let registry = load_registry(&path, seed)?;
        let resources = registry.get("resources").unwrap();
        match command {
            "validate" => {
                if options.contains_key("check-routes") {
                    let actual = fs::canonicalize(&path).map_err(|e| file_error(e, &path))?;
                    check_routes(&registry, actual.parent().unwrap(), seed)?;
                }
                println!(
                    "{}",
                    json_text::dumps(
                        &object([
                            ("valid", Json::Bool(true)),
                            (
                                "resources",
                                Json::Int(resources.as_object().unwrap().len() as i128)
                            ),
                            ("cache_payloads_read", Json::Bool(false))
                        ]),
                        false
                    )
                    .map_err(invalid)?
                );
            }
            "list" => {
                println!("resource\tlifecycle\tautomatic\trefresh_entrypoints");
                for (name, entry) in resources.as_object().unwrap() {
                    let refresh = entry.get("refresh").unwrap();
                    if options.contains_key("refreshable")
                        && matches!(refresh.get("keep_function"), Some(Json::Bool(false)))
                    {
                        continue;
                    }
                    println!(
                        "{name}\t{}\t{}\t{}",
                        entry.get("lifecycle").unwrap().as_str().unwrap(),
                        if matches!(refresh.get("automatic"), Some(Json::Bool(true))) {
                            "true"
                        } else {
                            "false"
                        },
                        strings(refresh.get("entrypoints").unwrap())
                            .unwrap()
                            .join(",")
                    );
                }
            }
            _ => {
                let name = &options["resource"];
                require(
                    resources.get(name).is_some(),
                    format!("Unknown resource: {name}"),
                )?;
                let result = if command == "show" {
                    resources.get(name).unwrap().clone()
                } else {
                    refresh_plan(resources, name, seed)
                };
                println!("{}", json_text::dumps_indented(&result).map_err(invalid)?);
            }
        }
        Ok(())
    })())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn installer_pattern_keeps_python_whitespace_and_multiline_anchor() {
        let shell = "function one() {}\n function hidden() {}\nfunction\n_two\u{2003}(\nfunction a-b() {}\nfunction\u{1c}four(\nfunction 9bad(\nfunctionlast()";
        assert_eq!(
            shell_functions(shell),
            HashSet::from(["one".into(), "_two".into(), "four".into()])
        );
    }
    #[test]
    fn identifiers_are_ascii_and_whole_string_matches() {
        assert!(id_valid("a_9"));
        for value in ["", "A", "a\n", "á", "9a", "a-"] {
            assert!(!id_valid(value));
        }
    }
    #[test]
    fn duplicate_hooks_wait_for_object_completion() {
        assert_eq!(
            parse_registry(r#"{"a":1,"a":2,"b":}"#).unwrap_err(),
            "Expecting value: line 1 column 18 (char 17)"
        );
        assert_eq!(
            parse_registry(r#"{"a":1,"a":2,"b":{"x":1,"x":2}}"#).unwrap_err(),
            "Duplicate manifest key: x"
        );
        assert_eq!(
            parse_registry(r#"{"x":{"a":1,"a":2},"b":}"#).unwrap_err(),
            "Duplicate manifest key: a"
        );
        assert_eq!(
            parse_registry(r#"{"x":{"a":1,"a":2},"b":0}"#).unwrap_err(),
            "Duplicate manifest key: a"
        );
    }
}
