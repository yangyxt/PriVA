//! `jsonschema.Draft202012Validator(schema).iter_errors(result)`, sorted by
//! error path, with jsonschema 4.23's messages for the keywords it checks.
use super::*;
use crate::python_text::is_python_space;
use std::cmp::Ordering;

#[derive(Clone, Debug)]
enum Step {
    Key(String),
    Index(usize),
}

struct Failure {
    path: Vec<Step>,
    message: String,
}

const SUPPORTED: [&str; 23] = [
    "$ref",
    "propertyNames",
    "type",
    "required",
    "properties",
    "patternProperties",
    "additionalProperties",
    "items",
    "const",
    "enum",
    "minLength",
    "maxLength",
    "minItems",
    "maxItems",
    "minProperties",
    "maxProperties",
    "pattern",
    "minimum",
    "maximum",
    "uniqueItems",
    "allOf",
    "anyOf",
    "if",
];
const REFUSED: [&str; 15] = [
    "$dynamicRef",
    "contains",
    "dependentRequired",
    "dependentSchemas",
    "exclusiveMaximum",
    "exclusiveMinimum",
    "format",
    "multipleOf",
    "not",
    "oneOf",
    "prefixItems",
    "unevaluatedItems",
    "unevaluatedProperties",
    "maxContains",
    "minContains",
];

/// jsonschema's `_utils.equal`: booleans never equal numbers.
fn same(left: &Json, right: &Json) -> bool {
    match (left, right) {
        (Json::Str(a), Json::Str(b)) => a == b,
        (Json::Str(_), _) | (_, Json::Str(_)) => false,
        (Json::Array(a), Json::Array(b)) => a.len() == b.len() && a.iter().zip(b).all(|(x, y)| same(x, y)),
        (Json::Object(a), Json::Object(b)) => {
            a.len() == b.len()
                && a.iter().all(|(key, value)| b.iter().any(|(other, item)| other == key && same(value, item)))
        }
        (Json::Bool(a), Json::Bool(b)) => a == b,
        (Json::Bool(_), _) | (_, Json::Bool(_)) => false,
        _ => equal(left, right),
    }
}

fn is_type(value: &Json, name: &str) -> bool {
    match name {
        "object" => matches!(value, Json::Object(_)),
        "array" => matches!(value, Json::Array(_)),
        "string" => matches!(value, Json::Str(_)),
        "integer" => match value {
            Json::Int(_) => true,
            Json::Float(number) => number.is_finite() && number.fract() == 0.0,
            _ => false,
        },
        "number" => matches!(value, Json::Int(_) | Json::Float(_)),
        "boolean" => matches!(value, Json::Bool(_)),
        "null" => matches!(value, Json::Null),
        _ => false,
    }
}

fn number(value: &Json) -> Option<f64> {
    match value {
        Json::Int(n) => Some(*n as f64),
        Json::Float(n) => Some(*n),
        _ => None,
    }
}

/// `a < b` between two JSON numbers.
fn less(a: &Json, b: &Json) -> bool {
    match (a, b) {
        (Json::Int(x), Json::Int(y)) => x < y,
        _ => match (number(a), number(b)) {
            (Some(x), Some(y)) => x < y,
            _ => false,
        },
    }
}

fn length(value: &Json) -> Option<usize> {
    match value {
        Json::Str(text) => Some(text.chars().count()),
        Json::Array(items) => Some(items.len()),
        Json::Object(items) => Some(items.len()),
        _ => None,
    }
}

struct Validator<'a> {
    root: &'a Json,
}

impl Validator<'_> {
    fn resolve(&self, reference: &str) -> Result<&Json> {
        let pointer = reference
            .strip_prefix('#')
            .ok_or_else(|| error("Unresolvable", format!("Unresolvable: {reference}")))?;
        let mut node = self.root;
        for part in pointer.split('/').skip(1) {
            let part = part.replace("~1", "/").replace("~0", "~");
            node = match node {
                Json::Object(_) => node.get(&part),
                Json::Array(items) => part.parse::<usize>().ok().and_then(|i| items.get(i)),
                _ => None,
            }
            .ok_or_else(|| error("Unresolvable", format!("Unresolvable: {reference}")))?;
        }
        Ok(node)
    }

    fn valid(&self, instance: &Json, schema: &Json) -> Result<bool> {
        let mut found = Vec::new();
        self.check(instance, schema, &[], &mut found)?;
        Ok(found.is_empty())
    }

    fn check(&self, instance: &Json, schema: &Json, path: &[Step], out: &mut Vec<Failure>) -> Result<()> {
        let fail = |out: &mut Vec<Failure>, message: String| {
            out.push(Failure {
                path: path.to_vec(),
                message,
            })
        };
        let entries = match schema {
            Json::Bool(true) => return Ok(()),
            Json::Bool(false) => {
                fail(out, format!("False schema does not allow {}", repr(instance)));
                return Ok(());
            }
            Json::Object(entries) => entries,
            _ => return Ok(()),
        };
        let child = |step: Step| {
            let mut next = path.to_vec();
            next.push(step);
            next
        };
        for (keyword, value) in entries {
            let keyword = keyword.as_str();
            if REFUSED.contains(&keyword) {
                return Err(error(
                    "ValueError",
                    format!("Unsupported JSON Schema keyword in the native validator: {keyword}"),
                ));
            }
            if !SUPPORTED.contains(&keyword) {
                continue;
            }
            match keyword {
                "$ref" => {
                    let target = self.resolve(value.as_str().unwrap_or(""))?;
                    self.check(instance, target, path, out)?;
                }
                "type" => {
                    let types: Vec<String> = match value {
                        Json::Array(items) => items.iter().map(py_str).collect(),
                        other => vec![py_str(other)],
                    };
                    if !types.iter().any(|t| is_type(instance, t)) {
                        let names: Vec<String> = types.iter().map(|t| crate::python_path::repr(t)).collect();
                        fail(out, format!("{} is not of type {}", repr(instance), names.join(", ")));
                    }
                }
                "required" => {
                    if let (Json::Object(_), Json::Array(names)) = (instance, value) {
                        for name in names {
                            if !instance.get(&py_str(name)).is_some() {
                                fail(out, format!("{} is a required property", repr(name)));
                            }
                        }
                    }
                }
                "properties" => {
                    if let (Json::Object(_), Json::Object(properties)) = (instance, value) {
                        for (name, subschema) in properties {
                            if let Some(item) = instance.get(name) {
                                self.check(item, subschema, &child(Step::Key(name.clone())), out)?;
                            }
                        }
                    }
                }
                "patternProperties" => {
                    if let (Json::Object(fields), Json::Object(patterns)) = (instance, value) {
                        for (pattern, subschema) in patterns {
                            for (name, item) in fields {
                                if search(pattern, name)? {
                                    self.check(item, subschema, &child(Step::Key(name.clone())), out)?;
                                }
                            }
                        }
                    }
                }
                "additionalProperties" => {
                    let Json::Object(fields) = instance else { continue };
                    let properties = schema.get("properties");
                    let patterns: Vec<String> = match schema.get("patternProperties") {
                        Some(Json::Object(items)) => items.iter().map(|(k, _)| k.clone()).collect(),
                        _ => vec![],
                    };
                    let mut extras = Vec::new();
                    for (name, _) in fields {
                        if properties.is_some_and(|p| p.get(name).is_some()) {
                            continue;
                        }
                        let mut matched = false;
                        for pattern in &patterns {
                            if search(pattern, name)? {
                                matched = true;
                                break;
                            }
                        }
                        if !matched {
                            extras.push(name.clone());
                        }
                    }
                    if let Json::Object(_) = value {
                        for name in &extras {
                            self.check(instance.get(name).unwrap(), value, &child(Step::Key(name.clone())), out)?;
                        }
                    } else if !truthy(value) && !extras.is_empty() {
                        extras.sort();
                        let joined = extras.iter().map(|e| crate::python_path::repr(e)).collect::<Vec<_>>().join(", ");
                        if schema.get("patternProperties").is_some() {
                            let mut sorted = patterns.clone();
                            sorted.sort();
                            let shown = sorted.iter().map(|p| crate::python_path::repr(p)).collect::<Vec<_>>().join(", ");
                            let verb = if extras.len() == 1 { "does" } else { "do" };
                            fail(out, format!("{joined} {verb} not match any of the regexes: {shown}"));
                        } else {
                            let verb = if extras.len() == 1 { "was" } else { "were" };
                            fail(out, format!("Additional properties are not allowed ({joined} {verb} unexpected)"));
                        }
                    }
                }
                "items" => {
                    if let Json::Array(items) = instance {
                        for (index, item) in items.iter().enumerate() {
                            self.check(item, value, &child(Step::Index(index)), out)?;
                        }
                    }
                }
                "const" => {
                    if !same(instance, value) {
                        fail(out, format!("{} was expected", repr(value)));
                    }
                }
                "enum" => {
                    let options: Vec<Json> = match value {
                        Json::Array(items) => items.clone(),
                        other => vec![other.clone()],
                    };
                    if options.iter().all(|option| !same(option, instance)) {
                        fail(out, format!("{} is not one of {}", repr(instance), repr(value)));
                    }
                }
                "minLength" | "maxLength" => {
                    if let (Json::Str(_), Some(limit)) = (instance, number(value)) {
                        let size = length(instance).unwrap() as f64;
                        if keyword == "minLength" && size < limit {
                            let words = if limit == 1.0 { "should be non-empty" } else { "is too short" };
                            fail(out, format!("{} {words}", repr(instance)));
                        } else if keyword == "maxLength" && size > limit {
                            let words = if limit == 0.0 { "is expected to be empty" } else { "is too long" };
                            fail(out, format!("{} {words}", repr(instance)));
                        }
                    }
                }
                "minItems" | "maxItems" => {
                    if let (Json::Array(_), Some(limit)) = (instance, number(value)) {
                        let size = length(instance).unwrap() as f64;
                        if keyword == "minItems" && size < limit {
                            let words = if limit == 1.0 { "should be non-empty" } else { "is too short" };
                            fail(out, format!("{} {words}", repr(instance)));
                        } else if keyword == "maxItems" && size > limit {
                            let words = if limit == 0.0 { "is expected to be empty" } else { "is too long" };
                            fail(out, format!("{} {words}", repr(instance)));
                        }
                    }
                }
                "minProperties" | "maxProperties" => {
                    if let (Json::Object(_), Some(limit)) = (instance, number(value)) {
                        let size = length(instance).unwrap() as f64;
                        if keyword == "minProperties" && size < limit {
                            let words = if limit == 1.0 { "should be non-empty" } else { "does not have enough properties" };
                            fail(out, format!("{} {words}", repr(instance)));
                        } else if keyword == "maxProperties" && size > limit {
                            let words = if limit == 0.0 { "is expected to be empty" } else { "has too many properties" };
                            fail(out, format!("{} {words}", repr(instance)));
                        }
                    }
                }
                "pattern" => {
                    if let Json::Str(body) = instance {
                        let pattern = py_str(value);
                        if !search(&pattern, body)? {
                            fail(out, format!("{} does not match {}", repr(instance), crate::python_path::repr(&pattern)));
                        }
                    }
                }
                "minimum" | "maximum" => {
                    if is_type(instance, "number") {
                        if keyword == "minimum" && less(instance, value) {
                            fail(out, format!("{} is less than the minimum of {}", repr(instance), repr(value)));
                        } else if keyword == "maximum" && less(value, instance) {
                            fail(out, format!("{} is greater than the maximum of {}", repr(instance), repr(value)));
                        }
                    }
                }
                "uniqueItems" => {
                    if let Json::Array(items) = instance {
                        if truthy(value) {
                            let repeated = items
                                .iter()
                                .enumerate()
                                .any(|(i, a)| items[i + 1..].iter().any(|b| same(a, b)));
                            if repeated {
                                fail(out, format!("{} has non-unique elements", repr(instance)));
                            }
                        }
                    }
                }
                "allOf" => {
                    if let Json::Array(subschemas) = value {
                        for subschema in subschemas {
                            self.check(instance, subschema, path, out)?;
                        }
                    }
                }
                "anyOf" => {
                    if let Json::Array(subschemas) = value {
                        let mut passed = false;
                        for subschema in subschemas {
                            if self.valid(instance, subschema)? {
                                passed = true;
                                break;
                            }
                        }
                        if !passed {
                            fail(out, format!("{} is not valid under any of the given schemas", repr(instance)));
                        }
                    }
                }
                "if" => {
                    if self.valid(instance, value)? {
                        if let Some(then) = schema.get("then") {
                            self.check(instance, then, path, out)?;
                        }
                    } else if let Some(otherwise) = schema.get("else") {
                        self.check(instance, otherwise, path, out)?;
                    }
                }
                "propertyNames" => {
                    if let Json::Object(fields) = instance {
                        for (name, _) in fields {
                            self.check(&text(name.clone()), value, path, out)?;
                        }
                    }
                }
                _ => {}
            }
        }
        Ok(())
    }
}

fn compare(a: &[Step], b: &[Step]) -> Ordering {
    for (x, y) in a.iter().zip(b) {
        let order = match (x, y) {
            (Step::Key(x), Step::Key(y)) => x.cmp(y),
            (Step::Index(x), Step::Index(y)) => x.cmp(y),
            (Step::Key(_), Step::Index(_)) => Ordering::Greater,
            (Step::Index(_), Step::Key(_)) => Ordering::Less,
        };
        if order != Ordering::Equal {
            return order;
        }
    }
    a.len().cmp(&b.len())
}

/// The script's `validate_canonical_json`.
pub(super) fn validate_canonical(result: &Json, schema_path: &Path) -> Result<()> {
    let schema = read_json(schema_path)?;
    if !truthy(&schema) {
        return Err(error(
            "FileNotFoundError",
            format!("Canonical JSON Schema is missing: {}", schema_path.display()),
        ));
    }
    let validator = Validator { root: &schema };
    let mut failures = Vec::new();
    validator.check(result, &schema, &[], &mut failures)?;
    if failures.is_empty() {
        return Ok(());
    }
    failures.sort_by(|a, b| compare(&a.path, &b.path));
    let examples: Vec<String> = failures
        .iter()
        .take(10)
        .map(|failure| {
            let location = failure
                .path
                .iter()
                .map(|step| match step {
                    Step::Key(key) => key.clone(),
                    Step::Index(index) => index.to_string(),
                })
                .collect::<Vec<_>>()
                .join("/");
            format!(
                "{}: {}",
                if location.is_empty() { "<root>" } else { &location },
                failure.message
            )
        })
        .collect();
    Err(error(
        "ValueError",
        format!(
            "Canonical JSON failed schema validation ({} errors): {}",
            failures.len(),
            examples.join(" | ")
        ),
    ))
}

#[derive(Clone, Debug)]
enum Atom {
    Literal(char),
    Any,
    Class(Vec<(char, char)>, Vec<char>, bool),
    Start,
    End,
}

#[derive(Clone, Debug)]
struct Piece {
    atom: Atom,
    min: usize,
    max: usize,
}

fn unsupported(pattern: &str) -> Box<dyn std::error::Error> {
    error(
        "ValueError",
        format!("Unsupported schema pattern in the native validator: {pattern}"),
    )
}

/// Class shorthands as single-letter markers: d, D, w, W, s, S.
fn shorthand_matches(kind: char, c: char) -> bool {
    match kind {
        'd' => crate::python_text::decimal_digit_value(c).is_some(),
        'D' => crate::python_text::decimal_digit_value(c).is_none(),
        'w' => crate::gene_mechanism_builder::is_word(c),
        'W' => !crate::gene_mechanism_builder::is_word(c),
        's' => is_python_space(c),
        'S' => !is_python_space(c),
        _ => false,
    }
}

fn parse_pattern(pattern: &str) -> Result<Vec<Piece>> {
    let chars: Vec<char> = pattern.chars().collect();
    let mut pieces: Vec<Piece> = Vec::new();
    let mut at = 0;
    while at < chars.len() {
        let c = chars[at];
        let atom = match c {
            '^' => Atom::Start,
            '$' => Atom::End,
            '.' => Atom::Any,
            '\\' => {
                at += 1;
                let next = *chars.get(at).ok_or_else(|| unsupported(pattern))?;
                if "dDwWsS".contains(next) {
                    Atom::Class(vec![], vec![next], false)
                } else if next.is_ascii_alphanumeric() {
                    return Err(unsupported(pattern));
                } else {
                    Atom::Literal(next)
                }
            }
            '[' => {
                at += 1;
                let negated = chars.get(at) == Some(&'^');
                if negated {
                    at += 1;
                }
                let mut ranges = Vec::new();
                let mut shorthands = Vec::new();
                let mut first = true;
                loop {
                    let c = *chars.get(at).ok_or_else(|| unsupported(pattern))?;
                    if c == ']' && !first {
                        break;
                    }
                    first = false;
                    let low = if c == '\\' {
                        at += 1;
                        let next = *chars.get(at).ok_or_else(|| unsupported(pattern))?;
                        if "dDwWsS".contains(next) {
                            shorthands.push(next);
                            at += 1;
                            continue;
                        }
                        if next.is_ascii_alphanumeric() {
                            return Err(unsupported(pattern));
                        }
                        next
                    } else {
                        c
                    };
                    if chars.get(at + 1) == Some(&'-') && chars.get(at + 2).is_some_and(|c| *c != ']') {
                        let high = chars[at + 2];
                        if high == '\\' {
                            return Err(unsupported(pattern));
                        }
                        ranges.push((low, high));
                        at += 3;
                    } else {
                        ranges.push((low, low));
                        at += 1;
                    }
                }
                Atom::Class(ranges, shorthands, negated)
            }
            '(' | ')' | '|' | '{' | '}' | '*' | '+' | '?' => return Err(unsupported(pattern)),
            other => Atom::Literal(other),
        };
        at += 1;
        let (min, max) = match chars.get(at) {
            Some('*') => (0, usize::MAX),
            Some('+') => (1, usize::MAX),
            Some('?') => (0, 1),
            Some('{') => return Err(unsupported(pattern)),
            _ => (1, 1),
        };
        if (min, max) != (1, 1) {
            if matches!(atom, Atom::Start | Atom::End) {
                return Err(unsupported(pattern));
            }
            at += 1;
            if matches!(chars.get(at), Some('?' | '+')) {
                return Err(unsupported(pattern));
            }
        }
        pieces.push(Piece { atom, min, max });
    }
    Ok(pieces)
}

fn atom_matches(atom: &Atom, c: char) -> bool {
    match atom {
        Atom::Literal(want) => *want == c,
        Atom::Any => c != '\n',
        Atom::Class(ranges, shorthands, negated) => {
            let hit = ranges.iter().any(|(low, high)| *low <= c && c <= *high)
                || shorthands.iter().any(|kind| shorthand_matches(*kind, c));
            hit != *negated
        }
        Atom::Start | Atom::End => false,
    }
}

fn match_here(pieces: &[Piece], text: &[char], at: usize) -> bool {
    let Some((piece, rest)) = pieces.split_first() else {
        return true;
    };
    match piece.atom {
        Atom::Start => at == 0 && match_here(rest, text, at),
        Atom::End => {
            (at == text.len() || (at + 1 == text.len() && text[at] == '\n')) && match_here(rest, text, at)
        }
        _ => {
            let mut count = 0;
            while count < piece.max && at + count < text.len() && atom_matches(&piece.atom, text[at + count]) {
                count += 1;
            }
            loop {
                if count >= piece.min && match_here(rest, text, at + count) {
                    return true;
                }
                if count == 0 || count <= piece.min {
                    return false;
                }
                count -= 1;
            }
        }
    }
}

/// `re.search(pattern, text) is not None` for the supported syntax.
fn search(pattern: &str, text: &str) -> Result<bool> {
    let pieces = parse_pattern(pattern)?;
    let chars: Vec<char> = text.chars().collect();
    Ok((0..=chars.len()).any(|start| match_here(&pieces, &chars, start)))
}
