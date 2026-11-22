//! Rust version of the lookup half of `scripts/protein_domain_mapping.py`
//! (class `DomainNormalizer`): turn a VEP domain annotation such as
//! `Pfam:PF00533` into InterPro entries and decide whether a domain is
//! functional. The splicing step (and later the PVS1 and domain steps) call it.
//!
//! The script's own command line (`--vcf --xml --report --mapping_output`,
//! which builds the InterPro pickle at installation from VEP's CSQ header and
//! InterPro's XML) is not converted here; it stays Python.
//!
//! Two data files are read when a `DomainNormalizer` is made, from the
//! repository's `data/InterPro` folder (the folder given by the environment
//! variable `PRIVA_DATA_DIR` if it is set, else `data` next to the `rust`
//! folder this binary was built from):
//! - `InterPro_domain_norm_map.json`: VEP database names and accession
//!   formatting rules. Rules of type `regex_capture` are matched with a small
//!   regular-expression subset (see `Pattern`); anything outside it fails.
//! - `curated_InterPro_func_domains.tsv.gz`: column 1 is the InterPro id and
//!   column 9 its curated functional status (`yes` or `no`), read the way
//!   Python's `csv.reader` reads a tab-separated file.
//!
//! `normalize_interpro_to_vep` is not ported: nothing in the steps converted
//! so far calls it.
//!
//! Where Python walks a set of normalised domains (whose order changes from
//! run to run), the entries are visited in sorted order here; the answer only
//! depends on the order when an unusable entry sits before a functional one.
use crate::cli::{Result, invalid};
use crate::interpro_mapping_io::{GoTerm, GoTerms, InterproEntry, InterproMapping};
use flate2::read::MultiGzDecoder;
use std::collections::{BTreeSet, HashMap};
use std::fs::{self, File};
use std::io::Read;
use std::path::{Path, PathBuf};

/// A failure at lookup time, which stops the run as Python's exception does.
pub(crate) type Lookup<T> = std::result::Result<T, String>;

/// Domains VEP predicts itself rather than takes from a member database.
const PREDICTED_DOMAINS: &[&str] = &[
    "Cleavage_site_(Signalp)",
    "Coiled-coils_(Ncoils)",
    "Low_complexity_(Seg)",
    "Transmembrane_helices",
];
/// The predicted domains minus low complexity, which count as functional.
const FUNCTIONAL_PREDICTED_DOMAINS: &[&str] = &[
    "Cleavage_site_(Signalp)",
    "Coiled-coils_(Ncoils)",
    "Transmembrane_helices",
];
const FUNCTIONAL_PROCESSES: &[&str] = &[
    "catalysis",
    "transport",
    "signaling",
    "regulation",
    "biosynthetic",
    "metabolic",
    "synthesis",
    "activity",
    "binding",
    "interaction",
    "transcription",
    "replication",
    "repair",
    "modification",
];
const FUNCTIONAL_KEYWORDS: &[&str] = &[
    "enzyme",
    "activity",
    "binding",
    "receptor",
    "kinase",
    "transferase",
    "hydrolase",
    "isomerase",
    "ligase",
    "lyase",
    "reductase",
    "oxidase",
    "synthase",
    "synthetase",
    "phosphatase",
    "channel",
    "transporter",
    "carrier",
    "permease",
    "catalytic",
    "regulatory",
    "effector",
    "interaction",
    "recognition",
    "transcription",
    "replication",
    "repair",
    "modification",
    "inhibitor",
    "activator",
    "homeobox",
    "helix-loop-helix",
    "zinc finger",
    "zn finger",
    "sh2",
    "sh3",
    "wd40",
    "ankyrin",
    "bzip",
    "bromodomain",
    "chromodomain",
    "active site",
    "binding site",
    "catalytic site",
    "substrate",
    "ligand",
];
/// Patterns of the form `<text>$`: the name ends with the text (or with the
/// text and one final line break, as a regular expression's `$` allows).
const FUNCTIONAL_SUFFIXES: &[&str] = &["ase", "_dom", "domain", "motif", "finger", "box"];
/// Patterns that are plain text searched anywhere in the name.
const FUNCTIONAL_PATTERNS: &[&str] = &[
    "active_site",
    "binding_site",
    "catalytic_site",
    "dna-binding",
    "rna-binding",
    "protein-binding",
];

/// The result of `interpret_functionality`; Python returns the strings
/// `"Functional"` and `"Unknown"`, or `nan` when nothing is known.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Functionality {
    Functional,
    Unknown,
    Missing,
}

/// What `query_interpro_entry_vep_anno` finds for one annotation.
pub(crate) struct DomainInfo<'a> {
    pub(crate) vep_domains: Vec<String>,
    /// `None` when no domain of the annotation maps to an InterPro database.
    pub(crate) interpro_entries: Option<Vec<&'a InterproEntry>>,
}

pub(crate) fn data_directory() -> PathBuf {
    match std::env::var_os("PRIVA_DATA_DIR") {
        Some(folder) if !folder.is_empty() => PathBuf::from(folder),
        _ => Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("..")
            .join("data"),
    }
}

// ---------------------------------------------------------------- JSON

#[derive(Clone, Debug, PartialEq)]
enum Json {
    Null,
    Bool(bool),
    Number(f64),
    Text(String),
    Array(Vec<Json>),
    Object(Vec<(String, Json)>),
}

impl Json {
    /// The last value written for a key, as Python's dict keeps it.
    fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Object(pairs) => pairs
                .iter()
                .rev()
                .find(|(name, _)| name.as_str() == key)
                .map(|(_, value)| value),
            _ => None,
        }
    }

    fn as_text(&self) -> Option<&str> {
        match self {
            Json::Text(text) => Some(text.as_str()),
            _ => None,
        }
    }

    /// Python's truth value.
    fn truthy(&self) -> bool {
        match self {
            Json::Null => false,
            Json::Bool(flag) => *flag,
            Json::Number(number) => *number != 0.0,
            Json::Text(text) => !text.is_empty(),
            Json::Array(items) => !items.is_empty(),
            Json::Object(pairs) => !pairs.is_empty(),
        }
    }

    fn into_map(self) -> Option<HashMap<String, Json>> {
        match self {
            Json::Object(pairs) => Some(pairs.into_iter().collect()),
            _ => None,
        }
    }
}

struct JsonParser<'a> {
    text: &'a [u8],
    position: usize,
}

impl JsonParser<'_> {
    fn peek(&self) -> Option<u8> {
        self.text.get(self.position).copied()
    }

    fn skip_whitespace(&mut self) {
        while matches!(self.peek(), Some(b' ' | b'\t' | b'\n' | b'\r')) {
            self.position += 1;
        }
    }

    fn expect(&mut self, byte: u8) -> Result<()> {
        if self.peek() == Some(byte) {
            self.position += 1;
            Ok(())
        } else {
            Err(invalid(format!(
                "JSON: expected {:?} at byte {}",
                byte as char, self.position
            )))
        }
    }

    fn keyword(&mut self, word: &str, value: Json) -> Result<Json> {
        if self.text[self.position..].starts_with(word.as_bytes()) {
            self.position += word.len();
            Ok(value)
        } else {
            Err(invalid(format!(
                "JSON: unexpected text at byte {}",
                self.position
            )))
        }
    }

    fn value(&mut self) -> Result<Json> {
        match self.peek() {
            Some(b'{') => self.object(),
            Some(b'[') => self.array(),
            Some(b'"') => Ok(Json::Text(self.string()?)),
            Some(b't') => self.keyword("true", Json::Bool(true)),
            Some(b'f') => self.keyword("false", Json::Bool(false)),
            Some(b'n') => self.keyword("null", Json::Null),
            Some(b'-' | b'0'..=b'9') => self.number(),
            _ => Err(invalid(format!(
                "JSON: unexpected text at byte {}",
                self.position
            ))),
        }
    }

    fn object(&mut self) -> Result<Json> {
        self.expect(b'{')?;
        let mut pairs = Vec::new();
        self.skip_whitespace();
        if self.peek() == Some(b'}') {
            self.position += 1;
            return Ok(Json::Object(pairs));
        }
        loop {
            self.skip_whitespace();
            let key = self.string()?;
            self.skip_whitespace();
            self.expect(b':')?;
            self.skip_whitespace();
            let value = self.value()?;
            pairs.push((key, value));
            self.skip_whitespace();
            match self.peek() {
                Some(b',') => self.position += 1,
                Some(b'}') => {
                    self.position += 1;
                    return Ok(Json::Object(pairs));
                }
                _ => {
                    return Err(invalid(format!(
                        "JSON: expected , or }} at byte {}",
                        self.position
                    )));
                }
            }
        }
    }

    fn array(&mut self) -> Result<Json> {
        self.expect(b'[')?;
        let mut items = Vec::new();
        self.skip_whitespace();
        if self.peek() == Some(b']') {
            self.position += 1;
            return Ok(Json::Array(items));
        }
        loop {
            self.skip_whitespace();
            items.push(self.value()?);
            self.skip_whitespace();
            match self.peek() {
                Some(b',') => self.position += 1,
                Some(b']') => {
                    self.position += 1;
                    return Ok(Json::Array(items));
                }
                _ => {
                    return Err(invalid(format!(
                        "JSON: expected , or ] at byte {}",
                        self.position
                    )));
                }
            }
        }
    }

    fn hex4(&mut self) -> Result<u32> {
        let digits = self
            .text
            .get(self.position..self.position + 4)
            .ok_or_else(|| invalid("JSON: a \\u escape is cut short"))?;
        let text = std::str::from_utf8(digits).map_err(|_| invalid("JSON: bad \\u escape"))?;
        let value = u32::from_str_radix(text, 16).map_err(|_| invalid("JSON: bad \\u escape"))?;
        self.position += 4;
        Ok(value)
    }

    fn string(&mut self) -> Result<String> {
        self.expect(b'"')?;
        let mut bytes: Vec<u8> = Vec::new();
        loop {
            let byte = self
                .peek()
                .ok_or_else(|| invalid("JSON: a string is not closed"))?;
            self.position += 1;
            match byte {
                b'"' => break,
                b'\\' => {
                    let escape = self
                        .peek()
                        .ok_or_else(|| invalid("JSON: a string ends inside an escape"))?;
                    self.position += 1;
                    match escape {
                        b'"' => bytes.push(b'"'),
                        b'\\' => bytes.push(b'\\'),
                        b'/' => bytes.push(b'/'),
                        b'b' => bytes.push(8),
                        b'f' => bytes.push(12),
                        b'n' => bytes.push(b'\n'),
                        b'r' => bytes.push(b'\r'),
                        b't' => bytes.push(b'\t'),
                        b'u' => {
                            let mut code = self.hex4()?;
                            if (0xd800..0xdc00).contains(&code)
                                && self.text[self.position..].starts_with(b"\\u")
                            {
                                self.position += 2;
                                let low = self.hex4()?;
                                if (0xdc00..0xe000).contains(&low) {
                                    code = 0x10000 + ((code - 0xd800) << 10) + (low - 0xdc00);
                                } else {
                                    return Err(invalid("JSON: a lone surrogate escape"));
                                }
                            }
                            let character = char::from_u32(code)
                                .ok_or_else(|| invalid("JSON: a lone surrogate escape"))?;
                            let mut buffer = [0_u8; 4];
                            bytes.extend_from_slice(character.encode_utf8(&mut buffer).as_bytes());
                        }
                        _ => return Err(invalid("JSON: an unknown escape")),
                    }
                }
                0..=0x1f => return Err(invalid("JSON: a control character inside a string")),
                _ => bytes.push(byte),
            }
        }
        String::from_utf8(bytes).map_err(|_| invalid("JSON: a string is not UTF-8"))
    }

    fn number(&mut self) -> Result<Json> {
        let start = self.position;
        if self.peek() == Some(b'-') {
            self.position += 1;
        }
        let digits_from = self.position;
        while matches!(self.peek(), Some(b'0'..=b'9')) {
            self.position += 1;
        }
        let whole = &self.text[digits_from..self.position];
        if whole.is_empty() || (whole.len() > 1 && whole[0] == b'0') {
            return Err(invalid(format!("JSON: bad number at byte {start}")));
        }
        if self.peek() == Some(b'.') {
            self.position += 1;
            let fraction_from = self.position;
            while matches!(self.peek(), Some(b'0'..=b'9')) {
                self.position += 1;
            }
            if self.position == fraction_from {
                return Err(invalid(format!("JSON: bad number at byte {start}")));
            }
        }
        if matches!(self.peek(), Some(b'e' | b'E')) {
            self.position += 1;
            if matches!(self.peek(), Some(b'+' | b'-')) {
                self.position += 1;
            }
            let exponent_from = self.position;
            while matches!(self.peek(), Some(b'0'..=b'9')) {
                self.position += 1;
            }
            if self.position == exponent_from {
                return Err(invalid(format!("JSON: bad number at byte {start}")));
            }
        }
        let text = std::str::from_utf8(&self.text[start..self.position])
            .map_err(|_| invalid("JSON: bad number"))?;
        let number: f64 = text
            .parse()
            .map_err(|_| invalid(format!("JSON: bad number at byte {start}")))?;
        Ok(Json::Number(number))
    }
}

fn parse_json(text: &str) -> Result<Json> {
    if text.starts_with('\u{feff}') {
        return Err(invalid("JSON: the text starts with a byte order mark"));
    }
    let mut parser = JsonParser {
        text: text.as_bytes(),
        position: 0,
    };
    parser.skip_whitespace();
    let value = parser.value()?;
    parser.skip_whitespace();
    if parser.position != parser.text.len() {
        return Err(invalid(format!(
            "JSON: extra data at byte {}",
            parser.position
        )));
    }
    Ok(value)
}

// ------------------------------------------------- regular expression subset

enum Atom {
    Digit,
    Any,
    Literal(char),
}

enum Item {
    Open(usize),
    Close(usize),
    Atom(Atom, usize, usize),
    End,
}

/// The part of Python's regular expressions that the curated accession rules
/// use: literals, `\d`, `.`, the quantifiers `*`, `+` and `?` after one of
/// those, capturing groups and a final `$`. Matching is greedy with
/// backtracking, as `re.match` is. Anything else is refused when the pattern
/// is compiled.
struct Pattern {
    items: Vec<Item>,
    groups: usize,
}

impl Pattern {
    fn compile(pattern: &str) -> Lookup<Pattern> {
        let unsupported = |why: &str| {
            format!("the regular expression {pattern:?} is outside the supported subset ({why})")
        };
        let chars: Vec<char> = pattern.chars().collect();
        let mut items = Vec::new();
        let mut groups = 0;
        let mut open: Vec<usize> = Vec::new();
        let mut index = 0;
        if chars.first() == Some(&'^') {
            index = 1;
        }
        while index < chars.len() {
            let character = chars[index];
            index += 1;
            let atom = match character {
                '(' => {
                    if chars.get(index) == Some(&'?') {
                        return Err(unsupported("group extension"));
                    }
                    groups += 1;
                    open.push(groups);
                    items.push(Item::Open(groups));
                    continue;
                }
                ')' => {
                    let Some(number) = open.pop() else {
                        return Err(unsupported("unbalanced parenthesis"));
                    };
                    items.push(Item::Close(number));
                    continue;
                }
                '$' => {
                    if index != chars.len() {
                        return Err(unsupported("$ before the end"));
                    }
                    items.push(Item::End);
                    continue;
                }
                '\\' => {
                    let Some(&escaped) = chars.get(index) else {
                        return Err(unsupported("trailing backslash"));
                    };
                    index += 1;
                    if escaped == 'd' {
                        Atom::Digit
                    } else if escaped.is_ascii_alphanumeric() {
                        return Err(unsupported("escape"));
                    } else {
                        Atom::Literal(escaped)
                    }
                }
                '.' => Atom::Any,
                '[' | ']' | '{' | '}' | '|' | '^' | '*' | '+' | '?' => {
                    return Err(unsupported("operator"));
                }
                other => Atom::Literal(other),
            };
            let (minimum, maximum) = match chars.get(index) {
                Some('*') => (0, usize::MAX),
                Some('+') => (1, usize::MAX),
                Some('?') => (0, 1),
                _ => (1, 1),
            };
            if (minimum, maximum) != (1, 1) {
                index += 1;
                if matches!(chars.get(index), Some('?' | '+' | '*' | '{')) {
                    return Err(unsupported("lazy or stacked quantifier"));
                }
            } else if chars.get(index) == Some(&'{') {
                return Err(unsupported("counted repetition"));
            }
            items.push(Item::Atom(atom, minimum, maximum));
        }
        if !open.is_empty() {
            return Err(unsupported("unbalanced parenthesis"));
        }
        Ok(Pattern { items, groups })
    }

    /// `re.match(pattern, text)`: the whole match and each group (`None` for a
    /// group that took no part), or `None` when the pattern does not match at
    /// the start of the text.
    fn match_start(&self, text: &str) -> Option<Vec<Option<String>>> {
        let chars: Vec<char> = text.chars().collect();
        let mut spans: Vec<Option<(usize, usize)>> = vec![None; self.groups + 1];
        let mut starts = vec![0; self.groups + 1];
        let mut end = 0;
        if !run(&self.items, &chars, 0, &mut spans, &mut starts, &mut end) {
            return None;
        }
        let mut found = vec![Some(chars[..end].iter().collect::<String>())];
        for number in 1..=self.groups {
            found.push(spans[number].map(|(from, to)| chars[from..to].iter().collect::<String>()));
        }
        Some(found)
    }
}

fn atom_matches(atom: &Atom, character: char) -> bool {
    match atom {
        Atom::Digit => character.is_ascii_digit(),
        Atom::Any => character != '\n',
        Atom::Literal(expected) => *expected == character,
    }
}

fn run(
    items: &[Item],
    text: &[char],
    position: usize,
    spans: &mut Vec<Option<(usize, usize)>>,
    starts: &mut Vec<usize>,
    end: &mut usize,
) -> bool {
    let Some((first, rest)) = items.split_first() else {
        *end = position;
        return true;
    };
    match first {
        Item::Open(number) => {
            let saved = starts[*number];
            starts[*number] = position;
            if run(rest, text, position, spans, starts, end) {
                true
            } else {
                starts[*number] = saved;
                false
            }
        }
        Item::Close(number) => {
            let saved = spans[*number];
            spans[*number] = Some((starts[*number], position));
            if run(rest, text, position, spans, starts, end) {
                true
            } else {
                spans[*number] = saved;
                false
            }
        }
        Item::End => {
            // `$` matches at the end, or before one final line break.
            let at_end =
                position == text.len() || (position + 1 == text.len() && text[position] == '\n');
            at_end && run(rest, text, position, spans, starts, end)
        }
        Item::Atom(atom, minimum, maximum) => {
            let mut count = 0;
            while count < *maximum
                && position + count < text.len()
                && atom_matches(atom, text[position + count])
            {
                count += 1;
            }
            if count < *minimum {
                return false;
            }
            loop {
                if run(rest, text, position + count, spans, starts, end) {
                    return true;
                }
                if count == *minimum {
                    return false;
                }
                count -= 1;
            }
        }
    }
}

// ------------------------------------------------------------ csv.reader

/// Python's `csv.reader(f, delimiter="\t")` over text that was read in
/// universal-newline mode (every line break is `\n`): a field starting with a
/// quote runs to the next lone quote, `""` is one quote, a quote inside a
/// field is a plain character, and an empty line is an empty row.
fn csv_rows(text: &str) -> Vec<Vec<String>> {
    #[derive(Clone, Copy, PartialEq)]
    enum State {
        StartRecord,
        StartField,
        InField,
        InQuotedField,
        QuoteInQuotedField,
    }
    let mut rows: Vec<Vec<String>> = Vec::new();
    let mut fields: Vec<String> = Vec::new();
    let mut field = String::new();
    let mut state = State::StartRecord;
    for character in text.chars() {
        loop {
            match state {
                State::StartRecord => {
                    if character == '\n' {
                        rows.push(Vec::new());
                        break;
                    }
                    state = State::StartField;
                }
                State::StartField => {
                    match character {
                        '\n' => {
                            fields.push(std::mem::take(&mut field));
                            rows.push(std::mem::take(&mut fields));
                            state = State::StartRecord;
                        }
                        '"' => state = State::InQuotedField,
                        '\t' => fields.push(std::mem::take(&mut field)),
                        _ => {
                            field.push(character);
                            state = State::InField;
                        }
                    }
                    break;
                }
                State::InField => {
                    match character {
                        '\n' => {
                            fields.push(std::mem::take(&mut field));
                            rows.push(std::mem::take(&mut fields));
                            state = State::StartRecord;
                        }
                        '\t' => {
                            fields.push(std::mem::take(&mut field));
                            state = State::StartField;
                        }
                        _ => field.push(character),
                    }
                    break;
                }
                State::InQuotedField => {
                    if character == '"' {
                        state = State::QuoteInQuotedField;
                    } else {
                        field.push(character);
                    }
                    break;
                }
                State::QuoteInQuotedField => {
                    match character {
                        '"' => {
                            field.push('"');
                            state = State::InQuotedField;
                        }
                        '\t' => {
                            fields.push(std::mem::take(&mut field));
                            state = State::StartField;
                        }
                        '\n' => {
                            fields.push(std::mem::take(&mut field));
                            rows.push(std::mem::take(&mut fields));
                            state = State::StartRecord;
                        }
                        _ => {
                            field.push(character);
                            state = State::InField;
                        }
                    }
                    break;
                }
            }
        }
    }
    // The end of the text ends the last row, even inside a quoted field.
    if state != State::StartRecord {
        fields.push(field);
        rows.push(fields);
    }
    rows
}

// ------------------------------------------------------- DomainNormalizer

pub(crate) struct DomainNormalizer {
    vep_to_interpro: HashMap<String, Json>,
    vep_to_interpro_rules: HashMap<String, Json>,
    curated_functionality: HashMap<String, String>,
}

impl DomainNormalizer {
    /// `DomainNormalizer()`: reads the two data files from the data directory.
    pub(crate) fn new() -> Result<DomainNormalizer> {
        let folder = data_directory().join("InterPro");
        let norm_path = folder.join("InterPro_domain_norm_map.json");
        let norm_text = fs::read_to_string(&norm_path)
            .map_err(|error| invalid(format!("{}: {error}", norm_path.display())))?;
        let curated_path = folder.join("curated_InterPro_func_domains.tsv.gz");
        let curated = read_curated_functional_domains(&curated_path)
            .map_err(|error| invalid(format!("{}: {error}", curated_path.display())))?;
        DomainNormalizer::from_parts(&norm_text, curated)
    }

    pub(crate) fn from_parts(
        norm_json: &str,
        curated_functionality: HashMap<String, String>,
    ) -> Result<DomainNormalizer> {
        let document = parse_json(norm_json)?;
        let section = |outer: &str, inner: &str| -> Result<HashMap<String, Json>> {
            document
                .get(outer)
                .and_then(|value| value.get(inner))
                .cloned()
                .and_then(Json::into_map)
                .ok_or_else(|| invalid(format!("the domain mapping JSON lacks {outer}/{inner}")))
        };
        let vep_to_interpro = section("database_name_mappings", "vep_to_interpro")?;
        section("database_name_mappings", "interpro_to_vep")?;
        let vep_to_interpro_rules = section("accession_id_normalization", "vep_to_interpro_key")?;
        section("accession_id_normalization", "interpro_key_to_vep_id")?;
        Ok(DomainNormalizer {
            vep_to_interpro,
            vep_to_interpro_rules,
            curated_functionality,
        })
    }

    /// `normalize_vep_to_interpro`: the InterPro database name and accession
    /// for a VEP database name and accession, or `None`.
    pub(crate) fn normalize_vep_to_interpro(
        &self,
        database: &str,
        accession: &str,
    ) -> Lookup<Option<(String, String)>> {
        if self.vep_to_interpro.is_empty() || self.vep_to_interpro_rules.is_empty() {
            return Ok(None);
        }
        let Some(mapped) = self.vep_to_interpro.get(database) else {
            return Ok(None);
        };
        if !mapped.truthy() {
            return Ok(None);
        }
        let Some(interpro_database) = mapped.as_text() else {
            return Err(format!("the InterPro name for {database} is not text"));
        };
        let default = self
            .vep_to_interpro_rules
            .get("__default__")
            .ok_or_else(|| "KeyError: '__default__'".to_owned())?;
        let rule = self.vep_to_interpro_rules.get(database).unwrap_or(default);
        if !matches!(rule, Json::Object(_)) {
            return Err(format!(
                "the accession rule for {database} is not an object"
            ));
        }
        let rule_type = rule.get("rule_type").and_then(Json::as_text);
        let mut normalized = accession.to_owned();
        match rule_type {
            Some("regex_capture") => {
                let pattern = match rule.get("pattern") {
                    Some(value) if value.truthy() => value
                        .as_text()
                        .ok_or_else(|| format!("the pattern for {database} is not text"))?,
                    _ => "",
                };
                let group = match rule.get("group") {
                    None => 1,
                    Some(Json::Number(number)) if number.fract() == 0.0 => *number as i64,
                    Some(_) => return Err(format!("the group for {database} is not an integer")),
                };
                if !pattern.is_empty() {
                    let compiled = Pattern::compile(pattern)?;
                    let Some(found) = compiled.match_start(accession) else {
                        return Ok(None);
                    };
                    if group < 0 {
                        return Err("IndexError: no such group".to_owned());
                    }
                    if (found.len() as i64) - 1 < group {
                        return Ok(None);
                    }
                    normalized = found[group as usize].clone().unwrap_or_default();
                }
            }
            Some("prefix_add") => {
                normalized = format!("{}{accession}", rule_prefix(rule, database)?);
            }
            Some("prefix_check") => {
                let prefix = rule_prefix(rule, database)?;
                let action = rule.get("action_if_missing").and_then(Json::as_text);
                if !accession.starts_with(prefix) && action == Some("prepend") {
                    normalized = format!("{prefix}{accession}");
                }
            }
            _ => {}
        }
        if normalized.is_empty() {
            return Ok(None);
        }
        Ok(Some((interpro_database.to_owned(), normalized)))
    }

    /// `query_interpro_entry_vep_anno` with the default predicted-domain set.
    pub(crate) fn query_interpro_entry_vep_anno<'a>(
        &self,
        annotation: &str,
        mapping: &'a InterproMapping,
    ) -> Lookup<DomainInfo<'a>> {
        let parts: Vec<&str> = annotation.split('&').collect();
        let mut normalized: BTreeSet<(String, String)> = BTreeSet::new();
        for part in &parts {
            let Some((database, accession)) = part.split_once(':') else {
                continue;
            };
            if let Some(pair) = self.normalize_vep_to_interpro(database, accession)? {
                normalized.insert(pair);
            }
        }
        let vep_domains: Vec<String> = parts
            .iter()
            .filter(|part| PREDICTED_DOMAINS.contains(&part.split(':').next().unwrap_or("")))
            .map(|part| (*part).to_owned())
            .collect();
        if normalized.is_empty() {
            return Ok(DomainInfo {
                vep_domains,
                interpro_entries: None,
            });
        }
        let mut entries: Vec<&InterproEntry> = Vec::new();
        for (database, accession) in &normalized {
            if let Some(found) = mapping.get(database, accession) {
                entries.extend(found.iter());
            }
        }
        Ok(DomainInfo {
            vep_domains,
            interpro_entries: Some(entries),
        })
    }

    /// `interpret_functionality` with the default predicted-domain set.
    pub(crate) fn interpret_functionality(
        &self,
        annotation: &str,
        mapping: &InterproMapping,
    ) -> Lookup<Functionality> {
        let info = self.query_interpro_entry_vep_anno(annotation, mapping)?;
        if info.vep_domains.iter().any(|domain| {
            FUNCTIONAL_PREDICTED_DOMAINS.contains(&domain.split(':').next().unwrap_or(""))
        }) {
            return Ok(Functionality::Functional);
        }
        match info.interpro_entries {
            Some(entries) if !entries.is_empty() => {
                for entry in entries {
                    if is_functional_domain(entry, Some(&self.curated_functionality))? {
                        return Ok(Functionality::Functional);
                    }
                }
                Ok(Functionality::Unknown)
            }
            _ => Ok(Functionality::Missing),
        }
    }
}

fn rule_prefix<'a>(rule: &'a Json, database: &str) -> Lookup<&'a str> {
    match rule.get("prefix") {
        None => Ok(""),
        Some(value) => value
            .as_text()
            .ok_or_else(|| format!("the prefix for {database} is not text")),
    }
}

/// Python's `csv.reader` over the gzip file, keeping column 1 (id) and
/// column 9 (functional status) of every row; the last row of an id wins.
fn read_curated_functional_domains(path: &Path) -> Result<HashMap<String, String>> {
    let mut bytes = Vec::new();
    MultiGzDecoder::new(File::open(path)?).read_to_end(&mut bytes)?;
    let text = String::from_utf8(bytes).map_err(|_| invalid("not UTF-8 text"))?;
    // Text mode turns "\r\n" and a lone "\r" into "\n".
    let text = text.replace("\r\n", "\n").replace('\r', "\n");
    parse_curated_rows(&text)
}

fn parse_curated_rows(text: &str) -> Result<HashMap<String, String>> {
    let mut curated = HashMap::new();
    for (number, mut row) in csv_rows(text).into_iter().enumerate() {
        if row.len() < 9 {
            return Err(invalid(format!(
                "row {} has {} fields, but the id and the status are fields 1 and 9",
                number + 1,
                row.len()
            )));
        }
        let status = std::mem::take(&mut row[8]);
        let id = std::mem::take(&mut row[0]);
        curated.insert(id, status);
    }
    Ok(curated)
}

fn ends_like_dollar(name: &str, suffix: &str) -> bool {
    name.ends_with(suffix)
        || name
            .strip_suffix('\n')
            .is_some_and(|stripped| stripped.ends_with(suffix))
}

/// `DomainNormalizer.is_functional_domain`: the curated status when the entry
/// has one, else the GO terms, then the name, then the entry type.
pub(crate) fn is_functional_domain(
    entry: &InterproEntry,
    curated: Option<&HashMap<String, String>>,
) -> Lookup<bool> {
    let InterproEntry::Usable(entry) = entry else {
        return Err("ValueError: Invalid entry data format".to_owned());
    };
    if let Some(status) = curated.and_then(|table| table.get(&entry.id)) {
        return Ok(status == "yes");
    }
    // Rule 1: GO terms that are lists (a tuple is not a list, so tuples never count).
    if let GoTerms::List(terms) = &entry.go_terms {
        for term in terms {
            let GoTerm::Annotated {
                category,
                description,
            } = term
            else {
                continue;
            };
            if category.as_deref() == Some("molecular_function") {
                return Ok(true);
            }
            if category.as_deref() == Some("biological_process") {
                let Some(description) = description else {
                    return Err("AttributeError: a GO description is not text".to_owned());
                };
                let lowered = description.to_lowercase();
                if FUNCTIONAL_PROCESSES
                    .iter()
                    .any(|process| lowered.contains(*process))
                {
                    return Ok(true);
                }
            }
        }
    }
    // Rule 2: the names.
    let name_text = format!("{} {}", entry.short_name, entry.full_name).to_lowercase();
    if FUNCTIONAL_KEYWORDS
        .iter()
        .any(|keyword| name_text.contains(*keyword))
    {
        return Ok(true);
    }
    if FUNCTIONAL_SUFFIXES
        .iter()
        .any(|suffix| ends_like_dollar(&name_text, *suffix))
        || FUNCTIONAL_PATTERNS
            .iter()
            .any(|pattern| name_text.contains(*pattern))
    {
        return Ok(true);
    }
    // Rule 3: the entry type.
    Ok(matches!(
        entry.kind.as_str(),
        "Active_site" | "Binding_site" | "Conserved_site" | "PTM"
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::interpro_mapping_io::UsableEntry;

    const NORM_MAP: &str = include_str!("../../data/InterPro/InterPro_domain_norm_map.json");

    fn normalizer(curated: &[(&str, &str)]) -> DomainNormalizer {
        let curated = curated
            .iter()
            .map(|(id, status)| ((*id).to_owned(), (*status).to_owned()))
            .collect();
        DomainNormalizer::from_parts(NORM_MAP, curated).unwrap()
    }

    fn entry(id: &str, kind: &str, short: &str, full: &str, go_terms: GoTerms) -> InterproEntry {
        InterproEntry::Usable(UsableEntry {
            id: id.into(),
            kind: kind.into(),
            short_name: short.into(),
            full_name: full.into(),
            go_terms,
        })
    }

    #[test]
    fn accessions_are_normalised_by_the_curated_rules() {
        let norm = normalizer(&[]);
        let pair = |database: &str, accession: &str| {
            norm.normalize_vep_to_interpro(database, accession).unwrap()
        };
        assert_eq!(
            pair("Pfam", "PF00533"),
            Some(("PFAM".into(), "PF00533".into()))
        );
        assert_eq!(
            pair("PANTHER", "PTHR10000:SF12"),
            Some(("PANTHER".into(), "PTHR10000".into()))
        );
        assert_eq!(pair("PANTHER", "other"), None);
        assert_eq!(
            pair("Gene3D", "1.10.8.10"),
            Some(("CATHGENE3D".into(), "G3DSA:1.10.8.10".into()))
        );
        assert_eq!(
            pair("Superfamily", "SSF46955"),
            Some(("SSF".into(), "SSF46955".into()))
        );
        assert_eq!(
            pair("TIGRFAM", "TIGR00001"),
            Some(("NCBIFAM".into(), "TIGR00001".into()))
        );
        assert_eq!(pair("Unheard", "x"), None);
        assert_eq!(pair("Pfam", ""), None);
        // The JSON's "comment" entries are real keys of the mapping.
        assert!(norm.normalize_vep_to_interpro("comment", "x").is_err());
    }

    #[test]
    fn the_regular_expression_subset_matches_like_python() {
        let pattern = Pattern::compile("^(PTHR\\d+)").unwrap();
        assert_eq!(
            pattern.match_start("PTHR10000:SF1"),
            Some(vec![Some("PTHR10000".into()), Some("PTHR10000".into())])
        );
        assert_eq!(pattern.match_start("xPTHR1"), None);
        assert_eq!(pattern.match_start("PTHR"), None);
        let greedy = Pattern::compile("a(.*)b(.?)$").unwrap();
        assert_eq!(
            greedy.match_start("axxbxb"),
            Some(vec![
                Some("axxbxb".into()),
                Some("xxbx".into()),
                Some("".into())
            ])
        );
        assert_eq!(
            greedy.match_start("ab"),
            Some(vec![Some("ab".into()), Some("".into()), Some("".into())])
        );
        for unsupported in ["[a-z]+", "a|b", "(?:a)", "a{2}", "a+?", "\\w", "$a", "(a"] {
            assert!(Pattern::compile(unsupported).is_err(), "{unsupported}");
        }
    }

    #[test]
    fn json_is_read_like_python() {
        let value =
            parse_json("{\"a\": [1, 2.5e1, \"x\\u00e9\\ud83d\\ude00\\n\", null, true], \"a\": {}}")
                .unwrap();
        assert_eq!(value.get("a"), Some(&Json::Object(Vec::new())));
        assert!(parse_json("{\"a\": 01}").is_err());
        assert!(parse_json("{\"a\": 1,}").is_err());
        assert!(parse_json("[1] x").is_err());
        assert!(parse_json("\"tab\there\"").is_err());
        assert!(parse_json("\u{feff}{}").is_err());
        let array = parse_json("[1, 2.5e1, \"x\\u00e9\\ud83d\\ude00\\n\", null, true]").unwrap();
        assert_eq!(
            array,
            Json::Array(vec![
                Json::Number(1.0),
                Json::Number(25.0),
                Json::Text("xé\u{1f600}\n".into()),
                Json::Null,
                Json::Bool(true),
            ])
        );
    }

    #[test]
    fn csv_rows_follow_python() {
        assert_eq!(
            csv_rows("a\tb\"c\t\"d\te\"\"f\"\tx\n\n\"two\nlines\"\ty\nlast"),
            vec![
                vec!["a".to_owned(), "b\"c".into(), "d\te\"f".into(), "x".into()],
                vec![],
                vec!["two\nlines".into(), "y".into()],
                vec!["last".into()],
            ]
        );
        assert_eq!(
            csv_rows("\"ab\"cd\te\t"),
            vec![vec!["abcd".to_owned(), "e".into(), String::new()]]
        );
        // The end of the text closes an unfinished quoted field.
        assert_eq!(
            csv_rows("a\t\"open"),
            vec![vec!["a".to_owned(), "open".into()]]
        );
        let rows = "id\ta\tb\tc\td\te\tf\tg\tstatus\nIPR1\t1\t2\t3\t4\t5\t6\t7\tyes\nIPR1\t1\t2\t3\t4\t5\t6\t7\tno\n";
        let curated = parse_curated_rows(rows).unwrap();
        assert_eq!(curated.get("IPR1").map(String::as_str), Some("no"));
        assert!(parse_curated_rows("a\tb\n").is_err());
        assert!(parse_curated_rows("\n").is_err());
    }

    #[test]
    fn functional_rules_apply_in_order() {
        let go = |category: &str, description: &str| GoTerm::Annotated {
            category: Some(category.into()),
            description: Some(description.into()),
        };
        let plain = |short: &str, full: &str, kind: &str, go_terms: GoTerms| {
            is_functional_domain(&entry("IPR9", kind, short, full, go_terms), None).unwrap()
        };
        // Curated status comes first.
        let curated: HashMap<String, String> = [("IPR9".to_owned(), "no".to_owned())].into();
        assert!(
            !is_functional_domain(
                &entry(
                    "IPR9",
                    "Domain",
                    "Kinase",
                    "Protein kinase",
                    GoTerms::NotAList
                ),
                Some(&curated)
            )
            .unwrap()
        );
        // Rule 1: GO terms count only when each is a list.
        assert!(plain(
            "a",
            "b",
            "Domain",
            GoTerms::List(vec![go("molecular_function", "x")])
        ));
        assert!(plain(
            "a",
            "b",
            "Domain",
            GoTerms::List(vec![go("biological_process", "Signaling pathway")])
        ));
        assert!(!plain(
            "a",
            "b",
            "Domain",
            GoTerms::List(vec![go("biological_process", "pigment")])
        ));
        assert!(!plain(
            "a",
            "b",
            "Domain",
            GoTerms::List(vec![GoTerm::Skipped])
        ));
        // Rule 2: names.
        assert!(plain(
            "Kringle",
            "Protein kinase-like",
            "Domain",
            GoTerms::NotAList
        ));
        assert!(plain("x", "Serine hydrolase", "Domain", GoTerms::NotAList));
        assert!(plain("x", "Homeobox", "Domain", GoTerms::NotAList));
        assert!(plain(
            "x",
            "RNA helicase domain",
            "Domain",
            GoTerms::NotAList
        ));
        assert!(plain(
            "x",
            "a DNA-binding thing",
            "Domain",
            GoTerms::NotAList
        ));
        assert!(!plain("x", "Plain thing", "Domain", GoTerms::NotAList));
        // Rule 3: type.
        assert!(plain("x", "Plain thing", "PTM", GoTerms::NotAList));
        assert!(!plain("x", "Plain thing", "Family", GoTerms::NotAList));
        // An entry of the wrong shape stops the run.
        assert!(is_functional_domain(&InterproEntry::Unusable, None).is_err());
    }

    #[test]
    fn functionality_follows_vep_predictions_then_interpro_entries() {
        let norm = normalizer(&[]);
        let mapping = InterproMapping::from_entries(vec![
            (
                "PFAM",
                "PF1",
                vec![entry(
                    "IPR1",
                    "Domain",
                    "x",
                    "Protein kinase",
                    GoTerms::NotAList,
                )],
            ),
            (
                "PFAM",
                "PF2",
                vec![entry(
                    "IPR2",
                    "Family",
                    "x",
                    "Plain thing",
                    GoTerms::NotAList,
                )],
            ),
            ("PFAM", "PF3", vec![InterproEntry::Unusable]),
        ]);
        let judge = |annotation: &str| norm.interpret_functionality(annotation, &mapping);
        assert_eq!(judge("Pfam:PF1"), Ok(Functionality::Functional));
        assert_eq!(judge("Pfam:PF2"), Ok(Functionality::Unknown));
        assert_eq!(judge("Pfam:PF9"), Ok(Functionality::Missing));
        assert_eq!(judge("something odd"), Ok(Functionality::Missing));
        assert_eq!(
            judge("Transmembrane_helices:TMhelix"),
            Ok(Functionality::Functional)
        );
        assert_eq!(
            judge("Low_complexity_(Seg):seg"),
            Ok(Functionality::Missing)
        );
        assert_eq!(
            judge("Low_complexity_(Seg):seg&Pfam:PF2"),
            Ok(Functionality::Unknown)
        );
        assert_eq!(
            judge("Pfam:PF3"),
            Err("ValueError: Invalid entry data format".to_owned())
        );
        // Entries are visited in sorted order, so the functional one answers first.
        assert_eq!(judge("Pfam:PF1&Pfam:PF3"), Ok(Functionality::Functional));
    }

    #[test]
    fn the_installed_data_files_are_read() {
        let norm = DomainNormalizer::new().unwrap();
        assert_eq!(norm.curated_functionality.len(), 47678);
        assert_eq!(
            norm.curated_functionality
                .get("IPR000001")
                .map(String::as_str),
            Some("yes")
        );
        assert!(
            norm.normalize_vep_to_interpro("Pfam", "PF00001")
                .unwrap()
                .is_some()
        );
    }
}
