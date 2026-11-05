//! `pd.read_excel(path, sheet_name=..., engine="openpyxl", dtype=str).fillna("")`.
//!
//! openpyxl's read-only parser decides each cell's value (shared and inline
//! strings, `_cast_number`, booleans, errors, date-styled serials through
//! `from_excel`); pandas' `_convert_cell` and `get_sheet_data` then make
//! integral numbers integers, trim trailing empty cells and rows and pad the
//! rows; the Python text parser names the columns from the first row and
//! turns pandas' missing-value spellings into empty strings. ISO `t="d"`
//! cells are refused.
use super::*;
use crate::interpro_mapping_builder::Element;
use std::cell::RefCell;
use std::rc::Rc;

#[derive(Clone, Debug, Default)]
struct Node {
    name: String,
    attributes: Vec<(String, String)>,
    text: String,
    children: Vec<Node>,
}

impl Node {
    fn local(&self) -> &str {
        self.name.rsplit('}').next().unwrap_or(&self.name)
    }
    fn plain(&self, name: &str) -> Option<&str> {
        self.attributes
            .iter()
            .find(|(key, _)| key == name)
            .map(|(_, value)| value.as_str())
    }
    fn children(&self, name: &'static str) -> impl Iterator<Item = &Node> {
        self.children.iter().filter(move |child| child.local() == name)
    }
    fn child(&self, name: &'static str) -> Option<&Node> {
        self.children(name).next()
    }
}

fn copy(element: &Element) -> Node {
    Node {
        name: element.name.clone(),
        attributes: element.ordered_attributes.clone(),
        text: element.text.clone(),
        children: element
            .children
            .iter()
            .map(|child| copy(&child.borrow()))
            .collect(),
    }
}

fn bad(message: impl Into<String>) -> Box<dyn std::error::Error> {
    error("ValueError", message)
}

struct Book {
    archive: zip::ZipArchive<File>,
}

impl Book {
    fn names(&self) -> Vec<String> {
        self.archive.file_names().map(str::to_owned).collect()
    }
    fn part(&mut self, name: &str) -> Result<Option<Node>> {
        let mut member = match self.archive.by_name(name) {
            Ok(member) => member,
            Err(zip::result::ZipError::FileNotFound) => return Ok(None),
            Err(other) => return Err(error("BadZipFile", other.to_string())),
        };
        let mut bytes = Vec::new();
        member
            .read_to_end(&mut bytes)
            .map_err(|e| error("BadZipFile", e.to_string()))?;
        let root = Rc::new(RefCell::new(None));
        let sink = root.clone();
        crate::interpro_mapping_builder::visit_document(&bytes[..], false, move |element| {
            *sink.borrow_mut() = Some(copy(element));
        })?;
        let node = root.borrow_mut().take();
        Ok(node)
    }
}

/// `posixpath.normpath(posixpath.join(base_dir, target))` without a leading `/`.
fn join_part(base_dir: &str, target: &str) -> String {
    let joined = if let Some(absolute) = target.strip_prefix('/') {
        absolute.to_owned()
    } else if base_dir.is_empty() {
        target.to_owned()
    } else {
        format!("{base_dir}/{target}")
    };
    let mut parts: Vec<&str> = Vec::new();
    for part in joined.split('/') {
        match part {
            "" | "." => {}
            ".." => {
                parts.pop();
            }
            part => parts.push(part),
        }
    }
    parts.join("/")
}

fn rels_path(part: &str) -> String {
    match part.rsplit_once('/') {
        Some((folder, name)) => format!("{folder}/_rels/{name}.rels"),
        None => format!("_rels/{part}.rels"),
    }
}

const WORKBOOK_TYPES: [&str; 4] = [
    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml",
    "application/vnd.ms-excel.sheet.macroEnabled.main+xml",
    "application/vnd.openxmlformats-officedocument.spreadsheetml.template.main+xml",
    "application/vnd.ms-excel.template.macroEnabled.main+xml",
];
const SHARED_STRINGS: &str =
    "application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml";

/// openpyxl's `Text.from_tree(node).content`.
fn rich_text(node: &Node) -> String {
    let mut out = String::new();
    if let Some(plain) = node.child("t") {
        out.push_str(&plain.text);
    }
    for run in node.children("r") {
        if let Some(part) = run.child("t") {
            out.push_str(&part.text);
        }
    }
    out
}

const BUILTIN: [(u32, &str); 43] = [
    (0, "General"),
    (1, "0"),
    (2, "0.00"),
    (3, "#,##0"),
    (4, "#,##0.00"),
    (5, "\"$\"#,##0_);(\"$\"#,##0)"),
    (6, "\"$\"#,##0_);[Red](\"$\"#,##0)"),
    (7, "\"$\"#,##0.00_);(\"$\"#,##0.00)"),
    (8, "\"$\"#,##0.00_);[Red](\"$\"#,##0.00)"),
    (9, "0%"),
    (10, "0.00%"),
    (11, "0.00E+00"),
    (12, "# ?/?"),
    (13, "# ??/??"),
    (14, "mm-dd-yy"),
    (15, "d-mmm-yy"),
    (16, "d-mmm"),
    (17, "mmm-yy"),
    (18, "h:mm AM/PM"),
    (19, "h:mm:ss AM/PM"),
    (20, "h:mm"),
    (21, "h:mm:ss"),
    (22, "m/d/yy h:mm"),
    (37, "#,##0_);(#,##0)"),
    (38, "#,##0_);[Red](#,##0)"),
    (39, "#,##0.00_);(#,##0.00)"),
    (40, "#,##0.00_);[Red](#,##0.00)"),
    (41, r#"_(* #,##0_);_(* \(#,##0\);_(* "-"_);_(@_)"#),
    (42, r#"_("$"* #,##0_);_("$"* \(#,##0\);_("$"* "-"_);_(@_)"#),
    (43, r#"_(* #,##0.00_);_(* \(#,##0.00\);_(* "-"??_);_(@_)"#),
    (44, r#"_("$"* #,##0.00_)_("$"* \(#,##0.00\)_("$"* "-"??_)_(@_)"#),
    (45, "mm:ss"),
    (46, "[h]:mm:ss"),
    (47, "mmss.0"),
    (48, "##0.0E+0"),
    (49, "@"),
    (23, "General"),
    (24, "General"),
    (25, "General"),
    (26, "General"),
    (27, "General"),
    (28, "General"),
    (29, "General"),
];

/// openpyxl's `is_date_format`.
fn is_date_format(code: &str) -> bool {
    let first = code.split(';').next().unwrap_or("");
    let chars: Vec<char> = first.chars().collect();
    let mut kept = Vec::new();
    let mut at = 0;
    while at < chars.len() {
        if chars[at] == '"' {
            if let Some(end) = chars[at + 1..].iter().position(|c| *c == '"') {
                at += end + 2;
                continue;
            }
        }
        if chars[at] == '[' {
            if let Some(end) = chars[at + 1..].iter().position(|c| *c == ']') {
                let inner: String = chars[at + 1..at + 1 + end].iter().collect();
                if !["h", "hh", "m", "mm", "s", "ss"].contains(&inner.as_str()) {
                    at += end + 2;
                    continue;
                }
            }
        }
        kept.push(chars[at]);
        at += 1;
    }
    kept.iter().enumerate().any(|(index, c)| {
        "dmhysDMHYS".contains(*c) && !(index > 0 && matches!(kept[index - 1], '_' | '\\'))
    })
}

/// openpyxl's `is_timedelta_format`.
fn is_timedelta_format(code: &str) -> bool {
    let first = code.split(';').next().unwrap_or("").to_lowercase();
    ["[h]", "[hh]", "[m]", "[mm]", "[s]", "[ss]"]
        .iter()
        .any(|token| first.contains(token))
}

/// A cell after openpyxl and pandas' `_convert_cell`.
#[derive(Clone, Debug)]
enum Value {
    Empty,
    Missing,
    Text(String),
    Int(i128),
    Float(f64),
    Bool(bool),
    Rendered(String),
}

impl Value {
    fn blank(&self) -> bool {
        matches!(self, Value::Empty) || matches!(self, Value::Text(text) if text.is_empty())
    }
}

fn time_text(hours: i64, minutes: i64, seconds: i64, micros: i64) -> String {
    if micros != 0 {
        format!("{hours:02}:{minutes:02}:{seconds:02}.{micros:06}")
    } else {
        format!("{hours:02}:{minutes:02}:{seconds:02}")
    }
}

/// `str(from_excel(value, epoch, timedelta))`; `None` for an overflow.
fn excel_date(value: f64, mac: bool, timedelta: bool) -> Option<String> {
    if timedelta {
        let mut micros = (value * 86_400_000_000.0).round_ties_even() as i128;
        if micros.rem_euclid(1_000_000) != 0 {
            let seconds = (micros as f64 / 1e6).floor() as i128;
            let rest = micros - seconds * 1_000_000;
            micros = seconds * 1_000_000 + ((rest as f64 / 1000.0).round_ties_even() as i128) * 1000;
        }
        let days = micros.div_euclid(86_400_000_000);
        let within = micros.rem_euclid(86_400_000_000);
        let clock = time_text(
            (within / 3_600_000_000) as i64,
            (within / 60_000_000 % 60) as i64,
            (within / 1_000_000 % 60) as i64,
            (within % 1_000_000) as i64,
        );
        let clock = clock.trim_start_matches('0');
        let clock = if clock.starts_with(':') { format!("0{clock}") } else { clock.to_owned() };
        return Some(if days != 0 {
            format!("{days} day{}, {clock}", if days.abs() != 1 { "s" } else { "" })
        } else {
            clock
        });
    }
    let day = value.floor();
    let fraction = value - day;
    let millis = (fraction * 86_400.0 * 1000.0).round_ties_even() as i64;
    if (0.0..1.0).contains(&value) && millis < 86_400_000 {
        return Some(time_text(
            millis / 3_600_000,
            millis / 60_000 % 60,
            millis / 1000 % 60,
            millis % 1000 * 1000,
        ));
    }
    let mut day = day as i64;
    if !mac && value > 0.0 && value < 60.0 {
        day += 1;
    }
    let epoch = if mac {
        days_from_civil(1904, 1, 1)
    } else {
        days_from_civil(1899, 12, 30)
    };
    let total = (epoch + day) as i128 * 86_400_000 + i128::from(millis);
    let days = total.div_euclid(86_400_000) as i64;
    let within = total.rem_euclid(86_400_000) as i64;
    let (year, month, date) = {
        let z = days + 719_468;
        let era = z.div_euclid(146_097);
        let doe = z - era * 146_097;
        let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
        let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
        let mp = (5 * doy + 2) / 153;
        let d = doy - (153 * mp + 2) / 5 + 1;
        let m = if mp < 10 { mp + 3 } else { mp - 9 };
        (yoe + era * 400 + i64::from(m <= 2), m, d)
    };
    if !(1..=9999).contains(&year) {
        return None;
    }
    Some(format!(
        "{year:04}-{month:02}-{date:02} {}",
        time_text(
            within / 3_600_000,
            within / 60_000 % 60,
            within / 1000 % 60,
            within % 1000 * 1000
        )
    ))
}

/// openpyxl's `coordinate_to_tuple`.
fn coordinate(text: &str) -> Result<(i128, i128)> {
    let body = text.replace('$', "");
    let split = body.find(|c: char| c.is_ascii_digit()).unwrap_or(body.len());
    let (letters, digits) = body.split_at(split);
    if letters.is_empty()
        || letters.len() > 3
        || !letters.chars().all(|c| c.is_ascii_alphabetic())
        || digits.is_empty()
        || !digits.chars().all(|c| c.is_ascii_digit())
    {
        return Err(error(
            "CellCoordinatesException",
            format!("Invalid cell coordinates ({text})"),
        ));
    }
    let column = letters
        .chars()
        .fold(0i128, |total, c| total * 26 + i128::from(c.to_ascii_uppercase() as u8 - b'A' + 1));
    Ok((digits.parse().unwrap_or(0), column))
}

fn python_int(text: &str) -> Result<i128> {
    crate::python_text::python_int(text).ok_or_else(|| {
        bad(format!(
            "invalid literal for int() with base 10: {}",
            crate::python_path::repr(text)
        ))
    })
}

struct Styles {
    dates: std::collections::HashSet<usize>,
    timedeltas: std::collections::HashSet<usize>,
}

fn styles(book: &mut Book) -> Result<Styles> {
    let mut result = Styles {
        dates: Default::default(),
        timedeltas: Default::default(),
    };
    let Some(root) = book.part("xl/styles.xml")? else {
        return Ok(result);
    };
    let mut custom: HashMap<u32, String> = HashMap::new();
    if let Some(formats) = root.child("numFmts") {
        for format in formats.children("numFmt") {
            if let (Some(id), Some(code)) = (format.plain("numFmtId"), format.plain("formatCode")) {
                custom.insert(id.parse().unwrap_or(0), code.to_owned());
            }
        }
    }
    if let Some(xfs) = root.child("cellXfs") {
        for (index, xf) in xfs.children("xf").enumerate() {
            let id: u32 = xf.plain("numFmtId").and_then(|v| v.parse().ok()).unwrap_or(0);
            let code = custom.get(&id).cloned().unwrap_or_else(|| {
                BUILTIN
                    .iter()
                    .find(|(builtin, _)| *builtin == id)
                    .map_or("General", |(_, code)| *code)
                    .to_owned()
            });
            if is_date_format(&code) {
                result.dates.insert(index);
            }
            if is_timedelta_format(&code) {
                result.timedeltas.insert(index);
            }
        }
    }
    Ok(result)
}

/// The DataFrame's column names and its text cells after `fillna("")`.
pub(super) fn read_sheet(path: &Path, sheet: &str) -> Result<(Vec<String>, Vec<Vec<String>>)> {
    let file = File::open(path).map_err(|e| os_error(e, path, None))?;
    let archive = zip::ZipArchive::new(file)
        .map_err(|_| error("BadZipFile", "File is not a zip file"))?;
    let mut book = Book { archive };
    let manifest = book
        .part("[Content_Types].xml")?
        .ok_or_else(|| error("KeyError", "\"There is no item named '[Content_Types].xml' in the archive\""))?;
    let override_part = |types: &[&str]| {
        manifest
            .children("Override")
            .find(|node| node.plain("ContentType").is_some_and(|t| types.contains(&t)))
            .and_then(|node| node.plain("PartName"))
            .map(|name| name.trim_start_matches('/').to_owned())
    };
    let workbook_part = override_part(&WORKBOOK_TYPES)
        .ok_or_else(|| error("OSError", "File contains no valid workbook part"))?;
    let strings: Vec<String> = match override_part(&[SHARED_STRINGS]) {
        Some(part) => match book.part(&part)? {
            Some(root) => root
                .children("si")
                .map(|item| rich_text(item).replace("x005F_", ""))
                .collect(),
            None => vec![],
        },
        None => vec![],
    };
    let workbook = book
        .part(&workbook_part)?
        .ok_or_else(|| error("KeyError", format!("There is no item named '{workbook_part}' in the archive")))?;
    let mac = workbook
        .child("workbookPr")
        .and_then(|p| p.plain("date1904"))
        .is_some_and(|v| matches!(v, "1" | "true"));
    let base_dir = workbook_part.rsplit_once('/').map_or("", |(folder, _)| folder).to_owned();
    let rels = book.part(&rels_path(&workbook_part))?;
    let names = book.names();
    let target = workbook
        .child("sheets")
        .into_iter()
        .flat_map(|sheets| sheets.children("sheet"))
        .filter_map(|node| {
            let id = node
                .attributes
                .iter()
                .find(|(key, _)| key.ends_with("}id"))
                .map(|(_, value)| value.clone())?;
            let target = rels.as_ref()?.children("Relationship").find(|r| r.plain("Id") == Some(id.as_str()))?;
            let part = join_part(&base_dir, target.plain("Target")?);
            names.contains(&part).then(|| (node.plain("name").unwrap_or("").to_owned(), part))
        })
        .find(|(name, _)| name == sheet)
        .map(|(_, part)| part)
        .ok_or_else(|| bad(format!("Worksheet named '{sheet}' not found")))?;
    let style = styles(&mut book)?;
    let root = book
        .part(&target)?
        .ok_or_else(|| error("KeyError", format!("There is no item named '{target}' in the archive")))?;
    let mut rows: Vec<Vec<Value>> = Vec::new();
    let mut counter: i128 = 1;
    let mut row_counter: i128 = 0;
    if let Some(data) = root.child("sheetData") {
        for row in data.children("row") {
            row_counter = match row.plain("r") {
                Some(text) => match crate::python_text::python_int(text) {
                    Some(number) => number,
                    None => {
                        let value = crate::python_text::parse_float(text)
                            .ok_or_else(|| bad(format!("could not convert string to float: '{text}'")))?;
                        if value.fract() != 0.0 {
                            return Err(bad(format!("{text} is not a valid row number")));
                        }
                        value as i128
                    }
                },
                None => row_counter + 1,
            };
            let mut column_counter: i128 = 0;
            let mut cells: Vec<(i128, Value)> = Vec::new();
            for cell in &row.children {
                let (_, column) = match cell.plain("r") {
                    Some(text) => {
                        let found = coordinate(text)?;
                        column_counter = found.1;
                        found
                    }
                    None => {
                        column_counter += 1;
                        (row_counter, column_counter)
                    }
                };
                let kind = cell.plain("t").unwrap_or("n");
                let style_id: usize = cell.plain("s").and_then(|v| v.parse().ok()).unwrap_or(0);
                let raw = if kind == "inlineStr" {
                    None
                } else {
                    cell.child("v").map(|v| v.text.clone()).filter(|v| !v.is_empty())
                };
                let value = match (kind, raw) {
                    ("inlineStr", _) => cell.child("is").map_or(Value::Empty, |node| Value::Text(rich_text(node))),
                    (_, None) => Value::Empty,
                    ("n", Some(raw)) => {
                        let number = if raw.contains(['.', 'E', 'e']) {
                            Value::Float(crate::python_text::parse_float(&raw).ok_or_else(|| {
                                bad(format!("could not convert string to float: {}", crate::python_path::repr(&raw)))
                            })?)
                        } else {
                            Value::Int(python_int(&raw)?)
                        };
                        if style.dates.contains(&style_id) {
                            let serial = match number {
                                Value::Int(n) => n as f64,
                                Value::Float(f) => f,
                                _ => 0.0,
                            };
                            match excel_date(serial, mac, style.timedeltas.contains(&style_id)) {
                                Some(text) => Value::Rendered(text),
                                None => Value::Missing,
                            }
                        } else {
                            match number {
                                Value::Float(f) if !f.is_finite() => {
                                    return Err(error(
                                        if f.is_nan() { "ValueError" } else { "OverflowError" },
                                        "cannot convert float to integer",
                                    ));
                                }
                                Value::Float(f) if f.fract() == 0.0 => {
                                    Value::Rendered(format!("{f:.0}").replace("-0", "0"))
                                }
                                other => other,
                            }
                        }
                    }
                    ("s", Some(raw)) => {
                        let index = python_int(&raw)?;
                        Value::Text(
                            usize::try_from(index)
                                .ok()
                                .and_then(|i| strings.get(i))
                                .cloned()
                                .ok_or_else(|| error("IndexError", "list index out of range"))?,
                        )
                    }
                    ("b", Some(raw)) => Value::Bool(python_int(&raw)? != 0),
                    ("e", Some(_)) => Value::Missing,
                    ("d", Some(raw)) => {
                        return Err(bad(format!("ISO date cells are not supported: {raw}")));
                    }
                    (_, Some(raw)) => Value::Text(raw),
                };
                cells.push((column, value));
            }
            for _ in counter..row_counter {
                counter += 1;
                rows.push(vec![]);
            }
            if counter <= row_counter {
                counter += 1;
                let width = cells.last().map_or(0, |(column, _)| *column);
                let mut values = vec![Value::Empty; usize::try_from(width.max(0)).unwrap_or(0)];
                for (column, value) in cells {
                    if (1..=width).contains(&column) {
                        values[column as usize - 1] = value;
                    }
                }
                rows.push(values);
            }
        }
    }
    let mut last = None;
    for (index, row) in rows.iter_mut().enumerate() {
        while row.last().is_some_and(Value::blank) {
            row.pop();
        }
        if !row.is_empty() {
            last = Some(index);
        }
    }
    rows.truncate(last.map_or(0, |index| index + 1));
    let Some(width) = rows.iter().map(Vec::len).max() else {
        return Ok((vec![], vec![]));
    };
    for row in &mut rows {
        row.resize(width, Value::Empty);
    }
    let header: Vec<String> = rows[0]
        .iter()
        .map(|value| match value {
            Value::Empty | Value::Missing => String::new(),
            Value::Text(text) | Value::Rendered(text) => text.clone(),
            Value::Int(number) => number.to_string(),
            Value::Float(number) => crate::python_text::float_repr(*number),
            Value::Bool(flag) => (if *flag { "True" } else { "False" }).into(),
        })
        .collect();
    let names = crate::pandas_table::column_names(header);
    let mut cells: Vec<Vec<String>> = rows[1..]
        .iter()
        .map(|row| {
            row.iter()
                .map(|value| match value {
                    Value::Empty | Value::Missing => String::new(),
                    Value::Text(text) if crate::pandas_table::is_na(text) => String::new(),
                    Value::Text(text) | Value::Rendered(text) => text.clone(),
                    Value::Int(number) => number.to_string(),
                    Value::Float(number) => crate::python_text::float_repr(*number),
                    Value::Bool(flag) => (if *flag { "True" } else { "False" }).into(),
                })
                .collect()
        })
        .collect();
    if width == 1 {
        let original = rows[1..].to_vec();
        cells = cells
            .into_iter()
            .zip(original)
            .filter(|(_, row)| !matches!(&row[0], Value::Empty) && !matches!(&row[0], Value::Text(text) if crate::python_text::python_strip(text).is_empty()))
            .map(|(cells, _)| cells)
            .collect();
    }
    Ok((names, cells))
}
