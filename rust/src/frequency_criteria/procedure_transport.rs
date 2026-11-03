//! Length-framed TSV for the native procedure/effect interface.
//!
//! There is one named procedure per owned process. A malformed request ends
//! the session. Neither transport nor host errors cause a procedure restart.
//! Objects are owner-local handles; tables enter only after a visited native
//! operation requests them. The protocol does not carry executable programs.

use std::io::{Read, Write};
use super::{Cell, Value};
use super::procedure::{Action, Effect, Invocation, KernelCall, Mode, Object, Procedure, Program, Reply, TypeTest};

pub const PROFILE: &str = "frequency-native-procedure-ef2e701-v1";
const MAX_FRAME: usize = 256 * 1024 * 1024;
type Row = Vec<String>;

fn row(fields: &[&str]) -> Row { fields.iter().map(|value| (*value).into()).collect() }

fn natural(text: &str) -> Result<u64, String> {
    let value = text.parse::<u64>().map_err(|_| "invalid unsigned integer")?;
    if value.to_string() != text { return Err("noncanonical unsigned integer".into()); }
    Ok(value)
}

fn bits(text: &str) -> Result<u64, String> {
    if text.len() != 16 || !text.bytes().all(|value| value.is_ascii_digit() || (b'a'..=b'f').contains(&value)) {
        return Err("invalid canonical float bits".into());
    }
    u64::from_str_radix(text, 16).map_err(|_| "invalid float bits".into())
}

fn records(bytes: &[u8]) -> Result<Vec<Row>, String> {
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').has_headers(false).flexible(false).from_reader(bytes);
    reader.records().map(|record| {
        let record = record.map_err(|error| error.to_string())?;
        if record.len() != 6 { return Err("procedure record width differs".into()); }
        Ok(record.iter().map(str::to_owned).collect())
    }).collect()
}

fn encode(rows: &[Row]) -> Result<Vec<u8>, String> {
    let mut writer = csv::WriterBuilder::new().delimiter(b'\t').has_headers(false).from_writer(Vec::new());
    for row in rows {
        if row.len() != 6 { return Err("internal procedure record width differs".into()); }
        writer.write_record(row).map_err(|error| error.to_string())?;
    }
    writer.into_inner().map_err(|error| error.to_string())
}

fn object(rows: &mut Vec<Row>, role: &str, ordinal: usize, value: &Object) {
    let (tag, value) = match value {
        Object::Handle(value) => ("handle", value.to_string()),
        Object::None => ("none", String::new()),
        Object::Integer(value) => ("int", value.to_string()),
        Object::Float(value) => ("float", format!("{value:016x}")),
        Object::Text(value) => ("str", (*value).into()),
    };
    rows.push(row(&["object", role, &ordinal.to_string(), tag, &value, ""]));
}

fn meta(rows: &mut Vec<Row>, name: &str, value: &str) { rows.push(row(&["meta", name, value, "", "", ""])); }

fn invocation(rows: &mut Vec<Row>, value: &Invocation) {
    match value {
        Invocation::Bound { callable, args } => {
            meta(rows, "invoke", "bound"); object(rows, "callable", 0, callable);
            for (ordinal, value) in args.iter().enumerate() { object(rows, "argument", ordinal, value); }
        }
        Invocation::Binary { operation, left, right } => {
            meta(rows, "invoke", "binary"); meta(rows, "operator", operation.name());
            object(rows, "left", 0, left); object(rows, "right", 0, right);
        }
    }
}

fn call(rows: &mut Vec<Row>, value: &KernelCall) {
    meta(rows, "kernel", value.kernel.name()); invocation(rows, &value.invocation);
    for (ordinal, value) in value.operands.iter().enumerate() { object(rows, "operand", ordinal, value); }
}

fn cell(value: &Cell) -> (&'static str, String) {
    match value {
        Cell::None => ("none", String::new()), Cell::Na => ("na", String::new()),
        Cell::Bool(value) => ("bool", u8::from(*value).to_string()),
        Cell::Int(value) => ("int", value.to_string()),
        Cell::Float(value) => ("float", format!("{:016x}", value.to_bits())),
        Cell::Text(value) => ("str", value.clone()),
    }
}

fn native_value(rows: &mut Vec<Row>, value: &Value) {
    rows.push(row(&["native", &value.kind, &value.dtype, &value.cells.len().to_string(), "", ""]));
    for (ordinal, label) in value.index.iter().enumerate() {
        rows.push(row(&["index", "0", &ordinal.to_string(), "label", label, ""]));
    }
    for (ordinal, value) in value.cells.iter().enumerate() {
        let (tag, value) = cell(value);
        rows.push(row(&["cell", "0", &ordinal.to_string(), tag, &value, ""]));
    }
}

fn action_rows(sequence: u64, action: &Action, native: u64, mature: u64) -> Vec<Row> {
    let mut rows = vec![row(&["action", &sequence.to_string(), "", "", "", PROFILE])];
    match action {
        Action::Effect { token, site, effect } => {
            rows[0][3] = token.to_string();
            rows[0][4] = format!("{}@{}:{}", site.procedure.name(), site.line, site.operation);
            let kind = match effect {
                Effect::IsInstance { object: value, class } => {
                    object(&mut rows, "object", 0, value);
                    meta(&mut rows, "class", match class { TypeTest::Float => "float", TypeTest::String => "str" });
                    "isinstance"
                }
                Effect::BindNumpyWhere => "bind-where",
                Effect::GetItem { owner, key } => { object(&mut rows, "owner", 0, owner); object(&mut rows, "key", 0, key); "getitem" }
                Effect::GetAttribute { owner, name } => { object(&mut rows, "owner", 0, owner); meta(&mut rows, "name", name); "getattr" }
                Effect::PrepareKernel(value) => { call(&mut rows, value); "prepare-kernel" }
                Effect::Materialize { call: value, value: result } => { call(&mut rows, value); native_value(&mut rows, result); "materialize" }
                Effect::InspectBuiltinFloat(value) => { object(&mut rows, "object", 0, value); "inspect-float" }
                Effect::InvokeMature { invocation: value, reason } => { invocation(&mut rows, value); meta(&mut rows, "reason", reason); "mature" }
            };
            rows[0][2] = kind.into();
        }
        Action::Return(value) => { rows[0][2] = "return".into(); object(&mut rows, "result", 0, value); }
        Action::Raise { class, message } => { rows[0][2] = "raise".into(); meta(&mut rows, "class", class); meta(&mut rows, "message", message); }
        Action::HostRaised { exception } => { rows[0][2] = "host-raised".into(); meta(&mut rows, "exception", &exception.to_string()); }
        Action::Unsupported { site, reason } => {
            rows[0][2] = "unsupported".into(); rows[0][4] = format!("{}@{}:{}", site.procedure.name(), site.line, site.operation);
            meta(&mut rows, "reason", reason);
        }
    }
    meta(&mut rows, "native-operations", &native.to_string());
    meta(&mut rows, "mature-operations", &mature.to_string());
    rows
}

fn decode_cell(tag: &str, text: &str) -> Result<Cell, String> {
    match tag {
        "none" if text.is_empty() => Ok(Cell::None), "na" if text.is_empty() => Ok(Cell::Na),
        "bool" if text == "0" || text == "1" => Ok(Cell::Bool(text == "1")),
        "int" => { let value = text.parse::<i64>().map_err(|_| "invalid typed integer")?;
            if value.to_string() != text { return Err("noncanonical typed integer".into()); } Ok(Cell::Int(value)) }
        "float" => Ok(Cell::Float(f64::from_bits(bits(text)?))), "str" => Ok(Cell::Text(text.into())),
        _ => Err("invalid typed cell".into()),
    }
}

fn operands(rows: &[Row], count: usize) -> Result<Vec<Value>, String> {
    if count > 3 { return Err("procedure operand count exceeds named operation".into()); }
    let mut values: Vec<Value> = Vec::new();
    for row in rows {
        let argument = usize::try_from(natural(&row[1])?).map_err(|_| "argument overflow")?;
        match row[0].as_str() {
            "arg" => {
                if argument != values.len() || argument >= count || !row[4].is_empty() || !row[5].is_empty()
                    || !matches!(row[2].as_str(), "scalar" | "series" | "array" | "list") {
                    return Err("invalid operand declaration".into());
                }
                values.push(Value { kind: row[2].clone(), dtype: row[3].clone(), index: Vec::new(), cells: Vec::new() });
            }
            "index" | "cell" => {
                if argument.checked_add(1) != Some(values.len()) || !row[5].is_empty() {
                    return Err("nonconsecutive/interleaved operand".into());
                }
                let value = values.get_mut(argument).ok_or("operand before declaration")?;
                let ordinal = usize::try_from(natural(&row[2])?).map_err(|_| "element overflow")?;
                if row[0] == "index" {
                    if value.kind != "series" || !value.cells.is_empty() || ordinal != value.index.len() || row[3] != "label" {
                        return Err("invalid operand axis".into());
                    }
                    value.index.push(row[4].clone());
                } else {
                    if ordinal != value.cells.len() { return Err("nonconsecutive operand cell".into()); }
                    value.cells.push(decode_cell(&row[3], &row[4])?);
                }
            }
            _ => return Err("unknown typed operand record".into()),
        }
    }
    if values.len() != count || values.iter().any(|value| {
        value.kind == "scalar" && value.cells.len() != 1
            || value.kind == "series" && value.index.len() != value.cells.len()
            || value.kind != "series" && !value.index.is_empty()
    }) { return Err("incomplete typed operand".into()); }
    Ok(values)
}

fn reply(rows: &[Row]) -> Result<Reply, String> {
    let first = rows.first().ok_or("missing procedure reply")?;
    if first[0] != "reply" || first[3..].iter().any(|value| !value.is_empty()) {
        return Err("invalid procedure reply declaration".into());
    }
    if first[1] == "prepared" {
        return Ok(Reply::Prepared(operands(&rows[1..], usize::try_from(natural(&first[2])?).map_err(|_| "operand count overflow")?)?));
    }
    if rows.len() != 1 { return Err("extra procedure reply records".into()); }
    match (first[1].as_str(), first[2].as_str()) {
        ("truth", "0") => Ok(Reply::Truth(false)), ("truth", "1") => Ok(Reply::Truth(true)),
        ("object", value) => Ok(Reply::Object(natural(value)?)),
        ("float", value) => Ok(Reply::BuiltinFloat(bits(value)?)),
        ("other-scalar", "") => Ok(Reply::OtherScalar),
        ("raised", value) => Ok(Reply::Raised(natural(value)?)),
        ("not-represented", value) if !value.is_empty() && value.len() <= 4096 => Ok(Reply::NotRepresented(value.into())),
        _ => Err("unknown or malformed procedure reply".into()),
    }
}

pub struct Session { next: u64, program: Option<Program>, terminal: bool, closed: bool }

impl Default for Session {
    fn default() -> Self { Self { next: 0, program: None, terminal: false, closed: false } }
}

impl Session {
    fn request_inner(&mut self, body: &[u8]) -> Result<Vec<u8>, String> {
        if self.closed { return Err("procedure session already closed".into()); }
        let rows = records(body)?;
        let header = rows.first().ok_or("missing procedure request")?;
        if header[0] != "request" || header[4] != PROFILE || !header[5].is_empty() || natural(&header[1])? != self.next {
            return Err("procedure request profile/sequence differs".into());
        }
        let sequence = self.next;
        self.next = sequence.checked_add(1).ok_or("procedure sequence overflow")?;
        if header[2] == "close" {
            if !self.terminal || !header[3].is_empty() || rows.len() != 1 { return Err("procedure close before terminal action".into()); }
            self.closed = true;
            return encode(&[row(&["action", &sequence.to_string(), "closed", "", "", PROFILE])]);
        }
        let action = match header[2].as_str() {
            "start" => {
                if self.program.is_some() || !header[3].is_empty() || rows.len() != 2 { return Err("duplicate/invalid procedure start".into()); }
                let arguments = &rows[1];
                if arguments[0] != "start" || !arguments[5].is_empty() { return Err("invalid procedure start fields".into()); }
                let mode = match arguments[2].as_str() { "native_math" => Mode::NativeMath, "mixed" => Mode::Mixed,
                    _ => return Err("explicit named execution mode required".into()) };
                let mut program = Program::new(Procedure::named(&arguments[1])?, mode, natural(&arguments[3])?, natural(&arguments[4])?)?;
                let action = program.start()?; self.program = Some(program); action
            }
            "resume" => {
                if self.terminal { return Err("reply after terminal procedure action".into()); }
                self.program.as_mut().ok_or("reply before procedure start")?.resume(natural(&header[3])?, reply(&rows[1..])?)?
            }
            _ => return Err("unknown procedure command".into()),
        };
        self.terminal = !matches!(action, Action::Effect { .. });
        let program = self.program.as_ref().ok_or("missing active procedure")?;
        encode(&action_rows(sequence, &action, program.native_operations, program.mature_operations))
    }

    pub fn request(&mut self, body: &[u8]) -> Result<Vec<u8>, String> {
        let result = self.request_inner(body);
        if result.is_err() { self.closed = true; }
        result
    }
}

fn read_frame(reader: &mut impl Read) -> Result<Option<Vec<u8>>, String> {
    let mut size = [0u8; 4];
    match reader.read(&mut size[..1]) { Ok(0) => return Ok(None), Ok(1) => (), _ => return Err("failed procedure frame prefix".into()) }
    reader.read_exact(&mut size[1..]).map_err(|_| "truncated procedure frame prefix")?;
    let length = u32::from_be_bytes(size) as usize;
    if length == 0 || length > MAX_FRAME { return Err("invalid procedure frame length".into()); }
    let mut body = vec![0; length]; reader.read_exact(&mut body).map_err(|_| "truncated procedure frame body")?;
    Ok(Some(body))
}

pub fn session(reader: &mut impl Read, writer: &mut impl Write) -> Result<(), String> {
    let mut session = Session::default();
    while let Some(body) = read_frame(reader)? {
        let answer = session.request(&body)?;
        if answer.is_empty() || answer.len() > MAX_FRAME { return Err("invalid procedure response length".into()); }
        writer.write_all(&(answer.len() as u32).to_be_bytes()).map_err(|error| error.to_string())?;
        writer.write_all(&answer).map_err(|error| error.to_string())?;
        writer.flush().map_err(|error| error.to_string())?;
        if session.closed { return Ok(()); }
    }
    Err("procedure session ended without terminal close".into())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn request(sequence: u64, command: &str, token: &str, body: &[Row]) -> Vec<u8> {
        let mut rows = vec![row(&["request", &sequence.to_string(), command, token, PROFILE, ""])];
        rows.extend_from_slice(body); encode(&rows).unwrap()
    }
    #[test]
    fn protocol_returns_other_type_none_then_closes() {
        let mut session = Session::default();
        let start = request(0, "start", "", &[row(&["start", "gnomAD_rare_AF", "mixed", "1", "2", ""])]);
        assert_eq!(records(&session.request(&start).unwrap()).unwrap()[0][2], "isinstance");
        let answer = session.request(&request(1, "resume", "0", &[row(&["reply", "truth", "0", "", "", ""])])).unwrap();
        assert_eq!(records(&answer).unwrap()[0][2], "isinstance");
        let answer = session.request(&request(2, "resume", "1", &[row(&["reply", "truth", "0", "", "", ""])])).unwrap();
        assert_eq!(records(&answer).unwrap()[0][2], "return");
        assert_eq!(records(&session.request(&request(3, "close", "", &[])).unwrap()).unwrap()[0][2], "closed");
        assert!(session.request(&request(4, "close", "", &[])).is_err());
    }
    #[test]
    fn invalid_start_and_early_close_are_terminal_faults() {
        let mut session = Session::default();
        assert!(session.request(&request(0, "close", "", &[])).is_err());
        assert!(session.request(&request(0, "start", "", &[row(&["start", "gnomAD_rare_AF", "mixed", "1", "2", ""])])).is_err());
        let mut session = Session::default();
        assert!(session.request(&request(0, "start", "", &[row(&["start", "BS1_criteria", "mixed", "1", "2", ""])])).is_err());
    }
    #[test]
    fn malformed_typed_operands_reject_before_control_resume() {
        assert!(operands(&[row(&["arg", "0", "scalar", "scalar", "", ""])], 1).is_err());
        assert!(operands(&[row(&["arg", "00", "array", "float64", "", ""])], 1).is_err());
        assert!(operands(&[row(&["arg", "0", "series", "float64", "", ""]),
            row(&["cell", "0", "0", "float", "3ff0000000000000", ""])], 1).is_err());
        assert!(bits("7FF8000000000000").is_err());
    }
}
