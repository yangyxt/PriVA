//! Rust-owned named procedures with a sequential mature-effect host.
//!
//! The single child drives the procedure. Every visited host effect has one
//! sequence number; scientific exceptions abort with their host-held identity.
//! Bad framing/foreign handles/dead streams are transport faults, never a
//! request to execute an original criterion or replay completed work.
use super::procedures::{self, Call, Host, Object, Site, Target};
use super::transport::{encode_cell, read_frame, request, response};
use super::{Cell, Value};
use std::io::{Read, Write};

pub const PROFILE: &str = "insilico-procedures-v2";
const MAX_FRAME: usize = 256 * 1024 * 1024;

#[derive(Debug)]
enum Failure {
    Source,
    Transport(String),
}

type Rows = Vec<Vec<String>>;

fn row(fields: [&str; 6]) -> Vec<String> {
    fields.into_iter().map(str::to_owned).collect()
}

fn parse(bytes: &[u8]) -> Result<Rows, Failure> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .flexible(false)
        .from_reader(bytes);
    let rows: Vec<csv::StringRecord> = reader
        .records()
        .collect::<Result<_, _>>()
        .map_err(|error| Failure::Transport(error.to_string()))?;
    if rows.is_empty() || rows.iter().any(|record| record.len() != 6) {
        return Err(Failure::Transport("procedure frame schema".into()));
    }
    Ok(rows
        .iter()
        .map(|record| record.iter().map(str::to_owned).collect())
        .collect())
}

fn write_frame(writer: &mut impl Write, rows: Rows) -> Result<(), Failure> {
    let mut csv = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(Vec::new());
    for record in rows {
        csv.write_record(record)
            .map_err(|error| Failure::Transport(error.to_string()))?;
    }
    let bytes = csv
        .into_inner()
        .map_err(|error| Failure::Transport(error.to_string()))?;
    if bytes.is_empty() || bytes.len() > MAX_FRAME {
        return Err(Failure::Transport("procedure frame length".into()));
    }
    writer
        .write_all(&(bytes.len() as u32).to_be_bytes())
        .and_then(|()| writer.write_all(&bytes))
        .and_then(|()| writer.flush())
        .map_err(|error| Failure::Transport(error.to_string()))
}

fn handle(value: &str) -> Result<Object, Failure> {
    if value.is_empty() || !value.bytes().all(|byte| byte.is_ascii_digit()) {
        return Err(Failure::Transport("invalid object handle".into()));
    }
    value
        .parse()
        .map_err(|_| Failure::Transport("object handle overflow".into()))
}

struct Remote<'a, R: Read, W: Write> {
    reader: &'a mut R,
    writer: &'a mut W,
    sequence: u64,
}

impl<R: Read, W: Write> Remote<'_, R, W> {
    fn exchange(&mut self, site: Site, action: &str, mut details: Rows) -> Result<Rows, Failure> {
        let sequence = self.sequence.to_string();
        let line = site.line.to_string();
        let mut rows = vec![row([
            "host",
            &sequence,
            action,
            site.function,
            &line,
            PROFILE,
        ])];
        rows.append(&mut details);
        write_frame(self.writer, rows)?;
        let bytes = read_frame(self.reader)
            .map_err(Failure::Transport)?
            .ok_or_else(|| Failure::Transport("host ended before reply".into()))?;
        let rows = parse(&bytes)?;
        let header = &rows[0];
        if header[0] != "host-result"
            || header[1] != sequence
            || header[3] != PROFILE
            || !header[4].is_empty()
            || !header[5].is_empty()
        {
            return Err(Failure::Transport("host reply sequence/profile".into()));
        }
        self.sequence = self
            .sequence
            .checked_add(1)
            .ok_or_else(|| Failure::Transport("sequence overflow".into()))?;
        if header[2] == "source-error" {
            if rows.len() != 1 {
                return Err(Failure::Transport(
                    "source exception payload forbidden".into(),
                ));
            }
            return Err(Failure::Source);
        }
        Ok(rows)
    }

    fn object(&mut self, site: Site, action: &str, details: Rows) -> Result<Object, Failure> {
        let rows = self.exchange(site, action, details)?;
        if rows.len() != 2
            || rows[0][2] != "object"
            || rows[1][0] != "object"
            || rows[1][1] != "0"
            || rows[1][3..].iter().any(|value| !value.is_empty())
        {
            return Err(Failure::Transport("object reply schema".into()));
        }
        handle(&rows[1][2])
    }

    fn boolean(&mut self, site: Site, action: &str, object: Object) -> Result<bool, Failure> {
        let rows = self.exchange(
            site,
            action,
            vec![row(["object", "0", &object.to_string(), "", "", ""])],
        )?;
        if rows.len() != 2
            || rows[0][2] != "boolean"
            || rows[1][0] != "boolean"
            || rows[1][1] != "0"
            || !["0", "1"].contains(&rows[1][2].as_str())
            || rows[1][3..].iter().any(|value| !value.is_empty())
        {
            return Err(Failure::Transport("boolean reply schema".into()));
        }
        Ok(rows[1][2] == "1")
    }

    fn call_rows(call: &Call) -> Rows {
        let (kind, function, receiver) = match &call.target {
            Target::Binary(op) => (*op, String::new(), String::new()),
            Target::Invert => ("invert", String::new(), String::new()),
            Target::Assignment => ("assignment", String::new(), String::new()),
            Target::Bound { function, receiver } => (
                "bound",
                function.to_string(),
                receiver.map(|value| value.to_string()).unwrap_or_default(),
            ),
        };
        let mut rows = vec![row([
            "call",
            call.operation,
            kind,
            &function,
            &receiver,
            "",
        ])];
        for (ordinal, value) in call.args.iter().enumerate() {
            rows.push(row([
                "argument",
                &ordinal.to_string(),
                &value.to_string(),
                "",
                "",
                "",
            ]));
        }
        for (ordinal, (name, value)) in call.keywords.iter().enumerate() {
            rows.push(row([
                "keyword",
                &ordinal.to_string(),
                name,
                &value.to_string(),
                "",
                "",
            ]));
        }
        rows
    }
}

impl<R: Read, W: Write> Host for Remote<'_, R, W> {
    type Error = Failure;

    fn literal(&mut self, site: Site, value: Cell) -> Result<Object, Failure> {
        let (tag, raw) = encode_cell(&value);
        self.object(
            site,
            "literal",
            vec![row(["literal", tag, &raw, "", "", ""])],
        )
    }
    fn column(&mut self, site: Site, frame: Object, key: &'static str) -> Result<Object, Failure> {
        self.object(
            site,
            "column",
            vec![
                row(["object", "0", &frame.to_string(), "", "", ""]),
                row(["text", "0", key, "", "", ""]),
            ],
        )
    }
    fn attribute(
        &mut self,
        site: Site,
        object: Object,
        name: &'static str,
    ) -> Result<Object, Failure> {
        self.object(
            site,
            "attribute",
            vec![
                row(["object", "0", &object.to_string(), "", "", ""]),
                row(["text", "0", name, "", "", ""]),
            ],
        )
    }
    fn global(
        &mut self,
        site: Site,
        module: Option<&'static str>,
        name: &'static str,
    ) -> Result<Object, Failure> {
        self.object(
            site,
            "global",
            vec![
                row(["text", "0", module.unwrap_or(""), "", "", ""]),
                row(["text", "1", name, "", "", ""]),
            ],
        )
    }
    fn length(&mut self, site: Site, object: Object) -> Result<Object, Failure> {
        self.object(
            site,
            "length",
            vec![row(["object", "0", &object.to_string(), "", "", ""])],
        )
    }
    fn effect(
        &mut self,
        site: Site,
        action: &'static str,
        args: &[Object],
        text: &[&str],
    ) -> Result<Object, Failure> {
        let mut rows = vec![row(["effect", action, "", "", "", ""])];
        for (ordinal, object) in args.iter().enumerate() {
            rows.push(row([
                "object",
                &ordinal.to_string(),
                &object.to_string(),
                "",
                "",
                "",
            ]));
        }
        for (ordinal, value) in text.iter().enumerate() {
            rows.push(row(["text", &ordinal.to_string(), value, "", "", ""]));
        }
        self.object(site, "effect", rows)
    }
    fn truth(&mut self, site: Site, object: Object) -> Result<bool, Failure> {
        self.boolean(site, "truth", object)
    }
    fn is_none(&mut self, site: Site, object: Object) -> Result<bool, Failure> {
        self.boolean(site, "is_none", object)
    }
    fn release(&mut self, site: Site, objects: &[Object]) -> Result<(), Failure> {
        let rows = objects
            .iter()
            .enumerate()
            .map(|(ordinal, object)| {
                row([
                    "object",
                    &ordinal.to_string(),
                    &object.to_string(),
                    "",
                    "",
                    "",
                ])
            })
            .collect();
        let reply = self.exchange(site, "release", rows)?;
        if reply.len() != 1 || reply[0][2] != "released" {
            return Err(Failure::Transport("release acknowledgment differs".into()));
        }
        Ok(())
    }
    fn prepare(&mut self, call: &Call) -> Result<Option<Vec<Value>>, Failure> {
        let sequence = self.sequence;
        let rows = self.exchange(call.site, "prepare", Self::call_rows(call))?;
        if rows[0][2] == "mature" && rows.len() == 1 {
            return Ok(None);
        }
        if rows[0][2] != "prepared"
            || rows.len() != 2
            || rows[1][0] != "operands"
            || rows[1][1] != "0"
            || rows[1][3..].iter().any(|value| !value.is_empty())
        {
            return Err(Failure::Transport("prepared operands schema".into()));
        }
        let (actual_sequence, operation, values) =
            request(rows[1][2].as_bytes()).map_err(Failure::Transport)?;
        if actual_sequence != sequence || operation != call.operation {
            return Err(Failure::Transport(
                "prepared operands operation/sequence".into(),
            ));
        }
        Ok(Some(values))
    }
    fn commit_native(&mut self, call: &Call, value: Value) -> Result<Object, Failure> {
        let mut rows = Self::call_rows(call);
        let bytes = response(self.sequence, Ok(value)).map_err(Failure::Transport)?;
        let text = String::from_utf8(bytes)
            .map_err(|_| Failure::Transport("native result UTF-8".into()))?;
        rows.push(row(["native", "0", &text, "", "", ""]));
        self.object(call.site, "commit", rows)
    }
    fn mature(&mut self, call: &Call) -> Result<Object, Failure> {
        self.object(call.site, "mature", Self::call_rows(call))
    }
}

fn execute(reader: &mut impl Read, writer: &mut impl Write) -> Result<(), Failure> {
    let bytes = read_frame(reader)
        .map_err(Failure::Transport)?
        .ok_or_else(|| Failure::Transport("missing named procedure request".into()))?;
    let rows = parse(&bytes)?;
    let header = &rows[0];
    if rows.len() != 1 || header[0] != "run" || header[5] != PROFILE {
        return Err(Failure::Transport("named procedure request schema".into()));
    }
    let frame = handle(&header[2])?;
    let pvs1 = handle(&header[3])?;
    let confidence = handle(&header[4])?;
    let mut host = Remote {
        reader,
        writer,
        sequence: 0,
    };
    let result = match header[1].as_str() {
        "BP7_criteria" => procedures::bp7_criteria(&mut host, frame).map(|value| vec![value]),
        "PP3_BP4_criteria" => procedures::pp3_bp4_criteria(&mut host, frame, pvs1, confidence)
            .map(|(left, right)| vec![left, right]),
        _ => return Err(Failure::Transport("unknown named procedure".into())),
    };
    let mut rows = match result {
        Ok(values) => {
            let mut rows = vec![row([
                "done",
                &host.sequence.to_string(),
                "ok",
                PROFILE,
                procedures::SOURCE_SHA256,
                "",
            ])];
            for (ordinal, value) in values.iter().enumerate() {
                rows.push(row([
                    "object",
                    &ordinal.to_string(),
                    &value.to_string(),
                    "",
                    "",
                    "",
                ]));
            }
            rows
        }
        Err(Failure::Source) => vec![row([
            "done",
            &host.sequence.to_string(),
            "source-error",
            PROFILE,
            procedures::SOURCE_SHA256,
            "",
        ])],
        Err(failure) => return Err(failure),
    };
    write_frame(host.writer, std::mem::take(&mut rows))
}

pub(crate) fn session(reader: &mut impl Read, writer: &mut impl Write) -> Result<(), String> {
    match execute(reader, writer) {
        Ok(()) => Ok(()),
        Err(Failure::Transport(message)) => Err(message),
        Err(Failure::Source) => Err("unacknowledged host exception".into()),
    }
}

#[cfg(test)]
mod tests {
    use super::super::transport::decode_cell;
    use super::*;

    #[test]
    fn named_procedure_rejects_foreign_or_incomplete_start_before_host_effects() {
        let mut output = Vec::new();
        assert!(session(&mut &[0, 0, 0, 5, 1][..], &mut output).is_err());
        assert!(output.is_empty());
        let mut input = Vec::new();
        write_frame(
            &mut input,
            vec![row(["run", "BP7_criteria", "0", "1", "2", "foreign"])],
        )
        .unwrap();
        assert!(session(&mut input.as_slice(), &mut output).is_err());
        assert!(output.is_empty());
    }

    #[test]
    fn procedure_literal_wire_preserves_bits_and_never_accepts_noncanonical_handles() {
        assert!(handle("-1").is_err());
        assert!(handle("true").is_err());
        let cell = Cell::Float(f64::from_bits(0x8000000000000000));
        let (tag, raw) = encode_cell(&cell);
        let Cell::Float(decoded) = decode_cell(tag, &raw).unwrap() else {
            panic!()
        };
        assert_eq!(decoded.to_bits(), 0x8000000000000000);
    }
}
