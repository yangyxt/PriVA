//! One persistent, framed process controls all current-record filter decisions.
use super::family;
use super::genotype::{Allele, Genotype};
use super::variant_source;
use super::{SourceError, SourceResult};
use crate::cli::{Result, invalid};
use std::collections::BTreeMap;
use std::io::{self, Read, Write};

// 2026-10-27: explicit v2 is separate from the preserved stateful v1 control.
#[path = "transition.rs"]
mod transition;

const VERSION: &str = "priva-pedigree-v1";
const MAX_FRAME: usize = 64 * 1024 * 1024;
type Fields = BTreeMap<String, String>;

fn take(fields: &mut Fields, name: &str) -> Result<String> {
    fields
        .remove(name)
        .ok_or_else(|| invalid(format!("missing protocol field {name}")))
}

fn number(fields: &mut Fields, name: &str) -> Result<usize> {
    let value = take(fields, name)?;
    let parsed = value
        .parse::<usize>()
        .map_err(|_| invalid(format!("invalid integer {name}")))?;
    if value != parsed.to_string() {
        return Err(invalid(format!("noncanonical integer {name}")));
    }
    Ok(parsed)
}

fn finished(fields: Fields) -> Result<()> {
    if fields.is_empty() {
        Ok(())
    } else {
        Err(invalid("unknown or extra protocol fields"))
    }
}

fn decimal(value: &str) -> bool {
    if value == "0" {
        return true;
    }
    let body = value.strip_prefix('-').unwrap_or(value);
    !body.is_empty() && !body.starts_with('0') && body.bytes().all(|value| value.is_ascii_digit())
}

fn genotype(fields: &mut Fields, prefix: &str) -> Result<Genotype> {
    let kind = take(fields, &format!("{prefix}.kind"))?;
    match kind.as_str() {
        "none" => Ok(Genotype::None),
        "int" | "bool" => {
            let value = take(fields, &format!("{prefix}.value"))?;
            if (kind == "int" && !decimal(&value))
                || (kind == "bool" && !["0", "1"].contains(&value.as_str()))
            {
                return Err(invalid("noncanonical scalar genotype"));
            }
            Ok(Genotype::Scalar(if kind == "int" { "int" } else { "bool" }))
        }
        "str" => Ok(Genotype::String(take(fields, &format!("{prefix}.value"))?)),
        "tuple" | "list" => {
            let count = number(fields, &format!("{prefix}.count"))?;
            let mut values = Vec::new();
            for index in 0..count {
                let stem = format!("{prefix}.{index}");
                let kind = take(fields, &format!("{stem}.kind"))?;
                values.push(match kind.as_str() {
                    "none" => Allele::None,
                    "str" => Allele::String(take(fields, &format!("{stem}.value"))?),
                    "int" => {
                        let value = take(fields, &format!("{stem}.value"))?;
                        if !decimal(&value) {
                            return Err(invalid("noncanonical integer allele"));
                        }
                        Allele::Integer(value)
                    }
                    "bool" => match take(fields, &format!("{stem}.value"))?.as_str() {
                        "0" => Allele::Bool(false),
                        "1" => Allele::Bool(true),
                        _ => return Err(invalid("invalid Boolean allele")),
                    },
                    _ => return Err(invalid("unqualified allele kind")),
                });
            }
            Ok(Genotype::Sequence(values))
        }
        _ => Err(invalid("unqualified genotype kind")),
    }
}

fn vector(fields: &mut Fields, prefix: &str) -> Result<Vec<String>> {
    let count = number(fields, &format!("{prefix}.count"))?;
    (0..count)
        .map(|index| take(fields, &format!("{prefix}.{index}")))
        .collect()
}

fn add_vector(output: &mut Fields, prefix: &str, values: &[String]) {
    output.insert(format!("{prefix}.count"), values.len().to_string());
    for (index, value) in values.iter().enumerate() {
        output.insert(format!("{prefix}.{index}"), value.clone());
    }
}

fn outcome<T: ToString>(value: SourceResult<T>) -> Fields {
    match value {
        Ok(value) => BTreeMap::from([
            ("status".into(), "ok".into()),
            ("value".into(), value.to_string()),
        ]),
        Err(SourceError { class, message }) => BTreeMap::from([
            ("status".into(), "source-error".into()),
            ("class".into(), class.into()),
            ("message".into(), message),
        ]),
    }
}

fn action(name: &str) -> Fields {
    BTreeMap::from([
        ("status".into(), "ok".into()),
        ("action".into(), name.into()),
    ])
}

#[derive(Clone, Copy)]
enum Phase {
    Idle,
    Patient(usize),
    Control(usize),
    Chrom,
    Proband,
    Father,
    Mother,
    DebugAuto,
    DebugDrop,
    Done,
}

struct Filter {
    patients: Vec<String>,
    controls: Vec<String>,
    proband: String,
    father: String,
    mother: String,
    sex: String,
    phase: Phase,
    ordinal: usize,
    chrom_type: String,
    proband_gt: Genotype,
    father_gt: Genotype,
    mother_gt: Genotype,
}

impl Filter {
    fn new(fields: &mut Fields) -> Result<Self> {
        Ok(Self {
            patients: vector(fields, "patients")?,
            controls: vector(fields, "controls")?,
            proband: take(fields, "proband")?,
            father: take(fields, "father")?,
            mother: take(fields, "mother")?,
            sex: take(fields, "sex")?,
            phase: Phase::Idle,
            ordinal: 0,
            chrom_type: String::new(),
            proband_gt: Genotype::None,
            father_gt: Genotype::None,
            mother_gt: Genotype::None,
        })
    }

    fn patient(&mut self, index: usize) -> Fields {
        if let Some(sample) = self.patients.get(index) {
            self.phase = Phase::Patient(index);
            let mut result = action("patient-gt");
            result.insert("sample".into(), sample.clone());
            result
        } else {
            self.phase = Phase::DebugDrop;
            action("debug-no-patient")
        }
    }

    fn control(&mut self, index: usize) -> Fields {
        if let Some(sample) = self.controls.get(index) {
            self.phase = Phase::Control(index);
            let mut result = action("control-gt");
            result.insert("sample".into(), sample.clone());
            result
        } else {
            self.phase = Phase::Chrom;
            action("constructor-chrom")
        }
    }

    fn source_gt(&mut self, phase: Phase, sample: &str) -> Fields {
        self.phase = phase;
        let mut result = action("source-gt");
        result.insert("sample".into(), sample.into());
        result
    }

    fn done(&mut self, name: &str, source: Option<&str>) -> Fields {
        self.phase = Phase::Done;
        let mut result = action(name);
        if let Some(value) = source {
            result.insert("value".into(), value.into());
        }
        result
    }

    fn finish_source(&mut self) -> SourceResult<Fields> {
        let value = variant_source::source(
            &self.proband_gt,
            &self.chrom_type,
            &self.sex,
            &self.father_gt,
            &self.mother_gt,
        )?;
        Ok(self.done("write", Some(value)))
    }

    fn reply(&mut self, fields: &mut Fields) -> Result<Fields> {
        let mut source_error = None;
        let result = match self.phase {
            Phase::Patient(index) | Phase::Control(index) => {
                let patient = matches!(self.phase, Phase::Patient(_));
                let present = take(fields, "present")?;
                if !["0", "1"].contains(&present.as_str()) {
                    return Err(invalid("invalid membership flag"));
                }
                let kept = if present == "1" {
                    let gt = genotype(fields, "gt")?;
                    if patient {
                        gt.carrier()
                    } else {
                        gt.homozygous_or_hemizygous()
                    }
                } else {
                    Ok(false)
                };
                match kept {
                    Err(error) => {
                        source_error = Some(error);
                        action("unreachable")
                    }
                    Ok(true) if patient => self.control(0),
                    Ok(true) => {
                        self.phase = Phase::DebugDrop;
                        let mut value = action("debug-control-reject");
                        value.insert("sample".into(), self.controls[index].clone());
                        value
                    }
                    Ok(false) if patient => self.patient(index + 1),
                    Ok(false) => self.control(index + 1),
                }
            }
            Phase::Chrom => {
                if take(fields, "ack")? != "1" || self.chrom_type.is_empty() {
                    return Err(invalid("exact native constructor classification required"));
                }
                self.source_gt(Phase::Proband, &self.proband.clone())
            }
            Phase::Proband => {
                self.proband_gt = genotype(fields, "gt")?;
                match self.proband_gt.has_variant() {
                    Ok(false) => self.done("write", Some("UNKNOWN")),
                    Ok(true) => self.source_gt(Phase::Father, &self.father.clone()),
                    Err(error) => {
                        source_error = Some(error);
                        action("unreachable")
                    }
                }
            }
            Phase::Father => {
                self.father_gt = genotype(fields, "gt")?;
                self.source_gt(Phase::Mother, &self.mother.clone())
            }
            Phase::Mother => {
                self.mother_gt = genotype(fields, "gt")?;
                if self.chrom_type == "AUTOSOMAL" || (self.chrom_type == "X" && self.sex == "2") {
                    self.phase = Phase::DebugAuto;
                    action("debug-autosomal")
                } else {
                    match self.finish_source() {
                        Ok(value) => value,
                        Err(error) => {
                            source_error = Some(error);
                            action("unreachable")
                        }
                    }
                }
            }
            Phase::DebugAuto => {
                if take(fields, "ack")? != "1" {
                    return Err(invalid("debug acknowledgment required"));
                }
                match self.finish_source() {
                    Ok(value) => value,
                    Err(error) => {
                        source_error = Some(error);
                        action("unreachable")
                    }
                }
            }
            Phase::DebugDrop => {
                if take(fields, "ack")? != "1" {
                    return Err(invalid("debug acknowledgment required"));
                }
                self.done("drop", None)
            }
            _ => return Err(invalid("filter reply outside current pending action")),
        };
        Ok(source_error.map_or(result, |error| outcome::<String>(Err(error))))
    }
}

fn dispatch_request(mut fields: Fields, filter: &mut Option<Filter>) -> Result<(Fields, bool)> {
    let op = take(&mut fields, "op")?;
    let sequence = number(&mut fields, "sequence")?;
    let mut stop = false;
    let mut output = match op.as_str() {
        "hello" => {
            if take(&mut fields, "protocol")? != VERSION {
                return Err(invalid("unqualified pedigree protocol"));
            }
            outcome(Ok(VERSION))
        }
        "quit" => {
            stop = true;
            outcome(Ok("closed"))
        }
        "genotype" => {
            let predicate = take(&mut fields, "predicate")?;
            let gt = genotype(&mut fields, "gt")?;
            outcome(match predicate.as_str() {
                "homo" => gt.homozygous_or_hemizygous(),
                "carrier" => gt.carrier(),
                "exists" => gt.exists(),
                "has" => gt.has_variant(),
                _ => return Err(invalid("unknown genotype predicate")),
            })
        }
        "chromosome" => {
            let value = variant_source::chromosome_type(&take(&mut fields, "chrom")?);
            if let Some(state) = filter.as_mut() {
                if matches!(state.phase, Phase::Chrom) {
                    if !state.chrom_type.is_empty() {
                        return Err(invalid("duplicate constructor classification"));
                    }
                    // Reuse this constructor's actual native classification;
                    // do not evaluate the original string operation twice.
                    state.chrom_type = value.into();
                }
            }
            outcome(Ok(value))
        }
        "route" => {
            let chrom = take(&mut fields, "chrom_type")?;
            let sex = take(&mut fields, "sex")?;
            outcome(Ok(match chrom.as_str() {
                "AUTOSOMAL" => "auto",
                "X" if sex == "2" => "auto",
                "X" => "X",
                "Y" => "Y",
                _ => "unknown",
            }))
        }
        "source" => {
            let method = take(&mut fields, "method")?;
            let sex = take(&mut fields, "sex")?;
            let proband = genotype(&mut fields, "proband")?;
            let father = genotype(&mut fields, "father")?;
            let mother = genotype(&mut fields, "mother")?;
            outcome(match method.as_str() {
                "auto" => variant_source::autosomal(&proband, &father, &mother),
                "X" => variant_source::chr_x(&proband, &sex, &father, &mother),
                "Y" => variant_source::chr_y(&sex, &father),
                _ => return Err(invalid("unknown direct source method")),
            })
        }
        "family" => {
            let count = number(&mut fields, "rows")?;
            let target_kind = take(&mut fields, "target.kind")?;
            if !["str", "other"].contains(&target_kind.as_str()) {
                return Err(invalid("unqualified target comparison kind"));
            }
            let target = take(&mut fields, "target.value")?;
            let mut rows = Vec::new();
            for row in 0..count {
                rows.push(
                    (0..6)
                        .map(|column| take(&mut fields, &format!("row.{row}.{column}")))
                        .collect::<Result<Vec<_>>>()?,
                );
            }
            let plan = family::select(&rows, (target_kind == "str").then_some(target.as_str()))?;
            let mut value = action("family");
            add_vector(
                &mut value,
                "positions",
                &plan
                    .positions
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>(),
            );
            add_vector(&mut value, "patients", &plan.patients);
            add_vector(&mut value, "controls", &plan.controls);
            value.insert("proband".into(), plan.proband.unwrap_or_default());
            value.insert("sex".into(), plan.sex.unwrap_or_default());
            value
        }
        "parents" => {
            let mut value = action("parents");
            for role in ["father", "mother"] {
                let count = number(&mut fields, &format!("{role}.count"))?;
                let mut rows = Vec::new();
                for index in 0..count {
                    let prefix = format!("{role}.{index}");
                    rows.push(
                        match take(&mut fields, &format!("{prefix}.kind"))?.as_str() {
                            "none" => None,
                            "str" => Some(take(&mut fields, &format!("{prefix}.value"))?),
                            _ => return Err(invalid("unqualified normalized parent cell")),
                        },
                    );
                }
                let selected = family::first_parent(&rows);
                value.insert(
                    format!("{role}.present"),
                    if selected.is_some() { "1" } else { "0" }.into(),
                );
                if let Some((position, text)) = selected {
                    value.insert(format!("{role}.position"), position.to_string());
                    value.insert(format!("{role}.value"), text);
                }
            }
            value
        }
        "filter-init" => {
            if filter.is_some() {
                return Err(invalid("filter context already initialized"));
            }
            *filter = Some(Filter::new(&mut fields)?);
            outcome(Ok("ready"))
        }
        "filter-begin" => {
            let state = filter
                .as_mut()
                .ok_or_else(|| invalid("filter not initialized"))?;
            if !matches!(state.phase, Phase::Idle | Phase::Done) {
                return Err(invalid("previous record action is incomplete"));
            }
            let ordinal = number(&mut fields, "ordinal")?;
            if ordinal != state.ordinal + 1 {
                return Err(invalid("record source ordinal changed"));
            }
            state.ordinal = ordinal;
            state.chrom_type.clear();
            state.patient(0)
        }
        "filter-reply" => filter
            .as_mut()
            .ok_or_else(|| invalid("filter not initialized"))?
            .reply(&mut fields)?,
        _ => return Err(invalid("unknown pedigree operation")),
    };
    finished(fields)?;
    output.insert("sequence".into(), sequence.to_string());
    Ok((output, stop))
}

fn read_frame(input: &mut impl Read) -> Result<Option<Fields>> {
    let mut length = [0_u8; 8];
    if input.read(&mut length[..1])? == 0 {
        return Ok(None);
    }
    input.read_exact(&mut length[1..])?;
    let length = u64::from_be_bytes(length);
    if length == 0 || length > MAX_FRAME as u64 {
        return Err(invalid("invalid pedigree frame length"));
    }
    let mut bytes = vec![0; length as usize];
    input.read_exact(&mut bytes)?;
    let mut reader = csv::ReaderBuilder::new().from_reader(bytes.as_slice());
    if reader.headers()?.iter().collect::<Vec<_>>() != ["name", "value"] {
        return Err(invalid("invalid pedigree frame columns"));
    }
    let mut fields = Fields::new();
    for row in reader.records() {
        let row = row?;
        if row.len() != 2
            || row[0].is_empty()
            || fields.insert(row[0].into(), row[1].into()).is_some()
        {
            return Err(invalid("empty or duplicate pedigree frame field"));
        }
    }
    Ok(Some(fields))
}

fn write_frame(output: &mut impl Write, fields: Fields) -> Result<()> {
    let mut writer = csv::Writer::from_writer(Vec::new());
    writer.write_record(["name", "value"])?;
    for (name, value) in fields {
        writer.write_record([name, value])?;
    }
    let bytes = writer
        .into_inner()
        .map_err(|error| invalid(error.to_string()))?;
    if bytes.is_empty() || bytes.len() > MAX_FRAME {
        return Err(invalid("invalid response length"));
    }
    output.write_all(&(bytes.len() as u64).to_be_bytes())?;
    output.write_all(&bytes)?;
    output.flush()?;
    Ok(())
}

pub fn dispatch(args: &[String]) -> Result<()> {
    if args == ["session", "--protocol", transition::VERSION] {
        return transition::dispatch();
    }
    if args == ["--help"] {
        println!(
            "priva pedigree-filtration session --protocol {VERSION}\nOne persistent process; genotype/chromosome/source/family/parents and whole-filter procedures use exact framed requests. Host I/O remains pandas/pysam. Explicit opt-in only."
        );
        return Ok(());
    }
    if args != ["session", "--protocol", VERSION] {
        return Err(invalid(
            "usage: priva pedigree-filtration session --protocol priva-pedigree-v1",
        ));
    }
    serve(&mut io::stdin().lock(), &mut io::stdout().lock())
}

fn serve(input: &mut impl Read, output: &mut impl Write) -> Result<()> {
    let mut filter = None;
    let mut expected_sequence = 1;
    while let Some(fields) = read_frame(input)? {
        if fields.get("sequence") != Some(&expected_sequence.to_string()) {
            return Err(invalid("pedigree request sequence changed"));
        }
        if expected_sequence == 1 && fields.get("op").map(String::as_str) != Some("hello") {
            return Err(invalid("pedigree hello required"));
        }
        if expected_sequence != 1 && fields.get("op").map(String::as_str) == Some("hello") {
            return Err(invalid("pedigree protocol cannot be renegotiated"));
        }
        let (response, stop) = dispatch_request(fields, &mut filter)?;
        write_frame(output, response)?;
        if stop {
            return Ok(());
        }
        expected_sequence += 1;
    }
    Err(invalid("pedigree session ended before explicit close"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn fields(values: &[(&str, &str)]) -> Fields {
        values
            .iter()
            .map(|(key, value)| ((*key).into(), (*value).into()))
            .collect()
    }

    fn hello() -> Fields {
        fields(&[("op", "hello"), ("sequence", "1"), ("protocol", VERSION)])
    }

    fn encoded(values: &[Fields]) -> Vec<u8> {
        let mut bytes = Vec::new();
        for value in values {
            write_frame(&mut bytes, value.clone()).unwrap();
        }
        bytes
    }

    fn raw_frame(body: &[u8]) -> Vec<u8> {
        let mut value = (body.len() as u64).to_be_bytes().to_vec();
        value.extend_from_slice(body);
        value
    }

    #[test]
    fn framed_session_requires_close_and_preserves_source_errors() {
        let requests = [
            hello(),
            fields(&[
                ("op", "genotype"),
                ("sequence", "2"),
                ("predicate", "carrier"),
                ("gt.kind", "str"),
                ("gt.value", "0/1"),
            ]),
            fields(&[("op", "quit"), ("sequence", "3")]),
        ];
        let mut output = Vec::new();
        serve(&mut Cursor::new(encoded(&requests)), &mut output).unwrap();
        let mut output = Cursor::new(output);
        assert_eq!(read_frame(&mut output).unwrap().unwrap()["value"], VERSION);
        let error = read_frame(&mut output).unwrap().unwrap();
        assert_eq!(error["status"], "source-error");
        assert_eq!(error["class"], "TypeError");
        assert_eq!(error["sequence"], "2");
        assert_eq!(read_frame(&mut output).unwrap().unwrap()["value"], "closed");
        assert!(read_frame(&mut output).unwrap().is_none());
        assert!(serve(&mut Cursor::new(encoded(&[hello()])), &mut Vec::new()).is_err());
    }

    #[test]
    fn actual_reader_rejects_duplicate_missing_utf8_and_truncated_fields() {
        for body in [
            b"name,value\nop,hello\nop,quit\n".as_slice(),
            b"name,value\n,empty\n",
            b"other,value\nop,hello\n",
            b"name,value\nop,\xff\n",
        ] {
            assert!(read_frame(&mut Cursor::new(raw_frame(body))).is_err());
        }
        let valid = encoded(&[hello()]);
        for length in [1, 7, valid.len() - 1] {
            assert!(read_frame(&mut Cursor::new(valid[..length].to_vec())).is_err());
        }
        for length in [0, MAX_FRAME as u64 + 1] {
            assert!(read_frame(&mut Cursor::new(length.to_be_bytes().to_vec())).is_err());
        }
    }

    #[test]
    fn session_rejects_wrong_order_renegotiation_and_extra_fields() {
        let mut bad_hello = hello();
        bad_hello.insert("foreign".into(), "extra".into());
        for requests in [
            vec![fields(&[("op", "quit"), ("sequence", "1")])],
            vec![bad_hello],
            vec![hello(), fields(&[("op", "quit"), ("sequence", "2.0")])],
            vec![hello(), fields(&[("op", "quit"), ("sequence", "3")])],
            vec![
                hello(),
                fields(&[("op", "hello"), ("sequence", "2"), ("protocol", VERSION)]),
            ],
        ] {
            assert!(serve(&mut Cursor::new(encoded(&requests)), &mut Vec::new()).is_err());
        }
        let malformed = fields(&[
            ("op", "genotype"),
            ("sequence", "2"),
            ("predicate", "has"),
            ("gt.kind", "tuple"),
            ("gt.count", "01"),
        ]);
        assert!(
            serve(
                &mut Cursor::new(encoded(&[hello(), malformed])),
                &mut Vec::new()
            )
            .is_err()
        );
    }

    fn filter() -> Filter {
        Filter::new(&mut fields(&[
            ("patients.count", "1"),
            ("patients.0", "p"),
            ("controls.count", "0"),
            ("proband", "p"),
            ("father", ""),
            ("mother", ""),
            ("sex", "1"),
        ]))
        .unwrap()
    }

    #[test]
    fn filter_transaction_guards_actual_lazy_phase_and_source_ordinal() {
        let mut state = Some(filter());
        let begin = fields(&[("op", "filter-begin"), ("sequence", "1"), ("ordinal", "1")]);
        let (reply, _) = dispatch_request(begin.clone(), &mut state).unwrap();
        assert_eq!(reply["action"], "patient-gt");
        assert!(dispatch_request(begin, &mut state).is_err());
        let (reply, _) = dispatch_request(
            fields(&[
                ("op", "filter-reply"),
                ("sequence", "2"),
                ("present", "1"),
                ("gt.kind", "tuple"),
                ("gt.count", "1"),
                ("gt.0.kind", "int"),
                ("gt.0.value", "1"),
            ]),
            &mut state,
        )
        .unwrap();
        assert_eq!(reply["action"], "constructor-chrom");
        assert!(
            state
                .as_mut()
                .unwrap()
                .reply(&mut fields(&[("ack", "1")]))
                .is_err()
        );
        dispatch_request(
            fields(&[("op", "chromosome"), ("sequence", "3"), ("chrom", "chrX")]),
            &mut state,
        )
        .unwrap();
        assert!(
            dispatch_request(
                fields(&[("op", "chromosome"), ("sequence", "4"), ("chrom", "chrY")]),
                &mut state
            )
            .is_err()
        );
        let (reply, _) = dispatch_request(
            fields(&[("op", "filter-reply"), ("sequence", "5"), ("ack", "1")]),
            &mut state,
        )
        .unwrap();
        assert_eq!(reply["action"], "source-gt");
        let (reply, _) = dispatch_request(
            fields(&[
                ("op", "filter-reply"),
                ("sequence", "6"),
                ("gt.kind", "none"),
            ]),
            &mut state,
        )
        .unwrap();
        assert_eq!(reply["value"], "UNKNOWN");
        assert!(
            dispatch_request(
                fields(&[("op", "filter-begin"), ("sequence", "7"), ("ordinal", "1")]),
                &mut state
            )
            .is_err()
        );
        assert!(
            dispatch_request(
                fields(&[("op", "filter-reply"), ("sequence", "8"), ("ack", "1")]),
                &mut state
            )
            .is_err()
        );
    }

    #[test]
    fn failed_writer_never_reports_a_completed_response() {
        struct Fail;
        impl Write for Fail {
            fn write(&mut self, _: &[u8]) -> io::Result<usize> {
                Err(io::Error::new(io::ErrorKind::BrokenPipe, "fixture"))
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        assert!(serve(&mut Cursor::new(encoded(&[hello()])), &mut Fail).is_err());
    }
}
