//! Strict, ordinal TSV session transport for source-ordered effect frontiers.
//!
//! Text payloads use UTF-8 hex and f64 values use sixteen raw IEEE hex digits.
//! Each request has one response. A malformed request closes the session; the
//! caller must not replay a public function after any host effect has started.

use std::io::{BufRead, Read, Write};
use std::collections::BTreeMap;

use super::evidence::{self, Evidence, Side};
use super::interpreter::{self, Interpretation, RawEffect, Spec, Stage};
use super::{PvalueRoute, combine_mask, contains, high_confidence_mask, pvalue_route, strength_table};

pub const PROTOCOL: &str = "reported-criteria-v1";
const MAX_LINE: usize = 16 * 1024 * 1024;

fn text(value: &str) -> String {
    value.as_bytes().iter().map(|byte| format!("{byte:02x}")).collect()
}

fn untext(value: &str) -> Result<String, String> {
    if !value.len().is_multiple_of(2) || !value.bytes().all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase()) {
        return Err("noncanonical UTF-8 hex".into());
    }
    let bytes: Result<Vec<u8>, _> = (0..value.len()).step_by(2).map(|index| u8::from_str_radix(&value[index..index + 2], 16)).collect();
    String::from_utf8(bytes.map_err(|_| "invalid UTF-8 hex")?).map_err(|_| "invalid UTF-8".into())
}

fn bit(value: &str) -> Result<bool, String> {
    match value { "0" => Ok(false), "1" => Ok(true), _ => Err("invalid Boolean bit".into()) }
}

fn number(value: &str) -> Result<f64, String> {
    if value.len() != 16 || !value.bytes().all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase()) {
        return Err("noncanonical f64 bits".into());
    }
    Ok(f64::from_bits(u64::from_str_radix(value, 16).map_err(|_| "invalid f64 bits")?))
}

fn optional(value: &str) -> Result<Option<f64>, String> {
    if value == "-" { Ok(None) } else { number(value).map(Some) }
}

fn bits(value: Option<f64>) -> String {
    value.map_or_else(|| "-".into(), |value| format!("{:016x}", value.to_bits()))
}

fn integer(value: &str) -> Result<usize, String> {
    let result: usize = value.parse().map_err(|_| "invalid source ordinal")?;
    if result.to_string() != value { return Err("noncanonical source ordinal".into()); }
    Ok(result)
}

fn exact(fields: &[&str], size: usize) -> Result<(), String> {
    if fields.len() == size { Ok(()) } else { Err("invalid reported request width".into()) }
}

fn spec(fields: &[&str]) -> Result<Spec, String> {
    exact(fields, 9)?;
    Ok(Spec { family: fields[0].into(), subtype: fields[1].into(), directionality: fields[2].into(),
        scale: fields[3].into(), lo: number(fields[4])?, hi: number(fields[5])?,
        band_text: untext(fields[6])?, units: untext(fields[7])?, flexibility: bit(fields[8])? })
}

fn result(value: Interpretation) -> Vec<String> {
    vec!["return".into(), text(&value.call), bits(value.standardized), text(&value.reason), text(&value.notes), u8::from(value.make_higher).to_string(),
        if value.notes_are_units { "units" } else { "text" }.into()]
}

fn stage(value: Stage) -> Vec<String> {
    match value {
        Stage::Return(value) => result(value),
        Stage::Standardize => vec!["standardize".into()],
        Stage::BindingLog10 => vec!["binding-log10".into()],
        Stage::Raw(effect) => vec!["raw".into(), match effect {
            RawEffect::DgNegate => "dg-negate", RawEffect::DgIdentity => "dg-identity",
            RawEffect::BindingKd => "binding-kd", RawEffect::BindingIdentity => "binding-identity",
        }.into()],
    }
}

fn bools(fields: &[&str]) -> Result<Vec<bool>, String> { fields.iter().map(|value| bit(value)).collect() }
fn list(values: &[usize]) -> String { values.iter().map(usize::to_string).collect::<Vec<_>>().join(",") }
fn ordinals(value: &str) -> Result<Vec<usize>, String> {
    if value.is_empty() { Ok(Vec::new()) } else { value.split(',').map(integer).collect() }
}

pub fn encode_text(value: &str) -> String { text(value) }
pub fn decode_text(value: &str) -> Result<String, String> { untext(value) }
pub fn encode_number(value: Option<f64>) -> String { bits(value) }
pub fn decode_number(value: &str) -> Result<f64, String> { number(value) }
pub fn decode_optional(value: &str) -> Result<Option<f64>, String> { optional(value) }
pub fn ordinal_list(values: &[usize]) -> String { list(values) }

fn records(fields: &[&str]) -> Result<Vec<Evidence>, String> {
    if !fields.len().is_multiple_of(3) { return Err("invalid evidence record width".into()); }
    fields.chunks_exact(3).enumerate().map(|(ordinal, row)| Ok(Evidence {
        ordinal, family: integer(row[0])?, side: match row[1] { "PS3" => Side::Ps3, "BS3" => Side::Bs3, _ => return Err("invalid evidence side".into()) },
        formatted: untext(row[2])?,
    })).collect()
}

pub fn request(operation: &str, fields: &[&str]) -> Result<Vec<String>, String> {
    match operation {
        "hello" => { exact(fields, 1)?; if fields[0] != PROTOCOL { return Err("unsupported reported protocol".into()); } Ok(vec![PROTOCOL.into()]) }
        "begin" => { exact(fields, 10)?; Ok(stage(interpreter::begin(&spec(&fields[..9])?, bit(fields[9])?))) }
        "after-standardize" => { exact(fields, 11)?; Ok(stage(interpreter::after_standardize(&spec(&fields[..9])?, optional(fields[9])?, &untext(fields[10])?))) }
        "gof-direction-needed" => {
            exact(fields, 2)?;
            Ok(vec![u8::from(fields[0] == "Electrophysiology" && optional(fields[1])?.is_some()).to_string()])
        }
        "after-raw" => { exact(fields, 11)?; Ok(result(interpreter::after_raw(&spec(&fields[..9])?, number(fields[9])?, &untext(fields[10])?))) }
        "binding-float" => { exact(fields, 1)?; Ok(stage(interpreter::after_binding_float(number(fields[0])?))) }
        "validate" => { Ok(vec![u8::from(interpreter::valid_spec(&spec(fields)?)).to_string()]) }
        "band-call" => { exact(fields, 10)?; Ok(vec![interpreter::call_from_band(optional(fields[9])?, &spec(&fields[..9])?, false).into()]) }
        "standardize" => {
            exact(fields, 2)?;
            Ok(vec![match interpreter::standardize_effect(bit(fields[0])?, fields[1]) {
                interpreter::StandardizeEffect::None => "none", interpreter::StandardizeEffect::Float => "float", interpreter::StandardizeEffect::Log2 => "log2",
            }.into()])
        }
        "log2-needed" => { exact(fields, 1)?; Ok(vec![u8::from(interpreter::log2_needed(optional(fields[0])?)).to_string()]) }
        "explicit-band" | "percentile-band" => {
            let value = if operation == "explicit-band" { exact(fields, 3)?; interpreter::explicit_band(optional(fields[0])?, optional(fields[1])?, bit(fields[2])?) }
                else { exact(fields, 2)?; interpreter::percentile_band(number(fields[0])?, number(fields[1])?) };
            Ok(value.map_or_else(|| vec!["none".into()], |(lo, hi)| vec!["band".into(), bits(Some(lo)), bits(Some(hi))]))
        }
        "quality" => { exact(fields, 3)?; Ok(vec![match interpreter::quality_route(bit(fields[0])?, bit(fields[1])?, bit(fields[2])?) {
            interpreter::QualityRoute::Ambiguous => "ambiguous", interpreter::QualityRoute::GeneRules => "gene-rules", interpreter::QualityRoute::Band => "band",
        }.into()]) }
        "metadata-criteria" => { exact(fields, 1)?; let (ps3, bs3) = interpreter::metadata_criteria(fields[0]); Ok(vec![u8::from(ps3).to_string(), u8::from(bs3).to_string()]) }
        "metadata-enums" => {
            exact(fields, 4)?;
            let route = match interpreter::FAMILIES.iter().find(|(family, _)| *family == fields[0]) {
                None => "unknown-family",
                Some((_, subtypes)) if !subtypes.contains(&fields[1]) => "unknown-subtype",
                _ if !["HIGHER_IS_FUNCTION", "LOWER_IS_FUNCTION", "NON_MONOTONIC"].contains(&fields[2]) => "unknown-directionality",
                _ if !["CENTERED_0", "WT_EQ_1", "RAW", "PERCENT", "FRACTION"].contains(&fields[3]) => "unknown-scale",
                _ => "spec",
            };
            Ok(vec![route.into()])
        }
        "pvalue" => { exact(fields, 1)?; Ok(vec![match pvalue_route(number(fields[0])?) { PvalueRoute::Interpreter => "interpreter", PvalueRoute::Ps3 => "PS3", PvalueRoute::Bs3 => "BS3" }.into()]) }
        "contains" => {
            if fields.is_empty() { return Err("missing pattern".into()); }
            let values: Result<Vec<_>, _> = fields[1..].iter().map(|value| untext(value)).collect();
            Ok(contains(&values?, fields[0])?.iter().map(|&value| u8::from(value).to_string()).collect())
        }
        "urn-mask" => {
            if fields.is_empty() || !(fields.len() - 1).is_multiple_of(2) { return Err("invalid typed URN column".into()); }
            let urn = untext(fields[0])?;
            let mut hits = Vec::new();
            for row in fields[1..].chunks_exact(2) {
                hits.push(match row[0] {
                    "str" => untext(row[1])? == urn,
                    "missing" if row[1].is_empty() => false,
                    _ => return Err("unsupported typed URN value".into()),
                });
            }
            let first = hits.iter().position(|hit| *hit);
            let mut result = vec![u8::from(first.is_none()).to_string(), first.map_or_else(|| "-".into(), |value| value.to_string())];
            result.extend(hits.into_iter().map(|hit| u8::from(hit).to_string()));
            Ok(result)
        }
        "high-confidence" => {
            if !fields.len().is_multiple_of(2) { return Err("invalid rank pairs".into()); }
            let found = fields.chunks_exact(2).map(|row| bit(row[0])).collect::<Result<Vec<_>, _>>()?;
            let ranks = fields.chunks_exact(2).map(|row| optional(row[1])).collect::<Result<Vec<_>, _>>()?;
            Ok(high_confidence_mask(&found, &ranks)?.iter().map(|&value| u8::from(value).to_string()).collect())
        }
        "logical" => {
            if fields.is_empty() || !(fields.len() - 1).is_multiple_of(2) { return Err("invalid logical mask pairs".into()); }
            let left = fields[1..].chunks_exact(2).map(|row| bit(row[0])).collect::<Result<Vec<_>, _>>()?;
            let right = fields[1..].chunks_exact(2).map(|row| bit(row[1])).collect::<Result<Vec<_>, _>>()?;
            Ok(super::logical_mask(&left, &right, fields[0])?.iter().map(|&value| u8::from(value).to_string()).collect())
        }
        "combine" | "strengths" => {
            if fields.is_empty() || !(fields.len() - 1).is_multiple_of(2) { return Err("invalid mask pairs".into()); }
            let mode = bit(fields[0])?;
            let left = bools(&fields[1..].chunks_exact(2).map(|row| row[0]).collect::<Vec<_>>())?;
            let right = bools(&fields[1..].chunks_exact(2).map(|row| row[1]).collect::<Vec<_>>())?;
            if operation == "combine" { return Ok(combine_mask(&left, &right, mode)?.iter().map(|&value| u8::from(value).to_string()).collect()); }
            let table = strength_table(&left, &right, mode)?;
            let first = table.get_columns()[1].i64().map_err(|error| error.to_string())?;
            let second = table.get_columns()[2].i64().map_err(|error| error.to_string())?;
            let mut output = Vec::with_capacity(left.len() * 2);
            for index in 0..left.len() {
                output.push(first.get(index).ok_or("null native strength")?.to_string());
                output.push(second.get(index).ok_or("null native strength")?.to_string());
            }
            Ok(output)
        }
        "resolve" => {
            if fields.is_empty() { return Err("missing conflict order".into()); }
            let output = evidence::resolve(&records(&fields[1..])?, &ordinals(fields[0])?)?;
            let mut result = vec![list(&output.final_ps3), list(&output.final_bs3), list(&output.clean_ps3), list(&output.clean_bs3),
                list(&output.ps3_families), list(&output.bs3_families), output.conflicts.len().to_string()];
            for item in output.conflicts {
                result.extend([item.family.to_string(), item.ps3.to_string(), item.bs3.to_string(), match item.winner { Some(Side::Ps3) => "PS3", Some(Side::Bs3) => "BS3", None => "tie" }.into(),
                    list(&item.ps3_ordinals), list(&item.bs3_ordinals)]);
            }
            Ok(result)
        }
        "format-evidence" => {
            if fields.is_empty() { return Err("missing formatting order".into()); }
            Ok(vec![text(&evidence::format_evidence(&records(&fields[1..])?, &ordinals(fields[0])?)?)])
        }
        "render-evidence" => {
            if !fields.len().is_multiple_of(5) { return Err("invalid evidence formatting width".into()); }
            let mut values = Vec::new();
            for row in fields.chunks_exact(5) {
                let row = row.iter().map(|value| untext(value)).collect::<Result<Vec<_>, _>>()?;
                values.push(format!("{}({}/{}, score={}, call={})", row[0], row[1], row[2], row[3], row[4]));
            }
            Ok(vec![text(&values.join("; "))])
        }
        "render-score-info" => {
            exact(fields, 3)?;
            if !["PS3", "BS3"].contains(&fields[0]) { return Err("invalid score INFO side".into()); }
            Ok(vec![text(&format!("There are {} variants with MaveDB determined {} criteria, they are determined by these functional assays: \n{}", untext(fields[1])?, fields[0], untext(fields[2])?))])
        }
        "render-metadata-info" => {
            exact(fields, 2)?;
            Ok(vec![text(&format!("There are {} unique URNs in the MaveDB metadata, which looks like: \n{}", untext(fields[0])?, untext(fields[1])?))])
        }
        "render-missing-warning" => {
            exact(fields, 1)?;
            Ok(vec![text(&format!("The URN {} is not in the MaveDB metadata, it is ignored", untext(fields[0])?))])
        }
        "render-pvalue-evidence" => {
            exact(fields, 2)?;
            let call = match fields[0] { "PS3" => "pvalue<0.05", "BS3" => "pvalue>=0.05", _ => return Err("invalid pvalue evidence side".into()) };
            Ok(vec![text(call), text(&format!("pvalue={}", untext(fields[1])?))])
        }
        "join-interpretations" => {
            let values = fields.iter().map(|value| untext(value)).collect::<Result<Vec<_>, _>>()?;
            Ok(vec![text(&values.join("&"))])
        }
        _ => Err("unknown reported procedure".into()),
    }
}

pub fn serve<R: BufRead, W: Write>(mut input: R, mut output: W) -> Result<(), String> {
    let mut ordinal = 0_usize;
    let mut next_context = 0_usize;
    let mut programs = BTreeMap::<usize, super::Program>::new();
    loop {
        let mut line = String::new();
        let bytes = Read::by_ref(&mut input).take(MAX_LINE as u64 + 1).read_line(&mut line).map_err(|error| error.to_string())?;
        if bytes == 0 { return Ok(()); }
        if bytes > MAX_LINE || !line.ends_with('\n') || !line.is_ascii() || line.contains('\r') {
            return Err("invalid reported request framing".into());
        }
        let fields: Vec<_> = line[..line.len() - 1].split('\t').collect();
        if fields.len() < 2 || integer(fields[0])? != ordinal || ((ordinal == 0) != (fields[1] == "hello")) {
            return Err("invalid reported request ordinal/handshake".into());
        }
        let response = if fields[1] == "start" {
            if fields.len() != 3 || programs.len() >= 32 { Err("invalid reported program start".into()) }
            else {
                let mut program = super::Program::new(fields[2])?;
                let effect = program.advance(&[])?;
                let context = next_context;
                next_context = next_context.checked_add(1).ok_or("reported context overflow")?;
                programs.insert(context, program);
                let mut response = vec![context.to_string()]; response.extend(effect.0); Ok(response)
            }
        } else if fields[1] == "next" {
            if fields.len() < 3 { Err("missing reported context".into()) }
            else {
                let context = integer(fields[2])?;
                let program = programs.get_mut(&context).ok_or("unknown or completed reported context")?;
                let values = fields[3..].iter().map(|value| (*value).to_owned()).collect::<Vec<_>>();
                let effect = program.resume(&values)?;
                if effect.0.first().is_some_and(|value| value == "return") { programs.remove(&context); }
                let mut response = vec![context.to_string()]; response.extend(effect.0); Ok(response)
            }
        } else { request(fields[1], &fields[2..]) };
        let (status, reply) = match response {
            Ok(reply) => ("ok", reply), Err(error) => ("error", vec![text(&error)]),
        };
        let mut row = vec![ordinal.to_string(), status.into()]; row.extend(reply);
        output.write_all(row.join("\t").as_bytes()).map_err(|error| error.to_string())?;
        output.write_all(b"\n").map_err(|error| error.to_string())?;
        output.flush().map_err(|error| error.to_string())?;
        if status == "error" { return Err("reported session rejected request".into()); }
        ordinal = ordinal.checked_add(1).ok_or("reported request ordinal overflow")?;
    }
}

/// The shared CLI can map `reported-criteria serve` here without changing any
/// scientific procedure or accepting a different transport profile.
pub fn dispatch(args: &[String]) -> Result<(), String> {
    if args == ["serve"] { return serve(std::io::stdin().lock(), std::io::stdout().lock()); }
    Err("expected reported-criteria serve".into())
}
