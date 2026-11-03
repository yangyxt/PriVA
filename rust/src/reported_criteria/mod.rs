//! Source-ordered PS3/BS3 and PP5/BP6 decisions.
//!
//! The Python facade retains pandas parsing/boxing, CPython formatting/set
//! iteration and `math.log` effects. These procedures own the PriVA rules and
//! return explicit next effects or ordered publications. No procedure performs
//! a file read, logs, or changes a caller-owned object.

pub mod evidence;
pub mod interpreter;
pub mod transport;

use polars_core::prelude::*;

/// Array writes are separate source operations; their order is part of the API.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Publication {
    MavedbPs3,
    MavedbBs3,
    MavedbInterpretation,
    Ps3,
    Bs3,
    Pp5,
    Bp6,
}

pub const ROW_PUBLICATIONS: [Publication; 3] = [
    Publication::MavedbPs3,
    Publication::MavedbBs3,
    Publication::MavedbInterpretation,
];
pub const PS3_PUBLICATIONS: [Publication; 2] = [Publication::Ps3, Publication::Bs3];
pub const PP5_PUBLICATIONS: [Publication; 2] = [Publication::Pp5, Publication::Bp6];

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PvalueRoute {
    Interpreter,
    Ps3,
    Bs3,
}

/// The original converts each token with Python float before this procedure.
/// In particular, exactly .05 belongs to BS3 and NaN requests interpretation.
pub fn pvalue_route(value: f64) -> PvalueRoute {
    if value.is_nan() {
        PvalueRoute::Interpreter
    } else if value < 0.05 {
        PvalueRoute::Ps3
    } else {
        PvalueRoute::Bs3
    }
}

/// Missing values have already passed the source's fillna/string seam. The
/// caller requests each pattern at the original getter position, not eagerly.
pub fn contains(values: &[String], pattern: &str) -> Result<Vec<bool>, String> {
    if !matches!(pattern, "Pathogenic" | "Benign" | "athogenic" | "enign") {
        return Err("unsupported reported-evidence pattern".into());
    }
    Ok(values.iter().map(|value| value.contains(pattern)).collect())
}

pub fn high_confidence_mask(found: &[bool], ranks: &[Option<f64>]) -> Result<Vec<bool>, String> {
    if found.len() != ranks.len() {
        return Err("confidence-mask length mismatch".into());
    }
    Ok(found.iter().zip(ranks).map(|(&hit, rank)| {
        hit && rank.is_some_and(|value| value >= 2.0)
    }).collect())
}

pub fn combine_mask(left: &[bool], right: &[bool], exclude: bool) -> Result<Vec<bool>, String> {
    if left.len() != right.len() {
        return Err("reported-mask length mismatch".into());
    }
    Ok(left.iter().zip(right).map(|(&a, &b)| if exclude { a && !b } else { a || b }).collect())
}

pub fn logical_mask(left: &[bool], right: &[bool], operation: &str) -> Result<Vec<bool>, String> {
    if left.len() != right.len() { return Err("logical-mask length mismatch".into()); }
    match operation {
        "and" => Ok(left.iter().zip(right).map(|(&left, &right)| left && right).collect()),
        "or" => Ok(left.iter().zip(right).map(|(&left, &right)| left || right).collect()),
        _ => Err("unsupported source mask operation".into()),
    }
}

/// A Polars result has stable source ordinals. Host NumPy assignments remain
/// separate so their genuine indexing/dtype failures retain the source prefix.
pub fn strength_table(first: &[bool], second: &[bool], reported: bool) -> Result<DataFrame, String> {
    if first.len() != second.len() {
        return Err("reported-strength length mismatch".into());
    }
    let strength = if reported { 1_i64 } else { 3_i64 };
    let ordinal: Vec<u64> = (0..first.len()).map(|value| value as u64).collect();
    let left: Vec<i64> = first.iter().map(|&hit| if hit { strength } else { 0 }).collect();
    let right: Vec<i64> = second.iter().map(|&hit| if hit { strength } else { 0 }).collect();
    DataFrame::new(vec![
        Series::new("source_position".into(), ordinal).into(),
        Series::new((if reported { "PP5" } else { "PS3" }).into(), left).into(),
        Series::new((if reported { "BP6" } else { "BS3" }).into(), right).into(),
    ]).map_err(|error| error.to_string())
}

/// A host effect is a single mature operation or a single source publication.
/// The host never chooses a PriVA branch or the next effect.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Effect(pub Vec<String>);

impl Effect {
    pub fn host(name: &str, arguments: &[&str]) -> Self {
        Self(std::iter::once("effect").chain(std::iter::once(name)).chain(arguments.iter().copied()).map(str::to_owned).collect())
    }
    pub fn owned(name: &str, arguments: Vec<String>) -> Self {
        let mut fields = vec!["effect".into(), name.into()]; fields.extend(arguments); Self(fields)
    }
    pub fn finish(names: &[&str]) -> Self {
        Self(std::iter::once("return").chain(names.iter().copied()).map(str::to_owned).collect())
    }
}

pub fn reply<'a>(fields: &'a [String], expected: &str, width: usize) -> Result<&'a [String], String> {
    if fields.len() == width + 1 && fields[0] == expected { Ok(&fields[1..]) }
    else { Err(format!("invalid {expected} effect reply")) }
}

pub fn truth(fields: &[String]) -> Result<bool, String> {
    match reply(fields, "bool", 1)?[0].as_str() { "0" => Ok(false), "1" => Ok(true), _ => Err("invalid effect truth bit".into()) }
}

fn done(fields: &[String]) -> Result<(), String> { reply(fields, "done", 0).map(|_| ()) }

/// Complete public-unit control program. All getters, dataframe effects, logs
/// and writes are selected here, including effects that later raise in pandas.
/// Their source expressions are the host's finite, source-pinned effect map.
pub struct Program {
    kind: String,
    pc: u16,
    row: Option<RowProgram>,
    metadata: Option<interpreter::MetadataProgram>,
}

impl Program {
    pub fn new(kind: &str) -> Result<Self, String> {
        if !["mavedb_interpretation_per_row", "mavedb_score_interpretation", "PS3_BS3_criteria", "PP5_BP6_criteria",
            "parse_neutral_band_from_metadata", "pick_interpreter_from_metadata"].contains(&kind) {
            return Err("unknown reported control program".into());
        }
        Ok(Self { kind: kind.into(), pc: 0,
            row: (kind == "mavedb_interpretation_per_row").then(RowProgram::default),
            metadata: if kind == "parse_neutral_band_from_metadata" || kind == "pick_interpreter_from_metadata" {
                Some(interpreter::MetadataProgram::new(kind))
            } else { None },
        })
    }

    pub fn advance(&mut self, fields: &[String]) -> Result<Effect, String> {
        if let Some(row) = &mut self.row { return row.advance(fields); }
        if let Some(metadata) = &mut self.metadata { return metadata.advance(fields); }
        let pc = self.pc;
        self.pc = self.pc.checked_add(1).ok_or("reported program counter overflow")?;
        if pc > 0 { done(fields)?; } else if !fields.is_empty() { return Err("initial effect has a reply".into()); }
        let action = match (self.kind.as_str(), pc) {
            ("mavedb_score_interpretation", 0) => Effect::host("score.apply", &[]),
            ("mavedb_score_interpretation", 1) => Effect::host("score.select", &["PS3"]),
            ("mavedb_score_interpretation", 2) => Effect::host("score.info", &["PS3"]),
            ("mavedb_score_interpretation", 3) => Effect::host("score.select", &["BS3"]),
            ("mavedb_score_interpretation", 4) => Effect::host("score.info", &["BS3"]),
            ("mavedb_score_interpretation", 5) => Effect::finish(&["df"]),
            ("PS3_BS3_criteria", 0) => Effect::host("criteria.result-init", &[]),
            ("PS3_BS3_criteria", 1) => Effect::host("criteria.clinvar", &["Pathogenic", "clinvar_lof"]),
            ("PS3_BS3_criteria", 2) => Effect::host("criteria.result-mask", &["clinvar_patho", "clinvar_lof"]),
            ("PS3_BS3_criteria", 3) => Effect::host("criteria.clinvar", &["Benign", "high_conf_benign"]),
            ("PS3_BS3_criteria", 4) => Effect::host("criteria.result-mask", &["clinvar_benign", "high_conf_benign"]),
            ("PS3_BS3_criteria", 5) => Effect::host("criteria.zeros", &["ps3_array"]),
            ("PS3_BS3_criteria", 6) => Effect::host("criteria.zeros", &["bs3_array"]),
            ("PS3_BS3_criteria", 7) => { self.pc = 100; Effect::host("criteria.metadata-truth", &[]) }
            ("PS3_BS3_criteria", 101) => Effect::host("criteria.read-table", &[]),
            ("PS3_BS3_criteria", 102) => Effect::host("criteria.drop-duplicates", &[]),
            ("PS3_BS3_criteria", 103) => Effect::host("criteria.metadata-info", &[]),
            ("PS3_BS3_criteria", 104) => Effect::host("criteria.score", &[]),
            ("PS3_BS3_criteria", 105) => Effect::host("criteria.combine", &["PS3", "ps3_criteria"]),
            ("PS3_BS3_criteria", 106) => { self.pc = 200; Effect::host("criteria.combine", &["BS3", "bs3_criteria"]) }
            ("PS3_BS3_criteria", 151) => Effect::host("criteria.alias", &["clinvar_lof", "ps3_criteria"]),
            ("PS3_BS3_criteria", 152) => { self.pc = 200; Effect::host("criteria.alias", &["high_conf_benign", "bs3_criteria"]) }
            ("PS3_BS3_criteria", 201) => Effect::host("criteria.assign", &["ps3_array", "ps3_criteria", "3"]),
            ("PS3_BS3_criteria", 202) => Effect::host("criteria.assign", &["bs3_array", "bs3_criteria", "3"]),
            ("PS3_BS3_criteria", 203) => Effect::host("criteria.result-array", &["PS3", "ps3_array"]),
            ("PS3_BS3_criteria", 204) => Effect::host("criteria.result-array", &["BS3", "bs3_array"]),
            ("PS3_BS3_criteria", 205) => Effect::finish(&["result_dict", "df"]),
            ("PP5_BP6_criteria", 0) => Effect::host("reported.mask", &["athogenic", "clinvar_patho", "pp5_criteria"]),
            ("PP5_BP6_criteria", 1) => Effect::host("reported.mask", &["enign", "clinvar_benign", "bp6_criteria"]),
            ("PP5_BP6_criteria", 2) => Effect::host("criteria.zeros", &["pp5_array"]),
            ("PP5_BP6_criteria", 3) => Effect::host("criteria.zeros", &["bp6_array"]),
            ("PP5_BP6_criteria", 4) => Effect::host("criteria.assign", &["pp5_array", "pp5_criteria", "1"]),
            ("PP5_BP6_criteria", 5) => Effect::host("criteria.assign", &["bp6_array", "bp6_criteria", "1"]),
            ("PP5_BP6_criteria", 6) => Effect::finish(&["pp5_array", "bp6_array"]),
            _ => return Err("invalid reported program counter".into()),
        };
        Ok(action)
    }

    pub fn resume(&mut self, fields: &[String]) -> Result<Effect, String> {
        if self.kind == "PS3_BS3_criteria" && self.pc == 100 {
            self.pc = if truth(fields)? { 101 } else { 151 };
            return self.advance(&["done".into()]);
        }
        // 200 is an internal join, with no additional host operation.
        if self.kind == "PS3_BS3_criteria" && self.pc == 200 { self.pc = 201; }
        self.advance(fields)
    }
}

#[derive(Default)]
struct RowProgram {
    pc: u16,
    index: usize,
    urns: Vec<String>,
    scores: Vec<String>,
    pvalues: Vec<String>,
    high_confs: Vec<String>,
    ps3: Vec<String>,
    bs3: Vec<String>,
    family: String,
    variant: String,
    plan: Option<evidence::Resolution>,
    conflict_index: usize,
}

fn decoded_reply(fields: &[String], tag: &str) -> Result<String, String> {
    transport::decode_text(&reply(fields, tag, 1)?[0])
}

impl RowProgram {
    fn effect(&mut self, pc: u16, name: &str, args: &[&str]) -> Effect {
        self.pc = pc;
        Effect::host(name, args)
    }
    fn next_assay(&mut self) -> Effect {
        if self.index == self.urns.len() { return self.effect(30, "row.variant", &[]); }
        self.pc = 10;
        Effect::owned("row.urn", vec![transport::encode_text(&self.urns[self.index]), self.index.to_string()])
    }
    fn append(&mut self, side: &str, interpreted: bool) -> Effect {
        if side == "PS3" { self.ps3.push(self.family.clone()); } else { self.bs3.push(self.family.clone()); }
        self.pc = if interpreted { if side == "PS3" { 23 } else { 24 } } else { 24 };
        Effect::host("row.append-evidence", &[side, if interpreted { "interpreted" } else { "pvalue" }])
    }
    fn publications(&mut self, missing: bool) -> Effect {
        self.pc = if missing { 40 } else { 50 };
        Effect::host("row.publish", &["MaveDB_PS3", if !missing && self.plan.as_ref().is_some_and(|plan| !plan.final_ps3.is_empty()) { "true" } else { "false" }])
    }
    fn after_resolution(&mut self) -> Effect {
        let plan = self.plan.as_ref().expect("source resolution exists");
        if self.conflict_index < plan.conflicts.len() {
            let conflict = &plan.conflicts[self.conflict_index];
            self.pc = 33;
            return Effect::owned("row.conflict", vec![conflict.family.to_string(), transport::ordinal_list(&conflict.ps3_ordinals),
                transport::ordinal_list(&conflict.bs3_ordinals), match conflict.winner { Some(evidence::Side::Ps3) => "PS3", Some(evidence::Side::Bs3) => "BS3", None => "tie" }.into()]);
        }
        self.pc = 35;
        Effect::owned("row.final-evidence", vec![transport::ordinal_list(&plan.final_ps3), transport::ordinal_list(&plan.final_bs3)])
    }
    fn advance(&mut self, fields: &[String]) -> Result<Effect, String> {
        match self.pc {
            0 => { if !fields.is_empty() { return Err("initial row has a reply".into()); } Ok(self.effect(1, "row.missing", &[])) }
            1 => if truth(fields)? { Ok(self.publications(true)) } else { Ok(self.effect(2, "row.text", &["MaveDB_high_conf"])) },
            2 => {
                self.high_confs = decoded_reply(fields, "text")?.split('&').map(str::to_owned).collect();
                self.pc = 3; let mut args = vec!["MaveDB_pvalue".into(), "high_confs".into()]; args.extend(self.high_confs.iter().map(|value| transport::encode_text(value)));
                Ok(Effect::owned("row.text", args))
            }
            3 => {
                self.pvalues = decoded_reply(fields, "text")?.split('&').map(str::to_owned).collect();
                self.pc = 4; let mut args = vec!["MaveDB_score".into(), "pvalues".into()]; args.extend(self.pvalues.iter().map(|value| transport::encode_text(value)));
                Ok(Effect::owned("row.text", args))
            }
            4 => {
                self.scores = decoded_reply(fields, "text")?.split('&').map(str::to_owned).collect();
                self.pc = 5; let mut args = vec!["MaveDB_urn".into(), "scores".into()]; args.extend(self.scores.iter().map(|value| transport::encode_text(value)));
                Ok(Effect::owned("row.text", args))
            }
            5 => { self.urns = decoded_reply(fields, "text")?.split('&').map(str::to_owned).collect(); self.pc = 6; Ok(Effect::owned("row.init-evidence", self.urns.iter().map(|value| transport::encode_text(value)).collect())) }
            6 => { done(fields)?; Ok(self.next_assay()) }
            10 => if truth(fields)? { Ok(self.effect(25, "row.missing-warning", &[])) }
                else { Ok(self.effect(11, "row.metadata", &[])) },
            11 => { done(fields)?; Ok(self.effect(12, "row.rationale", &[])) }
            12 => { done(fields)?; Ok(self.effect(13, "row.family", &[])) }
            13 => { self.family = decoded_reply(fields, "text")?; Ok(self.effect(14, "row.subtype", &[])) }
            14 => { done(fields)?; Ok(self.effect(15, "row.append-rationale", &[])) }
            15 => {
                done(fields)?;
                if self.index >= self.scores.len() { return Ok(self.effect(99, "raise.list-index", &["scores"])); }
                self.pc = 16; Ok(Effect::owned("row.index-score", vec![transport::encode_text(&self.scores[self.index])]))
            }
            16 => {
                done(fields)?;
                if self.index >= self.high_confs.len() { return Ok(self.effect(99, "raise.list-index", &["high_confs"])); }
                self.pc = 17; Ok(Effect::owned("row.index-high-conf", vec![transport::encode_text(&self.high_confs[self.index])]))
            }
            17 => {
                done(fields)?;
                if self.index >= self.pvalues.len() { return Ok(self.effect(99, "raise.list-index", &["pvalues"])); }
                self.pc = 18; Ok(Effect::owned("row.index-pvalue", vec![transport::encode_text(&self.pvalues[self.index])]))
            }
            18 => { done(fields)?; Ok(self.effect(19, "row.float-score", &[])) }
            19 => {
                if fields == ["ValueError"] { self.index += 1; return Ok(self.next_assay()); }
                transport::decode_number(&reply(fields, "number", 1)?[0])?;
                Ok(self.effect(20, "row.float-pvalue", &[]))
            }
            20 => {
                let value = if fields == ["ValueError"] { f64::NAN } else { transport::decode_number(&reply(fields, "number", 1)?[0])? };
                match pvalue_route(value) {
                    PvalueRoute::Interpreter => Ok(self.effect(21, "row.interpret", &[])),
                    PvalueRoute::Ps3 => Ok(self.append("PS3", false)),
                    PvalueRoute::Bs3 => Ok(self.append("BS3", false)),
                }
            }
            21 => { done(fields)?; Ok(self.effect(22, "row.interpreter-ps3", &[])) }
            22 => if truth(fields)? { Ok(self.append("PS3", true)) } else { Ok(self.effect(24, "row.interpreter-bs3", &[])) },
            23 => { done(fields)?; Ok(self.effect(24, "row.interpreter-bs3", &[])) }
            24 => {
                if fields.first().is_some_and(|value| value == "bool") {
                    if truth(fields)? { return Ok(self.append("BS3", true)); }
                } else { done(fields)?; }
                self.index += 1; Ok(self.next_assay())
            }
            25 => { done(fields)?; self.index += 1; Ok(self.next_assay()) }
            30 => {
                let parts = reply(fields, "texts", 4)?.iter().map(|value| transport::decode_text(value)).collect::<Result<Vec<_>, _>>()?;
                self.variant = parts.join(":"); Ok(self.effect(31, "row.family-sets", &[]))
            }
            31 => {
                if fields.is_empty() || fields[0] != "conflicts" { return Err("missing actual Python conflict order".into()); }
                let families: Vec<String> = self.ps3.iter().chain(&self.bs3).fold(Vec::new(), |mut result, family| {
                    if !result.contains(family) { result.push(family.clone()); } result
                });
                let order = fields[1..].iter().map(|value| transport::decode_text(value).and_then(|family|
                    families.iter().position(|item| item == &family).ok_or_else(|| "unknown conflict family".into()))).collect::<Result<Vec<_>, String>>()?;
                let rows = self.ps3.iter().map(|family| (family, evidence::Side::Ps3)).chain(self.bs3.iter().map(|family| (family, evidence::Side::Bs3)))
                    .enumerate().map(|(ordinal, (family, side))| evidence::Evidence { ordinal, family: families.iter().position(|value| value == family).expect("known family"), side, formatted: String::new() }).collect::<Vec<_>>();
                let plan = evidence::resolve(&rows, &order)?;
                let names = families.iter().map(|value| transport::encode_text(value)).collect::<Vec<_>>();
                let mut args = vec![transport::ordinal_list(&plan.clean_ps3), transport::ordinal_list(&plan.clean_bs3)]; args.extend(names);
                self.plan = Some(plan); self.pc = 32; Ok(Effect::owned("row.clean-evidence", args))
            }
            32 => { done(fields)?; Ok(self.after_resolution()) }
            33 => {
                done(fields)?;
                let plan = self.plan.as_ref().expect("resolution exists"); let conflict = &plan.conflicts[self.conflict_index];
                let family = self.ps3.iter().chain(&self.bs3).fold(Vec::<String>::new(), |mut result, family| { if !result.contains(family) { result.push(family.clone()); } result })[conflict.family].clone();
                let message = match conflict.winner {
                    Some(side) => format!("Variant {}: Family '{}' conflict resolved to {} by majority (PS3:{} vs BS3:{})", self.variant, family, if side == evidence::Side::Ps3 { "PS3" } else { "BS3" }, conflict.ps3, conflict.bs3),
                    None => format!("Variant {}: Family '{}' has tied conflict (PS3:{} vs BS3:{}), excluded from final decision", self.variant, family, conflict.ps3, conflict.bs3),
                };
                self.pc = 34; Ok(Effect::owned("log", vec!["DEBUG".into(), transport::encode_text(&message)]))
            }
            34 => { done(fields)?; self.conflict_index += 1; Ok(self.after_resolution()) }
            35 => {
                done(fields)?; let plan = self.plan.as_ref().expect("resolution exists");
                if !plan.final_ps3.is_empty() && !plan.final_bs3.is_empty() { Ok(self.effect(36, "row.final-log-parts", &["both"])) }
                else if !plan.conflicts.is_empty() { Ok(self.effect(36, "row.final-log-parts", &["conflict"])) }
                else { Ok(self.publications(false)) }
            }
            36 => {
                let parts = reply(fields, "log-parts", 5)?;
                let parts = parts.iter().map(|value| transport::decode_text(value)).collect::<Result<Vec<_>, _>>()?;
                let plan = self.plan.as_ref().expect("resolution exists");
                let both = !plan.final_ps3.is_empty() && !plan.final_bs3.is_empty();
                let message = if both { format!("Variant {} has PS3 and BS3 from DIFFERENT assay families (both retained). PS3 families: {}, BS3 families: {}. PS3: [{}] | BS3: [{}]", self.variant, parts[0], parts[1], parts[3], parts[4]) }
                else { format!("Variant {} had conflicts in families {}, resolved. Final: PS3={}, BS3={}. PS3: [{}] | BS3: [{}]", self.variant, parts[2], if plan.final_ps3.is_empty() { "False" } else { "True" }, if plan.final_bs3.is_empty() { "False" } else { "True" }, parts[3], parts[4]) };
                self.pc = 37; Ok(Effect::owned("log", vec![if both { "INFO" } else { "DEBUG" }.into(), transport::encode_text(&message)]))
            }
            37 => { done(fields)?; Ok(self.publications(false)) }
            40 => { done(fields)?; Ok(self.effect(41, "row.publish", &["MaveDB_BS3", "false"])) }
            41 => { done(fields)?; Ok(self.effect(42, "row.publish", &["MaveDB_score_interpretation", "nan"])) }
            42 | 53 => { done(fields)?; Ok(Effect::finish(&["row"])) }
            50 => {
                done(fields)?;
                let value = if self.plan.as_ref().is_some_and(|plan| !plan.final_bs3.is_empty()) { "true" } else { "false" };
                Ok(self.effect(51, "row.publish", &["MaveDB_BS3", value]))
            }
            51 => { done(fields)?; Ok(self.effect(52, "row.rationale-texts", &[])) }
            52 => {
                if fields.is_empty() || fields[0] != "texts" { return Err("invalid rationale formatting reply".into()); }
                let values = fields[1..].iter().map(|value| transport::decode_text(value)).collect::<Result<Vec<_>, _>>()?;
                self.pc = 53; Ok(Effect::owned("row.publish-text", vec!["MaveDB_score_interpretation".into(), transport::encode_text(&values.join("&"))]))
            }
            _ => Err("invalid row source frontier".into()),
        }
    }
}

#[cfg(test)]
mod tests;
