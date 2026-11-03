//! All thirteen source assay families and twenty-five categorical subtypes.
//!
//! Standardization is a suspended effect: Python's float conversion and exact
//! `math.log(x, 2)` / `math.log(1.0 / kd, 10)` run at their original positions.
//! A failed effect never produces a Binding directionality mutation.

pub const FAMILIES: [(&str, &[&str]); 13] = [
    ("Stability_K50", &["LOG10_K50_STANDARDIZED", "LOG10_K50_RELATIVE", "K50_RAW"]),
    ("Stability_DG", &["DDG_KCAL_MOL", "DG_KCAL_MOL"]),
    ("Abundance", &["VAMP_SEQ", "CELL_SURFACE_TRAFF", "TOTAL_EXPRESSION"]),
    ("Activity", &["ENZYME_ACTIVITY_REL", "REPORTER_ACTIVITY_REL"]),
    ("Binding", &["KD", "PKD", "KA", "KD_APPARENT"]),
    ("Fitness", &["LOG2_ENRICHMENT"]),
    ("Splicing", &["DELTA_PSI_PERCENT", "DELTA_PSI_FRACTION"]),
    ("Electrophysiology", &["V1_2", "CURRENT_DENSITY"]),
    ("Complementation", &["YEAST", "MAMMALIAN_CELL"]),
    ("Localization", &["CELLULAR_LOCALIZATION"]),
    ("Aggregation", &["AGGREGATION_METRIC"]),
    ("Interaction", &["TWO_HYBRID"]),
    ("MPRA", &["REPORTER_ACTIVITY"]),
];

#[derive(Clone, Debug)]
pub struct Spec {
    pub family: String,
    pub subtype: String,
    pub directionality: String,
    pub scale: String,
    pub lo: f64,
    pub hi: f64,
    pub band_text: String,
    pub units: String,
    pub flexibility: bool,
}

#[derive(Clone, Debug)]
pub struct Interpretation {
    pub call: String,
    pub standardized: Option<f64>,
    pub reason: String,
    pub notes: String,
    /// Some source returns publish the raw units object, not its string form.
    pub notes_are_units: bool,
    /// The host must publish this to the SAME AssaySpec before returning.
    pub make_higher: bool,
}

impl Interpretation {
    fn na(reason: &str, notes: &str) -> Self {
        Self { call: "NA".into(), standardized: None, reason: reason.into(), notes: notes.into(), notes_are_units: false, make_higher: false }
    }

    fn units(reason: &str) -> Self {
        let mut result = Self::na(reason, "");
        result.notes_are_units = true;
        result
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RawEffect { DgNegate, DgIdentity, BindingKd, BindingIdentity }

#[derive(Clone, Debug)]
pub enum Stage {
    Return(Interpretation),
    Standardize,
    Raw(RawEffect),
    BindingLog10,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum StandardizeEffect { None, Float, Log2 }

pub fn standardize_effect(score_is_none: bool, scale: &str) -> StandardizeEffect {
    if score_is_none { return StandardizeEffect::None; }
    match scale {
        "CENTERED_0" | "PERCENT" | "FRACTION" => StandardizeEffect::Float,
        "WT_EQ_1" => StandardizeEffect::Log2,
        _ => StandardizeEffect::None,
    }
}

pub fn log2_needed(value: Option<f64>) -> bool {
    value.is_some_and(|value| !(value <= 0.0))
}

pub fn valid_spec(spec: &Spec) -> bool {
    FAMILIES.iter().find(|(family, _)| *family == spec.family)
        .is_some_and(|(_, subtypes)| subtypes.contains(&spec.subtype.as_str())) && spec.lo < spec.hi
}

pub fn call_from_band(value: Option<f64>, spec: &Spec, make_higher: bool) -> &'static str {
    let Some(value) = value else { return "Unknown"; };
    if spec.lo <= value && value <= spec.hi { return "Neutral"; }
    match if make_higher { "HIGHER_IS_FUNCTION" } else { spec.directionality.as_str() } {
        "HIGHER_IS_FUNCTION" => if value > spec.hi { "GoF" } else { "LoF" },
        "LOWER_IS_FUNCTION" if spec.flexibility => if value < spec.lo { "GoF" } else { "LoF" },
        "LOWER_IS_FUNCTION" => if value < spec.lo { "LoF" } else { "GoF" },
        _ => "Unknown",
    }
}

/// Registry lookup precedes per-family validation in interpret_score.
pub fn begin(spec: &Spec, gene_rules_present: bool) -> Stage {
    if !FAMILIES.iter().any(|(family, _)| *family == spec.family) {
        return Stage::Return(Interpretation::na("No interpreter", &spec.family));
    }
    if !valid_spec(spec) { return Stage::Return(Interpretation::na("Invalid spec", "")); }
    if spec.family == "Stability_K50" && spec.subtype == "K50_RAW" {
        return Stage::Return(Interpretation::na("K50 RAW requires transform", &spec.subtype));
    }
    if spec.family == "Electrophysiology" && !gene_rules_present {
        return Stage::Return(Interpretation::na("Gene-specific rule required", ""));
    }
    Stage::Standardize
}

pub fn after_standardize(spec: &Spec, value: Option<f64>, gof_dir: &str) -> Stage {
    if value.is_none() {
        match spec.family.as_str() {
            "Stability_DG" => return Stage::Raw(if spec.subtype == "DDG_KCAL_MOL" { RawEffect::DgNegate } else { RawEffect::DgIdentity }),
            "Binding" => return Stage::Raw(if matches!(spec.subtype.as_str(), "KD" | "KD_APPARENT") { RawEffect::BindingKd } else { RawEffect::BindingIdentity }),
            _ => {
                let label = match spec.family.as_str() {
                    "Stability_K50" => "K50", "Electrophysiology" => "Ephys", family => family,
                };
                return Stage::Return(Interpretation::units(&format!("{label} not standardized")));
            }
        }
    }
    Stage::Return(finish(spec, value, false, false, gof_dir))
}

pub fn after_binding_float(value: f64) -> Stage {
    if value <= 0.0 { Stage::Return(Interpretation::na("Kd <= 0 invalid", "")) }
    else { Stage::BindingLog10 }
}

pub fn after_raw(spec: &Spec, value: f64, gof_dir: &str) -> Interpretation {
    let binding = spec.family == "Binding";
    finish(spec, Some(value), true, binding && spec.directionality == "NON_MONOTONIC", gof_dir)
}

fn finish(spec: &Spec, value: Option<f64>, raw_used: bool, make_higher: bool, gof_dir: &str) -> Interpretation {
    let band = &spec.band_text;
    let subtype = &spec.subtype;
    let flex = if spec.flexibility { " | flexibility_req=True" } else { "" };
    let (call, reason, notes) = match spec.family.as_str() {
        "Stability_K50" => (call_from_band(value, spec, false), "K50 heuristic band".into(), format!("{band} | {}{flex}", spec.units)),
        "Stability_DG" => (call_from_band(value, spec, false), if raw_used { "ΔG/ΔΔG transform" } else { "ΔG/ΔΔG standardized" }.into(), format!("{band} | {}{flex}", spec.units)),
        "Binding" if spec.directionality == "NON_MONOTONIC" && !make_higher => ("NA", "Binding directionality unknown".into(), spec.units.clone()),
        "Binding" => (call_from_band(value, spec, make_higher), "Binding band".into(), format!("{band} | {subtype}")),
        "Fitness" => (call_from_band(value, spec, false), "Fitness band".into(), band.clone()),
        "Splicing" => (call_from_band(value, spec, false), "Splicing band".into(), format!("{band} | ΔPSI mapping is gene/context-specific")),
        "Electrophysiology" => {
            let value = value.expect("standardized Ephys stage");
            if spec.lo <= value && value <= spec.hi {
                ("Neutral", "Ephys neutral band".into(), subtype.clone())
            } else {
                let call = match gof_dir {
                    "more_negative" => if value < spec.lo { "GoF" } else { "LoF" },
                    "more_positive" => if value > spec.hi { "GoF" } else { "LoF" },
                    _ => "Unknown",
                };
                (call, "Ephys rule".into(), format!("{band} | {subtype} | {gof_dir}"))
            }
        }
        family => (call_from_band(value, spec, false), format!("{family} band"), format!("{band} | {subtype}")),
    };
    let notes_are_units = spec.family == "Binding" && spec.directionality == "NON_MONOTONIC" && !make_higher;
    Interpretation { call: call.into(), standardized: value, reason, notes, notes_are_units, make_higher }
}

/// Selection after source-ordered Python token conversion. None means the
/// source must continue to q05/q95; this is not a synthesized neutral band.
pub fn explicit_band(lo: Option<f64>, hi: Option<f64>, flexibility: bool) -> Option<(f64, f64)> {
    match (lo, hi) {
        (None, Some(hi)) => Some(if flexibility { (f64::NEG_INFINITY, hi) } else { (hi, f64::INFINITY) }),
        (Some(lo), None) => Some(if flexibility { (lo, f64::INFINITY) } else { (f64::NEG_INFINITY, lo) }),
        (Some(lo), Some(hi)) if lo < hi => Some((lo, hi)),
        _ => None,
    }
}

pub fn percentile_band(lo: f64, hi: f64) -> Option<(f64, f64)> {
    if !lo.is_nan() && !hi.is_nan() && lo < hi { Some((lo, hi)) } else { None }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum QualityRoute { Ambiguous, GeneRules, Band }

pub fn quality_route(ambiguous: bool, nonmonotonic: bool, needs_gene_rules: bool) -> QualityRoute {
    if ambiguous { QualityRoute::Ambiguous }
    else if nonmonotonic || needs_gene_rules { QualityRoute::GeneRules }
    else { QualityRoute::Band }
}

pub fn metadata_criteria(call: &str) -> (bool, bool) {
    (matches!(call, "LoF" | "GoF"), call == "Neutral")
}

/// Metadata selection and parsing have native control, but mature pandas
/// selection, CPython parsing/conversion and actual Enum/AssaySpec effects.
pub struct MetadataProgram {
    parser: bool,
    pc: u16,
    lo: Option<f64>,
    hi: Option<f64>,
    family: String,
    subtype: String,
    pending_reason: String,
    publish_index: usize,
}

impl MetadataProgram {
    pub fn new(kind: &str) -> Self {
        Self { parser: kind == "parse_neutral_band_from_metadata", pc: 0, lo: None, hi: None,
            family: String::new(), subtype: String::new(), pending_reason: String::new(), publish_index: 0 }
    }
    fn effect(&mut self, pc: u16, name: &str, args: &[&str]) -> super::Effect {
        self.pc = pc; super::Effect::host(name, args)
    }
    fn band(&mut self, result: Option<(f64, f64)>) -> super::Effect {
        self.pc = 90;
        let (lo, hi) = result.map_or((None, None), |(lo, hi)| (Some(lo), Some(hi)));
        super::Effect::owned("band.result", vec![super::transport::encode_number(lo), super::transport::encode_number(hi)])
    }
    fn fail(&mut self, reason: &str) -> super::Effect {
        self.pending_reason = reason.into();
        self.effect(80, "pick.failure-parts", &[reason])
    }
    fn publish(&mut self) -> super::Effect {
        let keys = ["metadata_flags", "urn", "n_syn_used", "band_source", "assay_family", "subtype", "neutral_band_parsed", "MaveDB_BS3", "MaveDB_PS3"];
        if self.publish_index < keys.len() {
            let key = keys[self.publish_index]; self.publish_index += 1;
            self.effect(30, "pick.publish", &[key, "source"])
        } else { self.effect(31, "pick.call", &[]) }
    }
    pub fn advance(&mut self, fields: &[String]) -> Result<super::Effect, String> {
        if self.parser { return self.parse(fields); }
        use super::{Effect, reply, truth};
        match self.pc {
            0 => { if !fields.is_empty() { return Err("initial metadata selector has a reply".into()); } Ok(self.effect(1, "pick.mode", &[])) }
            1 => if truth(fields)? { Ok(self.effect(2, "pick.frame-lookup", &[])) } else { Ok(self.effect(3, "pick.dict-lookup", &[])) },
            2 => {
                if fields == ["not-found"] { Ok(self.fail("missing")) }
                else { reply(fields, "done", 0)?; Ok(self.effect(5, "pick.flags", &[])) }
            }
            3 => { reply(fields, "done", 0)?; Ok(self.effect(4, "pick.meta-truth", &[])) }
            4 => if truth(fields)? { Ok(self.effect(5, "pick.flags", &[])) } else { Ok(self.fail("missing")) },
            5 => { reply(fields, "done", 0)?; Ok(self.effect(6, "pick.flag-truth", &["ambiguous"])) }
            6 => if truth(fields)? { Ok(self.fail("ambiguous")) } else { Ok(self.effect(7, "pick.flag-truth", &["nonmonotonic"])) },
            7 => if truth(fields)? { Ok(self.fail("rules")) } else { Ok(self.effect(8, "pick.flag-truth", &["needs_gene_rules"])) },
            8 => if truth(fields)? { Ok(self.fail("rules")) } else { Ok(self.effect(9, "pick.parse-band", &[])) },
            9 => {
                let values = reply(fields, "band", 2)?;
                self.lo = super::transport::decode_optional(&values[0])?; self.hi = super::transport::decode_optional(&values[1])?;
                if self.lo.is_none() || self.hi.is_none() { Ok(self.fail("band")) }
                else { Ok(self.effect(10, "pick.band-ready", &[])) }
            }
            10 => { reply(fields, "done", 0)?; Ok(self.effect(11, "pick.family-subtype", &[])) }
            11 => {
                if fields.first().is_some_and(|value| value == "construct-error") { return Ok(self.fail("construct")); }
                let values = reply(fields, "family-subtype", 4)?;
                if values[0] != "str" { return Ok(self.fail("family")); }
                self.family = super::transport::decode_text(&values[1])?;
                if !FAMILIES.iter().any(|(family, _)| *family == self.family) { return Ok(self.fail("family")); }
                self.subtype = if values[2] == "str" { super::transport::decode_text(&values[3])? } else { String::new() };
                let family = self.family.clone(); Ok(self.effect(12, "pick.subtype", &[&family]))
            }
            12 => {
                if fields == ["KeyError"] { return Ok(self.fail("subtype")); }
                if fields == ["construct-error"] { return Ok(self.fail("construct")); }
                reply(fields, "done", 0)?;
                let expected = FAMILIES.iter().find(|(family, _)| *family == self.family).expect("family checked");
                if !expected.1.contains(&self.subtype.as_str()) { return Err("host Enum disagrees with native subtype inventory".into()); }
                Ok(self.effect(13, "pick.enum-strings", &[]))
            }
            13 => {
                if fields == ["construct-error"] { return Ok(self.fail("construct")); }
                reply(fields, "done", 0)?; Ok(self.effect(14, "pick.enums", &[]))
            }
            14 => {
                if fields == ["KeyError"] { return Ok(self.fail("enum")); }
                if fields == ["construct-error"] { return Ok(self.fail("construct")); }
                reply(fields, "done", 0)?; Ok(self.effect(15, "pick.construct", &[]))
            }
            15 => {
                if fields == ["construct-error"] { return Ok(self.fail("construct")); }
                reply(fields, "done", 0)?; Ok(self.effect(16, "pick.interpret", &[]))
            }
            16 => { reply(fields, "done", 0)?; Ok(self.publish()) }
            30 => { reply(fields, "done", 0)?; Ok(self.publish()) }
            31 => {
                let call = super::transport::decode_text(&reply(fields, "text", 1)?[0])?;
                let (ps3, bs3) = metadata_criteria(&call);
                if ps3 { Ok(self.effect(32, "pick.publish", &["MaveDB_PS3", "true"])) }
                else if bs3 { Ok(self.effect(32, "pick.publish", &["MaveDB_BS3", "true"])) }
                else { Ok(Effect::finish(&["result"])) }
            }
            32 => { reply(fields, "done", 0)?; Ok(Effect::finish(&["result"])) }
            80 => {
                if fields.is_empty() || fields[0] != "failure-parts" { return Err("invalid metadata error formatting parts".into()); }
                let parts = fields[1..].iter().map(|value| super::transport::decode_text(value)).collect::<Result<Vec<_>, _>>()?;
                let expected = match self.pending_reason.as_str() { "family" => 1, "subtype" => 2, "enum" => 3, "construct" => 2, "band" => 3, _ => 0 };
                if parts.len() != expected { return Err("invalid metadata error part count".into()); }
                let reason = match self.pending_reason.as_str() {
                    "missing" => "URN not found in metadata".into(), "ambiguous" => "Dataset marked ambiguous".into(),
                    "rules" => "Gene-specific rules required (not supported)".into(), "band" => "Unable to determine neutral band".into(),
                    "family" => format!("Unknown assay family: {}", parts[0]),
                    "subtype" => format!("Unknown subtype {} for {}", parts[0], parts[1]),
                    "enum" => format!("Invalid enum value: {}", parts[0]),
                    "construct" => format!("Failed to construct AssaySpec: {}", parts[0]),
                    _ => return Err("unknown source metadata failure".into()),
                };
                let note = match self.pending_reason.as_str() {
                    "band" => Some(format!("neutral_bands={}, q05={}, q95={}", parts[0], parts[1], parts[2])),
                    "enum" => Some(format!("directionality={}, scale={}", parts[1], parts[2])),
                    "construct" => Some(parts[1].clone()), _ => None,
                };
                self.pc = 81;
                Ok(Effect::owned("pick.na", vec![self.pending_reason.clone(), super::transport::encode_text(&reason), note.map_or_else(|| "-".into(), |value| super::transport::encode_text(&value))]))
            }
            81 => { reply(fields, "done", 0)?; Ok(Effect::finish(&["result"])) }
            _ => Err("invalid metadata selector source frontier".into()),
        }
    }

    fn parse(&mut self, fields: &[String]) -> Result<super::Effect, String> {
        use super::{Effect, reply, truth};
        match self.pc {
            0 => { if !fields.is_empty() { return Err("initial band parser has a reply".into()); } Ok(self.effect(1, "band.initial-text", &[])) }
            1 => {
                let text = super::transport::decode_text(&reply(fields, "text", 1)?[0])?;
                if ["[NA, NA]", "NA", "", "None", "nan"].contains(&text.as_str()) { return Ok(self.effect(20, "band.percentiles", &[])); }
                self.pc = 2; Ok(Effect::owned("band.strip-tokens", text.split(',').map(super::transport::encode_text).collect()))
            }
            2 => {
                if fields.is_empty() || fields[0] != "texts" { return Err("invalid stripped band tokens".into()); }
                let mut arguments = Vec::new();
                for value in &fields[1..] {
                    let token = super::transport::decode_text(value)?;
                    arguments.extend([if ["NA", "nan", "None", ""].contains(&token.as_str()) { "nan" } else { "float" }.into(), value.clone()]);
                }
                self.pc = 3; Ok(Effect::owned("band.float-tokens", arguments))
            }
            3 => {
                if fields == ["parse-error"] { return Ok(self.effect(20, "band.percentiles", &[])); }
                if fields.is_empty() || fields[0] != "numbers" { return Err("invalid parsed band numbers".into()); }
                if fields.len() != 3 { return Ok(self.effect(20, "band.percentiles", &[])); }
                self.lo = super::transport::decode_optional(&fields[1])?.filter(|value| !value.is_nan());
                self.hi = super::transport::decode_optional(&fields[2])?.filter(|value| !value.is_nan());
                Ok(self.effect(4, "band.flexibility", &[]))
            }
            4 => {
                if fields == ["parse-error"] { return Ok(self.effect(20, "band.percentiles", &[])); }
                reply(fields, "done", 0)?;
                match (self.lo, self.hi) {
                    (None, Some(_)) | (Some(_), None) => Ok(self.effect(5, "band.flex-truth", &[])),
                    (Some(lo), Some(hi)) if lo < hi => Ok(self.band(Some((lo, hi)))),
                    _ => Ok(self.effect(20, "band.percentiles", &[])),
                }
            }
            5 => {
                if fields == ["parse-error"] { return Ok(self.effect(20, "band.percentiles", &[])); }
                let flexibility = truth(fields)?; Ok(self.band(explicit_band(self.lo, self.hi, flexibility)))
            }
            20 => {
                let values = reply(fields, "none-bits", 2)?;
                if values == ["0", "0"] { Ok(self.effect(21, "band.float-percentiles", &[])) }
                else if values.iter().all(|value| value == "0" || value == "1") { Ok(self.band(None)) }
                else { Err("invalid percentile None flags".into()) }
            }
            21 => {
                if fields == ["parse-error"] { return Ok(self.band(None)); }
                let values = reply(fields, "numbers", 2)?;
                Ok(self.band(percentile_band(super::transport::decode_number(&values[0])?, super::transport::decode_number(&values[1])?)))
            }
            90 => { reply(fields, "done", 0)?; Ok(Effect::finish(&["_band_result"])) }
            _ => Err("invalid neutral-band parser source frontier".into()),
        }
    }
}
