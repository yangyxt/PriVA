use super::evidence::{Evidence, Side, resolve};
use super::interpreter::{self, RawEffect, Spec, Stage};
use super::{PvalueRoute, combine_mask, pvalue_route, strength_table};

fn spec(family: &str, subtype: &str) -> Spec {
    Spec { family: family.into(), subtype: subtype.into(), directionality: "HIGHER_IS_FUNCTION".into(),
        scale: "CENTERED_0".into(), lo: -0.5, hi: 0.5, band_text: "(-0.5, 0.5)".into(), units: "unit".into(), flexibility: false }
}

#[test]
fn complete_subtype_inventory() {
    assert_eq!(interpreter::FAMILIES.len(), 13);
    assert_eq!(interpreter::FAMILIES.iter().map(|(_, values)| values.len()).sum::<usize>(), 25);
    for (family, subtypes) in interpreter::FAMILIES {
        for subtype in subtypes { assert!(interpreter::valid_spec(&spec(family, subtype))); }
    }
}

#[test]
fn inclusive_bands_and_nan_branching() {
    let mut value = spec("Abundance", "VAMP_SEQ");
    assert_eq!(interpreter::call_from_band(Some(-0.5), &value, false), "Neutral");
    assert_eq!(interpreter::call_from_band(Some(0.5), &value, false), "Neutral");
    assert_eq!(interpreter::call_from_band(Some(f64::NAN), &value, false), "LoF");
    value.directionality = "LOWER_IS_FUNCTION".into();
    assert_eq!(interpreter::call_from_band(Some(f64::NAN), &value, false), "GoF");
    value.flexibility = true;
    assert_eq!(interpreter::call_from_band(Some(f64::NAN), &value, false), "LoF");
}

#[test]
fn binding_mutation_only_after_completed_raw_effect() {
    let mut value = spec("Binding", "KD"); value.scale = "RAW".into(); value.directionality = "NON_MONOTONIC".into();
    assert!(matches!(interpreter::after_standardize(&value, None, ""), Stage::Raw(RawEffect::BindingKd)));
    assert!(matches!(interpreter::after_binding_float(0.0), Stage::Return(ref result) if !result.make_higher && result.reason == "Kd <= 0 invalid"));
    assert!(matches!(interpreter::after_binding_float(1.0), Stage::BindingLog10));
    let result = interpreter::after_raw(&value, -0.0, "");
    assert!(result.make_higher); assert_eq!(result.call, "Neutral");
    assert_eq!(result.standardized.unwrap().to_bits(), (-0.0_f64).to_bits());
    assert_eq!(value.directionality, "NON_MONOTONIC");
}

#[test]
fn ephys_requires_rules_before_standardization() {
    let value = spec("Electrophysiology", "V1_2");
    assert!(matches!(interpreter::begin(&value, false), Stage::Return(ref result) if result.reason == "Gene-specific rule required"));
    assert!(matches!(interpreter::begin(&value, true), Stage::Standardize));
    assert!(matches!(interpreter::after_standardize(&value, Some(-1.0), "more_negative"), Stage::Return(ref result) if result.call == "GoF"));
}

#[test]
fn one_sided_band_direction_is_source_defined() {
    assert_eq!(interpreter::explicit_band(None, Some(0.5), false), Some((0.5, f64::INFINITY)));
    assert_eq!(interpreter::explicit_band(None, Some(0.5), true), Some((f64::NEG_INFINITY, 0.5)));
    assert_eq!(interpreter::explicit_band(Some(0.0), None, false), Some((f64::NEG_INFINITY, 0.0)));
    assert_eq!(interpreter::explicit_band(Some(0.0), None, true), Some((0.0, f64::INFINITY)));
    assert_eq!(interpreter::percentile_band(f64::NAN, 1.0), None);
}

#[test]
fn pvalue_boundary_and_infinities() {
    assert_eq!(pvalue_route(0.05), PvalueRoute::Bs3);
    assert_eq!(pvalue_route(f64::NEG_INFINITY), PvalueRoute::Ps3);
    assert_eq!(pvalue_route(f64::INFINITY), PvalueRoute::Bs3);
    assert_eq!(pvalue_route(f64::NAN), PvalueRoute::Interpreter);
}

#[test]
fn family_resolution_uses_given_conflict_order_and_original_list_order() {
    let rows = [(0, Side::Ps3), (1, Side::Ps3), (0, Side::Bs3), (1, Side::Bs3), (1, Side::Ps3), (2, Side::Bs3)]
        .into_iter().enumerate().map(|(ordinal, (family, side))| Evidence { ordinal, family, side, formatted: ordinal.to_string() }).collect::<Vec<_>>();
    let value = resolve(&rows, &[1, 0]).unwrap();
    assert_eq!(value.final_ps3, [1, 4]); assert_eq!(value.final_bs3, [5]);
    assert_eq!(value.conflicts[0].family, 1); assert_eq!(value.conflicts[1].winner, None);
    assert!(resolve(&rows, &[1]).is_err()); assert!(resolve(&rows, &[1, 1, 0]).is_err());
}

#[test]
fn pp5_excludes_only_the_given_clinvar_mask() {
    assert_eq!(combine_mask(&[true, true, false], &[false, true, false], true).unwrap(), [true, false, false]);
    let table = strength_table(&[true, false], &[false, true], true).unwrap();
    assert_eq!(table.column("PP5").unwrap().i64().unwrap().get(0), Some(1));
    assert_eq!(table.column("BP6").unwrap().i64().unwrap().get(1), Some(1));
}

#[test]
fn session_rejects_duplicate_ordinals_and_malformed_payloads() {
    let mut output = Vec::new();
    assert!(super::transport::serve(&b"0\thello\treported-criteria-v1\n0\tpvalue\t3fa999999999999a\n"[..], &mut output).is_err());
    assert_eq!(String::from_utf8(output).unwrap(), "0\tok\treported-criteria-v1\n");
    assert!(super::transport::request("pvalue", &["NaN"]).is_err());
    assert!(super::transport::request("contains", &["Pathogenic", "ff"]).is_err());
}

fn response(values: &[&str]) -> Vec<String> { values.iter().map(|value| (*value).into()).collect() }
fn action(value: &super::Effect) -> &str { &value.0[1] }

#[test]
fn complete_ps3_without_metadata_keeps_each_publication_order() {
    let mut program = super::Program::new("PS3_BS3_criteria").unwrap();
    let mut result = program.advance(&[]).unwrap();
    let mut effects = Vec::new();
    loop {
        if result.0[0] == "return" { break; }
        effects.push(result.0[1..].join(":"));
        let reply = if action(&result) == "criteria.metadata-truth" { response(&["bool", "0"]) } else { response(&["done"]) };
        result = program.resume(&reply).unwrap();
    }
    assert_eq!(effects, ["criteria.result-init", "criteria.clinvar:Pathogenic:clinvar_lof", "criteria.result-mask:clinvar_patho:clinvar_lof",
        "criteria.clinvar:Benign:high_conf_benign", "criteria.result-mask:clinvar_benign:high_conf_benign", "criteria.zeros:ps3_array", "criteria.zeros:bs3_array",
        "criteria.metadata-truth", "criteria.alias:clinvar_lof:ps3_criteria", "criteria.alias:high_conf_benign:bs3_criteria", "criteria.assign:ps3_array:ps3_criteria:3",
        "criteria.assign:bs3_array:bs3_criteria:3", "criteria.result-array:PS3:ps3_array", "criteria.result-array:BS3:bs3_array"]);
    assert_eq!(result.0, response(&["return", "result_dict", "df"]));
}

#[test]
fn complete_ps3_metadata_route_loads_deduplicates_logs_then_scores() {
    let mut program = super::Program::new("PS3_BS3_criteria").unwrap();
    let mut result = program.advance(&[]).unwrap();
    let mut effects = Vec::new();
    while result.0[0] != "return" {
        effects.push(action(&result).to_owned());
        let reply = if action(&result) == "criteria.metadata-truth" { response(&["bool", "1"]) } else { response(&["done"]) };
        result = program.resume(&reply).unwrap();
    }
    assert_eq!(&effects[8..14], ["criteria.read-table", "criteria.drop-duplicates", "criteria.metadata-info", "criteria.score", "criteria.combine", "criteria.combine"]);
}

#[test]
fn missing_row_returns_only_after_three_ordered_writes() {
    let mut program = super::Program::new("mavedb_interpretation_per_row").unwrap();
    assert_eq!(action(&program.advance(&[]).unwrap()), "row.missing");
    let first = program.resume(&response(&["bool", "1"])).unwrap();
    assert_eq!(first.0[1..], response(&["row.publish", "MaveDB_PS3", "false"]));
    let second = program.resume(&response(&["done"])).unwrap();
    assert_eq!(second.0[1..], response(&["row.publish", "MaveDB_BS3", "false"]));
    let third = program.resume(&response(&["done"])).unwrap();
    assert_eq!(third.0[1..], response(&["row.publish", "MaveDB_score_interpretation", "nan"]));
    assert_eq!(program.resume(&response(&["done"])).unwrap().0, response(&["return", "row"]));
}

#[test]
fn neutral_band_native_priority_retains_unused_flexibility_get() {
    let mut program = super::Program::new("parse_neutral_band_from_metadata").unwrap();
    assert_eq!(action(&program.advance(&[]).unwrap()), "band.initial-text");
    let text = super::transport::encode_text("NA, NA");
    assert_eq!(action(&program.resume(&response(&["text", &text])).unwrap()), "band.strip-tokens");
    let na = super::transport::encode_text("NA");
    assert_eq!(action(&program.resume(&response(&["texts", &na, &na])).unwrap()), "band.float-tokens");
    let nan = super::transport::encode_number(Some(f64::NAN));
    assert_eq!(action(&program.resume(&response(&["numbers", &nan, &nan])).unwrap()), "band.flexibility");
    assert_eq!(action(&program.resume(&response(&["done"])).unwrap()), "band.percentiles");
    let result = program.resume(&response(&["none-bits", "1", "1"])).unwrap();
    assert_eq!(result.0, response(&["effect", "band.result", "-", "-"]));
}

#[test]
fn metadata_ambiguous_gate_precedes_every_band_enum_and_score_effect() {
    let mut program = super::Program::new("pick_interpreter_from_metadata").unwrap();
    assert_eq!(action(&program.advance(&[]).unwrap()), "pick.mode");
    assert_eq!(action(&program.resume(&response(&["bool", "1"])).unwrap()), "pick.frame-lookup");
    assert_eq!(action(&program.resume(&response(&["done"])).unwrap()), "pick.flags");
    assert_eq!(program.resume(&response(&["done"])).unwrap().0, response(&["effect", "pick.flag-truth", "ambiguous"]));
    assert_eq!(program.resume(&response(&["bool", "1"])).unwrap().0, response(&["effect", "pick.failure-parts", "ambiguous"]));
    let result = program.resume(&response(&["failure-parts"])).unwrap();
    assert_eq!(&result.0[..3], response(&["effect", "pick.na", "ambiguous"]));
    assert_eq!(super::transport::decode_text(&result.0[3]).unwrap(), "Dataset marked ambiguous");
}

#[test]
fn public_programs_reject_wrong_reply_types_and_unknown_names() {
    assert!(super::Program::new("mask-only").is_err());
    let mut program = super::Program::new("mavedb_interpretation_per_row").unwrap();
    program.advance(&[]).unwrap();
    assert!(program.resume(&response(&["done"])).is_err());
    let mut output = Vec::new();
    assert!(super::transport::serve(&b"0\thello\treported-criteria-v1\n1\tnext\t0\tdone\n"[..], &mut output).is_err());
}

#[test]
fn typed_urn_selection_retains_first_duplicate_and_all_match_ordinals() {
    let one = super::transport::encode_text("one");
    let two = super::transport::encode_text("two");
    assert_eq!(super::transport::request("urn-mask", &[&one, "str", &two, "missing", "", "str", &one, "str", &one]).unwrap(),
        response(&["0", "2", "0", "0", "1", "1"]));
    assert_eq!(super::transport::request("urn-mask", &[&one, "missing", ""]).unwrap(), response(&["1", "-", "0"]));
    assert!(super::transport::request("urn-mask", &[&one, "missing", "00"]).is_err());
    assert!(super::transport::request("urn-mask", &[&one, "object", ""]).is_err());
}

#[test]
fn fixed_patterns_confidence_nan_and_logical_masks_have_native_values() {
    assert_eq!(super::contains(&["Likely_pathogenic".into(), "Pathogenic".into(), "Benign".into()], "Pathogenic").unwrap(), [false, true, false]);
    assert_eq!(super::high_confidence_mask(&[true, true, true, true], &[None, Some(f64::NAN), Some(2.0), Some(1.0)]).unwrap(), [false, false, true, false]);
    assert_eq!(super::logical_mask(&[false, true], &[true, false], "or").unwrap(), [true, true]);
    assert!(super::logical_mask(&[true], &[], "and").is_err());
}

#[test]
fn ephys_direction_getter_is_selected_only_after_a_real_standardized_value() {
    let none = "-";
    let zero = super::transport::encode_number(Some(0.0));
    assert_eq!(super::transport::request("gof-direction-needed", &["Electrophysiology", none]).unwrap(), response(&["0"]));
    assert_eq!(super::transport::request("gof-direction-needed", &["Electrophysiology", &zero]).unwrap(), response(&["1"]));
    assert_eq!(super::transport::request("gof-direction-needed", &["Binding", &zero]).unwrap(), response(&["0"]));
}

#[test]
fn unused_high_conf_index_failure_remains_after_rationale_and_score_index() {
    let mut program = super::Program::new("mavedb_interpretation_per_row").unwrap();
    let mut effect = program.advance(&[]).unwrap();
    let mut effects = Vec::new();
    let mut urn_count = 0;
    loop {
        effects.push(effect.0[1..].join(":"));
        if action(&effect) == "raise.list-index" { break; }
        let value = match action(&effect) {
            "row.missing" => response(&["bool", "0"]),
            "row.text" => {
                let token = match effect.0[2].as_str() { "MaveDB_high_conf" => "", "MaveDB_pvalue" => "nan&nan", "MaveDB_score" => "bad&1", "MaveDB_urn" => "one&two", _ => panic!("unexpected key") };
                vec!["text".into(), super::transport::encode_text(token)]
            }
            "row.urn" => { urn_count += 1; response(&["bool", "0"]) }
            "row.family" => vec!["text".into(), super::transport::encode_text("Activity")],
            "row.float-score" => response(&["ValueError"]),
            _ => response(&["done"]),
        };
        effect = program.resume(&value).unwrap();
    }
    assert_eq!(urn_count, 2);
    assert_eq!(&effects[effects.len() - 3..], response(&["row.append-rationale", &format!("row.index-score:{}", super::transport::encode_text("1")), "raise.list-index:high_confs"]));
}

#[test]
fn interleaved_native_contexts_preserve_unique_ids_and_reject_completed_reuse() {
    let input = b"0\thello\treported-criteria-v1\n1\tstart\tPP5_BP6_criteria\n2\tstart\tPS3_BS3_criteria\n3\tnext\t999\tdone\n";
    let mut output = Vec::new();
    assert!(super::transport::serve(&input[..], &mut output).is_err());
    let output = String::from_utf8(output).unwrap();
    assert!(output.contains("1\tok\t0\teffect\treported.mask\tathogenic\tclinvar_patho\tpp5_criteria\n"));
    assert!(output.contains("2\tok\t1\teffect\tcriteria.result-init\n"));
}
