use priva_tabular::acmg_scoring::{
    calculate_posterior_probability::{Probability, calculate_posterior_probability},
    create_criteria_summary::create_criteria_summary,
    fix_clingen_term::fix_clingen_term,
    translate_strength_level::translate_strength_level,
};

#[test]
fn preserves_clingen_suffix_splitting_quirks() {
    assert_eq!(fix_clingen_term("PVS1_Very_extra"), "PVS1_Very Strong");
    assert_eq!(fix_clingen_term("BA1_Standalone"), "BA1_Stand Alone");
    assert_eq!(fix_clingen_term("PM2_Supporting_extra"), "PM2_Supporting");
}

#[test]
fn bare_criterion_takes_priority_over_an_earlier_suffixed_token() {
    assert_eq!(
        translate_strength_level("PS1", Some("PS1_Supporting;PS1")).unwrap(),
        3
    );
    assert_eq!(
        translate_strength_level("PS1", Some("PS1_Supporting;PS1_Strong")).unwrap(),
        1
    );
}

#[test]
fn substring_false_positive_still_has_source_failure() {
    assert_eq!(
        translate_strength_level("PS1", Some("XPS1")).unwrap_err().0,
        "IndexError"
    );
    assert_eq!(
        translate_strength_level("PS1", Some("PS1_Unknown"))
            .unwrap_err()
            .0,
        "KeyError"
    );
    assert_eq!(translate_strength_level("PS1", None).unwrap(), 0);
}

#[test]
fn clingen_replaces_without_merging_or_reordering() {
    assert_eq!(
        create_criteria_summary(&[("PP3", 1), ("PM2", 1)], Some(" BP4 , PS1 ")).unwrap(),
        Some("BP4;PS1".into())
    );
    assert_eq!(
        create_criteria_summary(&[("PP3", 1), ("PM2", 1)], None).unwrap(),
        Some("PP3;PM2_Supporting".into())
    );
}

#[test]
fn whitespace_only_is_silence_but_commas_replace_with_no_codes() {
    assert_eq!(
        create_criteria_summary(&[("PP3", 1)], Some("\u{001c}\u{001f}")).unwrap(),
        Some("PP3".into())
    );
    assert_eq!(
        create_criteria_summary(&[("PP3", 1)], Some(", ,")).unwrap(),
        None
    );
}

#[test]
fn malformed_source_strength_fails_before_valid_clingen_override() {
    assert_eq!(
        create_criteria_summary(&[("PP3", 6)], Some("PS1"))
            .unwrap_err()
            .0,
        "TypeError"
    );
}

#[test]
fn original_special_case_discards_extra_very_strong_counts() {
    assert_eq!(
        calculate_posterior_probability(&[("PVS1", 4), ("PS1", 4), ("PM1", 2)]),
        calculate_posterior_probability(&[("PVS1", 4), ("PM1", 2)])
    );
}

#[test]
fn source_ignores_pathogenic_five_and_benign_four() {
    assert_eq!(
        calculate_posterior_probability(&[("PS1", 5), ("BS1", 4)]),
        (Probability::Float(0.1), "Uncertain Significance")
    );
}

#[test]
fn standalone_benign_overrides_pathogenic_evidence_last() {
    assert_eq!(
        calculate_posterior_probability(&[("PVS1", 4), ("PS1", 3), ("BA1", 5)]),
        (Probability::Integer(0), "Benign")
    );
}
