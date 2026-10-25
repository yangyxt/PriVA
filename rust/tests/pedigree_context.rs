use priva_tabular::acmg_scoring::pedigree_context::{Family, Warning, pedigree_context};
use priva_tabular::acmg_scoring::terminal_order::{Input, InputColumn};
use priva_tabular::acmg_variant_mechanism::variant_mechanism_masks::Cell;

fn text(value: &str) -> Cell {
    Cell::Text(value.into())
}
fn column(name: &str, dtype: &str, values: Vec<Cell>) -> InputColumn {
    InputColumn {
        name: name.into(),
        dtype: dtype.into(),
        values,
    }
}
fn frame(values: Vec<Cell>, dtype: &str) -> Input<usize> {
    Input {
        index: vec![7; values.len()],
        columns: vec![column("P", dtype, values)],
        axis_metadata: "opaque index/attrs/flags".into(),
    }
}
fn ped() -> Input<usize> {
    Input {
        index: vec![9, 9],
        columns: vec![
            column("#FamilyID", "object", vec![text("F"), text("F")]),
            column("Phenotype", "object", vec![text("2"), Cell::Integer(2)]),
            column("IndividualID", "object", vec![text("P"), text("OTHER")]),
            column("Sex", "object", vec![text("1"), text("2")]),
        ],
        axis_metadata: String::new(),
    }
}
fn family() -> Family {
    Family::Value(text("F"))
}
fn data(context: &priva_tabular::acmg_scoring::pedigree_context::Context<usize>) -> Vec<bool> {
    context
        .mask
        .column("value")
        .unwrap()
        .bool()
        .unwrap()
        .into_no_null_iter()
        .collect()
}

#[test]
fn first_affected_row_uses_source_position_despite_duplicate_ped_index() {
    let context =
        pedigree_context(&frame(vec![text("1")], "object"), Some(&ped()), &family()).unwrap();
    assert_eq!(context.pedigree_position, Some(0));
    assert_eq!(context.proband.as_deref(), Some("P"));
    assert!(context.is_male);
    assert!(!context.is_female);
    assert_eq!(context.index, vec![7]);
}
#[test]
fn counts_entire_sample_string_and_ascii_one_without_gt_split() {
    let input = frame(
        vec![
            text("0/1:10"),
            text("0/0:1"),
            text("1|1"),
            text("１"),
            text("\0\n1"),
            Cell::Null,
        ],
        "object",
    );
    let context = pedigree_context(&input, Some(&ped()), &family()).unwrap();
    assert_eq!(data(&context), vec![false, true, false, false, true, false]);
    assert_eq!(context.mask_dtype, "bool");
    assert_eq!(context.axis_metadata, input.axis_metadata);
    assert!(matches!(input.columns[0].values[5], Cell::Null));
}
#[test]
fn nullable_string_retains_true_data_bits_beneath_missing_mask() {
    let context = pedigree_context(
        &frame(vec![text("11"), Cell::Null, text("1")], "string"),
        Some(&ped()),
        &family(),
    )
    .unwrap();
    assert_eq!(context.mask_dtype, "boolean");
    assert_eq!(data(&context), vec![false, true, true]);
    assert_eq!(
        context
            .mask
            .column("is_null")
            .unwrap()
            .bool()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<_>>(),
        vec![false, true, false]
    );
    assert!(context.inherit_frame_metadata);
}
#[test]
fn disabled_pedigree_branch_does_not_validate_unused_ped_or_numeric_sample() {
    let input = frame(vec![], "UInt64");
    let mut pedigree = ped();
    pedigree.columns.clear();
    let context = pedigree_context(&input, Some(&pedigree), &Family::None).unwrap();
    assert_eq!(context.mask_dtype, "object");
    assert_eq!(context.warning, Some(Warning::NoPedigree));
    assert!(!context.inherit_frame_metadata);
    assert!(context.proband.is_none());
    let context =
        pedigree_context(&frame(vec![Cell::Integer(1)], "Int64"), None, &family()).unwrap();
    assert_eq!(data(&context), vec![false]);
    assert_eq!(context.mask_dtype, "bool");
}
#[test]
fn pd_na_family_is_enabled_and_cannot_become_python_none() {
    let error = pedigree_context(
        &frame(vec![], "object"),
        Some(&ped()),
        &Family::Value(Cell::Null),
    )
    .err()
    .unwrap();
    assert_eq!(error.error.0, "IndexError");
    assert_eq!(error.source_line, Some(502));
    assert!(error.warning.is_none());
}
#[test]
fn required_lookup_and_no_affected_failure_order_is_unchanged() {
    let input = frame(vec![], "object");
    let mut pedigree = ped();
    pedigree.columns.clear();
    let error = pedigree_context(&input, Some(&pedigree), &family())
        .err()
        .unwrap();
    assert_eq!(error.error.1, "'#FamilyID'");
    assert_eq!(error.source_line, Some(503));
    pedigree.columns = vec![column(
        "#FamilyID",
        "object",
        vec![text("other"), text("other")],
    )];
    assert_eq!(
        pedigree_context(&input, Some(&pedigree), &family())
            .err()
            .unwrap()
            .error
            .1,
        "'Phenotype'"
    );
    pedigree
        .columns
        .push(column("Phenotype", "object", vec![text("1"), text("1")]));
    assert_eq!(
        pedigree_context(&input, Some(&pedigree), &family())
            .err()
            .unwrap()
            .source_line,
        Some(502)
    );
}
#[test]
fn original_unknown_sex_warning_precedes_numeric_accessor_failure() {
    let mut pedigree = ped();
    pedigree.columns[3].values[0] = text("unknown");
    let error = pedigree_context(
        &frame(vec![Cell::Integer(1)], "int64"),
        Some(&pedigree),
        &family(),
    )
    .err()
    .unwrap();
    assert_eq!(error.error.0, "AttributeError");
    assert_eq!(error.source_line, Some(518));
    assert_eq!(error.warning, Some(Warning::UnknownSex));
    assert_eq!(error.pedigree_position, Some(0));
}
#[test]
fn mixed_object_accessor_inference_is_column_wide() {
    let context = pedigree_context(
        &frame(
            vec![Cell::Bool(true), Cell::Integer(1), Cell::Null],
            "object",
        ),
        Some(&ped()),
        &family(),
    )
    .unwrap();
    assert_eq!(data(&context), vec![false, false, false]);
    let error = pedigree_context(
        &frame(vec![Cell::Float(2.5), Cell::Integer(1)], "object"),
        Some(&ped()),
        &family(),
    )
    .err()
    .unwrap();
    assert_eq!(error.error.0, "AttributeError");
    assert_eq!(
        pedigree_context(&frame(vec![], "float64"), Some(&ped()), &family())
            .err()
            .unwrap()
            .source_line,
        Some(518)
    );
}
#[test]
fn sex_predicates_keep_python_strip_numeric_and_bool_distinctions() {
    for (value, male, female) in [
        (text("\u{2003}MaLe\u{a0}"), true, false),
        (text("\u{1c}F\u{1f}"), false, true),
        (Cell::Integer(1), true, false),
        (Cell::Float(2.), false, true),
        (Cell::Bool(true), false, false),
        (Cell::Null, false, false),
        (text("ΟΣ"), false, false),
    ] {
        let mut pedigree = ped();
        pedigree.columns[3].values[0] = value;
        let context = pedigree_context(
            &frame(vec![text("1")], "object"),
            Some(&pedigree),
            &family(),
        )
        .unwrap();
        assert_eq!((context.is_male, context.is_female), (male, female));
        assert_eq!(context.warning.is_some(), !male && !female);
    }
}
#[test]
fn exact_numeric_family_comparison_does_not_round_large_integers() {
    let mut pedigree = ped();
    pedigree.columns[0].values = vec![
        Cell::Float(9007199254740992.),
        Cell::Integer(9007199254740993),
    ];
    pedigree.columns[2].values[1] = text("P");
    let context = pedigree_context(
        &frame(vec![text("1")], "object"),
        Some(&pedigree),
        &Family::Value(Cell::Integer(9007199254740993)),
    )
    .unwrap();
    assert_eq!(context.pedigree_position, Some(1));
    pedigree.columns[0].values = vec![Cell::Float(18446744073709551616.), Cell::Unsigned(u64::MAX)];
    assert_eq!(
        pedigree_context(
            &frame(vec![], "object"),
            Some(&pedigree),
            &Family::Value(Cell::Unsigned(u64::MAX))
        )
        .unwrap()
        .pedigree_position,
        Some(1)
    );
}
#[test]
fn missing_and_duplicate_sample_fail_at_original_accessor_line() {
    let mut input = frame(vec![text("1")], "object");
    input.columns[0].name = "ABSENT".into();
    let error = pedigree_context(&input, Some(&ped()), &family())
        .err()
        .unwrap();
    assert_eq!(error.error.1, "'P'");
    input.columns[0].name = "P".into();
    input.columns.push(input.columns[0].clone());
    let error = pedigree_context(&input, Some(&ped()), &family())
        .err()
        .unwrap();
    assert_eq!(error.error.1, "'DataFrame' object has no attribute 'str'");
    assert_eq!(error.source_line, Some(518));
}
#[test]
fn unsupported_nullable_accessor_and_ped_objects_explicitly_defer() {
    let error = pedigree_context(&frame(vec![], "UInt64"), Some(&ped()), &family())
        .err()
        .unwrap();
    assert_eq!(error.error.0, "UnsupportedInput");
    assert!(error.error.1.contains("cache effects"));
    let mut pedigree = ped();
    pedigree.columns[2].values[0] = Cell::Integer(1);
    assert_eq!(
        pedigree_context(&frame(vec![], "object"), Some(&pedigree), &family())
            .err()
            .unwrap()
            .error
            .0,
        "UnsupportedInput"
    );
}

#[test]
fn native_family_column_numpy_comparison_remains_python() {
    let mut pedigree = ped();
    pedigree.columns[0] = column(
        "#FamilyID",
        "float64",
        vec![Cell::Float(9007199254740992.), Cell::Float(0.)],
    );
    let error = pedigree_context(
        &frame(vec![], "object"),
        Some(&pedigree),
        &Family::Value(Cell::Integer(9007199254740993)),
    )
    .err()
    .unwrap();
    assert_eq!(error.error.0, "UnsupportedInput");
    assert!(error.error.1.contains("NumPy comparison"));
}
