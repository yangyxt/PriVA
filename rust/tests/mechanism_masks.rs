use polars_core::prelude::*;
use priva_tabular::acmg_variant_mechanism::variant_mechanism_masks::*;
fn col(name: &str, values: Vec<Cell>) -> InputColumn {
    InputColumn {
        name: name.into(),
        dtype: Dtype::Object,
        values,
    }
}
#[test]
fn required_errors_precede_nullable_error() {
    let input = Input {
        index_names: vec![Some("row_label".into())],
        index: vec![Cell::Integer(9)],
        columns: vec![InputColumn {
            name: "variant_effect".into(),
            dtype: Dtype::Nullable("Int64"),
            values: vec![Cell::Null],
        }],
    };
    let err = variant_mechanism_masks(&input).err().unwrap();
    assert_eq!(err.0, "KeyError");
    assert!(err.1.contains("var_plausible_patho_mechs"));
}
#[test]
fn profile_error_precedes_effect_and_scores() {
    let input = Input {
        index_names: vec![Some("row_label".into())],
        index: vec![Cell::Integer(9)],
        columns: vec![
            InputColumn {
                name: "var_plausible_patho_mechs".into(),
                dtype: Dtype::Nullable("Int64"),
                values: vec![Cell::Null],
            },
            InputColumn {
                name: "variant_effect".into(),
                dtype: Dtype::Nullable("boolean"),
                values: vec![Cell::Null],
            },
        ],
    };
    assert_eq!(
        variant_mechanism_masks(&input).err().unwrap().1,
        "Invalid value '' for dtype Int64"
    );
}
#[test]
fn preserves_gaps_duplicates_and_contradictory_masks() {
    let labels = vec![Cell::Integer(9), Cell::Integer(-3), Cell::Integer(9)];
    let input = Input {
        index_names: vec![Some("row_label".into())],
        index: labels.clone(),
        columns: vec![
            col(
                "var_plausible_patho_mechs",
                vec![
                    Cell::Text("dominant_unknown".into()),
                    Cell::Text("\u{1c}".into()),
                    Cell::Text("recessive_LOF\n".into()),
                ],
            ),
            col(
                "variant_effect",
                vec![
                    Cell::Text("x+GOF".into()),
                    Cell::Text("uncertain".into()),
                    Cell::Text("".into()),
                ],
            ),
            col("variant_lof_score", vec![Cell::Integer(2); 3]),
        ],
    };
    let before = format!("{input:?}");
    let result = variant_mechanism_masks(&input).unwrap();
    assert_eq!(before, format!("{input:?}"));
    assert_eq!(result.index, labels);
    assert_eq!(result.frame.get_column_names_str(), NAMES);
    for (name, expected) in [
        ("has_dominant_compatible", vec![true, false, false]),
        ("has_dom_lof_history", vec![false, false, false]),
        ("is_exact_gof", vec![true, false, false]),
        ("is_predicted_lof", vec![true, false, true]),
        ("modern_profile", vec![true, false, true]),
        ("has_rec_lof_history", vec![false, false, true]),
    ] {
        assert_eq!(
            result
                .frame
                .column(name)
                .unwrap()
                .bool()
                .unwrap()
                .into_no_null_iter()
                .collect::<Vec<_>>(),
            expected
        );
    }
}
#[test]
fn polars_input_is_not_mutated() {
    let df = DataFrame::new(vec![
        Column::new("var_plausible_patho_mechs".into(), [Some("dominant"), None]),
        Column::new("variant_effect".into(), ["uncertain", ""]),
        Column::new("variant_gof_score".into(), [2.0, f64::NAN]),
    ])
    .unwrap();
    let before = df.clone();
    let result = from_polars(
        &df,
        vec![Cell::Integer(5); 2],
        vec![None],
        [
            Dtype::Object,
            Dtype::Object,
            Dtype::Numpy,
            Dtype::Numpy,
            Dtype::Numpy,
        ],
    )
    .unwrap();
    assert!(df.equals_missing(&before));
    assert_eq!(result.frame.height(), 2);
}
#[test]
fn numeric_threshold_and_newline_semantics() {
    assert_eq!(numeric_text("0.99999999999999999"), 1.0);
    assert_eq!(numeric_text("\u{a0}2"), 0.0);
    assert_eq!(numeric_text("2\0junk"), 0.0);
    assert_eq!(numeric_text("2.0\0junk"), 2.0);
    assert_eq!(numeric_text("inf"), f64::INFINITY);
    assert_eq!(numeric_text(" inf "), 0.0);
}

#[test]
fn tuple_labels_and_series_names_are_preserved() {
    let input = Input {
        index: vec![("chr1", 7), ("chr1", 7)],
        index_names: vec![Some("chrom".into()), Some("offset".into())],
        columns: vec![
            col("var_plausible_patho_mechs", vec![Cell::Null; 2]),
            col("variant_effect", vec![Cell::Text("uncertain".into()); 2]),
        ],
    };
    let out = variant_mechanism_masks(&input).unwrap();
    assert_eq!(out.index, input.index);
    assert_eq!(out.index_names, input.index_names);
    assert_eq!(
        out.source_series_names[0],
        Some("var_plausible_patho_mechs")
    );
    assert_eq!(out.source_series_names[1], None);
    assert_eq!(out.source_series_names[5], Some("variant_effect"));
}

#[test]
fn polars_nullable_phase_precedes_later_unsupported_binary() {
    for failing in 0..2 {
        let mut columns = vec![
            Column::new("var_plausible_patho_mechs".into(), [Some(1_i64)]),
            Column::new("variant_effect".into(), [Some(1_i64)]),
            Column::new("variant_lof_score".into(), [b"unsupported".as_slice()]),
        ];
        columns[failing] = Column::new(
            ["var_plausible_patho_mechs", "variant_effect"][failing].into(),
            [None::<i64>],
        );
        let df = DataFrame::new(columns).unwrap();
        let error = from_polars(
            &df,
            vec![7_i64],
            vec![None],
            [
                Dtype::Nullable("Int64"),
                Dtype::Nullable("Int64"),
                Dtype::Object,
                Dtype::Object,
                Dtype::Object,
            ],
        )
        .err()
        .unwrap();
        assert_eq!(error.0, "TypeError");
        assert_eq!(error.1, "Invalid value '' for dtype Int64");
    }
}

#[test]
fn nullable_score_failure_precedes_later_unsupported_column() {
    for failing in 2..4 {
        let mut columns = vec![
            Column::new("var_plausible_patho_mechs".into(), [""]),
            Column::new("variant_effect".into(), [""]),
            Column::new("variant_lof_score".into(), [Some(true)]),
            Column::new("variant_gof_score".into(), [Some(true)]),
            Column::new("variant_dn_score".into(), [Some(true)]),
        ];
        let names = [
            "var_plausible_patho_mechs",
            "variant_effect",
            "variant_lof_score",
            "variant_gof_score",
            "variant_dn_score",
        ];
        columns[failing] = Column::new(names[failing].into(), [None::<bool>]);
        columns[failing + 1] = Column::new(names[failing + 1].into(), [b"unsupported".as_slice()]);
        let df = DataFrame::new(columns).unwrap();
        let error = from_polars(
            &df,
            vec![1_i64],
            vec![None],
            [
                Dtype::Object,
                Dtype::Object,
                Dtype::Nullable("boolean"),
                Dtype::Nullable("boolean"),
                Dtype::Nullable("boolean"),
            ],
        )
        .err()
        .unwrap();
        assert_eq!(error.0, "TypeError");
        assert_eq!(error.1, "Invalid value '0' for dtype boolean");
    }
}

#[test]
fn empty_and_all_null_unsupported_dtype_are_explicitly_unfinished() {
    for len in [0_usize, 1] {
        let df = DataFrame::new(vec![
            Column::new("var_plausible_patho_mechs".into(), vec![""; len]),
            Column::new("variant_effect".into(), vec![""; len]),
            Series::full_null("variant_lof_score".into(), len, &DataType::Binary).into(),
        ])
        .unwrap();
        let error = from_polars(
            &df,
            vec![1_i64; len],
            vec![None],
            [
                Dtype::Object,
                Dtype::Object,
                Dtype::Object,
                Dtype::Object,
                Dtype::Object,
            ],
        )
        .err()
        .unwrap();
        assert_eq!(error.0, "UnfinishedParity");
        assert!(error.1.contains("variant_lof_score"));
    }
}
