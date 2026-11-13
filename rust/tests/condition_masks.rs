use polars_core::prelude::*;
use priva_tabular::acmg_variant_mechanism::variant_condition_masks::*;
use priva_tabular::acmg_variant_mechanism::variant_mechanism_masks::{
    Cell, Dtype, Input, InputColumn,
};

fn column(name: &str, dtype: Dtype, values: Vec<Cell>) -> InputColumn {
    InputColumn {
        name: name.into(),
        dtype,
        values,
    }
}
fn text(name: &str, values: &[&str]) -> InputColumn {
    column(
        name,
        Dtype::Object,
        values.iter().map(|s| Cell::Text((*s).into())).collect(),
    )
}
fn input(columns: Vec<InputColumn>, len: usize) -> Input {
    Input {
        columns,
        index: (0..len).map(|i| Cell::Integer(i as i64)).collect(),
        index_names: vec![Some("labels".into())],
    }
}
fn values(
    masks: &priva_tabular::acmg_variant_mechanism::variant_mechanism_masks::Masks,
    name: &str,
) -> Vec<bool> {
    masks
        .frame
        .column(name)
        .unwrap()
        .bool()
        .unwrap()
        .into_no_null_iter()
        .collect()
}

#[test]
fn missing_columns_are_sorted_and_precede_conversion() {
    let empty = input(vec![], 1);
    let error = variant_condition_masks(&empty).err().unwrap();
    assert_eq!(error.0, "KeyError");
    assert_eq!(
        error.1,
        "'variant-level condition annotations are required; missing columns: variant_inheritance, variant_penetrance'"
    );
    let missing = input(
        vec![column(
            "variant_penetrance",
            Dtype::Nullable("boolean"),
            vec![Cell::Null],
        )],
        1,
    );
    assert_eq!(
        variant_condition_masks(&missing).err().unwrap().1,
        "'variant-level condition annotations are required; missing columns: variant_inheritance'"
    );
}

#[test]
fn inheritance_conversion_precedes_penetrance_despite_column_order() {
    let df = input(
        vec![
            column(
                "variant_penetrance",
                Dtype::Nullable("boolean"),
                vec![Cell::Null],
            ),
            column(
                "variant_inheritance",
                Dtype::Nullable("Int64"),
                vec![Cell::Null],
            ),
        ],
        1,
    );
    let error = variant_condition_masks(&df).err().unwrap();
    assert_eq!(error.0, "TypeError");
    assert_eq!(error.1, "Invalid value '' for dtype Int64");
}

#[test]
fn no_inheritance_tags_still_converts_penetrance() {
    let df = input(
        vec![
            text("variant_inheritance", &[""]),
            column(
                "variant_penetrance",
                Dtype::Nullable("Float64"),
                vec![Cell::Null],
            ),
        ],
        1,
    );
    assert_eq!(
        variant_condition_masks(&df).err().unwrap().1,
        "Invalid value '' for dtype Float64"
    );
}

#[test]
fn duplicate_inheritance_fails_before_later_penetrance_conversion() {
    let df = input(
        vec![
            text("variant_inheritance", &[""]),
            text("variant_inheritance", &[""]),
            column(
                "variant_penetrance",
                Dtype::Nullable("boolean"),
                vec![Cell::Null],
            ),
        ],
        1,
    );
    let error = variant_condition_masks(&df).err().unwrap();
    assert_eq!(error.0, "AttributeError");
    assert_eq!(error.1, "'DataFrame' object has no attribute 'str'");
}

#[test]
fn duplicate_selected_columns_fill_before_dataframe_str_failure() {
    for name in ["variant_inheritance", "variant_penetrance"] {
        let other = if name == "variant_inheritance" {
            "variant_penetrance"
        } else {
            "variant_inheritance"
        };
        let df = input(
            vec![
                text(name, &[""]),
                column(name, Dtype::Nullable("Int64"), vec![Cell::Null]),
                text(other, &[""]),
            ],
            1,
        );
        assert_eq!(
            variant_condition_masks(&df).err().unwrap().1,
            "Invalid value '' for dtype Int64"
        );
    }
}

#[test]
fn thirteen_masks_keep_source_order_names_and_independent_histories() {
    let df = input(
        vec![
            text(
                "variant_inheritance",
                &[
                    "recessive;dominant;x_linked_recessive;x_linked_dominant;x_linked_unspecified;y_linked;mitochondrial;polygenic;digenic;oligogenic;non_mendelian",
                    "x_linked_unspecified",
                    "non_mendelian",
                    "digenic",
                    "unresolved",
                ],
            ),
            text(
                "variant_penetrance",
                &["complete;incomplete", "complete", "incomplete", "high", ""],
            ),
        ],
        5,
    );
    let out = variant_condition_masks(&df).unwrap();
    assert_eq!(out.frame.get_column_names_str(), NAMES);
    for c in out.frame.get_columns() {
        assert_eq!(c.dtype(), &DataType::Boolean);
        assert_eq!(c.null_count(), 0);
    }
    for (name, expected) in [
        ("has_recessive", vec![true, false, false, false, false]),
        ("has_dominant_like", vec![true, false, false, false, false]),
        ("has_mendelian", vec![true, true, false, false, false]),
        ("has_non_monogenic", vec![true, false, false, true, false]),
        ("has_non_mendelian", vec![true, false, true, false, false]),
        (
            "has_incomplete_penetrance",
            vec![true, false, true, false, false],
        ),
        (
            "has_complete_penetrance",
            vec![true, true, false, false, false],
        ),
    ] {
        assert_eq!(values(&out, name), expected);
    }
    assert_eq!(
        out.source_series_names,
        [
            vec![Some("variant_inheritance"); 11],
            vec![Some("variant_penetrance"); 2]
        ]
        .concat()
    );
}

#[test]
fn python_final_lf_and_exact_semicolon_boundaries() {
    let texts = [
        "recessive\n",
        "recessive\n\n",
        "x;recessive\n",
        "recessive\n;",
        "recessive;\n",
        "recessive\r\n",
        " recessive",
        "recessive ",
        "recessive;recessive\njunk",
        ";RECESSIVE;",
    ];
    let df = input(
        vec![
            text("variant_inheritance", &texts),
            text("variant_penetrance", &["incomplete\n"; 10]),
        ],
        10,
    );
    let out = variant_condition_masks(&df).unwrap();
    assert_eq!(
        values(&out, "has_recessive"),
        vec![
            true, false, true, false, true, false, false, false, true, true
        ]
    );
    assert_eq!(values(&out, "has_incomplete_penetrance"), vec![true; 10]);
    assert_eq!(values(&out, "has_complete_penetrance"), vec![false; 10]);
}

#[test]
fn lowercase_is_unicode_aware_without_trimming_or_substring_matching() {
    let df = input(
        vec![
            text(
                "variant_inheritance",
                &[
                    "X_LIN\u{212a}ED_DOMINANT",
                    "DOMINANTLY_INHERITED",
                    "\u{a0}recessive",
                    "M\u{130}TOCHONDRIAL",
                ],
            ),
            text(
                "variant_penetrance",
                &[
                    "COMPLETE",
                    "INCOMPLETE",
                    "incomplete_like",
                    "incomplete;complete",
                ],
            ),
        ],
        4,
    );
    let out = variant_condition_masks(&df).unwrap();
    assert_eq!(
        values(&out, "has_x_linked_dominant"),
        vec![true, false, false, false]
    );
    assert_eq!(values(&out, "has_dominant"), vec![false; 4]);
    assert_eq!(values(&out, "has_mitochondrial"), vec![false; 4]);
    assert_eq!(
        values(&out, "has_complete_penetrance"),
        vec![true, false, false, true]
    );
}

#[test]
fn object_nulls_nan_and_primitive_scalars_have_no_tokens() {
    let scalars = vec![
        Cell::Null,
        Cell::Float(f64::NAN),
        Cell::Float(f64::INFINITY),
        Cell::Bool(true),
        Cell::Integer(-1),
        Cell::Unsigned(u64::MAX),
    ];
    let df = input(
        vec![
            column("variant_inheritance", Dtype::Object, scalars.clone()),
            column("variant_penetrance", Dtype::Object, scalars),
        ],
        6,
    );
    let out = variant_condition_masks(&df).unwrap();
    assert!(
        out.frame
            .get_columns()
            .iter()
            .all(|c| c.bool().unwrap().into_no_null_iter().all(|b| !b))
    );
}

#[test]
fn nullable_dtype_missingness_is_preserved() {
    for dtype in ["Int64", "Float64", "boolean"] {
        let df = input(
            vec![
                column(
                    "variant_inheritance",
                    Dtype::Nullable(dtype),
                    vec![Cell::Null],
                ),
                text("variant_penetrance", &[""]),
            ],
            1,
        );
        let error = variant_condition_masks(&df).err().unwrap();
        assert_eq!(error.0, "TypeError");
        assert_eq!(error.1, format!("Invalid value '' for dtype {dtype}"));
    }
    let df = input(
        vec![
            column("variant_inheritance", Dtype::String, vec![Cell::Null]),
            column(
                "variant_penetrance",
                Dtype::Numpy,
                vec![Cell::Float(f64::NAN)],
            ),
        ],
        1,
    );
    assert!(!values(&variant_condition_masks(&df).unwrap(), "has_mendelian")[0]);
}

#[test]
fn generic_tuple_indices_names_duplicates_and_input_are_retained() {
    let df = Input {
        columns: vec![
            text(
                "variant_inheritance",
                &["dominant", "recessive", "dominant"],
            ),
            text("variant_penetrance", &["complete", "incomplete", ""]),
        ],
        index: vec![("chr1", 9), ("chr2", -3), ("chr1", 9)],
        index_names: vec![Some("chrom".into()), Some("position".into())],
    };
    let before = format!("{df:?}");
    let out = variant_condition_masks(&df).unwrap();
    assert_eq!(out.index, df.index);
    assert_eq!(out.index_names, df.index_names);
    assert_eq!(before, format!("{df:?}"));
}

#[test]
fn empty_masks_still_have_all_boolean_columns() {
    let df = input(
        vec![
            text("variant_inheritance", &[]),
            text("variant_penetrance", &[]),
        ],
        0,
    );
    let out = variant_condition_masks(&df).unwrap();
    assert_eq!(out.frame.shape(), (0, 13));
    assert_eq!(out.source_series_names.len(), 13);
    assert!(
        out.frame
            .get_columns()
            .iter()
            .all(|c| c.dtype() == &DataType::Boolean)
    );
}

#[test]
fn unrelated_duplicate_columns_are_not_evaluated() {
    let df = input(
        vec![
            text("variant_inheritance", &["dominant"]),
            text("variant_penetrance", &["complete"]),
            column("ignored", Dtype::Nullable("category"), vec![Cell::Null]),
            column("ignored", Dtype::Nullable("category"), vec![Cell::Null]),
        ],
        1,
    );
    assert_eq!(
        values(&variant_condition_masks(&df).unwrap(), "has_dominant"),
        vec![true]
    );
}

#[test]
fn polars_adapter_preserves_input_and_does_not_join_duplicate_indices() {
    let df = DataFrame::new(vec![
        Column::new("variant_inheritance".into(), [Some("dominant"), None]),
        Column::new("variant_penetrance".into(), ["complete", "incomplete"]),
    ])
    .unwrap();
    let before = df.clone();
    let out = from_polars(
        &df,
        vec![9_i64, 9],
        vec![Some("labels".into())],
        [Dtype::Object, Dtype::String],
    )
    .unwrap();
    assert!(df.equals_missing(&before));
    assert_eq!(out.index, vec![9_i64, 9]);
    assert_eq!(out.source_series_names[12], Some("variant_penetrance"));
    assert_eq!(
        out.frame
            .column("has_complete_penetrance")
            .unwrap()
            .bool()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<_>>(),
        vec![true, false]
    );
}

#[test]
fn polars_required_and_earlier_fill_errors_precede_unsupported_dtype() {
    let df = DataFrame::new(vec![Column::new(
        "variant_inheritance".into(),
        [b"unsupported".as_slice()],
    )])
    .unwrap();
    assert_eq!(
        from_polars(&df, vec![1_i64], vec![None], [Dtype::Object, Dtype::Object])
            .err()
            .unwrap()
            .0,
        "KeyError"
    );
    let df = DataFrame::new(vec![
        Column::new("variant_inheritance".into(), [None::<i64>]),
        Column::new("variant_penetrance".into(), [b"unsupported".as_slice()]),
    ])
    .unwrap();
    let error = from_polars(
        &df,
        vec![1_i64],
        vec![None],
        [Dtype::Nullable("Int64"), Dtype::Object],
    )
    .err()
    .unwrap();
    assert_eq!(error.0, "TypeError");
    assert_eq!(error.1, "Invalid value '' for dtype Int64");
}

#[test]
fn empty_and_all_null_unsupported_polars_dtypes_require_python() {
    for len in [0, 1] {
        let df = DataFrame::new(vec![
            Column::new("variant_inheritance".into(), vec![""; len]),
            Series::full_null("variant_penetrance".into(), len, &DataType::Binary).into(),
        ])
        .unwrap();
        let error = from_polars(
            &df,
            vec![1_i64; len],
            vec![None],
            [Dtype::Object, Dtype::Object],
        )
        .err()
        .unwrap();
        assert_eq!(error.0, "UnfinishedParity");
        assert!(error.1.contains("variant_penetrance"));
    }
}
