use priva_tabular::acmg_scoring::control_multiply::RowIndex;
use priva_tabular::acmg_scoring::pedigree_context::{Family, pedigree_context};
use priva_tabular::acmg_scoring::ranking_compatibility::ranking_compatibility;
use priva_tabular::acmg_scoring::ranking_gene_weights::{
    COMPATIBILITY, DISPENSABLE, GeneList, Penalty, RELEVANT, parse_gene_text, ranking_gene_weights,
};
use priva_tabular::acmg_scoring::terminal_order::{Input, InputColumn};
use priva_tabular::acmg_variant_mechanism::variant_mechanism_masks::Cell;

fn column(name: &str, dtype: &str, values: Vec<Cell>) -> InputColumn {
    InputColumn {
        name: name.into(),
        dtype: dtype.into(),
        values,
    }
}
fn frame() -> Input<i64> {
    Input {
        columns: vec![
            column(
                "SYMBOL",
                "object",
                ["A", "B", "C", "A"].map(|s| Cell::Text(s.into())).into(),
            ),
            column(
                "sort_index",
                "float64",
                [0.95, 0.475, -0., 0.12345678901234568]
                    .map(Cell::Float)
                    .into(),
            ),
            column(
                COMPATIBILITY,
                "float64",
                [0.9, 1., 0.9, 0.9].map(Cell::Float).into(),
            ),
            column(
                "untouched",
                "object",
                vec![Cell::Text("original".into()); 4],
            ),
        ],
        index: vec![7, 7, -2, 90],
        axis_metadata: "caller-owned opaque axis/attrs/flags".into(),
    }
}
fn penalty(frame: &Input<i64>) -> Penalty<i64> {
    Penalty {
        index: frame.index.clone(),
        dtype: "bool",
        data: vec![true, false, true, false],
        null: vec![false; 4],
    }
}
fn bits(frame: &Input<i64>, name: &str) -> Vec<u64> {
    frame
        .columns
        .iter()
        .find(|c| c.name == name)
        .unwrap()
        .values
        .iter()
        .map(|c| match c {
            Cell::Float(x) => x.to_bits(),
            _ => panic!("native float64 expected"),
        })
        .collect()
}

#[test]
fn actual_python_overlap_witness_preserves_three_separate_multiplications() {
    let mut data = frame();
    let mask = penalty(&data);
    let result = ranking_gene_weights(
        &mut data,
        &mask,
        GeneList::PlainText("A\nB\n"),
        GeneList::PlainText("A\nB\n"),
        RowIndex::Flat,
    )
    .unwrap();
    assert_eq!(
        result
            .steps
            .iter()
            .map(|s| s.source_line)
            .collect::<Vec<_>>(),
        [604, 611, 613]
    );
    // Captured actual unchanged Python extra_same_file, manifest a5851667… .
    assert_eq!(
        bits(&data, "sort_index"),
        [
            0x3fed8c7e28240b78,
            0x3fe06a7ef9db22d1,
            0x8000000000000000,
            0x3fc111110e7d5c14
        ]
    );
    assert_eq!(result.dispensable_genes.unwrap(), ["A", "B"]);
    assert_eq!(result.index, [7, 7, -2, 90]);
}

#[test]
fn none_disables_but_empty_file_resets_existing_index_column() {
    let mut data = frame();
    data.columns.insert(
        0,
        column(DISPENSABLE, "string", vec![Cell::Text("old".into()); 4]),
    );
    let mask = penalty(&data);
    let disabled = ranking_gene_weights(
        &mut data,
        &mask,
        GeneList::Disabled,
        GeneList::Disabled,
        RowIndex::Flat,
    )
    .unwrap();
    assert!(disabled.dispensable_genes.is_none());
    assert_eq!(data.columns[0].dtype, "string");
    let enabled = ranking_gene_weights(
        &mut data,
        &mask,
        GeneList::PlainText(""),
        GeneList::Disabled,
        RowIndex::Flat,
    )
    .unwrap();
    assert_eq!(enabled.dispensable_genes, Some(Vec::new()));
    assert_eq!(enabled.steps[0].multiplier_column_position, 0);
    assert_eq!(data.columns[0].dtype, "float64");
    assert_eq!(bits(&data, DISPENSABLE), [1f64.to_bits(); 4]);
}

#[test]
fn existing_target_position_unrelated_columns_and_axes_are_retained() {
    let mut data = frame();
    data.columns
        .insert(1, column(RELEVANT, "int64", vec![Cell::Integer(7); 4]));
    data.columns.push(data.columns[4].clone());
    let before = data
        .columns
        .iter()
        .map(|c| c.name.clone())
        .collect::<Vec<_>>();
    let original_untouched = format!("{:?}", data.columns[4]);
    let axis = data.axis_metadata.clone();
    let mask = penalty(&data);
    let result = ranking_gene_weights(
        &mut data,
        &mask,
        GeneList::Disabled,
        GeneList::PlainText("B"),
        RowIndex::Flat,
    )
    .unwrap();
    assert_eq!(result.steps[0].multiplier_column_position, 1);
    assert_eq!(
        data.columns
            .iter()
            .map(|c| c.name.clone())
            .collect::<Vec<_>>(),
        before
    );
    assert_eq!(format!("{:?}", data.columns[4]), original_untouched);
    assert_eq!(data.axis_metadata, axis);
    assert_eq!(result.axis_metadata, axis);
}

#[test]
fn nullable_mask_does_not_select_na_or_change_hidden_data() {
    let mut data = frame();
    let mask = Penalty {
        index: data.index.clone(),
        dtype: "boolean",
        data: vec![true; 4],
        null: vec![true, false, true, false],
    };
    ranking_gene_weights(
        &mut data,
        &mask,
        GeneList::Disabled,
        GeneList::Disabled,
        RowIndex::Flat,
    )
    .unwrap();
    assert_eq!(mask.data, [true; 4]);
    assert_eq!(mask.null, [true, false, true, false]);
    assert_eq!(bits(&data, "sort_index")[0], 0.95f64.to_bits());
    assert_eq!(bits(&data, "sort_index")[2], (-0.0f64).to_bits());
}

#[test]
fn literal_parser_retains_case_duplicates_comments_and_last_line() {
    assert_eq!(
        parse_gene_text("# header\nA\na\nA# comment\nB-2_3.4").unwrap(),
        ["A", "a", "A", "B-2_3.4"]
    );
    assert!(parse_gene_text("# only\n\n").unwrap().is_empty());
}

#[test]
fn inference_sensitive_parser_routes_are_explicit_python_boundaries() {
    for text in [
        "NA",
        "NaN",
        "nAn",
        "TRUE",
        "fAlSe",
        "iNf",
        "1",
        "01",
        "A\tB",
        " A",
        "A ",
        "\"A\"",
        "\u{feff}A",
        "A\r\n",
        "β",
        "A\0",
    ] {
        assert_eq!(
            parse_gene_text(text).unwrap_err().0,
            "UnsupportedInput",
            "{text:?}"
        );
    }
}

#[test]
fn second_parser_rejection_does_not_apply_first_multiplier() {
    let mut data = frame();
    let before = format!("{data:?}");
    let mask = penalty(&data);
    let error = ranking_gene_weights(
        &mut data,
        &mask,
        GeneList::PlainText("A"),
        GeneList::PlainText("NA"),
        RowIndex::Flat,
    )
    .err()
    .unwrap();
    assert_eq!(error.0, "UnsupportedInput");
    assert_eq!(format!("{data:?}"), before);
}

#[test]
fn two_nan_dtype_axis_and_missing_routes_reject_before_mutation() {
    for mode in 0..5 {
        let mut data = frame();
        let mut mask = penalty(&data);
        let row_index = if mode == 3 {
            RowIndex::MultiIndex
        } else {
            RowIndex::Flat
        };
        match mode {
            0 => {
                data.columns[1].values[0] = Cell::Float(f64::NAN);
                data.columns[2].values[0] = Cell::Float(f64::NAN);
            }
            1 => data.columns[1].dtype = "Float64".into(),
            2 => mask.index[0] = 99,
            4 => {
                data.columns.remove(2);
            }
            _ => (),
        }
        let before = format!("{data:?}");
        let error = ranking_gene_weights(
            &mut data,
            &mask,
            GeneList::PlainText("A"),
            GeneList::Disabled,
            row_index,
        )
        .err()
        .unwrap();
        assert_eq!(error.0, "UnsupportedInput");
        assert_eq!(format!("{data:?}"), before);
    }
}

#[test]
fn nonfinite_scores_and_signed_zero_are_not_normalized() {
    let mut data = frame();
    data.columns[1].values = [
        f64::INFINITY,
        f64::NEG_INFINITY,
        -0.,
        f64::from_bits(0xfff8000000001234),
    ]
    .map(Cell::Float)
    .into();
    let mask = penalty(&data);
    ranking_gene_weights(
        &mut data,
        &mask,
        GeneList::PlainText("A\nB"),
        GeneList::PlainText("A"),
        RowIndex::Flat,
    )
    .unwrap();
    assert_eq!(
        bits(&data, "sort_index"),
        [
            f64::INFINITY.to_bits(),
            f64::NEG_INFINITY.to_bits(),
            (-0f64).to_bits(),
            0xfff8000000001234
        ]
    );
}

#[test]
fn empty_frame_keeps_empty_axis_and_all_enabled_step_positions() {
    let mut data = frame();
    data.index.clear();
    for column in &mut data.columns {
        column.values.clear();
    }
    let mask = Penalty {
        index: Vec::new(),
        dtype: "object",
        data: Vec::new(),
        null: Vec::new(),
    };
    let result = ranking_gene_weights(
        &mut data,
        &mask,
        GeneList::PlainText("A"),
        GeneList::PlainText("B"),
        RowIndex::Flat,
    )
    .unwrap();
    assert_eq!(result.steps.len(), 3);
    assert!(result.steps.iter().all(|step| step.frame.height() == 0));
    assert!(data.index.is_empty());
}

#[test]
fn compose_r2_r4_r5_r6_in_memory_without_replaying_warning_or_prior_math() {
    let mut data = frame();
    for (name, dtype, value) in [
        ("var_plausible_patho_mechs", "object", Cell::Text("".into())),
        (
            "variant_effect",
            "object",
            Cell::Text("predicted_LOF_high_confidence".into()),
        ),
        ("variant_lof_score", "int64", Cell::Integer(1)),
        ("ACMG_quant_score", "float64", Cell::Float(0.95)),
        ("PROBAND", "string", Cell::Text("0/1".into())),
    ] {
        data.columns.push(column(name, dtype, vec![value; 4]));
    }
    let pedigree = Input {
        columns: vec![
            column("#FamilyID", "object", vec![Cell::Text("F".into())]),
            column("IndividualID", "object", vec![Cell::Text("PROBAND".into())]),
            column("Phenotype", "object", vec![Cell::Text("2".into())]),
            column("Sex", "object", vec![Cell::Text("unknown".into())]),
        ],
        index: vec![0],
        axis_metadata: String::new(),
    };
    let context = pedigree_context(
        &data,
        Some(&pedigree),
        &Family::Value(Cell::Text("F".into())),
    )
    .unwrap();
    let warning = context.warning;
    assert!(warning.is_some());
    let r4 = ranking_compatibility(&mut data, &context).unwrap();
    let mask = Penalty::from_compatibility(&r4).unwrap();
    let result = ranking_gene_weights(
        &mut data,
        &mask,
        GeneList::PlainText("A"),
        GeneList::PlainText("A"),
        RowIndex::Flat,
    )
    .unwrap();
    assert_eq!(context.warning, warning);
    assert_eq!(r4.log_counts[4], 4);
    assert_eq!(result.steps[2].source_line, 613);
    assert_eq!(bits(&data, "sort_index")[0], 0x3fed8c7e28240b78);
}
