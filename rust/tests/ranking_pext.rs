use priva_tabular::acmg_consequence::coerce_pext_value::Value;
use priva_tabular::acmg_scoring::control_multiply::RowIndex;
use priva_tabular::acmg_scoring::ranking_pext::{
    CoercionColumn, Modulation, Options, ranking_pext, select_column,
};
use priva_tabular::acmg_scoring::terminal_order::{Input, InputColumn};
use priva_tabular::acmg_variant_mechanism::variant_mechanism_masks::Cell;

fn column(name: &str, values: &[f64]) -> InputColumn {
    InputColumn {
        name: name.into(),
        dtype: "float64".into(),
        values: values.iter().copied().map(Cell::Float).collect(),
    }
}
fn frame() -> Input<i64> {
    Input {
        columns: vec![
            column("sort_index", &[0.475, 0.475, -0.0, 0.475]),
            column("PEXT_MEAN", &[0.02, 0.2, 0., 0.02]),
        ],
        index: vec![7, 7, -2, 91],
        axis_metadata: "original duplicate/gapped axis/attrs".into(),
    }
}
fn options() -> Options {
    Options {
        tissues: Cell::Null,
        cutoff: Cell::Float(0.1),
        floor: Cell::Float(0.8),
        shape: Cell::Float(0.5),
    }
}
fn mask(data: &Input<i64>) -> Modulation<i64> {
    Modulation {
        index: data.index.clone(),
        dtype: "bool",
        data: vec![true; 4],
        null: vec![false; 4],
    }
}
fn values(data: &Input<i64>) -> CoercionColumn<i64> {
    CoercionColumn {
        position: 1,
        index: data.index.clone(),
        values: data.columns[1]
            .values
            .iter()
            .cloned()
            .map(Value::Primitive)
            .collect(),
    }
}
fn bits(data: &Input<i64>, name: &str) -> Vec<u64> {
    data.columns
        .iter()
        .find(|c| c.name == name)
        .unwrap()
        .values
        .iter()
        .map(|v| match v {
            Cell::Float(v) => v.to_bits(),
            _ => panic!(),
        })
        .collect()
}

#[test]
fn source_default_operation_order_and_duplicate_axis_are_retained() {
    let mut data = frame();
    let result = ranking_pext(
        &mut data,
        &options(),
        Some(&values(&frame())),
        &mask(&frame()),
        RowIndex::Flat,
    )
    .unwrap();
    assert_eq!(result.index, [7, 7, -2, 91]);
    assert_eq!(result.axis_metadata, "original duplicate/gapped axis/attrs");
    assert_eq!(result.target_column_position, Some(2));
    assert_eq!(bits(&data, "pext_sort_index")[1], 1f64.to_bits());
    assert_eq!(bits(&data, "sort_index")[2], (-0f64).to_bits());
    assert_eq!(
        result
            .values
            .unwrap()
            .column("scaled_pext")
            .unwrap()
            .f64()
            .unwrap()
            .get(0)
            .unwrap()
            .to_bits(),
        (0.02f64 / 0.1).to_bits()
    );
}
#[test]
fn first_missing_tissue_falls_back_to_mean_without_searching_later_tissues() {
    let mut data = frame();
    data.columns.push(column("PEXT_Liver", &[0.01; 4]));
    let result = select_column(&data, &Cell::Text(" Missing , Liver ".into())).unwrap();
    assert_eq!(result.column.as_deref(), Some("PEXT_MEAN"));
    assert_eq!(result.first_tissue_column.as_deref(), Some("PEXT_Missing"));
    assert_eq!(result.tissue_list.unwrap(), ["Missing", "Liver"]);
    assert!(!result.used_tissue);
}
#[test]
fn python_strip_and_hyphen_conversion_preserve_source_selection_locals() {
    let mut data = frame();
    data.columns.push(column("PEXT_Brain_Cortex", &[0.01; 4]));
    let result = select_column(
        &data,
        &Cell::Text("\u{001c}Brain-Cortex\u{001f},Liver".into()),
    )
    .unwrap();
    assert_eq!(result.column.as_deref(), Some("PEXT_Brain_Cortex"));
    assert!(result.used_tissue);
    assert!(
        select_column(&data, &Cell::Text("".into()))
            .unwrap()
            .tissue_list
            .is_none()
    );
    assert_eq!(
        select_column(&data, &Cell::Text(" , , ".into()))
            .unwrap()
            .tissue_list,
        Some(Vec::new())
    );
}
#[test]
fn no_pext_column_skips_invalid_parameters_mask_score_and_axis() {
    let mut data = frame();
    data.columns.remove(1);
    data.columns[0].dtype = "object".into();
    let before = format!("{data:?}");
    let bad = Options {
        tissues: Cell::Null,
        cutoff: Cell::Null,
        floor: Cell::Text("invalid".into()),
        shape: Cell::Float(-1.),
    };
    let mut mask = mask(&data);
    mask.dtype = "object";
    mask.index.clear();
    let result = ranking_pext(&mut data, &bad, None, &mask, RowIndex::MultiIndex).unwrap();
    assert!(result.values.is_none());
    assert_eq!(format!("{data:?}"), before);
}
#[test]
fn nullable_mask_hidden_true_na_bits_are_unselected_and_unchanged() {
    let mut data = frame();
    let mut selected = mask(&data);
    selected.dtype = "boolean";
    selected.null = vec![true, false, true, false];
    let before = selected.data.clone();
    let old = bits(&data, "sort_index");
    ranking_pext(
        &mut data,
        &options(),
        Some(&values(&frame())),
        &selected,
        RowIndex::Flat,
    )
    .unwrap();
    assert_eq!(bits(&data, "sort_index")[0], old[0]);
    assert_eq!(bits(&data, "sort_index")[2], old[2]);
    assert_eq!(selected.data, before);
}
#[test]
fn numpy_bool_origin_is_distinct_from_python_bool_before_fill() {
    let mut data = frame();
    let mut boxed = values(&data);
    boxed.values = vec![
        Value::NumpyBool(false),
        Value::Primitive(Cell::Bool(false)),
        Value::Primitive(Cell::Text("NAN&0.02".into())),
        Value::Primitive(Cell::Text("0.02&NAN".into())),
    ];
    let result = ranking_pext(
        &mut data,
        &options(),
        Some(&boxed),
        &mask(&frame()),
        RowIndex::Flat,
    )
    .unwrap();
    let result = result.values.unwrap();
    let pext = result.column("pext_numeric").unwrap().f64().unwrap();
    assert_eq!(
        pext.into_no_null_iter().collect::<Vec<_>>(),
        [1., 0., 1., 0.02]
    );
}
#[test]
fn in_range_negative_zero_survives_clip_and_division() {
    let mut data = frame();
    let mut boxed = values(&data);
    boxed.values[0] = Value::Primitive(Cell::Text("-0&0".into()));
    let result = ranking_pext(
        &mut data,
        &options(),
        Some(&boxed),
        &mask(&frame()),
        RowIndex::Flat,
    )
    .unwrap();
    let result = result.values.unwrap();
    for column in ["pext_numeric", "scaled_pext", "power"] {
        assert_eq!(
            result
                .column(column)
                .unwrap()
                .f64()
                .unwrap()
                .get(0)
                .unwrap()
                .to_bits(),
            (-0f64).to_bits()
        );
    }
}
#[test]
fn actual_series_half_power_keeps_the_strict_cutoff_neighbor_witness() {
    let mut data = frame();
    let mut boxed = values(&data);
    boxed.values[0] = Value::Primitive(Cell::Float(f64::from_bits(0x3fb9999999999999)));
    let result = ranking_pext(
        &mut data,
        &options(),
        Some(&boxed),
        &mask(&frame()),
        RowIndex::Flat,
    )
    .unwrap();
    let result = result.values.unwrap();
    // Actual Series expression: old generic powf returned 3ff0000000000000.
    // Preserve the original source's lower adjacent float at this boundary.
    assert_eq!(
        result
            .column("scaled_pext")
            .unwrap()
            .f64()
            .unwrap()
            .get(0)
            .unwrap()
            .to_bits(),
        0x3fefffffffffffff
    );
    assert_eq!(
        result
            .column("power")
            .unwrap()
            .f64()
            .unwrap()
            .get(0)
            .unwrap()
            .to_bits(),
        0x3fefffffffffffff
    );
}

#[test]
fn existing_pext_column_is_replaced_at_its_original_position() {
    let mut data = frame();
    data.columns.insert(0, column("pext_sort_index", &[99.; 4]));
    let mut boxed = values(&frame());
    boxed.position = 2;
    let result = ranking_pext(
        &mut data,
        &options(),
        Some(&boxed),
        &mask(&frame()),
        RowIndex::Flat,
    )
    .unwrap();
    assert_eq!(result.target_column_position, Some(0));
    assert_eq!(result.score_column_position, Some(1));
    assert_eq!(data.columns.len(), 3);
}
#[test]
fn parameter_failures_and_nonfinite_routes_reject_atomically() {
    for option in 0..3 {
        for invalid in [0., -1., f64::NAN, f64::INFINITY] {
            let mut data = frame();
            let before = format!("{data:?}");
            let mut config = options();
            match option {
                0 => config.cutoff = Cell::Float(invalid),
                1 => config.floor = Cell::Float(if invalid == 0. { 2. } else { invalid }),
                _ => config.shape = Cell::Float(invalid),
            }
            assert!(
                ranking_pext(
                    &mut data,
                    &config,
                    Some(&values(&frame())),
                    &mask(&frame()),
                    RowIndex::Flat
                )
                .is_err()
            );
            assert_eq!(format!("{data:?}"), before);
        }
    }
}
#[test]
fn unsupported_alignment_storage_or_duplicate_labels_reject_before_assignments() {
    for route in 0..5 {
        let mut data = frame();
        let mut boxed = values(&data);
        let mut selected = mask(&data);
        match route {
            0 => data.columns[0].dtype = "Float64".into(),
            1 => selected.index.reverse(),
            2 => boxed.position = 0,
            3 => data.columns.push(data.columns[0].clone()),
            _ => boxed.values.pop().map(|_| ()).unwrap(),
        }
        let before = format!("{data:?}");
        assert!(
            ranking_pext(
                &mut data,
                &options(),
                Some(&boxed),
                &selected,
                RowIndex::Flat
            )
            .is_err()
        );
        assert_eq!(format!("{data:?}"), before);
    }
}
#[test]
fn scalar_overflow_after_other_rows_still_rejects_the_whole_plan_atomically() {
    let mut data = frame();
    let mut boxed = values(&data);
    boxed.values[3] = Value::BigInteger("9".repeat(400));
    let before = format!("{data:?}");
    assert!(
        ranking_pext(
            &mut data,
            &options(),
            Some(&boxed),
            &mask(&frame()),
            RowIndex::Flat
        )
        .is_err()
    );
    assert_eq!(format!("{data:?}"), before);
}
#[test]
fn empty_selected_column_and_multiindex_require_original_pandas() {
    let mut data = frame();
    let before = format!("{data:?}");
    assert!(
        ranking_pext(
            &mut data,
            &options(),
            Some(&values(&frame())),
            &mask(&frame()),
            RowIndex::MultiIndex
        )
        .is_err()
    );
    assert_eq!(format!("{data:?}"), before);
    data.index.clear();
    for column in &mut data.columns {
        column.values.clear();
    }
    let mask = Modulation {
        index: vec![],
        dtype: "bool",
        data: vec![],
        null: vec![],
    };
    assert!(ranking_pext(&mut data, &options(), None, &mask, RowIndex::Flat).is_err());
}

#[test]
fn current_r4_mask_composes_in_memory_without_replaying_r2_or_other_math() {
    use priva_tabular::acmg_scoring::pedigree_context::{Family, pedigree_context};
    use priva_tabular::acmg_scoring::ranking_compatibility::ranking_compatibility;
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
    ] {
        data.columns.push(InputColumn {
            name: name.into(),
            dtype: dtype.into(),
            values: vec![value; 4],
        });
    }
    let context = pedigree_context(&data, None, &Family::None).unwrap();
    let r4 = ranking_compatibility(&mut data, &context).unwrap();
    let selected = Modulation::from_compatibility(&r4).unwrap();
    assert_eq!(selected.data, [true; 4]);
    let boxed = values(&data);
    let result = ranking_pext(
        &mut data,
        &options(),
        Some(&boxed),
        &selected,
        RowIndex::Flat,
    )
    .unwrap();
    assert_eq!(result.selected, [true; 4]);
    // Without a pedigree, R2 has no heterozygous rows for this distinct R4
    // compatibility penalty. Its high-score PEXT mask above still selects all.
    assert_eq!(r4.log_counts[4], 0);
}

#[test]
fn native_nonfinite_scores_and_signed_nan_payload_are_not_normalized() {
    let mut data = frame();
    let original = [
        f64::from_bits(0xfff8000000001234),
        f64::INFINITY,
        f64::NEG_INFINITY,
        -0.,
    ];
    data.columns[0].values = original.map(Cell::Float).into();
    ranking_pext(
        &mut data,
        &options(),
        Some(&values(&frame())),
        &mask(&frame()),
        RowIndex::Flat,
    )
    .unwrap();
    assert_eq!(bits(&data, "sort_index"), original.map(f64::to_bits));
}
