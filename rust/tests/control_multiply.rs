use priva_tabular::acmg_scoring::control_multiply::{RowIndex, control_multiply};
use priva_tabular::acmg_scoring::terminal_order::{Input, InputColumn};
use priva_tabular::acmg_variant_mechanism::variant_mechanism_masks::Cell;

fn column(name: &str, values: &[f64]) -> InputColumn {
    InputColumn {
        name: name.into(),
        dtype: "float64".into(),
        values: values.iter().copied().map(Cell::Float).collect(),
    }
}
fn frame(left: &[f64], right: &[f64]) -> Input<i32> {
    Input {
        index: vec![7; left.len()],
        columns: vec![
            column("untouched", left),
            column("sort_index", left),
            column("control_common_index", right),
            column("untouched", right),
        ],
        axis_metadata: "opaque duplicate/gapped axes and original metadata".into(),
    }
}
fn bits(input: &Input<i32>, position: usize) -> Vec<u64> {
    input.columns[position]
        .values
        .iter()
        .map(|cell| match cell {
            Cell::Float(value) => value.to_bits(),
            _ => panic!("expected source float64"),
        })
        .collect()
}

#[test]
fn mutates_existing_score_position_and_preserves_all_other_buffers() {
    let mut input = frame(&[0.95, 0.5, 1.], &[0.5, 1., 0.9]);
    let before = input.clone();
    assert_eq!(control_multiply(&mut input, RowIndex::Flat).unwrap(), 1);
    assert_eq!(
        bits(&input, 1),
        vec![0.475f64.to_bits(), 0.5f64.to_bits(), 0.9f64.to_bits()]
    );
    assert_eq!(input.index, before.index);
    assert_eq!(input.axis_metadata, before.axis_metadata);
    assert_eq!(input.columns.len(), 4);
    for position in [0, 2, 3] {
        assert_eq!(bits(&input, position), bits(&before, position));
    }
    for (actual, original) in input.columns.iter().zip(&before.columns) {
        assert_eq!(actual.name, original.name);
        assert_eq!(actual.dtype, original.dtype);
    }
}

#[test]
fn preserves_negative_zero_and_gradual_underflow() {
    let mut input = frame(
        &[0., -0., -0., f64::from_bits(1), -f64::from_bits(1)],
        &[0.5, 0.5, -1., 0.5, 0.5],
    );
    control_multiply(&mut input, RowIndex::Flat).unwrap();
    assert_eq!(
        bits(&input, 1),
        vec![0, (-0.0f64).to_bits(), 0, 0, (-0.0f64).to_bits()]
    );
}

#[test]
fn retains_nonfinite_source_results_instead_of_clipping_or_null_conversion() {
    let mut input = frame(
        &[f64::MAX, f64::INFINITY, f64::NEG_INFINITY, f64::NAN],
        &[2., 0., 1., 0.5],
    );
    control_multiply(&mut input, RowIndex::Flat).unwrap();
    let values = bits(&input, 1)
        .into_iter()
        .map(f64::from_bits)
        .collect::<Vec<_>>();
    assert_eq!(values[0], f64::INFINITY);
    assert!(values[1].is_nan());
    assert_eq!(values[2], f64::NEG_INFINITY);
    assert!(values[3].is_nan());
    assert!(
        input.columns[1]
            .values
            .iter()
            .all(|cell| matches!(cell, Cell::Float(_)))
    );
}

#[test]
fn two_nan_payload_route_defers_before_any_finite_or_signalling_value_changes() {
    let mut input = frame(
        &[0.95, f64::from_bits(0x7ff0000000000789)],
        &[0.5, f64::from_bits(0xfff8000000000456)],
    );
    let before = input.clone();
    let error = control_multiply(&mut input, RowIndex::Flat).unwrap_err();
    assert!(error.1.contains("two-NaN operand payload"));
    for position in 0..input.columns.len() {
        assert_eq!(bits(&input, position), bits(&before, position));
    }
}

#[test]
fn missing_left_lookup_precedes_missing_right_even_on_empty_input() {
    let mut input = frame(&[], &[]);
    input.columns.clear();
    assert_eq!(
        control_multiply(&mut input, RowIndex::Flat).unwrap_err().1,
        "'sort_index'"
    );
    input.columns.push(column("sort_index", &[]));
    assert_eq!(
        control_multiply(&mut input, RowIndex::Flat).unwrap_err().1,
        "'control_common_index'"
    );
    input.columns.push(column("control_common_index", &[]));
    assert_eq!(control_multiply(&mut input, RowIndex::Flat).unwrap(), 0);
    assert_eq!(input.columns[0].dtype, "float64");
}

#[test]
fn duplicate_required_labels_defer_before_mutation() {
    let mut input = frame(&[0.95], &[0.5]);
    input.columns.push(column("sort_index", &[9.]));
    let before = bits(&input, 1);
    let error = control_multiply(&mut input, RowIndex::Flat).unwrap_err();
    assert_eq!(error.0, "UnsupportedInput");
    assert_eq!(bits(&input, 1), before);
}

#[test]
fn dtype_promotion_and_nullable_cells_retain_original_current_buffer_stage() {
    for dtype in ["int64", "Int64", "Float64", "object", "string"] {
        let mut input = frame(&[1.], &[0.5]);
        input.columns[1].dtype = dtype.into();
        let before = bits(&input, 1);
        assert_eq!(
            control_multiply(&mut input, RowIndex::Flat).unwrap_err().0,
            "UnsupportedInput"
        );
        assert_eq!(bits(&input, 1), before);
    }
    let mut input = frame(&[1.], &[0.5]);
    input.columns[2].values[0] = Cell::Null;
    assert_eq!(
        control_multiply(&mut input, RowIndex::Flat).unwrap_err().0,
        "UnsupportedInput"
    );
    assert_eq!(bits(&input, 1), vec![1f64.to_bits()]);
}

#[test]
fn unqualified_axis_warning_and_native_size_routes_defer_without_mutation() {
    let mut input = frame(&[0.95], &[0.5]);
    assert!(
        control_multiply(&mut input, RowIndex::MultiIndex)
            .unwrap_err()
            .1
            .contains("MultiIndex")
    );
    assert_eq!(bits(&input, 1), vec![0.95f64.to_bits()]);
    // Check size before allocating product buffers or accepting inconsistent data.
    input.index.resize(1_000_001, 0);
    assert!(
        control_multiply(&mut input, RowIndex::Flat)
            .unwrap_err()
            .1
            .contains("numexpr-size")
    );
    assert_eq!(bits(&input, 1), vec![0.95f64.to_bits()]);
}

#[test]
fn never_invokes_control_preparation_or_changes_prepared_control_fields() {
    let mut input = frame(&[0.95], &[0.5]);
    input.columns.push(InputColumn {
        name: "control_AC".into(),
        dtype: "object".into(),
        values: vec![Cell::Text("invalid original control text".into())],
    });
    control_multiply(&mut input, RowIndex::Flat).unwrap();
    assert!(
        matches!(&input.columns[4].values[0], Cell::Text(value) if value=="invalid original control text")
    );
    assert_eq!(bits(&input, 1), vec![0.475f64.to_bits()]);
}
