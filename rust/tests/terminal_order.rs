use priva_tabular::acmg_scoring::terminal_order::*;
use priva_tabular::acmg_variant_mechanism::variant_mechanism_masks::Cell;

fn column(name: &str, dtype: &str, values: Vec<Cell>) -> InputColumn {
    InputColumn {
        name: name.into(),
        dtype: dtype.into(),
        values,
    }
}
fn text(values: &[&str]) -> Vec<Cell> {
    values.iter().map(|x| Cell::Text((*x).into())).collect()
}
fn frame(scores: &[f64]) -> Input<&'static str> {
    let n = scores.len();
    Input {
        columns: vec![
            column("chrom", "object", vec![Cell::Text("1".into()); n]),
            column("pos", "int64", vec![Cell::Integer(100); n]),
            column("ref", "object", vec![Cell::Text("A".into()); n]),
            column("alt", "object", vec![Cell::Text("G".into()); n]),
            column(
                "sort_index",
                "float64",
                scores.iter().copied().map(Cell::Float).collect(),
            ),
            column(
                "ACMG_quant_score",
                "float64",
                scores.iter().copied().map(Cell::Float).collect(),
            ),
            column("variant_id", "object", vec![Cell::Text("same".into()); n]),
        ],
        index: vec!["opaque duplicate"; n],
        axis_metadata: "caller-owned index/attrs".into(),
    }
}

#[test]
fn first_tied_zero_maximum_and_source_tie_order() {
    let input = frame(&[-0., 0.]);
    let output = terminal_order(&input).unwrap();
    assert_eq!(output.source_positions, vec![0, 1]);
    assert!(
        output
            .maxima
            .iter()
            .all(|x| matches!(x, Cell::Float(x) if x.to_bits() == (-0.0f64).to_bits()))
    );
    assert_eq!(output.index, input.index);
    assert_eq!(output.axis_metadata, input.axis_metadata);
    assert_eq!(input.columns[4].values[0], Cell::Float(-0.));
}

#[test]
fn nonfinite_group_and_missing_coordinate_behavior() {
    let mut input = frame(&[f64::NAN, f64::INFINITY, f64::NEG_INFINITY, f64::NAN, 9.]);
    input.columns[1].values = vec![
        Cell::Integer(1),
        Cell::Integer(1),
        Cell::Integer(2),
        Cell::Integer(2),
        Cell::Null,
    ];
    input.columns[1].dtype = "Int64".into();
    let output = terminal_order(&input).unwrap();
    assert_eq!(output.source_positions, vec![1, 0, 2, 3, 4]);
    assert_eq!(
        &output.maxima[..4],
        &[
            Cell::Float(f64::INFINITY),
            Cell::Float(f64::INFINITY),
            Cell::Float(f64::NEG_INFINITY),
            Cell::Float(f64::NEG_INFINITY)
        ]
    );
    assert!(matches!(output.maxima[4], Cell::Float(x) if x.is_nan()));
    assert!(matches!(input.columns[4].values[0], Cell::Float(x) if x.is_nan()));
}

#[test]
fn raw_score_precedes_coordinate_sort_and_rank_identity_is_separate() {
    let mut input = frame(&[0.4, 0.9, 0.8]);
    input.columns[1].values = vec![Cell::Integer(1), Cell::Integer(1), Cell::Integer(2)];
    input.columns[6].values = text(&["b", "a", "a"]);
    let output = terminal_order(&input).unwrap();
    assert_eq!(output.source_positions, vec![1, 0, 2]);
    assert_eq!(output.ranks, vec![1, 2, 1]);
}

#[test]
fn missing_ids_rank_zero_and_do_not_consume_a_rank() {
    let mut input = frame(&[1., 1., 1., 1.]);
    input.columns[6].values = vec![
        Cell::Null,
        Cell::Text("a".into()),
        Cell::Float(f64::NAN),
        Cell::Text("b".into()),
    ];
    let output = terminal_order(&input).unwrap();
    assert_eq!(output.source_positions, vec![0, 1, 2, 3]);
    assert_eq!(output.ranks, vec![0, 1, 0, 2]);
}

#[test]
fn absent_id_actual_colon_join_collision_is_preserved() {
    let mut input = frame(&[1., 1.]);
    input.columns.pop();
    input.columns[0].values = text(&["a:b", "a"]);
    input.columns[1] = column("pos", "object", text(&["c", "b:c"]));
    let output = terminal_order(&input).unwrap();
    assert_eq!(output.ranks, vec![1, 1]);
}

#[test]
fn absent_float_ids_distinguish_zero_signs_but_explicit_ids_do_not() {
    let mut input = frame(&[1., 1.]);
    input.columns.pop();
    input.columns[1] = column("pos", "float64", vec![Cell::Float(0.), Cell::Float(-0.)]);
    assert_eq!(terminal_order(&input).unwrap().ranks, vec![1, 2]);
    input.columns.push(column(
        "variant_id",
        "float64",
        vec![Cell::Float(0.), Cell::Float(-0.)],
    ));
    assert_eq!(terminal_order(&input).unwrap().ranks, vec![1, 1]);
}

#[test]
fn column_overwrite_drop_positions_and_empty_frame() {
    let mut input = frame(&[]);
    input
        .columns
        .insert(0, column("variant_rank", "object", vec![]));
    input
        .columns
        .insert(2, column("max_variant_score", "object", vec![]));
    let output = terminal_order(&input).unwrap();
    assert!(output.source_positions.is_empty());
    assert_eq!(output.temporary_column_position, 2);
    assert_eq!(
        output.output_columns[0],
        ("variant_rank".into(), "int64".into())
    );
    assert!(
        !output
            .output_columns
            .iter()
            .any(|x| x.0 == "max_variant_score")
    );
}

#[test]
fn unsupported_mixed_objects_and_ambiguous_float_text_join_fail_closed() {
    let mut input = frame(&[1., 1.]);
    input.columns[0].values = vec![Cell::Integer(2), Cell::Text("10".into())];
    assert_eq!(terminal_order(&input).err().unwrap().0, "UnsupportedInput");
    input.columns.pop();
    input.columns[0].values = text(&["a:b", "a"]);
    input.columns[1] = column("pos", "float64", vec![Cell::Float(1.), Cell::Float(2.)]);
    assert!(
        terminal_order(&input)
            .err()
            .unwrap()
            .1
            .contains("colon-bearing")
    );
}

#[test]
fn required_lookup_order_precedes_unsupported_dtype() {
    let mut input = frame(&[1.]);
    input
        .columns
        .retain(|c| c.name != "chrom" && c.name != "pos");
    input.columns[0].dtype = "category".into();
    let error = terminal_order(&input).err().unwrap();
    assert_eq!((error.0, error.1.as_str()), ("KeyError", "'chrom'"));
}

#[test]
fn unsigned_ids_preserve_large_integer_identity() {
    let mut input = frame(&[1., 1., 1.]);
    input.columns[6] = column(
        "variant_id",
        "uint64",
        vec![
            Cell::Unsigned(u64::MAX),
            Cell::Unsigned(1 << 63),
            Cell::Unsigned(u64::MAX),
        ],
    );
    assert_eq!(terminal_order(&input).unwrap().ranks, vec![1, 2, 1]);
}

#[test]
fn impossible_typed_missing_cells_and_empty_unsupported_ids_are_deferred() {
    let mut input = frame(&[1.]);
    input.columns[1].values[0] = Cell::Float(f64::NAN);
    input.columns[1].dtype = "Int64".into();
    assert_eq!(terminal_order(&input).err().unwrap().0, "UnsupportedInput");
    input.columns[1].values[0] = Cell::Integer(1);
    input.columns[4].values[0] = Cell::Null;
    assert_eq!(terminal_order(&input).err().unwrap().0, "UnsupportedInput");
    let mut empty = frame(&[]);
    empty.columns[6].dtype = "category".into();
    assert_eq!(terminal_order(&empty).err().unwrap().0, "UnsupportedInput");
}

#[test]
fn nullable_integer_maximum_keeps_values_above_f64_precision_and_null_groups() {
    let mut input = frame(&[1., 1., 1.]);
    input.columns[4] = column(
        "sort_index",
        "Int64",
        vec![
            Cell::Integer((1 << 53) + 1),
            Cell::Integer(1 << 53),
            Cell::Null,
        ],
    );
    input.columns[1] = column(
        "pos",
        "Int64",
        vec![Cell::Integer(1), Cell::Integer(1), Cell::Null],
    );
    let output = terminal_order(&input).unwrap();
    assert_eq!(output.maximum_dtype, "Int64");
    assert_eq!(
        output.maxima,
        vec![
            Cell::Integer((1 << 53) + 1),
            Cell::Integer((1 << 53) + 1),
            Cell::Null
        ]
    );
}

#[test]
fn numpy_integer_maximum_promotes_only_when_source_has_excluded_groups() {
    let mut input = frame(&[1., 1.]);
    input.columns[4] = column(
        "sort_index",
        "uint64",
        vec![Cell::Unsigned(u64::MAX), Cell::Unsigned(0)],
    );
    input.columns[1] = column("pos", "Int64", vec![Cell::Integer(1), Cell::Integer(2)]);
    let output = terminal_order(&input).unwrap();
    assert_eq!(output.maximum_dtype, "uint64");
    assert_eq!(
        output.maxima,
        vec![Cell::Unsigned(u64::MAX), Cell::Unsigned(0)]
    );
    input.columns[1].values[1] = Cell::Null;
    let output = terminal_order(&input).unwrap();
    assert_eq!(output.maximum_dtype, "float64");
    assert!(
        matches!(output.maxima[0], Cell::Float(x) if x.to_bits() == (u64::MAX as f64).to_bits())
    );
    assert!(matches!(output.maxima[1], Cell::Float(x) if x.is_nan()));
}
