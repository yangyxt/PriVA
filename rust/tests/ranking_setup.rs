use priva_tabular::acmg_scoring::ranking_setup::ranking_setup;
use priva_tabular::acmg_scoring::terminal_order::{Input, InputColumn};
use priva_tabular::acmg_variant_mechanism::variant_mechanism_masks::Cell;

fn col(name: &str, dtype: &str, values: Vec<Cell>) -> InputColumn {
    InputColumn {
        name: name.into(),
        dtype: dtype.into(),
        values,
    }
}
fn frame(values: Vec<Cell>) -> Input<(&'static str, i64)> {
    let n = values.len();
    Input {
        columns: vec![
            col("ACMG_quant_score", "float64", vec![Cell::Float(-0.); n]),
            col("BIOTYPE", "object", values),
        ],
        index: vec![("duplicate", 7); n],
        axis_metadata: "caller axes/attrs".into(),
    }
}

#[test]
fn exact_protein_coding_token_only_and_source_order() {
    let input = frame(vec![
        Cell::Text("protein_coding".into()),
        Cell::Text("protein_coding ".into()),
        Cell::Text("Protein_coding".into()),
        Cell::Null,
        Cell::Float(f64::NAN),
        Cell::Integer(1),
        Cell::Text("protein_coding".into()),
    ]);
    let plan = ranking_setup(&input).unwrap();
    assert_eq!(plan.source_positions, vec![0, 6]);
    assert_eq!(plan.index, vec![("duplicate", 7); 2]);
    assert_eq!(plan.axis_metadata, input.axis_metadata);
    assert!(
        matches!(input.columns[0].values[0],Cell::Float(x) if x.to_bits()==(-0.0f64).to_bits())
    );
    assert_eq!(plan.columns.last().unwrap().source_column, 0);
}

#[test]
fn legacy_drop_order_and_existing_score_position_are_preserved() {
    let mut input = frame(vec![Cell::Text("protein_coding".into())]);
    input.columns.insert(
        0,
        col(
            "zygosity_inheritance_mechanism_index",
            "category",
            vec![Cell::Null],
        ),
    );
    input.columns.insert(
        1,
        col("sort_index", "object", vec![Cell::Text("stale".into())]),
    );
    input.columns.push(col(
        "haplo_insuf_index",
        "object",
        vec![Cell::Text("legacy".into())],
    ));
    let plan = ranking_setup(&input).unwrap();
    assert_eq!(plan.dropped_columns, vec![4, 0]);
    assert_eq!(
        plan.columns
            .iter()
            .map(|x| x.name.as_str())
            .collect::<Vec<_>>(),
        vec!["sort_index", "ACMG_quant_score", "BIOTYPE"]
    );
    assert!(plan.columns[0].copy_score);
    assert_eq!(plan.columns[0].source_column, 2);
    assert_eq!(plan.columns[0].dtype, "float64");
}

#[test]
fn raw_dtype_and_bits_are_copied_without_numerical_operations() {
    let mut input = frame(vec![Cell::Text("protein_coding".into()); 2]);
    input.columns[0] = col(
        "ACMG_quant_score",
        "UInt64",
        vec![Cell::Unsigned(u64::MAX), Cell::Null],
    );
    let plan = ranking_setup(&input).unwrap();
    assert_eq!(plan.columns.last().unwrap().dtype, "UInt64");
    assert_eq!(
        input.columns[0].values,
        vec![Cell::Unsigned(u64::MAX), Cell::Null]
    );
}

#[test]
fn required_lookup_order_also_applies_to_empty_frames() {
    let input: Input<()> = Input {
        columns: vec![],
        index: vec![],
        axis_metadata: String::new(),
    };
    let error = ranking_setup(&input).err().unwrap();
    assert_eq!(
        (error.0, error.1.as_str()),
        ("KeyError", "'ACMG_quant_score'")
    );
    let input = Input {
        columns: vec![col("ACMG_quant_score", "float64", vec![])],
        ..input
    };
    assert_eq!(ranking_setup(&input).err().unwrap().1, "'BIOTYPE'");
}

#[test]
fn empty_frame_still_appends_sort_column_and_keeps_opaque_axis() {
    let input = frame(vec![]);
    let plan = ranking_setup(&input).unwrap();
    assert!(plan.source_positions.is_empty());
    assert_eq!(plan.columns.len(), 3);
    assert_eq!(plan.columns[2].name, "sort_index");
    assert_eq!(plan.axis_metadata, "caller axes/attrs");
}

#[test]
fn unrelated_duplicate_labels_are_retained_but_affected_labels_defer() {
    let mut input = frame(vec![Cell::Text("protein_coding".into())]);
    input
        .columns
        .push(col("unused", "object", vec![Cell::Text("first".into())]));
    input
        .columns
        .push(col("unused", "object", vec![Cell::Text("second".into())]));
    let plan = ranking_setup(&input).unwrap();
    assert_eq!(plan.columns[2].source_column, 2);
    assert_eq!(plan.columns[3].source_column, 3);
    input.columns.push(input.columns[1].clone());
    assert_eq!(ranking_setup(&input).err().unwrap().0, "UnsupportedInput");
}

#[test]
fn nullable_boolean_biotype_missingness_never_becomes_a_retained_row() {
    let mut input = frame(vec![Cell::Bool(true), Cell::Null, Cell::Bool(false)]);
    input.columns[1].dtype = "boolean".into();
    assert!(ranking_setup(&input).unwrap().source_positions.is_empty());
    assert!(matches!(input.columns[1].values[1], Cell::Null));
}

#[test]
fn unsupported_dtype_or_inconsistent_primitive_cell_fails_closed() {
    let mut input = frame(vec![Cell::Text("protein_coding".into())]);
    input.columns[1].dtype = "category".into();
    assert_eq!(ranking_setup(&input).err().unwrap().0, "UnsupportedInput");
    input.columns[1].dtype = "int64".into();
    assert_eq!(ranking_setup(&input).err().unwrap().0, "UnsupportedInput");
}
