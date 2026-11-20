use polars_core::prelude::*;
use priva_tabular::acmg_scoring::pedigree_context::{Context, Family, pedigree_context};
use priva_tabular::acmg_scoring::ranking_compatibility::{
    MetadataOrigin, TARGET, ranking_compatibility,
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
fn frame(profiles: &[&str]) -> Input<i64> {
    let n = profiles.len();
    Input {
        columns: vec![
            column(
                "var_plausible_patho_mechs",
                "object",
                profiles.iter().map(|x| Cell::Text((*x).into())).collect(),
            ),
            column(
                "variant_effect",
                "object",
                vec![Cell::Text("predicted_LOF_high_confidence".into()); n],
            ),
            column("variant_lof_score", "int64", vec![Cell::Integer(1); n]),
            column("ACMG_quant_score", "float64", vec![Cell::Float(0.95); n]),
            column("sort_index", "float64", vec![Cell::Float(0.475); n]),
            column("chrom", "object", vec![Cell::Text("X".into()); n]),
        ],
        index: (0..n).map(|i| i as i64 * 3 - 2).collect(),
        axis_metadata: "opaque frame attrs/index/flags".into(),
    }
}
fn context(frame: &Input<i64>, male: bool, female: bool) -> Context<i64> {
    Context {
        pedigree_position: Some(0),
        proband: Some("PROBAND".into()),
        is_male: male,
        is_female: female,
        mask_dtype: "bool",
        mask: DataFrame::new(vec![
            Column::new("value".into(), vec![true; frame.index.len()]),
            Column::new("is_null".into(), vec![false; frame.index.len()]),
        ])
        .unwrap(),
        index: frame.index.clone(),
        axis_metadata: frame.axis_metadata.clone(),
        inherit_frame_metadata: true,
        warning: None,
    }
}
fn bits(column: &Column) -> Vec<bool> {
    column.bool().unwrap().into_no_null_iter().collect()
}

#[test]
fn sex_and_x_history_keep_original_predicate_order() {
    let profiles = [
        "",
        "dominant_LOF",
        "dominant",
        "dominant_DN",
        "x_linked_recessive_LOF",
        "x_linked_dominant_LOF",
        "x_linked_dominant",
    ];
    for (male, female, expected) in [
        (
            false,
            false,
            vec![false, true, true, false, false, false, false],
        ),
        (
            true,
            false,
            vec![false, true, true, false, true, true, true],
        ),
        (
            false,
            true,
            vec![false, true, true, false, false, true, true],
        ),
    ] {
        let mut frame = frame(&profiles);
        let context = context(&frame, male, female);
        let result = ranking_compatibility(&mut frame, &context).unwrap();
        assert_eq!(
            bits(result.role("has_single_alt_lof_compatible").unwrap()),
            expected
        );
        assert_eq!(
            bits(result.role("dn_without_hi_or_ambiguous").unwrap()),
            vec![false, false, false, true, false, false, false]
        );
        assert_eq!(
            frame.columns[4]
                .values
                .iter()
                .map(|cell| match cell {
                    Cell::Float(x) => x.to_bits(),
                    _ => panic!(),
                })
                .collect::<Vec<_>>(),
            vec![0.475f64.to_bits(); 7]
        );
    }
}

#[test]
fn strict_raw_score_neighbors_do_not_use_penalized_sort_score() {
    let mut frame = frame(&["", "", "", "", ""]);
    frame.columns[3].values = [
        f64::from_bits(0.89f64.to_bits() - 1),
        0.89,
        f64::from_bits(0.89f64.to_bits() + 1),
        f64::INFINITY,
        f64::NAN,
    ]
    .into_iter()
    .map(Cell::Float)
    .collect();
    let context = context(&frame, false, false);
    let result = ranking_compatibility(&mut frame, &context).unwrap();
    assert_eq!(
        bits(result.role("high_acmg_lof").unwrap()),
        vec![false, false, true, true, false]
    );
    assert_eq!(result.log_counts, [0, 5, 0, 5, 2]);
}

#[test]
fn nullable_r2_mask_retains_hidden_raw_bits_and_unmasks_false() {
    let mut frame = frame(&["", "", "dominant_LOF", ""]);
    let mut context = context(&frame, false, false);
    context.mask_dtype = "boolean";
    context.mask = DataFrame::new(vec![
        Column::new("value".into(), [true, false, true, true]),
        Column::new("is_null".into(), [true, true, true, false]),
    ])
    .unwrap();
    let result = ranking_compatibility(&mut frame, &context).unwrap();
    assert_eq!(result.penalty_dtype, "boolean");
    assert_eq!(
        bits(result.role("compatibility_penalty_rows").unwrap()),
        vec![true, false, false, true]
    );
    assert_eq!(
        bits(result.role("compatibility_penalty_rows_null").unwrap()),
        vec![true, true, false, false]
    );
    assert_eq!(result.log_counts[4], 1);
    let values = &frame.columns[result.compatibility_column_position].values;
    assert!(matches!(values[0], Cell::Float(1.0)));
    assert!(matches!(values[3], Cell::Float(0.9)));
}

#[test]
fn aliases_and_metadata_origins_are_explicit() {
    let mut frame = frame(&[""]);
    let context = context(&frame, true, false);
    let result = ranking_compatibility(&mut frame, &context).unwrap();
    assert!(std::ptr::eq(
        result.role("pext_modulation_rows").unwrap(),
        result.role("high_acmg_lof").unwrap()
    ));
    assert!(std::ptr::eq(
        result.role("exact_gof").unwrap(),
        result.mechanism.frame.column("is_exact_gof").unwrap()
    ));
    assert!(std::ptr::eq(
        result.role("asserted_lof_effect").unwrap(),
        result.mechanism.frame.column("is_predicted_lof").unwrap()
    ));
    assert_eq!(result.metadata_origin("exact_gof"), MetadataOrigin::Fresh);
    assert_eq!(
        result.metadata_origin("asserted_lof_effect"),
        MetadataOrigin::SourceColumn(1)
    );
    assert_eq!(
        result.source_name("high_acmg_lof"),
        Some("ACMG_quant_score")
    );
    assert_eq!(result.source_name("has_single_alt_lof_compatible"), None);
}

#[test]
fn existing_numeric_gof_and_dn_still_have_fresh_metadata() {
    let mut frame = frame(&[""]);
    frame
        .columns
        .push(column("variant_gof_score", "int64", vec![Cell::Integer(2)]));
    frame
        .columns
        .push(column("variant_dn_score", "int64", vec![Cell::Integer(2)]));
    let context = context(&frame, false, false);
    let result = ranking_compatibility(&mut frame, &context).unwrap();
    for name in [
        "is_exact_gof",
        "is_exact_nonlof",
        "exact_gof",
        "is_exact_dn",
    ] {
        assert_eq!(result.metadata_origin(name), MetadataOrigin::Fresh);
    }
    assert!(bits(result.role("exact_gof").unwrap())[0]);
}

#[test]
fn overwrite_preserves_column_position_and_opaque_duplicate_indices() {
    let mut frame = frame(&["", "dominant_LOF"]);
    frame.index = vec![7, 7];
    frame.columns.insert(
        1,
        column(TARGET, "string", vec![Cell::Text("old".into()); 2]),
    );
    let before = frame
        .columns
        .iter()
        .map(|x| x.name.clone())
        .collect::<Vec<_>>();
    let context = context(&frame, false, false);
    let result = ranking_compatibility(&mut frame, &context).unwrap();
    assert_eq!(result.compatibility_column_position, 1);
    assert_eq!(frame.columns[1].dtype, "float64");
    assert_eq!(
        frame
            .columns
            .iter()
            .map(|x| x.name.clone())
            .collect::<Vec<_>>(),
        before
    );
    assert_eq!(frame.index, vec![7, 7]);
    assert_eq!(result.axis_metadata, "opaque frame attrs/index/flags");
}

#[test]
fn helper_failure_precedes_missing_score_and_does_not_mutate_input() {
    let mut first = frame(&[""]);
    first
        .columns
        .retain(|x| x.name != "variant_effect" && x.name != "ACMG_quant_score");
    let before = format!("{:?}", first.columns);
    let ctx = context(&first, false, false);
    let error = ranking_compatibility(&mut first, &ctx).err().unwrap();
    assert_eq!(error.source_line, Some(537));
    assert_eq!(error.error.0, "KeyError");
    assert!(error.error.1.contains("variant_effect"));
    assert_eq!(format!("{:?}", first.columns), before);
    let mut frame = frame(&[""]);
    frame.columns.retain(|x| x.name != "ACMG_quant_score");
    let context = context(&frame, false, false);
    let error = ranking_compatibility(&mut frame, &context).err().unwrap();
    assert_eq!(error.source_line, Some(584));
}

#[test]
fn unsupported_context_or_cache_routes_reject_atomically() {
    for mode in [0, 1, 2] {
        let mut frame = frame(&[""]);
        let mut context = context(&frame, false, false);
        match mode {
            0 => context.index[0] = 99,
            1 => frame.columns[3].dtype = "Float64".into(),
            _ => frame.columns[2].dtype = "Int64".into(),
        }
        let before = format!("{:?}", frame.columns);
        let error = ranking_compatibility(&mut frame, &context).err().unwrap();
        assert_eq!(error.error.0, "UnsupportedInput");
        assert_eq!(format!("{:?}", frame.columns), before);
    }
}

#[test]
fn consume_existing_r2_context_in_memory_without_replaying_warning() {
    let mut frame = frame(&["", ""]);
    frame.columns.push(column(
        "PROBAND",
        "string",
        vec![Cell::Text("0/1".into()), Cell::Null],
    ));
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
        &frame,
        Some(&pedigree),
        &Family::Value(Cell::Text("F".into())),
    )
    .unwrap();
    assert!(context.warning.is_some());
    let warning = context.warning;
    let result = ranking_compatibility(&mut frame, &context).unwrap();
    assert_eq!(context.warning, warning);
    assert_eq!(result.penalty_dtype, "boolean");
    assert_eq!(result.log_counts[4], 1);
}
