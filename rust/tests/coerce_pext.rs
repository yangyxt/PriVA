use polars_core::prelude::*;
use priva_tabular::acmg_consequence::coerce_pext_value::{
    BooleanOrigin, Value, coerce_pext_value, coerce_values, from_polars,
};
use priva_tabular::acmg_variant_mechanism::variant_mechanism_masks::Cell;

fn text(value: &str) -> f64 {
    coerce_pext_value(&Value::Primitive(Cell::Text(value.into()))).unwrap()
}

#[test]
fn missing_first_resets_numeric_nan_payload_and_sign() {
    for bits in [
        0x7ff8_0000_0000_0123,
        0xfff8_0000_0000_0456,
        0x7ff0_0000_0000_0001,
    ] {
        let input = Value::Primitive(Cell::Float(f64::from_bits(bits)));
        assert_eq!(
            coerce_pext_value(&input).unwrap().to_bits(),
            0x7ff8_0000_0000_0000
        );
        let Value::Primitive(Cell::Float(unchanged)) = input else {
            unreachable!()
        };
        assert_eq!(unchanged.to_bits(), bits);
    }
    assert!(
        coerce_pext_value(&Value::Primitive(Cell::Null))
            .unwrap()
            .is_nan()
    );
}

#[test]
fn python_bool_and_numpy_bool_are_different_scalar_types() {
    assert_eq!(
        coerce_pext_value(&Value::Primitive(Cell::Bool(true))).unwrap(),
        1.0
    );
    assert_eq!(
        coerce_pext_value(&Value::Primitive(Cell::Bool(false))).unwrap(),
        0.0
    );
    for value in [true, false] {
        assert!(
            coerce_pext_value(&Value::NumpyBool(value))
                .unwrap()
                .is_nan()
        );
    }
}

#[test]
fn numeric_branch_retains_nonfinite_zero_and_float32_conversion() {
    for number in [
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0,
        -0.0,
        f64::MIN_POSITIVE,
    ] {
        assert_eq!(
            coerce_pext_value(&Value::Primitive(Cell::Float(number)))
                .unwrap()
                .to_bits(),
            number.to_bits()
        );
    }
    let small = f32::from_bits(0x3dcccccd);
    assert_eq!(
        coerce_pext_value(&Value::Float32(small)).unwrap(),
        f64::from(small)
    );
    assert_eq!(
        coerce_pext_value(&Value::Float32(f32::from_bits(0xffc00123)))
            .unwrap()
            .to_bits(),
        0x7ff8_0000_0000_0000
    );
}

#[test]
fn source_skip_tokens_are_case_sensitive_and_nan_order_is_retained() {
    assert_eq!(text("NaN&nan&NA&.&1"), 1.0);
    assert!(text("NAN&1&2").is_nan());
    assert_eq!(text("1&NAN&2"), 2.0);
    assert_eq!(text("-nan&inf").to_bits(), 0xfff8_0000_0000_0000);
    assert_eq!(text("1&-nan&inf"), f64::INFINITY);
    assert_eq!(text("nan&-nan&1").to_bits(), 0xfff8_0000_0000_0000);
}

#[test]
fn python_max_retains_the_first_equal_signed_zero() {
    assert_eq!(text("-0&0&-0").to_bits(), (-0.0f64).to_bits());
    assert_eq!(text("0&-0&0").to_bits(), 0.0f64.to_bits());
    assert_eq!(text("-0,0.0").to_bits(), (-0.0f64).to_bits());
}

#[test]
fn underscores_are_allowed_only_between_decimal_digits() {
    assert_eq!(text("1_2.3_4e+0_1"), 123.4);
    for value in [
        "_1", "1_", "1__2", "1_.0", "1._0", "1_e2", "1e_2", "1e+_2", "1e2_",
    ] {
        assert!(text(value).is_nan(), "{value}");
    }
    assert_eq!(text("bad&1e_2&.25"), 0.25);
}

#[test]
fn unicode_decimal_digits_use_the_frozen_python_unicode_version() {
    assert_eq!(text("١.٢_٣"), 1.23);
    assert_eq!(text("𝟙.𝟚e𝟛"), 1200.0);
    assert_eq!(text("\u{11f51}.\u{11f52}"), 1.2);
    assert_eq!(text("\u{1e4f1}.\u{1e4f2}"), 1.2);
    for value in ["²", "Ⅷ", "½", "①", "−1", "\u{10d41}"] {
        assert!(text(value).is_nan(), "{value}");
    }
}

#[test]
fn python_strip_and_full_float_consumption_are_preserved() {
    assert_eq!(text("\u{001c}\u{00a0}.125\u{202f}\u{001f}"), 0.125);
    assert_eq!(text("1,2&.5"), 2.0);
    for value in [
        "1\0",
        "1\0x",
        "1 2",
        "1;2",
        "1|2",
        "\u{200b}.5",
        "\u{feff}.5",
        "0x1.8p2",
    ] {
        assert!(text(value).is_nan(), "{value}");
    }
}

#[test]
fn integer_types_round_and_overflow_in_the_numeric_branch() {
    assert_eq!(
        coerce_pext_value(&Value::Primitive(Cell::Unsigned(u64::MAX))).unwrap(),
        18446744073709551616.0
    );
    assert_eq!(
        coerce_pext_value(&Value::BigInteger("9007199254740993".into())).unwrap(),
        9007199254740992.0
    );
    assert_eq!(
        coerce_pext_value(&Value::BigInteger("-0".into()))
            .unwrap()
            .to_bits(),
        0.0f64.to_bits()
    );
    let error = coerce_pext_value(&Value::BigInteger(format!("1{}", "0".repeat(309)))).unwrap_err();
    assert_eq!(error.0, "OverflowError");
    assert_eq!(error.1, "int too large to convert to float");
    assert_eq!(text("1e309"), f64::INFINITY); // A numeric integer overflow is not string-float overflow.
}

#[test]
fn invalid_big_integer_transport_is_explicit_fallback_and_first_error_wins() {
    assert_eq!(
        coerce_pext_value(&Value::BigInteger("1_000".into()))
            .unwrap_err()
            .0,
        "UnfinishedParity"
    );
    let inputs = [
        Value::BigInteger(format!("1{}", "0".repeat(309))),
        Value::BigInteger("bad".into()),
    ];
    assert_eq!(
        coerce_values(&inputs, "value").unwrap_err().0,
        "OverflowError"
    );
}

#[test]
fn polars_nulls_become_nan_and_numeric_input_bits_are_not_mutated() {
    let input = Column::new(
        "retained".into(),
        [
            Some(-0.0),
            None,
            Some(f64::from_bits(0xfff8_0000_0000_0123)),
            Some(f64::INFINITY),
        ],
    );
    let before: Vec<_> = input
        .f64()
        .unwrap()
        .into_iter()
        .map(|v| v.map(f64::to_bits))
        .collect();
    let output = from_polars(&input, BooleanOrigin::Python).unwrap();
    assert_eq!(output.name(), input.name());
    assert_eq!(output.dtype(), &DataType::Float64);
    assert_eq!(output.null_count(), 0);
    assert_eq!(
        output
            .f64()
            .unwrap()
            .into_no_null_iter()
            .map(f64::to_bits)
            .collect::<Vec<_>>(),
        vec![
            (-0.0f64).to_bits(),
            0x7ff8_0000_0000_0000,
            0x7ff8_0000_0000_0000,
            f64::INFINITY.to_bits()
        ]
    );
    assert_eq!(
        input
            .f64()
            .unwrap()
            .into_iter()
            .map(|v| v.map(f64::to_bits))
            .collect::<Vec<_>>(),
        before
    );
}

#[test]
fn polars_strings_preserve_row_order_and_boolean_origin_is_explicit() {
    let input = Column::new("source_name".into(), ["-nan&1", "1&NAN&2", "NaN&1", "-0&0"]);
    let output = from_polars(&input, BooleanOrigin::Python).unwrap();
    assert_eq!(
        output
            .f64()
            .unwrap()
            .into_no_null_iter()
            .map(f64::to_bits)
            .collect::<Vec<_>>(),
        vec![
            0xfff8_0000_0000_0000,
            2.0f64.to_bits(),
            1.0f64.to_bits(),
            (-0.0f64).to_bits()
        ]
    );
    let booleans = Column::new("bools".into(), [true, false]);
    assert_eq!(
        from_polars(&booleans, BooleanOrigin::Python)
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<_>>(),
        vec![1.0, 0.0]
    );
    assert!(
        from_polars(&booleans, BooleanOrigin::Numpy)
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .all(f64::is_nan)
    );
}

#[test]
fn empty_and_null_polars_columns_keep_output_dtype_and_name() {
    for length in [0, 3] {
        let input = Column::from(Series::full_null("missing".into(), length, &DataType::Null));
        let output = from_polars(&input, BooleanOrigin::Python).unwrap();
        assert_eq!(output.name().as_str(), "missing");
        assert_eq!(output.dtype(), &DataType::Float64);
        assert_eq!(output.len(), length);
        assert_eq!(output.null_count(), 0);
    }
}

#[test]
fn unsupported_polars_dtype_is_rejected_even_when_empty() {
    let input = Column::from(Series::full_null(
        "unsupported".into(),
        0,
        &DataType::List(Box::new(DataType::Int64)),
    ));
    assert_eq!(
        from_polars(&input, BooleanOrigin::Python).unwrap_err().0,
        "UnfinishedParity"
    );
}
