//! The narrow scalar forms emitted by the ClinVar gene-statistics producer.
//! These are descriptions interpreted by Rust, never executable Python globals.
//! The ordinary pickle reader does not enable any of these opcode handlers.
use super::{Result, Value, invalid};
use std::cell::RefCell;
use std::collections::HashSet;
use std::rc::Rc;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Scalar {
    PandasNA,
    NumpyI32(i32),
    /// IEEE-754 bits, without promotion or NaN canonicalization.
    NumpyF32(u32),
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum Global {
    Dtype,
    Scalar,
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum Kind {
    I32,
    F32,
}

#[derive(Clone, Debug)]
pub(crate) struct Dtype {
    kind: Kind,
    built: bool,
}

#[derive(Clone, Debug)]
pub(crate) enum Object {
    Scalar(Scalar),
    Global(Global),
    Dtype(Rc<RefCell<Dtype>>),
    Bytes([u8; 4]),
}

impl Object {
    pub(super) fn type_name(&self) -> &'static str {
        match self {
            Self::Scalar(Scalar::PandasNA) => "NAType",
            Self::Scalar(Scalar::NumpyI32(_)) => "int32",
            Self::Scalar(Scalar::NumpyF32(_)) => "float32",
            _ => "restricted reducer intermediate",
        }
    }
    pub(super) fn scalar(&self) -> Option<Scalar> {
        if let Self::Scalar(scalar) = self {
            Some(*scalar)
        } else {
            None
        }
    }
}

pub(super) fn global(module: Value, name: Value) -> Result<Value> {
    let value = match (module.as_str(), name.as_str()) {
        (Some("pandas._libs.missing"), Some("NA")) => Object::Scalar(Scalar::PandasNA),
        (Some("numpy"), Some("dtype")) => Object::Global(Global::Dtype),
        (Some("numpy.core.multiarray"), Some("scalar")) => Object::Global(Global::Scalar),
        _ => return Err(invalid("schema-5 pickle names an unapproved global")),
    };
    Ok(Value::Clinvar(value))
}

pub(super) fn bytes(value: [u8; 4]) -> Value {
    Value::Clinvar(Object::Bytes(value))
}

pub(super) fn reduce(target: Value, args: Value) -> Result<Value> {
    let Value::Tuple(args) = args else {
        return Err(invalid("schema-5 reducer arguments must be a tuple"));
    };
    let value = match (target, args.as_slice()) {
        (
            Value::Clinvar(Object::Global(Global::Dtype)),
            [Value::Str(name), Value::Bool(false), Value::Bool(true)],
        ) => {
            let kind = match name.as_ref() {
                "i4" => Kind::I32,
                "f4" => Kind::F32,
                _ => return Err(invalid("schema-5 accepts only i4/f4 dtypes")),
            };
            Object::Dtype(Rc::new(RefCell::new(Dtype { kind, built: false })))
        }
        (
            Value::Clinvar(Object::Global(Global::Scalar)),
            [
                Value::Clinvar(Object::Dtype(dtype)),
                Value::Clinvar(Object::Bytes(bytes)),
            ],
        ) => {
            let dtype = dtype.borrow();
            if !dtype.built {
                return Err(invalid("schema-5 scalar dtype has no validated state"));
            }
            Object::Scalar(match dtype.kind {
                Kind::I32 => Scalar::NumpyI32(i32::from_le_bytes(*bytes)),
                Kind::F32 => Scalar::NumpyF32(u32::from_le_bytes(*bytes)),
            })
        }
        _ => {
            return Err(invalid(
                "schema-5 pickle has unapproved reducer arguments or target",
            ));
        }
    };
    Ok(Value::Clinvar(value))
}

pub(super) fn build(target: &Value, state: Value) -> Result<()> {
    let (Value::Clinvar(Object::Dtype(dtype)), Value::Tuple(state)) = (target, state) else {
        return Err(invalid(
            "schema-5 BUILD needs a dtype and its exact state tuple",
        ));
    };
    if !matches!(state.as_slice(),
        [Value::Int(3), Value::Str(endian), Value::None, Value::None, Value::None,
         Value::Int(-1), Value::Int(-1), Value::Int(0)] if endian.as_ref() == "<")
    {
        return Err(invalid(
            "schema-5 dtype state is not the approved little-endian i4/f4 form",
        ));
    }
    let mut dtype = dtype.borrow_mut();
    if dtype.built {
        return Err(invalid("schema-5 dtype state may be assigned only once"));
    }
    dtype.built = true;
    Ok(())
}

/// Intermediates are allowed only while decoding a scalar, never as output.
/// Visit aliases once so a malicious cycle cannot recurse without a bound.
pub(super) fn validate(root: &Value) -> Result<()> {
    let mut pending = vec![root.clone()];
    let mut seen = HashSet::new();
    while let Some(value) = pending.pop() {
        match value {
            Value::Clinvar(Object::Scalar(_)) => {}
            Value::Clinvar(_) => {
                return Err(invalid(
                    "schema-5 pickle retains an incomplete reducer value",
                ));
            }
            Value::Tuple(items) => {
                if seen.insert((0, Rc::as_ptr(&items) as usize)) {
                    pending.extend(items.iter().cloned());
                }
            }
            Value::List(items) => {
                if seen.insert((1, Rc::as_ptr(&items) as usize)) {
                    pending.extend(items.borrow().iter().cloned());
                }
            }
            Value::Dict(dict) => {
                if seen.insert((2, Rc::as_ptr(&dict) as usize)) {
                    for (key, value) in &dict.borrow().entries {
                        pending.push(key.clone());
                        pending.push(value.clone());
                    }
                }
            }
            Value::Set(set) => {
                if seen.insert((3, Rc::as_ptr(&set) as usize)) {
                    pending.extend(set.borrow().items.iter().cloned());
                }
            }
            Value::FrozenSet(set) => {
                if seen.insert((4, Rc::as_ptr(&set) as usize)) {
                    pending.extend(set.items.iter().cloned());
                }
            }
            _ => {}
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pickle_reader::{load, load_clinvar_gene};

    fn string(value: &str) -> Vec<u8> {
        let mut out = vec![0x8c, value.len() as u8];
        out.extend_from_slice(value.as_bytes());
        out
    }
    fn named(module: &str, name: &str) -> Vec<u8> {
        let mut out = string(module);
        out.extend(string(name));
        out.push(0x93);
        out
    }
    fn dtype(kind: &str, endian: &str) -> Vec<u8> {
        let mut out = named("numpy", "dtype");
        out.extend(string(kind));
        out.extend_from_slice(b"\x89\x88\x87R(K\x03");
        out.extend(string(endian));
        out.extend_from_slice(b"NNNJ\xff\xff\xff\xffJ\xff\xff\xff\xffK\0tb");
        out
    }
    fn scalar(kind: &str, bits: u32) -> Vec<u8> {
        let mut out = vec![0x80, 5];
        out.extend(named("numpy.core.multiarray", "scalar"));
        out.extend(dtype(kind, "<"));
        out.extend_from_slice(b"C\x04");
        out.extend_from_slice(&bits.to_le_bytes());
        out.extend_from_slice(b"\x86R.");
        out
    }

    #[test]
    fn explicitly_selected_scalars_keep_kind_and_bits() {
        let mut na = named("pandas._libs.missing", "NA");
        na.push(b'.');
        assert_eq!(
            load_clinvar_gene(&na).unwrap().clinvar_scalar(),
            Some(Scalar::PandasNA)
        );
        assert!(load(&na).is_err());
        for bits in [
            0,
            1,
            0x8000_0000,
            0x7f80_0000,
            0xff80_0000,
            0x7fc0_0042,
            0x7f80_0042,
        ] {
            let bytes = scalar("f4", bits);
            assert_eq!(
                load_clinvar_gene(&bytes).unwrap().clinvar_scalar(),
                Some(Scalar::NumpyF32(bits))
            );
            assert!(load(&bytes).is_err());
        }
        for value in [0, 1, -1, i32::MIN, i32::MAX] {
            let bytes = scalar("i4", value as u32);
            assert_eq!(
                load_clinvar_gene(&bytes).unwrap().clinvar_scalar(),
                Some(Scalar::NumpyI32(value))
            );
            assert!(load(&bytes).is_err());
        }
    }

    #[test]
    fn refuses_other_globals_reducers_dtypes_states_and_incomplete_values() {
        for (module, name) in [
            ("os", "system"),
            ("builtins", "eval"),
            ("numpy", "ndarray"),
            ("numpy._core.multiarray", "scalar"),
            ("pandas._libs.missing", "NAType"),
        ] {
            let mut bytes = named(module, name);
            bytes.extend_from_slice(b")R.");
            assert!(load_clinvar_gene(&bytes).is_err());
        }
        for kind in ["i8", "f8", "u4", "O", "<i4", ">i4"] {
            assert!(load_clinvar_gene(&scalar(kind, 0)).is_err());
        }
        for endian in [">", "=", "|"] {
            let mut bytes = dtype("i4", endian);
            bytes.push(b'.');
            assert!(load_clinvar_gene(&bytes).is_err());
        }
        let base = scalar("i4", 0);
        // Reducer booleans and every field of the dtype state are exact, not
        // Python-equal substitutes or values passed into a Python callable.
        let patterns: &[(&[u8], &[u8])] = &[
            (b"\x89\x88\x87R", b"\x88\x88\x87R"),
            (b"\x89\x88\x87R", b"K\0\x88\x87R"),
            (b"(K\x03", b"(K\x02"),
            (b"NNN", b"NN}"),
            (b"K\0tb", b"K\x01tb"),
            (b"C\x04", b"C\x03"),
            (b"tb", b"t"),
        ];
        for (old, new) in patterns {
            let at = base
                .windows(old.len())
                .position(|part| part == *old)
                .unwrap();
            let mut bytes = base[..at].to_vec();
            bytes.extend_from_slice(new);
            bytes.extend_from_slice(&base[at + old.len()..]);
            assert!(load_clinvar_gene(&bytes).is_err(), "{old:?}");
        }
        for mut bytes in [
            named("numpy", "dtype"),
            dtype("i4", "<"),
            b"C\x04abcd".to_vec(),
        ] {
            bytes.push(b'.');
            assert!(load_clinvar_gene(&bytes).is_err());
        }
        // A completed dtype remains shared through the memo, but cannot be
        // changed by a second BUILD or called with arbitrary arguments.
        let mut repeated = dtype("i4", "<");
        repeated.extend_from_slice(b")b.");
        assert!(load_clinvar_gene(&repeated).is_err());
        let mut na_call = named("pandas._libs.missing", "NA");
        na_call.extend_from_slice(b")R.");
        assert!(load_clinvar_gene(&na_call).is_err());
    }
}
