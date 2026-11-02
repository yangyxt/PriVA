//! Common-missense filtering over ordered residue AF lists or the legacy array.
//! Fallback is model-wide, not per domain. Python min/max retain a first NaN,
//! while a later NaN cannot replace an already finite accumulator.
use crate::das_reference::{count_ratio, round_four};

#[derive(Debug)]
pub(crate) enum Fraction {
    IntegerZero,
    Float(f64),
}

#[derive(Debug)]
pub(crate) struct Summary {
    pub(crate) fully_common: u64,
    pub(crate) total: u64,
    pub(crate) fraction: Fraction,
    pub(crate) removed: bool,
    pub(crate) with_data: bool,
}

impl Summary {
    fn empty(with_data: bool) -> Self {
        Self {
            fully_common: 0,
            total: 0,
            fraction: Fraction::IntegerZero,
            removed: false,
            with_data,
        }
    }
}

/// Counts include empty residue lists in the denominator. Zero totals model
/// missing/empty residue data; they never evaluate the removal threshold.
pub(crate) fn residue_counts(total: u64, fully_common: u64, threshold: f64) -> Summary {
    if total == 0 {
        return Summary::empty(false);
    }
    let fraction = count_ratio(fully_common, total).expect("nonzero residue count");
    Summary {
        fully_common,
        total,
        fraction: Fraction::Float(round_four(fraction)),
        removed: fraction > threshold,
        with_data: true,
    }
}

pub(crate) fn filter<R, A, F>(
    residues: Option<R>,
    flat: Option<F>,
    any_residue_data: bool,
    af_common: f64,
    af_rare: f64,
    residue_fraction: f64,
) -> Summary
where
    R: IntoIterator<Item = A>,
    A: IntoIterator<Item = f64>,
    F: IntoIterator<Item = f64>,
{
    let mut total = 0;
    let mut fully_common = 0;
    if let Some(residues) = residues {
        for afs in residues {
            total += 1;
            let mut values = afs.into_iter();
            if let Some(first) = values.next() {
                let mut minimum = first;
                let mut maximum = first;
                for value in values {
                    if value < minimum {
                        minimum = value;
                    }
                    if value > maximum {
                        maximum = value;
                    }
                }
                if maximum > af_common && minimum > af_rare {
                    fully_common += 1;
                }
            }
        }
    }
    if total != 0 {
        return residue_counts(total, fully_common, residue_fraction);
    }
    if !any_residue_data {
        if let Some(flat) = flat {
            for value in flat {
                total += 1;
                fully_common += u64::from(value > af_common);
            }
            if total == 0 {
                return Summary::empty(true);
            }
            let fraction = count_ratio(fully_common, total).expect("nonzero flat count");
            return Summary {
                fully_common,
                total,
                fraction: Fraction::Float(round_four(fraction)),
                removed: fully_common >= 1,
                with_data: true,
            };
        }
    }
    Summary::empty(false)
}
