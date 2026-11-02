//! Sorted rescue task planning and ordered FDR results. Only these two rescue
//! iterators use the approved sort; upstream KS/reference order is untouched.
use crate::das_rescue::{self, Action, Entry};
use crate::fdr;
use std::collections::{BTreeSet, HashMap};

#[derive(Clone, Copy)]
pub(crate) enum Mode {
    GoMf,
    NonGoMf,
}

pub(crate) struct Task<'a> {
    pub(crate) domain: &'a str,
    pub(crate) query: &'a [f64],
    pub(crate) composite: &'a [f64],
    /// None means the global composite. An existing empty typed array remains
    /// selected, so the numerical worker can raise its original failure.
    pub(crate) composite_type: Option<&'a str>,
    pub(crate) seed: u32,
}

pub(crate) struct Plan<'a> {
    pub(crate) eligible: usize,
    pub(crate) tasks: Vec<Task<'a>>,
    pub(crate) direct: Vec<&'a str>,
}

pub(crate) struct Inputs<'a> {
    pub(crate) scores: &'a HashMap<String, Vec<f64>>,
    pub(crate) global: &'a [f64],
    pub(crate) by_type: &'a HashMap<String, Vec<f64>>,
    pub(crate) minimum_median: f64,
    pub(crate) seed: u32,
}

pub(crate) fn plan<'a, E>(
    mode: Mode,
    tolerant: impl IntoIterator<Item = &'a str>,
    inputs: Inputs<'a>,
    mut lookup: impl FnMut(&str) -> Result<Option<Entry>, E>,
) -> Result<Plan<'a>, E> {
    let mut plan = Plan {
        eligible: 0,
        tasks: Vec::new(),
        direct: Vec::new(),
    };
    // Python string order equals UTF-8 lexicographic order for scalar Unicode.
    // Deduplication is the original input set's membership, not a score merge.
    for domain in tolerant.into_iter().collect::<BTreeSet<_>>() {
        let Some(query) = inputs.scores.get(domain) else {
            continue;
        };
        let entry = lookup(domain)?;
        let selected = match mode {
            Mode::GoMf => das_rescue::go_mf(Some(query), entry.as_ref(), inputs.minimum_median),
            Mode::NonGoMf => {
                das_rescue::non_go_mf(domain, Some(query), entry.as_ref(), inputs.minimum_median)
            }
        };
        plan.eligible += usize::from(selected.eligible);
        match selected.action {
            Action::Skip => {}
            Action::Direct => plan.direct.push(domain),
            Action::Test(kind) => {
                let (composite_type, composite) = inputs
                    .by_type
                    .get_key_value(&kind)
                    .map_or((None, inputs.global), |(name, values)| {
                        (Some(name.as_str()), values.as_slice())
                    });
                plan.tasks.push(Task {
                    domain,
                    query,
                    composite,
                    composite_type,
                    seed: inputs.seed.wrapping_add(plan.tasks.len() as u32),
                });
            }
        }
    }
    Ok(plan)
}

pub(crate) struct Decision {
    pub(crate) task_index: usize,
    pub(crate) adjusted: f64,
    pub(crate) rescued: bool,
}

pub(crate) struct Decisions {
    /// Keep candidate order regardless of the workers' completion order.
    pub(crate) tested: Vec<Decision>,
    /// Python discards direct rescues on either early return, even when direct
    /// candidates exist. The caller unions them only after a valid FDR result.
    pub(crate) include_direct: bool,
}

impl Plan<'_> {
    pub(crate) fn workers(&self, threads: i64) -> Result<usize, &'static str> {
        if self.tasks.is_empty() {
            return Ok(0);
        }
        if threads <= 0 {
            return Err("Number of processes must be at least 1");
        }
        Ok((threads as usize).min(self.tasks.len()))
    }

    /// Probabilities are joined back to task order by the caller. Worker field
    /// values/order stay in their original records; append these two fields.
    pub(crate) fn finish(
        &self,
        probabilities: &[f64],
        alpha: f64,
    ) -> Result<Decisions, &'static str> {
        if probabilities.len() != self.tasks.len() {
            return Err("incomplete rescue task results");
        }
        let indices: Vec<_> = probabilities
            .iter()
            .enumerate()
            .filter_map(|(index, value)| (!value.is_nan()).then_some(index))
            .collect();
        if indices.is_empty() {
            return Ok(Decisions {
                tested: Vec::new(),
                include_direct: false,
            });
        }
        let selected: Vec<_> = indices.iter().map(|&index| probabilities[index]).collect();
        let (rejected, adjusted) = fdr::statsmodels_bh(&selected, alpha)?;
        Ok(Decisions {
            tested: indices
                .into_iter()
                .zip(adjusted)
                .zip(rejected)
                .map(|((task_index, adjusted), rejected)| Decision {
                    task_index,
                    adjusted,
                    rescued: !rejected,
                })
                .collect(),
            include_direct: true,
        })
    }
}
