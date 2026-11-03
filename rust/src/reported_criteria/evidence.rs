//! Source-order evidence, within-family majority and exact formatting plans.
//!
//! Family identifiers are assigned by the host's ordinary Python string set.
//! The actual CPython intersection iteration order is an explicit input. Rust
//! never substitutes a hash iteration order or sorts a source set.

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Side { Ps3, Bs3 }

#[derive(Clone, Debug)]
pub struct Evidence {
    pub ordinal: usize,
    pub family: usize,
    pub side: Side,
    /// Python evaluates primitive formatting at the source formatting seam.
    pub formatted: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Conflict {
    pub family: usize,
    pub ps3: usize,
    pub bs3: usize,
    pub winner: Option<Side>,
    pub ps3_ordinals: Vec<usize>,
    pub bs3_ordinals: Vec<usize>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Resolution {
    pub conflicts: Vec<Conflict>,
    pub final_ps3: Vec<usize>,
    pub final_bs3: Vec<usize>,
    pub clean_ps3: Vec<usize>,
    pub clean_bs3: Vec<usize>,
    pub ps3_families: Vec<usize>,
    pub bs3_families: Vec<usize>,
    pub publications: [&'static str; 3],
}

pub fn resolve(evidence: &[Evidence], conflict_order: &[usize]) -> Result<Resolution, String> {
    if evidence.iter().enumerate().any(|(i, item)| item.ordinal != i) {
        return Err("noncontiguous evidence ordinal".into());
    }
    let size = evidence.iter().map(|item| item.family + 1).max().unwrap_or(0);
    let mut ps3 = vec![0_usize; size];
    let mut bs3 = vec![0_usize; size];
    for item in evidence {
        match item.side { Side::Ps3 => ps3[item.family] += 1, Side::Bs3 => bs3[item.family] += 1 }
    }
    let mut visited = vec![false; size];
    for &family in conflict_order {
        if family >= size || visited[family] || ps3[family] == 0 || bs3[family] == 0 {
            return Err("invalid source conflict order".into());
        }
        visited[family] = true;
    }
    if (0..size).any(|family| (ps3[family] > 0 && bs3[family] > 0) != visited[family]) {
        return Err("incomplete source conflict order".into());
    }
    let clean_ps3: Vec<usize> = evidence.iter().filter(|item| item.side == Side::Ps3 && bs3[item.family] == 0).map(|item| item.ordinal).collect();
    let clean_bs3: Vec<usize> = evidence.iter().filter(|item| item.side == Side::Bs3 && ps3[item.family] == 0).map(|item| item.ordinal).collect();
    let mut final_ps3 = clean_ps3.clone();
    let mut final_bs3 = clean_bs3.clone();
    let mut conflicts = Vec::with_capacity(conflict_order.len());
    for &family in conflict_order {
        let winner = if ps3[family] > bs3[family] { Some(Side::Ps3) }
            else if bs3[family] > ps3[family] { Some(Side::Bs3) } else { None };
        let chosen = match winner { Some(Side::Ps3) => Some(&mut final_ps3), Some(Side::Bs3) => Some(&mut final_bs3), None => None };
        if let Some(chosen) = chosen {
            chosen.extend(evidence.iter().filter(|item| item.family == family && Some(item.side) == winner).map(|item| item.ordinal));
        }
        conflicts.push(Conflict { family, ps3: ps3[family], bs3: bs3[family], winner,
            ps3_ordinals: evidence.iter().filter(|item| item.family == family && item.side == Side::Ps3).map(|item| item.ordinal).collect(),
            bs3_ordinals: evidence.iter().filter(|item| item.family == family && item.side == Side::Bs3).map(|item| item.ordinal).collect(),
        });
    }
    let families = |ordinals: &[usize]| {
        let mut result = Vec::new();
        for &ordinal in ordinals {
            let family = evidence[ordinal].family;
            if !result.contains(&family) { result.push(family); }
        }
        result
    };
    Ok(Resolution { conflicts, ps3_families: families(&final_ps3), bs3_families: families(&final_bs3),
        final_ps3, final_bs3, clean_ps3, clean_bs3, publications: ["MaveDB_PS3", "MaveDB_BS3", "MaveDB_score_interpretation"] })
}

pub fn format_evidence(evidence: &[Evidence], ordinals: &[usize]) -> Result<String, String> {
    let mut tokens = Vec::with_capacity(ordinals.len());
    for &ordinal in ordinals {
        let item = evidence.get(ordinal).ok_or("invalid evidence formatting ordinal")?;
        tokens.push(item.formatted.as_str());
    }
    Ok(tokens.join("; "))
}
