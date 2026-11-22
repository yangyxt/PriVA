//! Candidate decisions from the approved domain-sorted-rescue-v1 reference.
//!
//! The caller resolves the FIRST InterPro entry in the source's entry order.
//! GO molecular-function membership includes tuple GO terms; the existing
//! runtime functionality interpreter has a different list-only rule.
//! These procedures do not sort domains, choose task seeds, or activate DAS.
use crate::exon_af_stats;

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Entry {
    pub(crate) kind: String,
    pub(crate) has_go_molecular_function: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Action {
    Skip,
    Direct,
    /// The caller uses this type's composite, falling back to the global one.
    Test(String),
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Selection {
    /// Source counter before the median threshold: n_go_mf / n_eligible_type.
    pub(crate) eligible: bool,
    pub(crate) action: Action,
}

fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        // NumPy's empty median is NaN; NaN < the threshold is false.
        f64::NAN
    } else {
        exon_af_stats::calculate(values)[0]
    }
}

fn skipped(eligible: bool) -> Selection {
    Selection {
        eligible,
        action: Action::Skip,
    }
}

/// The reference's GO-MF branch increments its counter before the median test.
pub(crate) fn go_mf(
    scores: Option<&[f64]>,
    first_entry: Option<&Entry>,
    minimum_median: f64,
) -> Selection {
    let (Some(scores), Some(entry)) = (scores, first_entry) else {
        return skipped(false);
    };
    if !entry.has_go_molecular_function {
        return skipped(false);
    }
    if median(scores) < minimum_median {
        return skipped(true);
    }
    Selection {
        eligible: true,
        action: if [
            "Domain",
            "Conserved_site",
            "Binding_site",
            "Active_site",
            "PTM",
        ]
        .contains(&entry.kind.as_str())
        {
            Action::Direct
        } else {
            Action::Test(entry.kind.clone())
        },
    }
}

/// Preserve the implemented Homologous_superfamily inclusion even though the
/// original docstring says it is excluded. Orphan eligibility uses field 1 of
/// the full colon-separated domain path, only when at least 3 fields exist.
pub(crate) fn non_go_mf(
    domain_path: &str,
    scores: Option<&[f64]>,
    first_entry: Option<&Entry>,
    minimum_median: f64,
) -> Selection {
    let Some(scores) = scores else {
        return skipped(false);
    };
    let composite = if let Some(entry) = first_entry {
        if ![
            "Domain",
            "Homologous_superfamily",
            "Repeat",
            "Conserved_site",
            "Binding_site",
            "Active_site",
            "PTM",
        ]
        .contains(&entry.kind.as_str())
            || entry.has_go_molecular_function
        {
            return skipped(false);
        }
        entry.kind.clone()
    } else {
        let parts: Vec<_> = domain_path.split(':').collect();
        let source = if parts.len() >= 3 { parts[1] } else { "" };
        if ![
            "Gene3D",
            "Superfamily",
            "PANTHER",
            "Pfam",
            "PIRSF",
            "SMART",
            "PROSITE_profiles",
            "Prints",
        ]
        .contains(&source)
            || [
                "Coiled-coils_(Ncoils)",
                "Low_complexity_(Seg)",
                "Transmembrane_helices",
                "Cleavage_site_(Signalp)",
                "PROSITE_patterns",
            ]
            .contains(&source)
        {
            return skipped(false);
        }
        "Domain".to_owned()
    };
    if median(scores) < minimum_median {
        return skipped(true);
    }
    Selection {
        eligible: true,
        action: Action::Test(composite),
    }
}
