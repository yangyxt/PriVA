//! Source-faithful PM4/BP3 dataframe procedures. Interval intersection stays in
//! the mature bedtools/pybedtools caller; no genomic intersection is reimplemented.
pub(crate) mod transport;

use polars_core::prelude::*;
use std::collections::HashSet;

pub type Result<T> = std::result::Result<T, String>;

fn error(value: impl std::fmt::Display) -> String {
    value.to_string()
}

fn text(frame: &DataFrame, name: &str) -> Result<Vec<String>> {
    frame
        .column(name)
        .map_err(error)?
        .str()
        .map_err(error)?
        .into_iter()
        .map(|value| {
            value
                .map(str::to_owned)
                .ok_or_else(|| format!("null {name}"))
        })
        .collect()
}

fn integer(frame: &DataFrame, name: &str) -> Result<Vec<i64>> {
    frame
        .column(name)
        .map_err(error)?
        .i64()
        .map_err(error)?
        .into_iter()
        .map(|value| value.ok_or_else(|| format!("null {name}")))
        .collect()
}

fn boolean(frame: &DataFrame, name: &str) -> Result<Vec<bool>> {
    frame
        .column(name)
        .map_err(error)?
        .bool()
        .map_err(error)?
        .into_iter()
        .map(|value| value.ok_or_else(|| format!("null {name}")))
        .collect()
}

fn bcolumn(name: &str, values: Vec<bool>) -> Column {
    Series::new(name.into(), values).into()
}
fn icolumn(name: &str, values: Vec<i64>) -> Column {
    Series::new(name.into(), values).into()
}
fn scolumn(name: &str, values: Vec<String>) -> Column {
    Series::new(name.into(), values).into()
}

/// The original BED uses the complete REF length, including Unicode codepoints,
/// with start=pos-1 and end=pos-1+len(ref). Negative/zero intervals are left for
/// actual bedtools to accept/reject. Overflow is rejected before publication.
pub fn prepare_bed(frame: &DataFrame) -> Result<DataFrame> {
    let positions = integer(frame, "pos")?;
    let refs = text(frame, "ref")?;
    let mut starts = Vec::with_capacity(frame.height());
    let mut ends = Vec::with_capacity(frame.height());
    for (position, reference) in positions.into_iter().zip(refs) {
        let start = position.checked_sub(1).ok_or("BED start overflow")?;
        let length = i64::try_from(reference.chars().count()).map_err(error)?;
        let end = start.checked_add(length).ok_or("BED end overflow")?;
        starts.push(start);
        ends.push(end);
    }
    DataFrame::new(vec![
        icolumn("ordinal", integer(frame, "ordinal")?),
        scolumn("chrom", text(frame, "chrom")?),
        icolumn("start", starts),
        icolumn("end", ends),
        scolumn("variant_id", text(frame, "variant_id")?),
    ])
    .map_err(error)
}

/// Preserve all source masks and the two final source assignments. In particular,
/// PM4 eligibility suppresses BP3 BEFORE the final PVS1<2 assignment; this is not
/// interchangeable with suppressing BP3 only when a PM4 score was assigned.
pub fn assign(frame: &DataFrame, overlapping_ids: &HashSet<String>) -> Result<DataFrame> {
    let ids = text(frame, "variant_id")?;
    let domains = text(frame, "DOMAINS")?;
    let consequence = text(frame, "Consequence")?;
    let nmd = text(frame, "NMD")?;
    let lof = text(frame, "LoF_filter")?;
    let splice_len = boolean(frame, "splicing_len_changing")?;
    let splice_frame = boolean(frame, "splicing_frameshift")?;
    let utr_len = boolean(frame, "5UTR_len_changing")?;
    let utr_frame = boolean(frame, "5UTR_frameshift")?;
    let vep_len = boolean(frame, "vep_consq_len_changing")?;
    let splice_domain = boolean(frame, "splicing_span_intol_domain")?;
    let utr_domain = boolean(frame, "5UTR_span_intol_domain")?;
    let span_domain = boolean(frame, "span_functional_domains")?;
    let pvs = integer(frame, "pvs1")?;
    let intolerant = boolean(frame, "loc_intol_domain")?;
    let mut repeated = Vec::new();
    let mut low_complexity = Vec::new();
    let mut repetitive = Vec::new();
    let mut inframe_indels = Vec::new();
    let mut splicing_inframe = Vec::new();
    let mut utr_inframe = Vec::new();
    let mut vep_inframe = Vec::new();
    let mut all_inframe = Vec::new();
    let mut stop_loss = Vec::new();
    let mut escaping_trunc = Vec::new();
    let mut functional = Vec::new();
    let mut frameshift_pvs = Vec::new();
    let mut strong_nmd = Vec::new();
    let mut excluded = Vec::new();
    let mut pm4_inframe = Vec::new();
    let mut pm4_eligible = Vec::new();
    let mut bp3_eligible = Vec::new();
    let mut pm4 = Vec::new();
    let mut bp3 = Vec::new();
    for i in 0..frame.height() {
        // contains() deliberately preserves the source's unanchored,
        // case-sensitive literal substrings, including misleading prefixes.
        let r = overlapping_ids.contains(&ids[i]);
        let l = domains[i].contains("Low_complexity");
        let repetitive_row = r | l;
        let indel = consequence[i].contains("inframe_deletion")
            | consequence[i].contains("inframe_insertion");
        let splice = splice_len[i] & !splice_frame[i];
        let utr = utr_len[i] & !utr_frame[i];
        let vep = vep_len[i]
            & !consequence[i].contains("frameshift")
            & !consequence[i].contains("stop_gained");
        let inframe = indel | splice | utr | vep;
        let stop = consequence[i].contains("stop_lost");
        let escaping = (consequence[i].contains("frameshift")
            | consequence[i].contains("stop_gained"))
            & (nmd[i].contains("escaping") | lof[i].contains("END_TRUNC"));
        let domain = intolerant[i] | splice_domain[i] | utr_domain[i] | span_domain[i];
        let fs_pvs = consequence[i].contains("frameshift") & (pvs[i] > 0) & !stop;
        let strong = (pvs[i] >= 3) & !stop;
        let exclude = fs_pvs | strong;
        let pm4_if = inframe & !exclude & (!repetitive_row | domain);
        let eligible = pm4_if | stop | escaping;
        let bp3_if = inframe & repetitive_row & !domain;
        let bp3_final_eligible = bp3_if & !eligible;
        let pm4_score = if eligible & (pvs[i] < 2) { 2 } else { 0 };
        let bp3_score = if bp3_final_eligible & (pm4_score == 0) {
            1
        } else {
            0
        };
        repeated.push(r);
        low_complexity.push(l);
        repetitive.push(repetitive_row);
        inframe_indels.push(indel);
        splicing_inframe.push(splice);
        utr_inframe.push(utr);
        vep_inframe.push(vep);
        all_inframe.push(inframe);
        stop_loss.push(stop);
        escaping_trunc.push(escaping);
        functional.push(domain);
        frameshift_pvs.push(fs_pvs);
        strong_nmd.push(strong);
        excluded.push(exclude);
        pm4_inframe.push(pm4_if);
        pm4_eligible.push(eligible);
        bp3_eligible.push(bp3_final_eligible);
        pm4.push(pm4_score);
        bp3.push(bp3_score);
    }
    DataFrame::new(vec![
        icolumn("ordinal", integer(frame, "ordinal")?),
        bcolumn("in_repeat", repeated),
        bcolumn("low_complexity_domain", low_complexity),
        bcolumn("in_repetitive_region", repetitive),
        bcolumn("inframe_indels", inframe_indels),
        bcolumn("splicing_inframe", splicing_inframe),
        bcolumn("utr_inframe", utr_inframe),
        bcolumn("vep_inframe", vep_inframe),
        bcolumn("all_inframe", all_inframe),
        bcolumn("stop_loss_variants", stop_loss),
        bcolumn("nmd_escaping_trunc", escaping_trunc),
        bcolumn("has_functional_domain", functional),
        bcolumn("frameshift_with_pvs1", frameshift_pvs),
        bcolumn("strong_nmd", strong_nmd),
        bcolumn("exclude_from_pm4", excluded),
        bcolumn("pm4_inframe_eligible", pm4_inframe),
        bcolumn("pm4_eligible", pm4_eligible),
        bcolumn("bp3_eligible", bp3_eligible),
        icolumn("pm4_array", pm4),
        icolumn("bp3_array", bp3),
    ])
    .map_err(error)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(consequence: &[&str], pvs: Vec<i64>) -> DataFrame {
        let n = consequence.len();
        let mut columns = vec![
            icolumn("ordinal", (0..i64::try_from(n).unwrap()).collect()),
            scolumn("variant_id", vec!["shared".into(); n]),
            scolumn("DOMAINS", vec![".".into(); n]),
            scolumn(
                "Consequence",
                consequence.iter().map(|s| (*s).into()).collect(),
            ),
            scolumn("NMD", vec!["escaping".into(); n]),
            scolumn("LoF_filter", vec![".".into(); n]),
            icolumn("pvs1", pvs),
            bcolumn("loc_intol_domain", vec![false; n]),
        ];
        for name in [
            "splicing_len_changing",
            "splicing_frameshift",
            "5UTR_len_changing",
            "5UTR_frameshift",
            "vep_consq_len_changing",
            "splicing_span_intol_domain",
            "5UTR_span_intol_domain",
            "span_functional_domains",
        ] {
            columns.push(bcolumn(name, vec![false; n]));
        }
        DataFrame::new(columns).unwrap()
    }

    #[test]
    fn full_reference_and_unicode_lengths() {
        let input = DataFrame::new(vec![
            icolumn("ordinal", vec![0, 1]),
            scolumn("chrom", vec!["chr1".into(); 2]),
            icolumn("pos", vec![11, 20]),
            scolumn("variant_id", vec!["a".into(), "a".into()]),
            scolumn("ref", vec!["😀AA".into(), "A\u{301}".into()]),
        ])
        .unwrap();
        let got = prepare_bed(&input).unwrap();
        assert_eq!(integer(&got, "start").unwrap(), vec![10, 19]);
        assert_eq!(integer(&got, "end").unwrap(), vec![13, 21]);
    }

    #[test]
    fn eligibility_suppresses_bp3_before_pvs1_gate() {
        let input = fixture(
            &["inframe_deletion&stop_lost", "inframe_deletion&stop_gained"],
            vec![2, 4],
        );
        let got = assign(&input, &HashSet::from(["shared".into()])).unwrap();
        assert_eq!(integer(&got, "pm4_array").unwrap(), vec![0, 0]);
        assert_eq!(integer(&got, "bp3_array").unwrap(), vec![0, 0]);
        assert_eq!(boolean(&got, "pm4_eligible").unwrap(), vec![true, true]);
    }

    #[test]
    fn shared_string_id_has_global_repeat_membership() {
        let input = fixture(&["inframe_deletion", "inframe_insertion"], vec![0, 0]);
        let got = assign(&input, &HashSet::from(["shared".into()])).unwrap();
        assert_eq!(integer(&got, "pm4_array").unwrap(), vec![0, 0]);
        assert_eq!(integer(&got, "bp3_array").unwrap(), vec![1, 1]);
    }

    #[test]
    fn source_substrings_are_not_corrected() {
        let input = fixture(
            &[
                "not_inframe_deletion",
                "not_stop_lost",
                "frameshift_variant",
            ],
            vec![0, 0, 0],
        );
        let got = assign(&input, &HashSet::new()).unwrap();
        assert_eq!(integer(&got, "pm4_array").unwrap(), vec![2, 2, 2]);
    }

    #[test]
    fn bed_integer_overflow_is_rejected() {
        let input = DataFrame::new(vec![
            icolumn("ordinal", vec![0]),
            scolumn("chrom", vec!["chr1".into()]),
            icolumn("pos", vec![i64::MIN]),
            scolumn("variant_id", vec!["a".into()]),
            scolumn("ref", vec!["A".into()]),
        ])
        .unwrap();
        assert!(prepare_bed(&input).is_err());
    }
}
