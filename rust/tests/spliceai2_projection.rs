use priva_tabular::spliceai2::{Statistics, VariantKind, native_columns, project};
use std::io::{Cursor, Write};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};

const CONTEXT_HEADER: &str = "gene_id\tchrom\tgene_start\tgene_end\tstrand\tgene_name\n";
const DEFAULT_GENES: &str =
    "ENSG001.1\tchr1\t10\t100\t+\tPLUS\nENSG002.2\tchr1\t10\t100\t-\tMINUS\n";

fn genes(body: &str) -> String {
    format!("{CONTEXT_HEADER}{body}")
}

fn row() -> Vec<String> {
    native_columns()
        .iter()
        .map(|column| {
            match column.as_str() {
                "chrom" => "chr1",
                "pos" => "20",
                "ref" => "A",
                "alt" => "C",
                "strand" => "+",
                "spliceai2_summary_score" => "0.40",
                "donor_gain_delta_score_0" => "0.30",
                "donor_loss_delta_score_0" => "0.40",
                "acceptor_gain_delta_score_0" => "0.10",
                "acceptor_loss_delta_score_0" => "0.20",
                "donor_gain_dist_0" => "-11",
                "donor_loss_dist_0" => "12",
                "acceptor_gain_dist_0" => "-13",
                "acceptor_loss_dist_0" => "14",
                name if name.contains("_dist_") => "0",
                _ => "0.0",
            }
            .to_owned()
        })
        .collect()
}

fn set(row: &mut [String], column: &str, value: &str) {
    let index = native_columns()
        .iter()
        .position(|name| name == column)
        .unwrap();
    row[index] = value.to_owned();
}

fn table(rows: &[Vec<String>]) -> String {
    let mut text = native_columns().join("\t") + "\n";
    for row in rows {
        text.push_str(&row.join("\t"));
        text.push('\n');
    }
    text
}

fn run(
    context: &str,
    text: &str,
    kind: VariantKind,
    header: bool,
) -> Result<(String, Statistics), String> {
    let mut output = Vec::new();
    let stats = project(
        Cursor::new(context),
        Cursor::new(text),
        &mut output,
        "chr1",
        kind,
        header,
    )
    .map_err(|error| error.to_string())?;
    Ok((String::from_utf8(output).unwrap(), stats))
}

fn rejects(row: Vec<String>, message: &str) {
    let error = run(
        &genes(DEFAULT_GENES),
        &table(&[row]),
        VariantKind::Snv,
        false,
    )
    .unwrap_err();
    assert!(
        error.contains(message),
        "{error:?} does not contain {message:?}"
    );
}

#[test]
fn published_header_exact_order() {
    // Actual top-three source header: 84 fields, SHA256 including newline
    // 1feb799bed12dae1030e64b03930eab470340848aa4a7e0166072a140e7b0d66.
    let expected = "chrom,pos,ref,alt,strand,donor_gain_delta_score_0,donor_gain_delta_score_1,donor_gain_delta_score_2,donor_gain_ref_score_0,donor_gain_ref_score_1,donor_gain_ref_score_2,donor_gain_alt_score_0,donor_gain_alt_score_1,donor_gain_alt_score_2,donor_gain_dist_0,donor_gain_dist_1,donor_gain_dist_2,donor_loss_delta_score_0,donor_loss_delta_score_1,donor_loss_delta_score_2,donor_loss_ref_score_0,donor_loss_ref_score_1,donor_loss_ref_score_2,donor_loss_alt_score_0,donor_loss_alt_score_1,donor_loss_alt_score_2,donor_loss_dist_0,donor_loss_dist_1,donor_loss_dist_2,acceptor_gain_delta_score_0,acceptor_gain_delta_score_1,acceptor_gain_delta_score_2,acceptor_gain_ref_score_0,acceptor_gain_ref_score_1,acceptor_gain_ref_score_2,acceptor_gain_alt_score_0,acceptor_gain_alt_score_1,acceptor_gain_alt_score_2,acceptor_gain_dist_0,acceptor_gain_dist_1,acceptor_gain_dist_2,acceptor_loss_delta_score_0,acceptor_loss_delta_score_1,acceptor_loss_delta_score_2,acceptor_loss_ref_score_0,acceptor_loss_ref_score_1,acceptor_loss_ref_score_2,acceptor_loss_alt_score_0,acceptor_loss_alt_score_1,acceptor_loss_alt_score_2,acceptor_loss_dist_0,acceptor_loss_dist_1,acceptor_loss_dist_2,jxn_gain_delta_score_0,jxn_gain_delta_score_1,jxn_gain_delta_score_2,jxn_gain_ref_score_0,jxn_gain_ref_score_1,jxn_gain_ref_score_2,jxn_gain_alt_score_0,jxn_gain_alt_score_1,jxn_gain_alt_score_2,jxn_gain_donor_dist_0,jxn_gain_donor_dist_1,jxn_gain_donor_dist_2,jxn_gain_acceptor_dist_0,jxn_gain_acceptor_dist_1,jxn_gain_acceptor_dist_2,jxn_loss_delta_score_0,jxn_loss_delta_score_1,jxn_loss_delta_score_2,jxn_loss_ref_score_0,jxn_loss_ref_score_1,jxn_loss_ref_score_2,jxn_loss_alt_score_0,jxn_loss_alt_score_1,jxn_loss_alt_score_2,jxn_loss_donor_dist_0,jxn_loss_donor_dist_1,jxn_loss_donor_dist_2,jxn_loss_acceptor_dist_0,jxn_loss_acceptor_dist_1,jxn_loss_acceptor_dist_2,spliceai2_summary_score";
    assert_eq!(native_columns().join(","), expected);
    assert_eq!(native_columns().len(), 84);
}

#[test]
fn exact_projection_and_only_info_field() {
    let (out, stats) = run(
        &genes(DEFAULT_GENES),
        &table(&[row()]),
        VariantKind::Snv,
        false,
    )
    .unwrap();
    assert_eq!(
        out,
        "chr1\t20\t.\tA\tC\t.\t.\tSpliceAI=C|PLUS|0.10|0.20|0.30|0.40|-13|14|-11|12\n"
    );
    assert_eq!(
        stats,
        Statistics {
            input_rows: 1,
            output_records: 1
        }
    );
    let info = out.trim_end().split('\t').nth(7).unwrap();
    assert_eq!(info.split('|').count(), 10);
    assert!(!info.contains(';') && !info.contains(','));
}

#[test]
fn minus_strand_genomic_distances_are_not_reflected() {
    let mut r = row();
    set(&mut r, "strand", "-");
    let (out, _) = run(&genes(DEFAULT_GENES), &table(&[r]), VariantKind::Snv, false).unwrap();
    assert!(out.ends_with("C|MINUS|0.10|0.20|0.30|0.40|-13|14|-11|12\n"));
}

#[test]
fn overlapping_genes_emit_separate_records_in_deterministic_symbol_order() {
    let context = genes(
        "ENSG003\tchr1\t10\t100\t+\tZETA\nENSG001\tchr1\t10\t100\t+\tALPHA\nENSG002\tchr1\t10\t100\t-\tOPPOSITE\n",
    );
    let (out, stats) = run(&context, &table(&[row()]), VariantKind::Snv, false).unwrap();
    let lines: Vec<_> = out.lines().collect();
    assert_eq!(
        stats,
        Statistics {
            input_rows: 1,
            output_records: 2
        }
    );
    assert!(lines[0].contains("|ALPHA|"));
    assert!(lines[1].contains("|ZETA|"));
    assert!(!out.contains("OPPOSITE"));
}

#[test]
fn gene_interval_bounds_are_inclusive_and_expired_intervals_removed() {
    let context = genes("G1\tchr1\t10\t20\t+\tFIRST\nG2\tchr1\t21\t30\t+\tSECOND\n");
    let rows: Vec<_> = [10, 20, 21, 30]
        .into_iter()
        .map(|pos| {
            let mut r = row();
            set(&mut r, "pos", &pos.to_string());
            r
        })
        .collect();
    let (out, stats) = run(&context, &table(&rows), VariantKind::Snv, false).unwrap();
    assert_eq!(stats.output_records, 4);
    assert_eq!(out.matches("|FIRST|").count(), 2);
    assert_eq!(out.matches("|SECOND|").count(), 2);
}

#[test]
fn no_matching_gene_or_opposite_strand_only_fails() {
    for context in ["G1\tchr1\t1\t19\t+\tOTHER\n", "G1\tchr1\t1\t30\t-\tOTHER\n"] {
        assert!(
            run(&genes(context), &table(&[row()]), VariantKind::Snv, false)
                .unwrap_err()
                .contains("no verified gene interval")
        );
    }
}

#[test]
fn same_symbol_same_strand_overlap_keeps_all_intervals_and_emits_once() {
    let context = genes("G1\tchr1\t1\t30\t+\tSAME\nG2\tchr1\t15\t50\t+\tSAME\n");
    let rows: Vec<_> = [5, 20, 40]
        .into_iter()
        .map(|pos| {
            let mut r = row();
            set(&mut r, "pos", &pos.to_string());
            r
        })
        .collect();
    let (out, stats) = run(&context, &table(&rows), VariantKind::Snv, false).unwrap();
    assert_eq!(
        stats,
        Statistics {
            input_rows: 3,
            output_records: 3
        }
    );
    assert_eq!(out.matches("|SAME|").count(), 3);
    assert_eq!(
        out.lines()
            .map(|line| line.split('\t').nth(1).unwrap())
            .collect::<Vec<_>>(),
        ["5", "20", "40"]
    );
}

#[test]
fn same_chromosome_symbol_opposite_strands_rejected_even_without_overlap() {
    for second in ["G2\tchr1\t1\t30\t-\tSAME\n", "G2\tchr1\t50\t90\t-\tSAME\n"] {
        let context = genes(&("G1\tchr1\t1\t30\t+\tSAME\n".to_owned() + second));
        assert!(
            run(&context, &table(&[row()]), VariantKind::Snv, false)
                .unwrap_err()
                .contains("ambiguous symbol SAME on chr1 has opposite-strand context IDs")
        );
    }
    let context =
        genes("G1\tchr1\t1\t30\t+\tSAME\nG2\tchrX\t50\t90\t-\tSAME\nG3\tchrY\t50\t90\t+\tSAME\n");
    assert!(run(&context, &table(&[row()]), VariantKind::Snv, false).is_ok());
}

#[test]
fn same_symbol_same_strand_nonoverlapping_intervals_are_retained() {
    let context = genes("G1\tchr1\t1\t30\t+\tSAME\nG2\tchr1\t50\t90\t+\tSAME\n");
    let mut next = row();
    set(&mut next, "pos", "60");
    let (out, stats) = run(&context, &table(&[row(), next]), VariantKind::Snv, false).unwrap();
    assert_eq!(
        stats,
        Statistics {
            input_rows: 2,
            output_records: 2
        }
    );
    assert_eq!(out.matches("|SAME|").count(), 2);
    let mut gap = row();
    set(&mut gap, "pos", "40");
    assert!(
        run(&context, &table(&[gap]), VariantKind::Snv, false)
            .unwrap_err()
            .contains("no verified gene interval")
    );
}

#[test]
fn coalescing_context_labels_preserves_prediction_rows_order_and_duplicate_guard() {
    let context =
        genes("G1\tchr1\t1\t30\t+\tZETA\nG2\tchr1\t15\t50\t+\tZETA\nG3\tchr1\t1\t30\t+\tALPHA\n");
    let mut alternate = row();
    set(&mut alternate, "alt", "G");
    let (out, stats) = run(
        &context,
        &table(&[row(), alternate]),
        VariantKind::Snv,
        false,
    )
    .unwrap();
    assert_eq!(
        stats,
        Statistics {
            input_rows: 2,
            output_records: 4
        }
    );
    assert_eq!(
        out.lines()
            .map(|line| line
                .split('\t')
                .nth(7)
                .unwrap()
                .split('|')
                .take(2)
                .collect::<Vec<_>>())
            .collect::<Vec<_>>(),
        [
            vec!["SpliceAI=C", "ALPHA"],
            vec!["SpliceAI=C", "ZETA"],
            vec!["SpliceAI=G", "ALPHA"],
            vec!["SpliceAI=G", "ZETA"]
        ]
    );
    assert!(
        run(&context, &table(&[row(), row()]), VariantKind::Snv, false)
            .unwrap_err()
            .contains("duplicate allele + model strand")
    );
}

#[test]
fn duplicate_context_ids_fail_even_if_rows_identical() {
    for second in ["G1\tchr1\t1\t30\t+\tSAME\n", "G1\tchr2\t1\t40\t-\tOTHER\n"] {
        let context = genes(&("G1\tchr1\t1\t30\t+\tSAME\n".to_owned() + second));
        assert!(
            run(&context, &table(&[row()]), VariantKind::Snv, false)
                .unwrap_err()
                .contains("duplicate gene context ID")
        );
    }
}

#[test]
fn malformed_context_rejected() {
    for body in [
        "G1\tchr1\t0\t30\t+\tGENE\n",
        "G1\tchr1\t30\t1\t+\tGENE\n",
        "G1\tchr1\t1\t30\t.\tGENE\n",
        "G1\tchr1\t1\t30\t+\tBAD;GENE\n",
        "G1\tchr01\t1\t30\t+\tGENE\n",
        "G1\tchr1\t1\t30\t+\n",
        "\n",
    ] {
        assert!(
            run(&genes(body), &table(&[row()]), VariantKind::Snv, false).is_err(),
            "accepted {body:?}"
        );
    }
    assert!(run(DEFAULT_GENES, &table(&[row()]), VariantKind::Snv, false).is_err());
}

#[test]
fn no_score_filtering_at_any_recommended_threshold() {
    for score in [
        "-0.0", "0.00", "0.01", "0.099", "0.1", "0.1001", "0.2499", "0.25", "0.2501", "0.4999",
        "0.5", "0.5001", "1.0",
    ] {
        let mut r = row();
        for column in [
            "acceptor_gain_delta_score_0",
            "acceptor_loss_delta_score_0",
            "donor_gain_delta_score_0",
            "donor_loss_delta_score_0",
            "spliceai2_summary_score",
        ] {
            set(&mut r, column, score);
        }
        let (out, stats) =
            run(&genes(DEFAULT_GENES), &table(&[r]), VariantKind::Snv, false).unwrap();
        assert_eq!(stats.output_records, 1);
        assert!(out.contains(&format!("|{score}|{score}|{score}|{score}|")));
    }
}

#[test]
fn signed_deltas_preserved_without_abs_clamping_or_ref_alt_recalculation() {
    let mut r = row();
    for column in [
        "acceptor_gain_delta_score_0",
        "acceptor_loss_delta_score_0",
        "donor_gain_delta_score_0",
        "donor_loss_delta_score_0",
        "spliceai2_summary_score",
    ] {
        set(&mut r, column, "-0.25");
    }
    set(&mut r, "jxn_loss_delta_score_2", "-1.0");
    // Independent rounded ref/alt values deliberately do not reconstruct delta.
    set(&mut r, "acceptor_gain_ref_score_0", "0.0");
    set(&mut r, "acceptor_gain_alt_score_0", "0.0");
    let (out, _) = run(&genes(DEFAULT_GENES), &table(&[r]), VariantKind::Snv, false).unwrap();
    assert!(out.contains("|-0.25|-0.25|-0.25|-0.25|"));
}

#[test]
fn summary_uses_top_zero_sites_only_and_ties_keep_original_offsets() {
    let mut r = row();
    set(&mut r, "acceptor_gain_delta_score_0", "0.40");
    set(&mut r, "jxn_gain_delta_score_0", "0.99");
    set(&mut r, "donor_loss_dist_1", "32767");
    let (out, _) = run(&genes(DEFAULT_GENES), &table(&[r]), VariantKind::Snv, false).unwrap();
    assert!(out.ends_with("|0.40|0.20|0.30|0.40|-13|14|-11|12\n"));
    let mut bad = row();
    set(&mut bad, "spliceai2_summary_score", "0.4000001");
    rejects(bad, "summary_score differs");
}

#[test]
fn validates_all_probability_and_signed_delta_fields() {
    let cases = [
        ("acceptor_gain_delta_score_2", "NaN"),
        ("donor_loss_ref_score_2", "-0.1"),
        ("jxn_gain_alt_score_1", "1.01"),
        ("jxn_loss_delta_score_2", "-1.01"),
        ("jxn_loss_ref_score_0", "inf"),
        ("donor_gain_delta_score_1", ""),
        ("donor_gain_delta_score_1", "NA"),
        ("acceptor_loss_alt_score_2", " 0.1"),
        ("spliceai2_summary_score", "-inf"),
    ];
    for (name, value) in cases {
        let mut r = row();
        set(&mut r, name, value);
        assert!(
            run(&genes(DEFAULT_GENES), &table(&[r]), VariantKind::Snv, false).is_err(),
            "accepted {name}={value:?}"
        );
    }
}

#[test]
fn distance_integer_lexemes_preserved_and_unsupported_forms_rejected() {
    let mut r = row();
    set(&mut r, "acceptor_gain_dist_0", "-32768");
    set(&mut r, "acceptor_loss_dist_0", "32767");
    set(&mut r, "donor_gain_dist_0", "-0");
    set(&mut r, "donor_loss_dist_0", "00012");
    let (out, _) = run(&genes(DEFAULT_GENES), &table(&[r]), VariantKind::Snv, false).unwrap();
    assert!(out.ends_with("|-32768|32767|-0|00012\n"));
    for value in [
        "-0.0",
        "1.0",
        "1e3",
        "+2",
        "NaN",
        "",
        "2147483648",
        "-2147483649",
    ] {
        let mut r = row();
        set(&mut r, "jxn_gain_acceptor_dist_2", value);
        rejects(r, "distance");
    }
}

#[test]
fn input_must_have_exact_native_header() {
    let text = table(&[row()]);
    for bad in [
        text.replacen("chrom\tpos", "pos\tchrom", 1),
        text.replacen("donor_gain_delta_score_0", "unknown", 1),
        text.replacen("chrom", "#chrom", 1),
        row().join("\t") + "\n",
        String::new(),
    ] {
        assert!(
            run(&genes(DEFAULT_GENES), &bad, VariantKind::Snv, false)
                .unwrap_err()
                .contains("exact native 84-column")
        );
    }
}

#[test]
fn missing_extra_and_blank_rows_fail() {
    let mut missing = row();
    missing.pop();
    rejects(missing, "expected 84 columns");
    let mut extra = row();
    extra.push("x".into());
    rejects(extra, "expected 84 columns");
    assert!(
        run(
            &genes(DEFAULT_GENES),
            &(table(&[]) + "\n"),
            VariantKind::Snv,
            false
        )
        .is_err()
    );
}

#[test]
fn invalid_coordinates_chromosomes_strands_and_alleles_fail() {
    for (column, value) in [
        ("pos", "0"),
        ("pos", "020"),
        ("pos", "-2"),
        ("pos", "2.0"),
        ("chrom", "1"),
        ("chrom", "chr2"),
        ("strand", "1"),
        ("strand", "."),
        ("ref", "a"),
        ("ref", "N"),
        ("alt", "*"),
        ("alt", "<DEL>"),
        ("alt", "C,G"),
        ("alt", "A"),
    ] {
        let mut r = row();
        set(&mut r, column, value);
        assert!(
            run(&genes(DEFAULT_GENES), &table(&[r]), VariantKind::Snv, false).is_err(),
            "accepted {column}={value}"
        );
    }
}

#[test]
fn duplicate_variant_model_strand_fails_but_distinct_alt_and_strands_pass() {
    assert!(
        run(
            &genes(DEFAULT_GENES),
            &table(&[row(), row()]),
            VariantKind::Snv,
            false
        )
        .unwrap_err()
        .contains("duplicate allele + model strand")
    );
    let mut different = row();
    set(&mut different, "alt", "G");
    let mut other_strand = row();
    set(&mut other_strand, "strand", "-");
    let (_, stats) = run(
        &genes(DEFAULT_GENES),
        &table(&[row(), different, other_strand]),
        VariantKind::Snv,
        false,
    )
    .unwrap();
    assert_eq!(stats.input_rows, 3);
}

#[test]
fn duplicate_tracking_resets_only_after_coordinate_advance() {
    let mut next = row();
    set(&mut next, "pos", "21");
    assert!(
        run(
            &genes(DEFAULT_GENES),
            &table(&[row(), next.clone()]),
            VariantKind::Snv,
            false
        )
        .is_ok()
    );
    assert!(
        run(
            &genes(DEFAULT_GENES),
            &table(&[next, row()]),
            VariantKind::Snv,
            false
        )
        .unwrap_err()
        .contains("not ascending")
    );
}

#[test]
fn insertions_deletions_and_complex_indels_remain_exactly_anchored() {
    for (reference, alternate) in [("A", "AT"), ("AT", "A"), ("AC", "AGT"), ("TC", "A")] {
        let mut r = row();
        set(&mut r, "ref", reference);
        set(&mut r, "alt", alternate);
        let (out, _) = run(
            &genes(DEFAULT_GENES),
            &table(&[r]),
            VariantKind::Indel,
            false,
        )
        .unwrap();
        assert!(out.starts_with(&format!(
            "chr1\t20\t.\t{reference}\t{alternate}\t.\t.\tSpliceAI={alternate}|"
        )));
    }
}

#[test]
fn variant_kind_mismatch_rejected() {
    let mut r = row();
    set(&mut r, "alt", "CT");
    rejects(r, "snv input requires");
    assert!(
        run(
            &genes(DEFAULT_GENES),
            &table(&[row()]),
            VariantKind::Indel,
            false
        )
        .unwrap_err()
        .contains("indel input requires")
    );
    let mut r = row();
    set(&mut r, "ref", "AC");
    set(&mut r, "alt", "GT");
    assert!(
        run(
            &genes(DEFAULT_GENES),
            &table(&[r]),
            VariantKind::Indel,
            false
        )
        .is_err()
    );
}

#[test]
fn output_header_is_complete_and_marks_lossy_annotation_only_view() {
    let (out, _) = run(
        &genes(DEFAULT_GENES),
        &table(&[row()]),
        VariantKind::Snv,
        true,
    )
    .unwrap();
    assert_eq!(
        out.lines()
            .filter(|line| line.starts_with("##contig="))
            .count(),
        24
    );
    for text in [
        "##reference=GRCh38",
        "##SpliceAI2Version=2",
        "##SpliceAI2ProjectionLossy=true",
        "##SpliceAI2RecommendedThresholds=0.1,0.25,0.5",
        "##SpliceAI2Policy=annotation_only_no_ACMG_or_CADD_activation",
        "##SpliceAI2PublisherGeneSpecificPrediction=false",
        "##contig=<ID=chrX>",
        "##contig=<ID=chrY>",
    ] {
        assert!(out.contains(text));
    }
    assert_eq!(
        out.lines()
            .filter(|line| line.starts_with("##INFO="))
            .count(),
        1
    );
    assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"));
    let (body, _) = run(
        &genes(DEFAULT_GENES),
        &table(&[row()]),
        VariantKind::Snv,
        false,
    )
    .unwrap();
    assert!(!body.contains('#'));
    assert!(out.ends_with(&body));
}

#[test]
fn crlf_and_unterminated_final_row_do_not_change_field_strings() {
    let text = table(&[row()]);
    let (expected, _) = run(&genes(DEFAULT_GENES), &text, VariantKind::Snv, false).unwrap();
    for other in [
        text.replace('\n', "\r\n"),
        text.trim_end_matches('\n').to_owned(),
    ] {
        assert_eq!(
            run(&genes(DEFAULT_GENES), &other, VariantKind::Snv, false)
                .unwrap()
                .0,
            expected
        );
    }
}

#[test]
fn empty_shard_is_explicit_zero_count_and_still_validates_header_context() {
    let (out, stats) = run(&genes(DEFAULT_GENES), &table(&[]), VariantKind::Snv, false).unwrap();
    assert!(out.is_empty());
    assert_eq!(
        stats,
        Statistics {
            input_rows: 0,
            output_records: 0
        }
    );
}

#[test]
fn output_write_failure_is_propagated() {
    struct BrokenWriter;
    impl Write for BrokenWriter {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::new(
                std::io::ErrorKind::BrokenPipe,
                "fixture consumer stopped",
            ))
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let error = project(
        Cursor::new(genes(DEFAULT_GENES)),
        Cursor::new(table(&[row()])),
        BrokenWriter,
        "chr1",
        VariantKind::Snv,
        false,
    )
    .unwrap_err();
    assert!(error.to_string().contains("fixture consumer stopped"));
}

static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);
struct TempDir(std::path::PathBuf);
impl TempDir {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "priva-spliceai2-test-{}-{}",
            std::process::id(),
            TEMP_COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        Self(path)
    }
}
impl Drop for TempDir {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).unwrap();
    }
}

#[test]
fn cli_stdin_stdout_counts_and_no_header_contract() {
    let temp = TempDir::new();
    let context = temp.0.join("genes.tsv");
    std::fs::write(&context, genes(DEFAULT_GENES)).unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_priva"))
        .args(["spliceai2-vcf", "--genes"])
        .arg(&context)
        .args([
            "--chrom",
            "chr1",
            "--kind",
            "snv",
            "--input",
            "-",
            "--output",
            "-",
            "--no-header",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(table(&[row()]).as_bytes())
        .unwrap();
    let result = child.wait_with_output().unwrap();
    assert!(result.status.success());
    assert_eq!(
        String::from_utf8(result.stdout).unwrap(),
        run(
            &genes(DEFAULT_GENES),
            &table(&[row()]),
            VariantKind::Snv,
            false
        )
        .unwrap()
        .0
    );
    assert_eq!(
        String::from_utf8(result.stderr).unwrap(),
        "{\"module\":\"spliceai2-vcf\",\"chrom\":\"chr1\",\"kind\":\"snv\",\"input_rows\":1,\"output_records\":1}\n"
    );
}

#[test]
fn cli_refuses_overwrite_and_same_input_output() {
    let temp = TempDir::new();
    let context = temp.0.join("genes.tsv");
    let input = temp.0.join("native.tsv");
    let output = temp.0.join("existing.txt");
    std::fs::write(&context, genes(DEFAULT_GENES)).unwrap();
    std::fs::write(&input, table(&[row()])).unwrap();
    std::fs::write(&output, b"retain me").unwrap();
    for path in [&input, &output] {
        let result = Command::new(env!("CARGO_BIN_EXE_priva"))
            .args(["spliceai2-vcf", "--genes"])
            .arg(&context)
            .args(["--chrom", "chr1", "--kind", "snv", "--input"])
            .arg(&input)
            .arg("--output")
            .arg(path)
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(2));
    }
    assert_eq!(std::fs::read_to_string(&input).unwrap(), table(&[row()]));
    assert_eq!(std::fs::read(&output).unwrap(), b"retain me");
}

#[test]
fn late_input_error_returns_failure_without_success_counts() {
    let temp = TempDir::new();
    let context = temp.0.join("genes.tsv");
    let input = temp.0.join("native.tsv");
    std::fs::write(&context, genes(DEFAULT_GENES)).unwrap();
    let mut bad = row();
    set(&mut bad, "pos", "21");
    set(&mut bad, "jxn_loss_ref_score_2", "NaN");
    std::fs::write(&input, table(&[row(), bad])).unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_priva"))
        .args(["spliceai2-vcf", "--genes"])
        .arg(&context)
        .args(["--chrom", "chr1", "--kind", "snv", "--input"])
        .arg(&input)
        .args(["--output", "-", "--no-header"])
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(2));
    // A validated prefix can already be streamed. It must never be published as
    // completed merely because bytes exist; the driver must check every exit.
    assert_eq!(
        String::from_utf8(result.stdout).unwrap(),
        run(
            &genes(DEFAULT_GENES),
            &table(&[row()]),
            VariantKind::Snv,
            false
        )
        .unwrap()
        .0
    );
    let stderr = String::from_utf8(result.stderr).unwrap();
    assert!(stderr.contains("native TSV line 3"));
    assert!(stderr.contains("jxn_loss_ref_score_2"));
    assert!(!stderr.contains("input_rows"));
}

#[test]
fn cli_bad_options_and_help() {
    let binary = env!("CARGO_BIN_EXE_priva");
    for args in [
        vec!["spliceai2-vcf"],
        vec!["spliceai2-vcf", "--no-header", "--no-header"],
        vec!["spliceai2-vcf", "--input", "x", "--input", "y"],
        vec!["spliceai2-vcf", "--unknown"],
    ] {
        let result = Command::new(binary).args(args).output().unwrap();
        assert_eq!(result.status.code(), Some(2));
    }
    let result = Command::new(binary)
        .args(["spliceai2-vcf", "--help"])
        .output()
        .unwrap();
    assert!(result.status.success());
    assert!(
        String::from_utf8(result.stdout)
            .unwrap()
            .contains("annotation-only")
    );
}
