use priva_tabular::acmg_scoring::{CRITERIA_ORDER, criteria_summary::criteria_summary};
use std::{
    fs,
    path::PathBuf,
    process::Command,
    sync::atomic::{AtomicUsize, Ordering},
};

static NEXT: AtomicUsize = AtomicUsize::new(0);

struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "priva-summary-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn input(&self, strength: &str, kind: &str, evidence: &str, position: &str) -> PathBuf {
        let path = self.0.join("input.csv");
        let mut w = csv::Writer::from_path(&path).unwrap();
        let header: Vec<&str> = ["source_position"]
            .into_iter()
            .chain(CRITERIA_ORDER)
            .chain(["evidence_kind", "evidence"])
            .collect();
        w.write_record(header).unwrap();
        let mut row = vec![position.to_owned()];
        row.extend(vec![strength.to_owned(); 28]);
        row.extend([kind.to_owned(), evidence.to_owned()]);
        w.write_record(row).unwrap();
        w.flush().unwrap();
        path
    }
    fn invoke(&self, input: &PathBuf, output: &PathBuf) -> std::process::Output {
        Command::new(env!("CARGO_BIN_EXE_priva"))
            .args(["criteria-summary", "--input"])
            .arg(input)
            .arg("--output")
            .arg(output)
            .output()
            .unwrap()
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}

fn columns(rows: usize) -> Vec<Vec<i64>> {
    vec![vec![0; rows]; CRITERIA_ORDER.len()]
}

#[test]
fn summaries_replace_without_merging_and_keep_duplicate_row_order() {
    let mut values = columns(3);
    values[5] = vec![2, 1, 2];
    let plan =
        criteria_summary(&values, &[Some("PVS1_Very, BA1_Stand".into()), None, None]).unwrap();
    assert_eq!(
        plan.summaries,
        vec![
            Some("PVS1_Very Strong;BA1_Stand Alone".into()),
            Some("PM1_Supporting".into()),
            Some("PM1".into())
        ]
    );
    assert_eq!(plan.overridden, vec![true, false, false]);
    assert_eq!(
        plan.strengths.column("PVS1").unwrap().i64().unwrap().get(0),
        Some(4)
    );
    assert_eq!(
        plan.strengths.column("BA1").unwrap().i64().unwrap().get(0),
        Some(5)
    );
    assert_eq!(
        plan.strengths.column("PM1").unwrap().i64().unwrap().get(0),
        Some(0)
    );
}

#[test]
fn empty_and_nonstring_resolved_evidence_do_not_erase_priva_criteria() {
    let mut values = columns(3);
    values[13] = vec![1, 1, 1];
    let plan = criteria_summary(
        &values,
        &[None, Some("\u{001c}\u{001f}".into()), Some(", ,".into())],
    )
    .unwrap();
    assert_eq!(
        plan.summaries,
        vec![Some("PP3".into()), Some("PP3".into()), None]
    );
    assert_eq!(plan.overridden, vec![false, false, true]);
    assert_eq!(
        plan.strengths.column("PP3").unwrap().i64().unwrap().get(2),
        Some(0)
    );
}

#[test]
fn malformed_second_row_pvs1_precedes_first_row_ps1_failure() {
    let error = criteria_summary(
        &columns(2),
        &[Some("XPS1".into()), Some("PVS1_Unknown".into())],
    );
    let message = match error {
        Ok(_) => panic!("malformed evidence accepted"),
        Err(e) => e.to_string(),
    };
    assert_eq!(message, "KeyError: Unknown");
}

#[test]
fn invalid_strengths_empty_or_misaligned_requests_are_rejected() {
    let mut values = columns(1);
    values[27][0] = 6;
    assert!(criteria_summary(&values, &[Some("PS1".into())]).is_err());
    assert!(criteria_summary(&columns(0), &[]).is_err());
    assert!(criteria_summary(&columns(1), &[]).is_err());
    assert!(criteria_summary(&[], &[None]).is_err());
}

#[test]
fn actual_cli_rejects_bad_typed_rows_before_output_publication() {
    for (strength, kind, evidence, position) in [
        ("6", "none", "", "0"),
        ("1.0", "none", "", "0"),
        ("0", "none", "unexpected", "0"),
        ("0", "str", "PS1_Unknown", "0"),
        ("0", "none", "", "1"),
    ] {
        let f = Fixture::new();
        let input = f.input(strength, kind, evidence, position);
        let output = f.0.join("new.csv");
        let before = fs::read(&input).unwrap();
        let result = f.invoke(&input, &output);
        assert!(!result.status.success());
        assert!(!output.exists());
        assert_eq!(fs::read(input).unwrap(), before);
    }
}

#[test]
fn actual_cli_preserves_existing_output_and_literal_or_hardlink_input_aliases() {
    let f = Fixture::new();
    let input = f.input("0", "str", "PS1", "0");
    let output = f.0.join("existing.csv");
    fs::write(&output, b"keep existing").unwrap();
    assert!(!f.invoke(&input, &output).status.success());
    assert_eq!(fs::read(&output).unwrap(), b"keep existing");
    let before = fs::read(&input).unwrap();
    assert!(!f.invoke(&input, &input).status.success());
    let alias = f.0.join("alias.csv");
    fs::hard_link(&input, &alias).unwrap();
    assert!(!f.invoke(&input, &alias).status.success());
    assert_eq!(fs::read(&input).unwrap(), before);
    assert_eq!(fs::read(alias).unwrap(), before);
}

#[test]
fn actual_cli_success_keeps_opaque_source_positions_and_summary_only_schema() {
    let f = Fixture::new();
    let input = f.input("0", "str", "PS1, PM2_Supporting", "0");
    let output = f.0.join("new.csv");
    let result = f.invoke(&input, &output);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let mut reader = csv::Reader::from_path(output).unwrap();
    let header = reader.headers().unwrap().clone();
    assert_eq!(&header[0], "source_position");
    assert_eq!(header.len(), 32);
    assert!(
        !header
            .iter()
            .any(|name| ["ACMG_quant_score", "ACMG_class"].contains(&name))
    );
    let rows: Vec<_> = reader.records().map(Result::unwrap).collect();
    assert_eq!(rows.len(), 1);
    assert_eq!(&rows[0][0], "0");
    assert_eq!(&rows[0][2], "PS1;PM2_Supporting");
    assert_eq!(&rows[0][3], "1");
}
