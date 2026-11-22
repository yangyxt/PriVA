use priva_tabular::acmg_scoring::{
    CRITERIA_ORDER,
    calculate_posterior_probability::{Probability, calculate_posterior_probability},
};
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
            "priva-posterior-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn input(&self, rows: &[Vec<String>]) -> PathBuf {
        let path = self.0.join("input.csv");
        let mut writer = csv::Writer::from_path(&path).unwrap();
        writer
            .write_record(["source_position"].into_iter().chain(CRITERIA_ORDER))
            .unwrap();
        for row in rows {
            writer.write_record(row).unwrap();
        }
        writer.flush().unwrap();
        path
    }
    fn invoke(&self, input: &PathBuf, output: &PathBuf) -> std::process::Output {
        Command::new(env!("CARGO_BIN_EXE_priva"))
            .args(["posterior-probability", "--input"])
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

fn row(position: usize) -> Vec<String> {
    let mut result = vec![position.to_string()];
    result.extend(vec!["0".to_owned(); 28]);
    result
}

#[test]
fn actual_cli_preserves_scalar_kind_raw_bits_classes_and_duplicate_value_order() {
    let fixture = Fixture::new();
    let mut benign = row(1);
    benign[17] = "5".into();
    let input = fixture.input(&[row(0), benign, row(2)]);
    let output = fixture.0.join("new.csv");
    let before = fs::read(&input).unwrap();
    let process = fixture.invoke(&input, &output);
    assert!(
        process.status.success(),
        "{}",
        String::from_utf8_lossy(&process.stderr)
    );
    let mut reader = csv::Reader::from_path(output).unwrap();
    assert_eq!(
        reader.headers().unwrap().iter().collect::<Vec<_>>(),
        [
            "source_position",
            "probability_kind",
            "probability_value",
            "ACMG_class"
        ]
    );
    let rows: Vec<_> = reader.records().map(Result::unwrap).collect();
    assert_eq!(rows.len(), 3);
    let (probability, class) = calculate_posterior_probability(&[]);
    let Probability::Float(value) = probability else {
        panic!("non-BA scalar kind changed")
    };
    for position in [0, 2] {
        assert_eq!(&rows[position][0], position.to_string());
        assert_eq!(&rows[position][1], "float64");
        assert_eq!(&rows[position][2], format!("{:016x}", value.to_bits()));
        assert_eq!(&rows[position][3], class);
    }
    assert_eq!(
        rows[1].iter().collect::<Vec<_>>(),
        ["1", "int", "0", "Benign"]
    );
    assert_eq!(fs::read(input).unwrap(), before);
}

#[test]
fn actual_cli_rejects_every_noncanonical_strength_before_publication() {
    for value in ["", "NaN", "1.0", "1.5", "-1", "6", "01", "+1", "true"] {
        let fixture = Fixture::new();
        let mut bad = row(0);
        bad[1] = value.into();
        let input = fixture.input(&[bad]);
        let output = fixture.0.join("new.csv");
        assert!(!fixture.invoke(&input, &output).status.success(), "{value}");
        assert!(!output.exists());
    }
}

#[test]
fn actual_cli_rejects_empty_and_incomplete_or_reordered_positions() {
    for rows in [
        vec![],
        vec![row(1)],
        vec![row(0), row(0)],
        vec![row(0), row(2)],
    ] {
        let fixture = Fixture::new();
        let input = fixture.input(&rows);
        let output = fixture.0.join("new.csv");
        assert!(!fixture.invoke(&input, &output).status.success());
        assert!(!output.exists());
    }
}

#[test]
fn actual_cli_rejects_changed_missing_or_duplicate_headers() {
    for header in [
        "source_position,PVS1\n0,0\n",
        "source_position,PVS1,PVS1\n0,0,0\n",
        "\n",
    ] {
        let fixture = Fixture::new();
        let input = fixture.0.join("input.csv");
        fs::write(&input, header).unwrap();
        let output = fixture.0.join("new.csv");
        assert!(!fixture.invoke(&input, &output).status.success());
        assert!(!output.exists());
    }
}

#[test]
fn later_malformed_row_cannot_publish_an_earlier_success() {
    let fixture = Fixture::new();
    let mut bad = row(1);
    bad[28] = "6".into();
    let input = fixture.input(&[row(0), bad]);
    let output = fixture.0.join("new.csv");
    assert!(!fixture.invoke(&input, &output).status.success());
    assert!(!output.exists());
}

#[test]
fn existing_output_and_literal_or_hardlink_aliases_stay_unchanged() {
    let fixture = Fixture::new();
    let input = fixture.input(&[row(0)]);
    let before = fs::read(&input).unwrap();
    let output = fixture.0.join("existing.csv");
    fs::write(&output, b"retain").unwrap();
    assert!(!fixture.invoke(&input, &output).status.success());
    assert_eq!(fs::read(&output).unwrap(), b"retain");
    assert!(!fixture.invoke(&input, &input).status.success());
    let alias = fixture.0.join("alias.csv");
    fs::hard_link(&input, &alias).unwrap();
    assert!(!fixture.invoke(&input, &alias).status.success());
    assert_eq!(fs::read(&input).unwrap(), before);
    assert_eq!(fs::read(alias).unwrap(), before);
}

#[test]
fn named_options_reject_duplicates_unknown_model_options_and_missing_values() {
    let fixture = Fixture::new();
    let input = fixture.input(&[row(0)]);
    for tail in [
        vec!["--prior", "0.1"],
        vec!["--input", "duplicate"],
        vec!["--output"],
    ] {
        let output = fixture.0.join("new.csv");
        let process = Command::new(env!("CARGO_BIN_EXE_priva"))
            .args(["posterior-probability", "--input"])
            .arg(&input)
            .arg("--output")
            .arg(&output)
            .args(tail)
            .output()
            .unwrap();
        assert!(!process.status.success());
        assert!(!output.exists());
    }
}

#[test]
fn help_exposes_the_default_only_procedure_without_input_or_publication() {
    let process = Command::new(env!("CARGO_BIN_EXE_priva"))
        .args(["posterior-probability", "--help"])
        .output()
        .unwrap();
    assert!(process.status.success());
    assert!(
        String::from_utf8(process.stdout)
            .unwrap()
            .contains("default prior0.1")
    );
}
