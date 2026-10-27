//! Strict named transport for the PM4/BP3 procedures.
use super::{assign, prepare_bed};
use polars_core::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::OpenOptions;
use std::path::Path;

const FLAGS: [&str; 8] = [
    "splicing_len_changing",
    "splicing_frameshift",
    "5UTR_len_changing",
    "5UTR_frameshift",
    "vep_consq_len_changing",
    "splicing_span_intol_domain",
    "5UTR_span_intol_domain",
    "span_functional_domains",
];

fn err(value: impl std::fmt::Display) -> String {
    value.to_string()
}

fn arguments(arguments: &[String]) -> Result<(String, HashMap<String, String>), String> {
    let mut args = arguments.iter().cloned();
    let procedure = args.next().ok_or("expected prepare-bed or assign")?;
    let mut options = HashMap::new();
    while let Some(key) = args.next() {
        if !["--input", "--output", "--overlap"].contains(&key.as_str()) {
            return Err(format!("unknown option {key}"));
        }
        let value = args
            .next()
            .ok_or_else(|| format!("missing value for {key}"))?;
        if options.insert(key.clone(), value).is_some() {
            return Err(format!("duplicate option {key}"));
        }
    }
    let expected = if procedure == "prepare-bed" {
        vec!["--input", "--output"]
    } else if procedure == "assign" {
        vec!["--input", "--output", "--overlap"]
    } else {
        return Err(format!("unknown procedure {procedure}"));
    };
    if options.len() != expected.len() || expected.iter().any(|key| !options.contains_key(*key)) {
        return Err("incorrect procedure options".into());
    }
    Ok((procedure, options))
}

fn read(
    path: &str,
    expected: &[&str],
    integers: &[&str],
    booleans: &[&str],
) -> Result<DataFrame, String> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .map_err(err)?;
    let header = reader.headers().map_err(err)?.clone();
    if header.iter().collect::<Vec<_>>() != expected {
        return Err("input header mismatch".into());
    }
    let mut cells: Vec<Vec<String>> = vec![Vec::new(); expected.len()];
    for (ordinal, row) in reader.records().enumerate() {
        let row = row.map_err(err)?;
        if row.get(0) != Some(ordinal.to_string().as_str()) {
            return Err("nonconsecutive input ordinal".into());
        }
        for (column, value) in cells.iter_mut().zip(row.iter()) {
            column.push(value.into());
        }
    }
    let mut columns = Vec::new();
    for (name, values) in expected.iter().zip(cells) {
        let series = if integers.contains(name) {
            let numbers = values
                .into_iter()
                .map(|value| value.parse::<i64>().map_err(err))
                .collect::<Result<Vec<_>, _>>()?;
            Series::new((*name).into(), numbers)
        } else if booleans.contains(name) {
            let numbers = values
                .into_iter()
                .map(|value| match value.as_str() {
                    "0" => Ok(false),
                    "1" => Ok(true),
                    _ => Err("non-boolean input field".into()),
                })
                .collect::<Result<Vec<_>, String>>()?;
            Series::new((*name).into(), numbers)
        } else {
            Series::new((*name).into(), values)
        };
        columns.push(series.into());
    }
    DataFrame::new(columns).map_err(err)
}

fn write(path: &str, frame: &DataFrame) -> Result<(), String> {
    // Results are private transport files. A failed process never exposes a
    // result plan; create_new also prevents accidental reuse/overwrite.
    let handle = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(err)?;
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(handle);
    writer
        .write_record(frame.get_column_names().iter().map(|name| name.as_str()))
        .map_err(err)?;
    for row in 0..frame.height() {
        let mut values = Vec::new();
        for column in frame.get_columns() {
            values.push(match column.get(row).map_err(err)? {
                AnyValue::String(value) => value.into(),
                // Polars may expose a singleton string column as an owned scalar.
                AnyValue::StringOwned(value) => value.as_str().into(),
                AnyValue::Int64(value) => value.to_string(),
                AnyValue::Boolean(value) => {
                    if value {
                        "1".into()
                    } else {
                        "0".into()
                    }
                }
                other => return Err(format!("unsupported result dtype {other:?}")),
            });
        }
        writer.write_record(values).map_err(err)?;
    }
    writer.flush().map_err(err)?;
    writer.into_inner().map_err(err)?.sync_all().map_err(err)
}

fn execute(args: &[String]) -> Result<(), String> {
    let (procedure, options) = arguments(args)?;
    let input = &options["--input"];
    let output = &options["--output"];
    if Path::new(output).exists() {
        return Err("output already exists".into());
    }
    if procedure == "prepare-bed" {
        let input = read(
            input,
            &["ordinal", "chrom", "pos", "variant_id", "ref"],
            &["ordinal", "pos"],
            &[],
        )?;
        write(output, &prepare_bed(&input)?)
    } else {
        let mut expected = vec![
            "ordinal",
            "variant_id",
            "DOMAINS",
            "Consequence",
            "NMD",
            "LoF_filter",
        ];
        expected.extend(FLAGS);
        expected.extend(["pvs1", "loc_intol_domain"]);
        let mut booleans = FLAGS.to_vec();
        booleans.push("loc_intol_domain");
        let input = read(input, &expected, &["ordinal", "pvs1"], &booleans)?;
        let mut overlap = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(&options["--overlap"])
            .map_err(err)?;
        if overlap.headers().map_err(err)?.iter().collect::<Vec<_>>() != vec!["variant_id"] {
            return Err("overlap header mismatch".into());
        }
        let mut ids = HashSet::new();
        for row in overlap.records() {
            ids.insert(row.map_err(err)?[0].to_owned());
        }
        write(output, &assign(&input, &ids)?)
    }
}

pub(crate) fn dispatch(args: &[String]) -> crate::cli::Result<()> {
    if args == ["--help"] {
        println!(
            "priva pm4-bp3 prepare-bed --input INPUT.tsv --output NEW.tsv\npriva pm4-bp3 assign --input INPUT.tsv --overlap OVERLAP.tsv --output NEW.tsv\nOpt-in typed procedures; interval intersection remains in pybedtools/bedtools."
        );
        return Ok(());
    }
    execute(args).map_err(crate::cli::invalid)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering};

    static NEXT_FIXTURE: AtomicU64 = AtomicU64::new(0);

    struct Fixture(PathBuf);

    impl Fixture {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "priva-pm4-transport-{}-{}",
                std::process::id(),
                NEXT_FIXTURE.fetch_add(1, Ordering::Relaxed)
            ));
            std::fs::create_dir(&path).unwrap();
            Self(path)
        }

        fn prepare(&self) -> Vec<String> {
            vec![
                "prepare-bed".into(),
                "--input".into(),
                self.0.join("input.tsv").to_str().unwrap().into(),
                "--output".into(),
                self.0.join("output.tsv").to_str().unwrap().into(),
            ]
        }
    }

    impl Drop for Fixture {
        fn drop(&mut self) {
            // Keep the exact failed transport files for diagnosis.
            if !std::thread::panicking() {
                std::fs::remove_dir_all(&self.0).unwrap();
            }
        }
    }

    #[test]
    fn singleton_unicode_prepare_preserves_tsv_and_ref_length() {
        let fixture = Fixture::new();
        let input = "ordinal\tchrom\tpos\tvariant_id\tref\n0\tchr1\t10\t\"v\tλ\"\tAλ\n";
        std::fs::write(fixture.0.join("input.tsv"), input).unwrap();

        execute(&fixture.prepare()).unwrap();

        let expected = "ordinal\tchrom\tstart\tend\tvariant_id\n0\tchr1\t9\t11\t\"v\tλ\"\n";
        assert_eq!(
            std::fs::read_to_string(fixture.0.join("output.tsv")).unwrap(),
            expected
        );
        assert_eq!(
            std::fs::read_to_string(fixture.0.join("input.tsv")).unwrap(),
            input
        );
        assert_eq!(
            execute(&fixture.prepare()).unwrap_err(),
            "output already exists"
        );
        assert_eq!(
            std::fs::read_to_string(fixture.0.join("output.tsv")).unwrap(),
            expected
        );
    }

    #[test]
    fn empty_typed_prepare_emits_only_header() {
        let fixture = Fixture::new();
        let input = "ordinal\tchrom\tpos\tvariant_id\tref\n";
        std::fs::write(fixture.0.join("input.tsv"), input).unwrap();

        execute(&fixture.prepare()).unwrap();

        assert_eq!(
            std::fs::read_to_string(fixture.0.join("output.tsv")).unwrap(),
            "ordinal\tchrom\tstart\tend\tvariant_id\n"
        );
        assert_eq!(
            std::fs::read_to_string(fixture.0.join("input.tsv")).unwrap(),
            input
        );
    }
}
