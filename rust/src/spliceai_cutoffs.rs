//! Rust version of `scripts/spliceai_cutoffs.py`: the score cutoffs PriVA's
//! splicing and ACMG rules use, chosen by the SpliceAI version of the
//! configured score files (v1 or SpliceAI2).
//!
//! Command line, as the script takes it: `priva spliceai-cutoffs CONFIG
//! [SNV_VCF INDEL_VCF]` prints the version and its cutoffs on one line, or
//! exits with code 1 and a message when the files differ or cannot be read.
//! Any other number of arguments prints the usage text (the script's
//! docstring) and exits with code 1. Like the script, it has no `--help`.
//!
//! Library functions for later conversions: `cutoffs_for_config` and the
//! `SPLICEAI_V1` and `SPLICEAI2` tables.
//!
//! The configuration is read with the flat YAML reader of `yaml_config`: a
//! construct it does not support stops the run with exit code 1 (the script
//! also exits with 1 for YAML it cannot parse).
use crate::cli::Result;
use crate::spliceai_cache_input::{Done, Failure, canonical_assembly, inspect_header};
use crate::yaml_config::{self, ConfigValue, Parsed, Scalar, YamlError};
use std::path::Path;

const USAGE_TEXT: &str = r##"SpliceAI cutoffs for PriVA's splicing and ACMG rules, by score version.

PriVA's rules were written for SpliceAI v1 scores. SpliceAI2 (Jaganathan et al.,
"A unified framework for quantitative splicing and transcript prediction", 2026)
reports the same four delta scores, and its authors recommend their own
thresholds on the largest of them (their summary score). Their README maps each
threshold to the SpliceAI v1 tier it replaces:

  tier                      SpliceAI v1   SpliceAI2
  high recall                  0.2          0.1
  balance of both              0.5          0.25
  high precision               0.8          0.5

The manuscript calls variants scoring above 0.1 predicted splice-altering and
prioritises disease variants at 0.25. Each PriVA rule takes the SpliceAI2
threshold of the tier its SpliceAI v1 cutoff belongs to:

  rule (where)                                       v1    SpliceAI2
  splice_altering (SpliceAI_interpretation, PVS1/PP3) 0.8   0.5   high precision
  nearby_acceptor_rescue (SpliceAI_interpretation)    0.5   0.25  balance
  known_pathogenic_comparison (PS1)                   0.5   0.25  balance
  no_splicing_impact (BP4; score below the cutoff)    0.1   0.1   not predicted
                                                                  splice-altering
  donor_loss_putative_nmd (ClinVar gene statistics)   0.5   0.25  balance

The version comes from the headers of the configured spliceai_snv_prescore and
spliceai_indel_prescore files, which must be the same version and match the
configured assembly (SpliceAI2 is GRCh38 only). Indels use the same cutoffs.

Command line: spliceai_cutoffs.py CONFIG [SNV_VCF INDEL_VCF]
prints the version and its cutoffs, or exits with an error when the files differ.
"##;

/// The cutoffs for one SpliceAI score version.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Cutoffs {
    pub(crate) splice_altering: f64,
    pub(crate) nearby_acceptor_rescue: f64,
    pub(crate) known_pathogenic_comparison: f64,
    pub(crate) no_splicing_impact: f64,
    pub(crate) donor_loss_putative_nmd: f64,
}

pub(crate) const SPLICEAI_V1: Cutoffs = Cutoffs {
    splice_altering: 0.8,
    nearby_acceptor_rescue: 0.5,
    known_pathogenic_comparison: 0.5,
    no_splicing_impact: 0.1,
    donor_loss_putative_nmd: 0.5,
};

pub(crate) const SPLICEAI2: Cutoffs = Cutoffs {
    splice_altering: 0.5,
    nearby_acceptor_rescue: 0.25,
    known_pathogenic_comparison: 0.25,
    no_splicing_impact: 0.1,
    donor_loss_putative_nmd: 0.25,
};

impl Cutoffs {
    fn entries(&self) -> [(&'static str, f64); 5] {
        [
            ("splice_altering", self.splice_altering),
            ("nearby_acceptor_rescue", self.nearby_acceptor_rescue),
            (
                "known_pathogenic_comparison",
                self.known_pathogenic_comparison,
            ),
            ("no_splicing_impact", self.no_splicing_impact),
            ("donor_loss_putative_nmd", self.donor_loss_putative_nmd),
        ]
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum ScoreVersion {
    V1,
    V2,
}

impl ScoreVersion {
    fn name(self) -> &'static str {
        match self {
            ScoreVersion::V1 => "SpliceAI v1",
            ScoreVersion::V2 => "SpliceAI2",
        }
    }

    pub(crate) fn cutoffs(self) -> Cutoffs {
        match self {
            ScoreVersion::V1 => SPLICEAI_V1,
            ScoreVersion::V2 => SPLICEAI2,
        }
    }
}

/// A setting as the script receives it: missing or null, a string, or a value of another type.
#[derive(Debug, Clone)]
enum Setting {
    Absent,
    Text(String),
    Other(ConfigValue),
}

impl Setting {
    fn from_config(value: Option<&ConfigValue>) -> Setting {
        match value {
            None | Some(ConfigValue::Scalar(Scalar::Null)) => Setting::Absent,
            Some(ConfigValue::Scalar(Scalar::Str(text))) => Setting::Text(text.clone()),
            Some(other) => Setting::Other(other.clone()),
        }
    }

    /// `value in (None, "", "null")`
    fn is_unset(&self) -> bool {
        match self {
            Setting::Absent => true,
            Setting::Text(text) => text.is_empty() || text == "null",
            Setting::Other(_) => false,
        }
    }
}

fn setting_repr(value: &ConfigValue) -> String {
    match value {
        ConfigValue::Scalar(scalar) => scalar.python_repr(),
        ConfigValue::Opaque(text) => text.clone(),
    }
}

/// The SpliceAI version of the given score files. Files that are unset or
/// missing are skipped so the annotation step can report them; with neither
/// present the original SpliceAI v1 cutoffs apply.
fn score_version(snv: &Setting, indel: &Setting, assembly: &Setting) -> Done<ScoreVersion> {
    let expected: Option<&'static str> = if assembly.is_unset() {
        None
    } else {
        match assembly {
            Setting::Text(text) => Some(canonical_assembly(text)?),
            Setting::Other(value) => {
                return Err(Failure::Rejected(format!(
                    "Unsupported assembly: {}",
                    setting_repr(value)
                )));
            }
            Setting::Absent => None,
        }
    };
    let mut found: Vec<(&str, &'static str)> = Vec::new();
    for (kind, setting) in [("SNV", snv), ("indel", indel)] {
        if setting.is_unset() {
            continue;
        }
        let path = match setting {
            Setting::Text(text) => text,
            _ => {
                return Err(Failure::Crash(
                    "TypeError: expected str, bytes or os.PathLike object".to_owned(),
                ));
            }
        };
        if !Path::new(path).is_file() {
            continue;
        }
        let header = inspect_header(Path::new(path))?;
        if let (Some(wanted), Some(actual)) = (expected, header.assembly) {
            if actual != wanted {
                return Err(Failure::Rejected(format!(
                    "The SpliceAI {kind} score file is {actual}, but the configured assembly is {wanted}: {path}"
                )));
            }
        }
        found.push((kind, header.format));
    }
    let mut versions: Vec<&'static str> = Vec::new();
    for (_, format) in &found {
        if !versions.contains(format) {
            versions.push(*format);
        }
    }
    if versions.len() > 1 {
        let described: Vec<String> = found
            .iter()
            .map(|(kind, format)| {
                let name = match *format {
                    "vcf_v1" => "SpliceAI v1",
                    "vcf_spliceai2" => "SpliceAI2",
                    other => other,
                };
                format!("{kind} {name}")
            })
            .collect();
        return Err(Failure::Rejected(format!(
            "The SpliceAI SNV and indel score files are different versions: {}",
            described.join(", ")
        )));
    }
    match versions.first().copied().unwrap_or("vcf_v1") {
        "vcf_v1" => Ok(ScoreVersion::V1),
        "vcf_spliceai2" => Ok(ScoreVersion::V2),
        _ => Err(Failure::Rejected(
            "SpliceAI2 publisher tables must be converted to indexed VCF by the installer first"
                .to_owned(),
        )),
    }
}

/// The three settings the cutoffs depend on, from a configuration file.
struct Settings {
    snv: Setting,
    indel: Setting,
    assembly: Setting,
}

fn read_config(config_path: Option<&str>) -> Done<Settings> {
    let unset = Settings {
        snv: Setting::Absent,
        indel: Setting::Absent,
        assembly: Setting::Absent,
    };
    let path = match config_path {
        Some(path) if !path.is_empty() && Path::new(path).is_file() => path,
        _ => return Ok(unset),
    };
    let bytes = std::fs::read(path)?;
    let text = String::from_utf8(bytes).map_err(|error| {
        Failure::Rejected(format!("'utf-8' codec can't decode the file: {error}"))
    })?;
    match yaml_config::parse_config(&text) {
        Ok(Parsed::Mapping(document)) => Ok(Settings {
            snv: Setting::from_config(document.get("spliceai_snv_prescore")),
            indel: Setting::from_config(document.get("spliceai_indel_prescore")),
            assembly: Setting::from_config(document.get("assembly")),
        }),
        Ok(Parsed::Empty) => Ok(unset),
        Ok(Parsed::List) => Err(Failure::Rejected(
            "PriVA config must be a mapping".to_owned(),
        )),
        Err(error @ YamlError::Invalid(_)) | Err(error @ YamlError::Unsupported(_)) => {
            Err(Failure::Crash(format!("YAML error: {error}")))
        }
    }
}

/// The cutoffs for the SpliceAI score version configured in `config_path`.
pub(crate) fn cutoffs_for_config(config_path: Option<&Path>) -> Done<Cutoffs> {
    let text = config_path.map(|path| path.to_string_lossy().into_owned());
    let settings = read_config(text.as_deref())?;
    Ok(score_version(&settings.snv, &settings.indel, &settings.assembly)?.cutoffs())
}

fn describe(version: ScoreVersion) -> String {
    let entries: Vec<String> = version
        .cutoffs()
        .entries()
        .iter()
        .map(|(name, value)| format!("{name} {value}"))
        .collect();
    format!("{} cutoffs: {}", version.name(), entries.join(", "))
}

fn run(args: &[String]) -> Done<ScoreVersion> {
    let settings = read_config(Some(args[0].as_str()))?;
    let (snv, indel) = if args.len() == 3 {
        (
            Setting::Text(args[1].clone()),
            Setting::Text(args[2].clone()),
        )
    } else {
        (settings.snv.clone(), settings.indel.clone())
    };
    score_version(&snv, &indel, &settings.assembly)
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args.len() != 1 && args.len() != 3 {
        eprintln!("{USAGE_TEXT}");
        std::process::exit(1);
    }
    match run(args) {
        Ok(version) => {
            println!("{}", describe(version));
            Ok(())
        }
        Err(Failure::Rejected(message)) => {
            eprintln!("Cannot use the SpliceAI score files: {message}");
            std::process::exit(1)
        }
        Err(Failure::Crash(message)) => {
            eprintln!("priva spliceai-cutoffs: unhandled error: {message}");
            std::process::exit(1)
        }
        Err(Failure::Helper(code)) => std::process::exit(code),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;
    use std::path::PathBuf;

    const HEADER: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    const V2_HEADER: &str = "##fileformat=VCFv4.2\n##SpliceAI2Version=2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    fn scratch(name: &str) -> PathBuf {
        let folder = std::env::temp_dir().join(format!(
            "priva_spliceai_cutoffs_{}_{name}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&folder);
        std::fs::create_dir_all(&folder).unwrap();
        folder
    }

    fn config(folder: &Path, snv: &str, indel: &str, assembly: &str) -> PathBuf {
        let mut paths = Vec::new();
        for (name, header) in [("snv.vcf.gz", snv), ("indel.vcf.gz", indel)] {
            let path = folder.join(name);
            let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
            encoder.write_all(header.as_bytes()).unwrap();
            std::fs::write(&path, encoder.finish().unwrap()).unwrap();
            paths.push(path);
        }
        let config = folder.join("config.yaml");
        std::fs::write(
            &config,
            format!(
                "assembly: {assembly}\nspliceai_snv_prescore: {}\nspliceai_indel_prescore: {}\n",
                paths[0].display(),
                paths[1].display()
            ),
        )
        .unwrap();
        config
    }

    #[test]
    fn each_version_selects_its_own_cutoffs() {
        let folder = scratch("versions");
        let v1 = config(&folder, HEADER, HEADER, "hg38");
        assert_eq!(cutoffs_for_config(Some(&v1)).unwrap(), SPLICEAI_V1);
        let v2 = config(&folder, V2_HEADER, V2_HEADER, "hg38");
        assert_eq!(cutoffs_for_config(Some(&v2)).unwrap(), SPLICEAI2);
        let _ = std::fs::remove_dir_all(&folder);
    }

    #[test]
    fn version_one_applies_without_score_files() {
        assert_eq!(cutoffs_for_config(None).unwrap(), SPLICEAI_V1);
        let missing = Path::new("/definitely/not/here.yaml");
        assert_eq!(cutoffs_for_config(Some(missing)).unwrap(), SPLICEAI_V1);
    }

    #[test]
    fn mixed_versions_and_a_wrong_assembly_are_refused() {
        let folder = scratch("refused");
        let mixed = config(&folder, V2_HEADER, HEADER, "hg38");
        match cutoffs_for_config(Some(&mixed)) {
            Err(Failure::Rejected(message)) => {
                assert!(message.contains("different versions"), "{message}")
            }
            other => panic!("{other:?}"),
        }
        let wrong = config(&folder, V2_HEADER, V2_HEADER, "hg19");
        match cutoffs_for_config(Some(&wrong)) {
            Err(Failure::Rejected(message)) => {
                assert!(
                    message.contains("configured assembly is GRCh37"),
                    "{message}"
                )
            }
            other => panic!("{other:?}"),
        }
        let _ = std::fs::remove_dir_all(&folder);
    }

    #[test]
    fn the_printed_line_lists_every_rule() {
        assert_eq!(
            describe(ScoreVersion::V1),
            "SpliceAI v1 cutoffs: splice_altering 0.8, nearby_acceptor_rescue 0.5, known_pathogenic_comparison 0.5, no_splicing_impact 0.1, donor_loss_putative_nmd 0.5"
        );
        assert_eq!(
            describe(ScoreVersion::V2),
            "SpliceAI2 cutoffs: splice_altering 0.5, nearby_acceptor_rescue 0.25, known_pathogenic_comparison 0.25, no_splicing_impact 0.1, donor_loss_putative_nmd 0.25"
        );
    }

    #[test]
    fn unset_settings_are_skipped() {
        for unset in [
            Setting::Absent,
            Setting::Text(String::new()),
            Setting::Text("null".into()),
        ] {
            assert!(unset.is_unset());
            assert_eq!(
                score_version(&unset, &unset, &unset).unwrap(),
                ScoreVersion::V1
            );
        }
        let bad = Setting::Text("hg37".into());
        assert!(matches!(
            score_version(&Setting::Absent, &Setting::Absent, &bad),
            Err(Failure::Rejected(message)) if message == "Unsupported assembly: 'hg37'"
        ));
    }
}
