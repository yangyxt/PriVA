//! Rust version of `scripts/utr_anno_interpret.py` (`priva utr-anno-interpret`).
//!
//! Reads the combined annotation table (the file given by `--variants_table`),
//! interprets the UTRAnnotator columns `5UTR_consequence` and `5UTR_annotation`
//! row by row, and rewrites the same file with four true/false columns added
//! after the existing ones: `5UTR_lof`, `5UTR_len_changing`, `5UTR_frameshift`
//! and `5UTR_span_intol_domain` (a column of that name that exists is
//! replaced where it stands).
//!
//! The table is read and written the way pandas does it (`pandas_table`), so
//! the other cells are rewritten as pandas prints them: missing-value
//! spellings become empty, `0.10` becomes `0.1`, `007` becomes `7`, an integer
//! column with a gap becomes floats.
//!
//! Behaviour reproduced from the Python script:
//! - a row whose consequence or annotation is not text (a missing value, or a
//!   column that read as numbers) has all four flags false; a column that does
//!   not exist reads as the empty text, as `row.get(name, '')` gives;
//! - consequences and annotations are paired after splitting at `&`, the pair
//!   count being the shorter of the two; each pair is read by the first label
//!   found among premature start codon gain, uORF stop codon loss and uORF
//!   frameshift; an annotation is a list of `name=value` items split at `:`,
//!   and an item without exactly one `=` stops the run with exit code 1
//!   before anything is written, as does a `DistanceToCDS` that Python's
//!   `int()` rejects (in start-gain annotations);
//! - a length-changing effect on a transcript whose exon 1 lists a domain of
//!   the intolerant set also counts as loss of function and sets
//!   `5UTR_span_intol_domain`; the two checked native caches are loaded only when both
//!   paths are non-empty, and a failure to read them stops the run;
//! - exit code 2 for command-line errors, 1 for any other failure.
//!
//! Differences from the Python script: log lines are fewer and carry no
//! timings. Refused with exit code 1 and a message, the table left untouched,
//! where the script succeeds:
//! - a table whose name says it is compressed (pandas reads and writes it
//!   compressed);
//! - what the table reader refuses (`pandas_table`): a first data row longer
//!   than the header, and numbers past the signed 64-bit range that pandas
//!   reads as unsigned integers or keeps as text;
//! - a resource without the declared native schema or a damaged archive is
//!   refused when loaded, before any row is judged.
//!
//! All output is prepared before the table is opened for writing, so a
//! failing run leaves the table as it was.
use crate::cache_convert::intolerant_domains::IntolerantDomainCache;
use crate::cache_convert::transcript_domains::TranscriptDomainCache;
use crate::cli::Result;
use crate::pandas_table::{Cell, Column, Table, read_table, write_table_to};
use crate::python_args::{self, Parsed};
use crate::python_text::python_int_is_valid;
use crate::yaml_config::python_str_repr;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

const USAGE: &str = "usage: priva utr-anno-interpret --variants_table VARIANTS_TABLE --domain_map DOMAIN_MAP --intolerant_domains INTOLERANT_DOMAINS";
const COMPRESSED_SUFFIXES: [&str; 6] = [".gz", ".bz2", ".zip", ".xz", ".zst", ".tar"];

/// A failure while reading the caches or judging a row, which stops the run
/// like a Python exception.
type Outcome<T> = std::result::Result<T, String>;

/// Both resources stay mapped, and are consulted only when both are nonempty.
struct Domains {
    map: TranscriptDomainCache,
    intolerant: IntolerantDomainCache,
}

impl Domains {
    fn load(map_path: &str, intolerant_path: &str) -> Outcome<Option<Domains>> {
        // Keep Python's read order, including reading the second resource when
        // the first map is empty, before testing the pair's truthiness.
        let map = TranscriptDomainCache::open(map_path).map_err(|error| error.to_string())?;
        let intolerant =
            IntolerantDomainCache::open(intolerant_path).map_err(|error| error.to_string())?;
        Ok(if map.len() != 0 && intolerant.len() != 0 {
            Some(Domains { map, intolerant })
        } else {
            None
        })
    }

    /// The native schema's lookup keys are text. A missing or nontext table
    /// cell cannot equal one, and absent/empty exon lists still give false.
    fn exon_one_spans_intolerant_domain(&self, transcript: &Cell) -> Outcome<bool> {
        let Cell::Text(transcript) = transcript else {
            return Ok(false);
        };
        let Some(entry) = self.map.get(transcript) else {
            return Ok(false);
        };
        let Some(exon) = entry
            .exon_domains
            .iter()
            .find(|exon| self.map.text(&exon.exon) == "1")
        else {
            return Ok(false);
        };
        Ok(exon
            .domains
            .iter()
            .any(|domain| self.intolerant.contains(self.map.text(domain))))
    }
}

/// The `name=value` items of an annotation; a later repeat of a name wins.
fn parse_fields(annotation: &str) -> Outcome<HashMap<&str, &str>> {
    let mut fields = HashMap::new();
    for (index, item) in annotation.split(':').enumerate() {
        let parts: Vec<&str> = item.split('=').collect();
        if parts.len() != 2 {
            return Err(format!(
                "ValueError: dictionary update sequence element #{index} has length {}; 2 is required",
                parts.len()
            ));
        }
        fields.insert(parts[0], parts[1]);
    }
    Ok(fields)
}

fn value<'a>(fields: &HashMap<&'a str, &'a str>, name: &str) -> Option<&'a str> {
    fields.get(name).copied()
}

fn is_confident(fields: &HashMap<&str, &str>) -> bool {
    value(fields, "KozakStrength") == Some("Strong")
        || (value(fields, "KozakStrength") == Some("Moderate")
            && value(fields, "Evidence") == Some("True"))
}

/// (loss of function, frameshift, length changing) of a premature start codon gain.
fn interpret_start_gained(annotation: &str) -> Outcome<(bool, bool, bool)> {
    let fields = parse_fields(annotation)?;
    let confidence = is_confident(&fields);
    if let Some(distance) = value(&fields, "DistanceToCDS") {
        if !python_int_is_valid(distance) {
            return Err(format!(
                "ValueError: invalid literal for int() with base 10: {}",
                python_str_repr(distance)
            ));
        }
    }
    let kind = value(&fields, "type");
    let is_lof = kind == Some("OutOfFrame_oORF") && confidence;
    let length_changing = kind == Some("inFrame_oORF") && confidence;
    Ok((is_lof, is_lof, length_changing || is_lof))
}

/// The same for the loss of a uORF stop codon.
fn interpret_stop_lost(annotation: &str) -> Outcome<(bool, bool, bool)> {
    let fields = parse_fields(annotation)?;
    let confidence = is_confident(&fields);
    let alternative_stop = value(&fields, "AltStop") == Some("False");
    let frame = value(&fields, "FrameWithCDS");
    let is_lof = alternative_stop && frame == Some("outOfFrame") && confidence;
    let length_changing = confidence && alternative_stop && frame == Some("inFrame");
    Ok((is_lof, is_lof, length_changing || is_lof))
}

/// The same for a uORF frameshift.
fn interpret_frameshift(annotation: &str) -> Outcome<(bool, bool, bool)> {
    let fields = parse_fields(annotation)?;
    let confidence = is_confident(&fields);
    let no_length = value(&fields, "alt_type_length") == Some("NA");
    let kind = value(&fields, "alt_type");
    let is_lof = no_length && kind == Some("OutOfFrame_oORF") && confidence;
    let length_changing = kind == Some("InFrame_oORF") && no_length && confidence;
    Ok((is_lof, is_lof, length_changing || is_lof))
}

/// (loss of function, frameshift, length changing, spans an intolerant domain) of one row.
/// The consequence and the annotation are `None` when the cell is not text.
fn interpret_utr_annotation(
    consequence: Option<&str>,
    annotation: Option<&str>,
    transcript: &Cell,
    domains: Option<&Domains>,
) -> Outcome<(bool, bool, bool, bool)> {
    let (Some(consequence), Some(annotation)) = (consequence, annotation) else {
        return Ok((false, false, false, false));
    };
    let (mut any_lof, mut any_frameshift, mut any_length, mut any_span) =
        (false, false, false, false);
    for (label, details) in consequence.split('&').zip(annotation.split('&')) {
        let (mut lof, mut frameshift, mut length_changing, mut span) = (false, false, false, false);
        if label.contains("5_prime_UTR_premature_start_codon_gain_variant") {
            (lof, frameshift, length_changing) = interpret_start_gained(details)?;
        } else if label.contains("5_prime_UTR_uORF_stop_codon_loss_variant") {
            (lof, frameshift, length_changing) = interpret_stop_lost(details)?;
        } else if label.contains("5_prime_UTR_uORF_frameshift_variant") {
            (lof, frameshift, length_changing) = interpret_frameshift(details)?;
        }
        if length_changing {
            if let Some(tables) = domains {
                // The variant is in the 5'UTR, so it must be in the first exon.
                if tables.exon_one_spans_intolerant_domain(transcript)? {
                    lof = true;
                    span = true;
                }
            }
        }
        any_lof |= lof;
        any_frameshift |= frameshift;
        any_length |= length_changing;
        any_span |= span;
    }
    Ok((any_lof, any_frameshift, any_length, any_span))
}

/// The text of a cell as the script sees it: only text columns hold text, a
/// missing value is not text, and a column that is not there reads as `''`.
fn text_of<'a>(table: &'a Table, column: Option<usize>, row: usize) -> Option<&'a str> {
    match column {
        None => Some(""),
        Some(column) => match table.cell(column, row) {
            Cell::Text(text) => Some(text),
            _ => None,
        },
    }
}

/// The four result columns in the order the script assigns them.
struct Flags {
    lof: Vec<bool>,
    len_changing: Vec<bool>,
    frameshift: Vec<bool>,
    span_intol_domain: Vec<bool>,
}

fn judge_rows(table: &Table, domains: Option<&Domains>) -> Outcome<Flags> {
    let consequence_column = table.column_index("5UTR_consequence");
    let annotation_column = table.column_index("5UTR_annotation");
    let feature_column = table.column_index("Feature");
    let mut flags = Flags {
        lof: Vec::with_capacity(table.rows),
        len_changing: Vec::with_capacity(table.rows),
        frameshift: Vec::with_capacity(table.rows),
        span_intol_domain: Vec::with_capacity(table.rows),
    };
    for row in 0..table.rows {
        let consequence = text_of(table, consequence_column, row);
        let annotation = text_of(table, annotation_column, row);
        // `row.get('Feature', '')`
        let transcript = match feature_column {
            None => Cell::Text(""),
            Some(column) => table.cell(column, row),
        };
        let (lof, frameshift, length_changing, span) =
            interpret_utr_annotation(consequence, annotation, &transcript, domains)?;
        flags.lof.push(lof);
        flags.len_changing.push(length_changing);
        flags.frameshift.push(frameshift);
        flags.span_intol_domain.push(span);
    }
    Ok(flags)
}

fn process_variants_table(table_path: &str, pickle_paths: Option<(&str, &str)>) -> Outcome<()> {
    eprintln!("INFO: Loading variants table from {table_path}");
    let lower = table_path.to_lowercase();
    if COMPRESSED_SUFFIXES
        .iter()
        .any(|suffix| lower.ends_with(suffix))
    {
        return Err(
            "a compressed table is not supported (pandas would read and write it compressed)"
                .to_owned(),
        );
    }
    let mut table =
        read_table(Path::new(table_path)).map_err(|error| format!("{table_path}: {error}"))?;
    eprintln!("INFO: Loaded {} variants", table.rows);
    let domains = match pickle_paths {
        Some((map, intolerant)) => {
            eprintln!("INFO: Loading domain data from {map} and {intolerant}");
            Domains::load(map, intolerant)?
        }
        None => None,
    };
    eprintln!("INFO: Processing {} rows", table.rows);
    let judged = judge_rows(&table, domains.as_ref());
    let flags = judged?;
    table.set_column("5UTR_lof", Column::Bool(flags.lof));
    table.set_column("5UTR_len_changing", Column::Bool(flags.len_changing));
    table.set_column("5UTR_frameshift", Column::Bool(flags.frameshift));
    table.set_column(
        "5UTR_span_intol_domain",
        Column::Bool(flags.span_intol_domain),
    );
    // Like `to_csv`, the file is opened (and emptied) only when everything is ready.
    let mut output: Vec<u8> = Vec::new();
    write_table_to(&mut output, &table).map_err(|error| error.to_string())?;
    eprintln!("INFO: Saving results to {table_path}");
    fs::write(table_path, &output).map_err(|error| format!("{table_path}: {error}"))
}

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let names = ["variants_table", "domain_map", "intolerant_domains"];
    let values = match python_args::parse_long_options(args, &names) {
        Ok(Parsed::Help) => {
            println!("{USAGE}");
            return Ok(());
        }
        Ok(Parsed::Values(values)) => values,
        Err(message) => {
            eprintln!("{USAGE}");
            eprintln!("priva utr-anno-interpret: error: {message}");
            std::process::exit(2);
        }
    };
    if let Err(message) = python_args::require_options(&values, &names) {
        eprintln!("{USAGE}");
        eprintln!("priva utr-anno-interpret: error: {message}");
        std::process::exit(2);
    }
    // Like the script, the pickles are used only when both paths are non-empty.
    let map = values["domain_map"].as_str();
    let intolerant = values["intolerant_domains"].as_str();
    let pickle_paths = if map.is_empty() || intolerant.is_empty() {
        None
    } else {
        Some((map, intolerant))
    };
    if let Err(message) = process_variants_table(&values["variants_table"], pickle_paths) {
        eprintln!("priva utr-anno-interpret: {message}");
        std::process::exit(1);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const START: &str = "5_prime_UTR_premature_start_codon_gain_variant";
    const STOP: &str = "5_prime_UTR_uORF_stop_codon_loss_variant";
    const SHIFT: &str = "5_prime_UTR_uORF_frameshift_variant";
    const WEAK_START: &str = "KozakStrength=Weak:type=OutOfFrame_oORF";
    const STRONG_START: &str = "KozakStrength=Strong:type=OutOfFrame_oORF";
    const STRONG_STOP: &str = "KozakStrength=Strong:AltStop=False:FrameWithCDS=outOfFrame";
    const STRONG_SHIFT: &str = "KozakStrength=Strong:alt_type=OutOfFrame_oORF:alt_type_length=NA";
    const IN_FRAME_STOP: &str = "KozakStrength=Strong:AltStop=False:FrameWithCDS=inFrame";
    const NO_TRANSCRIPT: Cell<'static> = Cell::Missing;

    // pickle.dumps({"ENST1": {"pp_domains": [1, 5, ("dom:a",)], "exon_domains":
    //   {"1": ["dom:b"], "2": []}}}, protocol=4) as the test bytes, with the
    //   intolerant set {"dom:b", "dom:z"}
    const DOMAIN_MAP: &[u8] = &[
        0x80, 0x04, 0x95, 0x59, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x7d, 0x94, 0x8c, 0x05,
        0x45, 0x4e, 0x53, 0x54, 0x31, 0x94, 0x7d, 0x94, 0x28, 0x8c, 0x0a, 0x70, 0x70, 0x5f, 0x64,
        0x6f, 0x6d, 0x61, 0x69, 0x6e, 0x73, 0x94, 0x5d, 0x94, 0x4b, 0x01, 0x4b, 0x05, 0x8c, 0x05,
        0x64, 0x6f, 0x6d, 0x3a, 0x61, 0x94, 0x87, 0x94, 0x61, 0x8c, 0x0c, 0x65, 0x78, 0x6f, 0x6e,
        0x5f, 0x64, 0x6f, 0x6d, 0x61, 0x69, 0x6e, 0x73, 0x94, 0x7d, 0x94, 0x28, 0x8c, 0x01, 0x31,
        0x94, 0x5d, 0x94, 0x28, 0x68, 0x05, 0x8c, 0x05, 0x64, 0x6f, 0x6d, 0x3a, 0x62, 0x94, 0x65,
        0x8c, 0x01, 0x32, 0x94, 0x5d, 0x94, 0x75, 0x75, 0x73, 0x2e,
    ];
    const INTOLERANT: &[u8] = &[
        0x80, 0x04, 0x95, 0x15, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x8f, 0x94, 0x28, 0x8c,
        0x05, 0x64, 0x6f, 0x6d, 0x3a, 0x62, 0x94, 0x8c, 0x05, 0x64, 0x6f, 0x6d, 0x3a, 0x7a, 0x94,
        0x90, 0x2e,
    ];

    fn interpret(effects: &[(&str, &str)]) -> (bool, bool, bool, bool) {
        let consequence: Vec<&str> = effects.iter().map(|effect| effect.0).collect();
        let annotation: Vec<&str> = effects.iter().map(|effect| effect.1).collect();
        let (joined_consequence, joined_annotation) = (consequence.join("&"), annotation.join("&"));
        interpret_utr_annotation(
            Some(joined_consequence.as_str()),
            Some(joined_annotation.as_str()),
            &NO_TRANSCRIPT,
            None,
        )
        .unwrap()
    }

    #[test]
    fn a_weak_start_does_not_hide_a_strong_stop_loss_or_frameshift() {
        for (label, details) in [(STOP, STRONG_STOP), (SHIFT, STRONG_SHIFT)] {
            assert_eq!(
                interpret(&[(START, WEAK_START), (label, details)]),
                (true, true, true, false)
            );
            assert_eq!(
                interpret(&[(label, details), (START, WEAK_START)]),
                (true, true, true, false)
            );
        }
    }

    #[test]
    fn a_repeated_label_keeps_every_orf() {
        assert_eq!(
            interpret(&[(START, WEAK_START), (START, STRONG_START)]),
            (true, true, true, false)
        );
    }

    #[test]
    fn single_effect_rules_and_other_effects() {
        for (label, details, expected) in [
            (START, WEAK_START, false),
            (START, STRONG_START, true),
            (STOP, STRONG_STOP, true),
            (SHIFT, STRONG_SHIFT, true),
            (
                "5_prime_UTR_uORF_stop_codon_gain_variant",
                STRONG_STOP,
                false,
            ),
            ("5_prime_UTR_start_codon_loss_variant", STRONG_START, false),
            ("-", "-", false),
        ] {
            assert_eq!(
                interpret(&[(label, details)]),
                (expected, expected, expected, false),
                "{label}"
            );
        }
    }

    #[test]
    fn in_frame_effects_change_length_without_losing_function() {
        let moderate = "KozakStrength=Moderate:Evidence=True:type=inFrame_oORF";
        assert_eq!(interpret(&[(START, moderate)]), (false, false, true, false));
        let unproven = "KozakStrength=Moderate:Evidence=False:type=inFrame_oORF";
        assert_eq!(
            interpret(&[(START, unproven)]),
            (false, false, false, false)
        );
        assert_eq!(
            interpret(&[(STOP, IN_FRAME_STOP)]),
            (false, false, true, false)
        );
        let in_frame_shift = "KozakStrength=Strong:alt_type=InFrame_oORF:alt_type_length=NA";
        assert_eq!(
            interpret(&[(SHIFT, in_frame_shift)]),
            (false, false, true, false)
        );
    }

    #[test]
    fn rows_without_text_and_unpaired_effects() {
        let nothing = (false, false, false, false);
        assert_eq!(
            interpret_utr_annotation(None, Some(STRONG_START), &NO_TRANSCRIPT, None).unwrap(),
            nothing
        );
        assert_eq!(
            interpret_utr_annotation(Some(START), None, &NO_TRANSCRIPT, None).unwrap(),
            nothing
        );
        // The shorter list decides how many pairs there are.
        let two = format!("{START}&{STOP}");
        assert_eq!(
            interpret_utr_annotation(Some(two.as_str()), Some(STRONG_START), &NO_TRANSCRIPT, None)
                .unwrap(),
            (true, true, true, false)
        );
    }

    #[test]
    fn malformed_annotations_and_distances_stop_the_run() {
        let run = |label: &str, annotation: &str| {
            interpret_utr_annotation(Some(label), Some(annotation), &NO_TRANSCRIPT, None)
        };
        for annotation in ["", "KozakStrength", "a=b=c", "KozakStrength=Strong:type"] {
            assert!(run(START, annotation).is_err(), "{annotation:?}");
        }
        // Labels that no rule reads never look at their annotation.
        assert!(run("-", "not valid").is_ok());
        assert!(run(START, "DistanceToCDS=12:type=x").is_ok());
        assert!(run(START, "DistanceToCDS= -3_0 ").is_ok());
        assert!(run(START, "DistanceToCDS=1.5").is_err());
        assert!(run(START, "DistanceToCDS=").is_err());
        // Digits of other scripts are digits to int(); the size limit applies.
        assert!(run(START, "DistanceToCDS=\u{661}\u{662}").is_ok());
        assert!(run(START, "DistanceToCDS=\u{ff11}").is_ok());
        assert!(run(START, "DistanceToCDS=\u{b2}").is_err());
        assert!(run(START, &format!("DistanceToCDS={}", "9".repeat(4300))).is_ok());
        assert!(run(START, &format!("DistanceToCDS={}", "9".repeat(4301))).is_err());
        // The distance is not read for the other two kinds.
        assert!(run(STOP, "DistanceToCDS=x:KozakStrength=Weak").is_ok());
    }

    fn native_tables(map: &[u8], intolerant: &[u8]) -> Outcome<Option<Domains>> {
        use std::sync::atomic::{AtomicUsize, Ordering};
        static SERIAL: AtomicUsize = AtomicUsize::new(0);
        let folder = std::env::temp_dir().join(format!(
            "priva_utr_native_{}_{}",
            std::process::id(),
            SERIAL.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&folder).unwrap();
        let source_map = folder.join("map.pkl");
        let source_set = folder.join("set.pkl");
        let native_map = folder.join("map.rkyv");
        let native_set = folder.join("set.rkyv");
        fs::write(&source_map, map).unwrap();
        fs::write(&source_set, intolerant).unwrap();
        let outcome = (|| {
            crate::cache_convert::transcript_domains::convert(
                source_map.to_str().unwrap(),
                &native_map,
            )
            .map_err(|error| error.to_string())?;
            crate::cache_convert::intolerant_domains::convert(
                source_set.to_str().unwrap(),
                &native_set,
            )
            .map_err(|error| error.to_string())?;
            Domains::load(native_map.to_str().unwrap(), native_set.to_str().unwrap())
        })();
        fs::remove_dir_all(folder).unwrap();
        outcome
    }

    fn tables() -> Domains {
        native_tables(DOMAIN_MAP, INTOLERANT).unwrap().unwrap()
    }

    #[test]
    fn a_length_change_in_an_intolerant_exon_one_domain_is_loss_of_function() {
        let domains = tables();
        let with = |transcript: Cell| {
            interpret_utr_annotation(Some(STOP), Some(IN_FRAME_STOP), &transcript, Some(&domains))
                .unwrap()
        };
        assert_eq!(with(Cell::Text("ENST1")), (true, false, true, true));
        assert_eq!(with(Cell::Text("ENST2")), (false, false, true, false));
        assert_eq!(with(Cell::Missing), (false, false, true, false));
        // A number is a different key from any text.
        assert_eq!(with(Cell::Int(1)), (false, false, true, false));
        // A loss of function also counts as length changing, so the check runs here too.
        let out_of_frame = interpret_utr_annotation(
            Some(STOP),
            Some(STRONG_STOP),
            &Cell::Text("ENST1"),
            Some(&domains),
        );
        assert_eq!(out_of_frame.unwrap(), (true, true, true, true));
    }

    #[test]
    fn empty_domain_tables_switch_the_check_off() {
        // pickle.dumps({}) and pickle.dumps(set()) at protocol 4
        let empty_dict: &[u8] = &[
            0x80, 0x04, 0x95, 0x03, 0, 0, 0, 0, 0, 0, 0, 0x7d, 0x94, 0x2e,
        ];
        let empty_set: &[u8] = &[
            0x80, 0x04, 0x95, 0x03, 0, 0, 0, 0, 0, 0, 0, 0x8f, 0x94, 0x2e,
        ];
        assert!(native_tables(empty_dict, empty_set).unwrap().is_none());
        assert!(native_tables(DOMAIN_MAP, empty_set).unwrap().is_none());
        assert!(native_tables(empty_dict, INTOLERANT).unwrap().is_none());
    }

    #[test]
    fn foreign_python_container_shapes_are_refused_before_rows_are_judged() {
        // Runtime accepts the declared native schemas; unsupported Python
        // containers are rejected by conversion rather than rehydrated at runtime.
        for bytes in [b"\x80\x04].".as_slice(), b"\x80\x04).", b"\x80\x04N."] {
            assert!(native_tables(DOMAIN_MAP, bytes).is_err());
            assert!(native_tables(bytes, INTOLERANT).is_err());
        }
    }

    #[test]
    fn missing_columns_read_as_empty_text_and_other_cells_as_not_text() {
        let flags = |text: &str| {
            let table = crate::pandas_table::parse_table(text).unwrap();
            judge_rows(&table, None).map(|flags| flags.lof)
        };
        // No annotation column: the label meets an empty annotation, which is malformed.
        assert!(flags(&format!("5UTR_consequence\n{START}\n")).is_err());
        assert_eq!(flags("5UTR_consequence\n-\n").unwrap(), [false]);
        assert_eq!(flags("5UTR_annotation\nanything\n").unwrap(), [false]);
        assert_eq!(flags("a\n1\n").unwrap(), [false]);
        // Numbers and missing values are not text.
        assert_eq!(
            flags("5UTR_consequence\t5UTR_annotation\n1\t2\nNA\t\n").unwrap(),
            [false, false]
        );
        // A header with no rows.
        assert!(
            flags("5UTR_consequence\t5UTR_annotation\n")
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn a_table_gets_four_columns_and_is_rewritten_like_pandas() {
        let folder = std::env::temp_dir().join(format!("priva_utr_anno_{}", std::process::id()));
        let _ = fs::remove_dir_all(&folder);
        fs::create_dir_all(&folder).unwrap();
        let table = folder.join("table.tsv");
        let table_path = table.to_str().unwrap();
        let text = format!(
            "chrom\tpos\tFeature\t5UTR_consequence\t5UTR_annotation\n1\t007\tENST1\t{START}\t{STRONG_START}\n1\t8\tENST2\t-\tNA\n1\t9\t\t{STOP}&{SHIFT}\t{STRONG_STOP}&{STRONG_SHIFT}\n"
        );
        fs::write(&table, text).unwrap();
        process_variants_table(table_path, None).unwrap();
        let result = fs::read_to_string(&table).unwrap();
        let expected = format!(
            "chrom\tpos\tFeature\t5UTR_consequence\t5UTR_annotation\t5UTR_lof\t5UTR_len_changing\t5UTR_frameshift\t5UTR_span_intol_domain\n1\t7\tENST1\t{START}\t{STRONG_START}\tTrue\tTrue\tTrue\tFalse\n1\t8\tENST2\t-\t\tFalse\tFalse\tFalse\tFalse\n1\t9\t\t{STOP}&{SHIFT}\t{STRONG_STOP}&{STRONG_SHIFT}\tTrue\tTrue\tTrue\tFalse\n"
        );
        assert_eq!(result, expected);
        // Running it again replaces the four columns where they stand.
        process_variants_table(table_path, None).unwrap();
        assert_eq!(fs::read_to_string(&table).unwrap(), expected);
        // A bad annotation leaves the file as it was.
        let broken = format!("5UTR_consequence\t5UTR_annotation\n{START}\tbroken\n");
        fs::write(&table, &broken).unwrap();
        assert!(process_variants_table(table_path, None).is_err());
        assert_eq!(fs::read_to_string(&table).unwrap(), broken);
        // So does a table that cannot be read, and a compressed name.
        fs::write(&table, "a\tb\n1\t2\t3\n").unwrap();
        assert!(process_variants_table(table_path, None).is_err());
        assert_eq!(fs::read_to_string(&table).unwrap(), "a\tb\n1\t2\t3\n");
        let zipped = folder.join("table.tsv.gz");
        fs::write(&zipped, "x").unwrap();
        assert!(process_variants_table(zipped.to_str().unwrap(), None).is_err());
        assert_eq!(fs::read_to_string(&zipped).unwrap(), "x");
        let _ = fs::remove_dir_all(&folder);
    }

    #[test]
    fn the_native_caches_are_read_and_applied_to_the_table() {
        let folder = std::env::temp_dir().join(format!("priva_utr_domains_{}", std::process::id()));
        let _ = fs::remove_dir_all(&folder);
        fs::create_dir_all(&folder).unwrap();
        let table = folder.join("table.tsv");
        let source_map = folder.join("map.pkl");
        let source_set = folder.join("intolerant.pkl");
        let map = folder.join("map.rkyv");
        let intolerant = folder.join("intolerant.rkyv");
        fs::write(&source_map, DOMAIN_MAP).unwrap();
        fs::write(&source_set, INTOLERANT).unwrap();
        crate::cache_convert::transcript_domains::convert(source_map.to_str().unwrap(), &map)
            .unwrap();
        crate::cache_convert::intolerant_domains::convert(
            source_set.to_str().unwrap(),
            &intolerant,
        )
        .unwrap();
        let input = format!(
            "Feature\t5UTR_consequence\t5UTR_annotation\nENST1\t{STOP}\t{IN_FRAME_STOP}\nENST2\t{STOP}\t{IN_FRAME_STOP}\n"
        );
        fs::write(&table, &input).unwrap();
        let paths = Some((map.to_str().unwrap(), intolerant.to_str().unwrap()));
        process_variants_table(table.to_str().unwrap(), paths).unwrap();
        assert_eq!(
            fs::read_to_string(&table).unwrap(),
            format!(
                "Feature\t5UTR_consequence\t5UTR_annotation\t5UTR_lof\t5UTR_len_changing\t5UTR_frameshift\t5UTR_span_intol_domain\nENST1\t{STOP}\t{IN_FRAME_STOP}\tTrue\tTrue\tFalse\tTrue\nENST2\t{STOP}\t{IN_FRAME_STOP}\tFalse\tTrue\tFalse\tFalse\n"
            )
        );
        // A missing, foreign or damaged cache leaves the table untouched.
        fs::write(&table, &input).unwrap();
        let missing = folder.join("absent.pkl");
        let gone = Some((missing.to_str().unwrap(), intolerant.to_str().unwrap()));
        assert!(process_variants_table(table.to_str().unwrap(), gone).is_err());
        fs::write(&map, b"not a pickle").unwrap();
        assert!(process_variants_table(table.to_str().unwrap(), paths).is_err());
        assert_eq!(fs::read_to_string(&table).unwrap(), input);
        let _ = fs::remove_dir_all(&folder);
    }
}
