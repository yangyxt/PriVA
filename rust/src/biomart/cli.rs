use super::*;
use crate::python_args::{parse_long_option_occurrences, require_options};
use crate::python_path::path_value;
use std::collections::HashMap;
const HELP: &str = "usage: priva biomart PROCEDURE --helper-root SOURCE [named options]\n\
Procedures: identity, receipt-check, receipt-write, artifact-route, artifact-receipt,\n\
            changed-target, verify, validate-table, parse-metadata, source,\n\
            snapshot-identity, is-bundled, ftp-specification, ftp-table-identity,\n\
            ftp-verify-record, ftp-parse-xml, ftp-sqlite-info, ftp-schema,\n\
            ftp-load, ftp-audit, ftp-export, ftp-project, ftp-verify-raw, ftp-metadata,\n\
            ftp-make-provenance, ftp-portable-provenance, external-path,\n\
            atomic-json and runtime-receipt (--method service|ftp), fetch.\n\
Use --mapping, --artifact, repeated --input, --assembly and --kind where applicable.\n\
artifact-route and changed-target accept --artifact-format {pickle,rkyv} (default pickle;\n\
rkyv names the natively built PM1 regions archive).\n\
fetch --source SOURCE --output TSV [--raw-output RAW] [--resource-dir DIR] [--timeout S]\n\
[--attempts N] acquires one pinned source through curl; a completed output is reused offline.\n\
The explicit helper root supplies pinned source assets and canonical bundle context;\n\
the optional native implementation never runs Python or selects a scientific release implicitly.";
pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args
        .first()
        .is_some_and(|s| matches!(s.as_str(), "-h" | "--help"))
    {
        println!("{HELP}");
        return Ok(());
    }
    let command = args.first().ok_or_else(|| invalid(HELP))?.as_str();
    let (names, flags, required): (&[&str], &[&str], &[&str]) = match command {
        "identity" => (
            &["helper-root", "mapping", "input"],
            &[],
            &["helper-root", "mapping"],
        ),
        "receipt-check" => (
            &["helper-root", "artifact", "mapping", "input"],
            &[],
            &["helper-root", "artifact", "mapping"],
        ),
        "receipt-write" => (
            &[
                "helper-root",
                "artifact",
                "mapping",
                "input",
                "expected-identity",
            ],
            &[],
            &["helper-root", "artifact", "mapping", "expected-identity"],
        ),
        "artifact-route" => (
            &[
                "helper-root",
                "artifact",
                "mapping",
                "cache-root",
                "output-dir",
                "assembly",
                "kind",
                "artifact-format",
                "input",
            ],
            &["print-path"],
            &[
                "helper-root",
                "artifact",
                "mapping",
                "cache-root",
                "output-dir",
                "assembly",
                "kind",
            ],
        ),
        "changed-target" => (
            &[
                "helper-root",
                "identity",
                "cache-root",
                "output-dir",
                "assembly",
                "kind",
                "artifact-format",
            ],
            &[],
            &[
                "helper-root",
                "identity",
                "cache-root",
                "output-dir",
                "assembly",
                "kind",
            ],
        ),
        "artifact-receipt" => (
            &["helper-root", "artifact"],
            &[],
            &["helper-root", "artifact"],
        ),
        "verify" | "ftp-verify-record" => (
            &["helper-root", "mapping", "assembly", "record"],
            &[],
            &["helper-root", "mapping"],
        ),
        "validate-table" => (
            &["helper-root", "table", "output"],
            &["raw"],
            &["helper-root", "table"],
        ),
        "parse-metadata" => (
            &["helper-root", "registry", "datasets", "source"],
            &[],
            &["helper-root", "registry", "datasets"],
        ),
        "source" => (
            &[
                "helper-root",
                "source",
                "basename",
                "assembly",
                "vep-release",
            ],
            &[],
            &["helper-root"],
        ),
        "snapshot-identity" => (
            &["helper-root", "mapping", "source"],
            &[],
            &["helper-root", "mapping"],
        ),
        "is-bundled" | "ftp-table-identity" => (
            &["helper-root", "mapping"],
            &[],
            &["helper-root", "mapping"],
        ),
        "ftp-specification" | "ftp-sqlite-info" => (&["helper-root"], &[], &["helper-root"]),
        "ftp-make-provenance" => (
            &["helper-root", "mapping", "operation"],
            &[],
            &["helper-root", "mapping", "operation"],
        ),
        "ftp-portable-provenance" => (&["helper-root", "record"], &[], &["helper-root", "record"]),
        "external-path" => (
            &["helper-root", "path", "method"],
            &[],
            &["helper-root", "path", "method"],
        ),
        "atomic-json" => (
            &["helper-root", "output", "record", "method"],
            &[],
            &["helper-root", "output", "record", "method"],
        ),
        "runtime-receipt" => (
            &["helper-root", "workspace", "basename", "record", "method"],
            &[],
            &["helper-root", "workspace", "basename", "record", "method"],
        ),
        "fetch" => (
            &[
                "helper-root",
                "source",
                "output",
                "raw-output",
                "resource-dir",
                "timeout",
                "attempts",
            ],
            &[],
            &["helper-root", "source", "output"],
        ),
        "ftp-schema" => (&["helper-root", "schema"], &[], &["helper-root", "schema"]),
        "ftp-verify-raw" => (
            &["helper-root", "path", "pin", "name"],
            &[],
            &["helper-root", "path"],
        ),
        "ftp-metadata" => (
            &["helper-root", "raw", "schemas", "output-dir"],
            &[],
            &["helper-root", "raw", "schemas", "output-dir"],
        ),
        "ftp-load" => (
            &[
                "helper-root",
                "db",
                "raw",
                "alias",
                "table",
                "columns",
                "schema",
            ],
            &[],
            &[
                "helper-root",
                "db",
                "raw",
                "alias",
                "table",
                "columns",
                "schema",
            ],
        ),
        "ftp-audit" => (&["helper-root", "db"], &[], &["helper-root", "db"]),
        "ftp-export" => (
            &["helper-root", "db", "output", "audit"],
            &[],
            &["helper-root", "db", "output", "audit"],
        ),
        "ftp-project" => (
            &["helper-root", "db", "raw", "schema", "output"],
            &[],
            &["helper-root", "db", "raw", "schema", "output"],
        ),
        "ftp-parse-xml" => (
            &["helper-root", "xml", "dataset-row"],
            &[],
            &["helper-root", "xml", "dataset-row"],
        ),
        _ => {
            return Err(invalid(format!(
                "Unknown BioMart procedure: {command}\n{HELP}"
            )));
        }
    };
    let Some(occurrences) =
        parse_long_option_occurrences(&args[1..], names, flags).map_err(invalid)?
    else {
        println!("{HELP}");
        return Ok(());
    };
    let values: HashMap<String, String> = occurrences.iter().cloned().collect();
    require_options(&values, required).map_err(invalid)?;
    if command == "ftp-verify-record" && !values.contains_key("record") {
        return Err(invalid("the following arguments are required: --record"));
    }
    if command == "ftp-verify-raw" && option_count(&values, &["pin", "name"]) != 1 {
        return Err(invalid("provide exactly one of --pin or --name"));
    }
    let artifact_format = values
        .get("artifact-format")
        .map_or("pickle", String::as_str);
    if !matches!(artifact_format, "pickle" | "rkyv") {
        return Err(invalid(format!(
            "argument --artifact-format: invalid choice: '{artifact_format}' \
             (choose from 'pickle', 'rkyv')"
        )));
    }
    let option = |key: &str| values.get(key).map(String::as_str);
    if option("method").is_some_and(|m| !matches!(m, "service" | "ftp")) {
        return Err(invalid("--method must be service or ftp"));
    }
    let fetch_options = if command == "fetch" {
        Some(fetch::Options::admit(&values)?)
    } else {
        None
    };
    let p = |key: &str| path_value(option(key).expect("required option"));
    let inputs = occurrences
        .iter()
        .filter(|(k, _)| k == "input")
        .map(|(_, v)| path_value(v))
        .collect::<Vec<_>>();
    crate::builder_error::finish((|| {
        let context = Context::new(&p("helper-root"))?;
        let value = match command {
            "identity" => Some(mapping_identity(&context, &p("mapping"), &inputs)?),
            "receipt-check" => {
                check_receipt(&context, &p("artifact"), &p("mapping"), &inputs)?;
                None
            }
            "receipt-write" => {
                identity::write_receipt(
                    &context,
                    &p("artifact"),
                    &p("mapping"),
                    &inputs,
                    option("expected-identity").unwrap(),
                )?;
                None
            }
            "artifact-receipt" => {
                // This procedure exposes the original function's complete
                // JSON value. Keep wide integers, lone-surrogate escapes,
                // duplicate-key last-wins meaning and deeply nested extras.
                // Lexical whitespace is retained, not claimed to be dumps().
                println!("{}", receipt::artifact_receipt(&p("artifact"))?.json);
                None
            }
            "changed-target" => Some(string(
                identity::changed_target(
                    &parse_json(option("identity").unwrap())?,
                    &p("cache-root"),
                    &p("output-dir"),
                    option("assembly").unwrap(),
                    option("kind").unwrap(),
                    artifact_format,
                )?
                .display()
                .to_string(),
            )),
            "artifact-route" => {
                let path = artifact_route(
                    &context,
                    &p("artifact"),
                    &p("mapping"),
                    &p("cache-root"),
                    &p("output-dir"),
                    option("assembly").unwrap(),
                    option("kind").unwrap(),
                    artifact_format,
                    &inputs,
                )?;
                if values.contains_key("print-path") {
                    println!("{}", path.display());
                    None
                } else {
                    Some(object([("artifact", string(path.display().to_string()))]))
                }
            }
            "verify" => Some(verify::pinned(&p("mapping"), option("assembly"))?),
            "ftp-verify-record" => Some(verify::ftp_record(
                &p("mapping"),
                read_json(&p("record"))?,
                option("assembly"),
            )?),
            "validate-table" => Some(table::validate(
                &p("table"),
                values.contains_key("raw"),
                option("output").map(path_value).as_deref(),
            )?),
            "parse-metadata" => {
                let read = |path: &Path| -> Result<String> {
                    let bytes = fs::read(path).map_err(|e| os_error(e, path, None))?;
                    Ok(crate::python_io::decode_utf8(&bytes)
                        .map_err(|e| error("UnicodeDecodeError", e))?
                        .replace("\r\n", "\n")
                        .replace('\r', "\n"))
                };
                Some(verify::parse_metadata(
                    &read(&p("registry"))?,
                    &read(&p("datasets"))?,
                    &string(option("source").unwrap_or(GRCH37)),
                )?)
            }
            "source" => {
                let selected = if let Some(assembly) = option("assembly") {
                    source_for_assembly(assembly, option("vep-release").map(string).as_ref())?
                } else if let Some(basename) = option("basename") {
                    source_for_basename(basename).ok_or_else(|| {
                        error(
                            "ValueError",
                            format!("Unknown pinned BioMart filename: {basename}"),
                        )
                    })?
                } else {
                    source(&string(option("source").unwrap_or(GRCH37)))?
                };
                Some(selected.clone())
            }
            "snapshot-identity" => Some(snapshot_identity(
                &p("mapping"),
                &string(option("source").unwrap_or(GRCH37)),
            )?),
            "is-bundled" => Some(Json::Bool(identity::is_bundled(&context, &p("mapping"))?)),
            "ftp-table-identity" => Some(verify::ftp_table_identity(&p("mapping"))?),
            "ftp-specification" => Some(get(facts(), "ftp_specification").clone()),
            "ftp-make-provenance" => {
                let (path, sha256) = publication::helper_producer(&context)?;
                Some(publication::make_provenance(
                    &p("mapping"),
                    option("operation").unwrap(),
                    &path,
                    &sha256,
                )?)
            }
            "ftp-portable-provenance" => {
                let record =
                    publication_json::Record::parse(&read_json_text(&p("record"))?)?.portable()?;
                println!("{}", record.text()?);
                None
            }
            "external-path" => {
                publication::external_path(&context, &p("path"), option("method") == Some("ftp"))?;
                None
            }
            "atomic-json" => {
                let record = publication_json::Record::parse(&read_json_text(&p("record"))?)?;
                publication::atomic(&p("output"), &record, option("method") == Some("ftp"))?;
                None
            }
            "runtime-receipt" => {
                let record = publication_json::Record::parse(&read_json_text(&p("record"))?)?;
                Some(string(
                    publication::runtime(
                        &context,
                        &record,
                        &p("workspace"),
                        option("basename").unwrap(),
                        option("method") == Some("ftp"),
                    )?
                    .display()
                    .to_string(),
                ))
            }
            "fetch" => {
                fetch::acquire(&context, fetch_options.as_ref().unwrap())?;
                None
            }
            "ftp-sqlite-info" => Some(ftp::sqlite_information()),
            "ftp-schema" => Some(ftp::schema(&p("schema"))?),
            "ftp-verify-raw" => Some(ftp::verify_raw(
                &p("path"),
                &match option("pin") {
                    Some(value) => parse_json(value)?,
                    None => ftp::named_pin(option("name").unwrap())?,
                },
            )?),
            "ftp-metadata" => Some(ftp::publisher(
                &p("raw"),
                &read_json(&p("schemas"))?,
                &p("output-dir"),
            )?),
            "ftp-load" => Some(ftp::load(
                &p("db"),
                &p("raw"),
                option("alias").unwrap(),
                option("table").unwrap(),
                &parse_json(option("columns").unwrap())?,
                &read_json(&p("schema"))?,
            )?),
            "ftp-audit" => Some(ftp::audit(&p("db"))?),
            "ftp-export" => Some(ftp::export(
                &p("db"),
                &p("output"),
                &read_json(&p("audit"))?,
            )?),
            "ftp-project" => Some(ftp::project(
                &p("raw"),
                &p("schema"),
                &p("db"),
                &p("output"),
            )?),
            "ftp-parse-xml" => Some(verify::parse_xml(
                &fs::read(p("xml")).map_err(|e| os_error(e, &p("xml"), None))?,
                &parse_json(option("dataset-row").unwrap())?,
            )?),
            _ => unreachable!(),
        };
        if let Some(value) = value {
            let sorted = !matches!(
                command,
                "ftp-schema"
                    | "ftp-load"
                    | "ftp-audit"
                    | "ftp-export"
                    | "ftp-project"
                    | "ftp-verify-raw"
                    | "ftp-metadata"
            );
            println!("{}", json_text::dumps(&value, sorted).map_err(invalid)?)
        }
        Ok(())
    })())
}
fn option_count(values: &HashMap<String, String>, names: &[&str]) -> usize {
    names
        .iter()
        .filter(|name| values.contains_key(**name))
        .count()
}
pub(super) fn source_for_assembly(assembly: &str, release: Option<&Json>) -> Result<&'static Json> {
    if let Some(release) = release {
        require(
            !matches!(release, Json::Bool(_)) && matches!(py_str(release).as_str(), "113" | "116"),
            format!(
                "No pinned automatic BioMart mapping for VEP {}; supply an explicit assembly-compatible nm_to_enst_tsv",
                py_str(release)
            ),
        )?;
    }
    for selected in &array(get(facts(), "sources"))[..2] {
        if array(get(selected, "aliases"))
            .iter()
            .any(|a| a.as_str() == Some(assembly))
        {
            return if text(get(selected, "source_id")) != GRCH37
                && release.is_some_and(|v| py_str(v) == "116")
            {
                source(&string(FTP))
            } else {
                Ok(selected)
            };
        }
    }
    Err(error(
        "ValueError",
        format!("No explicitly pinned BioMart source for assembly: {assembly}"),
    ))
}
