//! Shared text-table operations and the two existing builder publication rules.
//! TSV bytes use the pandas-compatible writer; gzip metadata follows the actual
//! staged filename and current time, as `DataFrame.to_csv(path)` does. DEFLATE
//! bytes may differ between backends; decompressed TSV bytes must match.
use crate::cli::{Result, invalid};
use crate::gzip_writer::{GzipOptions, GzipWriter};
use crate::json_text::{self, Json};
use crate::pandas_table::{self, Cell, Column, Table};
use crate::python_text::python_strip;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

static NEXT: AtomicU64 = AtomicU64::new(0);

pub(crate) fn read(path: &Path, comments: bool) -> Result<Table> {
    pandas_table::read_builder_table(path, comments)
}

pub(crate) fn text<'a>(table: &'a Table, name: &str, row: usize) -> &'a str {
    match table
        .column(name)
        .expect("required column checked")
        .cell(row)
    {
        Cell::Text(value) => value,
        Cell::Missing => "nan", // pandas astype(str), including an unmatched merge.
        _ => unreachable!("builder inputs use dtype=str, keep_default_na=False"),
    }
}

pub(crate) fn optional(value: &str) -> String {
    let value = python_strip(value);
    if value.is_empty() {
        "-".to_owned()
    } else {
        value.to_owned()
    }
}

pub(crate) fn require(table: &Table, required: &[&str], label: &str) -> Result<()> {
    let mut missing: Vec<&str> = required
        .iter()
        .copied()
        .filter(|name| !table.has_column(name))
        .collect();
    missing.sort_unstable();
    if !missing.is_empty() {
        return Err(invalid(format!(
            "{label} is missing columns: {}",
            repr_list(&missing)
        )));
    }
    Ok(())
}

pub(crate) fn repr_list(values: &[&str]) -> String {
    // These diagnostics contain source identifiers, never serialized TSV data.
    format!(
        "[{}]",
        values
            .iter()
            .map(|value| pandas_table::python_text_repr(value))
            .collect::<Vec<_>>()
            .join(", ")
    )
}

pub(crate) fn select(table: &Table, names: &[&str]) -> Result<Table> {
    Ok(Table {
        names: names.iter().map(|name| (*name).to_owned()).collect(),
        columns: names
            .iter()
            .map(|name| {
                table
                    .column(name)
                    .cloned()
                    .map_err(|error| crate::builder_error::category("KeyError", error))
            })
            .collect::<Result<_>>()?,
        rows: table.rows,
    })
}

pub(crate) fn strings(names: &[&str], rows: Vec<Vec<String>>) -> Table {
    Table {
        names: names.iter().map(|name| (*name).to_owned()).collect(),
        rows: rows.len(),
        columns: (0..names.len())
            .map(|column| Column::Text(rows.iter().map(|row| Some(row[column].clone())).collect()))
            .collect(),
    }
}

pub(crate) fn json_counts(values: impl IntoIterator<Item = (&'static str, usize)>) -> Json {
    Json::Object(
        values
            .into_iter()
            .map(|(name, count)| (name.to_owned(), Json::Int(count as i128)))
            .collect(),
    )
}

pub(crate) fn print_summary(value: &Json) -> Result<()> {
    println!("{}", json_text::dumps(value, false).map_err(invalid)?);
    Ok(())
}

#[derive(Clone, Copy)]
pub(crate) enum Publication {
    /// collapse_HPO_anno: private directory, original basename, prior file mode.
    Hpo,
    /// build_mondo_disease_scope: sibling .tmp[.gz], ordinary pandas file mode.
    Mondo,
}

fn parent(path: &Path) -> &Path {
    path.parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or(Path::new("."))
}

fn gzip_options(path: &Path) -> Result<GzipOptions> {
    let name = path
        .file_name()
        .ok_or_else(|| invalid("output requires a filename"))?
        .to_string_lossy();
    // GzipFile encodes the basename as Latin-1, omitting FNAME if not encodable.
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    let filename: Option<Vec<u8>> = name
        .chars()
        .map(|ch| u8::try_from(ch as u32).ok())
        .collect();
    let mtime = u32::try_from(SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs())?;
    Ok(GzipOptions {
        mtime,
        filename,
        ..GzipOptions::default()
    })
}

fn write_staged(table: &Table, path: &Path) -> Result<()> {
    let name = path.to_string_lossy().to_lowercase();
    if [".bz2", ".xz", ".zip", ".zst", ".tar"]
        .iter()
        .any(|suffix| name.ends_with(suffix))
    {
        return Err(invalid(
            "this builder writer supports gzip and plain TSV destinations",
        ));
    }
    let file = BufWriter::new(File::create(path)?);
    if name.ends_with(".gz") {
        let mut writer = GzipWriter::new(file, gzip_options(path)?)?;
        pandas_table::write_table_to(&mut writer, table)?;
        writer.finish()?;
    } else {
        let mut writer = file;
        pandas_table::write_table_to(&mut writer, table)?;
        writer.flush()?;
    }
    Ok(())
}

pub(crate) fn publish(table: &Table, output: &Path, policy: Publication) -> Result<()> {
    publish_with(output, policy, |path| write_staged(table, path))
}

fn publish_with(
    output: &Path,
    policy: Publication,
    write: impl FnOnce(&Path) -> Result<()>,
) -> Result<()> {
    fs::create_dir_all(parent(output))?;
    let name = output
        .file_name()
        .ok_or_else(|| invalid("output requires a filename"))?
        .to_string_lossy();
    let (staged, directory, mode): (PathBuf, Option<PathBuf>, _) = match policy {
        Publication::Hpo => {
            let mode = match fs::metadata(output) {
                Ok(metadata) => Some(metadata.permissions()),
                Err(error) if error.kind() == io::ErrorKind::NotFound => None,
                Err(error) => return Err(error.into()),
            };
            let directory = loop {
                let nonce = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
                let directory = parent(output).join(format!(
                    ".{name}.write-{}-{nonce}-{}",
                    std::process::id(),
                    NEXT.fetch_add(1, Ordering::Relaxed)
                ));
                let mut builder = fs::DirBuilder::new();
                #[cfg(unix)]
                {
                    use std::os::unix::fs::DirBuilderExt;
                    builder.mode(0o700);
                }
                match builder.create(&directory) {
                    Ok(()) => break directory,
                    Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
                    Err(error) => return Err(error.into()),
                }
            };
            (
                directory.join(output.file_name().unwrap()),
                Some(directory),
                mode,
            )
        }
        Publication::Mondo => {
            let suffix = if name.ends_with(".gz") {
                ".tmp.gz"
            } else {
                ".tmp"
            };
            (
                parent(output).join(format!(".{name}.{}{suffix}", std::process::id())),
                None,
                None,
            )
        }
    };
    let result = (|| {
        write(&staged)?;
        if let Some(mode) = mode {
            fs::set_permissions(&staged, mode)?;
        }
        fs::rename(&staged, output)?;
        Ok(())
    })();
    // Clean only this invocation's staging path. Never touch earlier attempts.
    if let Some(directory) = directory {
        fs::remove_dir_all(directory)?;
    } else if let Err(error) = fs::remove_file(staged) {
        if error.kind() != io::ErrorKind::NotFound {
            return Err(error.into());
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn partial_writes_and_rename_failures_do_not_replace_complete_outputs() {
        let root = std::env::temp_dir().join(format!(
            "priva_builder_publish_{}_{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&root).unwrap();
        for policy in [Publication::Hpo, Publication::Mondo] {
            let target = root.join("output.tsv.gz");
            fs::write(&target, b"complete").unwrap();
            let result = publish_with(&target, policy, |path| {
                fs::write(path, b"partial")?;
                Err(io::Error::other("injected writer failure").into())
            });
            assert!(result.is_err());
            assert_eq!(fs::read(&target).unwrap(), b"complete");
            fs::remove_file(&target).unwrap();
            fs::create_dir(&target).unwrap();
            assert!(
                publish_with(&target, policy, |path| {
                    fs::write(path, b"new")?;
                    Ok(())
                })
                .is_err()
            );
            assert!(target.is_dir());
            fs::remove_dir(target).unwrap();
            assert_eq!(fs::read_dir(&root).unwrap().count(), 0);
        }
        fs::remove_dir(root).unwrap();
    }
}
