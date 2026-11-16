//! Joins the splicing and 5' UTR copies of a family table, as if the two
//! interpretation steps had run one after another on one table: the
//! combined table's columns, then the splicing step's new columns, then the
//! 5' UTR step's new columns.
//!
//! Both steps read and rewrite the table with the same reader and writer
//! (`pandas_table`), so their copies of the combined table's columns must
//! match cell for cell. The join refuses copies whose row counts differ,
//! whose shared columns differ, or that both added a column of one name.
//! The output is the table the same writer writes, as the 5' UTR step wrote
//! the table when it ran after the splicing step.
//!
//! Memory: the join reads the three files one record at a time, so it holds
//! one record of each file and the columns the two steps added. The combined
//! table's columns are copied as they are written: the splicing step read
//! and wrote them, and the reader gives such text back unchanged. The added
//! columns are read and written again as the whole-table join did, because
//! the reader may give them another type than the step held (a column of
//! exon numbers such as "4" is read as numbers and written as "4.0"). When a
//! header is not plain writer output (a name the reader would rename, or a
//! quoted name) or a record has another number of fields, the join reads the
//! whole tables instead (`join_in_memory`), as it did before.
use super::{Result, refuse};
use crate::pandas_table::{Table, parse_table, read_table, write_table, write_table_to};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Cursor, Write};
use std::path::{Path, PathBuf};

fn project(table: &Table, names: &[String]) -> Result<Vec<u8>> {
    let mut columns = Vec::new();
    for name in names {
        let Some(index) = table.column_index(name) else {
            return refuse(format!("A step dropped the column {name}"));
        };
        columns.push(table.columns[index].clone());
    }
    let projected = Table {
        names: names.to_vec(),
        columns,
        rows: table.rows,
    };
    let mut text = Vec::new();
    write_table_to(&mut text, &projected)?;
    Ok(text)
}

/// The join over whole tables in memory: about seven times each table's size.
fn join_in_memory(combined: &Path, splicing: &Path, utr: &Path, output: &Path) -> Result<()> {
    let base = read_table(combined)?;
    let mut first = read_table(splicing)?;
    let mut second = read_table(utr)?;
    if first.rows != base.rows || second.rows != base.rows {
        return refuse(format!(
            "Row counts differ: combined {}, splicing {}, 5' UTR {}",
            base.rows, first.rows, second.rows
        ));
    }
    let names = base.names;
    drop(base.columns);
    if project(&first, &names)? != project(&second, &names)? {
        return refuse("The splicing and 5' UTR copies differ in the combined table's columns");
    }
    for (index, name) in second.names.iter().enumerate() {
        if names.contains(name) {
            continue;
        }
        if first.column_index(name).is_some() {
            return refuse(format!(
                "Both the splicing and the 5' UTR step added the column {name}"
            ));
        }
        first.names.push(name.clone());
        first.columns.push(std::mem::replace(
            &mut second.columns[index],
            crate::pandas_table::Column::Int(Vec::new()),
        ));
    }
    write_table(output, &first)
}

/// One file read record by record. A record ends at a line end outside
/// quotes; the writer quotes a field with a tab, quote or line end and
/// doubles its quotes, so a field never holds an unquoted quote.
struct Records<R: BufRead> {
    reader: R,
    record: Vec<u8>,
    count: usize,
}

fn open(path: &Path) -> Result<Records<BufReader<File>>> {
    Ok(Records::new(BufReader::with_capacity(
        1 << 22,
        File::open(path)?,
    )))
}

impl<R: BufRead> Records<R> {
    fn new(reader: R) -> Self {
        Records {
            reader,
            record: Vec::new(),
            count: 0,
        }
    }

    /// The next record without its line end; `false` at the end of the file.
    fn next(&mut self) -> Result<bool> {
        self.record.clear();
        let mut quotes = 0usize;
        loop {
            let start = self.record.len();
            if self.reader.read_until(b'\n', &mut self.record)? == 0 {
                if self.record.is_empty() {
                    return Ok(false);
                }
                break;
            }
            quotes += self.record[start..].iter().filter(|b| **b == b'"').count();
            if quotes % 2 == 0 {
                break;
            }
        }
        if self.record.last() == Some(&b'\n') {
            self.record.pop();
        }
        if self.record.is_empty() {
            // pandas skips a blank line; the writer never writes one.
            return self.next();
        }
        self.count += 1;
        Ok(true)
    }
}

/// The fields of a record as they are written, split at tabs outside quotes.
fn fields(record: &[u8]) -> Vec<&[u8]> {
    let mut out = Vec::new();
    let (mut start, mut quoted) = (0, false);
    for (index, byte) in record.iter().enumerate() {
        match byte {
            b'"' => quoted = !quoted,
            b'\t' if !quoted => {
                out.push(&record[start..index]);
                start = index + 1;
            }
            _ => {}
        }
    }
    out.push(&record[start..]);
    out
}

/// The column names the reader gives a header record, when they are the
/// header's fields as written; `None` sends the join to `join_in_memory`.
fn plain_names(header: &[u8]) -> Result<Option<Vec<String>>> {
    let text = std::str::from_utf8(header).map_err(|_| "a table header is not UTF-8 text")?;
    let names = parse_table(&format!("{text}\n"))?.names;
    let written: Vec<&[u8]> = fields(header);
    let plain = names.len() == written.len()
        && names
            .iter()
            .zip(&written)
            .all(|(name, field)| name.as_bytes() == *field);
    Ok(plain.then_some(names))
}

/// The columns at `positions` (ascending) of a written table as the reader
/// and writer give them when the table is read again: the written text, one
/// record per row, header first. A constant first column keeps every line
/// non-blank and the table wider than one column, so the writer quotes an
/// empty value as it does among other columns; callers skip that field.
/// `None` when a record has another number of fields than the header.
fn rewritten(path: &Path, positions: &[usize], width: usize) -> Result<Option<Vec<u8>>> {
    let mut records = open(path)?;
    let mut text = Vec::new();
    while records.next()? {
        let fields = fields(&records.record);
        if fields.len() != width {
            return Ok(None);
        }
        text.extend_from_slice(if records.count == 1 {
            b"priva_join_row"
        } else {
            b"r"
        });
        for position in positions {
            text.push(b'\t');
            text.extend_from_slice(fields[*position]);
        }
        text.push(b'\n');
    }
    let text = String::from_utf8(text).map_err(|_| "a table is not UTF-8 text")?;
    let table = parse_table(&text)?;
    drop(text);
    let mut out = Vec::new();
    write_table_to(&mut out, &table)?;
    Ok(Some(out))
}

/// A temporary file next to `output`, renamed over it once the join is done.
fn staged(output: &Path) -> PathBuf {
    let mut name = output.file_name().unwrap_or_default().to_os_string();
    name.push(format!(".join.{}.partial", std::process::id()));
    output.with_file_name(name)
}

/// Joins `splicing` and `utr`, both made from `combined`, into `output`
/// (which may be `combined`). An empty combined table (no variants left
/// after filtration) stays empty.
pub(crate) fn join(combined: &Path, splicing: &Path, utr: &Path, output: &Path) -> Result<()> {
    if fs::metadata(combined)?.len() == 0 {
        if output != combined {
            fs::write(output, b"")?;
        }
        return Ok(());
    }
    let mut base = open(combined)?;
    let mut first = open(splicing)?;
    let mut second = open(utr)?;
    let headers = (base.next()?, first.next()?, second.next()?);
    if headers != (true, true, true) {
        // An empty copy: the whole-table reader gives pandas' error.
        return join_in_memory(combined, splicing, utr, output);
    }
    let (Some(base_names), Some(first_names), Some(second_names)) = (
        plain_names(&base.record)?,
        plain_names(&first.record)?,
        plain_names(&second.record)?,
    ) else {
        return join_in_memory(combined, splicing, utr, output);
    };
    // The refusals keep the earlier order: row counts, a dropped column,
    // differing copies, then a column both steps added.
    let position = |names: &[String], name: &String| names.iter().position(|n| n == name);
    // As `project`: the first name missing from the splicing copy, else the
    // first name missing from the 5' UTR copy.
    let missing = |names: &[String]| base_names.iter().find(|n| !names.contains(n)).cloned();
    let dropped = missing(&first_names).or_else(|| missing(&second_names));
    let shared: Vec<(usize, usize)> = base_names
        .iter()
        .filter_map(|name| {
            Some((
                position(&first_names, name)?,
                position(&second_names, name)?,
            ))
        })
        .collect();
    let mut both_added = None;
    let mut added = Vec::new();
    for (index, name) in second_names.iter().enumerate() {
        if base_names.contains(name) {
            continue;
        }
        if first_names.contains(name) {
            both_added.get_or_insert_with(|| name.clone());
        }
        added.push(index);
    }
    // The columns the splicing step added, at their places in its table.
    let first_added: Vec<usize> = (0..first_names.len())
        .filter(|index| !base_names.contains(&first_names[*index]))
        .collect();
    let (Some(first_new), Some(second_new)) = (
        rewritten(splicing, &first_added, first_names.len())?,
        rewritten(utr, &added, second_names.len())?,
    ) else {
        return join_in_memory(combined, splicing, utr, output);
    };
    let mut first_new = Records::new(Cursor::new(first_new));
    let mut second_new = Records::new(Cursor::new(second_new));
    first_new.next()?;
    second_new.next()?;
    let path = staged(output);
    let result = (|| -> Result<bool> {
        let mut out = BufWriter::with_capacity(1 << 22, File::create(&path)?);
        let mut differ = false;
        // One output record: the splicing copy's fields, its added ones as
        // read again, then the 5' UTR step's added fields as read again.
        let write = |first: &[u8],
                     first_new: &[u8],
                     second_new: &[u8],
                     out: &mut BufWriter<File>|
         -> Result<()> {
            let mut again = fields(first_new).into_iter().skip(1);
            for (index, field) in fields(first).into_iter().enumerate() {
                if index > 0 {
                    out.write_all(b"\t")?;
                }
                if first_added.binary_search(&index).is_ok() {
                    out.write_all(again.next().unwrap_or_default())?;
                } else {
                    out.write_all(field)?;
                }
            }
            for field in fields(second_new).into_iter().skip(1) {
                out.write_all(b"\t")?;
                out.write_all(field)?;
            }
            out.write_all(b"\n")?;
            Ok(())
        };
        // The header as written: the names are plain, so the reader keeps them.
        out.write_all(&first.record)?;
        for index in &added {
            out.write_all(b"\t")?;
            out.write_all(fields(&second.record)[*index])?;
        }
        out.write_all(b"\n")?;
        loop {
            let more = (base.next()?, first.next()?, second.next()?);
            match more {
                (false, false, false) => break,
                (_, true, true) if !differ && dropped.is_none() => {
                    let (a, b) = (fields(&first.record), fields(&second.record));
                    if shared.iter().any(|(i, j)| a[*i] != b[*j]) {
                        differ = true;
                    } else if first_new.next()? && second_new.next()? {
                        write(
                            &first.record,
                            &first_new.record,
                            &second_new.record,
                            &mut out,
                        )?;
                    } else {
                        return Ok(false);
                    }
                }
                _ => {}
            }
        }
        let counts = (base.count - 1, first.count - 1, second.count - 1);
        if counts.1 != counts.0 || counts.2 != counts.0 {
            return Err(format!(
                "Row counts differ: combined {}, splicing {}, 5' UTR {}",
                counts.0, counts.1, counts.2
            )
            .into());
        }
        if let Some(name) = &dropped {
            return Err(format!("A step dropped the column {name}").into());
        }
        if differ {
            return Err(
                "The splicing and 5' UTR copies differ in the combined table's columns".into(),
            );
        }
        if let Some(name) = &both_added {
            return Err(
                format!("Both the splicing and the 5' UTR step added the column {name}").into(),
            );
        }
        out.flush()?;
        out.get_ref().sync_all()?;
        Ok(true)
    })();
    match result {
        Ok(true) => {
            fs::rename(&path, output)?;
            Ok(())
        }
        Ok(false) => {
            let _ = fs::remove_file(&path);
            join_in_memory(combined, splicing, utr, output)
        }
        Err(error) => {
            let _ = fs::remove_file(&path);
            Err(error)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Folder(std::path::PathBuf);
    impl Drop for Folder {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }
    fn folder(name: &str) -> Folder {
        let path =
            std::env::temp_dir().join(format!("priva_table_join_{name}_{}", std::process::id()));
        let _ = fs::remove_dir_all(&path);
        fs::create_dir_all(&path).unwrap();
        Folder(path)
    }

    #[test]
    fn splicing_columns_come_before_utr_columns() {
        let f = folder("order");
        let p = |n: &str| f.0.join(n);
        fs::write(p("c.tsv"), "Chr\tPos\tScore\n1\t10\t0.5\n1\t20\t\n").unwrap();
        fs::write(
            p("s.tsv"),
            "Chr\tPos\tScore\tsplicing_lof\n1\t10\t0.5\tTrue\n1\t20\t\tFalse\n",
        )
        .unwrap();
        fs::write(
            p("u.tsv"),
            "Chr\tPos\tScore\t5UTR_lof\n1\t10\t0.5\tFalse\n1\t20\t\tTrue\n",
        )
        .unwrap();
        join(&p("c.tsv"), &p("s.tsv"), &p("u.tsv"), &p("out.tsv")).unwrap();
        assert_eq!(
            fs::read_to_string(p("out.tsv")).unwrap(),
            "Chr\tPos\tScore\tsplicing_lof\t5UTR_lof\n1\t10\t0.5\tTrue\tFalse\n1\t20\t\tFalse\tTrue\n"
        );
    }

    #[test]
    fn differing_copies_are_refused() {
        let f = folder("refuse");
        let p = |n: &str| f.0.join(n);
        fs::write(p("c.tsv"), "Chr\tPos\n1\t10\n").unwrap();
        fs::write(p("s.tsv"), "Chr\tPos\tx\n1\t10\t1\n").unwrap();
        fs::write(p("u.tsv"), "Chr\tPos\ty\n1\t11\t1\n").unwrap();
        let error = join(&p("c.tsv"), &p("s.tsv"), &p("u.tsv"), &p("o.tsv")).unwrap_err();
        assert!(error.to_string().contains("differ"), "{error}");
        fs::write(p("u.tsv"), "Chr\tPos\tx\n1\t10\t2\n").unwrap();
        let error = join(&p("c.tsv"), &p("s.tsv"), &p("u.tsv"), &p("o.tsv")).unwrap_err();
        assert!(error.to_string().contains("Both"), "{error}");
        fs::write(p("u.tsv"), "Chr\tPos\n1\t10\n1\t11\n").unwrap();
        assert!(join(&p("c.tsv"), &p("s.tsv"), &p("u.tsv"), &p("o.tsv")).is_err());
    }

    /// The table text the writer gives `text`, as a step's output would be.
    fn written(text: &str) -> String {
        let mut out = Vec::new();
        write_table_to(&mut out, &parse_table(text).unwrap()).unwrap();
        String::from_utf8(out).unwrap()
    }

    /// `base` with `extra` columns added, written by the writer. A value is
    /// one line of a one-column table, so a missing value is written "NA".
    fn with_columns(base: &str, extra: &[(&str, &[&str])]) -> String {
        let mut table = parse_table(base).unwrap();
        for (name, values) in extra {
            let column = parse_table(&format!("{name}\n{}\n", values.join("\n")))
                .unwrap()
                .columns
                .remove(0);
            table.names.push((*name).to_owned());
            table.columns.push(column);
        }
        let mut out = Vec::new();
        write_table_to(&mut out, &table).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn record_by_record_join_equals_the_whole_table_join() {
        let f = folder("stream");
        let p = |n: &str| f.0.join(n);
        // Integers with gaps, floats in several spellings, missing-value
        // spellings, bools, text with tabs, quotes and line ends, and a
        // column of empty values: the writer's output of each.
        let raw = "Chr\tPos\tAF\tN\tFlag\tNote\tEmpty\n\
                   1\t10\t1e-05\t3\tTrue\tplain\t\n\
                   1\t20\t0.10000000000000001\t\tFalse\t\"a\tb\"\t\n\
                   X\t30\tnan\t7\tTrue\t\"say \"\"hi\"\"\"\t\n\
                   Y\t40\t2.50\t-1\tFalse\t\"two\nlines\"\t\n\
                   MT\t50\tNA\t0\tTrue\tNA\t\n";
        let combined = written(raw);
        fs::write(p("c.tsv"), &combined).unwrap();
        fs::write(
            p("s.tsv"),
            with_columns(
                &combined,
                &[
                    ("splicing_lof", &["True", "False", "NA", "True", "False"]),
                    (
                        "splicing_affected_exons",
                        &["1,3", "NA", "\"4\t5\"", "2", "x"],
                    ),
                ],
            ),
        )
        .unwrap();
        fs::write(
            p("u.tsv"),
            with_columns(
                &combined,
                &[("5UTR_lof", &["False", "True", "False", "NA", "1.0"])],
            ),
        )
        .unwrap();
        join(&p("c.tsv"), &p("s.tsv"), &p("u.tsv"), &p("stream.tsv")).unwrap();
        join_in_memory(&p("c.tsv"), &p("s.tsv"), &p("u.tsv"), &p("memory.tsv")).unwrap();
        assert_eq!(
            fs::read(p("stream.tsv")).unwrap(),
            fs::read(p("memory.tsv")).unwrap()
        );
        // The output may be the combined table itself, as `priva task run` writes it.
        join(&p("c.tsv"), &p("s.tsv"), &p("u.tsv"), &p("c.tsv")).unwrap();
        assert_eq!(
            fs::read(p("c.tsv")).unwrap(),
            fs::read(p("memory.tsv")).unwrap()
        );
        assert!(fs::read_dir(&f.0).unwrap().all(|e| {
            !e.unwrap()
                .file_name()
                .to_string_lossy()
                .ends_with(".partial")
        }));
    }

    #[test]
    fn added_columns_are_read_again_like_the_whole_table_join() {
        let f = folder("again");
        let p = |n: &str| f.0.join(n);
        let combined = written("Chr\tPos\tAF\n1\t10\t0.5\n1\t20\t\n2\t30\t1e-05\n");
        fs::write(p("c.tsv"), &combined).unwrap();
        // The steps write their added values as they hold them: exon numbers
        // as text ("4"), which the reader takes as numbers ("4.0"), and an
        // offset column mixing "0" and "-3.0".
        let append = |extra: &[&str]| -> String {
            combined
                .lines()
                .zip(extra)
                .map(|(line, value)| format!("{line}\t{value}\n"))
                .collect()
        };
        fs::write(
            p("s.tsv"),
            append(&["splicing_affected_exons\tgain", "4\t0", "\t-3.0", "12\t0"]),
        )
        .unwrap();
        fs::write(p("u.tsv"), append(&["5UTR_exon", "7", "", "\"a\tb\""])).unwrap();
        join(&p("c.tsv"), &p("s.tsv"), &p("u.tsv"), &p("stream.tsv")).unwrap();
        join_in_memory(&p("c.tsv"), &p("s.tsv"), &p("u.tsv"), &p("memory.tsv")).unwrap();
        let memory = fs::read_to_string(p("memory.tsv")).unwrap();
        assert!(memory.contains("\t4.0\t0.0\t"), "{memory}");
        assert_eq!(fs::read_to_string(p("stream.tsv")).unwrap(), memory);
    }

    #[test]
    fn headers_the_reader_renames_use_the_whole_table_join() {
        let f = folder("rename");
        let p = |n: &str| f.0.join(n);
        // Two columns of one name: the reader names the second "Pos.1".
        fs::write(p("c.tsv"), "Chr\tPos\tPos\n1\t10\t11\n").unwrap();
        fs::write(p("s.tsv"), "Chr\tPos\tPos\tx\n1\t10\t11\t1\n").unwrap();
        fs::write(p("u.tsv"), "Chr\tPos\tPos\ty\n1\t10\t11\t2\n").unwrap();
        join(&p("c.tsv"), &p("s.tsv"), &p("u.tsv"), &p("stream.tsv")).unwrap();
        join_in_memory(&p("c.tsv"), &p("s.tsv"), &p("u.tsv"), &p("memory.tsv")).unwrap();
        assert_eq!(
            fs::read(p("stream.tsv")).unwrap(),
            fs::read(p("memory.tsv")).unwrap()
        );
    }

    #[test]
    fn refusals_keep_their_order() {
        let f = folder("order_of_refusals");
        let p = |n: &str| f.0.join(n);
        // Fewer rows and a dropped column: the row counts come first.
        fs::write(p("c.tsv"), "Chr\tPos\n1\t10\n1\t20\n").unwrap();
        fs::write(p("s.tsv"), "Chr\tx\n1\t1\n").unwrap();
        fs::write(p("u.tsv"), "Chr\tPos\ty\n1\t10\t2\n1\t20\t3\n").unwrap();
        let error = join(&p("c.tsv"), &p("s.tsv"), &p("u.tsv"), &p("o.tsv")).unwrap_err();
        assert!(
            error
                .to_string()
                .starts_with("Row counts differ: combined 2, splicing 1"),
            "{error}"
        );
        fs::write(p("s.tsv"), "Chr\tx\n1\t1\n1\t2\n").unwrap();
        let error = join(&p("c.tsv"), &p("s.tsv"), &p("u.tsv"), &p("o.tsv")).unwrap_err();
        assert_eq!(error.to_string(), "A step dropped the column Pos");
        assert!(!p("o.tsv").exists());
    }

    #[test]
    fn an_empty_table_stays_empty() {
        let f = folder("empty");
        let p = |n: &str| f.0.join(n);
        fs::write(p("c.tsv"), "").unwrap();
        join(&p("c.tsv"), &p("missing1"), &p("missing2"), &p("o.tsv")).unwrap();
        assert_eq!(fs::read(p("o.tsv")).unwrap(), b"");
    }
}
