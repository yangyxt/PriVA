//! Joins the splicing and 5' UTR copies of a family table, as if the two
//! interpretation steps had run one after another on one table: the
//! combined table's columns, then the splicing step's new columns, then the
//! 5' UTR step's new columns.
//!
//! Both steps read and rewrite the table with the same reader and writer
//! (`pandas_table`), so their copies of the combined table's columns must
//! match cell for cell. The join refuses copies whose row counts differ,
//! whose shared columns differ, or that both added a column of one name.
//! The output is written by the same writer, as the 5' UTR step wrote the
//! table when it ran after the splicing step.
use super::{Result, refuse};
use crate::pandas_table::{Table, read_table, write_table, write_table_to};
use std::fs;
use std::path::Path;

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
    let base = read_table(combined)?;
    let first = read_table(splicing)?;
    let second = read_table(utr)?;
    if first.rows != base.rows || second.rows != base.rows {
        return refuse(format!(
            "Row counts differ: combined {}, splicing {}, 5' UTR {}",
            base.rows, first.rows, second.rows
        ));
    }
    if project(&first, &base.names)? != project(&second, &base.names)? {
        return refuse("The splicing and 5' UTR copies differ in the combined table's columns");
    }
    let mut joined = first.clone();
    for (index, name) in second.names.iter().enumerate() {
        if base.column_index(name).is_some() {
            continue;
        }
        if first.column_index(name).is_some() {
            return refuse(format!(
                "Both the splicing and the 5' UTR step added the column {name}"
            ));
        }
        joined.names.push(name.clone());
        joined.columns.push(second.columns[index].clone());
    }
    write_table(output, &joined)
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

    #[test]
    fn an_empty_table_stays_empty() {
        let f = folder("empty");
        let p = |n: &str| f.0.join(n);
        fs::write(p("c.tsv"), "").unwrap();
        join(&p("c.tsv"), &p("missing1"), &p("missing2"), &p("o.tsv")).unwrap();
        assert_eq!(fs::read(p("o.tsv")).unwrap(), b"");
    }
}
