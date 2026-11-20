//! A tab-separated table read as Python's `csv.DictReader(handle,
//! delimiter="\t")` reads it, for the installer's table checks
//! (`mechanism-evidence validate`, `mondo-scope validate`): the first record
//! names the fields, later blank lines are skipped, and rows may be shorter
//! or longer than the header. Text that is not UTF-8 is an error, as Python's
//! decoder makes it one.
use std::io::Read;

pub(crate) struct DictRows<R: Read> {
    records: csv::StringRecordsIntoIter<R>,
    pub(crate) header: Vec<String>,
}

pub(crate) fn open<R: Read>(source: R) -> Result<DictRows<R>, String> {
    let mut records = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .flexible(true)
        .from_reader(source)
        .into_records();
    let header = match records.next() {
        None => Vec::new(),
        Some(record) => record
            .map_err(|error| error.to_string())?
            .iter()
            .map(str::to_owned)
            .collect(),
    };
    Ok(DictRows { records, header })
}

impl<R: Read> DictRows<R> {
    /// The position of a field: a later column of the same name wins, as in
    /// `dict(zip(fieldnames, row))`.
    pub(crate) fn column(&self, name: &str) -> Option<usize> {
        self.header.iter().rposition(|field| field == name)
    }
}

impl<R: Read> Iterator for DictRows<R> {
    type Item = Result<Vec<String>, String>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let record = match self.records.next()? {
                Ok(record) => record,
                Err(error) => return Some(Err(error.to_string())),
            };
            if record.len() == 1 && record[0].is_empty() {
                continue;
            }
            return Some(Ok(record.iter().map(str::to_owned).collect()));
        }
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn blank_lines_short_rows_and_quotes() {
        let text = "a\tb\n\n1\t2\n3\n\"x\ty\"\t4\t5\n\n";
        let mut rows = super::open(text.as_bytes()).unwrap();
        assert_eq!(rows.header, ["a", "b"]);
        let rows: Vec<Vec<String>> = rows.by_ref().map(Result::unwrap).collect();
        assert_eq!(
            rows,
            vec![vec!["1", "2"], vec!["3"], vec!["x\ty", "4", "5"]]
        );
    }
}
