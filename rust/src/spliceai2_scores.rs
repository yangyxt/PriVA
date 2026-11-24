//! SpliceAI2 scores for the variants of a family table that the published
//! SpliceAI2 tables do not hold. Those tables hold every single-base change of
//! the protein-coding genes, but only the insertions, deletions and multi-base
//! substitutions seen in human populations, so a rare one has empty SpliceAI
//! columns. The model itself (scripts/spliceai2_score_variants.py, run by
//! pipeline_steps/spliceai2.rs) can score it.
//!
//! - missing: the variants to score, one line per variant and strand: those of
//!   the rows whose SpliceAI columns are empty, which are not a single-base
//!   change, and whose gene is a scored gene with the variant inside its scored
//!   interval. Variants a store folder already holds are written with their
//!   scores instead.
//! - fill: the table with the scores written into those rows; every other byte
//!   stays as it is.
//!
//! The scored genes and their intervals are the gene context the SpliceAI2
//! files were prepared with (`priva spliceai2-prepare`). Genes match by Ensembl
//! gene ID without its version, so a gene renamed between releases keeps its
//! scores, and a row takes the scores of its own gene's strand. The columns are
//! written as VEP's SpliceAI plugin writes them from the published tables: the
//! four scores and distances, the gene symbol, and PASS when the largest score
//! reaches the plugin's cutoff of 0.5.
use crate::cli::{Result, invalid};
use crate::nextflow_task::table_join::{fields, open};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

const USAGE: &str = "usage: priva spliceai2-scores missing --table TABLE.tsv --gene-context CONTEXT.tsv --new NEW.tsv [--store FOLDER --known KNOWN.tsv]
       priva spliceai2-scores fill --table TABLE.tsv --gene-context CONTEXT.tsv --scores SCORES.tsv --output OUTPUT.tsv
--scores may be given more than once.";

const CONTEXT_HEADER: &str = "gene_id\tchrom\tgene_start\tgene_end\tstrand\tgene_name";
/// The header of a list of variants to score.
pub(crate) const VARIANT_HEADER: &str = "chrom\tpos\tref\talt\tstrand";
/// The header of a table of scores: what scripts/spliceai2_score_variants.py
/// writes and a store folder holds.
pub(crate) const SCORE_HEADER: &str =
    "chrom\tpos\tref\talt\tstrand\tDS_AG\tDS_AL\tDS_DG\tDS_DL\tDP_AG\tDP_AL\tDP_DG\tDP_DL";
const SCORES: [&str; 4] = ["DS_AG", "DS_AL", "DS_DG", "DS_DL"];
const DISTANCES: [&str; 4] = ["DP_AG", "DP_AL", "DP_DG", "DP_DL"];
/// The cutoff VEP's SpliceAI plugin runs with (pipeline_steps/vep.rs): the
/// score from which SpliceAI_cutoff is PASS.
const PLUGIN_CUTOFF: f64 = 0.5;

/// The scored intervals of each gene: chromosome, start, end (1-based, both
/// included) and strand, by Ensembl gene ID without version.
type Genes = HashMap<String, Vec<(String, u64, u64, char)>>;

fn read_genes(path: &Path) -> Result<Genes> {
    let mut lines = BufReader::new(File::open(path)?).lines();
    if lines.next().transpose()?.as_deref() != Some(CONTEXT_HEADER) {
        return Err(invalid(format!(
            "{} is not a SpliceAI2 gene context ({CONTEXT_HEADER})",
            path.display()
        )));
    }
    let mut genes: Genes = HashMap::new();
    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        let [id, chrom, start, end, strand @ ("+" | "-"), _] = parts.as_slice() else {
            return Err(invalid(format!(
                "{}: a line is not a gene of the context: {line:?}",
                path.display()
            )));
        };
        let id = id.split('.').next().unwrap_or(id);
        let interval = (
            (*chrom).to_owned(),
            start.parse::<u64>()?,
            end.parse::<u64>()?,
            if *strand == "+" { '+' } else { '-' },
        );
        let intervals = genes.entry(id.to_owned()).or_default();
        if !intervals.contains(&interval) {
            intervals.push(interval);
        }
    }
    Ok(genes)
}

/// The columns of a family table the step reads and writes.
struct Columns {
    chrom: usize,
    pos: usize,
    reference: usize,
    alternate: usize,
    gene: usize,
    symbol: usize,
    scores: [usize; 4],
    distances: [usize; 4],
    scored_symbol: usize,
    cutoff: Option<usize>,
    width: usize,
}

impl Columns {
    fn of(header: &[u8], path: &Path) -> Result<Self> {
        let names = fields(header);
        let find = |name: &str| names.iter().position(|field| *field == name.as_bytes());
        let need = |name: &str| {
            find(name).ok_or_else(|| invalid(format!("{} has no column {name}", path.display())))
        };
        let four = |kinds: [&str; 4]| -> Result<[usize; 4]> {
            let mut out = [0; 4];
            for (slot, kind) in out.iter_mut().zip(kinds) {
                *slot = need(&format!("SpliceAI_pred_{kind}"))?;
            }
            Ok(out)
        };
        Ok(Self {
            chrom: need("chrom")?,
            pos: need("pos")?,
            reference: need("ref")?,
            alternate: need("alt")?,
            gene: need("Gene")?,
            symbol: need("SYMBOL")?,
            scores: four(SCORES)?,
            distances: four(DISTANCES)?,
            scored_symbol: need("SpliceAI_pred_SYMBOL")?,
            cutoff: find("SpliceAI_cutoff"),
            width: names.len(),
        })
    }

    /// The variant and strand of a row to score, as a line of the variant
    /// list; `None` for a row the model's scores do not belong in.
    fn wanted(&self, row: &[&[u8]], genes: &Genes) -> Option<String> {
        let bases =
            |allele: &[u8]| !allele.is_empty() && allele.iter().all(|base| b"ACGT".contains(base));
        let (reference, alternate) = (row[self.reference], row[self.alternate]);
        if (reference.len() == 1 && alternate.len() == 1)
            || !bases(reference)
            || !bases(alternate)
            || self.scores.iter().any(|column| !row[*column].is_empty())
        {
            return None;
        }
        let text = |column: usize| std::str::from_utf8(row[column]).ok();
        let (chrom, pos) = (text(self.chrom)?, text(self.pos)?.parse::<u64>().ok()?);
        let (_, _, _, strand) = genes
            .get(text(self.gene)?)?
            .iter()
            .find(|(c, start, end, _)| c == chrom && *start <= pos && pos <= *end)?;
        Some(format!(
            "{chrom}\t{pos}\t{}\t{}\t{strand}",
            text(self.reference)?,
            text(self.alternate)?
        ))
    }
}

/// The scores of a table of scores, by variant and strand; only those in
/// `keys` when it is given. The first file that holds a variant wins.
fn read_scores(
    path: &Path,
    keys: Option<&HashSet<String>>,
    into: &mut HashMap<String, [String; 8]>,
) -> Result<()> {
    let mut lines = BufReader::new(File::open(path)?).lines();
    if lines.next().transpose()?.as_deref() != Some(SCORE_HEADER) {
        return Err(invalid(format!(
            "{} is not a table of SpliceAI2 scores ({SCORE_HEADER})",
            path.display()
        )));
    }
    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() != 13 {
            return Err(invalid(format!(
                "{}: a line has {} fields instead of 13",
                path.display(),
                parts.len()
            )));
        }
        let key = parts[..5].join("\t");
        if keys.is_some_and(|keys| !keys.contains(&key)) || into.contains_key(&key) {
            continue;
        }
        let mut values: [String; 8] = Default::default();
        for (slot, value) in values.iter_mut().zip(&parts[5..]) {
            *slot = (*value).to_owned();
        }
        if values[..4]
            .iter()
            .any(|value| value.parse::<f64>().is_err())
            || values[4..]
                .iter()
                .any(|value| value.parse::<i64>().is_err())
        {
            return Err(invalid(format!(
                "{}: the scores of {key:?} are not numbers",
                path.display()
            )));
        }
        into.insert(key, values);
    }
    Ok(())
}

/// The tables of scores in a store folder, by name.
fn store_files(store: &Path) -> Vec<PathBuf> {
    let Ok(entries) = fs::read_dir(store) else {
        return Vec::new();
    };
    let mut files: Vec<PathBuf> = entries
        .filter_map(|entry| entry.ok().map(|entry| entry.path()))
        .filter(|path| {
            let name = path
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or_default();
            name.starts_with("scores-") && name.ends_with(".tsv")
        })
        .collect();
    files.sort();
    files
}

/// Writes a whole file under another name first, then renames it.
fn write_file(path: &Path, text: &str) -> Result<()> {
    let staged = path.with_extension(format!("{}.partial", std::process::id()));
    fs::write(&staged, text)?;
    fs::rename(staged, path)?;
    Ok(())
}

/// `missing`: how many variants the table needs scores for and how many of
/// them the store held.
pub(crate) fn missing(
    table: &Path,
    context: &Path,
    new: &Path,
    store: Option<(&Path, &Path)>,
) -> Result<(usize, usize)> {
    let mut wanted = BTreeSet::new();
    if fs::metadata(table)?.len() > 0 {
        let genes = read_genes(context)?;
        let mut records = open(table)?;
        records.next()?;
        let columns = Columns::of(&records.record, table)?;
        while records.next()? {
            let row = fields(&records.record);
            if row.len() != columns.width {
                return Err(invalid(format!(
                    "{}: a row has {} fields, the header {}",
                    table.display(),
                    row.len(),
                    columns.width
                )));
            }
            wanted.extend(columns.wanted(&row, &genes));
        }
    }
    let mut known = HashMap::new();
    if let Some((store, known_path)) = store {
        let keys: HashSet<String> = wanted.iter().cloned().collect();
        if !keys.is_empty() {
            for file in store_files(store) {
                read_scores(&file, Some(&keys), &mut known)?;
            }
        }
        let mut text = format!("{SCORE_HEADER}\n");
        for key in wanted.iter().filter(|key| known.contains_key(*key)) {
            text.push_str(&format!("{key}\t{}\n", known[key].join("\t")));
        }
        write_file(known_path, &text)?;
    }
    let mut text = format!("{VARIANT_HEADER}\n");
    for key in wanted.iter().filter(|key| !known.contains_key(*key)) {
        text.push_str(key);
        text.push('\n');
    }
    write_file(new, &text)?;
    Ok((wanted.len(), known.len()))
}

/// `fill`: the number of rows whose scores were written.
pub(crate) fn fill(
    table: &Path,
    context: &Path,
    scores: &[PathBuf],
    output: &Path,
) -> Result<usize> {
    if fs::metadata(table)?.len() == 0 {
        fs::copy(table, output)?;
        return Ok(0);
    }
    let genes = read_genes(context)?;
    let mut known = HashMap::new();
    for path in scores {
        read_scores(path, None, &mut known)?;
    }
    let mut records = open(table)?;
    records.next()?;
    let columns = Columns::of(&records.record, table)?;
    let staged = output.with_extension(format!("{}.partial", std::process::id()));
    let mut out = BufWriter::with_capacity(1 << 22, File::create(&staged)?);
    out.write_all(&records.record)?;
    out.write_all(b"\n")?;
    let mut filled = 0;
    while records.next()? {
        let row = fields(&records.record);
        let values = (row.len() == columns.width)
            .then(|| columns.wanted(&row, &genes))
            .flatten()
            .and_then(|key| known.get(&key));
        let Some(values) = values else {
            out.write_all(&records.record)?;
            out.write_all(b"\n")?;
            continue;
        };
        let mut row: Vec<&[u8]> = row;
        let symbol = row[columns.symbol];
        for (column, value) in columns.scores.iter().chain(&columns.distances).zip(values) {
            row[*column] = value.as_bytes();
        }
        row[columns.scored_symbol] = symbol;
        if let Some(column) = columns.cutoff {
            let best = values[..4]
                .iter()
                .filter_map(|value| value.parse::<f64>().ok())
                .fold(f64::NEG_INFINITY, f64::max);
            row[column] = if best >= PLUGIN_CUTOFF {
                b"PASS"
            } else {
                b"FAIL"
            };
        }
        for (index, field) in row.iter().enumerate() {
            if index > 0 {
                out.write_all(b"\t")?;
            }
            out.write_all(field)?;
        }
        out.write_all(b"\n")?;
        filled += 1;
    }
    out.flush()?;
    drop(out);
    fs::rename(staged, output)?;
    Ok(filled)
}

/// `priva spliceai2-scores`.
pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    let Some(command) = args.first() else {
        return Err(invalid(USAGE));
    };
    if command == "--help" || command == "-h" {
        println!("{USAGE}");
        return Ok(());
    }
    let mut options: HashMap<&str, Vec<&str>> = HashMap::new();
    let mut rest = args[1..].iter();
    while let Some(option) = rest.next() {
        let name = option
            .strip_prefix("--")
            .ok_or_else(|| invalid(format!("unexpected argument {option}\n{USAGE}")))?;
        let value = rest
            .next()
            .ok_or_else(|| invalid(format!("{option} needs a value\n{USAGE}")))?;
        options.entry(name).or_default().push(value);
    }
    let get = |name: &str| options.get(name).and_then(|values| values.last()).copied();
    let need = |name: &str| {
        get(name)
            .map(Path::new)
            .ok_or_else(|| invalid(format!("--{name} is required\n{USAGE}")))
    };
    match command.as_str() {
        "missing" => {
            let store = match (get("store"), get("known")) {
                (Some(store), Some(known)) => Some((Path::new(store), Path::new(known))),
                (None, None) => None,
                _ => return Err(invalid(format!("--store and --known go together\n{USAGE}"))),
            };
            let (wanted, known) =
                missing(need("table")?, need("gene-context")?, need("new")?, store)?;
            println!("{wanted}\t{known}");
            Ok(())
        }
        "fill" => {
            let scores: Vec<PathBuf> = options
                .get("scores")
                .map(|paths| paths.iter().map(PathBuf::from).collect())
                .unwrap_or_default();
            if scores.is_empty() {
                return Err(invalid(format!("--scores is required\n{USAGE}")));
            }
            let filled = fill(
                need("table")?,
                need("gene-context")?,
                &scores,
                need("output")?,
            )?;
            println!("{filled}");
            Ok(())
        }
        _ => Err(invalid(format!("unknown command {command}\n{USAGE}"))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEADER: &str = "chrom\tpos\tref\talt\tNote\tSYMBOL\tGene\tSpliceAI_cutoff\tSpliceAI_pred_DP_AG\tSpliceAI_pred_DP_AL\tSpliceAI_pred_DP_DG\tSpliceAI_pred_DP_DL\tSpliceAI_pred_DS_AG\tSpliceAI_pred_DS_AL\tSpliceAI_pred_DS_DG\tSpliceAI_pred_DS_DL\tSpliceAI_pred_SYMBOL\tLast";

    /// A row with the given variant, symbol, gene and scores.
    fn row(variant: &str, symbol: &str, gene: &str, scores: &str) -> String {
        let note = "\"a \"\"quoted\"\"\ttab\"";
        let empty = "\t\t\t";
        let (scores, symbol_scored) = match scores {
            "" => (empty.to_owned(), ""),
            given => (given.replace(' ', "\t"), symbol),
        };
        format!(
            "{}\t{note}\t{symbol}\t{gene}\t\t{empty}\t{scores}\t{symbol_scored}\tend",
            variant.replace(':', "\t")
        )
    }

    struct Folder(PathBuf);

    impl Folder {
        fn new(name: &str) -> Self {
            let path =
                std::env::temp_dir().join(format!("priva_spliceai2_{name}_{}", std::process::id()));
            let _ = fs::remove_dir_all(&path);
            fs::create_dir_all(&path).unwrap();
            Folder(path)
        }

        fn file(&self, name: &str, text: &str) -> PathBuf {
            let path = self.0.join(name);
            fs::write(&path, text).unwrap();
            path
        }

        fn read(&self, name: &str) -> String {
            fs::read_to_string(self.0.join(name)).unwrap()
        }
    }

    impl Drop for Folder {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    fn context(folder: &Folder) -> PathBuf {
        folder.file(
            "context.tsv",
            &format!(
                "{CONTEXT_HEADER}\nENSG01.7\tchr1\t100\t200\t+\tONE\nENSG01.7\tchr1\t100\t200\t+\tOLD1\nENSG02.1\tchr1\t150\t300\t-\tTWO\n"
            ),
        )
    }

    fn table(folder: &Folder) -> PathBuf {
        let rows = [
            // An unscored deletion: once per transcript and once per strand's gene.
            row("chr1:160:CG:C", "ONE", "ENSG01", ""),
            row("chr1:160:CG:C", "ONE", "ENSG01", ""),
            row("chr1:160:CG:C", "TWO", "ENSG02", ""),
            // A scored one, a single-base change, a gene without scores, a variant
            // outside the interval, another chromosome and a symbolic allele.
            row("chr1:161:C:CA", "ONE", "ENSG01", "0.10 0.00 0.00 0.00"),
            row("chr1:162:C:T", "ONE", "ENSG01", ""),
            row("chr1:163:CG:C", "LNC", "ENSG09", ""),
            row("chr1:250:CG:C", "ONE", "ENSG01", ""),
            row("chr2:160:CG:C", "ONE", "ENSG01", ""),
            row("chr1:164:C:<DEL>", "ONE", "ENSG01", ""),
            // A multi-base substitution the model also scores.
            row("chr1:170:AC:TT", "ONE", "ENSG01", ""),
        ];
        folder.file("table.tsv", &format!("{HEADER}\n{}\n", rows.join("\n")))
    }

    #[test]
    fn missing_lists_each_unscored_variant_once_per_strand() {
        let folder = Folder::new("missing");
        let (context, table) = (context(&folder), table(&folder));
        let new = folder.0.join("new.tsv");
        assert_eq!(missing(&table, &context, &new, None).unwrap(), (3, 0));
        assert_eq!(
            folder.read("new.tsv"),
            format!(
                "{VARIANT_HEADER}\nchr1\t160\tCG\tC\t+\nchr1\t160\tCG\tC\t-\nchr1\t170\tAC\tTT\t+\n"
            )
        );
        // A store gives the variants it holds; other files of the folder are not read.
        let store = folder.0.join("store");
        fs::create_dir_all(&store).unwrap();
        fs::write(
            store.join("scores-a.tsv"),
            format!("{SCORE_HEADER}\nchr1\t160\tCG\tC\t-\t0.00\t0.62\t0.01\t0.00\t5\t-3\t12\t40\nchr9\t1\tA\tAT\t+\t0.00\t0.00\t0.00\t0.00\t1\t1\t1\t1\n"),
        )
        .unwrap();
        fs::write(store.join("notes.txt"), "not a table of scores").unwrap();
        let known = folder.0.join("known.tsv");
        assert_eq!(
            missing(&table, &context, &new, Some((&store, &known))).unwrap(),
            (3, 1)
        );
        assert_eq!(
            folder.read("known.tsv"),
            format!("{SCORE_HEADER}\nchr1\t160\tCG\tC\t-\t0.00\t0.62\t0.01\t0.00\t5\t-3\t12\t40\n")
        );
        assert_eq!(
            folder.read("new.tsv"),
            format!("{VARIANT_HEADER}\nchr1\t160\tCG\tC\t+\nchr1\t170\tAC\tTT\t+\n")
        );
        // An empty table (no variants left after filtration) needs nothing.
        let empty = folder.file("empty.tsv", "");
        assert_eq!(missing(&empty, &context, &new, None).unwrap(), (0, 0));
        assert_eq!(folder.read("new.tsv"), format!("{VARIANT_HEADER}\n"));
    }

    #[test]
    fn fill_writes_only_the_rows_of_the_scored_variants() {
        let folder = Folder::new("fill");
        let (context, table) = (context(&folder), table(&folder));
        let scores = folder.file(
            "scores.tsv",
            &format!(
                "{SCORE_HEADER}\nchr1\t160\tCG\tC\t+\t0.00\t0.49\t0.01\t0.00\t5\t-3\t12\t40\nchr1\t160\tCG\tC\t-\t0.50\t0.00\t0.00\t0.00\t7\t8\t9\t-10\n"
            ),
        );
        let output = folder.0.join("filled.tsv");
        assert_eq!(fill(&table, &context, &[scores], &output).unwrap(), 3);
        let (before, after) = (folder.read("table.tsv"), folder.read("filled.tsv"));
        let (before, after): (Vec<&str>, Vec<&str>) =
            (before.lines().collect(), after.lines().collect());
        assert_eq!(before.len(), after.len());
        let note = "\"a \"\"quoted\"\"\ttab\"";
        for line in [1, 2] {
            assert_eq!(
                after[line],
                format!(
                    "chr1\t160\tCG\tC\t{note}\tONE\tENSG01\tFAIL\t5\t-3\t12\t40\t0.00\t0.49\t0.01\t0.00\tONE\tend"
                )
            );
        }
        assert_eq!(
            after[3],
            format!(
                "chr1\t160\tCG\tC\t{note}\tTWO\tENSG02\tPASS\t7\t8\t9\t-10\t0.50\t0.00\t0.00\t0.00\tTWO\tend"
            )
        );
        for line in (0..before.len()).filter(|line| !(1..=3).contains(line)) {
            assert_eq!(before[line], after[line], "line {line}");
        }
        // A table that is already filled has nothing left for those variants.
        let new = folder.0.join("new.tsv");
        assert_eq!(missing(&output, &context, &new, None).unwrap(), (1, 0));
    }

    #[test]
    fn other_files_are_refused() {
        let folder = Folder::new("refused");
        let (context, table) = (context(&folder), table(&folder));
        let new = folder.0.join("new.tsv");
        let broken = folder.file("broken.tsv", "chrom\tpos\nchr1\t160\n");
        let error = |result: Result<(usize, usize)>| result.unwrap_err().to_string();
        assert!(error(missing(&broken, &context, &new, None)).contains("has no column ref"));
        assert!(
            error(missing(&table, &broken, &new, None)).contains("not a SpliceAI2 gene context")
        );
        let scores = folder.file(
            "scores.tsv",
            &format!("{SCORE_HEADER}\nchr1\t160\tCG\tC\t+\t\t0.49\t0.01\t0.00\t5\t-3\t12\t40\n"),
        );
        let output = folder.0.join("filled.tsv");
        let error = fill(&table, &context, &[scores], &output)
            .unwrap_err()
            .to_string();
        assert!(error.contains("are not numbers"), "{error}");
    }
}
