//! Opt-in DAS domain procedures. The raw archive's tuple GO terms are kept
//! distinct from the runtime functionality interpreter's list-only profile.
use crate::builder_error;
use crate::cache_convert::interpro_mapping::InterproCache;
use crate::cli::{Result, invalid};
use crate::das_rescue::Entry;
use crate::json_text::{self, Json};
use crate::protein_domain_mapping::{DomainNormalizer, Lookup};
use crate::python_set_order::PairSet;
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};

pub(crate) struct LocatedEntry {
    pub(crate) id: String,
    pub(crate) entry: Entry,
}

pub(crate) struct DomainLookup {
    normalizer: DomainNormalizer,
    mapping: InterproCache,
    hash_seed: u32,
}

impl DomainLookup {
    pub(crate) fn new(
        normalizer: DomainNormalizer,
        mapping: InterproCache,
        hash_seed: u32,
    ) -> Self {
        Self {
            normalizer,
            mapping,
            hash_seed,
        }
    }

    /// Source `domain_path.split(':', 1)[-1]`, then normalization into a newly
    /// constructed set of two-text tuples. First entry and GO tuple order are
    /// source behavior; neither a sorted pair set nor the runtime projection
    /// can substitute for this DAS-specific query.
    pub(crate) fn first(&self, domain_path: &str) -> Lookup<Option<LocatedEntry>> {
        let annotation = domain_path
            .split_once(':')
            .map_or(domain_path, |(_, rest)| rest);
        let mut normalized = PairSet::new(self.hash_seed);
        for part in annotation.split('&') {
            if let Some((database, accession)) = part.split_once(':') {
                if let Some(pair) = self
                    .normalizer
                    .normalize_vep_to_interpro(database, accession)?
                {
                    normalized.insert(pair);
                }
            }
        }
        for (database, accession) in normalized.into_order() {
            if let Some(entry) = self
                .mapping
                .get(&database, &accession)
                .and_then(|entries| entries.first())
            {
                return Ok(Some(LocatedEntry {
                    id: self.mapping.text(&entry.id).to_owned(),
                    entry: Entry {
                        kind: self.mapping.text(&entry.kind).to_owned(),
                        has_go_molecular_function: entry
                            .go_terms
                            .iter()
                            .any(|term| self.mapping.text(&term.category) == "molecular_function"),
                    },
                }));
            }
        }
        Ok(None)
    }
}

const HELP: &str = "usage: priva das-domains lookup --domains JSONL --interpro RKYV --normalization-json JSON --python-hash-seed UINT32\n\
       priva das-domains transcripts --statistics PARQUET --exon-history LINES --source-builder-sha256 SHA256 --python-hash-seed UINT32 --output RKYV [--history-builder-sha256 SHA256]\n\
       priva das-domains transcript-provenance --input RKYV\n\
       priva das-domains benign --domains JSONL --vcf VCF --min-benign INTEGER\n\
       priva das-domains benign-chromosome --domains JSONL --vcf INDEXED_VCF --chrom CONTIG --output NEW_JSON\n\
       priva das-domains benign-merge --domains JSONL --min-benign INTEGER --input CHR1_JSON --input CHR2_JSON ...\n\
       priva das-domains coldspot-table --spots TSV\n\
       priva das-domains coldspot --domains JSONL --metadata PARQUET --spots TSV --interpro RKYV --normalization-json JSON --python-hash-seed UINT32 --overlap-threshold FLOAT\n\
Each domain input line is a full domain-path JSON string. Output lines retain input\n\
order and contain the first source-ordered InterPro entry's identity, type and\n\
GO molecular-function status. The explicit hash seed describes the reference\n\
process; no global seed or installed cache is changed. The transcripts procedure\n\
requires retained serialized-set history bound to every typed transcript; it\n\
does not qualify fresh native Alpha builds that lack this history. Benign counts\n\
every qualifying ClinVar annotation occurrence and emits counts/removal flags\n\
in requested domain order. It reads sequentially without an index. Coldspot\n\
combines native metadata, source-ordered InterPro lookup and pandas TSV parsing;\n\
it emits exact fractions/removal flags in requested domain order. Neither filter\n\
procedure reconstructs the enclosing DAS candidate-set history. All are opt-in.";

pub(crate) fn dispatch(args: &[String]) -> Result<()> {
    if args
        .first()
        .is_some_and(|s| matches!(s.as_str(), "benign-chromosome" | "benign-merge"))
    {
        return crate::das_benign_shards::dispatch(args);
    }
    if args.is_empty() || args.iter().any(|arg| arg == "--help" || arg == "-h") {
        println!("{HELP}");
        return Ok(());
    }
    if args[0] == "benign" {
        if args.len() != 7 {
            return Err(invalid(HELP));
        }
        let mut options = HashMap::new();
        for pair in args[1..].chunks_exact(2) {
            if !["--domains", "--vcf", "--min-benign"].contains(&pair[0].as_str())
                || options.insert(pair[0].as_str(), pair[1].as_str()).is_some()
            {
                return Err(invalid(HELP));
            }
        }
        let minimum = crate::das_benign::minimum(options["--min-benign"])?;
        return builder_error::finish(benign_lines(&options, minimum));
    }
    if args[0] == "coldspot-table" {
        if args.len() != 3 || args[1] != "--spots" {
            return Err(invalid(HELP));
        }
        return builder_error::finish(coldspot_rows(std::path::Path::new(&args[2])));
    }
    if args[0] == "transcript-provenance" {
        if args.len() != 3 || args[1] != "--input" {
            return Err(invalid(HELP));
        }
        let cache =
            crate::cache_convert::transcript_domains::TranscriptDomainCache::open(&args[2])?;
        let rows = cache
            .provenance()
            .into_iter()
            .map(|(key, value)| Json::Array(vec![Json::Str(key), Json::Str(value)]))
            .collect();
        println!(
            "{}",
            json_text::dumps(&Json::Array(rows), false).map_err(invalid)?
        );
        return Ok(());
    }
    if args[0] == "coldspot" {
        if args.len() != 15 {
            return Err(invalid(HELP));
        }
        let mut options = HashMap::new();
        for pair in args[1..].chunks_exact(2) {
            if ![
                "--domains",
                "--metadata",
                "--spots",
                "--interpro",
                "--normalization-json",
                "--python-hash-seed",
                "--overlap-threshold",
            ]
            .contains(&pair[0].as_str())
                || options.insert(pair[0].as_str(), pair[1].as_str()).is_some()
            {
                return Err(invalid(HELP));
            }
        }
        let seed = options["--python-hash-seed"]
            .parse::<u32>()
            .map_err(|_| invalid("--python-hash-seed requires UINT32"))?;
        let threshold = crate::python_text::parse_float(options["--overlap-threshold"])
            .ok_or_else(|| invalid("--overlap-threshold requires a float"))?;
        return builder_error::finish(crate::das_coldspot_filter::lines(&options, seed, threshold));
    }
    if args[0] == "transcripts" {
        if ![11, 13].contains(&args.len()) {
            return Err(invalid(HELP));
        }
        let mut options = HashMap::new();
        for pair in args[1..].chunks_exact(2) {
            if ![
                "--statistics",
                "--exon-history",
                "--source-builder-sha256",
                "--history-builder-sha256",
                "--python-hash-seed",
                "--output",
            ]
            .contains(&pair[0].as_str())
                || options.insert(pair[0].as_str(), pair[1].as_str()).is_some()
            {
                return Err(invalid(HELP));
            }
        }
        if ![
            "--statistics",
            "--exon-history",
            "--source-builder-sha256",
            "--python-hash-seed",
            "--output",
        ]
        .iter()
        .all(|key| options.contains_key(key))
        {
            return Err(invalid(HELP));
        }
        let seed = options["--python-hash-seed"]
            .parse::<u32>()
            .map_err(|_| invalid("--python-hash-seed requires UINT32"))?;
        let source = options["--source-builder-sha256"];
        if source.len() != 64 || !source.bytes().all(|b| b.is_ascii_hexdigit()) {
            return Err(invalid("--source-builder-sha256 requires SHA256"));
        }
        let history_builder = options.get("--history-builder-sha256").copied();
        if history_builder
            .is_some_and(|value| value.len() != 64 || !value.bytes().all(|b| b.is_ascii_hexdigit()))
        {
            return Err(invalid("--history-builder-sha256 requires SHA256"));
        }
        let (transcripts, sets) =
            builder_error::finish(crate::alpha_domain::das_transcripts::build(
                std::path::Path::new(options["--statistics"]),
                std::path::Path::new(options["--exon-history"]),
                source,
                history_builder,
                seed,
                std::path::Path::new(options["--output"]),
            ))?;
        println!("{{\"transcripts\":{transcripts},\"exon_sets\":{sets}}}");
        return Ok(());
    }
    if args[0] != "lookup" || args.len() != 9 {
        return Err(invalid(HELP));
    }
    let mut options = HashMap::new();
    for pair in args[1..].chunks_exact(2) {
        if ![
            "--domains",
            "--interpro",
            "--normalization-json",
            "--python-hash-seed",
        ]
        .contains(&pair[0].as_str())
            || options.insert(pair[0].as_str(), pair[1].as_str()).is_some()
        {
            return Err(invalid(HELP));
        }
    }
    if options.len() != 4 {
        return Err(invalid(HELP));
    }
    let seed = options["--python-hash-seed"]
        .parse::<u32>()
        .map_err(|_| invalid("--python-hash-seed requires UINT32"))?;
    builder_error::finish(lookup_lines(&options, seed))
}

fn coldspot_rows(path: &std::path::Path) -> Result<()> {
    let rows = crate::das_coldspot_table::read(path)?;
    let mut output = std::io::BufWriter::new(std::io::stdout().lock());
    for row in rows {
        let value = Json::Object(vec![
            ("gene".into(), row.gene.json()),
            ("start".into(), Json::Str(row.start)),
            ("end".into(), Json::Str(row.end)),
        ]);
        writeln!(
            output,
            "{}",
            json_text::dumps(&value, false).map_err(invalid)?
        )?;
    }
    output.flush()?;
    Ok(())
}

fn benign_lines(options: &HashMap<&str, &str>, minimum: num_bigint::BigInt) -> Result<()> {
    let input = BufReader::new(std::fs::File::open(options["--domains"])?);
    let mut domains = Vec::new();
    for line in input.lines() {
        let Json::Str(domain) = json_text::parse(&line?, true).map_err(invalid)? else {
            return Err(invalid("each domain-path line must be a JSON string"));
        };
        domains.push(domain);
    }
    let wanted: HashSet<_> = domains.iter().cloned().collect();
    let counts = crate::das_benign::counts(std::path::Path::new(options["--vcf"]), &wanted)?;
    let mut output = std::io::BufWriter::new(std::io::stdout().lock());
    for domain in domains {
        let count = counts.get(&domain).copied().unwrap_or(0);
        let row = Json::Object(vec![
            ("domain".into(), Json::Str(domain)),
            ("n_benign_missense".into(), Json::Int(count as i128)),
            (
                "removed".into(),
                Json::Bool(num_bigint::BigInt::from(count) >= minimum),
            ),
        ]);
        writeln!(
            output,
            "{}",
            json_text::dumps(&row, false).map_err(invalid)?
        )?;
    }
    output.flush()?;
    Ok(())
}

fn lookup_lines(options: &HashMap<&str, &str>, seed: u32) -> Result<()> {
    let normalizer = DomainNormalizer::from_parts(
        &std::fs::read_to_string(options["--normalization-json"])?,
        HashMap::new(),
    )?;
    let lookup = DomainLookup::new(
        normalizer,
        InterproCache::open(options["--interpro"])?,
        seed,
    );
    let input = BufReader::new(std::fs::File::open(options["--domains"])?);
    let mut output = std::io::BufWriter::new(std::io::stdout().lock());
    for line in input.lines() {
        let Json::Str(domain) = json_text::parse(&line?, true).map_err(invalid)? else {
            return Err(invalid("each domain-path line must be a JSON string"));
        };
        let found = lookup.first(&domain).map_err(invalid)?;
        let entry = found.map_or(Json::Null, |found| {
            Json::Object(vec![
                ("id".into(), Json::Str(found.id)),
                ("type".into(), Json::Str(found.entry.kind)),
                (
                    "go_molecular_function".into(),
                    Json::Bool(found.entry.has_go_molecular_function),
                ),
            ])
        });
        let row = Json::Object(vec![
            ("domain".into(), Json::Str(domain)),
            ("entry".into(), entry),
        ]);
        writeln!(
            output,
            "{}",
            json_text::dumps(&row, false).map_err(invalid)?
        )?;
    }
    output.flush()?;
    Ok(())
}
