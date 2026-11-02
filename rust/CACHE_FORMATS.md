# Rust-native cache formats of the branch (frozen 2026-10-31)

The branch writes and reads its run-time caches in formats chosen for Rust's IO and memory, and
its `config.yaml` points at them; main keeps its Python pickles. A cache stays shared with main
only when its format suits both languages as it is (VCF, TSV, JSON, gzip JSON). This note fixes
the formats, the container, the rules for turning Python values into them, and how a reader
checks a file. The benchmark behind the choice is in section 5.

## 1. One format per shape

| Shape | Format | Caches |
| --- | --- | --- |
| A dict looked up by key (transcript, gene, variant, domain), or a set asked for membership | **rkyv 0.8 archive in the `PRIVARKV` container, memory-mapped** (section 2) | `clinvar_aa_stat`, `clinvar_splice_stat`, `alphamissense_tranx_domain_map`, `clinvar_patho_exon_af_stat`, `clinvar_gene_stat`, `clinvar_patho_af_stat`, `pm1_regions_pkl`, `all_intolerant_domains`, `interpro_mapping_pickle`, `clingen_map` |
| A table (rows and typed columns, read whole or scanned, worked on with polars) | **Parquet, zstd**, read and written with polars `=0.51.0` | builder intermediates of groups 4 and 5 (for example the AlphaMissense domain statistics that only DAS reads); no run-time pickle is a table today |
| Small text resources main also reads | unchanged: TSV, JSON, gzip JSON, VCF | `alphamissense_intolerant_domains` (TSV), `hpo_assertions`, the mechanism JSONs, `nm_to_enst_tsv`, the VCFs |

Not used: Arrow IPC. Plain IPC (the memory-mappable kind) stores each string in full, so its
files are 2.7 to 8 times larger than the pickles; random access through a memory map of it on
StorNext took 16 to 91 s cold for 11,000 keys. Compressed IPC cannot be mapped and is larger and
slower to open than Parquet.

## 2. The `PRIVARKV` container (container version 1)

A file is a 32-byte header followed by exactly one rkyv archive. Integers are little-endian.

| Bytes | Field |
| --- | --- |
| 0..8 | magic `PRIVARKV` |
| 8..12 | container version, `1` |
| 12..16 | schema id (section 4) |
| 16..20 | schema version |
| 20..24 | archive layout: `1` = rkyv 0.8 (crate pinned `=0.8.18`), little-endian, aligned, 32-bit relative pointers (archives up to 2 GiB) |
| 24..32 | archive length in bytes; the file is exactly header plus archive |

The archive starts at offset 32, so a page-aligned memory map keeps rkyv's 16-byte alignment.

**Reading** (`rust/src/rust_cache.rs`, `MappedCache::open`):

1. Map the file read-only.
2. Check magic, container version, schema id, schema version (the exact pair the reader was
   built for), layout and length. A Python pickle has no `PRIVARKV` magic and is refused with a
   message naming `priva cache-convert`.
3. Validate the whole archive once with `rkyv::access` (bytecheck): every pointer, length and
   UTF-8 string is checked, so a damaged or truncated file is an error, never undefined
   behaviour. Measured cost: 0.09 to 0.5 s warm, 0.35 to 1.7 s cold for 170 to 660 MB, and on
   StorNext it is also the fastest way to read a cold file (it reads the file in order; the
   unchecked map took 1.6 to 12 s cold for the same look-ups, section 5).
4. Check what bytecheck cannot: top-level keys sorted by bytes and unique, one entry per key,
   every pool index in range.
5. Hand out `&Archived...` references; look-ups binary-search the sorted key vector and touch
   only the pages of the entries asked for. A reader copies out only what it keeps.

**Writing** (`rust_cache::write_cache`): the archive is built in memory, written with its header
to a temporary name in the target folder, synced and renamed into place, so a reader never maps a
half-written file. Caches are never modified in place. The same input gives the same bytes.

**Versioning**: a change of any archived type bumps the schema version; a reader accepts only
explicitly supported versions and selects each version's unchanged archived type. Unsupported
versions are refused. A change of rkyv's binary format (a new rkyv version or
feature set) gets a new layout number. The container version changes only with the header.

**Compression**: none, so that the file can be mapped. Repeated text is stored once instead: low
cardinality texts (ClinVar significance and review status, consequences, exon numbers, domain
paths, SpliceAI scores) go into a `texts` pool (`Vec<String>`) and records hold `u32` indexes.
For files bundled in git, see the distribution policy in section 6.

**Provenance**: every root holds `provenance: Vec<{key, value}>`. Converters record `source`
(the pickle or input it was made from), `source_bytes` and `converter` (`priva cache-convert`
and the crate version). Direct builders record the identities of their prepared inputs;
the PM1 builder records a SHA256 for each of its four sources.

## 3. Python values in the archive

Every cache has a typed schema (Rust structs deriving `rkyv::Archive`); there is no generic
"Python value" type. The converter reads the pickle with the branch's own plain-data reader
(`pickle_dict_stream`) and refuses any value the schema does not describe, so a converted file
holds everything the pickle held, in the readers' terms:

| Python | Archived as | Rule |
| --- | --- | --- |
| top-level `dict` | `keys: Vec<String>` sorted by bytes plus a parallel `entries: Vec<V>` | Python's top-level order is not kept (no reader uses it); a key set twice in the pickle is refused |
| inner `dict` with free keys | `Vec<{key, value}>` in Python's insertion order | order kept |
| inner `dict` with fixed keys (a record) | a struct with one field per key | other keys refused; a key that may be absent is `Option<T>`; key order not kept (no reader iterates record keys) |
| `list` | `Vec<T>`, order kept | a `tuple` where the schema says list is refused |
| `tuple` of fixed length | a struct (or fixed fields) | a `list` where the schema says tuple is refused; the schema's doc comment names the Python type of each field, so the value can be rebuilt |
| `set`, `frozenset` | `Vec<T>` sorted by bytes (numbers by value), unique | Python's set order is hash-random and no reader may depend on it |
| `None` | `Option<T>::None` | where both `None` and an absent key occur, an enum keeps them apart |
| `float` | `f64`, bits kept exactly | NaN is stored as NaN; where a field mixes text and NaN (the ClinGen map) the schema says so (`None` = NaN) and refuses other kinds |
| `int` | `i64` | integers past signed 64 bits are refused except for the explicitly versioned ClinGen scalar extension below; `bool` is `bool`, never folded into an int |
| `str` | `String` (UTF-8) | low-cardinality fields go through the `texts` pool |
| numpy arrays, other classes | refused | (the pickle reader refuses them already) |

In Parquet: one row per innermost record, the top-level key as a column, rows sorted by it with
the inner order kept; `None` is null, NaN stays NaN in `Float64` (polars keeps them apart);
lists are `List` columns; sets are sorted `List` columns.

## 4. Schemas

| Id | Cache | Version 1 layout (root fields after `provenance`) | State |
| ---: | --- | --- | --- |
| 1 | `clinvar_aa_stat` | `texts`; `transcripts` (sorted); `entries[i].positions: Vec<{protein_pos, changes: Vec<{hgvsp, clnsig: Vec<u32>, clnrevstat: Vec<u32>}>}>` | **implemented** (`rust/src/cache_convert.rs`) |
| 2 | `clinvar_splice_stat` | `texts`; `splice_ai_names` (the nine `SpliceAI_pred_*` keys in pickle order); `transcripts`; `entries[i]: Vec<{chrom: u32, pos: i64, reference, alternative, exon: u32, intron: u32, hgvsc, consequence: u32, clinvar_sig: u32, clinvar_review: u32, splice_ai: Vec<u32>}>` (`splice_ai` empty for `{}`, else the values of `splice_ai_names` in order; any other key set or inconsistent nonempty dictionary order refused) | **implemented** (`cache_convert/clinvar_splice_stat.rs`) |
| 3 | `alphamissense_tranx_domain_map` | `texts` (domain paths and exon names); `transcripts`; `entries[i]: {exon_domains: Vec<{exon: u32, domains: Vec<u32>}>, pp_domains: Vec<{start: i64, end: i64, domain: u32}>}` (Python: `pp_domains` items are 3-tuples, `exon_domains` values lists; both record keys required; exon dictionary and list order/multiplicity retained) | **implemented** (`rust/src/cache_convert/transcript_domains.rs`) |
| 4 | `clinvar_patho_exon_af_stat` | `transcripts`; `entries[i]: Vec<{exon_key ("<exon>/<count>"), item0..item3: f64, item4: Vec<f64>}>` (the Python 5-tuple; item 4 a tuple of floats) | **implemented** (`cache_convert/clinvar_exon_af.rs`) |
| 5 | `clinvar_gene_stat` | `variants` pool (sorted, unique variant ids); `genes`; `entries[i]: {ensembl_id, gene_symbol, sets, domains: Vec<{domain, sets}>, transcripts: Vec<{transcript, sets}>}`, `sets` = the seven named Python `set` fields as sorted `Vec<u32>`; free-key domain/transcript dict order retained | **implemented** (`cache_convert/clinvar_gene_stat.rs`) |
| 6 | `clinvar_patho_af_stat` | `genes`; `entries[i]: Option<Vec<{key, value: Text/Int/Float}>>` (`None` for Python `None`; only the twelve builder field names are accepted, with missing keys and dict order retained) | **implemented** (`cache_convert/clinvar_gene_af.rs`) |
| 7 | `pm1_regions_pkl` | `das_intolerant_domains: Vec<String>` (sorted Python set); `hcseeker_hotspots` and `rmc_constrained: Vec<{transcript, residues: Vec<i64>}>` (sorted transcripts, sorted Python residue sets); `metadata: Vec<Pm1Metadata>` (closed variants for `assembly: str`, two `rmc_*_threshold: float` and five `n_*: int` fields; field presence, insertion order and float bits retained; undeclared fields/types refused) | **implemented** (`rust/src/cache_convert/pm1_regions.rs`) |
| 8 | `all_intolerant_domains` | `domains`: UTF-8 byte-sorted, unique (a Python set of text; lists and frozen sets refused) | **implemented** (`rust/src/cache_convert/intolerant_domains.rs`) |
| 9 | `interpro_mapping_pickle` | `texts` pool; keys `(database, accession)` sorted; `entries[i]: Vec<{id, kind, short_name, name, go_terms: Vec<{id, category, description}>, abstract_text}>`, every text field a pool index (Python: list of exact 6-tuples, GO terms a list of exact 3-tuples) | **implemented** (`cache_convert/interpro_mapping.rs`) |
| 10 | `clingen_map` | `variants` (sorted); `codes: Vec<Option<String>>` (`None` = the float NaN of curated variants without codes; Python `None` and every non-text/non-NaN value refused) | **implemented** (`cache_convert/clingen_map.rs`) |

All ten version-1 layouts remain frozen, with field names taken from the builder scripts.
File names: the pickle's name
with `.pkl`/`.pkl.gz` replaced by `.v<schema version>.rkyv`. The branch's configuration keeps
the key names and points them at the `.rkyv` files.

### ClinVar splice records (schema 2, version 2)

The splice converter and `priva clinvar-aa-splice` retain version 1 for records
with text fields and either empty SpliceAI dictionaries or the same ordered
nine-key dictionary. The original version-1 converter, archive layout and
reader remain unchanged. Conversion of compatible input retains identical
version-1 bytes; existing bundles are not rewritten.

Version 2 stores `provenance`, the unique `texts` pool, sorted `transcripts`,
and one ordered record list per transcript. Each record keeps the eleven
source fields. `alt`, `exon`, `intron`, `hgvsc` and `consequence` are nullable;
`None` and empty text remain distinct. `splice_ai` is an ordered list of
unique `SpliceAI_pred_*` key/value pairs, with text values. Key subsets,
additional names and differing per-record dictionary order are preserved.
`chrom`, clinical significance/review, exon/intron/consequence, and SpliceAI
keys/values use checked pool indices; reference, alternative and HGVSc remain
inline strings. Position is an `i64`. Fixed record keys export in the same
declared order as version 1; inner dictionary and record-list order remain intact.

The container stays `PRIVARKV` version 1. An older schema-2-v1-only reader
refuses version 2. Managed conversion selects `.v1.rkyv` or `.v2.rkyv` and
checks that version before publication, recovery or reuse. Prior archives
and receipts are retained. The splice consumer treats absent scores as NaN,
missing HGVSc as an unavailable parsed location, and evaluates a missing ALT
only if a matching record reaches its length calculation, where Python raises
`TypeError`. These are representation and source-behavior rules; the installer
does not activate a different VCF producer.

### ClinVar gene AF values (schema 6, version 2)

The gene-AF converter and `priva clinvar-gene-af` retain version 1 whenever
all surviving record values are text, integer or float. That path uses the
unchanged version-1 serializer, including its provenance, field order and float
bits. Existing archives are read with their declared schema version.

Version 2 retains `provenance`, sorted unique `genes`, and
`entries: Vec<Option<Vec<Entry>>>`. Each ordered `Entry` holds a known field
name and a value whose declared enum order is `Text(String)`, `Int(i64)`,
`Float(f64)`, `Bool(bool)`. Boolean values remain distinct from integers;
missing fields, empty records and `None` genes remain distinct. The container
version stays 1. An older reader that supports only schema version 1 refuses
version 2.

Both versions expose the same checked borrowed record view and typed export.
PM2 preserves pandas' Boolean or mixed-object gene-AF column; PM2 and BS1 use
Boolean zero/one only in their gene-AF arithmetic. Other numeric-column rules
are unchanged. Managed installation chooses `.v1.rkyv` or `.v2.rkyv` from the
actual values, checks the declared version before publication or reuse, and
keeps prior archives and receipts. Creating an archive does not activate the
native VCF builder in an installer caller.

### ClinGen scalar codes (schema 10, version 2)

The ClinGen converter and `priva clingen-map build` select version 1 whenever every
surviving code is text or missing. They use the original version-1 serializer, so
the same values and provenance keep identical archive bytes. Bundled version-1
files are unchanged. The ClinGen reader accepts both versions and validates the
selected version's archive; an older version-1-only reader refuses version 2.

Version 2 retains `provenance`, sorted unique `variants`, and one tagged `Code`
per variant in `codes`. The declared enum order is `Text(String)`, `Missing`,
`Integer(i128)`, `Float(u64)`, `Boolean(bool)`. Integers are limited to the range
of pandas int64 and uint64 values (`-2^63` through `2^64-1`). Float payloads are
their exact non-NaN IEEE-754 bits, including signed zero and infinities. Float
NaN keeps the existing `Missing` normalization; Python `None`, containers, and
integers outside that range remain refused. Only the ClinGen pickle visitor
opts into uint64; other cache schemas retain their previous integer boundary.

Export preserves scalar types and float bits. Runtime criteria override still
uses only nonempty text after Python-compatible whitespace stripping. Present
numeric, Boolean, or missing codes keep PriVA's assigned criteria, matching the
Python consumer; they remain distinct from absent keys in the cache. Creating
a version-2 archive does not select it in a configuration or change a caller.

## 5. Format selection

Keyed caches use checked memory-mapped rkyv archives so readers can access requested entries
without rebuilding whole Python object graphs. Shared strings are pooled; dictionaries keep
scientifically relevant inner order and exact values. Checked access rejects malformed archives
before any borrowed reference is exposed. Tabular intermediates use compressed Parquet through
Polars; the keyed-cache layout and table format serve different access patterns.

## 6. Distribution and installation

- Bundle the InterPro archive and both assembly-specific ClinGen archives alongside their
  existing pickles. These are portable snapshots; normal bundled reuse checks only existence
  and nonzero size.
- Keep the compact ClinVar gene-statistics pickle in Git. Convert it during installation
  into the configured `cache_root`, preserving all original sets and values. Its expanded
  archive is not a Git asset; Git LFS and a different archive layout are not required.
- Convert other acquired or derived pickle caches once into managed archives under
  `cache_root`. Keep the existing YAML key names and publish their archive paths only after
  conversion succeeds. Preserve the source pickles while installation-time Python builders
  still use them; changing storage format does not require repeating VEP annotation.
- `clinvar_gene_stat` holds 13 million variant-id strings in 7 sets per gene, domain and
  transcript; the PP2/BP1 reader uses set sizes and intersections. Version 1 keeps the sets
  (as pooled `u32` ids), not only their counts.

`priva cache-install prepare --config CONFIG --build-config PRIVATE --source-root REPO`
creates the builder view, restoring the recorded pickle paths when CONFIG already selects
native archives. Use a stable PRIVATE filename so interrupted builders can resume.
`priva cache-install publish` with the same options converts the completed resources and
atomically replaces CONFIG. It retains pending builder updates, refuses conflicting user
edits, and reuses completed conversions. It does not run the scientific builders itself.

Managed archive names contain the existing YAML key, assembly, VEP release, source digest
prefix and format version. They stay in the configured resource directories under
`cache_root`; each has a `.source.json` receipt identifying its original builder input.
ClinGen managed names use `.v1.rkyv` or `.v2.rkyv` according to their actual schema.
Version-2 receipts add `schema_version: "2"`; existing version-1 receipts remain
unchanged. Both use `format: "PRIVARKV-v1"` for the unchanged container. Publication
and reuse check that the archive's schema matches its managed generation.
Changed sources get a new archive path, preserving the previous configuration until the
final publication. Keep those source files and their scientific companion receipts while
the corresponding builders still use them.

## 7. Build PM1 regions from prepared inputs

The standalone `pm1-regions build` procedure writes schema 7 directly, without a pickle
intermediate. Supply the transcript mapping, existing schema-8 DAS archive, HCSeeker
hotspots and RMC table from the same assembly and scientific profile:

```bash
priva pm1-regions build --assembly hg38 \
  --nm-to-enst data/hgnc/ensembl_GRCh38_mart116_ftp_candidate.tsv \
  --intolerant-domains /path/to/cache/all_intolerant_domains.hg38.v1.rkyv \
  --hcseeker data/HCSeeker/HC_spots.hg38.tsv \
  --rmc data/gnomAD_RMC/gnomad_v2.1.1_rmc_flat.hg38.tsv \
  --output /path/to/cache/pm1_regions.hg38.v1.rkyv
```

This is the prepared-DAS route of `scripts/prepare_pm1_regions.py --skip_das`.
It reads mapping, DAS, HCSeeker and RMC in that order; the RMC arm applies float64
Benjamini-Hochberg adjustment and strict `oe < 0.4` and adjusted `p < 0.05` selection.
Named `--rmc-oe-threshold` and `--rmc-p-bh-threshold` options preserve explicit profile
values. Transcript lookup keys and residue sets follow the schema's canonical order;
all values, empty transcript sets and metadata float bits are retained. Residues must
fit the schema's signed 64-bit integers.

The procedure checks source identities before publishing the completed archive, retaining
an existing output on a source/parse failure. It does not calculate DAS or download data.
The shell installer's existing Python producer remains separate from this opt-in procedure.
