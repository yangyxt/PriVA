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

**Versioning**: a change of any archived type bumps the schema version; a reader accepts only its
own version and says which it found. A change of rkyv's binary format (a new rkyv version or
feature set) gets a new layout number. The container version changes only with the header.

**Compression**: none, so that the file can be mapped. Repeated text is stored once instead: low
cardinality texts (ClinVar significance and review status, consequences, exon numbers, domain
paths, SpliceAI scores) go into a `texts` pool (`Vec<String>`) and records hold `u32` indexes.
For files bundled in git, see the distribution policy in section 6.

**Provenance**: every root holds `provenance: Vec<{key, value}>` with at least `source` (the
pickle or input it was made from), `source_bytes` and `converter` (`priva cache-convert` and the
crate version).

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
| `int` | `i64` | integers past 64 bits are refused; `bool` is `bool`, never folded into an int |
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

All ten schemas are frozen as version 1, with field names taken from the builder scripts.
File names: the pickle's name
with `.pkl`/`.pkl.gz` replaced by `.v<schema version>.rkyv`. The branch's configuration keeps
the key names and points them at the `.rkyv` files.

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
Changed sources get a new archive path, preserving the previous configuration until the
final publication. Keep those source files and their scientific companion receipts while
the corresponding builders still use them.
