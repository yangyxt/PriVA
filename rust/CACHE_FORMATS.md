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
For files bundled in git, see the open question in section 6.

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
| 2 | `clinvar_splice_stat` | `texts`; `splice_ai_names` (the nine `SpliceAI_pred_*` keys in pickle order); `transcripts`; `entries[i]: Vec<{chrom: u32, pos: i64, reference, alternative, exon: u32, intron: u32, hgvsc, consequence: u32, clinvar_sig: u32, clinvar_review: u32, splice_ai: Vec<u32>}>` (`splice_ai` empty for `{}`, else the values of `splice_ai_names` in order; any other key set refused) | benchmarked, converter pending |
| 3 | `alphamissense_tranx_domain_map` | `texts` (domain paths and exon names); `transcripts`; `entries[i]: {exon_domains: Vec<{exon: u32, domains: Vec<u32>}>, pp_domains: Vec<{start: i64, end: i64, domain: u32}>}` (Python: `pp_domains` items are 3-tuples, `exon_domains` values lists) | benchmarked, converter pending |
| 4 | `clinvar_patho_exon_af_stat` | `transcripts`; `entries[i]: Vec<{exon_key ("<exon>/<count>"), item0..item3: f64, item4: Vec<f64>}>` (the Python 5-tuple; item 4 a tuple of floats) | pending |
| 5 | `clinvar_gene_stat` | `variants` pool (sorted, unique variant ids); `genes`; `entries[i]: {ensembl_id, gene_symbol, sets, domains: Vec<{domain, sets}>, transcripts: Vec<{transcript, sets}>}`, `sets` = the seven named Python `set` fields as sorted `Vec<u32>`; free-key domain/transcript dict order retained | **implemented** (`cache_convert/clinvar_gene_stat.rs`) |
| 6 | `clinvar_patho_af_stat` | `genes`; `entries[i]: Option<Vec<{key, value: Text/Int/Float}>>` (`None` for Python `None`) | pending |
| 7 | `pm1_regions_pkl` | `das_intolerant_domains` (sorted); `hcseeker_hotspots` and `rmc_constrained` as sorted transcripts with sorted `Vec<i64>` residues; `metadata: Vec<{key, value}>` | pending |
| 8 | `all_intolerant_domains` | `domains`: sorted, unique (a Python set) | pending |
| 9 | `interpro_mapping_pickle` | `texts` pool; keys `(database, accession)` sorted; `entries[i]: Vec<{id, kind, short_name, name, go_terms: Vec<{id, category, description}>, abstract_text}>`, every text field a pool index (Python: list of exact 6-tuples, GO terms a list of exact 3-tuples) | **implemented** (`cache_convert/interpro_mapping.rs`) |
| 10 | `clingen_map` | `variants` (sorted); `codes: Vec<Option<String>>` (`None` = the float NaN of curated variants without codes; Python `None` and every non-text/non-NaN value refused) | **implemented** (`cache_convert/clingen_map.rs`) |

Schemas 2 to 10 are fixed in outline here; each is frozen as version 1 in the commit that adds
its converter, with field names taken from the builder script. File names: the pickle's name
with `.pkl`/`.pkl.gz` replaced by `.v<schema version>.rkyv`. The branch's configuration keeps
the key names and points them at the `.rkyv` files.

## 5. Measurements behind the choice (2026-10-31, HGSC rackc81o, release builds)

Cluster folder: `PriVA_upgrade/validation/rust_conversion_20261029/cache_formats`
(`README_section.md` there has every number). Benchmark program: `source/rust/cache_format_bench`
in that folder (a separate crate that includes the branch's pickle readers; not committed). Each read is its own process: open the file, look up
10,000 random existing keys and 1,000 missing ones, copy the values out; peak memory is VmHWM
(it counts mapped file pages too). Every format gave the pickle reader's values for all 11,000
keys (equal digests).

Total seconds (cold = first read of the file on the node / warm, mean of two runs) and peak MB:

| Cache (pickle size) | pickle stream, wanted keys | pickle stream, all keys (today's PS1 reader) | Parquet zstd | IPC mapped | IPC zstd | rkyv checked | rkyv unchecked |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `clinvar_aa_stat` (217 MB) | 3.0 / 3.4 s, 251 MB | 12.1 s, 662 MB | 1.16 / 0.89 s, 692 MB | 16.0 / 0.51 s, 391 MB | 1.5 / 1.5 s, 688 MB | **0.79 / 0.48 s, 345 MB** | 1.6 / 0.27 s, 165 MB |
| `alphamissense_tranx_domain_map` (131 MB) | 2.4 / 3.2 s, 127 MB | 8.0 s, 135 MB (whole-pickle reader: 8.2 s, 2,444 MB) | 2.6 / 1.2 s, 989 MB | 30.2 / 0.47 s, 445 MB | 2.2 / 2.2 s, 1,063 MB | **0.91 / 0.25 s, 219 MB** | 2.0 / 0.21 s, 115 MB |
| `clinvar_splice_stat` (1,156 MB) | 14.0 / 15.9 s, 649 MB | 42.7 s, 2,097 MB | 4.1 / 2.4 s, 2,133 MB | 91.3 / 1.6 s, 1,226 MB | 5.2 / 5.0 s, 2,572 MB | **2.2 / 1.5 s, 1,109 MB** | 12.0 / 1.0 s, 563 MB |

File sizes (MB) and conversion time from the decoded values (s):

| Cache | pickle | Parquet zstd | IPC plain | IPC zstd | rkyv |
| --- | ---: | ---: | ---: | ---: | ---: |
| `clinvar_aa_stat` | 217 | 28 (4.2 s) | 595 (0.5 s) | 46 (2.1 s) | 245 (4.3 s) |
| `alphamissense_tranx_domain_map` | 131 | 51 (5.1 s) | 1,054 (0.9 s) | 88 (3.3 s) | 170 (4.0 s) |
| `clinvar_splice_stat` | 1,156 | 59 (5.4 s) | 2,194 (1.7 s) | 105 (4.5 s) | 660 (13.1 s) |

Why rkyv for keyed caches: fastest cold and warm on all three, 3 to 30 times faster than today's
pickle streaming, memory below Parquet (the mapped pages are page cache and can be dropped by the
kernel; the anonymous memory is only what the reader copies out), and no decoding step at all for
keys not asked for. Parquet is 4 to 20 times smaller on disk but must decode whole row groups
into memory (2.1 GB for the splice dictionary). Why checked access: it makes a damaged file an
error and, on StorNext, its sequential pass is faster cold than the random faults of an
unchecked map.

## 6. Open questions for the owner

- Bundled caches in git (`clinvar_gene_stat`, `interpro_mapping_pickle`, `clingen_map`): either
  the branch installer converts the bundled pickle into `cache_root` at installation (git
  unchanged, shared with main), or the branch bundles the `.rkyv` files (larger than the gzip
  pickles; a zstd-compressed copy would need unpacking before mapping). Recommended: convert at
  installation.
- `clinvar_gene_stat` holds 13 million variant-id strings in 7 sets per gene, domain and
  transcript; the PP2/BP1 reader uses set sizes and intersections. Version 1 keeps the sets
  (as pooled `u32` ids), not only their counts.
