# GoFCards cache: schema reference

Deployed artifact: `/paedyl01/disk1/yangyxt/PriVA/data/gofcards/gofcards_exact_gof.json.gz`

There is exactly one deployed GoFCards cache, and it already carries its ClinVar
conditions. It is built in two steps, the first of which writes a build
intermediate rather than a second deliverable:

| Step | Script | Writes |
| --- | --- | --- |
| `gofcards_exact_gof_cache_install` | `/paedyl01/disk1/yangyxt/PriVA/scripts/build_gofcards_exact_gof_cache.py` | `<gofcards_workdir>/gofcards_exact_gof.normalized.json.gz` (intermediate) |
| `gofcards_clinvar_injection_install` | `/paedyl01/disk1/yangyxt/PriVA/scripts/inject_clinvar_into_gofcards.py` | `data/gofcards/gofcards_exact_gof.json.gz` (deployed) |

Both are functions in `/paedyl01/disk1/yangyxt/PriVA/scripts/install_utils.sh`.
The split exists so a ClinVar refresh does not force the VEP pipeline to rerun,
not because two caches are published.

`mechanism_resource_install` runs both steps, and runs one more between them.
Injection reads the weekly ClinVar VCV XML, but the component that downloads
that XML is the non-LOF builder, which in turn cannot start until injection has
written the deployed cache it reads. The installer breaks that standstill by
calling the non-LOF builder twice: once as
`gene_nonlof_mechanism_cache_install <config> sources-only`, which downloads and
stops, and once normally after injection. The resulting order is:

1. `gofcards_exact_gof_cache_install` — normalize, write the intermediate
2. `gene_nonlof_mechanism_cache_install <config> sources-only` — place the ClinVar VCV XML
3. `gofcards_clinvar_injection_install` — write the deployed cache
4. `gene_nonlof_mechanism_cache_install <config>` — build the non-LOF cache
5. `hpo_condition_mechanism_cache_install` — build the condition cache

## When each step rebuilds

Neither step rebuilds on a timer alone. Both compare the cache they own against
the inputs that actually determine its contents.

**Step 1, normalization.** Rebuilds when the normalization script, the mechanism
review table, the HGNC table, the liftover chain, or either FASTA is newer than
the intermediate; when the intermediate records a different VEP annotation
cache, VEP cache version, or chain than the one now configured; or when it is
older than `gofcards_refresh_days` (180). The refresh interval is an upper bound
for picking up a new GoFCards release, not the whole test — polarity review is
applied only while building, so an allele newly marked LOF or mixed-mechanism in
the review table would otherwise keep its GOF eligibility until the interval
lapsed. `PRIVA_FORCE_GOFCARDS_CACHE=1` or `PRIVA_FORCE_ALL_CACHES=1` forces it.

The rebuild is affordable because the slow part is cached: the roughly 2,000
sequential per-allele calls to the GoFCards summary endpoint are appended to
`<gofcards_workdir>/gofcards_summary_cache.jsonl`, and a rerun fetches only the
alleles missing from it.

**Step 3, ClinVar injection.** Rereads the entire 5.8 GB VCV XML, so the only
question the gate asks is whether rereading it could change anything.

It answers that from what the cache records about its own making, not from file
timestamps. The injector writes three facts into `metadata.clinvar`: its own
`injector_sha256`, the `source_cache_sha256` of the normalized intermediate it
read, and the `min_review_stars` it applied. Step 3 reruns when any of those
differs from what is in use now, or when the cache no longer validates. A
content hash is used rather than a modification time because a timestamp cannot
distinguish a changed file from a copied or checked-out one, and getting that
wrong costs hours. When a hash does differ the rebuild happens immediately: a
correction to the injector should reach the cache at once, not after a waiting
period.

The XML is the one exception, because it is the one input that changes
constantly and almost never matters. A new VCV release lands every week and
rarely alters the conditions attached here, so the cache is allowed to trail it
and is refreshed once the gap reaches `gofcards_clinvar_reinjection_lag_days`,
which defaults to 90 days.

`PRIVA_FORCE_GOFCARDS_CLINVAR_INJECTION=1` or `PRIVA_FORCE_ALL_CACHES=1` forces it.

The deployed cache is published by writing a temporary file beside it and moving
that into place, so an injection that dies partway through leaves the previous
complete cache untouched rather than a truncated one where every consumer looks.

## Pipeline position

```
  normalization                       ClinVar VCV injection
  build_gofcards_exact_gof_cache      inject_clinvar_into_gofcards
        |                                       |
        v                                       v
  <workdir>/gofcards_exact_gof   ----->  gofcards_exact_gof.json.gz
           .normalized.json.gz                  |
        (build intermediate)          the ONE deployed cache; every
                                      consumer reads this file
                                             |
                        +--------------------+--------------------+
                        |                                         |
              DIRECT VARIANT MATCH                     CONDITION-BASED MATCH
              gene_mechanism_hub.py                    build_hpo_condition_
              reads record + assemblies                mechanism_cache.py
              IGNORES clinvar                          reads ONLY clinvar
                        |                                         |
                        v                                         v
              GOF / DN evidence for a                  hpo_condition_mechanism
              query variant                            _cache.json
```

The two consumers read **disjoint** parts of the same file. ClinVar contributes
conditions and never mechanism: ClinVar pathogenicity establishes neither gain of
function nor a dominant-negative effect, so the variant matcher does not look at
the `clinvar` block at all.

The injection is a separate step so that refreshing ClinVar does not force the
VEP pipeline to rerun.

## What establishes a record

A GoFCards record is usable for production variant matching when **either**
route is fully established:

| Route | Required |
| --- | --- |
| **1 — transcript route** | assembly + gene + transcript ID + HGVSc *or* HGVSp |
| **2 — coordinate route** | assembly + chromosome + position + reference + alternate |

Everything else is affiliated annotation. Both routes are build-bound, so
assembly sits above everything it qualifies.

## The variant identifier

```
loc_<chrom>:<start>:<ref>-><alt>_grch37
```

Examples across every variant shape:

```
loc_9:5073770:G->T_grch37             single substitution      (JAK2 V617F)
loc_7:140453136:AC->CT_grch37         compound / in-cis MNV    (BRAF V600R)
loc_6:91261902:-->TACTAC_grch37       insertion, no reference
loc_7:117199647:TTT->-_grch37         deletion, no alternate
```

### Why this construction

GoFCards publishes **no identifier that is both complete and unique**. Its
`Accession` field is a ClinVar cross-reference — verified against ClinVar, all
values are ClinVar VariationIDs — and it covers only 1,245 of 2,033 alleles
(61.2%), because 788 of the variants are not in ClinVar. It was never intended
to trace GoFCards' own records.

Identity therefore has to be constructed. Chromosome and position alone are not
enough:

| Candidate | Distinct | Verdict |
| --- | ---: | --- |
| chrom + start | 1,930 | **collides — 104 variants lost** |
| chrom + start + end | 1,938 | **collides — 96 variants lost** |
| chrom + start + ref + alt | **2,034** | unique |

The collisions are concentrated exactly on the hotspot codons that matter most
in a gain-of-function catalogue, where one base carries several substitutions:

```
chr17:7577120   TP53    C>T, C>A, C>G
chr3:178952085  PIK3CA  A>T, A>G, A>C
chr2:209113113  IDH1    G>A, G>T, G>C
chr9:5073770    JAK2    G>T (V617F)  and  G>A (V617I)
```

`end` is omitted because `ref` already encodes the span, and adding it recovers
only 8 of the 104 collisions.

**The `grch37` suffix is deliberate.** GoFCards states coordinates on GRCh37
only, so the identifier is minted from that build and says so. It is an opaque
label: never parse it, and never read a GRCh38 position out of it. The GRCh38
position lives under `assemblies.hg38.genomic`, where its build is named.

Missing alleles are written `-`, exactly as GoFCards writes them. The string
still parses unambiguously: `"-->TACTAC".split("->")` yields `["-", "TACTAC"]`.

As a side effect the identifier repairs a real source inconsistency — the `chr`
column mixes integer `4` and string `"4"` in the same file, which would
otherwise split 13 variants into phantom duplicates.

## Key hierarchy

```
gofcards_exact_gof.json.gz
│
├── metadata                                        provenance, stored once
│
└── genes
    └── "<GENE SYMBOL>"                             ── LEVEL 1
        ├── hgnc_id
        └── variants
            └── "loc_<chrom>:<start>:<ref>-><alt>_grch37"    ── LEVEL 2
                │
                ├── record                          build-independent, stored once
                │   ├── eligibility                 the gate verdict
                │   ├── source                      what GoFCards published, verbatim
                │   ├── liftover_status
                │   ├── annotations                 optional extras (rsid, …)
                │   └── evidence                    one entry per publication
                │
                └── assemblies
                    └── "hg19" | "hg38"             ── LEVEL 3
                        │
                        ├── genomic                         ROUTE 2
                        │   { chrom, pos, ref, alt, status }
                        │   one object, not a map: a variant has exactly
                        │   one position per build
                        │
                        └── transcripts                     ROUTE 1
                            └── "<TRANSCRIPT ID + VERSION>" ── LEVEL 4
                                ├── by_hgvsc
                                │   └── "<c.HGVS>"          ── LEVEL 5  unique
                                │       ├── hgvsp           essential (may be null)
                                │       ├── consequence     affiliated
                                │       ├── canonical       affiliated
                                │       └── mane_select     affiliated
                                └── by_hgvsp
                                    └── "<p.HGVS>" : [ "<c.HGVS>", … ]   NOT unique → list
```

### Why each key is what it is

| Level | Key | Reason |
| --- | --- | --- |
| 1 | gene symbol | What PriVA queries by. `hgnc_id` is carried as a value so a rename is detectable. |
| 2 | constructed variant ID | The variant is the thing that owns the evidence and the eligibility verdict, and both are build-independent. |
| 3 | assembly | Both establishment routes are build-dependent. Nothing build-derived sits above this. |
| 4 | versioned transcript ID | The version matters — the same accession has different coding sequence between GENCODE 19 and GENCODE 47. |
| 5 | HGVSc | Measured unique: 25,106 keys, **0 collisions**. |
| — | `by_hgvsp` → **list** | Measured **not** unique: 135 collisions. JAK2 `c.1605G>C` and `c.1605G>T` both give `p.Met535Ile`. |

Note that there is no coordinate-keyed map. A variant has one position per build,
so `genomic` is a plain object. A lookup by coordinate is built in memory by the
reader, which keeps the file normalized and free of duplicated indexes.

## Compound variants

44 variants have both a reference and an alternate longer than one base.
Evaluated on the MANE transcript:

| Protein consequence | Count |
| --- | ---: |
| single residue substitution | **36** |
| frameshift | 3 |
| protein-altering `delins` | 2 |
| in-frame insertion | 2 |
| missense `delins` | 1 |

**36 of 44 collapse to a single amino-acid substitution** — two in-cis
substitutions inside one codon giving one residue change:

```
gene      ref   alt    protein change on the MANE transcript
BRAF      AC -> CT     p.Val640Arg
BRAF      AC -> TT     p.Val640Lys
MPL       TG -> AA     p.Trp515Lys
MPL       TGG-> AAC    p.Trp515Asn
H3-3A     GGA-> TGG    p.Gly34Trp
JAK2      AA -> CT     p.Lys539Leu
CHRNA4    TG -> AA     p.Ser424Phe
FGFR3     GG -> AA     p.Gly382Lys
```

The remaining eight, in full:

| Gene | Variant ID | HGVSc | HGVSp | Consequence |
| --- | --- | --- | --- | --- |
| CDKN1C | `loc_11:2905326:ACTTC->CAGCTC_grch37` | `c.824_825delinsGCT` | `p.Lys275SerfsTer79` | frameshift |
| CSF1R | `loc_5:149433642:TGAT->TGATTGAT_grch37` | `c.2906_2909dup` | `p.Phe971SerfsTer7` | frameshift |
| KAT6B | `loc_10:76790205:CA->CACA_grch37` | `c.5623_5624dup` | `p.Gln1875HisfsTer5` | frameshift |
| JAK2 | `loc_9:5070023:CACAAA->CTT_grch37` | `c.1614_1617delinsT` | `p.Lys539del` | protein_altering |
| KIT | `loc_4:55593602:GTGGAAGGTTGTTGAG->ACAA_grch37` | `c.1669_1683delinsCAA` | `p.Trp557_Glu561delinsGln` | protein_altering |
| KIT | `loc_4:55595517:TAC->GAT_grch37` | `c.2007_2009delinsGAT` | `p.Ile669_Thr670delinsMetIle` | missense |
| KIT | `loc_4:55592180:GCCTAT->GCCTATGCCTAT_grch37` | `c.1504_1509dup` | `p.Ala502_Tyr503dup` | inframe_insertion |
| LRP5 | `loc_11:68125138:GGGGTG->GGGGTGGGGGTG_grch37` | `c.509_514dup` | `p.Gly171_Glu172insGlyGly` | inframe_insertion |

### Which compound variants can use the coordinate route

Classified on the **normalized** alleles, which is what matching actually sees:

| | Count |
| --- | ---: |
| single event on both builds — coordinate route usable | 5 |
| single event on one build only | 1 |
| irreducible MNV on both builds — coordinate route not usable | 38 |

The five that reduce cleanly are tandem duplications and deletions written out
longhand in the source; normalization recovers the single-event form, and both
builds agree:

```
gene    variant id                                    hg19            hg38
CSF1R   loc_5:149433642:TGAT->TGATTGAT_grch37         C->CTGAT        C->CTGAT
KAT6B   loc_10:76790205:CA->CACA_grch37               G->GCA          G->GCA
KIT     loc_4:55592180:GCCTAT->GCCTATGCCTAT_grch37    C->CTGCCTA      C->CTGCCTA
KIT     loc_4:55593603:TGGA->AA_grch37                TGG->A          TGG->T
LRP5    loc_11:68125138:GGGGTG->GGGGTGGGGGTG_grch37   T->TGGGGTG      T->TGGGGTG
```

The 38 irreducible ones stay as multi-base substitutions on both builds, for
example `7:140453136 AC->CT` and `7:140753336 AC->CT`. The key is well formed but
a standard caller would emit the substitutions separately, so it is unreachable
in practice. **For these, do not attempt the coordinate route at all.**

**Every one of the 44 carries a well-formed HGVSp**, so the transcript route
reaches all of them.

Note that three records which look like plain multi-base substitutions are
really tandem duplications written out longhand — `CSF1R TGAT->TGATTGAT`,
`KIT GCCTAT->GCCTATGCCTAT`, `LRP5 GGGGTG->GGGGTGGGGGTG` — which VEP renders as
`dup`. This is why the identifier writes alleles verbatim instead of trying to
classify them.

## Matching order

**Match by HGVSp first. A valid HGVSp match is sufficient to establish the
gain-of-function effect**, and no agreement between genomic representations is
required. HGVSc is the next route, and coordinates are the fallback.

This ordering is what makes compound variants work. BRAF V600R sits at
`7:140453136` on GRCh37 and `7:140753336` on GRCh38, yet is `p.Val600Arg` on
both: the protein change is stable across builds where the coordinate is not,
and it is stable across genomic representations of the same protein consequence.

## Worked example — BRAF V600R, a compound variant

```json
"BRAF": {
  "hgnc_id": "HGNC:1097",
  "variants": {
    "loc_7:140453136:AC->CT_grch37": {
      "record": {
        "eligibility": {
          "status": "eligible", "gene_match_status": "gene_concordant",
          "reason": null, "vep_symbol": "BRAF"
        },
        "source": {
          "gofcards_allele_key": "Indel|7|140453136|140453137|AC|CT",
          "variant_type_label": "Indel", "assembly": "hg19",
          "chrom": "7", "start": "140453136", "ref": "AC", "alt": "CT",
          "protein_change": "BRAF:NM_004333:exon15:c.1798_1799delinsAG:p.V600R,…"
        },
        "liftover_status": "mapped",
        "annotations": { "rsid": null },
        "evidence": [ { "pmid": "…", "pscore": "…", "…": "…" } ]
      },
      "assemblies": {
        "hg19": {
          "genomic": { "chrom": "7", "pos": 140453136, "ref": "AC", "alt": "CT",
                       "status": "raw_ref_alt" },
          "transcripts": { "…": "…" }
        },
        "hg38": {
          "genomic": { "chrom": "7", "pos": 140753336, "ref": "AC", "alt": "CT",
                       "status": "lifted_ref_match" },
          "transcripts": {
            "ENST00000646891.2": {
              "by_hgvsc": {
                "c.1798_1799delinsAG": {
                  "hgvsp": "p.Val600Arg",
                  "consequence": "missense_variant",
                  "canonical": true,
                  "mane_select": "NM_004333.6"
                }
              },
              "by_hgvsp": { "p.Val600Arg": ["c.1798_1799delinsAG"] }
            }
          }
        }
      }
    }
  }
}
```

## ClinVar annotation (injected file only)

Added by the injection step under the variant it matched. Matching is by exact
normalized allele on either build, so a variant that failed liftover is still
reachable through its GRCh37 coordinates.

```
"loc_9:5073770:G->T_grch37": {
  "record":     { ... },
  "assemblies": { ... },
  "clinvar": {                                    <-- present only when matched
    "vcv_accession": "VCV000014986",
    "variation_id": "14986",
    "variation_name": "NM_004972.4(JAK2):c.1849G>T (p.Val617Phe)",
    "hgvs": [ "NM_004972.4:c.1849G>T", ... ],
    "condition_assertions": [
      { "condition": "...", "identifiers": {...},
        "review_status": "...", "clinical_significance": "..." }
    ]
  },
  "clinvar_additional": [ { ... } ]               <-- only if several VCV records match
}
```

| Key | Notes |
| --- | --- |
| `vcv_accession` | ClinVar VCV accession |
| `variation_id` | ClinVar VariationID — the same identifier GoFCards stores in its own `Accession` field for the 61% of variants that are in ClinVar |
| `hgvs` | HGVS expressions ClinVar lists for the record |
| `condition_assertions` | what the HPO condition cache is built from |
| `clinvar_additional` | further VCV records for the same variant, kept rather than overwritten |

A variant with no `clinvar` key simply has no ClinVar record at the review-star
threshold. It remains fully eligible for variant matching, which does not depend
on ClinVar in any way.

Injection counters are recorded under `metadata.clinvar`.

## Field classification

### Essential — no match without these

gene symbol (level 1) · assembly (level 3) · transcript ID with version (level 4)
· HGVSc (level 5) · HGVSp (value, and key of `by_hgvsp`) · chromosome, position,
reference, alternate (`genomic`).

### Affiliated annotation — optional

`consequence`, `canonical`, `mane_select`, `genomic.status`,
`record.annotations.rsid`, and every field inside `record.evidence[]`.

### Gate

`record.eligibility.status`. **Only `eligible` variants may be used as
gain-of-function evidence.** Every other state is retained for audit.

## Complete value sets

### `record.eligibility.status`

| Value | Count |
| --- | ---: |
| `eligible` | 1,999 |
| `quarantined_gene_discordance` | 19 |
| `quarantined_reviewed_lof` | 5 |
| `quarantined_reviewed_mixed` | 3 |
| `quarantined_allele_gene_discordance` | 2 |

Permitted but unused in this build: `quarantined_missing_source_gene`,
`quarantined_reviewed_dominant_negative`, `quarantined_reviewed_uncertain`,
`quarantined_reviewed_exclude`, `quarantined_mechanism_review_required`.

### `record.eligibility.gene_match_status`

`gene_concordant` (2,007) · `gene_discordant` (19) · `source_gene_only` (2)

### `record.eligibility.reason`

| Value | Count |
| --- | ---: |
| `null` | 1,999 |
| `curated_gene_absent_from_locus` | 19 |
| `reduced_pump_capacity:LOF` | 3 |
| `article_mechanism_leakage:LOF` | 2 |
| `mixed_functional_effects:MIXED` | 2 |
| `sibling_of_gene_discordant_assertion` | 2 |
| `tissue_dependent_polarity:MIXED` | 1 |

### `record.liftover_status`

`mapped` (2,027) · `unmapped` (1)

### `genomic.status`

| Value | Applies to |
| --- | --- |
| `raw_ref_alt` | hg19 — both alleles given and matching the genome |
| `deletion_padded_ref_match` | hg19 — sparse deletion, anchor base from GRCh37 |
| `insertion_padded` | hg19 — sparse insertion, anchor base from GRCh37 |
| `lifted_ref_match` | hg38 — produced by the chain-file liftover |

Variants whose reference did not match GRCh37 never enter this file; they go to
`gofcards_reference_rejects.tsv` in the build working directory.

### `source.variant_type_label`

`SNV` (1,929) · `Indel` (99). **Provenance only** — the public workbook and the
backend disagree on this label for about 50 multi-base substitutions, which is
why it appears in no key.

### `by_hgvsc.<c.>.consequence`

53 distinct VEP values. Twelve most frequent: `missense_variant` (12,461),
`non_coding_transcript_exon_variant` (5,900),
`3_prime_UTR_variant,NMD_transcript_variant` (2,112), `downstream_gene_variant`
(1,225), `intron_variant` (1,113), `upstream_gene_variant` (919), `stop_gained`
(429), `missense_variant,NMD_transcript_variant` (411),
`intron_variant,non_coding_transcript_variant` (386),
`missense_variant,splice_region_variant` (368),
`intron_variant,NMD_transcript_variant` (313), `5_prime_UTR_variant` (300).

### `record.evidence[]`

One entry per publication, never merged — the basic unit of GoFCards is an
assertion (a variant plus one paper), not a variant. 3,161 workbook rows describe
2,034 variants; 472 variants carry several rows and one carries 45.

| Key | Notes |
| --- | --- |
| `pmid` | single PubMed identifier |
| `pscore` | GoFCards P score for this record |
| `disease` | GoFCards "Disorder involved" |
| `function` | free text — **not reliably scoped to this variant**; see the limitations section of `/paedyl01/disk1/yangyxt/PriVA/docs/GOFCARDS_NORMALIZATION_METHODS.md` |
| `pathway` | GoFCards "Pathways proteins involved" |
| `animal_model`, `cell_model` | `Y` / `N` |
| `source_order` | row number in the public workbook |
| `source_transcript` | RefSeq transcript GoFCards cited |
| `source_protein_change` | ANNOVAR protein string for this record |

## Scale

Measured on a full build (2026-07-27). 1.05 MB compressed, 15.1 MB uncompressed.

| Quantity | Count |
| --- | ---: |
| genes | 575 |
| variants | 2,028 |
| distinct variant IDs | 2,028 (unique) |
| variants present on hg19 | 2,026 |
| variants present on hg38 | 2,026 |
| transcript views (variant x assembly x transcript) | 65,105 |
| HGVSc keys | 65,105 |
| HGVSp keys | 51,482 |
| evidence entries | 3,154 |
| variants carrying an rsID | 1,438 |

Every transcript is retained. Collapsing transcripts that share an HGVS string
would delete a lookup path, because a MANE pair yields the same HGVSc from its
Ensembl and its RefSeq member and PriVA may annotate with either.

## Querying

**Protein-change route — try this first:**

```python
import gzip, json
cache = json.load(gzip.open(
    "/paedyl01/disk1/yangyxt/PriVA/data/gofcards/gofcards_exact_gof.json.gz", "rt"))

for vid, variant in cache["genes"]["BRAF"]["variants"].items():
    if variant["record"]["eligibility"]["status"] != "eligible":
        continue
    view = variant["assemblies"].get("hg38", {}).get("transcripts", {}).get("ENST00000646891.2")
    if view and "p.Val600Arg" in view["by_hgvsp"]:
        ...          # usable as gain-of-function evidence
```

**Coordinate route — fallback only:**

```python
for vid, variant in cache["genes"]["BRAF"]["variants"].items():
    g = variant["assemblies"].get("hg38", {}).get("genomic")
    if g and (g["chrom"], g["pos"], g["ref"], g["alt"]) == ("7", 140753336, "AC", "CT"):
        if variant["record"]["eligibility"]["status"] == "eligible":
            ...
```

Both routes reach the same `record`, so the gate and the literature are identical
whichever way the variant was found.

## Known defects

1. **A protein-change agreement field was removed.** It compared our three-letter
   form (`VAL617PHE`) against GoFCards' one-letter ANNOVAR form (`p.V617F`),
   which can never match, so it reported disagreement for every record. The raw
   source string is preserved at `record.source.protein_change` so the comparison
   can be done properly later.

2. **Six variants have normalized alleles that differ between the two builds.**
   Both builds are normalized against their own reference, so a genuine sequence
   difference between GRCh37 and GRCh38 produces a different minimal
   representation:

   | Gene | Source | hg19 | hg38 |
   | --- | --- | --- | --- |
   | ABCC9 | `TCTACAGTGA->G` | `TCTACAGTGA->G` | `TCTACAGTGA->T` |
   | BRAF | `C->TGTA` | `C->TGTA` | `C->CGTA` |
   | CDKN1C | `ACTTC->CAGCTC` | `ACT->CAGC` | `CT->AGC` |
   | JAK2 | `CACAAA->CTT` | `ACAAA->TT` | `CAAA->T` |
   | KIT | `GTGGAAGGTTGTTGAG->ACAA` | `GTGGAAGGTTGTTGAG->ACAA` | `TGGAAGGTTGTTGAG->CAA` |
   | KIT | `TGGA->AA` | `TGG->A` | `TGG->T` |

   Three of these sit in the KIT region already known to carry reference
   discrepancies. The protein route is unaffected, which is another reason it is
   tried first.
