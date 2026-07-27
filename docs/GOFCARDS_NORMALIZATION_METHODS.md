# GoFCards normalization: methods, quality control, and limitations

## Scope

This document describes how PriVA turns the GoFCards catalogue into variant-level
evidence. It deliberately separates two questions that must not be conflated:

1. **Allele identity:** are the gene, coordinates, alleles, transcript, and HGVS
   representation internally consistent and consistent with the genome?
2. **Mechanism polarity:** does the cited experiment establish gain of function
   (GOF), loss of function (LOF), a mixed effect, or no resolved mechanism?

The workflow establishes allele identity. It does not independently establish
mechanism polarity. Membership in the GoFCards catalogue is therefore a source
claim to be audited, not by itself a GOF classification.

## Source

- Public workbook: `https://download.genemed.tech/upload/GainFunCards/gofcards_data_download.xlsx`
- Backend tables: `GainFunCards_SNV` and `GainFunCards_Indel` under
  `https://java.genemed.tech/admin-api/backend/data/hg19`

The public workbook is the citable artifact. The backend tables are pulled
because they carry the ANNOVAR protein change (`AAChange_refGene`) that the
workbook omits; that field is what lets us detect alleles where GoFCards' own
protein numbering disagrees with ours. Three network requests in total.

**GoFCards publishes GRCh37 coordinates only.** The public workbook has `chr`,
`hg19start`, `hg19end`, `ref`, and `alt`, and no GRCh38 field. GRCh37 is
therefore the sole ground truth, and every other representation in the cache is
derived from it.

The GoFCards publication states that PubMed was searched for articles whose
title or abstract contained `gain of function` or `GOF` together with `variant`
or `variation`, that review articles were excluded, and that included variants
were manually curated (PMID: 39578693; DOI: 10.1093/nar/gkae1079). The workbook
has no explicit row-level mechanism field.

## Workflow

The whole pipeline is `scripts/build_gofcards_exact_gof_cache.py`, deployed by
`gofcards_exact_gof_cache_install` in `scripts/install_utils.sh`. Seven steps,
one function each:

| Step | Function | Job |
| --- | --- | --- |
| 1 | `fetch_sources` | download the public workbook and both backend tables |
| 2 | `build_source_vcf` | pad sparse indels, check REF against GRCh37, reject failures |
| 3 | `normalize_vcf` | `bcftools norm` against GRCh37 |
| 4 | `run_vep` | Ensembl VEP on GRCh37, all transcripts |
| 5 | `liftover_vcf` + `normalize_vcf` | CrossMap to GRCh38, then re-left-align |
| 6 | `run_vep` | Ensembl VEP on GRCh38, all transcripts |
| 7 | `build_cache` | join, resolve HGNC symbols, apply reviews, export |

`normalize_vcf` and `run_vep` are each a single function invoked twice, once per
assembly, rather than an assembly-specific pair.

### Why the two checks are separate

The **reference check** in step 2 decides which records are trustworthy. GoFCards
asserts a reference base at a GRCh37 position; only GRCh37 can adjudicate that
claim. An allele whose reference does not match describes a variant that does not
exist, so its HGVS is meaningless regardless of how confidently VEP reports it.
Those records are rejected before any annotation happens.

This check cannot be deferred to the liftover. CrossMap rewrites the reference
allele by reading the target genome, so it will silently convert a wrong GRCh37
reference into a plausible-looking GRCh38 one.

The **liftover** in step 5 decides only which records can offer a GRCh38
coordinate key. Failing it costs an allele one match route, not its eligibility:
its HGVS came from the GRCh37 run and is unaffected. Unmapped records are counted
and reported.

Re-running `bcftools norm` after the liftover is required because repeat
structure differs between assemblies, so an indel correctly left-aligned in
GRCh37 is not guaranteed to be left-aligned in GRCh38.

### Why VEP runs twice

HGVS is a property of the transcript sequence, not of the genome assembly. The
same variant on the same transcript version yields the same HGVSc and HGVSp from
either assembly; measured across this dataset, agreement was 9,701 to 1.

What differs is the transcript catalogue each VEP cache ships: GENCODE 19 with
the GRCh37 cache, GENCODE 47 with the GRCh38 cache. The same transcript accession
usually carries a different version between them, and a different version means a
different coding sequence and different residue numbering. Of transcripts present
in both runs, 623 gave a different protein change, and every one of those was a
version difference.

The two runs therefore exist so that PriVA can match against whichever catalogue
it annotated its own query variants with, using the same VEP cache directory that
`config.yaml` points at. Because both sides read the same cache, transcript
versions align by construction rather than by coincidence.

All transcripts are retained, not just a preferred one. Across this dataset, 66%
of variants have transcripts that disagree with each other about the protein
change, so keeping one would make a match depend on PriVA having chosen the same
transcript. Anchoring instead to GRCh38 MANE Select was tested and rejected: it
matches only 80.8% of GRCh37-annotated queries, against 99.6% when all
transcripts are kept.

### Evidence handling

Source records are aggregated per allele rather than reduced to one
representative. An allele can carry many independent literature records — KIT
p.Asp816Val has 33 — and retaining only one silently discards the rest, including
any that describe a different mechanism from the one retained.

## Eligibility

A row is eligible only if all of the following hold. Anything else is written to
the cache with an explicit `match_eligibility` and `reject_reason` rather than
being dropped silently.

1. The reference allele matched GRCh37 (step 2).
2. The resolved curated gene equals the resolved VEP gene. Both symbols pass
   through the HGNC complete set, so aliases and withdrawn symbols are normalized
   before comparison. A disagreement quarantines the whole curated gene-allele
   assertion, including concordant sibling transcript rows.
3. The allele has no reviewed polarity decision marking it non-GOF
   (`data/gofcards/gofcards_mechanism_reviews.tsv`).

Match routes, in order of preference: gene plus HGVSp, gene plus HGVSc, then
genomic coordinates for the assembly the query variant is mapped to. Coordinates
are an alternative route, not the default.

## Mechanism-polarity audit finding

The GoFCards P score measures the amount and type of supporting material: one
point for a literature report, two for a cellular model, three for an animal
model. It does not encode whether the observed effect is GOF or LOF. A
well-performed LOF experiment can therefore receive a nontrivial P score inside
the GOF catalogue.

A deterministic screen of all 3,161 public rows found 11 rows across eight genes
whose `Function` text explicitly contains LOF language:

- **Clear polarity conflicts:** CFTR p.Phe508del and TET2 p.Gln886Ter are LOF.
  CFTR is an in-frame deletion causing defective maturation and trafficking; TET2
  p.Gln886Ter is a stop-gain allele whose paper discusses GOF KIT D816V and LOF
  TET2 in the same neoplasm.
- **Variant-specific leakage within one paper:** the ATP2A2 study tested four
  alleles. p.Gly23Arg, p.Asp567Tyr, and p.Ile1014Val were LOF pumps; only
  p.Gly860Ser acquired a calcium-leak GOF.
- **Mixed or context-dependent alleles:** NLGN3 p.Arg451Cys, KCNQ1 p.Gln147Arg,
  and KCNA2 p.Phe302Leu have GOF and LOF language tied to different measured
  properties, tissues, or phenotypes.
- **GOF rows with background LOF discussion:** HNF1B p.Ser36Phe and ACTL6B
  p.Gly343Arg are GOF-relevant alleles whose text also discusses other mechanisms.

The likely failure mode is article-level mechanism leakage: the source search
selects papers containing GOF terminology, but a paper can discuss a GOF allele
alongside a different LOF allele, gene, or disease.

### This screen is a lower bound, not a census

The screen was a keyword search over the `Function` text, and that text is not
reliably scoped to the row's own allele. Of 1,896 eligible alleles with a
parseable protein change:

| | Alleles |
| --- | ---: |
| `Function` text names no protein change at all | 1,671 (88%) |
| text names the same change as the row | 163 (8.6%) |
| text names **only other** protein changes | 62 (3.3%) |

The 62 split into 29 isoform-numbering offsets, 4 same-position different
substitutions, and 29 genuinely unrelated alleles — for example VPS35 `G174S`
whose text describes `D620N`, and ERBB2 `L755S` whose text discusses `I767M`.

Because the text is not row-scoped, the keyword screen can both miss LOF alleles
whose text never uses the phrase and flag alleles because a different allele in
shared text is LOF. Eleven rows is therefore a lower bound produced by a method
that cannot establish an upper bound.

## Limitations and required remediation

1. Row-level polarity adjudication currently covers 11 alleles. It needs to cover
   every source allele, with `GOF`, `LOF`, `mixed`, `DN`, `uncertain`, and
   `exclude` states, and a recorded PMID, publication type, quoted evidence,
   curator, and review date for every non-GOF disposition.
2. The exact experimental endpoint is not stored, nor whether it measures the
   encoded protein directly or a downstream phenotype.
3. Polarity conflicts are quarantined rather than relabelled. That is deliberate:
   no replacement mechanism is asserted without a primary source.
4. Until polarity review is complete, manuscript language should state that
   GoFCards alleles were normalized and independently checked for gene and allele
   consistency, but that source GOF assertions were treated as candidates
   requiring polarity review.

## Key literature

- Zhao W, et al. GoFCards: an integrated database and analytic platform for gain
  of function variants in humans. *Nucleic Acids Res.* 2025;53:D976-D988.
  PMID: 39578693; DOI: 10.1093/nar/gkae1079.
- Zhao JH, et al. Chemical chaperone and inhibitor discovery: potential treatments
  for protein conformational diseases. *Perspect Medicin Chem.* 2007;1:39-48.
  PMID: 19812735; DOI: 10.4137/pmc.s212.
- Cheng SH, et al. Defective intracellular transport and processing of CFTR is the
  molecular basis of most cystic fibrosis. *Cell.* 1990;63:827-834.
  PMID: 1699669; DOI: 10.1016/0092-8674(90)90148-8.
- Jensen TJ, et al. Multiple proteolytic systems, including the proteasome,
  contribute to CFTR processing. *Cell.* 1995;83:129-135. PMID: 7553864;
  DOI: 10.1016/0092-8674(95)90241-4.
- Ward CL, et al. Degradation of CFTR by the ubiquitin-proteasome pathway.
  *Cell.* 1995;83:121-127. PMID: 7553863; DOI: 10.1016/0092-8674(95)90240-6.
- Zhang W, et al. Restoration of Sarco/Endoplasmic Reticulum Ca2+-ATPase activity
  functions as a pivotal therapeutic target. *Front Pharmacol.* 2022;13:877175.
  PMID: 35517826; DOI: 10.3389/fphar.2022.877175.
- Chow S, et al. Diagnosis of systemic mastocytosis with cryptic deletion of TET2
  and DNMT3A resulting from unbalanced translocation. *Br J Haematol.*
  2024;205:961-966. PMID: 38702998; DOI: 10.1111/bjh.19501.

## Frozen provenance for the superseded build

The cache deployed before this workflow was produced on 2026-07-07 by a separate
repository that queried a per-allele GoFCards summary endpoint for GRCh38
coordinates. That endpoint is no longer used: it was the source of eleven GRCh38
records that failed reference validation, and everything else it supplied was
either duplicated elsewhere (`gofcards_variant_id` was identical to the allele
key for all rows), recoverable from VEP (`rsID`), or never exported (ClinVar
fields). Its inputs, VEP outputs, and audit logs are retained unchanged under
`/paedyl01/disk1/yangyxt/gofcards_hg38_normalizer/work_20260707_full/` and
`/paedyl01/disk1/yangyxt/gofcards_hg38_normalizer/audit/` as the audit trail for
that superseded build.
