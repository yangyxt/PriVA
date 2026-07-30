# PriVA Gene Mechanism Hub Decision Flow

![PriVA gene mechanism and inheritance decision flow](GENE_MECHANISM_HUB_DECISION_FLOW.svg)

The diagram represents the current implementation in
`/paedyl01/disk1/yangyxt/PriVA/scripts/gene_mechanism_hub.py` and its call path
from `/paedyl01/disk1/yangyxt/PriVA/scripts/acmg_criteria_assign.py`.

## Source key

The most important distinction in this table is which sources reach the hub and
which do not. Only S1, S5, S6 and S7 can create a condition history. S2, S3, S4
and S8 are gene-wide: none of them says which condition a variant acts on, nor by
what mechanism, so none can create a history. They are read by PVS1's gene gate
alone and never enter the hub.

| ID | Source or resource | Current role |
|---|---|---|
| S1 | `/paedyl01/disk1/yangyxt/PriVA/data/hpo/genes_to_phenotype.assertions.tsv.gz` | Gene to condition, with the condition's inheritance, penetrance and onset assertions. The backbone of the cache. |
| S5 | `/paedyl01/disk1/yangyxt/PriVA/data/patho_mechanism/gene_pathogenic_mechanism_evidence.tsv` | Curated G2P/DDG2P mechanisms and their allelic requirements. |
| S6 | `/paedyl01/disk1/yangyxt/PriVA/data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.json.gz` | PriVA-local canonical schema-v2 non-LOF gene-condition, GoFCards, and exact matched ClinVar VCV evidence. This is the only mechanism cache the hub reads; the older combined cache it once fell back to is no longer published. |
| S7 | `/paedyl01/disk1/yangyxt/PriVA/data/gofcards/gofcards_exact_gof.json.gz` | Build-time GoFCards normalization and quarantine. The canonical S6 JSON embeds eligible exact alleles. |
| S9 | Orphadata | Curated mechanisms per condition, admitted alongside S5. |
| S2 | `/paedyl01/disk1/yangyxt/PriVA/data/clingen/gene_dosage_sensitivity.hg19.tsv` | Two separate roles. Haploinsufficiency scores 3, 2 and 1 are admitted **into the cache** as a curated LOF mechanism, keyed by MONDO. Scores 3 and 30/40 are also read **by PVS1's gate** as a gene-wide signal. |
| S3 | `/paedyl01/disk1/yangyxt/PriVA/data/loeuf/loeuf_dataset.tsv.gz` | `LOEUF < 0.35`. Read by PVS1's gate, and by the hub's third inheritance fallback tier only. |
| S4 | `/paedyl01/disk1/yangyxt/PriVA/data/alphamissense/alphamissense_mean_score.tsv` | Per-transcript `mean_am_pathogenicity`, mapped to genes. `> 0.6` is read by PVS1's gate. |
| S8 | `clinvar_pathogenic_genes`, computed at run time | Genes carrying a ClinVar pathogenic or likely-pathogenic variant at two stars or better. Read by PVS1's gate. |

`gene_nonlof_mechanism_json` selects S6. It is the only mechanism-cache
configuration key; the older `gene_mechanism_json` fallback has been removed
rather than repointed, because the non-LOF cache is tracked in this repository
and so that choice had one possible outcome.

## Implemented variant-level contract

```text
gene-condition assertion, from the cache
    mechanism: LOF / GOF / DOMINANT_NEGATIVE / TRIPLOSENSITIVITY
               (a history with no mechanism is labelled UNRESOLVED)
    assertion_basis: curated | deduced
    inheritance, penetrance, disease, source, confidence

query-variant mechanism, three scores of 0 / 1 / 2
    2  exclusively established -- an exact curated allele match, or a
       transcript destroyed by nonsense-mediated decay. When any mechanism
       scores 2 the other two are forced to 0.
    1  plausible from consequence or prediction only
    0  confidently not this mechanism

variant_effect
    exact_known_LOF                 decay triggered. Outranks even a curated
                                    call: there is no protein left to gain a
                                    function
    exact_known_GOF                 an exact curated allele
    exact_known_DOMINANT_NEGATIVE   "+"-joined when both apply
    predicted_LOF_high_confidence   2 with a LOFTEE HC rescue, otherwise 1
    uncertain                       no mechanism evidence at all

variant-to-condition applicability
    applicable  the variant's mechanism is established for that history
    uncertain   the mechanism is possible but not established
```

`assertion_basis` records whether a curator asserted the mechanism or the cache
build deduced it. The build writes `deduced` LOF for a recessive condition that
states no mechanism from any source. That deduction carries inheritance
everywhere, but PVS1's gate counts only `curated` entries, because the deduction
restates the inheritance the same record already supplies and so adds no
independent observation to the one criterion that can reach Very Strong alone.

Gene-wide signals do not promote a bare `recessive` tag to `recessive_LOF`. That
older behaviour is gone: those signals never enter the hub at all.

`var_plausible_patho_mechs` examples are `recessive`, `recessive_LOF`,
`recessive_GOF`, `dominant`, `dominant_LOF`, `dominant_GOF`, and `dominant_DN`.

## The fifteen columns attached to every row

```text
variant_effect                  variant_lof_score
variant_gof_score               variant_dn_score
variant_mechanism_exclusive     variant_exact_mechanisms
variant_mechanism_applicable    variant_mechanism_uncertain
variant_condition_ids           variant_condition_histories
variant_inheritance             variant_inheritance_basis
variant_x_linked                variant_penetrance
gene_lof_mechanism_history
```

`variant_condition_histories` keeps one entry per history as
`condition|mechanism|inheritance|penetrance`, so the pairing between a condition
and its inheritance and penetrance survives. De-duplicating those three as
separate lists would make them different lengths and destroy the pairing.

`variant_inheritance_basis` says which of three tiers answered:
`matched_history` (3,897 genes, 70.0%), `gene_consensus` (879, 15.8%), or
`gene_constraint` (789, 14.2%). No gene is left without an inheritance.

`gene_lof_mechanism_history` is a property of the gene rather than of the
variant, and is read only by PVS1.

PVS1, BS1, BS2, BS4, PP1, BP5 and BP2/PM3 require the variant-level fields.
Missing fields raise an error rather than activating a gene-level or raw-
annotation fallback.

## LOFTEE HC and OS

Both count as predicted loss of function, but they are no longer graded the same.
`OS` means "other splice" and covers predicted disruptions in extended donor or
acceptor regions as well as newly created donor sites.

`HC` scores 2, and is the one signal that rescues a truncating variant which
escapes nonsense-mediated decay. `OS` scores 1. The original LOFTEE category
remains in `variant_effect_evidence`, so OS is not collapsed into HC.

A curated allele still outranks a LOFTEE `HC` call: those are predictions, this
is a curator's verdict on this exact allele. Whatever the winner outranked is
recorded in `variant_effect_suppressed_evidence` rather than dropped.

## PVS1 asks two questions and never mixes them

PVS1 is the only criterion that can carry a variant to Likely Pathogenic on its
own, so its two halves are kept strictly apart.

**Is loss of function a disease mechanism of this gene?** Purely a question about
the gene. It reads none of the variant's mechanism scores. Any one of five
independent gene-level sources opens the gate, deliberately a union rather than a
consensus, because each is incomplete on its own and a gene missing from one is
routinely present in another:

1. `gene_lof_mechanism_history` — a **curated** LOF in the condition cache
2. a ClinVar pathogenic or likely-pathogenic variant at two stars or better
3. `LOEUF < 0.35`
4. mean AlphaMissense `> 0.6`
5. ClinGen haploinsufficiency score 3, or 30/40 for a recessive phenotype

**How convincing is this particular null variant?** Decided from the consequence
— whether decay is triggered, whether an intolerant domain is spanned, what
fraction of the protein is lost — and not from any mechanism score. Letting the
mechanism layer pre-compute a verdict here would decide for PVS1 a question that
is PVS1's own.

The single exception runs the other way: `is_exact_gof` withdraws PVS1 from a
variant curated as gain of function. It takes evidence away, never grants it.

### Source 2 reads both ClinVar structures

PriVA's ClinVar build writes two files in one pass, keyed on different things.
The amino-acid file is keyed on the protein consequence, and every entry has a
non-blank HGVSp. Frameshift and nonsense variants are there, carrying `p.*fs` and
`p.*Ter`, but a canonical splice-site variant such as `c.1234+1G>A` produces no
HGVSp in VEP and cannot appear at all — nor can deep-intronic variants or
whole-exon deletions. The splice file is keyed on having an `EXON` or `INTRON`
field, so those variants do appear. Reading only the first silently missed any
gene whose two-star pathogenic variants all lack a protein consequence.

`summarize_clinvar_gene_pathogenicity` reads both, applying one strict test to
each: `CLNSIG` split on `,` `/` `|` `;`, lowercased, and intersected with
`{pathogenic, likely_pathogenic}`, with review status at two stars or better.
Tokenising is what keeps `Conflicting_classifications_of_pathogenicity` out.

Measured on the hg19 build: 2,869 genes to 3,040, and nothing is lost. `VMA21` is
the clearest case recovered — X-linked myopathy with excessive autophagy, whose
pathogenic variants sit in the untranslated and first-intron region, so the
protein-keyed file could never see it.

## ClinVar hierarchy and review-tier audit

The exact allele join is made at the VCV level because VCV is ClinVar's stable
variation-centric record. Each matched VCV retains its condition-specific RCV
records as `condition_assertions`; linked submission-level SCV observations
remain nested under the corresponding RCV. The review threshold is applied to
the RCV aggregate classification, not to the VCV container.

For an exact query GoFCards match, PriVA now joins the row's GoFCards variant or
accession ID to the nested schema-v2 `ClinVar_VCV.match` records. This adds the
VCV accession, RCV condition names, maximum retained review stars, and ClinVar
HGVS expressions. ClinVar observed zygosity remains evidence only and is not
converted into a causal allelic requirement.

The 2026-07-23 full-release audit found:

| Review coverage among exact, gene-concordant GoFCards matches | VCV count |
|---|---:|
| Maximum RCV review level is 0 stars | 207 |
| Maximum RCV review level is 1 star | 626 |
| Maximum RCV review level is 2 stars | 445 |
| Maximum RCV review level is 3 stars | 80 |
| Total with at least one RCV at 2 or more stars | 525 |

The two-star-or-higher canonical subset covers 216 distinct genes. Lower-review
matched VCVs occur in 330 genes; 208 of those genes have no matched VCV reaching
two stars. Review stars measure review quality and consensus, not whether the
classification itself is pathogenic or benign.

The review-tier audit is generated on demand by
`/paedyl01/disk1/yangyxt/PriVA/scripts/audit_clinvar_gofcards_review_tiers.py` and is not kept in the
repository. It writes a summary and the full matched-RCV table into
`/paedyl01/disk1/yangyxt/PriVA/data/patho_mechanism/`, where build byproducts are ignored by git.
