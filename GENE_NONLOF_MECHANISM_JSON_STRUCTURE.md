# Canonical Non-LOF Mechanism JSON And ClinVar VCV Structure

## Authority And Scope

The only authoritative mechanism JSON remains:

`/paedyl01/disk1/yangyxt/PriVA/data/gene_pathogenic_mechanism/prepared/gene_nonlof_mechanism_curated_assertions.json`

ClinVar VCV evidence is injected into that file. The flattened TSV is only an
audit view, not a parallel authority:

- `/paedyl01/disk1/yangyxt/PriVA/data/gene_pathogenic_mechanism/prepared/clinvar_vcv_gofcards_matches.tsv`

The machine-readable output contract is:

`/paedyl01/disk1/yangyxt/PriVA/data/gene_pathogenic_mechanism/schema/gene_nonlof_mechanism_curated_assertions.schema.json`

## What VCV, RCV, And SCV Mean

ClinVar uses three accession levels:

| Accession | Plain-language unit | Meaning in this integration |
|---|---|---|
| `VCV` | variation record | One ClinVar variation, haplotype, or genotype record. It describes what the molecular change is and can contain several condition-specific RCV records. |
| `RCV` | variation-condition record | One aggregate relationship between that VCV and one condition or trait set. Review stars and aggregate germline classification are evaluated here. |
| `SCV` | submitter assertion | One laboratory, expert panel, or other submitter's classification and supporting observations. Several SCVs can contribute to one aggregate RCV classification. |

The conceptual relationship is:

```text
VCV: molecular variation
└── RCV: variation + condition A
│   ├── SCV: submission from laboratory 1
│   └── SCV: submission from laboratory 2
└── RCV: same variation + condition B
    └── SCV: submission from expert panel 3
```

An SCV is not assumed to belong to an RCV merely because both occur in the
same VCV XML record. The parser links them using condition identifiers/names
and ClinVar `TraitMapping` records. Ambiguous links to several eligible RCVs
are retained and marked explicitly.

## Eligibility Contract

A ClinVar condition is retained only when all of these statements are true:

1. The `VariationArchive` is `current` and contains a `ClassifiedRecord`.
2. A `SimpleAllele`, or a component allele of a haplotype/genotype, exactly
   matches a normalized GoFCards `hg19_vcf_key` or `hg38_vcf_key`.
3. The condition has an RCV-level `GermlineClassification`.
4. The condition-specific RCV `ReviewStatus` is at least two stars.

The allowed review statuses are exactly:

| Stars | Retained RCV `ReviewStatus` |
|---:|---|
| 2 | `criteria provided, multiple submitters, no conflicts` |
| 3 | `reviewed by expert panel` |
| 4 | `practice guideline` |

The filter is applied at
`RCVAccession/RCVClassifications/GermlineClassification/ReviewStatus`, not at a
VCV-wide classification. A single VCV can therefore retain one disease and
exclude another.

No clinical-significance filter is added after the star filter. Pathogenic,
benign, uncertain, and other eligible aggregate assertions remain visible as
condition-specific evidence. GoFCards supplies the GOF assertion; ClinVar does
not independently prove GOF.

## Complete Emitted JSON Shape

### Presence And Empty-Value Notation

The following notation is used in the trees and field tables:

- **Required**: the key must exist for that object to pass the JSON Schema.
- **Optional**: the key is omitted entirely when the source did not provide a
  usable value.
- **`[1..*]`**: a present array has at least one item.
- **`[0..*]`**: the array is present but is allowed to be empty.

The canonical serializer does not use JSON `null`, empty strings, or literal
placeholder strings such as `"NA"`, `"N/A"`, `"_"`, or `"-"` to represent
missing evidence. Missing scalar evidence means that the key is absent.
The sole source-notation exception is top-level `GoFCards.ref` or
`GoFCards.alt`: the public workbook uses `"-"` as a meaningful sparse empty
allele for an insertion or deletion. It is preserved there as evidence, not as
a missing-value marker. The nested normalized VCF REF/ALT fields are padded
and never use `"-"`.
`"not provided"` can still occur when it is an actual controlled ClinVar
value, for example a submitted zygosity category; in that case it is evidence
reported by ClinVar, not a JSON missing-value marker.

Two empty arrays have deliberate meanings in the current canonical output:

- `matched_scvs: []`: the RCV passed the review filter but no current,
  contributing germline SCV could be linked confidently to that condition.
- `discarded_gene_discordant_symbols: []`: exact matched GoFCards rows did not
  contain a gene symbol that needed to be discarded.

Optional collection fields such as `pmids`, `observations`, `trait_mappings`,
and `hgvs` are omitted when empty.

### Root And Gene Presence

| Object | Required keys | Optional keys | Collection rule |
|---|---|---|---|
| root | `_meta` | Any number of `HGNC:<number>` or fallback `SYMBOL:<symbol>` gene keys | At least zero gene keys are schema-valid; a production build normally has many. |
| `_meta` | `version`, `built_at`, `total_genes`, `sources` | `unmapped_symbols`, `mechanism_values`, `clinvar_vcv_integration` | The builder currently emits all three optional members; `unmapped_symbols` may be empty. |
| gene | `symbol`, `mechanisms` | `entrez_id`, `ensembl_id`, `gene_level`, `variant_level` | `mechanisms` is non-empty. `gene_level` and `variant_level` are omitted when they would be empty and contain at least one item when present. |
| `gene_level` item | exactly one of `G2P_DDG2P`, `PanelApp`, or `ClinGen_Dosage` | no second source key in the same item | One source assertion per array item. |
| `variant_level` item | exactly one of `GoFCards` or `ClinVar_VCV` | no second source key in the same item | One source assertion per array item. |

### `_meta.sources` Provenance Objects

These objects describe how to interpret each source. They do not themselves
assert that a gene or variant is pathogenic.

```text
_meta.sources
├── G2P_DDG2P
│   ├── level: "gene_level"
│   ├── keys: [mechanism, mechanism_raw, disease, inheritance, confidence, pmids]
│   ├── mechanism_raw: field-description string
│   └── confidence_values: field-description string
├── PanelApp
│   ├── level: "gene_level"
│   ├── keys: [mechanism, panel, panel_id, inheritance, confidence]
│   └── note: interpretation caveat
├── ClinGen_Dosage
│   ├── level: "gene_level"
│   ├── keys: [mechanism, score, score_description, pmids]
│   └── score: score-definition string
├── GoFCards
│   ├── level: "variant_level"
│   ├── keys: [mechanism, source_record_id, allele_key,
│   │          exact_normalization_status, exact_cache_gene_symbols,
│   │          exact_normalized_variants, disease, chr, pos, ref, alt,
│   │          transcript, consequence, pscore, pmids, function, pathway,
│   │          animal_model, cell_model]
│   ├── exact_normalized_variant_keys: [source, mechanism, build, HGNC_Symbol,
│   │   VEP_assembly, VEP_transcript, feature_type, consequence, HGVSc, HGVSp,
│   │   hgvsp_key, match_status, raw_GoFCards_HGVS, GoFCards_transcript,
│   │   canonical_transcript, hg19_chrom, hg19_pos, hg19_ref, hg19_alt,
│   │   hg19_vcf_status, hg38_chrom, hg38_pos, hg38_ref, hg38_alt,
│   │   hg38_refalt_status, gofcards_accession_id, gofcards_variant_id,
│   │   disease, pmids, pscore, function, pathway, derived_on, allele_key,
│   │   hg19_genomic_key, hg19_vcf_pos, hg19_vcf_ref, hg19_vcf_alt,
│   │   hg38_vcf_pos, hg38_vcf_ref, hg38_vcf_alt, hg19_vcf_key,
│   │   hg38_genomic_key, hg38_vcf_key, match_key_types]
│   ├── raw_public_allele_fields
│   │   ├── assembly: "hg19"
│   │   ├── keys: [chr, pos, ref, alt]
│   │   └── note: scope clarification; does not describe normalized records
│   ├── exact_match_methods: [gene + exact HGVSp,
│   │                         gene + exact hg19/hg38 genomic allele]
│   ├── pscore: score-definition string
│   ├── animal_model: flag-definition string
│   └── cell_model: flag-definition string
└── ClinVar_VCV
    ├── level: "variant_level"
    ├── scope: retained RCV/SCV evidence description
    ├── minimum_review_stars: 2
    ├── eligible_review_statuses: [three allowed status strings]
    ├── mechanism_caveat: ClinVar is not the GOF source
    └── inheritance_caveat: observations do not resolve allelic requirement
```

These metadata objects are emitted as a complete fixed documentation block by
the current builder. Their descriptive strings are not used as analysis
inputs; the actual assertions described below are the analysis inputs.

```text
root object
├── _meta
│   ├── version: "2.0"
│   ├── built_at: ISO-8601 string
│   ├── total_genes: integer
│   ├── unmapped_symbols: [string]
│   ├── mechanism_values: [string]
│   ├── sources
│   │   ├── G2P_DDG2P: source schema/provenance
│   │   ├── PanelApp: source schema/provenance
│   │   ├── ClinGen_Dosage: source schema/provenance
│   │   ├── GoFCards: source schema/provenance
│   │   └── ClinVar_VCV: review/matching/inheritance caveats
│   └── clinvar_vcv_integration: parser/injection counts
└── HGNC:<number> or SYMBOL:<symbol>
    ├── symbol: string
    ├── entrez_id: string, optional
    ├── ensembl_id: string, optional
    ├── mechanisms: [normalized mechanism string]
    ├── gene_level: [source assertion], optional [1..*]
    │   ├── {G2P_DDG2P
    │   │   ├── mechanism: normalized mechanism string, required
    │   │   ├── mechanism_raw: original G2P text, optional
    │   │   ├── disease: condition label, optional
    │   │   ├── inheritance: allelic requirement, optional
    │   │   ├── confidence: normalized confidence, optional
    │   │   └── pmids: [string], optional [1..*]
    │   ├── {PanelApp
    │   │   ├── mechanism: "PANELAPP_GREEN_NON_LOF_PATHO_HISTORY", required
    │   │   ├── confidence: "high", required
    │   │   ├── panel: panel name, optional
    │   │   ├── panel_id: panel identifier, optional
    │   │   └── inheritance: submitted mode of inheritance, optional
    │   └── {ClinGen_Dosage
    │       ├── mechanism: "TRIPLOSENSITIVITY", required
    │       ├── score: "2" | "3", required
    │       ├── score_description: source evidence description, optional
    │       └── pmids: [string], optional [1..*]
    └── variant_level: [source assertion], optional [1..*]
        ├── {GoFCards
        │   ├── mechanism: "GOF", required
        │   ├── source_record_id: public workbook row identifier, required
        │   ├── allele_key: stable raw hg19 source-allele key, required
        │   ├── exact_normalization_status, required
        │   │   ├── matched_gene_concordant
        │   │   ├── gene_discordant_coordinate_collision
        │   │   └── unmatched_public_source_allele
        │   ├── exact_cache_gene_symbols[], required only for a gene collision
        │   ├── exact_normalized_variants[], required for a matched assertion
        │   │   └── full 45-key normalized record; expanded below
        │   ├── disease: disorder label, optional
        │   ├── chr, pos, ref, alt: hg19 allele fields, optional individually
        │   ├── transcript: transcript accession, optional
        │   ├── consequence: missense_variant | indel, optional
        │   ├── pscore: number, optional
        │   ├── pmids: [string], optional [1..*]
        │   ├── function: curated functional text, optional
        │   ├── pathway: pathway text, optional
        │   ├── animal_model: true, optional
        │   └── cell_model: true, optional
        └── {ClinVar_VCV
            ├── mechanism_context: ["GOF"]
            ├── match
            │   ├── method: "exact_normalized_vcf_allele"
            │   ├── gene_concordance: boolean
            │   ├── clinvar_gene_symbols[]
            │   ├── exact_cache_gene_symbols[]
            │   ├── discarded_gene_discordant_symbols[]
            │   └── matched_gofcards_records[]
            │       ├── HGNC_Symbol
            │       ├── gofcards_accession_id, optional
            │       ├── gofcards_variant_id
            │       ├── allele_key
            │       ├── hg19_vcf_key
            │       ├── hg38_vcf_key
            │       ├── disease, optional
            │       ├── pmids, optional
            │       └── pscore, optional
            ├── variation
            │   ├── variation_id
            │   ├── vcv_accession
            │   ├── vcv_version
            │   ├── variation_name
            │   ├── variation_type
            │   ├── record_type: "classified"
            │   ├── record_status: "current"
            │   ├── date_created
            │   ├── date_last_updated, optional
            │   ├── classification_scope
            │   │   └── simple_allele | haplotype | genotype
            │   ├── matched_component_context
            │   │   └── allele | haplotype_component |
            │   │       genotype_component | genotype_haplotype_component
            │   ├── component_allele_id, optional
            │   ├── component_variation_id
            │   ├── component_name
            │   ├── component_variant_type
            │   ├── canonical_spdi, optional
            │   ├── genes[]
            │   │   └── symbol, hgnc_id, entrez_id, relationship_type
            │   ├── locations
            │   │   ├── GRCh37, optional
            │   │   │   └── assembly accession, chromosome, VCF pos/ref/alt/key
            │   │   └── GRCh38, optional
            │   │       └── assembly accession, chromosome, VCF pos/ref/alt/key
            │   └── hgvs[]
            │       └── type, expression, sequence_accession_version
            ├── condition_assertions[]
            │   ├── rcv_accession
            │   ├── rcv_version
            │   ├── title
            │   ├── trait_set_id
            │   ├── conditions[]
            │   │   └── name, database, id
            │   ├── germline_classification
            │   │   ├── clinical_significance
            │   │   ├── review_status
            │   │   ├── review_stars: 2 | 3 | 4
            │   │   ├── date_last_evaluated, optional
            │   │   └── submission_count, optional
            │   ├── matched_scvs[]
            │   │   ├── clinical_assertion_id
            │   │   ├── scv_accession and scv_version
            │   │   ├── record_status: "current"
            │   │   ├── contributes_to_aggregate_classification: true
            │   │   ├── submitter
            │   │   ├── classification
            │   │   ├── submitted_mode_of_inheritance[]
            │   │   ├── penetrance[]
            │   │   ├── trait_linkage_method
            │   │   ├── trait_linkage_ambiguous_across_eligible_rcvs: true,
            │   │   │   optional
            │   │   ├── linked_eligible_rcv_accessions[], optional
            │   │   ├── trait_mappings[]
            │   │   ├── observed_zygosity_counts
            │   │   └── observations[]
            │   │       ├── origin
            │   │       ├── affected_status
            │   │       ├── sex, optional
            │   │       ├── sample/chromosome counts, optional
            │   │       ├── observed_data[]
            │   │       └── cooccurrence_sets[]
            │   │           ├── asserted_variant_zygosity
            │   │           ├── count
            │   │           └── cooccurring_alleles[]
            │   │               └── name, relative_orientation,
            │   │                   zygosity, clinical_significance
            │   └── scv_linkage
            │       └── method and matched/unlinked counts
            └── allelic_requirement
                ├── value: "unresolved"
                └── reason: observation is not causal inheritance proof
```

### Required And Optional ClinVar Fields

| Object | Required keys | Optional keys | Empty collection allowed? |
|---|---|---|---|
| `ClinVar_VCV` | `mechanism_context`, `match`, `variation`, `condition_assertions`, `allelic_requirement` | none | `condition_assertions` must contain at least one RCV. |
| `match` | `method`, `gene_concordance`, `matched_gofcards_records` | `clinvar_gene_symbols`, `exact_cache_gene_symbols`, `discarded_gene_discordant_symbols` | The matched-record array is non-empty. The three symbol arrays may be empty; the builder currently emits them even when empty. |
| matched GoFCards record inside `match` | `HGNC_Symbol`, `gofcards_variant_id`, `allele_key`, and at least one of `hg19_vcf_key` or `hg38_vcf_key` | `gofcards_accession_id`, the other assembly key, `disease`, `pmids`, `pscore` | Not applicable. Optional scalar keys are omitted, not empty. |
| `variation` | `variation_id`, `vcv_accession`, `record_type`, `record_status`, `classification_scope`, `matched_component_context`, `locations` | `vcv_version`, `variation_name`, `variation_type`, creation/update dates, component identifiers/name/type, `canonical_spdi`, `genes`, `hgvs` | `locations` has at least one assembly. Optional `genes` and `hgvs` are omitted when empty. |
| `locations` | at least one of `GRCh37` or `GRCh38` | the other assembly | No: the object must contain at least one usable normalized VCF location. |
| assembly location | `chromosome`, `position_vcf`, `reference_allele_vcf`, `alternate_allele_vcf`, `vcf_key` | `assembly_accession_version` | Not applicable. |
| RCV `condition_assertion` | `rcv_accession`, `conditions`, `germline_classification`, `matched_scvs`, `scv_linkage` | `rcv_version`, `title`, `trait_set_id` | `conditions` and `matched_scvs` may be empty. An empty `matched_scvs` is a valid, explicit zero-link result. |
| RCV `germline_classification` | `review_status`, `review_stars` | `clinical_significance`, `date_last_evaluated`, `submission_count` | Not applicable. Review stars are always 2, 3, or 4 here. |
| SCV inside `matched_scvs` | `scv_accession`, `record_status`, `contributes_to_aggregate_classification`, `classification`, `trait_linkage_method` | `clinical_assertion_id`, `scv_version`, `submitter`, submitted MOI, penetrance, ambiguity fields, trait mappings, zygosity counts, observations | Optional SCV arrays are omitted when empty. `linked_eligible_rcv_accessions`, when present, has at least two RCVs. |
| SCV `classification` | object is required | `germline_significance`, SCV-level `review_status`, `date_last_evaluated` | Missing scalar members are omitted. The RCV review-star filter is not applied to this SCV-level review status. |
| `observation` | none | origin, affected status, sex, tested counts, `observed_data`, `cooccurrence_sets` | The observation is retained only when ClinVar supplied the branch; empty child arrays are omitted. |
| `cooccurrence_set` | none | asserted-variant zygosity, count, co-occurring alleles | Missing zygosity or phase is not inferred. |
| `allelic_requirement` | `value`, `reason` | none | `value` is always `"unresolved"`. |

The word **optional** therefore means “the key may not exist.” It does not
mean “the key exists with `null`, an empty string, or `NA`.” Required arrays
with `[0..*]` cardinality are the deliberate exception: the key exists and an
empty array communicates a real zero-result state.

## Complete Non-ClinVar Source Assertion Structures

These are evidence objects under a gene's `gene_level` or `variant_level`
array. They are distinct from `_meta.sources`, which only describes source
provenance and parser policy. Each array item contains exactly one source key.

### G2P / DDG2P

```text
{ "G2P_DDG2P": {
    "mechanism": string,          required
    "mechanism_raw": string,      optional
    "disease": string,            optional
    "inheritance": string,        optional
    "confidence": string,         optional
    "pmids": [string, ...]        optional, non-empty when present
} }
```

| Field | Emitted meaning and possible values |
|---|---|
| `mechanism` | Normalized, pipe-delimited mechanism label(s). Retained non-LOF labels can include `GOF`, `ACTIVATING`, `DOMINANT_NEGATIVE`, `TRIPLOSENSITIVITY`, or `INCREASED_DOSAGE`. The current production G2P assertions contain `GOF` or `DOMINANT_NEGATIVE`. |
| `mechanism_raw` | Original G2P molecular-mechanism text, without reinterpretation. |
| `disease` | G2P disease/condition label. |
| `inheritance` | Original G2P allelic-requirement text. It is not inferred from the mechanism. |
| `confidence` | Normalized `high`, `moderate`, `limited`, or `conflicting_or_refuted`. |
| `pmids` | Digit-only publication identifiers parsed from the source. |

The canonical cache excludes G2P rows whose normalized mechanism has no
curated non-LOF prompt-exception label. Thus a purely `LOF` G2P row can exist
in the raw download but does not appear in this canonical non-LOF JSON.

### PanelApp

```text
{ "PanelApp": {
    "mechanism": "PANELAPP_GREEN_NON_LOF_PATHO_HISTORY", required
    "confidence": "high",                              required
    "panel": string,                                   optional
    "panel_id": string,                                optional
    "inheritance": string                              optional
} }
```

Only PanelApp entries with source `confidence_level = 3` and a mode-of-
pathogenicity statement saying that loss-of-function variants do **not** cause
the phenotype are retained. This is a high-confidence non-LOF history marker;
it is not converted into a specific GOF or dominant-negative assertion.

### ClinGen Dosage

```text
{ "ClinGen_Dosage": {
    "mechanism": "TRIPLOSENSITIVITY", required
    "score": "2" | "3",              required
    "score_description": string,      optional
    "pmids": [string, ...]            optional, non-empty when present
} }
```

Score `2` means emerging/some evidence for dosage pathogenicity; score `3`
means sufficient evidence. The raw ClinGen file also contains
haploinsufficiency and non-positive dosage scores, but this canonical JSON is
a non-LOF mechanism cache, so only positive triplosensitivity assertions
survive its source filter.

### GoFCards

```text
{ "GoFCards": {
    "mechanism": "GOF",                    required
    "source_record_id": string,             required
    "allele_key": string,                   required
    "exact_normalization_status":           required; one of
        "matched_gene_concordant" |
        "gene_discordant_coordinate_collision" |
        "unmatched_public_source_allele",
    "exact_cache_gene_symbols": [string],   required only for gene collision
    "exact_normalized_variants": [          required only when matched [1..*]
        { full normalized record }
    ],
    "disease": string,                      optional
    "chr": string,                          optional, raw hg19
    "pos": string,                          optional, raw hg19
    "ref": string,                          optional, raw hg19; "-" may mean empty allele
    "alt": string,                          optional, raw hg19; "-" may mean empty allele
    "transcript": string,                   optional, raw GoFCards value
    "consequence": "missense_variant" |
                   "indel",                 optional, coarse source shape
    "pscore": number,                       optional
    "pmids": [string, ...],                 optional, non-empty when present
    "function": string,                     optional
    "pathway": string,                      optional
    "animal_model": true,                   optional
    "cell_model": true                      optional
} }
```

The public assertion is **not HGVS-only**. `allele_key` is derived from the
public hg19 chromosome/start/end/ref/alt tuple. The metadata object
`raw_public_allele_fields` scopes hg19 only to the legacy top-level
`chr`/`pos`/`ref`/`alt` source fields; it does not assign one assembly to the
assertion or to its normalized records. For a gene-concordant cache
join, `exact_normalized_variants` carries both source alleles and normalized
VCF alleles for hg19 and hg38, plus VEP-calibrated HGVSc/HGVSp when available.
All transcript/assembly rows are retained; the serializer does not select one
arbitrary transcript.

The three status branches prevent unsafe coordinate-only enrichment:

| `exact_normalization_status` | Emitted companion data | Meaning |
|---|---|---|
| `matched_gene_concordant` | `exact_normalized_variants[1..*]` | The stable source allele and resolved HGNC gene both agree. |
| `gene_discordant_coordinate_collision` | `exact_cache_gene_symbols[1..*]` | The allele key exists in the normalized cache but only under another gene. No normalized record is attached. |
| `unmatched_public_source_allele` | neither companion array | The public workbook allele is absent from the backend-derived normalized cache. No hg38/HGVS value is invented. |

The 2026-07-23 production-source audit contains 3,161 public assertions:
3,079 `matched_gene_concordant`, 28
`gene_discordant_coordinate_collision`, and 54
`unmatched_public_source_allele`. The 3,079 matched assertions emit 6,108
gene-concordant transcript/assembly records, and every matched assertion has
both an `hg19_vcf_key` and an `hg38_vcf_key`. The raw audit log is:

`/paedyl01/disk1/yangyxt/PriVA/data/gene_pathogenic_mechanism/audits/audit_gofcards_expanded_source_build_20260723.log`

Each object in `exact_normalized_variants` can emit the following complete
45-key structure. Optional means the key is absent when its normalized-cache
cell is blank or a placeholder; it never means JSON `null` or an empty string.

```text
exact_normalized_variants[]
├── source: "GoFCards"                                      required
├── mechanism: "GOF"                                       required
├── build: string                                           required
├── HGNC_Symbol: string                                     required
├── VEP_assembly: "hg19" | "hg38"                          optional
├── VEP_transcript: string                                  optional
├── feature_type: string                                    optional
├── consequence: string                                     optional
├── HGVSc: string                                           optional
├── HGVSp: string                                           optional
├── hgvsp_key: normalized protein-change key                optional
├── match_status: string                                    required
├── raw_GoFCards_HGVS: string                               optional
├── GoFCards_transcript: string                             optional
├── canonical_transcript: string                            optional
├── hg19_chrom: string                                      required
├── hg19_pos: string                                        required
├── hg19_ref: string                                        required
├── hg19_alt: string                                        required
├── hg19_vcf_status: string                                 required
├── hg38_chrom: string                                      required
├── hg38_pos: string                                        required
├── hg38_ref: string                                        required
├── hg38_alt: string                                        required
├── hg38_refalt_status: string                              required
├── gofcards_accession_id: dbSNP/other accession            optional
├── gofcards_variant_id: stable source allele identifier    required
├── disease: string                                         optional
├── pmids: [string, ...]                                    optional [1..*]
├── pscore: number | source string                          optional
├── function: string                                        optional
├── pathway: string                                         optional
├── derived_on: date string                                 required
├── allele_key: stable source allele key                    required
├── hg19_genomic_key: chrom|source-pos|source-ref|source-alt required
├── hg19_vcf_pos: string                                    required
├── hg19_vcf_ref: string                                    required
├── hg19_vcf_alt: string                                    required
├── hg38_vcf_pos: string                                    required
├── hg38_vcf_ref: string                                    required
├── hg38_vcf_alt: string                                    required
├── hg19_vcf_key: chrom|normalized-pos|ref|alt               required
├── hg38_genomic_key: chrom|source-pos|source-ref|source-alt required
├── hg38_vcf_key: chrom|normalized-pos|ref|alt               required
└── match_key_types: [hgvsp, hg19_genomic, hg19_vcf,
                      hg38_genomic, hg38_vcf]                required [1..*]
```

The difference between the genomic fields is deliberate:

- `hg19_pos/ref/alt` and `hg38_pos/ref/alt` are assembly-specific source-style
  alleles retained by the normalizer.
- `hg19_vcf_pos/ref/alt` and `hg38_vcf_pos/ref/alt` are padded and normalized
  VCF alleles used for caller/ClinVar-style exact comparison.
- `*_genomic_key` preserves the source-style tuple; `*_vcf_key` preserves the
  exact normalized VCF tuple.
- `HGVSc`, `HGVSp`, and `hgvsp_key` are absent for the 83 current rows without
  a concordant parseable protein annotation. Those rows remain usable through
  exact genomic matching.

`animal_model` and `cell_model` are presence-only positive flags: when the raw
value is `Y`, the key is emitted as `true`; otherwise the key is absent. An
absent key must not be interpreted as experimentally proven `false`.
`consequence` is the cache parser's coarse allele-shape classification: a
single-base substitution becomes `missense_variant`; other ref/alt shapes
become `indel`. Exact ClinVar matching does not use that label or HGVS. It uses
only `hg19_vcf_key` and `hg38_vcf_key`. PriVA additionally supports normalized
gene plus exact HGVSp and normalized gene plus exact hg19/hg38 genomic allele.

`ClinVar_VCV.match.matched_gofcards_records` remains a deliberately compact
9-key projection for the ClinVar join (`HGNC_Symbol`, optional
`gofcards_accession_id`, `gofcards_variant_id`, `allele_key`,
`hg19_vcf_key`, `hg38_vcf_key`, and optional `disease`, `pmids`, `pscore`). It
is not the full standalone GoFCards assertion shown above.

## ClinVar VCV XML Structure Used By The Parser

The authoritative structural metadata is the NCBI XSD, not an inferred table:

- `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/ClinVar_VCV_2.6.xsd`
- `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/_README`
- `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/`

The release root also embeds
`xsi:noNamespaceSchemaLocation="http://ftp.ncbi.nlm.nih.gov/pub/clinvar/xsd_public/ClinVar_VCV_2.6.xsd"`.
The downloaded weekly file reports `ReleaseDate="2026-07-20"`.

NCBI retired the old single-classification VCV XML format on 2025-08-07. This
parser targets the current split germline/somatic/oncogenic format only and
discovers the newest public `ClinVar_VCV_<version>.xsd` from the XSD index.

The PriVA installer stores the exact XSD and README beside the release under:

`/paedyl01/disk1/yangyxt/PriVA/data/gene_pathogenic_mechanism/raw/clinvar_vcv/format/`

Parser-relevant XSD branches are:

```text
ClinVarVariationRelease @ReleaseDate
└── VariationArchive*                                  VariationArchiveType
    ├── @VariationID, @Accession, @Version
    ├── @VariationName, @VariationType, @RecordType
    ├── @DateCreated, @DateLastUpdated
    ├── RecordStatus
    │   └── current | previous | replaced | deleted
    ├── ReplacedBy? / ReplacedList? / Comment? / Species
    └── choice
        ├── IncludedRecord                             skipped: no direct assertion
        │   ├── choice SimpleAllele | Haplotype
        │   ├── Classifications
        │   ├── SubmittedClassificationList/SCV+
        │   └── ClassifiedVariationList/ClassifiedVariation+
        └── ClassifiedRecord
            ├── choice: aggregate classification scope
            │   ├── SimpleAllele                       allele scope
            │   ├── Haplotype                          haplotype scope
            │   │   └── SimpleAllele+                  matched component(s)
            │   └── Genotype                           genotype scope
            │       ├── SimpleAllele+                  direct component(s)
            │       └── Haplotype+
            │           └── SimpleAllele+              nested component(s)
            ├── RCVList
            │   └── RCVAccession+
            │       ├── @Accession, @Version, @Title
            │       ├── ClassifiedConditionList
            │       │   └── ClassifiedCondition* @DB @ID
            │       └── RCVClassifications
            │           ├── GermlineClassification?    only this branch is eligible
            │           │   ├── ReviewStatus
            │           │   └── Description
            │           ├── SomaticClinicalImpact?     not used
            │           ├── OncogenicityClassification? not used
            │           └── NoClassification?          not used
            ├── Classifications                        VCV-wide aggregates; not filter
            │   ├── GermlineClassification?
            │   ├── SomaticClinicalImpactList?
            │   └── OncogenicityClassification?
            ├── ClinicalAssertionList
            │   └── ClinicalAssertion+                 SCV submission
            │       ├── @ID, @ContributesToAggregateClassification
            │       ├── ClinVarSubmissionID
            │       ├── ClinVarAccession @Accession=SCV... @Version
            │       ├── AdditionalSubmitters?
            │       ├── RecordStatus
            │       │   └── current | replaced | removed
            │       ├── Replaces* | ReplacedList?
            │       ├── Classification
            │       │   ├── ReviewStatus
            │       │   └── choice
            │       │       ├── GermlineClassification
            │       │       ├── SomaticClinicalImpact
            │       │       ├── OncogenicityClassification
            │       │       └── NoClassification
            │       ├── Assertion
            │       ├── AttributeSet*
            │       │   └── Attribute @Type
            │       │       ├── ModeOfInheritance
            │       │       ├── Penetrance
            │       │       ├── AgeOfOnset / Severity
            │       │       └── AssertionMethod / history
            │       ├── ObservedInList?
            │       │   └── ObservedIn+
            │       │       ├── Sample
            │       │       │   ├── Origin
            │       │       │   ├── AffectedStatus
            │       │       │   ├── Sex / age / family / counts?
            │       │       │   └── tissue/somatic fields?
            │       │       ├── Method*
            │       │       ├── choice ObservedData+ | FunctionalData
            │       │       ├── Co-occurrenceSet*
            │       │       │   ├── Zygosity?
            │       │       │   │   └── Homozygote | SingleHeterozygote |
            │       │       │   │       CompoundHeterozygote | Hemizygote |
            │       │       │   │       not provided
            │       │       │   ├── AlleleDescSet*
            │       │       │   │   ├── Name
            │       │       │   │   ├── RelativeOrientation?: cis | trans | unknown
            │       │       │   │   ├── Zygosity?
            │       │       │   │   └── ClinicalSignificance?
            │       │       │   └── Count?
            │       │       └── TraitSet? / citations / xrefs / comments
            │       ├── choice SimpleAllele | Haplotype | Genotype
            │       ├── TraitSet                       submitted condition
            │       └── citations/study/comments?
            └── TraitMappingList?
                └── TraitMapping+
                    ├── @ClinicalAssertionID           links mapping to SCV
                    ├── @TraitType, @MappingType, @MappingValue, @MappingRef
                    └── MedGen @Name @CUI
```

### SimpleAllele Subtree

Every aggregate `SimpleAllele` may contain:

```text
SimpleAllele @AlleleID @VariationID
├── GeneList?/Gene*                                   zero, one, or many genes
├── Name
├── CommonName?
├── CanonicalSPDI?
├── VariantType?
├── Location?
│   ├── CytogeneticLocation*
│   └── SequenceLocation*                             multiple assemblies
│       └── @Assembly @Chr @positionVCF
│           @referenceAlleleVCF @alternateAlleleVCF
├── OtherNameList?
├── ProteinChange*
├── HGVSlist?/HGVS*
│   └── NucleotideExpression and/or ProteinExpression
├── Classifications?
├── XRefList?
├── Comment*
├── AlleleFrequencyList?
└── GlobalMinorAlleleFrequency?
```

Only the normalized VCF tuple is used for exact matching. HGVS, SPDI, genes,
and names are retained as provenance but are not substitute match routes.

## XML Edge Decisions

| XML edge | Parser decision |
|---|---|
| `IncludedRecord` | Skip. It has no direct submitted classification. |
| non-current VCV | Skip `previous`, `replaced`, and `deleted`. |
| simple allele | RCV classification is directly about the matched allele. |
| haplotype/genotype | Retain with aggregate `classification_scope` and component-match label; never present the aggregate classification as allele-only. |
| several matching components | Emit one component-aware match per matching component. |
| several genes | Prefer intersection of ClinVar gene symbols and exact-cache symbols; preserve gene-concordance flag. |
| ClinVar supplies gene symbols but none match GoFCards | Reject the component as a coordinate collision. Fall back to the GoFCards symbol only when ClinVar supplies no gene symbol. |
| missing GRCh37 or GRCh38 | Match using whichever normalized assembly is available. |
| RCV has only somatic/oncogenic/no-classification branch | Exclude. |
| RCV has 0-1 star germline status | Exclude only that condition, not necessarily the whole VCV. |
| eligible RCV with no linkable SCV | Keep the RCV aggregate; emit an empty `matched_scvs` and linkage counts. |
| replaced/removed or non-contributing SCV | Do not use its observations. |
| SCV condition maps to another RCV | Do not attach it to the eligible RCV. |
| one SCV maps to multiple eligible RCVs | Retain it under each matching RCV and mark `trait_linkage_ambiguous_across_eligible_rcvs: true` with all `linked_eligible_rcv_accessions`; do not choose one condition. |
| no zygosity | Keep the SCV observation with no invented zygosity. |
| zygosity present | Store the observation/count exactly; do not infer dominance/recessivity. |
| cis/trans co-occurrence | Store relative orientation and co-occurring allele name. |
| affected/unaffected/unknown | Store the submitted status; do not treat an observation as segregation proof by itself. |
| submitted MOI/penetrance | Store as submitted evidence; leave causal allelic requirement unresolved. |

## Weekly Release And Structural Metadata

Release URL:

`https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/weekly_release/ClinVarVCVRelease_00-latest_weekly.xml.gz`

Companion checksum:

`https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/weekly_release/ClinVarVCVRelease_00-latest_weekly.xml.gz.md5`

The updater checks the small MD5 once per seven days. It downloads the XML only
when the remote MD5 differs, resumes an interrupted transfer into an
MD5-specific `.part` file, verifies the completed MD5, and atomically replaces
the previous valid release. The exact XSD version, XSD checksum, README
checksum, XML MD5, XML SHA256, file size, and timestamps are stored in:

`/paedyl01/disk1/yangyxt/PriVA/data/gene_pathogenic_mechanism/metadata/nonlof_source_manifest.json`

The XSD is the structural-format metadata file requested for this integration.
The README describes release naming/cadence; the companion MD5 describes file
integrity, not XML structure.
