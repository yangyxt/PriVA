#!/usr/bin/env python3
r"""Gene-condition and query-variant mechanism hub for PriVA.

THE LOGIC CHAIN
===============

This module exists to serve one chain, and nothing in it that fails to serve
that chain is justified. The purpose is to raise the resolution at which
inheritance and penetrance are known, from the gene to the individual
variant-transcript-consequence. A gene's inheritance history on its own does
not say which condition produced it, nor what pathogenic mechanism that
condition acts by, so it cannot be attached to the variant in front of us.
This chain can.

This filename remains PriVA's pipeline-facing entrance and command-line
orchestrator. The gene_mechanism_common, gene_mechanism_conditions,
gene_mechanism_variants, gene_mechanism_resources and
gene_mechanism_annotation modules hold definitions only. Pipeline callers
must import this facade, and none of those definition modules may import it.

    the variant                          the gene's condition histories
    (one row, one transcript)            (from the HPO condition cache)
          |                                          |
          |  STEP 1                                  |  STEP 2
          |  score LOF / GOF / DN                    |  germline-included
          |  each as 0, 1 or 2                       |  conditions only
          v                                          v
      +---------------------+              +----------------------------+
      | variant_lof_score   |              | condition                  |
      | variant_gof_score   |              | inheritance                |
      | variant_dn_score    |              | penetrance                 |
      | variant_effect      |              | mechanism (LOF/GOF/DN)     |
      +---------------------+              +----------------------------+
                    \                              /
                     \        STEP 3              /
                      \   keep the histories     /
                       \  this mechanism reaches
                        v                        v
                   +--------------------------------------+
                   | condition + inheritance + penetrance |
                   |   ...and nothing else travels        |
                   +--------------------------------------+
                                    |
                                    |  STEP 4  attach to the row
                                    v
                       every ACMG criterion now reads an
                       inheritance and a penetrance that
                       belong to THIS variant's mechanism

STEP 1  the variant's own mechanism, as three scores
====================================================

    2   exclusively established -- an exact curated allele match.
        When any mechanism scores 2 the other two are forced to 0.
    1   plausible -- either supported by prediction, or not excluded when no
        mechanism-specific variant evidence exists.
    0   excluded by a stronger, exclusive variant call.

The scale grades the QUERY VARIANT. It is not a statement about how well
attested a condition's mechanism is; that is a separate thing, recorded in
the cache as assertion_basis.

    match_curated_nonlof_variant   supplies the exact score (2) from the
                                   non-LOF cache, by HGVSp then HGVSc then
                                   normalized genomic allele
    infer_query_variant_effect     combines it with PriVA's own annotation

    variant_effect, one of:
        "exact_known_GOF"                  |  also _DOMINANT_NEGATIVE, and
        "exact_known_DOMINANT_NEGATIVE"    |  "+"-joined when both apply
        "predicted_LOF_high_confidence"
        "uncertain"

    Predicted loss of function is GRADED, not flat:

        score 2, the transcript is destroyed
            NMD_PREDICTED_LOF   stop_gained / frameshift / splicing_frameshift
                                that does NOT escape nonsense-mediated decay
                                (NMD does not say "escaping", and LoF_filter
                                is not END_TRUNC)
            LOFTEE_HC           LoF == "HC". A high-confidence call rescues a
                                variant that does escape decay.

        score 1, loss of function is plausible, the transcript may survive
            NMD_ESCAPING_TRUNCATION  truncating or splice-frameshift, but
                                     escaping decay, with no LOFTEE rescue
            LOFTEE_OS           LoF == "OS", LOFTEE's other-splice category
            VEP_LOF             vep_consq_lof is true
            PRIVA_SPLICE_LOF    splicing_lof is true, frame not shifted
            PRIVA_5UTR_LOF      5UTR_lof is true

        With no mechanism-specific evidence, LOF, GOF and DN all score 1. The
        three parallel hypotheses are the uncertainty representation; no
        variant-side UNRESOLVED mechanism is invented.

    THE PRECEDENCE, HIGHEST FIRST
        1. nonsense-mediated decay        effect = "exact_known_LOF"
           The transcript is destroyed, so there is no protein left to gain a
           function or to poison a complex. This is the one place a prediction
           outranks curation: a curated gain-of-function call for the same
           allele is moved into variant_effect_suppressed_evidence, not kept.
        2. a curated allele              effect = "exact_known_GOF" | "..._DN"
           Outranks every remaining loss-of-function signal, including a
           LOFTEE high-confidence rescue: those are predictions, this is a
           curator's verdict on this exact allele.
        3. other predicted loss           effect = "predicted_LOF_high_confidence"
           2 if LOFTEE rescues it, otherwise 1.
        4. nothing                        effect = "uncertain"; LOF=GOF=DN=1

    After an exact match these are kept in
    variant_effect_suppressed_evidence but cannot create a competing LOF
    call: an exact curated mechanism is exclusive.

STEP 2  the gene's condition histories, germline disease only
=============================================================

Read from the HPO condition cache built by
build_hpo_condition_mechanism_cache.py. The builder physically removes every
condition whose final merged scope decision is ``exclude`` and records a
per-gene exclusion count. Of the retained conditions, a history is admissible
only when

    condition.priva_scope.decision == "include"

enforced in condition_cache_context. Retained ``review`` or unscoped conditions
are audit only and can never influence a criterion.

STEP 3  match by mechanism, and carry three things back
=======================================================

select_condition_histories_for_variant filters histories whose curated or
deduced mechanism is incompatible with the variant. A variant with no LOF,
GOF, or dominant-negative evidence cannot prove any of those histories
incompatible, so all three remain plausible. A condition with no usable
mechanism assertion is retained as ``UNRESOLVED``, because absence of a
condition-mechanism assertion is not evidence that the condition is irrelevant:

    variant_effect                      histories kept
    ---------------------------------   ----------------------------
    exact_known_GOF                     GOF + UNRESOLVED
    exact_known_DOMINANT_NEGATIVE       DOMINANT_NEGATIVE + UNRESOLVED
    predicted_LOF_high_confidence       LOF + UNRESOLVED
    uncertain                           LOF + GOF + DOMINANT_NEGATIVE + UNRESOLVED

Each surviving history is reduced to exactly this, and nothing else travels:

    {
      "mechanism":   "LOF" | "GOF" | "DOMINANT_NEGATIVE" | "UNRESOLVED",
      "condition":   the OMIM / ORPHA / MONDO identifier,
      "inheritance": see the vocabulary below,
      "x_linked":    True | False,
      "penetrance":  "complete" | "incomplete" | "",
    }

THE INHERITANCE VOCABULARY, VALUE BY VALUE
==========================================

Three source vocabularies contribute different facts. G2P states an allelic
requirement, HPO states an inheritance mode, and ClinGen HI=3 establishes a
condition-linked loss-of-function dosage mechanism. ClinGen supplies an
autosomal monoallelic requirement. On X, condition-linked G2P/HPO inheritance
wins; an HI=3 condition with neither defaults to X-linked dominant. HI=3 on Y
remains Y-linked.

    delivered value   source values that fold to it
    ---------------   ---------------------------------------------------
    "recessive"       biallelic_autosomal          (G2P)
                      autosomal_recessive          (HPO)
                      pseudoautosomal_recessive    (HPO)
    "dominant"        monoallelic_autosomal        (G2P or autosomal HI=3)
                      autosomal_dominant           (HPO)
                      autosomal_dominant_maternal_imprinting  (HPO)
    "x_linked_recessive"
                      x_linked_recessive           (HPO)
    "x_linked_dominant"
                      monoallelic_X_heterozygous   (G2P)
                      x_linked_dominant            (HPO)
    "x_linked_unspecified"
                      monoallelic_X or
                      monoallelic_X_hemizygous     (G2P observation)
                      x_linked                     (HPO, unqualified)
    "y_linked"        monoallelic_Y_hemizygous (G2P), y_linked (HPO)
    "mitochondrial"   mitochondrial (both)
    "non_mendelian"   \
    "polygenic"        |  HPO only, and reported ONLY when no Mendelian
    "digenic"          |  mode accompanies them on the same condition
    "oligogenic"      /
    ""                nothing stated -- 17.3% of assertions

    x_linked is also returned as a compatibility flag. Hemizygosity is not an
    inheritance value; it remains an observation used locally by sex-aware
    criteria.

    Downstream reading of the values:
        DOMINANT_LIKE_INHERITANCE = {dominant,
                                     x_linked_dominant,
                                     y_linked, mitochondrial}
            one allele in the asserted carrier state is enough
        x_linked_recessive
            explicit X-linked recessive history; criteria use XY observations
        x_linked_unspecified
            no genotype-based benign shortcut is permitted
        NON_MENDELIAN_INHERITANCE = {non_mendelian, polygenic, digenic,
                                     oligogenic}
            germline but not single-gene. Delivered rather than discarded
            precisely so benign-supporting criteria are NOT assigned
            easily against them.

THE PENETRANCE VOCABULARY
=========================

    "incomplete"   literal incomplete/moderate/low penetrance, or late onset.
                   Other onset terms, age-dependent penetrance, and variable
                   expressivity remain source context but do not block an
                   apparently unaffected carrier observation.
    "complete"     explicit complete or high penetrance
    ""             nothing usable stated. No third ``unknown`` category is
                   emitted.

    Sex-limited expression (HP:0001470) and slowly progressive disease
    (HP:0003677) are deliberately excluded: the former is not penetrance and
    the latter describes change after onset rather than whether onset occurs.

WHEN A MECHANISM IS DEDUCED, UNKNOWN OR INCOMPATIBLE
===================================================

For a germline-included recessive condition without an explicit curated
mechanism, the cache supplies a high-likelihood LOF working inference with
``assertion_basis="deduced"``. The inference participates in
variant-to-condition matching but is not established LOF history: PVS1 accepts
only ``assertion_basis="curated"``. Other included conditions without a usable
mechanism emit an ``UNRESOLVED`` history, preserving their condition-linked
inheritance, penetrance, onset and provenance. ``variant_inheritance_basis``
makes the selected result auditable:

    matched_history                 compatible curated or deduced histories
    unresolved_condition_history   mechanism-free condition histories only
    matched_and_unresolved_history  both kinds contributed
    mechanism_mismatch              known histories existed, but all were
                                    incompatible with this variant
    condition_scope_blocked         retained review conditions or a build-time
                                    exclusion tombstone blocks fallback
    gene_constraint                 no condition record existed for the gene

The last fallback supplies inheritance only and never a mechanism or
penetrance. It is not used after a known mechanism mismatch and cannot revive
a retained review condition or a build-pruned exclusion.

STEP 4  what lands on the row
=============================

    var_plausible_patho_mechs        "<inheritance>_<MECHANISM>" tags,
                                     e.g. "dominant_GOF;recessive_LOF" or
                                     "x_linked_recessive_LOF"
                                     (DOMINANT_NEGATIVE abbreviates to DN)
    variant_effect                   see step 1
    variant_lof_score                0 | 1 | 2
    variant_gof_score                0 | 1 | 2
    variant_dn_score                 0 | 1 | 2
    variant_mechanism_exclusive      True when any score is 2
    variant_exact_mechanisms         ";"-joined mechanisms scoring 2
    variant_mechanism_applicable     tags whose mechanism is ESTABLISHED
                                     for this variant
    variant_mechanism_uncertain      tags whose mechanism is only POSSIBLE;
                                     includes all matched histories when the
                                     variant has no mechanism-specific evidence
    variant_condition_ids            ";"-joined condition identifiers
    variant_condition_histories      the same facts kept TOGETHER, one
                                     entry per history:

        <condition>|<mechanism>|<inheritance>|<penetrance>

        e.g. for a truncating variant in ABCB4, a gene with both a dominant
        and a recessive condition:

        OMIM:600803|LOF|dominant|;OMIM:602347|LOF|recessive|

        and for ATP1A2, where a penetrance is actually recorded:

        OMIM:602481|LOF|dominant|incomplete;OMIM:619602|LOF|recessive|

        This column exists because the three flat lists below are each
        de-duplicated separately, so they can have different lengths and a
        reader cannot tell which inheritance belongs to which condition.
        Empty when no history was reached; the basis column then says which
        fallback supplied the inheritance.

    variant_inheritance             ";"-joined inheritance values
    variant_inheritance_basis       one of the bases documented above
    variant_x_linked                "true" | "false"
    variant_penetrance              ";"-joined complete | incomplete, or empty.
                                    Stated for
                                    only 1.0% of assertions, because HPO
                                    rarely annotates penetrance at all.

WHAT IS OBSOLETE
================

Anything that does not serve the four steps. That explicitly includes
re-reading ClinVar per variant to discover a condition: the condition comes
from matching a mechanism to a history, never from a second pass over
ClinVar. Gene-wide signals are likewise absent by design -- a ClinVar
pathogenic history, a constrained LOEUF, a high average AlphaMissense
score, a ClinGen dosage call at the gene level -- because none of them says
which condition a variant acts on, or by what mechanism.
"""


from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from gene_mechanism_annotation import (
    HISTORY_IDENTITY_FIELDS,
    annotate_gene_mechanism_categories,
    select_condition_histories_for_variant,
)
from gene_mechanism_common import (
    CANONICAL_MECHANISMS,
    CONDITION_MECHANISM_EVIDENCE_COLUMNS,
    CONDITION_MECHANISM_SOURCES,
    CURATED_EXACT_MECHANISMS,
    DATA_DIR,
    DDG2P_DOMINANT_LOF_INHERITANCE,
    DDG2P_LOF_RAW_VALUES,
    DDG2P_RECESSIVE_LOF_INHERITANCE,
    DDG2P_USABLE_DISEASE_CONFIDENCE,
    DDG2P_USABLE_MECHANISM_CONFIDENCE,
    DEFAULT_CLINGEN_DOSAGE,
    DEFAULT_DDG2P_MECHANISM_EVIDENCE,
    DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
    DEFAULT_HGNC_TABLE,
    DEFAULT_HPO_ASSERTIONS,
    DEFAULT_HPO_COLLAPSED,
    DEFAULT_HPO_CONDITION_MECHANISM_CACHE,
    DEFAULT_LOEUF_TABLE,
    DEFAULT_MECHANISM_JSON,
    EXACT_SEQUENCE_MECHANISMS,
    HPO_CACHE_INHERITANCE_LABELS,
    HPO_CONDITION_MECHANISM_SCHEMA_VERSION,
    HPO_INHERITANCE_TERMS,
    HPO_ONSET_TERMS,
    HPO_PENETRANCE_TERMS,
    LOOKUP_FIELD_PRIORITY,
    PRIVA_ROOT,
    SCRIPT_DIR,
    UNRESOLVED_MECHANISM,
    VARIANT_MECHANISM_OUTPUT_COLUMNS,
    VARIANT_MECHANISM_SCORE_KEYS,
)
from gene_mechanism_conditions import (
    ASSERTION_IDENTITY_FIELDS,
    DOMINANT_LIKE_INHERITANCE,
    NON_MENDELIAN_INHERITANCE,
    build_hpo_condition_index,
    condition_cache_context,
    condition_cache_mechanism_assertions,
    condition_cache_mechanism_entries,
    enrich_condition_mechanism_assertion,
    gene_has_lof_mechanism_history,
    gene_inheritance_consensus,
    gene_inheritance_from_constraint,
    normalize_inheritance,
    normalize_inheritance_histories,
    normalize_penetrance,
)
from gene_mechanism_resources import GeneMechanismHub, _LocalHgncResolver
from gene_mechanism_variants import infer_query_variant_effect


__all__ = (
    "GeneMechanismHub",
    "annotate_gene_mechanism_categories",
    "build_hpo_condition_index",
    "condition_cache_context",
    "condition_cache_mechanism_assertions",
    "condition_cache_mechanism_entries",
    "enrich_condition_mechanism_assertion",
    "gene_has_lof_mechanism_history",
    "gene_inheritance_consensus",
    "gene_inheritance_from_constraint",
    "infer_query_variant_effect",
    "match_curated_nonlof_variant",
    "match_gofcards_variant_gof",
    "match_gofcards_variant_gof_by_genomic",
    "normalize_inheritance",
    "normalize_inheritance_histories",
    "normalize_penetrance",
    "resolve_gene_symbol",
    "select_condition_histories_for_variant",
    "DEFAULT_CLINGEN_DOSAGE",
    "DEFAULT_DDG2P_MECHANISM_EVIDENCE",
    "DEFAULT_GOFCARDS_EXACT_GOF_HGVSP",
    "DEFAULT_HGNC_TABLE",
    "DEFAULT_HPO_ASSERTIONS",
    "DEFAULT_HPO_COLLAPSED",
    "DEFAULT_HPO_CONDITION_MECHANISM_CACHE",
    "DEFAULT_LOEUF_TABLE",
    "DEFAULT_MECHANISM_JSON",
    "VARIANT_MECHANISM_OUTPUT_COLUMNS",
)


def resolve_gene_symbol(gene_symbol: Any) -> str:
    """Convenience function returning one current HGNC symbol."""
    return GeneMechanismHub().resolve_symbol(gene_symbol)




def match_curated_nonlof_variant(
    gene_symbol: Any,
    *,
    hgvsp: Any = "",
    chrom: Any = "",
    pos: Any = "",
    ref: Any = "",
    alt: Any = "",
    assembly: str = "hg38",
    key_type: str = "auto",
    mechanism_json: str | Path = DEFAULT_MECHANISM_JSON,
    use_hgnc_package: bool = False,
) -> dict[str, Any]:
    """Convenience API for canonical exact GoF/dominant-negative matching."""
    return GeneMechanismHub(
        mechanism_json=mechanism_json,
        use_hgnc_package=use_hgnc_package,
    ).match_curated_nonlof_variant(
        gene_symbol,
        hgvsp=hgvsp,
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        assembly=assembly,
        key_type=key_type,
    )


def match_gofcards_variant_gof(
    gene_symbol: Any,
    hgvsp: Any,
    *,
    gofcards_exact_hgvsp_tsv: str | Path = DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
    gofcards_step1_tsv: str | Path | None = None,
    gofcards_active_tsv: str | Path | None = None,
    gofcards_raw_xlsx: str | Path | None = None,
    gofcards_conversion_audit_tsv: str | Path | None = None,
    use_hgnc_package: bool = False,
) -> dict[str, Any]:
    """Convenience wrapper for exact GoFCards HGVSp variant-level GOF matching."""
    return GeneMechanismHub(use_hgnc_package=use_hgnc_package).match_gofcards_variant_gof(
        gene_symbol,
        hgvsp,
        gofcards_exact_hgvsp_tsv=gofcards_exact_hgvsp_tsv,
        gofcards_step1_tsv=gofcards_step1_tsv,
        gofcards_active_tsv=gofcards_active_tsv,
        gofcards_raw_xlsx=gofcards_raw_xlsx,
        gofcards_conversion_audit_tsv=gofcards_conversion_audit_tsv,
    )


def match_gofcards_variant_gof_by_genomic(
    gene_symbol: Any,
    chrom: Any,
    pos: Any,
    ref: Any,
    alt: Any,
    *,
    assembly: str = "hg38",
    key_type: str = "auto",
    gofcards_exact_hgvsp_tsv: str | Path = DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
    use_hgnc_package: bool = False,
) -> dict[str, Any]:
    """Convenience wrapper for exact GoFCards genomic variant-level GOF matching."""
    return GeneMechanismHub(use_hgnc_package=use_hgnc_package).match_gofcards_variant_gof_by_genomic(
        gene_symbol,
        chrom,
        pos,
        ref,
        alt,
        assembly=assembly,
        key_type=key_type,
        gofcards_exact_hgvsp_tsv=gofcards_exact_hgvsp_tsv,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gene", nargs="+", help="Gene symbol/alias/HGNC/Ensembl/Entrez query")
    parser.add_argument("--include-entries", action="store_true")
    parser.add_argument(
        "--condition-cache",
        default=str(DEFAULT_HPO_CONDITION_MECHANISM_CACHE),
        help="Integrated HPO condition-mechanism JSON cache",
    )
    parser.add_argument(
        "--mechanism-json",
        default=str(DEFAULT_MECHANISM_JSON),
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--ddg2p-evidence",
        default=str(DEFAULT_DDG2P_MECHANISM_EVIDENCE),
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--hpo-collapsed", default=str(DEFAULT_HPO_COLLAPSED))
    parser.add_argument("--clingen-dosage", default=str(DEFAULT_CLINGEN_DOSAGE))
    parser.add_argument("--loeuf-table", default=str(DEFAULT_LOEUF_TABLE))
    parser.add_argument("--hgnc-table", default=str(DEFAULT_HGNC_TABLE))
    args = parser.parse_args()

    hub = GeneMechanismHub(
        condition_cache=args.condition_cache,
        mechanism_json=args.mechanism_json,
        ddg2p_evidence=args.ddg2p_evidence,
        hpo_collapsed=args.hpo_collapsed,
        clingen_dosage=args.clingen_dosage,
        loeuf_table=args.loeuf_table,
        hgnc_table=args.hgnc_table,
    )
    result = {
        gene: hub.gene_summary(gene, include_entries=args.include_entries)
        for gene in args.gene
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
