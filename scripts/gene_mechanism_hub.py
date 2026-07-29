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
    1   plausible -- supported by consequence or prediction only.
    0   confidently not this mechanism.

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

        score 0 is reserved for no loss-of-function evidence at all.

    A loss-of-function score of 2 does NOT make the mechanism exclusive:
    EXACT_SEQUENCE_MECHANISMS is {GOF, DOMINANT_NEGATIVE}, so a curated
    gain-of-function or dominant-negative allele still overrides a predicted
    loss of function. Curation beats prediction.

    After an exact match these are kept in
    variant_effect_suppressed_evidence but cannot create a competing LOF
    call: an exact curated mechanism is exclusive.

STEP 2  the gene's condition histories, germline disease only
=============================================================

Read from the HPO condition cache built by
build_hpo_condition_mechanism_cache.py. A history is admissible only when

    condition.priva_scope.decision == "include"

enforced in condition_cache_context. Everything else -- "review",
"exclude", or unscoped -- is audit only and can never influence a criterion.

STEP 3  match by mechanism, and carry three things back
=======================================================

select_condition_histories_for_variant keeps only histories whose mechanism
the variant plausibly acts by:

    variant_effect                      histories kept
    ---------------------------------   ----------------------------
    exact_known_GOF                     GOF only
    exact_known_DOMINANT_NEGATIVE       DOMINANT_NEGATIVE only
    predicted_LOF_high_confidence       LOF only
    uncertain                           all three, marked uncertain

Each surviving history is reduced to exactly this, and nothing else travels:

    {
      "mechanism":   "LOF" | "GOF" | "DOMINANT_NEGATIVE",
      "condition":   the OMIM / ORPHA / MONDO identifier,
      "inheritance": see the vocabulary below,
      "x_linked":    True | False,
      "penetrance":  "complete" | "incomplete" | "unknown",
    }

THE INHERITANCE VOCABULARY, VALUE BY VALUE
==========================================

Two source vocabularies describe one fact. G2P states an allelic
requirement; HPO states an inheritance mode. normalize_inheritance folds
both, and accepts an HPO mode as either the cache's snake_case key or the
human-readable label that condition_cache_context substitutes.

    delivered value   source values that fold to it
    ---------------   ---------------------------------------------------
    "recessive"       biallelic_autosomal          (G2P)
                      monoallelic_X_hemizygous     (G2P)
                      monoallelic_X                (G2P)
                      autosomal_recessive          (HPO)
                      x_linked_recessive           (HPO)
                      pseudoautosomal_recessive    (HPO)
                      x_linked                     (HPO, bare -- read as
                                                    X-linked recessive)
    "dominant"        monoallelic_autosomal        (G2P)
                      monoallelic_X_heterozygous   (G2P)
                      autosomal_dominant           (HPO)
                      x_linked_dominant            (HPO)
                      autosomal_dominant_maternal_imprinting  (HPO)
    "y_linked"        monoallelic_Y_hemizygous (G2P), y_linked (HPO)
    "mitochondrial"   mitochondrial (both)
    "non_mendelian"   \
    "polygenic"        |  HPO only, and reported ONLY when no Mendelian
    "digenic"          |  mode accompanies them on the same condition
    "oligogenic"      /
    ""                nothing stated -- 17.3% of assertions

    x_linked is returned SEPARATELY, true when the requirement starts
    monoallelic_X or any mode starts x_linked. Folding it into the value
    would erase it, and a hemizygous male affected by one allele is still
    the recessive pattern.

    Downstream reading of the values:
        DOMINANT_LIKE_INHERITANCE = {dominant, y_linked, mitochondrial}
            one allele in a carrier is enough
        NON_MENDELIAN_INHERITANCE = {non_mendelian, polygenic, digenic,
                                     oligogenic}
            germline but not single-gene. Delivered rather than discarded
            precisely so benign-supporting criteria are NOT assigned
            easily against them.

THE PENETRANCE VOCABULARY
=========================

    "incomplete"   HP:0003829 incomplete, HP:4000159 moderate,
                   HP:4000160 low  -- all three are forms of the condition
                   failing to appear in some carriers, which is the only
                   distinction the criteria act on
    "complete"     HP:0034950
    "unknown"      nothing stated -- 98.9% of assertions, because HPO
                   rarely annotates penetrance

WHEN NO HISTORY STATES AN INHERITANCE: THE FALLBACK LADDER
==========================================================

A variant must not be left with an empty inheritance simply because the
gene's conditions carry no mechanism. Three tiers are tried in order, most
grounded first, and variant_inheritance_basis always says which one answered.

    basis = "matched_history"      3,897 genes   70.0%
        The variant's own mechanism reached at least one condition history,
        and that history states the inheritance. This is the resolution the
        whole chain exists to reach.

    basis = "gene_consensus"         879 genes   15.8%
        No history states one, but every germline-included condition of the
        gene agrees on a single inheritance, so whatever this variant does,
        that is what the disease requires. Unanimity is what makes it safe.
        Delivers 839 dominant, 23 mitochondrial, 17 y_linked. It delivers no
        recessive genes, because a recessive condition with no mechanism has
        already been given LOF by the cache build, so it resolves as a
        matched history instead.

    basis = "gene_constraint"        789 genes   14.2%
        HPO states no inheritance for this gene at all -- for 755 of these
        the annotation simply does not exist. The constraint data then
        decides how many copies must be affected:

            dominant    ClinGen haploinsufficiency score 3
                        OR LOEUF below 0.35
            recessive   otherwise, the default

        Either signal is enough: of these genes 196 carry only the LOEUF
        signal, 6 only the ClinGen one, and 9 both, so requiring both would
        reduce the rule to 9 genes. Delivers 580 recessive, 209 dominant.

    basis = ""                         0 genes
        Reserved for a gene whose germline conditions disagree with each
        other AND whose mechanism reached nothing. No such gene exists in
        the current data: the 501 genes with mixed inheritance all carry
        mechanisms, so they resolve as matched histories.

    No mechanism accompanies either fallback. For a gene answered by
    constraint we know how many copies must be affected, not what a variant
    has to do to the protein.

STEP 4  what lands on the row
=============================

    var_plausible_patho_mechs        "<inheritance>_<MECHANISM>" tags,
                                     e.g. "dominant_GOF;recessive_LOF"
                                     (DOMINANT_NEGATIVE abbreviates to DN)
    variant_effect                   see step 1
    variant_lof_score                0 | 1 | 2
    variant_gof_score                0 | 1 | 2
    variant_dn_score                 0 | 1 | 2
    variant_mechanism_exclusive      True when any score is 2
    variant_exact_mechanisms         ";"-joined mechanisms scoring 2
    variant_mechanism_applicable     tags whose mechanism is ESTABLISHED
                                     for this variant
    variant_mechanism_uncertain      tags whose mechanism is only POSSIBLE
                                     because the variant effect is
                                     unresolved
    variant_condition_ids            ";"-joined condition identifiers
    variant_condition_histories      the same facts kept TOGETHER, one
                                     entry per history:

        <condition>|<mechanism>|<inheritance>|<penetrance>

        e.g. for a truncating variant in ABCB4, a gene with both a dominant
        and a recessive condition:

        OMIM:600803|LOF|dominant|unknown;OMIM:602347|LOF|recessive|unknown

        and for ATP1A2, where a penetrance is actually recorded:

        OMIM:602481|LOF|dominant|incomplete;OMIM:619602|LOF|recessive|unknown

        This column exists because the three flat lists below are each
        de-duplicated separately, so they can have different lengths and a
        reader cannot tell which inheritance belongs to which condition.
        Empty when no history was reached; the basis column then says which
        fallback supplied the inheritance.

    variant_inheritance             ";"-joined inheritance values
    variant_inheritance_basis       matched_history | gene_consensus
                                    | gene_constraint | ""
    variant_x_linked                "true" | "false"
    variant_penetrance              ";"-joined, or "unknown". Stated for
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
import csv
import json
import logging
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from clinvar_vcv import open_text

SCRIPT_DIR = Path(__file__).resolve().parent
PRIVA_ROOT = SCRIPT_DIR.parent
DATA_DIR = PRIVA_ROOT / "data"

DEFAULT_HGNC_TABLE = DATA_DIR / "hgnc" / "non_alt_loci_set.tsv"
DEFAULT_HPO_ASSERTIONS = DATA_DIR / "hpo" / "genes_to_phenotype.assertions.tsv.gz"
# Keep the existing constructor argument name while migrating its file schema.
DEFAULT_HPO_COLLAPSED = DEFAULT_HPO_ASSERTIONS
DEFAULT_CLINGEN_DOSAGE = DATA_DIR / "clingen" / "gene_dosage_sensitivity.hg19.tsv"
DEFAULT_LOEUF_TABLE = DATA_DIR / "loeuf" / "loeuf_dataset.tsv.gz"
DEFAULT_HPO_CONDITION_MECHANISM_CACHE = (
    DATA_DIR
    / "gene_pathogenic_mechanism"
    / "prepared"
    / "hpo_condition_mechanism_cache.json.gz"
)
HPO_CONDITION_MECHANISM_SCHEMA_VERSION = "1.0"
PACKAGED_MECHANISM_JSON = (
    DATA_DIR / "gene_pathogenic_mechanism" / "prepared" / "gene_mechanism_curated_assertions.json"
)
CANONICAL_NONLOF_MECHANISM_JSON = (
    DATA_DIR
    / "gene_pathogenic_mechanism"
    / "prepared"
    / "gene_nonlof_mechanism_curated_assertions.json.gz"
)
DEFAULT_MECHANISM_JSON = (
    CANONICAL_NONLOF_MECHANISM_JSON
    if CANONICAL_NONLOF_MECHANISM_JSON.exists()
    else PACKAGED_MECHANISM_JSON
)
DEFAULT_DDG2P_MECHANISM_EVIDENCE = (
    DATA_DIR / "gene_pathogenic_mechanism" / "prepared" / "gene_pathogenic_mechanism_evidence.tsv"
)
DEFAULT_GOFCARDS_EXACT_GOF_HGVSP = DATA_DIR / "gofcards" / "gofcards_exact_gof.json.gz"
DEFAULT_GOFCARDS_RAW_XLSX = (
    DATA_DIR / "gene_pathogenic_mechanism" / "raw" / "gofcards" / "gofcards_data_download.xlsx"
)

LOOKUP_FIELD_PRIORITY = (
    "symbol",
    "hgnc_id",
    "ensembl_gene_id",
    "entrez_id",
    "prev_symbol",
    "alias_symbol",
    "refseq_accession",
    "uniprot_ids",
    "mane_select",
)

CANONICAL_MECHANISMS = {"LOF", "GOF", "DOMINANT_NEGATIVE", "TRIPLOSENSITIVITY"}
EXACT_SEQUENCE_MECHANISMS = {"GOF", "DOMINANT_NEGATIVE"}
VARIANT_MECHANISM_SCORE_KEYS = ("LOF", "GOF", "DOMINANT_NEGATIVE")
# Every source that may state a condition's mechanism. GoFCards is included:
# a curated gain-of-function allele in a gene IS that gene's curated history
# for the condition it was curated against, and 97 mechanism blocks have no
# other source at all. The old worry -- that one allele would give every
# variant in the gene a gain-of-function history -- does not arise, because
# select_condition_histories_for_variant keeps only the histories the query
# variant's own mechanism reaches. A predicted loss-of-function allele never
# sees the gain-of-function history; an unresolved one sees it as possible,
# not established, which is exactly what it is.
CONDITION_MECHANISM_SOURCES = {
    "G2P_DDG2P",
    "Orphadata",
    "ClinGen_haploinsufficiency",
    "GoFCards_exact+ClinVar_VCV",
    "deduced_from_inheritance",
}
CONDITION_MECHANISM_EVIDENCE_COLUMNS = {
    "gene_symbol",
    "source",
    "source_record_id",
    "source_condition_id",
    "mondo_id",
    "disease_scope",
    "priva_scope",
    "scope_review_status",
    "disease_label",
    "inheritance",
    "patho_mode_raw",
    "normalized_mechanisms",
    "mechanism_confidence",
    "disease_confidence",
    "pmids",
}
HPO_INHERITANCE_TERMS = {
    "HP:0000006": "Autosomal dominant inheritance",
    "HP:0000007": "Autosomal recessive inheritance",
    "HP:0001417": "X-linked inheritance",
    "HP:0001419": "X-linked recessive inheritance",
    "HP:0001423": "X-linked dominant inheritance",
    "HP:0001426": "Non-Mendelian inheritance",
    "HP:0001427": "Mitochondrial inheritance",
    "HP:0001450": "Y-linked inheritance",
    "HP:0010982": "Polygenic inheritance",
    "HP:0010983": "Oligogenic inheritance",
    "HP:0010984": "Digenic inheritance",
    "HP:0012275": "Autosomal dominant inheritance with maternal imprinting",
    "HP:0034341": "Pseudoautosomal recessive inheritance",
}
HPO_CACHE_INHERITANCE_LABELS = {
    "autosomal_dominant": "Autosomal dominant inheritance",
    "autosomal_recessive": "Autosomal recessive inheritance",
    "x_linked": "X-linked inheritance",
    "x_linked_recessive": "X-linked recessive inheritance",
    "x_linked_dominant": "X-linked dominant inheritance",
    "non_mendelian": "Non-Mendelian inheritance",
    "mitochondrial": "Mitochondrial inheritance",
    "y_linked": "Y-linked inheritance",
    "polygenic": "Polygenic inheritance",
    "oligogenic": "Oligogenic inheritance",
    "digenic": "Digenic inheritance",
    "autosomal_dominant_maternal_imprinting": (
        "Autosomal dominant inheritance with maternal imprinting"
    ),
    "pseudoautosomal_recessive": "Pseudoautosomal recessive inheritance",
}
HPO_PENETRANCE_TERMS = {
    "HP:0003829",  # Incomplete penetrance
    "HP:0034950",  # Complete penetrance
    "HP:4000159",  # Moderate penetrance
    "HP:4000160",  # Low penetrance
}
HPO_ONSET_TERMS = {
    "HP:0030674",  # Antenatal onset
    "HP:0011460",  # Embryonal onset
    "HP:0011461",  # Fetal onset
    "HP:0003577",  # Congenital onset
    "HP:0003623",  # Neonatal onset
    "HP:0003593",  # Infantile onset
    "HP:0011463",  # Childhood onset
    "HP:0003621",  # Juvenile onset
    "HP:0410280",  # Pediatric onset
    "HP:0003581",  # Adult onset
    "HP:0011462",  # Young adult onset
    "HP:0003596",  # Middle age onset
    "HP:0003584",  # Late onset
}
DDG2P_LOF_RAW_VALUES = {
    "loss of function",
    "absent gene product",
    "decreased gene product level",
    "reduced gene product level",
}
DDG2P_USABLE_MECHANISM_CONFIDENCE = {"high", "moderate"}
DDG2P_USABLE_DISEASE_CONFIDENCE = {"definitive", "strong", "moderate"}
DDG2P_RECESSIVE_LOF_INHERITANCE = {
    "biallelic_autosomal",
    "monoallelic_X",
    "monoallelic_X_hemizygous",
}
DDG2P_DOMINANT_LOF_INHERITANCE = {
    "monoallelic_autosomal",
    "monoallelic_X_heterozygous",
    "monoallelic_Y_hemizygous",
}
VARIANT_MECHANISM_OUTPUT_COLUMNS = (
    # Step 1 -- the variant's own mechanism, as three scores.
    "variant_effect",
    "variant_lof_score",
    "variant_gof_score",
    "variant_dn_score",
    "variant_mechanism_exclusive",
    "variant_exact_mechanisms",
    # Steps 2 and 3 -- the histories that mechanism reaches, split by whether
    # the mechanism is established for them or only possible.
    "variant_mechanism_applicable",
    "variant_mechanism_uncertain",
    # The three facts the chain exists to deliver, at variant resolution.
    "variant_condition_ids",
    # condition|mechanism|inheritance|penetrance, one entry per history, so the
    # pairing between a condition and its inheritance and penetrance survives.
    "variant_condition_histories",
    "variant_inheritance",
    # How that inheritance was arrived at: matched to this variant's own
    # mechanism, taken from a gene whose every germline condition agrees, or
    # inferred from the gene's constraint data when HPO states none.
    "variant_inheritance_basis",
    "variant_x_linked",
    "variant_penetrance",
)
logger = logging.getLogger(__name__)


def _clean(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "<na>"}:
        return ""
    return text


def _bool_value(value: Any) -> bool:
    """Coerce common dataframe boolean encodings without treating NaN as true."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return value != 0
    return _clean(value).lower() in {"1", "true", "t", "yes", "y"}


def _split_annotation_tokens(value: Any) -> set[str]:
    return {
        token.strip().upper()
        for token in re.split(r"[&,;|]", _clean(value))
        if token.strip()
    }


def _bounded_mechanism_score(value: Any) -> int:
    try:
        score = int(float(_clean(value)))
    except ValueError:
        return 0
    return score if score in {0, 1, 2} else 0


def _exact_mechanisms_from_effect(value: Any) -> set[str]:
    effect = _clean(value)
    prefix = "exact_known_"
    if not effect.startswith(prefix):
        return set()
    return {
        _normalize_sequence_mechanism(token)
        for token in effect[len(prefix):].split("+")
        if _normalize_sequence_mechanism(token) in EXACT_SEQUENCE_MECHANISMS
    }


def infer_query_variant_effect(row: dict[str, Any] | pd.Series) -> dict[str, Any]:
    """Infer one exclusive variant effect while retaining suppressed predictions.

    Exact curated GoF or dominant-negative evidence has score ``2`` and defines
    the allowed mechanism set for this allele. Generic consequence, NMD, or
    LOFTEE signals remain visible as audit evidence but cannot create a competing
    LoF call. Without an exact assertion, prediction-based LoF evidence receives
    score ``1``. LOFTEE ``OS`` is its other-splice category.
    """
    predicted_lof_evidence: list[str] = []
    loftee_tokens = _split_annotation_tokens(row.get("LoF"))
    if "HC" in loftee_tokens:
        predicted_lof_evidence.append("LOFTEE_HC")
    if "OS" in loftee_tokens:
        predicted_lof_evidence.append("LOFTEE_OS")

    consequence = _clean(row.get("Consequence")).lower()
    nmd = _clean(row.get("NMD")).lower()
    lof_filter = _clean(row.get("LoF_filter")).upper()
    truncating = "stop_gained" in consequence or "frameshift" in consequence
    escapes_nmd = "escaping" in nmd or "END_TRUNC" in lof_filter
    # A splice variant that shifts the reading frame destroys the transcript
    # the same way a coding frameshift does, so it is judged by the same
    # nonsense-mediated decay question rather than by being "a splice variant".
    splice_frameshift = _bool_value(row.get("splicing_frameshift"))
    nmd_lof = (truncating or splice_frameshift) and not escapes_nmd
    if nmd_lof:
        predicted_lof_evidence.append("NMD_PREDICTED_LOF")
    elif truncating or splice_frameshift:
        # Escaping decay weakens the claim, it does not withdraw it: the
        # protein is still truncated, so this remains plausible loss of
        # function. Without its own token the variant would fall through to
        # "uncertain" and score 0, which is reserved for no evidence at all.
        predicted_lof_evidence.append("NMD_ESCAPING_TRUNCATION")
    if _bool_value(row.get("vep_consq_lof")):
        predicted_lof_evidence.append("VEP_LOF")
    if _bool_value(row.get("splicing_lof")):
        predicted_lof_evidence.append("PRIVA_SPLICE_LOF")
    if _bool_value(row.get("5UTR_lof")):
        predicted_lof_evidence.append("PRIVA_5UTR_LOF")

    scores = {
        "LOF": _bounded_mechanism_score(row.get("variant_lof_score")),
        "GOF": _bounded_mechanism_score(row.get("variant_gof_score")),
        "DOMINANT_NEGATIVE": _bounded_mechanism_score(
            row.get("variant_dn_score")
        ),
    }
    gof_tokens = _split_annotation_tokens(row.get("variant_gof_tag"))
    if "GOF" in gof_tokens:
        scores["GOF"] = 2

    exact_mechanisms = {
        mechanism for mechanism, score in scores.items() if score == 2
    } & EXACT_SEQUENCE_MECHANISMS
    predicted_lof_evidence = list(dict.fromkeys(predicted_lof_evidence))
    suppressed_evidence: list[str] = []
    evidence: list[str] = []
    if exact_mechanisms:
        for mechanism in scores:
            scores[mechanism] = 2 if mechanism in exact_mechanisms else 0
        evidence.extend(
            f"CANONICAL_EXACT_{mechanism}" for mechanism in sorted(exact_mechanisms)
        )
        suppressed_evidence = predicted_lof_evidence
        effect = "exact_known_" + "+".join(sorted(exact_mechanisms))
    elif predicted_lof_evidence:
        # Predicted loss of function is graded, not flat. A transcript that is
        # destroyed by nonsense-mediated decay is established loss of function,
        # not merely a plausible one, so it scores 2. A transcript that escapes
        # decay may still make a partly functional protein, so it scores 1 --
        # unless LOFTEE's high-confidence call rescues it, which is a direct
        # statement that the loss holds despite the escape.
        established_lof = (
            "NMD_PREDICTED_LOF" in predicted_lof_evidence
            or "LOFTEE_HC" in predicted_lof_evidence
        )
        scores["LOF"] = max(scores["LOF"], 2 if established_lof else 1)
        evidence.extend(predicted_lof_evidence)
        effect = "predicted_LOF_high_confidence"
    else:
        effect = "uncertain"
    return {
        "variant_effect": effect,
        "variant_effect_evidence": ";".join(evidence),
        "variant_effect_suppressed_evidence": ";".join(suppressed_evidence),
        "variant_lof_score": scores["LOF"],
        "variant_gof_score": scores["GOF"],
        "variant_dn_score": scores["DOMINANT_NEGATIVE"],
        "variant_mechanism_exclusive": bool(exact_mechanisms),
        "variant_exact_mechanisms": ";".join(sorted(exact_mechanisms)),
    }


def _norm(value: Any) -> str:
    return _clean(value).upper()


def _norm_chrom(value: Any) -> str:
    text = _clean(value)
    if not text:
        return ""
    return text if text.startswith("chr") else f"chr{text}"


def _norm_chrom_key(value: Any) -> str:
    text = _clean(value)
    if not text:
        return ""
    if text.lower().startswith("chr"):
        text = text[3:]
    if text.upper() in {"M", "MT"}:
        return "MT"
    return text.upper()


def _norm_pos_key(value: Any) -> str:
    text = _clean(value)
    if not text:
        return ""
    try:
        as_float = float(text)
    except ValueError:
        return text
    if as_float.is_integer():
        return str(int(as_float))
    return text


def _norm_allele_key(value: Any) -> str:
    return _clean(value).upper()


def _genomic_match_key(chrom: Any, pos: Any, ref: Any, alt: Any) -> str:
    chrom_key = _norm_chrom_key(chrom)
    pos_key = _norm_pos_key(pos)
    if not chrom_key or not pos_key:
        return ""
    return "|".join(
        [
            chrom_key,
            pos_key,
            _norm_allele_key(ref),
            _norm_allele_key(alt),
        ]
    )


def _genomic_match_key_from_text(value: Any) -> str:
    text = _clean(value)
    if not text:
        return ""
    parts = text.split("|")
    if len(parts) != 4:
        return ""
    return _genomic_match_key(parts[0], parts[1], parts[2], parts[3])


def _normalize_assembly(value: Any) -> str:
    aliases = {
        "grch37": "hg19",
        "hg19": "hg19",
        "b37": "hg19",
        "grch38": "hg38",
        "hg38": "hg38",
        "b38": "hg38",
    }
    cleaned = _clean(value).lower()
    return aliases.get(cleaned, cleaned)


def _requested_genomic_key_types(value: Any) -> tuple[str, ...]:
    cleaned = _clean(value).lower() or "auto"
    if cleaned == "auto":
        return ("vcf", "genomic")
    if cleaned in {"vcf", "genomic"}:
        return (cleaned,)
    return ()


def _norm_hgvsp(value: Any) -> str:
    text = _clean(value)
    if not text:
        return ""
    protein = text.split(":")[-1].strip()
    if protein.startswith("p."):
        protein = protein[2:]
    return protein.upper()


def _normalize_sequence_mechanism(value: Any) -> str:
    """Normalize an explicitly curated sequence-variant mechanism label."""
    token = re.sub(r"[^A-Z0-9]+", "_", _clean(value).upper()).strip("_")
    aliases = {
        "GAIN_OF_FUNCTION": "GOF",
        "ACTIVATING": "GOF",
        "DN": "DOMINANT_NEGATIVE",
        "DOMINANTNEGATIVE": "DOMINANT_NEGATIVE",
    }
    return aliases.get(token, token)


def _split_multi(value: Any) -> list[str]:
    text = _clean(value)
    if not text:
        return []
    return [part.strip() for part in text.split("|") if part.strip()]


def _split_pmids(value: Any) -> list[str]:
    text = _clean(value)
    if not text:
        return []
    return [part.strip() for part in re.split(r"[;|,]", text) if part.strip()]


def _safe_int(value: Any) -> int | None:
    text = _clean(value)
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def _safe_float(value: Any) -> float:
    text = _clean(value)
    if not text:
        return np.nan
    try:
        return float(text)
    except ValueError:
        return np.nan


# The former gene-only category classifier is deliberately disabled. Mechanism
# and inheritance are now resolved together in variant-level assertions.


def _hpo_inheritance_flags(inheritance_modes: Any) -> dict[str, bool]:
    """Parse the HPO inheritance mode text used by the collapsed gene table."""
    modes = {_clean(part) for part in _clean(inheritance_modes).split(";") if _clean(part)}
    return {
        "recessive": bool(
            {
                "Autosomal recessive inheritance",
                "X-linked recessive inheritance",
                "X-linked inheritance",
            }
            & modes
        ),
        "dominant": bool(
            {
                "Autosomal dominant inheritance",
                "X-linked dominant inheritance",
                "Y-linked inheritance",
            }
            & modes
        ),
    }


def build_hpo_condition_index(
    assertions: pd.DataFrame,
) -> dict[tuple[str, str], dict[str, Any]]:
    """Index HPO evidence by the inseparable gene-plus-disease identity.

    The index deliberately requires both a gene symbol and a stable disease ID.
    A MONDO ID alone is insufficient because one disease can involve multiple
    genes with different mechanisms. A gene symbol alone is also insufficient
    because the same gene can cause several disorders with different
    inheritance, penetrance, or onset.

    Every original HPO assertion is retained with its frequency, evidence code,
    and reference. Compact inheritance, penetrance, and onset lists are derived
    only for convenient downstream decisions; they never replace the auditable
    source rows. Each condition can be looked up by its native OMIM/ORPHA ID or
    by MONDO, allowing G2P and Orphadata to use the same condition record.
    """
    required = {
        "gene_symbol",
        "disease_id",
        "hpo_id",
        "frequency",
        "evidence",
        "reference",
    }
    missing = sorted(required - set(assertions.columns))
    if missing:
        raise ValueError(f"HPO assertion table missing columns: {', '.join(missing)}")

    frame = assertions.copy().fillna("")
    for column in (
        "mondo_id",
        "disease_scope",
        "priva_scope",
        "scope_review_status",
    ):
        if column not in frame.columns:
            frame[column] = ""
    frame = frame.loc[
        frame["gene_symbol"].astype(str).str.strip().ne("")
        & frame["disease_id"].astype(str).str.strip().ne("")
    ]

    index: dict[tuple[str, str], dict[str, Any]] = {}
    group_columns = [
        "gene_symbol",
        "disease_id",
        "mondo_id",
        "disease_scope",
        "priva_scope",
        "scope_review_status",
    ]
    for group_key, group in frame.groupby(group_columns, sort=False, dropna=False):
        (
            gene_symbol,
            disease_id,
            mondo_id,
            disease_scope,
            priva_scope,
            scope_review_status,
        ) = map(_clean, group_key)
        hpo_rows = [
            {
                "hpo_id": _clean(row.get("hpo_id")),
                "frequency": _clean(row.get("frequency")),
                "evidence": _clean(row.get("evidence")),
                "reference": _clean(row.get("reference")),
            }
            for row in group.to_dict(orient="records")
        ]
        hpo_ids = list(
            dict.fromkeys(row["hpo_id"] for row in hpo_rows if row["hpo_id"])
        )
        record = {
            "disease_id": disease_id,
            "mondo_id": mondo_id,
            "disease_scope": disease_scope,
            "priva_scope": priva_scope,
            "scope_review_status": scope_review_status,
            "hpo_ids": hpo_ids,
            "inheritance_modes": [
                HPO_INHERITANCE_TERMS[hpo_id]
                for hpo_id in hpo_ids
                if hpo_id in HPO_INHERITANCE_TERMS
            ],
            "penetrance_hpo_ids": [
                hpo_id for hpo_id in hpo_ids if hpo_id in HPO_PENETRANCE_TERMS
            ],
            "onset_hpo_ids": [
                hpo_id for hpo_id in hpo_ids if hpo_id in HPO_ONSET_TERMS
            ],
            "hpo_assertions": hpo_rows,
        }
        for condition_id in (disease_id, mondo_id):
            if condition_id:
                index[(gene_symbol.upper(), condition_id.upper())] = record
    return index


class _LocalHgncResolver:
    """Small fallback resolver using PriVA's local HGNC table."""

    def __init__(self, hgnc_table: Path = DEFAULT_HGNC_TABLE) -> None:
        self.hgnc_table = Path(hgnc_table)
        self.rows: list[dict[str, str]] = []
        self.lookup: dict[str, list[tuple[int, dict[str, str]]]] = defaultdict(list)
        self._load()

    def _load(self) -> None:
        with open(self.hgnc_table, encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.rows = [dict(row) for row in reader]
        for row in self.rows:
            for priority, field in enumerate(LOOKUP_FIELD_PRIORITY):
                for value in _split_multi(row.get(field, "")):
                    key = _norm(value)
                    if key:
                        self.lookup[key].append((priority, row))

    def resolve(self, query: Any, *, warn_ambiguous: bool = True) -> str:
        cleaned = _clean(query)
        if not cleaned:
            return ""
        matches = self.lookup.get(_norm(cleaned), [])
        if not matches:
            return cleaned
        matches = sorted(matches, key=lambda item: (item[0], item[1].get("symbol", "")))
        best_priority = matches[0][0]
        best = [row for priority, row in matches if priority == best_priority]
        symbols = sorted({row.get("symbol", "") for row in best if row.get("symbol", "")})
        if len(symbols) > 1:
            if warn_ambiguous:
                logger.warning(
                    "Ambiguous HGNC query %s resolves to %s; keeping input symbol",
                    cleaned,
                    ",".join(symbols),
                )
            return cleaned
        return symbols[0] if symbols else cleaned


class GeneMechanismHub:
    """Normalize symbols, then summarize gene mechanism and inheritance history."""

    def __init__(
        self,
        *,
        condition_cache: str | Path = DEFAULT_HPO_CONDITION_MECHANISM_CACHE,
        mechanism_json: str | Path = DEFAULT_MECHANISM_JSON,
        ddg2p_evidence: str | Path = DEFAULT_DDG2P_MECHANISM_EVIDENCE,
        hpo_collapsed: str | Path = DEFAULT_HPO_COLLAPSED,
        clingen_dosage: str | Path = DEFAULT_CLINGEN_DOSAGE,
        loeuf_table: str | Path = DEFAULT_LOEUF_TABLE,
        hgnc_table: str | Path = DEFAULT_HGNC_TABLE,
        use_hgnc_package: bool = True,
    ) -> None:
        self.condition_cache = Path(condition_cache)
        self.mechanism_json = Path(mechanism_json)
        self.ddg2p_evidence = Path(ddg2p_evidence) if ddg2p_evidence else Path("")
        self.hpo_collapsed = Path(hpo_collapsed)
        self.clingen_dosage = Path(clingen_dosage)
        self.loeuf_table = Path(loeuf_table)
        self.hgnc_table = Path(hgnc_table)

        self._resolver = self._build_resolver(use_hgnc_package)
        self._condition_cache_by_symbol: dict[str, dict[str, Any]] | None = None
        self._condition_cache_meta: dict[str, Any] | None = None
        self._mechanism_by_symbol: dict[str, dict[str, Any]] | None = None
        self._condition_mechanism_by_symbol: dict[str, list[dict[str, Any]]] | None = None
        self._ddg2p_lof_by_symbol: dict[str, list[dict[str, Any]]] | None = None
        self._hpo_by_symbol: dict[str, dict[str, str]] | None = None
        self._hpo_condition_index: dict[tuple[str, str], dict[str, Any]] | None = None
        self._clingen_by_symbol: dict[str, dict[str, Any]] | None = None
        self._loeuf_by_symbol: dict[str, float] | None = None
        self._canonical_exact_nonlof_rows: list[dict[str, Any]] | None = None
        self._gofcards_by_symbol_hgvsp: dict[tuple[str, str], list[dict[str, str]]] | None = None
        self._gofcards_by_symbol_hgvsc: dict[tuple[str, str], list[dict[str, str]]] | None = None
        self._gofcards_by_symbol_genomic: dict[tuple[str, str, str, str], list[dict[str, str]]] | None = None

    def _build_resolver(self, use_hgnc_package: bool) -> Any:
        if use_hgnc_package:
            try:
                from hgnc_symbol_resolver import HgncResolver

                return HgncResolver()
            except Exception:
                pass
        return _LocalHgncResolver(self.hgnc_table)

    def resolve_symbol(self, gene_symbol: Any) -> str:
        """Return one current HGNC symbol for a query string."""
        if hasattr(self._resolver, "resolve"):
            return self._resolver.resolve(gene_symbol)
        raise TypeError("HGNC resolver object does not expose resolve()")

    def _try_resolve_symbol(self, gene_symbol: Any) -> str:
        """Resolve cache-side symbols without failing on ambiguous legacy aliases."""
        try:
            if isinstance(self._resolver, _LocalHgncResolver):
                return self._resolver.resolve(gene_symbol, warn_ambiguous=False)
            return self.resolve_symbol(gene_symbol)
        except ValueError:
            return _clean(gene_symbol)

    def _resolved_symbol_key(self, gene_symbol: Any) -> str:
        return self.resolve_symbol(gene_symbol)

    def _load_condition_cache(self) -> dict[str, dict[str, Any]]:
        """Load and normalize PriVA's integrated condition-mechanism cache.

        The builder writes one atomic JSON document with a versioned top-level
        contract. Runtime loading is deliberately strict: using a stale schema
        could silently move inheritance or penetrance between conditions, which
        is not acceptable for ACMG criteria. Gene keys are indexed under both
        their cache spelling and current HGNC symbol so historical aliases do
        not force callers to understand resource-specific naming.
        """
        if self._condition_cache_by_symbol is not None:
            return self._condition_cache_by_symbol
        if not self.condition_cache.is_file():
            raise FileNotFoundError(
                "Integrated HPO condition-mechanism cache not found: "
                f"{self.condition_cache}. Run "
                "`bash scripts/install_utils.sh "
                "hpo_condition_mechanism_cache_install config.yaml`."
            )

        with open_text(self.condition_cache) as handle:
            payload = json.load(handle)
        if not isinstance(payload, dict):
            raise ValueError(f"{self.condition_cache} must contain a JSON object")
        meta = payload.get("_meta")
        if not isinstance(meta, dict):
            raise ValueError(f"{self.condition_cache} is missing object _meta")
        schema_version = _clean(meta.get("schema_version"))
        if schema_version != HPO_CONDITION_MECHANISM_SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported condition cache schema {schema_version!r}; expected "
                f"{HPO_CONDITION_MECHANISM_SCHEMA_VERSION!r}"
            )
        genes = payload.get("genes")
        if not isinstance(genes, dict):
            raise ValueError(f"{self.condition_cache} is missing object genes")

        by_symbol: dict[str, dict[str, Any]] = {}
        for cache_symbol, gene in genes.items():
            if not isinstance(gene, dict):
                raise ValueError(
                    f"Condition cache gene {cache_symbol!r} must be an object"
                )
            symbol = _clean(cache_symbol)
            if not symbol:
                continue
            by_symbol[symbol] = gene
            resolved = self._try_resolve_symbol(symbol)
            if resolved:
                by_symbol.setdefault(resolved, gene)

        self._condition_cache_meta = meta
        self._condition_cache_by_symbol = by_symbol
        return by_symbol

    def _load_mechanisms(self) -> dict[str, dict[str, Any]]:
        if self._mechanism_by_symbol is not None:
            return self._mechanism_by_symbol
        with open_text(self.mechanism_json) as handle:
            raw = json.load(handle)
        by_symbol: dict[str, dict[str, Any]] = {}
        for info in raw.values():
            symbol = _clean(info.get("symbol"))
            if not symbol:
                continue
            by_symbol[symbol] = info
            resolved = self._try_resolve_symbol(symbol)
            if resolved and resolved not in by_symbol:
                by_symbol[resolved] = info
        self._mechanism_by_symbol = by_symbol
        return by_symbol

    def _load_condition_mechanism_evidence(
        self,
    ) -> dict[str, list[dict[str, Any]]]:
        """Load condition-resolved G2P and Orphadata mechanism assertions.

        The configured mechanism JSON and this TSV have deliberately different
        responsibilities. The shared JSON preserves exact ClinVar/GoFCards
        variant relationships, whereas the TSV is the canonical runtime source
        for G2P and Orphadata's gene-condition-mechanism assertions. Reading the
        TSV here lets PriVA keep both resources without replacing the richer
        shared JSON or copying condition records into it.

        A complete identity schema is mandatory. Older evidence tables that do
        not contain source disease IDs, MONDO mappings, or PriVA disease-scope
        fields are rejected rather than being interpreted gene-wide. Such old
        tables remain usable by ``_load_ddg2p_lof`` for the narrowly defined
        legacy PVS1/inheritance audit, but they cannot transfer condition-level
        inheritance, onset, or penetrance.

        One source row can encode more than one normalized mechanism. Each
        canonical mechanism is emitted as its own record so downstream variant
        matching can select GOF or LOF without carrying an unrelated mechanism.
        No disease-name, PMID, or phenotype-similarity matching occurs here.
        """
        if self._condition_mechanism_by_symbol is not None:
            return self._condition_mechanism_by_symbol

        by_symbol: dict[str, list[dict[str, Any]]] = defaultdict(list)
        if not self.ddg2p_evidence or not self.ddg2p_evidence.exists():
            self._condition_mechanism_by_symbol = {}
            return self._condition_mechanism_by_symbol

        try:
            with open(self.ddg2p_evidence, encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                missing = sorted(
                    CONDITION_MECHANISM_EVIDENCE_COLUMNS
                    - set(reader.fieldnames or [])
                )
                if missing:
                    logger.warning(
                        "Condition mechanism evidence table %s has a stale schema; "
                        "missing columns: %s",
                        self.ddg2p_evidence,
                        ", ".join(missing),
                    )
                    self._condition_mechanism_by_symbol = {}
                    return self._condition_mechanism_by_symbol

                seen_by_symbol: dict[str, set[tuple[str, ...]]] = defaultdict(set)
                for row in reader:
                    source = _clean(row.get("source"))
                    symbol = _clean(row.get("gene_symbol"))
                    if source not in CONDITION_MECHANISM_SOURCES or not symbol:
                        continue

                    mechanisms = {
                        token.strip().upper()
                        for token in re.split(
                            r"[|;,]", _clean(row.get("normalized_mechanisms"))
                        )
                        if token.strip()
                    } & CANONICAL_MECHANISMS
                    if not mechanisms:
                        continue

                    keys = {symbol}
                    resolved = self._try_resolve_symbol(symbol)
                    if resolved:
                        keys.add(resolved)
                    for mechanism in sorted(mechanisms):
                        record = {
                            "level": "gene_level",
                            "source": source,
                            "mechanism": mechanism,
                            "mechanism_raw": _clean(row.get("patho_mode_raw")),
                            "pmids": _split_pmids(row.get("pmids", "")),
                            "source_record_id": _clean(
                                row.get("source_record_id")
                            ),
                            "source_condition_id": _clean(
                                row.get("source_condition_id")
                            ),
                            "mondo_id": _clean(row.get("mondo_id")),
                            "disease": _clean(row.get("disease_label")),
                            "disease_scope": _clean(row.get("disease_scope")),
                            "priva_scope": _clean(row.get("priva_scope")),
                            "scope_review_status": _clean(
                                row.get("scope_review_status")
                            ),
                            "allelic_requirement": _clean(
                                row.get("inheritance")
                            ),
                            "confidence": _clean(
                                row.get("mechanism_confidence")
                            ),
                            "mechanism_confidence": _clean(
                                row.get("mechanism_confidence")
                            ),
                            "disease_confidence": _clean(
                                row.get("disease_confidence")
                            ),
                        }
                        identity = (
                            source,
                            record["source_record_id"],
                            record["source_condition_id"],
                            record["mondo_id"],
                            mechanism,
                            record["allelic_requirement"],
                        )
                        for key in keys:
                            if identity in seen_by_symbol[key]:
                                continue
                            seen_by_symbol[key].add(identity)
                            by_symbol[key].append(record)
        except Exception as exc:
            logger.warning(
                "Failed to load condition mechanism evidence table %s: %s",
                self.ddg2p_evidence,
                exc,
            )
            self._condition_mechanism_by_symbol = {}
            return self._condition_mechanism_by_symbol

        self._condition_mechanism_by_symbol = {
            symbol: records for symbol, records in by_symbol.items() if records
        }
        return self._condition_mechanism_by_symbol

    @staticmethod
    def _is_strict_ddg2p_lof(row: pd.Series) -> bool:
        raw = _clean(row.get("patho_mode_raw")).lower()
        if raw.startswith("undetermined") or "non-loss" in raw:
            return False
        norm_tokens = {
            token.strip().upper()
            for token in re.split(r"[;,]", _clean(row.get("normalized_mechanisms")))
            if token.strip()
        }
        return raw in DDG2P_LOF_RAW_VALUES or norm_tokens == {"LOF"}

    def _load_ddg2p_lof(self) -> dict[str, list[dict[str, Any]]]:
        if self._ddg2p_lof_by_symbol is not None:
            return self._ddg2p_lof_by_symbol
        by_symbol: dict[str, list[dict[str, Any]]] = defaultdict(list)
        if not self.ddg2p_evidence or not self.ddg2p_evidence.exists():
            self._ddg2p_lof_by_symbol = {}
            return self._ddg2p_lof_by_symbol

        try:
            df = pd.read_csv(
                self.ddg2p_evidence,
                sep="\t",
                dtype=str,
                low_memory=False,
            ).fillna("")
        except Exception as exc:
            logger.warning("Failed to load DDG2P/G2P evidence table %s: %s", self.ddg2p_evidence, exc)
            self._ddg2p_lof_by_symbol = {}
            return self._ddg2p_lof_by_symbol

        if "source" in df.columns:
            df = df.loc[df["source"].eq("G2P_DDG2P")].copy()
        if df.empty or "gene_symbol" not in df.columns:
            self._ddg2p_lof_by_symbol = {}
            return self._ddg2p_lof_by_symbol

        for _, row in df.iterrows():
            if not self._is_strict_ddg2p_lof(row):
                continue
            mechanism_confidence = _clean(row.get("mechanism_confidence")).lower()
            disease_confidence = _clean(row.get("disease_confidence")).lower()
            if mechanism_confidence not in DDG2P_USABLE_MECHANISM_CONFIDENCE:
                continue
            if disease_confidence not in DDG2P_USABLE_DISEASE_CONFIDENCE:
                continue

            symbol = _clean(row.get("gene_symbol"))
            if not symbol:
                continue
            inheritance = _clean(row.get("inheritance"))
            record = {
                "level": "gene_level",
                "source": "G2P_DDG2P",
                "mechanism": "LOF",
                "pmids": _split_pmids(row.get("pmids", "")),
                "disease": _clean(row.get("disease_label")),
                "inheritance": inheritance,
                "allelic_requirement": inheritance,
                "confidence": mechanism_confidence,
                "mechanism_confidence": mechanism_confidence,
                "disease_confidence": disease_confidence,
                "source_record_id": _clean(row.get("source_record_id")),
                "source_condition_id": _clean(row.get("source_condition_id")),
                "mondo_id": _clean(row.get("mondo_id")),
                "disease_scope": _clean(row.get("disease_scope")),
                "priva_scope": _clean(row.get("priva_scope")),
                "scope_review_status": _clean(row.get("scope_review_status")),
            }
            resolved = self._try_resolve_symbol(symbol)
            keys = {symbol}
            if resolved:
                keys.add(resolved)
            for key in keys:
                by_symbol[key].append(record)

        self._ddg2p_lof_by_symbol = {
            symbol: records for symbol, records in by_symbol.items() if records
        }
        return self._ddg2p_lof_by_symbol

    def _load_hpo(self) -> dict[str, dict[str, str]]:
        if self._hpo_by_symbol is not None:
            return self._hpo_by_symbol
        df = pd.read_csv(self.hpo_collapsed, sep="\t", dtype=str, low_memory=False).fillna("")
        by_symbol: dict[str, dict[str, str]] = {}
        assertion_columns = {
            "gene_symbol",
            "disease_id",
            "hpo_id",
            "frequency",
            "evidence",
            "reference",
        }
        if assertion_columns.issubset(df.columns):
            raw_condition_index = build_hpo_condition_index(df)
            condition_index: dict[tuple[str, str], dict[str, Any]] = {}
            for (symbol, condition_id), record in raw_condition_index.items():
                condition_index[(symbol, condition_id)] = record
                resolved = self._try_resolve_symbol(symbol)
                if resolved:
                    condition_index.setdefault(
                        (resolved.upper(), condition_id),
                        record,
                    )
            self._hpo_condition_index = condition_index
            hpo_records = []
            for symbol, group in df.groupby("gene_symbol", sort=False):
                if "priva_scope" in group.columns:
                    included_group = group.loc[group["priva_scope"].eq("include")]
                    review_diseases = list(
                        dict.fromkeys(
                            filter(
                                None,
                                map(
                                    _clean,
                                    group.loc[
                                        group["priva_scope"].eq("review"), "disease_id"
                                    ],
                                ),
                            )
                        )
                    )
                    excluded_diseases = list(
                        dict.fromkeys(
                            filter(
                                None,
                                map(
                                    _clean,
                                    group.loc[
                                        group["priva_scope"].eq("exclude"), "disease_id"
                                    ],
                                ),
                            )
                        )
                    )
                else:
                    included_group = group
                    review_diseases = []
                    excluded_diseases = []
                hpo_ids = list(
                    dict.fromkeys(filter(None, map(_clean, included_group["hpo_id"])))
                )
                inheritance_rows = included_group[
                    included_group["hpo_id"].isin(HPO_INHERITANCE_TERMS)
                ]
                inheritance_modes = list(
                    dict.fromkeys(
                        HPO_INHERITANCE_TERMS[hpo_id]
                        for hpo_id in inheritance_rows["hpo_id"]
                    )
                )
                inheritance_diseases = list(
                    dict.fromkeys(filter(None, map(_clean, inheritance_rows["disease_id"])))
                )
                hpo_records.append(
                    (
                        symbol,
                        {
                            "ncbi_gene_id": "",
                            "hpo_id": ";".join(hpo_ids),
                            "inheritance_modes": ";".join(inheritance_modes),
                            "inheritance_disease_ids": ";".join(inheritance_diseases),
                            "scope_review_required": bool(review_diseases),
                            "scope_review_disease_ids": ";".join(review_diseases),
                            "scope_excluded_disease_ids": ";".join(excluded_diseases),
                        },
                    )
                )
        else:
            self._hpo_condition_index = {}
            hpo_records = []
            for _, row in df.iterrows():
                symbol = _clean(row.get("gene_symbol"))
                if not symbol:
                    continue
                hpo_records.append(
                    (
                        symbol,
                        {
                            "ncbi_gene_id": _clean(row.get("ncbi_gene_id")),
                            "hpo_id": _clean(row.get("hpo_id")),
                            "inheritance_modes": _clean(row.get("inheritance_modes")),
                            "inheritance_disease_ids": _clean(
                                row.get("inheritance_disease_ids")
                            ),
                        },
                    )
                )

        for symbol, record in hpo_records:
            symbol = _clean(symbol)
            if not symbol:
                continue
            by_symbol.setdefault(symbol, record)
            resolved = self._try_resolve_symbol(symbol)
            if resolved:
                by_symbol.setdefault(resolved, record)
        self._hpo_by_symbol = by_symbol
        return by_symbol

    def _load_clingen(self) -> dict[str, dict[str, Any]]:
        if self._clingen_by_symbol is not None:
            return self._clingen_by_symbol
        df = pd.read_csv(self.clingen_dosage, sep="\t", dtype=str, low_memory=False).fillna("")
        by_symbol: dict[str, dict[str, Any]] = {}
        for _, row in df.iterrows():
            symbol = _clean(row.get("#Gene Symbol"))
            if not symbol:
                continue
            record = {
                "haploinsufficiency_score": _safe_int(row.get("Haploinsufficiency Score")),
                "haploinsufficiency_description": _clean(
                    row.get("Haploinsufficiency Description")
                ),
                "triplosensitivity_score": _safe_int(row.get("Triplosensitivity Score")),
                "triplosensitivity_description": _clean(
                    row.get("Triplosensitivity Description")
                ),
            }
            by_symbol.setdefault(symbol, record)
            resolved = self._try_resolve_symbol(symbol)
            if resolved:
                by_symbol.setdefault(resolved, record)
        self._clingen_by_symbol = by_symbol
        return by_symbol

    def _load_loeuf(self) -> dict[str, float]:
        if self._loeuf_by_symbol is not None:
            return self._loeuf_by_symbol
        df = pd.read_csv(
            self.loeuf_table,
            sep="\t",
            dtype=str,
            low_memory=False,
            usecols=["#gene", "canonical", "oe_lof_upper"],
        ).fillna("")
        df["__loeuf"] = df["oe_lof_upper"].map(_safe_float)
        out: dict[str, float] = {}
        canonical = df[df["canonical"].str.lower().eq("true")]
        for symbol, grp in canonical.groupby("#gene"):
            values = [v for v in grp["__loeuf"] if not math.isnan(v)]
            if values:
                out.setdefault(symbol, min(values))
        for symbol, grp in df.groupby("#gene"):
            if symbol in out:
                continue
            values = [v for v in grp["__loeuf"] if not math.isnan(v)]
            if values:
                out[symbol] = min(values)
        resolved_out = dict(out)
        for symbol, value in out.items():
            resolved = self._try_resolve_symbol(symbol)
            if resolved:
                resolved_out.setdefault(resolved, value)
        self._loeuf_by_symbol = resolved_out
        return resolved_out

    # ------------------------------------------------------------------
    # Exact sequence-variant mechanism evidence
    # ------------------------------------------------------------------

    def _load_canonical_exact_nonlof_rows(self) -> list[dict[str, Any]]:
        """Return runtime-eligible GoFCards variant blocks from the canonical JSON.

        ``build_gene_nonlof_mechanism_cache.py`` embeds each variant under
        ``exact_normalized_variants`` exactly as the GoFCards cache stores it,
        so the block still has its ``record`` and its per-assembly ``genomic``
        and ``transcripts``. Runtime reads that structure and never reinterprets
        the upstream cache independently.

        ClinVar VCV entries link an already curated variant to clinical
        assertions; ClinVar pathogenicity alone establishes neither gain of
        function nor a dominant-negative effect, so those links are excluded
        here on purpose.
        """
        if self._canonical_exact_nonlof_rows is not None:
            return self._canonical_exact_nonlof_rows

        entries: list[dict[str, Any]] = []
        seen: set[str] = set()
        for gene in self._load_mechanisms().values():
            symbol = self._try_resolve_symbol(gene.get("symbol"))
            if not symbol:
                continue
            for wrapped in gene.get("variant_level", []) or []:
                if not isinstance(wrapped, dict):
                    continue
                for source, assertion in wrapped.items():
                    if source == "ClinVar_VCV" or not isinstance(assertion, dict):
                        continue
                    status = assertion.get("exact_normalization_status")
                    if status and status != "matched_gene_concordant":
                        continue
                    for variant in assertion.get("exact_normalized_variants", []) or []:
                        if not isinstance(variant, dict):
                            continue
                        record = variant.get("record") or {}
                        if (record.get("eligibility") or {}).get("status") != "eligible":
                            continue
                        variant_id = _clean(variant.get("variant_id"))
                        # One allele can be asserted by two curated sources with
                        # different mechanisms, so the identity that de-duplicates
                        # must include the source. Keying on the allele alone
                        # would silently drop the second mechanism.
                        marker = f"{source}|{symbol}|{variant_id}"
                        if not variant_id or marker in seen:
                            continue
                        seen.add(marker)
                        entries.append({
                            **variant,
                            "symbol": symbol,
                            # Carried from the wrapping assertion so a curated
                            # dominant-negative source is not reported as GoFCards
                            # gain of function.
                            "assertion_source": source,
                            "assertion_mechanism": _clean(assertion.get("mechanism")),
                            # Which canonical file this evidence came from, so an
                            # audit row identifies its own provenance.
                            "canonical_mechanism_json": str(self.mechanism_json),
                        })

        self._canonical_exact_nonlof_rows = entries
        return entries

    def _load_gofcards_variant_hgvsp(
        self,
        *,
        gofcards_exact_hgvsp_tsv: str | Path = DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
    ) -> dict[tuple[str, str], list[dict[str, str]]]:
        """Index eligible variants by (HGNC symbol, normalized protein change).

        A variant registers under every protein change its transcripts present,
        because one variant is numbered differently on different isoforms --
        JAK2 V617F is also Val468Phe and Val212Phe. Indexing only one of them
        would miss a query annotated on any other transcript.
        """
        if self._gofcards_by_symbol_hgvsp is not None:
            return self._gofcards_by_symbol_hgvsp

        index: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for variant in self._load_canonical_exact_nonlof_rows():
            symbol = variant["symbol"]
            for assembly, block in (variant.get("assemblies") or {}).items():
                for transcript, view in (block.get("transcripts") or {}).items():
                    for hgvsp, coding in (view.get("by_hgvsp") or {}).items():
                        key = _norm_hgvsp(hgvsp)
                        if not key:
                            continue
                        index[(symbol, key)].append(
                            self._gofcards_match_entry(
                                variant, assembly=assembly, transcript=transcript,
                                hgvsp=hgvsp, hgvsc=coding[0] if coding else "",
                            )
                        )

        self._gofcards_by_symbol_hgvsp = {
            key: self._deduplicate_gofcards_matches(rows) for key, rows in index.items()
        }
        return self._gofcards_by_symbol_hgvsp

    def _load_gofcards_variant_hgvsc(
        self,
    ) -> dict[tuple[str, str], list[dict[str, str]]]:
        """Index eligible variants by (HGNC symbol, coding change).

        Measured unique across the catalogue at 25,106 keys with no collisions,
        so a coding change identifies a variant within its gene where a protein
        change may not.
        """
        if self._gofcards_by_symbol_hgvsc is not None:
            return self._gofcards_by_symbol_hgvsc

        index: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for variant in self._load_canonical_exact_nonlof_rows():
            symbol = variant["symbol"]
            for assembly, block in (variant.get("assemblies") or {}).items():
                for transcript, view in (block.get("transcripts") or {}).items():
                    for hgvsc, detail in (view.get("by_hgvsc") or {}).items():
                        key = _clean(hgvsc).upper()
                        if not key:
                            continue
                        index[(symbol, key)].append(
                            self._gofcards_match_entry(
                                variant, assembly=assembly, transcript=transcript,
                                hgvsp=detail.get("hgvsp") or "", hgvsc=hgvsc,
                            )
                        )

        self._gofcards_by_symbol_hgvsc = {
            key: self._deduplicate_gofcards_matches(rows) for key, rows in index.items()
        }
        return self._gofcards_by_symbol_hgvsc

    @staticmethod
    def _gofcards_match_entry(
        variant: dict[str, Any], *, assembly: str = "", transcript: str = "",
        hgvsp: str = "", hgvsc: str = "",
    ) -> dict[str, str]:
        """Flatten one variant into the shape the matcher reports back."""
        record = variant.get("record") or {}
        source = record.get("source") or {}
        genomic = ((variant.get("assemblies") or {}).get(assembly) or {}).get("genomic") or {}
        # GoFCards is the only variant-level source in the shipped cache, so it
        # is the default; a curated source that declares its own mechanism keeps
        # that mechanism rather than being reported as gain of function.
        assertion_source = _clean(variant.get("assertion_source")) or "GoFCards"
        assertion_mechanism = _clean(variant.get("assertion_mechanism")) or "GOF"
        return {
            "source": assertion_source,
            "mechanism": assertion_mechanism,
            "canonical_assertion_source": assertion_source,
            "canonical_mechanism_json": _clean(variant.get("canonical_mechanism_json")),
            "symbol": variant.get("symbol", ""),
            "gofcards_variant_id": _clean(variant.get("variant_id")),
            "gofcards_allele_key": _clean(source.get("gofcards_allele_key")),
            "gofcards_accession_id": _clean(
                (record.get("annotations") or {}).get("clinvar_variation_id")),
            "assembly": assembly,
            "vep_transcript": transcript,
            "HGVSp": hgvsp,
            "HGVSc": hgvsc,
            "chrom": _clean(genomic.get("chrom")),
            "pos": _clean(genomic.get("pos")),
            "ref": _clean(genomic.get("ref")),
            "alt": _clean(genomic.get("alt")),
            "match_eligibility": _clean((record.get("eligibility") or {}).get("status")),
            "pmids": ";".join(
                e.get("pmid", "") for e in (record.get("evidence") or []) if e.get("pmid")
            ),
        }

    @staticmethod
    def _deduplicate_gofcards_matches(rows: list[dict[str, str]]) -> list[dict[str, str]]:
        seen: set[tuple[str, ...]] = set()
        unique_rows: list[dict[str, str]] = []
        for row in sorted(
            rows,
            key=lambda item: (
                item.get("gofcards_variant_id", ""),
                item.get("assembly", ""),
                item.get("vep_transcript", ""),
                item.get("HGVSp", ""),
            ),
        ):
            # The variant identifier already distinguishes variants, so a match
            # reported twice through different transcripts of the same variant
            # is one match, not two.
            dedup_key = (
                row.get("mechanism", ""),
                row.get("canonical_assertion_source", ""),
                row.get("gofcards_variant_id", ""),
            )
            if dedup_key in seen:
                continue
            seen.add(dedup_key)
            unique_rows.append(row)
        return unique_rows

    def _load_gofcards_variant_genomic(
        self,
        *,
        gofcards_exact_hgvsp_tsv: str | Path = DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
    ) -> dict[tuple[str, str, str, str], list[dict[str, str]]]:
        """Index eligible variants by (symbol, assembly, key type, chrom|pos|ref|alt).

        The coordinate is read from the assembly block it belongs to, so a
        variant that failed liftover simply contributes no hg38 key rather than
        an empty one. ``key_type`` is retained for the existing call contract;
        the cache stores a single normalized representation per build, so both
        types resolve to the same key.
        """
        if self._gofcards_by_symbol_genomic is not None:
            return self._gofcards_by_symbol_genomic

        index: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
        for variant in self._load_canonical_exact_nonlof_rows():
            symbol = variant["symbol"]
            for assembly, block in (variant.get("assemblies") or {}).items():
                coords = block.get("genomic") or {}
                if not coords.get("chrom") or coords.get("pos") is None:
                    continue
                key = _genomic_match_key(
                    coords.get("chrom"), coords.get("pos"),
                    coords.get("ref"), coords.get("alt"),
                )
                if not key:
                    continue
                entry = self._gofcards_match_entry(variant, assembly=assembly)
                for key_type in ("vcf", "genomic"):
                    index[(symbol, assembly, key_type, key)].append(entry)

        self._gofcards_by_symbol_genomic = {
            key: self._deduplicate_gofcards_matches(rows) for key, rows in index.items()
        }
        return self._gofcards_by_symbol_genomic

    def match_curated_nonlof_variant(
        self,
        gene_symbol: Any,
        *,
        hgvsp: Any = "",
        hgvsc: Any = "",
        chrom: Any = "",
        pos: Any = "",
        ref: Any = "",
        alt: Any = "",
        assembly: str = "hg38",
        key_type: str = "auto",
    ) -> dict[str, Any]:
        """Score exact curated GoF/DN evidence for one query allele.

        Scores use the agreed variant-level contract: ``0`` means no exact
        support, and ``2`` means an exact curated assertion that excludes
        unasserted mechanisms. Score ``1`` is reserved for prediction-based
        plausibility and is added later by :func:`infer_query_variant_effect`;
        this exact matcher does not infer loss of function from ClinVar
        pathogenicity or consequence alone.

        ``current HGNC gene + normalized HGVSp`` is an exact match route.
        Normalized genomic alleles are the fallback when HGVSp is absent or not
        found. More than one mechanism receives score ``2`` only if separate
        exact records explicitly assert those mechanisms for the same allele.
        """
        symbol = self._resolved_symbol_key(gene_symbol)
        hgvsp_key = _norm_hgvsp(hgvsp)
        hgvsc_key = _clean(hgvsc).split(":")[-1].upper()
        assembly_key = _normalize_assembly(assembly)
        genomic_key = _genomic_match_key(chrom, pos, ref, alt)
        scores = {mechanism: 0 for mechanism in VARIANT_MECHANISM_SCORE_KEYS}
        matches: list[dict[str, str]] = []
        route = ""
        matched_key_type = ""

        if symbol and hgvsp_key:
            matches = list(
                self._load_gofcards_variant_hgvsp().get(
                    (symbol, hgvsp_key),
                    [],
                )
            )
            if matches:
                route = "hgvsp"

        # Coding change is the second route: it separates two different coding
        # changes that produce the same protein change.
        if symbol and not matches and hgvsc_key:
            matches = list(
                self._load_gofcards_variant_hgvsc().get((symbol, hgvsc_key), [])
            )
            if matches:
                route = "hgvsc"

        if symbol and not matches and genomic_key and assembly_key:
            genomic_lookup = self._load_gofcards_variant_genomic()
            for candidate_key_type in _requested_genomic_key_types(key_type):
                candidate = genomic_lookup.get(
                    (
                        symbol,
                        assembly_key,
                        candidate_key_type,
                        genomic_key,
                    ),
                    [],
                )
                if candidate:
                    matches.extend(candidate)
                    route = "genomic"
                    matched_key_type = candidate_key_type
            matches = self._deduplicate_gofcards_matches(matches)

        matches = [
            match
            for match in matches
            if _normalize_sequence_mechanism(match.get("mechanism"))
            in EXACT_SEQUENCE_MECHANISMS
        ]
        mechanisms = sorted(
            {
                _normalize_sequence_mechanism(match.get("mechanism"))
                for match in matches
            }
        )
        for mechanism in mechanisms:
            scores[mechanism] = 2

        accession_ids = sorted(
            {
                _clean(match.get("gofcards_accession_id"))
                for match in matches
                if _clean(match.get("gofcards_accession_id"))
            }
        )
        variant_ids = sorted(
            {
                _clean(match.get("gofcards_variant_id"))
                for match in matches
                if _clean(match.get("gofcards_variant_id"))
            }
        )
        return {
            "input_symbol": _clean(gene_symbol),
            "symbol": symbol,
            "input_hgvsp": _clean(hgvsp),
            "matched_hgvsp_key": hgvsp_key if route == "hgvsp" else "",
            "matched_hgvsc_key": hgvsc_key if route == "hgvsc" else "",
            "input_assembly": _clean(assembly),
            "assembly": assembly_key,
            "input_genomic_key": genomic_key,
            "matched_key_type": matched_key_type,
            "match_route": route,
            "mechanism_scores": scores,
            "lof_score": scores["LOF"],
            "gof_score": scores["GOF"],
            "dn_score": scores["DOMINANT_NEGATIVE"],
            "mechanism_exclusive": bool(mechanisms),
            "exclusive_mechanisms": mechanisms,
            "matched_mechanisms": mechanisms,
            "gofcards_accession_id": ";".join(accession_ids),
            "gofcards_variant_id": ";".join(variant_ids),
            "matches": matches,
        }

    def match_gofcards_variant_gof(
        self,
        gene_symbol: Any,
        hgvsp: Any,
        *,
        gofcards_exact_hgvsp_tsv: str | Path = DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
        gofcards_step1_tsv: str | Path | None = None,
        gofcards_active_tsv: str | Path | None = None,
        gofcards_raw_xlsx: str | Path | None = None,
        gofcards_conversion_audit_tsv: str | Path | None = None,
    ) -> dict[str, Any]:
        """Compatibility view of canonical exact GoF HGVSp evidence."""
        result = self.match_curated_nonlof_variant(gene_symbol, hgvsp=hgvsp)
        gof_matches = [
            match
            for match in result["matches"]
            if _normalize_sequence_mechanism(match.get("mechanism")) == "GOF"
        ]
        result["matches"] = gof_matches
        result["variant_gof_tag"] = "GOF" if gof_matches else ""
        return result

    def match_gofcards_variant_gof_by_genomic(
        self,
        gene_symbol: Any,
        chrom: Any,
        pos: Any,
        ref: Any,
        alt: Any,
        *,
        assembly: str = "hg38",
        key_type: str = "auto",
        gofcards_exact_hgvsp_tsv: str | Path = DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
    ) -> dict[str, Any]:
        """Compatibility view of canonical exact GoF genomic evidence."""
        result = self.match_curated_nonlof_variant(
            gene_symbol,
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            assembly=assembly,
            key_type=key_type,
        )
        gof_matches = [
            match
            for match in result["matches"]
            if _normalize_sequence_mechanism(match.get("mechanism")) == "GOF"
        ]
        result.update(
            {
                "input_chrom": _clean(chrom),
                "input_pos": _clean(pos),
                "input_ref": _clean(ref),
                "input_alt": _clean(alt),
                "matched_genomic_key": result["input_genomic_key"],
                "variant_gof_tag": "GOF" if gof_matches else "",
                "matches": gof_matches,
            }
        )
        return result

    @staticmethod
    def _iter_mechanism_entries(info: dict[str, Any]) -> list[dict[str, Any]]:
        entries: list[dict[str, Any]] = []
        for level in ("gene_level", "variant_level"):
            for entry in info.get(level, []) or []:
                for source, data in entry.items():
                    mechanism = _clean(data.get("mechanism"))
                    if not mechanism:
                        continue
                    pmids = [_clean(p) for p in data.get("pmids", []) if _clean(p)]
                    entries.append(
                        {
                            "level": level,
                            "source": source,
                            "mechanism": mechanism,
                            "pmids": pmids,
                            "source_record_id": _clean(data.get("source_record_id")),
                            "source_condition_id": _clean(
                                data.get("source_condition_id")
                            ),
                            "mondo_id": _clean(data.get("mondo_id")),
                            "disease": _clean(data.get("disease")),
                            "disease_scope": _clean(data.get("disease_scope")),
                            "priva_scope": _clean(data.get("priva_scope")),
                            "scope_review_status": _clean(
                                data.get("scope_review_status")
                            ),
                            "consequence": _clean(data.get("consequence")),
                            "chrom": _clean(data.get("chr")),
                            "pos": _clean(data.get("pos")),
                            "ref": _clean(data.get("ref")),
                            "alt": _clean(data.get("alt")),
                            "transcript": _clean(data.get("transcript")),
                            "allelic_requirement": _clean(
                                data.get("inheritance")
                                or data.get("allelic_requirement")
                            ),
                            "confidence": _clean(
                                data.get("confidence")
                                or data.get("mechanism_confidence")
                            ),
                        }
                    )
        return entries

    def mechanism_history(
        self,
        gene_symbol: Any,
        *,
        include_entries: bool = False,
    ) -> dict[str, Any]:
        """Return condition-resolved and unresolved integrated-cache history."""
        symbol = self._resolved_symbol_key(gene_symbol)
        gene = self._load_condition_cache().get(symbol, {})
        entries = condition_cache_mechanism_entries(gene)
        mechanism_counts = Counter(entry["mechanism"] for entry in entries)
        pmids_by_mechanism: dict[str, set[str]] = defaultdict(set)
        variant_counts = Counter()
        gene_counts = Counter()
        for entry in entries:
            mechanism = entry["mechanism"]
            pmids_by_mechanism[mechanism].update(entry["pmids"])
            if entry["level"] == "variant_level":
                variant_counts[mechanism] += 1
            elif entry["level"] == "gene_level":
                gene_counts[mechanism] += 1

        out: dict[str, Any] = {
            "input_symbol": _clean(gene_symbol),
            "symbol": symbol,
            "mechanisms": sorted(mechanism_counts),
            "mechanism_counts": dict(sorted(mechanism_counts.items())),
            "gene_level_counts": dict(sorted(gene_counts.items())),
            "variant_level_counts": dict(sorted(variant_counts.items())),
            "pmids_by_mechanism": {
                key: sorted(value) for key, value in sorted(pmids_by_mechanism.items())
            },
            "has_lof_history": mechanism_counts.get("LOF", 0) > 0,
            "has_gof_history": mechanism_counts.get("GOF", 0) > 0,
            "has_dn_history": mechanism_counts.get("DOMINANT_NEGATIVE", 0) > 0,
            "has_triplosensitivity_history": mechanism_counts.get("TRIPLOSENSITIVITY", 0) > 0,
        }
        if include_entries:
            out["entries"] = entries
        return out

    def ddg2p_lof_history(
        self,
        gene_symbol: Any,
        *,
        include_entries: bool = False,
    ) -> dict[str, Any]:
        """Return DDG2P/G2P LoF mechanism history and allelic requirement flags."""
        symbol = self._resolved_symbol_key(gene_symbol)
        entries = [
            entry
            for entry in condition_cache_mechanism_entries(
                self._load_condition_cache().get(symbol, {})
            )
            if entry.get("level") == "gene_level"
            and entry.get("source") == "G2P_DDG2P"
            and entry.get("mechanism") == "LOF"
            and _clean(entry.get("mechanism_confidence")).lower()
            in DDG2P_USABLE_MECHANISM_CONFIDENCE
            and _clean(entry.get("disease_confidence")).lower()
            in DDG2P_USABLE_DISEASE_CONFIDENCE
        ]
        inheritance_counts = Counter(_clean(entry.get("inheritance")) for entry in entries)
        recessive = any(
            inheritance in DDG2P_RECESSIVE_LOF_INHERITANCE
            for inheritance in inheritance_counts
        )
        dominant = any(
            inheritance in DDG2P_DOMINANT_LOF_INHERITANCE
            for inheritance in inheritance_counts
        )
        out: dict[str, Any] = {
            "input_symbol": _clean(gene_symbol),
            "symbol": symbol,
            "has_ddg2p_lof_history": bool(entries),
            "ddg2p_lof_recessive": recessive,
            "ddg2p_lof_dominant": dominant,
            "ddg2p_lof_inheritance_counts": dict(sorted(inheritance_counts.items())),
            "ddg2p_lof_disease_count": len({entry.get("disease", "") for entry in entries if entry.get("disease", "")}),
        }
        if include_entries:
            out["entries"] = entries
        return out

    def known_inheritance_mode(
        self,
        gene_symbol: Any,
        *,
        gene_mean_am_score: float = np.nan,
    ) -> dict[str, Any]:
        """Return inheritance/HI calls using PriVA's existing inheritance function."""
        from acmg_criteria_assign import identify_inheritance_mode_per_row

        symbol = self._resolved_symbol_key(gene_symbol)
        hpo_record = self._load_hpo().get(symbol, {})
        clingen_record = self._load_clingen().get(symbol, {})
        loeuf = self._load_loeuf().get(symbol, np.nan)
        clingen_hi_score = clingen_record.get("haploinsufficiency_score")
        ddg2p_lof = self.ddg2p_lof_history(symbol)

        row_dict = {
            "Gene": symbol,
            "SYMBOL": symbol,
            "LOEUF": loeuf,
            "HPO_IDs": hpo_record.get("hpo_id", ""),
            "HPO_gene_inheritance": hpo_record.get("inheritance_modes") or None,
        }
        result = identify_inheritance_mode_per_row(
            row_dict,
            gene_mean_am_score,
            clingen_hi_score,
        )
        recessive, dominant, non_monogenic, non_mendelian, haplo_insufficient, incomplete = result
        effective_recessive = bool(recessive) or bool(ddg2p_lof["ddg2p_lof_recessive"])
        effective_dominant = bool(dominant) or bool(ddg2p_lof["ddg2p_lof_dominant"])
        effective_haplo = bool(haplo_insufficient) or bool(ddg2p_lof["ddg2p_lof_dominant"])
        return {
            "input_symbol": _clean(gene_symbol),
            "symbol": symbol,
            "recessive": effective_recessive,
            "dominant": effective_dominant,
            "non_monogenic": bool(non_monogenic),
            "non_mendelian": bool(non_mendelian),
            "haplo_insufficient": effective_haplo,
            "incomplete_penetrance": bool(incomplete),
            "hpo_inheritance_modes": hpo_record.get("inheritance_modes", ""),
            "hpo_inheritance_disease_ids": hpo_record.get("inheritance_disease_ids", ""),
            "hpo_scope_review_required": bool(
                hpo_record.get("scope_review_required", False)
            ),
            "hpo_scope_review_disease_ids": hpo_record.get(
                "scope_review_disease_ids", ""
            ),
            "hpo_scope_excluded_disease_ids": hpo_record.get(
                "scope_excluded_disease_ids", ""
            ),
            "ddg2p_lof_history": ddg2p_lof,
            "clingen_haploinsufficiency_score": clingen_hi_score,
            "clingen_haploinsufficiency_description": clingen_record.get(
                "haploinsufficiency_description", ""
            ),
            "clingen_triplosensitivity_score": clingen_record.get("triplosensitivity_score"),
            "clingen_triplosensitivity_description": clingen_record.get(
                "triplosensitivity_description", ""
            ),
            "loeuf": None if isinstance(loeuf, float) and math.isnan(loeuf) else loeuf,
            "decision_function": "acmg_criteria_assign.identify_inheritance_mode_per_row",
        }

    def gene_summary(self, gene_symbol: Any, *, include_entries: bool = False) -> dict[str, Any]:
        """Return normalized symbol, mechanism history, and known inheritance summary."""
        symbol = self._resolved_symbol_key(gene_symbol)
        return {
            "input_symbol": _clean(gene_symbol),
            "symbol": symbol,
            "mechanism_history": self.mechanism_history(symbol, include_entries=include_entries),
            "condition_mechanism_assertions": self.condition_mechanism_assertions(symbol),
            "ddg2p_lof_history": self.ddg2p_lof_history(symbol, include_entries=include_entries),
            "known_inheritance_mode": self.known_inheritance_mode(symbol),
        }

    def condition_mechanism_assertions(self, gene_symbol: Any) -> list[dict[str, Any]]:
        """Return automatic condition histories from the integrated cache."""
        symbol = self._resolved_symbol_key(gene_symbol)
        gene = self._load_condition_cache().get(symbol, {})
        return condition_cache_mechanism_assertions(gene)

def resolve_gene_symbol(gene_symbol: Any) -> str:
    """Convenience function returning one current HGNC symbol."""
    return GeneMechanismHub().resolve_symbol(gene_symbol)


def _hpo_allelic_requirements(inheritance_modes: Any) -> set[str]:
    """Collapse HPO inheritance to the supported compact assertion vocabulary.

    HPO establishes only a broad dominant or recessive history here. It does
    not establish autosomal/X-linked dosage or a molecular mechanism.
    """
    modes = {
        _clean(part)
        for part in _clean(inheritance_modes).split(";")
        if _clean(part)
    }
    requirements: set[str] = set()
    if {
        "Autosomal recessive inheritance",
        "X-linked inheritance",
        "X-linked recessive inheritance",
        "Pseudoautosomal recessive inheritance",
    } & modes:
        requirements.add("recessive")
    if {
        "Autosomal dominant inheritance",
        "Autosomal dominant inheritance with maternal imprinting",
        "X-linked dominant inheritance",
        "Y-linked inheritance",
    } & modes:
        requirements.add("dominant")
    if "Mitochondrial inheritance" in modes:
        requirements.add("mitochondrial")
    return requirements


def condition_cache_context(
    condition_key: Any,
    condition: dict[str, Any],
) -> dict[str, Any]:
    """Return the shared PriVA context for one included cache condition.

    An empty result means the condition cannot influence automatic ACMG logic.
    The cache retains only the HPO axes used here; ``hpo_assertion_count`` keeps
    the size of the complete phenotype record visible without loading it.
    """
    scope = condition.get("priva_scope", {})
    if not isinstance(scope, dict) or _clean(scope.get("decision")).lower() != "include":
        return {}

    identifiers = condition.get("identifiers", {})
    if not isinstance(identifiers, dict):
        identifiers = {}
    mondo_ids = identifiers.get("MONDO", [])
    mondo_id = _clean(mondo_ids[0]) if isinstance(mondo_ids, list) and mondo_ids else ""

    inheritance = condition.get("inheritance", {})
    inheritance = inheritance if isinstance(inheritance, dict) else {}
    hpo_inheritance_modes = [
        HPO_CACHE_INHERITANCE_LABELS.get(_clean(mode), _clean(mode))
        for mode in inheritance.get("modes", [])
        if _clean(mode)
    ]
    penetrance = condition.get("penetrance", {})
    penetrance = penetrance if isinstance(penetrance, dict) else {}
    onset = condition.get("onset", {})
    onset = onset if isinstance(onset, dict) else {}

    axis_assertions: list[dict[str, str]] = []
    axis_seen: set[str] = set()
    for axis in (inheritance, penetrance, onset):
        for raw_assertion in axis.get("assertions", []) or []:
            if not isinstance(raw_assertion, dict):
                continue
            hpo_assertion = {
                "hpo_id": _clean(raw_assertion.get("hpo_id")),
                "frequency": _clean(raw_assertion.get("frequency")),
                "evidence": _clean(raw_assertion.get("evidence")),
                "reference": _clean(raw_assertion.get("reference")),
            }
            key = json.dumps(hpo_assertion, sort_keys=True, separators=(",", ":"))
            if key not in axis_seen:
                axis_seen.add(key)
                axis_assertions.append(hpo_assertion)

    return {
        "hpo_match_status": "matched_gene_and_condition",
        "hpo_disease_id": _clean(condition_key),
        "mondo_id": mondo_id,
        "disease_scope": _clean(scope.get("category")),
        "priva_scope": "include",
        "scope_review_status": _clean(scope.get("review_status")),
        "hpo_inheritance_modes": hpo_inheritance_modes,
        "penetrance_hpo_ids": [
            _clean(item.get("hpo_id"))
            for item in penetrance.get("assertions", []) or []
            if isinstance(item, dict) and _clean(item.get("hpo_id"))
        ],
        "onset_hpo_ids": [
            _clean(item.get("hpo_id"))
            for item in onset.get("assertions", []) or []
            if isinstance(item, dict) and _clean(item.get("hpo_id"))
        ],
        "hpo_assertions": axis_assertions,
        "hpo_assertion_count": condition.get("hpo_assertion_count", 0),
    }


def condition_cache_mechanism_entries(
    gene_record: dict[str, Any],
) -> list[dict[str, Any]]:
    """Flatten one integrated gene record for hub summaries and audits.

    Unlike automatic assertion selection, this audit view retains mechanisms
    from review/excluded conditions and evidence that could not be joined to an
    HPO condition. Scope and condition-link status remain explicit on every
    row. Exact GoFCards alleles are emitted at ``variant_level`` and never
    converted into gene-level mechanism evidence.
    """
    entries: list[dict[str, Any]] = []
    seen: set[str] = set()

    def append_entry(entry: dict[str, Any]) -> None:
        normalized = {
            **entry,
            "level": _clean(entry.get("level")),
            "source": _clean(entry.get("source")),
            "mechanism": _clean(entry.get("mechanism")).upper(),
            "pmids": list(entry.get("pmids", []) or []),
        }
        if normalized["mechanism"] not in CANONICAL_MECHANISMS:
            return
        identity = json.dumps(normalized, sort_keys=True, separators=(",", ":"))
        if identity not in seen:
            seen.add(identity)
            entries.append(normalized)

    conditions = gene_record.get("conditions", {})
    if isinstance(conditions, dict):
        for condition_key, condition in sorted(conditions.items()):
            if not isinstance(condition, dict):
                continue
            scope = condition.get("priva_scope", {})
            scope = scope if isinstance(scope, dict) else {}
            condition_common = {
                "source_condition_id": _clean(condition_key),
                "disease": _clean(condition.get("label")),
                "disease_scope": _clean(scope.get("category")),
                "priva_scope": _clean(scope.get("decision")),
                "scope_review_status": _clean(scope.get("review_status")),
                "condition_link_status": "exact",
            }
            mechanisms = condition.get("pathogenic_mechanisms", {})
            if not isinstance(mechanisms, dict):
                continue
            for mechanism, block in sorted(mechanisms.items()):
                mechanism = _clean(mechanism).upper()
                if mechanism not in CANONICAL_MECHANISMS or not isinstance(block, dict):
                    continue
                for evidence in block.get("evidence", []) or []:
                    if not isinstance(evidence, dict):
                        continue
                    source = _clean(evidence.get("source"))
                    if source == "GoFCards_exact+ClinVar_VCV":
                        continue
                    append_entry(
                        {
                            **condition_common,
                            "level": "gene_level",
                            "source": source,
                            "source_record_id": _clean(
                                evidence.get("source_record_id")
                            ),
                            "mechanism": mechanism,
                            "mechanism_raw": _clean(evidence.get("mechanism_raw")),
                            "allelic_requirement": _clean(
                                evidence.get("allelic_requirement")
                            ),
                            "confidence": _clean(
                                evidence.get("mechanism_confidence")
                            ),
                            "mechanism_confidence": _clean(
                                evidence.get("mechanism_confidence")
                            ),
                            "disease_confidence": _clean(
                                evidence.get("disease_confidence")
                            ),
                            "pmids": list(evidence.get("pmids", []) or []),
                        }
                    )
                for variant_key, variant in sorted(
                    (block.get("variants", {}) or {}).items()
                ):
                    if not isinstance(variant, dict):
                        continue
                    links = [
                        link
                        for link in variant.get("clinvar_links", []) or []
                        if isinstance(link, dict)
                    ]
                    append_entry(
                        {
                            **condition_common,
                            "level": "variant_level",
                            "source": "GoFCards_exact+ClinVar_VCV",
                            "source_record_id": _clean(variant_key),
                            "mechanism": mechanism,
                            "allelic_requirement": "",
                            "confidence": "exact_variant",
                            "mechanism_confidence": "exact_variant",
                            "disease_confidence": "ClinVar_germline_assertion",
                            "pmids": list(variant.get("pmids", []) or []),
                            "vcv_accessions": [
                                _clean(link.get("vcv_accession"))
                                for link in links
                                if _clean(link.get("vcv_accession"))
                            ],
                        }
                    )

    unmapped = gene_record.get("unmapped_evidence", {})
    unmapped = unmapped if isinstance(unmapped, dict) else {}
    for evidence in unmapped.get("mechanisms", []) or []:
        if not isinstance(evidence, dict):
            continue
        source_scope = evidence.get("source_scope", {})
        source_scope = source_scope if isinstance(source_scope, dict) else {}
        append_entry(
            {
                "level": "gene_level",
                "source": _clean(evidence.get("source")),
                "source_record_id": _clean(evidence.get("source_record_id")),
                "source_condition_id": next(
                    iter(evidence.get("condition_identifiers", []) or []), ""
                ),
                "disease": _clean(evidence.get("condition_label")),
                "mechanism": _clean(evidence.get("mechanism")),
                "mechanism_raw": _clean(evidence.get("mechanism_raw")),
                "allelic_requirement": _clean(evidence.get("allelic_requirement")),
                "confidence": _clean(evidence.get("mechanism_confidence")),
                "mechanism_confidence": _clean(
                    evidence.get("mechanism_confidence")
                ),
                "disease_confidence": _clean(evidence.get("disease_confidence")),
                "disease_scope": _clean(source_scope.get("category")),
                "priva_scope": _clean(source_scope.get("decision")),
                "scope_review_status": _clean(source_scope.get("review_status")),
                "condition_link_status": "unresolved",
                "unmapped_reason": _clean(evidence.get("reason")),
                "pmids": list(evidence.get("pmids", []) or []),
            }
        )
    variants = unmapped.get("variants", {})
    variants = variants if isinstance(variants, dict) else {}
    for variant_key, variant in sorted(variants.items()):
        if not isinstance(variant, dict):
            continue
        append_entry(
            {
                "level": "variant_level",
                "source": "GoFCards",
                "source_record_id": _clean(variant_key),
                "source_condition_id": "",
                # The disease name comes from the ClinVar records linked to
                # this variant. GoFCards' own free-text disease field is not
                # carried into the condition cache: it has no identifier, so
                # it can name a condition but never join to one.
                "disease": ";".join(
                    dict.fromkeys(
                        name
                        for link in variant.get("clinvar_links", []) or []
                        if isinstance(link, dict)
                        for name in link.get("condition_names", []) or []
                        if _clean(name)
                    )
                ),
                "mechanism": _clean(variant.get("mechanism")),
                "allelic_requirement": "",
                "confidence": "exact_variant_condition_unresolved",
                "mechanism_confidence": "exact_variant",
                "disease_confidence": "",
                "condition_link_status": "unresolved",
                "unmapped_reason": _clean(
                    (variant.get("condition_link", {}) or {}).get("reason")
                ),
                "pmids": list(variant.get("pmids", []) or []),
            }
        )
    return entries


def condition_cache_mechanism_assertions(
    gene_record: dict[str, Any],
) -> list[dict[str, Any]]:
    """Translate one integrated gene record to PriVA's assertion contract.

    The integrated cache has already performed the only permitted join:
    gene + exact condition identifier. This adapter therefore does not repeat
    disease matching. It exposes only conditions whose effective PriVA scope is
    ``include`` and retains the condition's inheritance, penetrance, onset, and
    assertion provenance beside every mechanism record.

    The cache also stores synthetic ``GoFCards_exact+ClinVar_VCV`` evidence in a
    condition block when a particular allele has an exact ClinVar condition
    link. That source is intentionally excluded here: it becomes an assertion
    only when the query allele matches the nested cache variant. Otherwise one
    exact GOF allele would incorrectly create gene-wide GOF history.
    """
    assertions: list[dict[str, Any]] = []
    conditions = gene_record.get("conditions", {})
    if not isinstance(conditions, dict):
        return assertions

    for condition_key, condition in sorted(conditions.items()):
        if not isinstance(condition, dict):
            continue
        common = condition_cache_context(condition_key, condition)
        if not common:
            continue

        mechanism_blocks = condition.get("pathogenic_mechanisms", {})
        if not isinstance(mechanism_blocks, dict):
            continue
        for mechanism, block in sorted(mechanism_blocks.items()):
            mechanism = _clean(mechanism).upper()
            if mechanism not in CANONICAL_MECHANISMS or not isinstance(block, dict):
                continue
            for evidence in block.get("evidence", []) or []:
                if not isinstance(evidence, dict):
                    continue
                source = _clean(evidence.get("source"))
                if source not in CONDITION_MECHANISM_SOURCES:
                    continue
                evidence_identifiers = evidence.get("condition_identifiers", [])
                evidence_identifiers = (
                    evidence_identifiers if isinstance(evidence_identifiers, list) else []
                )
                source_condition_id = next(
                    (
                        _clean(identifier)
                        for identifier in evidence_identifiers
                        if _clean(identifier) and not _clean(identifier).upper().startswith("MONDO:")
                    ),
                    _clean(condition_key),
                )
                requirement = _clean(evidence.get("allelic_requirement"))
                requirements = [requirement] if requirement else sorted(
                    _hpo_allelic_requirements(
                        ";".join(common["hpo_inheritance_modes"])
                    )
                )
                if not requirements:
                    requirements = [""]

                for effective_requirement in requirements:
                    assertions.append(
                        {
                            **common,
                            "source": source,
                            "source_record_id": _clean(evidence.get("source_record_id")),
                            "source_condition_id": source_condition_id,
                            "disease": _clean(evidence.get("condition_label"))
                            or _clean(condition.get("label")),
                            "mechanism": mechanism,
                            "mechanism_raw": _clean(evidence.get("mechanism_raw")),
                            "allelic_requirement": effective_requirement,
                            "confidence": _clean(evidence.get("mechanism_confidence")),
                            "mechanism_confidence": _clean(
                                evidence.get("mechanism_confidence")
                            ),
                            "disease_confidence": _clean(
                                evidence.get("disease_confidence")
                            ),
                            "pmids": list(evidence.get("pmids", []) or []),
                        }
                    )
    return _deduplicate_by(assertions, ASSERTION_IDENTITY_FIELDS)


def enrich_condition_mechanism_assertion(
    assertion: dict[str, Any],
    *,
    gene_symbol: Any,
    hpo_condition_index: dict[tuple[str, str], dict[str, Any]],
) -> list[dict[str, Any]]:
    """Bind one mechanism assertion to HPO for the same gene and disease.

    Native disease identity is tried before MONDO identity. This avoids linking
    through names and prevents another disease of the same gene from donating
    inheritance, penetrance, or onset. Automatic transfer requires an effective
    ``priva_scope=include``; review, excluded, and still-unscoped diseases remain
    audit-only.

    G2P allelic requirements remain authoritative when present. Orphadata does
    not encode allelic requirements, so an unambiguous HPO inheritance term can
    supply a compact dominant, recessive, or mitochondrial requirement. If HPO
    legitimately records more than one mode, separate records are returned so
    no mode is hidden in a lossy combined string.
    """
    record = dict(assertion)
    source_scope = _clean(record.get("priva_scope")).lower()
    if source_scope in {"review", "exclude"}:
        return []

    symbol_key = _clean(gene_symbol).upper()
    hpo_record: dict[str, Any] | None = None
    for condition_id in (
        record.get("source_condition_id"),
        record.get("mondo_id"),
    ):
        condition_key = _clean(condition_id).upper()
        if condition_key:
            hpo_record = hpo_condition_index.get((symbol_key, condition_key))
        if hpo_record is not None:
            break

    if hpo_record is None:
        # A source assertion without a registry-supported inherited-disease
        # scope is useful for manual review, but it cannot safely influence an
        # automated germline classification. Explicitly included G2P diseases
        # remain usable even when HPO has no phenotype rows for that condition.
        if source_scope != "include":
            return []
        record.update(
            {
                "hpo_match_status": "no_matching_gene_condition_hpo_record",
                "hpo_inheritance_modes": [],
                "penetrance_hpo_ids": [],
                "onset_hpo_ids": [],
                "hpo_assertions": [],
            }
        )
        return [record]
    if _clean(hpo_record.get("priva_scope")).lower() in {"review", "exclude"}:
        return []

    record.update(
        {
            "hpo_match_status": "matched_gene_and_condition",
            "hpo_disease_id": _clean(hpo_record.get("disease_id")),
            "mondo_id": _clean(record.get("mondo_id"))
            or _clean(hpo_record.get("mondo_id")),
            "hpo_inheritance_modes": list(
                hpo_record.get("inheritance_modes", [])
            ),
            "penetrance_hpo_ids": list(
                hpo_record.get("penetrance_hpo_ids", [])
            ),
            "onset_hpo_ids": list(hpo_record.get("onset_hpo_ids", [])),
            "hpo_assertions": list(hpo_record.get("hpo_assertions", [])),
            "disease_scope": _clean(hpo_record.get("disease_scope"))
            or _clean(record.get("disease_scope")),
            "priva_scope": _clean(hpo_record.get("priva_scope"))
            or _clean(record.get("priva_scope")),
            "scope_review_status": _clean(
                hpo_record.get("scope_review_status")
            )
            or _clean(record.get("scope_review_status")),
        }
    )
    if _clean(record.get("priva_scope")).lower() != "include":
        return []
    if _clean(record.get("allelic_requirement")):
        return [record]

    requirements = sorted(
        _hpo_allelic_requirements(";".join(record["hpo_inheritance_modes"]))
    )
    if not requirements:
        return [record]
    return [
        {**record, "allelic_requirement": requirement}
        for requirement in requirements
    ]


# The reverse of HPO_CACHE_INHERITANCE_LABELS, so a mode is recognized whether
# it arrives as the cache's key or as the label substituted for display.
_INHERITANCE_LABEL_TO_KEY = {
    label.lower(): key for key, label in HPO_CACHE_INHERITANCE_LABELS.items()
}
# Downstream treats these three the same way: one allele in a carrier is
# enough for the condition to appear.
DOMINANT_LIKE_INHERITANCE = {"dominant", "y_linked", "mitochondrial"}
# Germline, but not single-gene. Delivered rather than discarded, so that
# benign-supporting criteria are not assigned easily against them.
NON_MENDELIAN_INHERITANCE = {
    "non_mendelian",
    "polygenic",
    "digenic",
    "oligogenic",
}


def normalize_inheritance(
    allelic_requirement: Any = "",
    hpo_inheritance_modes: Any = (),
) -> tuple[str, bool]:
    """Reduce the two source vocabularies to one value plus an X-linked flag.

    Two sources describe the same fact in different words: G2P states an
    allelic requirement (``biallelic_autosomal``, ``monoallelic_X_hemizygous``)
    and HPO states an inheritance mode (``autosomal_recessive``,
    ``x_linked_dominant``). Both fold to the same answer.

    The delivered value is ``recessive`` or ``dominant`` wherever that question
    has an answer, because that is what the criteria reason about. Being on the
    X chromosome is a separate fact and is returned separately: folding it into
    the value would erase it, and a hemizygous male affected by one allele is
    still the recessive pattern.

    ``y_linked`` and ``mitochondrial`` are delivered as themselves rather than
    forced into dominant; downstream treats them the same as dominant.
    ``non_mendelian``, ``polygenic``, ``digenic`` and ``oligogenic`` are also
    delivered as themselves. They are not discarded -- they exist so that
    benign-supporting criteria are not assigned easily against them.
    """
    requirement = _clean(allelic_requirement).lower()
    # Modes arrive in either form: the cache's own snake_case key, or the
    # human-readable HPO label that condition_cache_context substitutes. Both
    # name the same mode, so both are accepted here rather than requiring every
    # caller to know which one it is holding.
    modes = [
        _INHERITANCE_LABEL_TO_KEY.get(_clean(mode).lower(), _clean(mode).lower())
        for mode in (
            [hpo_inheritance_modes]
            if isinstance(hpo_inheritance_modes, str)
            else list(hpo_inheritance_modes or [])
        )
        if _clean(mode)
    ]
    x_linked = requirement.startswith("monoallelic_x") or any(
        mode.startswith("x_linked") for mode in modes
    )

    # A non-Mendelian mode is reported only when no Mendelian mode accompanies
    # it. Both surviving genes in the current cache carry the two together, and
    # the Mendelian one is the actionable statement.
    mendelian_modes = [mode for mode in modes if mode not in NON_MENDELIAN_INHERITANCE]
    if not mendelian_modes:
        for mode in modes:
            if mode in NON_MENDELIAN_INHERITANCE:
                return mode, x_linked

    if requirement.startswith("biallelic_") or requirement in {
        "monoallelic_x",
        "monoallelic_x_hemizygous",
    }:
        return "recessive", x_linked
    if requirement == "mitochondrial":
        return "mitochondrial", False
    if requirement.startswith("monoallelic_y"):
        return "y_linked", False
    if requirement.startswith("monoallelic_"):
        return "dominant", x_linked

    for mode in mendelian_modes:
        if mode in {"autosomal_recessive", "x_linked_recessive", "pseudoautosomal_recessive"}:
            return "recessive", x_linked
        # A bare "x_linked" with no recessive or dominant qualifier reads as
        # X-linked recessive.
        if mode == "x_linked":
            return "recessive", True
        if mode in {
            "autosomal_dominant",
            "x_linked_dominant",
            "autosomal_dominant_maternal_imprinting",
        }:
            return "dominant", x_linked
        if mode == "mitochondrial":
            return "mitochondrial", False
        if mode == "y_linked":
            return "y_linked", False
    return "", x_linked


def gene_inheritance_consensus(gene_record: dict[str, Any]) -> tuple[str, bool, int]:
    """Return the one inheritance every germline condition of a gene agrees on.

    This is the fallback for a variant whose mechanism reaches no history --
    because no condition of that gene states a mechanism at all. The inheritance
    is still knowable: if every germline-included condition of the gene says
    dominant, then whatever this variant does, one allele is what the disease
    requires. Leaving it empty would discard a fact we hold.

    Unanimity is what makes it safe. A gene carrying both a dominant and a
    recessive condition is deliberately given nothing here, because there the
    question genuinely has two answers and only the per-history match can choose
    between them. Measured on the current cache: 2,772 genes are unanimously
    recessive, 1,453 unanimously dominant, and 501 are mixed.

    Conditions that are not germline-inherited disease never contribute, so a
    review-scoped or excluded condition cannot donate its inheritance here.

    The third return value is how many distinct inheritances the gene's germline
    conditions stated, which separates the two ways this can come back empty:
    ``0`` means nothing was stated at all and the caller may fall back further,
    while ``2`` or more means the gene genuinely disagrees with itself and no
    fallback is legitimate.
    """
    values: set[tuple[str, bool]] = set()
    for condition in (gene_record.get("conditions") or {}).values():
        if (condition.get("priva_scope") or {}).get("decision") != "include":
            continue
        inheritance, x_linked = normalize_inheritance(
            "", (condition.get("inheritance") or {}).get("modes") or []
        )
        if inheritance:
            values.add((inheritance, x_linked))
    distinct = len({inheritance for inheritance, _x in values})
    if distinct != 1:
        return "", False, distinct
    inheritance, x_linked = next(iter(values))
    return inheritance, x_linked, 1


def gene_inheritance_from_constraint(
    symbol: str,
    *,
    clingen: dict[str, dict[str, Any]],
    loeuf: dict[str, float],
) -> str:
    """Last resort for a gene HPO says nothing about: read the constraint data.

    Reached only when no germline condition of the gene states an inheritance at
    all -- 798 genes, for 755 of which HPO holds no inheritance annotation
    anywhere. Rather than deliver nothing, the inheritance is inferred from
    whether the gene tolerates losing one copy:

        dominant   ClinGen haploinsufficiency score 3, or LOEUF below 0.35.
                   Either is a statement that one lost copy already causes
                   disease, which is the dominant pattern.
        recessive  otherwise, as the default. Most disease genes are recessive,
                   and a gene with no haploinsufficiency signal is far more
                   likely to need both copies disabled.

    The two signals barely overlap -- of these genes 196 have only the LOEUF
    signal, 6 only the ClinGen one, 9 both -- so requiring both would reduce
    the rule to 9 genes. Either is therefore enough.

    No mechanism accompanies this. We do not know what a variant must do to
    cause disease here, only how many copies must be affected.
    """
    score = clingen.get(symbol, {}).get("haploinsufficiency_score")
    try:
        haploinsufficient = int(str(score)) == 3
    except (TypeError, ValueError):
        haploinsufficient = False
    constraint = loeuf.get(symbol, float("nan"))
    constrained = (
        isinstance(constraint, (int, float))
        and not math.isnan(constraint)
        and constraint < 0.35
    )
    return "dominant" if haploinsufficient or constrained else "recessive"


def normalize_penetrance(penetrance_hpo_ids: Any = ()) -> str:
    """Reduce the penetrance HPO terms to complete, incomplete, or unknown.

    Moderate and low penetrance are forms of incomplete penetrance: in each the
    condition fails to appear in some carriers, which is the only distinction
    the criteria act on.
    """
    ids = {
        _clean(value).upper()
        for value in (penetrance_hpo_ids or [])
        if _clean(value)
    }
    if ids & {"HP:0003829", "HP:4000159", "HP:4000160"}:
        return "incomplete"
    if "HP:0034950" in ids:
        return "complete"
    return "unknown"


def _compact_inheritance(allelic_requirement: Any) -> str:
    """Inheritance value alone, for the compact per-row mechanism tags."""
    return normalize_inheritance(allelic_requirement)[0]


def _mechanism_profile_tag(history: dict[str, Any]) -> str:
    """One compact `<inheritance>_<MECHANISM>` tag for a selected history.

    The inheritance is read straight off the history, which
    select_condition_histories_for_variant already normalized. Re-deriving it
    here from the raw allelic requirement would be a second copy of that rule.
    """
    inheritance = _clean(history.get("inheritance"))
    mechanism = _clean(history.get("mechanism")).upper() or "UNRESOLVED"
    if mechanism == "UNRESOLVED":
        return inheritance or "uncertain"
    mechanism = "DN" if mechanism == "DOMINANT_NEGATIVE" else mechanism
    return f"{inheritance}_{mechanism}" if inheritance else mechanism


def _compact_profile_tags(assertions: list[dict[str, Any]]) -> list[str]:
    tags = list(dict.fromkeys(_mechanism_profile_tag(assertion) for assertion in assertions))
    inheritance_with_mechanism = {
        tag.split("_", 1)[0]
        for tag in tags
        if "_" in tag and tag.split("_", 1)[0] in {"recessive", "dominant"}
    }
    qualified_mechanisms = {
        tag.split("_", 1)[1]
        for tag in tags
        if "_" in tag and tag.split("_", 1)[0] in {"recessive", "dominant"}
    }
    return [
        tag
        for tag in tags
        if tag not in inheritance_with_mechanism
        and tag not in qualified_mechanisms
        and tag != "uncertain"
    ]


def select_condition_histories_for_variant(
    assertions: list[dict[str, Any]],
    *,
    variant_effect: str,
) -> list[dict[str, Any]]:
    """Step 3: keep the histories this variant's mechanism reaches, and carry
    back the condition, the inheritance and the penetrance.

    A curated mechanism elsewhere in a gene is background history, not evidence
    that every variant acts through that mechanism. Selection therefore happens
    before inheritance or penetrance can influence any criterion:

    * an exact known GOF or DN allele selects that mechanism's histories only;
    * a high-confidence predicted LOF allele selects LOF histories; and
    * an unresolved allele keeps every history, marked uncertain.

    Exact curated mechanisms are exclusive: a consequence-based loss-of-function
    prediction cannot reintroduce a different history after an exact match.

    Each surviving history is reduced to the three facts the chain exists to
    deliver, plus the mechanism that selected it. Nothing else travels, because
    nothing else is read.
    """
    effect = _clean(variant_effect)
    exact_mechanisms = _exact_mechanisms_from_effect(effect)
    if exact_mechanisms:
        allowed = exact_mechanisms
    elif effect == "predicted_LOF_high_confidence":
        allowed = {"LOF"}
    else:
        allowed = CANONICAL_MECHANISMS

    histories: list[dict[str, Any]] = []
    for assertion in assertions:
        mechanism = _clean(assertion.get("mechanism")).upper()
        if mechanism not in allowed:
            continue
        inheritance, x_linked = normalize_inheritance(
            assertion.get("allelic_requirement"),
            assertion.get("hpo_inheritance_modes"),
        )
        histories.append(
            {
                "mechanism": mechanism,
                "condition": _clean(
                    assertion.get("hpo_disease_id")
                    or assertion.get("source_condition_id")
                ),
                "inheritance": inheritance,
                "x_linked": x_linked,
                "penetrance": normalize_penetrance(
                    assertion.get("penetrance_hpo_ids")
                ),
            }
        )
    return histories


HISTORY_IDENTITY_FIELDS = ("mechanism", "condition", "inheritance", "x_linked", "penetrance")
ASSERTION_IDENTITY_FIELDS = (
    "source",
    "source_record_id",
    "source_condition_id",
    "mondo_id",
    "disease",
    "mechanism",
    "allelic_requirement",
    "confidence",
)


def _deduplicate_by(
    records: list[dict[str, Any]],
    fields: tuple[str, ...],
) -> list[dict[str, Any]]:
    """Drop repeats, comparing only the named fields.

    Used twice at different stages: on the gene's raw assertions, where two
    sources can state the same thing, and on the selected histories, where the
    identity is only the facts the chain delivers. One rule, two field lists,
    rather than two functions that would drift apart.
    """
    seen: set[tuple[Any, ...]] = set()
    output: list[dict[str, Any]] = []
    for record in records:
        key = tuple(
            bool(record.get(field))
            if isinstance(record.get(field), bool)
            else _clean(record.get(field))
            for field in fields
        )
        if key in seen:
            continue
        seen.add(key)
        output.append(record)
    return output


def _classify_variant_applicability(
    histories: list[dict[str, Any]],
    effect: str,
) -> dict[str, Any]:
    """Split the selected histories into established and merely possible.

    Step 3 has already removed the histories this variant's mechanism cannot
    reach, so nothing here is incompatible. What remains is one distinction the
    criteria act on: whether the variant's mechanism is *established* for that
    history (an exact curated allele, or a high-confidence predicted loss of
    function) or only *possible* because the variant's own effect is unresolved.
    """
    exact_mechanisms = _exact_mechanisms_from_effect(effect)
    applicable: list[str] = []
    uncertain: list[str] = []
    for history in histories:
        mechanism = _clean(history.get("mechanism")).upper()
        established = (
            mechanism in exact_mechanisms
            or (mechanism == "LOF" and effect == "predicted_LOF_high_confidence")
        )
        tag = _mechanism_profile_tag(history)
        target = applicable if established else uncertain
        if tag not in target:
            target.append(tag)
    return {
        "plausible": ";".join(_compact_profile_tags(histories)),
        "applicable": ";".join(applicable),
        "uncertain": ";".join(uncertain),
    }


def annotate_gene_mechanism_categories(
    df: pd.DataFrame,
    *,
    condition_cache: str | Path = DEFAULT_HPO_CONDITION_MECHANISM_CACHE,
    symbol_col: str = "SYMBOL",
    output_col: str = "var_plausible_patho_mechs",
    use_hgnc_package: bool = False,
    hpo_collapsed: str | Path = DEFAULT_HPO_COLLAPSED,
    clingen_dosage: str | Path = DEFAULT_CLINGEN_DOSAGE,
    loeuf_table: str | Path = DEFAULT_LOEUF_TABLE,
    hgnc_table: str | Path = DEFAULT_HGNC_TABLE,
) -> pd.DataFrame:
    """Step 4: run the chain over every row and attach the result.

    The only evidence sources are the variant's own annotation and the HPO
    condition cache. Gene-wide signals -- a ClinVar pathogenic history, a
    constrained LOEUF, a high average AlphaMissense score, a ClinGen dosage
    call -- are deliberately absent: none of them says which condition a
    variant acts on or by what mechanism, so none can create a history.
    """
    if symbol_col not in df.columns:
        raise KeyError(f"missing symbol column: {symbol_col}")

    hub = GeneMechanismHub(
        condition_cache=condition_cache,
        hpo_collapsed=hpo_collapsed,
        clingen_dosage=clingen_dosage,
        loeuf_table=loeuf_table,
        hgnc_table=hgnc_table,
        use_hgnc_package=use_hgnc_package,
    )
    out = df.copy()
    plausible_mechanism_values: list[str] = []
    variant_outputs: dict[str, list[str]] = {
        column: [] for column in VARIANT_MECHANISM_OUTPUT_COLUMNS
    }
    assertion_cache: dict[str, list[dict[str, Any]]] = {}
    for _, row in out.iterrows():
        gene = row[symbol_col]
        symbol = hub.resolve_symbol(gene)
        if symbol not in assertion_cache:
            assertion_cache[symbol] = hub.condition_mechanism_assertions(gene)
        assertions = list(assertion_cache[symbol])

        # STEP 1: what mechanism does this variant plausibly act by?
        effect_call = infer_query_variant_effect(row)
        # STEPS 2 and 3: the gene's germline condition histories this
        # mechanism reaches, reduced to condition, inheritance and penetrance.
        histories = _deduplicate_by(
            select_condition_histories_for_variant(
                assertions,
                variant_effect=effect_call["variant_effect"],
            ),
            HISTORY_IDENTITY_FIELDS,
        )
        applicability = _classify_variant_applicability(
            histories,
            effect_call["variant_effect"],
        )

        # STEP 4: attach at variant-transcript resolution.
        plausible_mechanism_values.append(applicability["plausible"])
        variant_outputs["variant_effect"].append(effect_call["variant_effect"])
        variant_outputs["variant_lof_score"].append(effect_call["variant_lof_score"])
        variant_outputs["variant_gof_score"].append(effect_call["variant_gof_score"])
        variant_outputs["variant_dn_score"].append(effect_call["variant_dn_score"])
        variant_outputs["variant_mechanism_exclusive"].append(
            effect_call["variant_mechanism_exclusive"]
        )
        variant_outputs["variant_exact_mechanisms"].append(
            effect_call["variant_exact_mechanisms"]
        )
        variant_outputs["variant_mechanism_applicable"].append(applicability["applicable"])
        variant_outputs["variant_mechanism_uncertain"].append(applicability["uncertain"])
        # The three facts the chain exists to deliver, at this variant's own
        # resolution. One entry per selected history, in a stable order.
        matched_inheritance = list(
            dict.fromkeys(h["inheritance"] for h in histories if h["inheritance"])
        )
        matched_x_linked = any(h["x_linked"] for h in histories)
        if matched_inheritance:
            basis = "matched_history"
        else:
            # No history stated an inheritance, either because the variant's
            # mechanism reached none or because none of this gene's conditions
            # records a mechanism at all. Two fallbacks remain, in order of how
            # much they are grounded in the gene's own disease record.
            consensus, consensus_x, stated = gene_inheritance_consensus(
                hub._load_condition_cache().get(symbol, {})
            )
            if consensus:
                matched_inheritance = [consensus]
                matched_x_linked = consensus_x
                basis = "gene_consensus"
            elif stated == 0:
                # HPO says nothing about this gene's inheritance at all, so the
                # constraint data decides. A gene that disagrees with itself
                # (stated > 1) deliberately falls through to nothing.
                matched_inheritance = [
                    gene_inheritance_from_constraint(
                        symbol, clingen=hub._load_clingen(), loeuf=hub._load_loeuf()
                    )
                ]
                basis = "gene_constraint"
            else:
                matched_inheritance = []
                basis = ""

        variant_outputs["variant_condition_ids"].append(
            ";".join(dict.fromkeys(h["condition"] for h in histories if h["condition"]))
        )
        # The same facts kept together rather than as three parallel lists.
        # De-duplicating each list separately makes them different lengths, so
        # a reader cannot tell which inheritance belongs to which condition;
        # here each entry is one whole history and the pairing survives.
        #
        #   <condition>|<mechanism>|<inheritance>|<penetrance>
        #
        # Empty when no history was reached at all: the basis column then says
        # whether the inheritance came from the gene's consensus or from its
        # constraint data, neither of which belongs to a named condition.
        variant_outputs["variant_condition_histories"].append(
            ";".join(
                dict.fromkeys(
                    "|".join(
                        (
                            h["condition"],
                            "DN" if h["mechanism"] == "DOMINANT_NEGATIVE"
                            else h["mechanism"],
                            h["inheritance"] or "unknown",
                            h["penetrance"],
                        )
                    )
                    for h in histories
                    if h["condition"]
                )
            )
        )
        variant_outputs["variant_inheritance"].append(";".join(matched_inheritance))
        variant_outputs["variant_inheritance_basis"].append(basis)
        variant_outputs["variant_x_linked"].append(
            "true" if matched_x_linked else "false"
        )
        variant_outputs["variant_penetrance"].append(
            ";".join(
                dict.fromkeys(
                    h["penetrance"] for h in histories if h["penetrance"] != "unknown"
                )
            )
            or "unknown"
        )

    out[output_col] = plausible_mechanism_values
    for column, values in variant_outputs.items():
        out[column] = values
    return out


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
