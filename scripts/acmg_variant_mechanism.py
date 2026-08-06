#!/usr/bin/env python3
"""Shared foundation: the hub's mechanism output, made ACMG-usable.

No ACMG criterion lives here. This module owns the two jobs every
mechanism-reading criterion depends on, so that all of them read one
interpretation rather than each rolling its own.

ATTACH AN EXACT CURATED ALLELE MATCH
    annotate_exact_nonlof_variants looks the query allele up in the non-LOF
    cache -- by HGVSp, then HGVSc, then the normalized genomic allele -- and
    writes the resulting gain-of-function or dominant-negative score of 2.
    It runs before the mechanism hub, because the hub's precedence ladder
    needs to know whether a curated call exists for this exact allele.

TURN THE HUB'S COLUMNS INTO BOOLEAN MASKS
    _variant_mechanism_masks takes two columns and returns fifteen masks in
    two families: what the GENE's condition history says, split by
    inheritance, and what THIS VARIANT does. There is deliberately no third
    family conjoining them -- whether the two line up is each criterion's own
    judgement, and pre-computing one verdict here would decide it for all of
    them at once.

    _variant_condition_masks parses the hub's categorical inheritance and
    penetrance columns. It is the one adapter used by every criterion that
    previously received the six gene-level arrays, so exact semicolon-token
    matching and the incomplete-penetrance gate cannot drift between criteria.

Read by PVS1, BS1, BS2, BS4, BP2/PM3, PP1 and the ranking step.
"""

import logging
import json
import re
import pandas as pd
import numpy as np

from gene_mechanism_hub import DEFAULT_MECHANISM_JSON, GeneMechanismHub


logger = logging.getLogger(__name__)


EXACT_NONLOF_OUTPUT_COLS = (
    "variant_lof_score",
    "variant_gof_score",
    "variant_dn_score",
    "variant_mechanism_exclusive",
    "variant_exact_mechanisms",
    "variant_exact_mechanism_evidence",
    "variant_mechanism_match_route",
    "variant_gof_tag",
    "gofcards_accession_id",
    "gofcards_variant_id",
    "gofcards_match_route",
)


def _clean_text(value):
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"", "nan", "none", "<na>"} else text


def _row_assembly_for_gof_lookup(row):
    for col in ("assembly", "genome_build", "reference_build", "build"):
        value = _clean_text(row.get(col))
        if value:
            return value
    return "hg19"


def _merge_semicolon_values(*values):
    seen = set()
    merged = []
    for value in values:
        for token in _clean_text(value).replace(",", ";").split(";"):
            token = token.strip()
            if token and token not in seen:
                seen.add(token)
                merged.append(token)
    return ";".join(merged)


def _ensure_exact_nonlof_columns(df):
    score_columns = ("variant_lof_score", "variant_gof_score", "variant_dn_score")
    for col in score_columns:
        if col not in df.columns:
            df[col] = 0
        else:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)
    if "variant_mechanism_exclusive" not in df.columns:
        df["variant_mechanism_exclusive"] = False
    else:
        df["variant_mechanism_exclusive"] = (
            df["variant_mechanism_exclusive"]
            .fillna(False)
            .astype(str)
            .str.lower()
            .isin({"1", "true", "yes"})
        )
    for col in EXACT_NONLOF_OUTPUT_COLS:
        if col in score_columns or col == "variant_mechanism_exclusive":
            continue
        if col not in df.columns:
            df[col] = ""
        else:
            df[col] = df[col].fillna("").astype(str)
    if "gofcards_accession" in df.columns:
        df["gofcards_accession_id"] = [
            _merge_semicolon_values(current, legacy)
            for current, legacy in zip(df["gofcards_accession_id"], df["gofcards_accession"])
        ]


def _record_exact_nonlof_match(df, indices, match):
    if not match.get("mechanism_exclusive"):
        return
    route = _clean_text(match.get("match_route"))
    mechanisms = list(match.get("exclusive_mechanisms", []) or [])
    scores = match.get("mechanism_scores", {}) or {}
    evidence = json.dumps(
        match.get("matches", []),
        sort_keys=True,
        separators=(",", ":"),
    )
    score_columns = {
        "LOF": "variant_lof_score",
        "GOF": "variant_gof_score",
        "DOMINANT_NEGATIVE": "variant_dn_score",
    }
    for mechanism, column in score_columns.items():
        score = int(scores.get(mechanism, 0) or 0)
        df.loc[indices, column] = np.maximum(
            pd.to_numeric(df.loc[indices, column], errors="coerce").fillna(0),
            score,
        )
    df.loc[indices, "variant_mechanism_exclusive"] = True
    df.loc[indices, "variant_exact_mechanisms"] = [
        _merge_semicolon_values(value, ";".join(mechanisms))
        for value in df.loc[indices, "variant_exact_mechanisms"]
    ]
    df.loc[indices, "variant_exact_mechanism_evidence"] = evidence
    if route:
        df.loc[indices, "variant_mechanism_match_route"] = [
            _merge_semicolon_values(value, route)
            for value in df.loc[indices, "variant_mechanism_match_route"]
        ]

    if "GOF" not in mechanisms:
        return
    accession_id = _clean_text(match.get("gofcards_accession_id"))
    variant_id = _clean_text(match.get("gofcards_variant_id"))
    df.loc[indices, "variant_gof_tag"] = [
        _merge_semicolon_values(value, "GOF") for value in df.loc[indices, "variant_gof_tag"]
    ]
    if accession_id:
        df.loc[indices, "gofcards_accession_id"] = [
            _merge_semicolon_values(value, accession_id)
            for value in df.loc[indices, "gofcards_accession_id"]
        ]
    if variant_id:
        df.loc[indices, "gofcards_variant_id"] = [
            _merge_semicolon_values(value, variant_id)
            for value in df.loc[indices, "gofcards_variant_id"]
        ]
    if route:
        df.loc[indices, "gofcards_match_route"] = [
            _merge_semicolon_values(value, route)
            for value in df.loc[indices, "gofcards_match_route"]
        ]


def annotate_exact_nonlof_variants(
    df,
    row_mask=None,
    *,
    context="",
    mechanism_json=str(DEFAULT_MECHANISM_JSON),
):
    """Annotate exclusive exact GoF/DN matches from the canonical JSON.

    Matching is variant-level only: current HGNC symbol plus normalized HGVSp is
    attempted first, followed by current HGNC symbol plus normalized genomic
    allele. Gene-level mechanism history never creates score ``2``.
    """
    if "SYMBOL" not in df.columns:
        return df
    _ensure_exact_nonlof_columns(df)
    if row_mask is None:
        row_mask = pd.Series(True, index=df.index)
    else:
        row_mask = pd.Series(row_mask, index=df.index).fillna(False).astype(bool)
    row_mask = row_mask & df["SYMBOL"].map(_clean_text).ne("")
    if not row_mask.any():
        return df

    hub = GeneMechanismHub(
        mechanism_json=mechanism_json,
        use_hgnc_package=False,
    )
    matched_rows = 0
    hgvsp_matches = 0
    genomic_matches = 0

    if "HGVSp" in df.columns:
        hgvsp_frame = df.loc[row_mask & df["HGVSp"].map(_clean_text).ne(""), ["SYMBOL", "HGVSp"]].copy()
        if not hgvsp_frame.empty:
            hgvsp_frame["_symbol"] = hgvsp_frame["SYMBOL"].map(_clean_text)
            hgvsp_frame["_hgvsp"] = hgvsp_frame["HGVSp"].map(_clean_text)
            for (symbol, hgvsp), idx in hgvsp_frame.groupby(["_symbol", "_hgvsp"], sort=False).groups.items():
                try:
                    match = hub.match_curated_nonlof_variant(symbol, hgvsp=hgvsp)
                except Exception as exc:
                    logger.debug(f"Canonical non-LOF HGVSp lookup failed for {context} {symbol} {hgvsp}: {exc}")
                    continue
                if match.get("mechanism_exclusive"):
                    _record_exact_nonlof_match(df, idx, match)
                    hgvsp_matches += 1
                    matched_rows += len(idx)

    has_genomic = {"chrom", "pos", "ref", "alt"}.issubset(df.columns)
    if has_genomic:
        already_matched = df["variant_mechanism_exclusive"].astype(bool)
        genomic_mask = row_mask & ~already_matched
        for col in ("chrom", "pos", "ref", "alt"):
            genomic_mask = genomic_mask & df[col].map(_clean_text).ne("")
        if genomic_mask.any():
            genomic_cols = ["SYMBOL", "chrom", "pos", "ref", "alt"]
            genomic_frame = df.loc[genomic_mask, genomic_cols].copy()
            genomic_frame["_assembly"] = df.loc[genomic_mask].apply(_row_assembly_for_gof_lookup, axis=1)
            for key, idx in genomic_frame.groupby(genomic_cols + ["_assembly"], sort=False).groups.items():
                symbol, chrom, pos, ref, alt, assembly = key
                try:
                    match = hub.match_curated_nonlof_variant(
                        symbol,
                        chrom=chrom,
                        pos=pos,
                        ref=ref,
                        alt=alt,
                        assembly=assembly,
                        key_type="auto",
                    )
                except Exception as exc:
                    logger.debug(
                        f"Canonical non-LOF genomic lookup failed for {context} {symbol} "
                        f"{chrom}:{pos}:{ref}>{alt} ({assembly}): {exc}"
                    )
                    continue
                if match.get("mechanism_exclusive"):
                    _record_exact_nonlof_match(df, idx, match)
                    genomic_matches += 1
                    matched_rows += len(idx)

    if matched_rows:
        logger.info(
            "Exact canonical non-LOF variant annotations for %s: rows=%s, "
            "unique_hgvsp_matches=%s, unique_genomic_matches=%s",
            context,
            matched_rows,
            hgvsp_matches,
            genomic_matches,
        )
    return df


def _variant_condition_masks(df: pd.DataFrame) -> dict[str, pd.Series]:
    """Parse condition-linked inheritance and penetrance categories.

    The mechanism hub has already selected the histories applicable to this
    variant: compatible known mechanisms plus every included condition whose
    mechanism is unresolved. Criteria must therefore consume these final
    categorical columns rather than recomputing a gene-wide answer from HPO.

    Values are semicolon-delimited sets. Matching is token-exact so, for
    example, ``non_mendelian`` cannot be confused with another substring. Any
    linked ``incomplete`` value is conservative and wins over simultaneous
    ``complete`` histories; consumers implement that by gating on
    ``has_incomplete_penetrance`` first. Upstream normalizes high penetrance to
    complete and emits an empty value when no usable assertion exists.

    ``non_monogenic`` preserves the old distinction for polygenic, digenic and
    oligogenic histories, while ``non_mendelian`` represents the literal HPO
    non-Mendelian category. Both remain visible so each criterion can retain
    its existing gates.
    """
    required = {"variant_inheritance", "variant_penetrance"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise KeyError(
            "variant-level condition annotations are required; missing columns: "
            + ", ".join(missing)
        )

    inheritance = df["variant_inheritance"].fillna("").astype(str).str.lower()
    penetrance = df["variant_penetrance"].fillna("").astype(str).str.lower()

    def has_token(series: pd.Series, token: str) -> pd.Series:
        return series.str.contains(
            rf"(?:^|;){re.escape(token)}(?:;|$)",
            regex=True,
        ).astype(bool)

    recessive = has_token(inheritance, "recessive")
    dominant = has_token(inheritance, "dominant")
    x_linked_recessive = has_token(inheritance, "x_linked_recessive")
    x_linked_dominant = has_token(inheritance, "x_linked_dominant")
    x_linked_unspecified = has_token(inheritance, "x_linked_unspecified")
    y_linked = has_token(inheritance, "y_linked")
    mitochondrial = has_token(inheritance, "mitochondrial")
    polygenic = has_token(inheritance, "polygenic")
    digenic = has_token(inheritance, "digenic")
    oligogenic = has_token(inheritance, "oligogenic")
    non_mendelian = has_token(inheritance, "non_mendelian")

    return {
        "has_recessive": recessive,
        "has_dominant": dominant,
        "has_x_linked_recessive": x_linked_recessive,
        "has_x_linked_dominant": x_linked_dominant,
        "has_x_linked_unspecified": x_linked_unspecified,
        "has_y_linked": y_linked,
        "has_mitochondrial": mitochondrial,
        "has_dominant_like": (
            dominant | x_linked_dominant | y_linked | mitochondrial
        ).astype(bool),
        "has_mendelian": (
            recessive
            | dominant
            | x_linked_recessive
            | x_linked_dominant
            | x_linked_unspecified
            | y_linked
            | mitochondrial
        ).astype(bool),
        "has_non_monogenic": (polygenic | digenic | oligogenic).astype(bool),
        "has_non_mendelian": non_mendelian,
        "has_incomplete_penetrance": has_token(penetrance, "incomplete"),
        "has_complete_penetrance": has_token(penetrance, "complete"),
    }


def _variant_mechanism_masks(
    df: pd.DataFrame,
) -> dict[str, pd.Series]:
    """Return masks from the required variant-level mechanism contract.

    There is intentionally no fallback to ``gene_mech_inher_history``, raw
    HPO inheritance arrays, consequences, or GoFCards tags. The upstream hub
    owns those interpretations and ACMG consumers use its final assertions.

    THE DECISION TREE
    =================

    Every criterion that reasons about mechanism reads its inputs from here,
    so this is the single place the hub's output becomes ACMG-usable. Two
    columns come in; fifteen boolean masks go out, in two families.

        var_plausible_patho_mechs      "dominant_GOF;recessive_LOF"
        variant_effect                 exact_known_* | predicted_LOF_high_
                                       confidence | uncertain
              |
              |  (a KeyError is raised if either is absent -- there is
              |   deliberately no silent default)
              v

    FAMILY 1  what the GENE's condition history says, split by inheritance
    ---------------------------------------------------------------------
    Read out of var_plausible_patho_mechs by matching the tag text.

        has_recessive_compatible     any tag starting "recessive"
        has_dominant_compatible      any tag starting "dominant"

        has_rec_lof_history          tag == "recessive_LOF"
        has_rec_gof_history          tag == "recessive_GOF"
        has_rec_dn_history           tag == "recessive_DN"
        has_rec_unresolved_history   bare "recessive", no mechanism suffix
        has_dom_lof_history          tag == "dominant_LOF"
        has_dom_gof_history          tag == "dominant_GOF"
        has_dom_dn_history           tag == "dominant_DN"
        has_dom_unresolved_history   bare "dominant", no mechanism suffix

    FAMILY 2  what THIS VARIANT does, from variant_effect and the scores
    -------------------------------------------------------------------
        is_exact_gof        variant_gof_score == 2, or effect names GOF
        is_exact_dn         variant_dn_score  == 2, or effect names DN
        is_exact_nonlof     either of the two above
        is_uncertain        effect == "uncertain"
        is_predicted_lof    not uncertain, and variant_lof_score >= 1 (or
                            effect == "predicted_LOF_high_confidence"). True
                            for BOTH actual LOF-prediction grades; read the
                            score itself to tell them apart.
        modern_profile      the row carries any history at all

    There is deliberately no third family conjoining the two. Whether the
    gene's history and this variant's effect line up is each criterion's own
    judgement, made from the two families above, and pre-computing a single
    verdict here would decide it for all of them at once.

    HOW LOSS OF FUNCTION IS GRADED
    ==============================

    Predicted loss of function is graded, and is_predicted_lof deliberately
    does NOT carry that grade -- it is true for both evidence-backed levels. A
    consumer that wants the distinction reads variant_lof_score together with
    variant_effect:

        2   the transcript is destroyed. Either nonsense-mediated decay is
            triggered -- a truncating variant, or a splice variant that
            shifts the reading frame, in neither case escaping decay -- or
            LOFTEE's high-confidence call rescues a variant that does escape.
        1   loss of function is plausible but the transcript may survive:
            escapes decay without a LOFTEE rescue, LOFTEE_OS, a splice
            change that does not shift the frame, a 5'UTR change.
        1   also appears when variant_effect is ``uncertain``. In that case LOF,
            GOF and DN all score 1 as parallel hypotheses, and
            is_predicted_lof remains false because no LOF observation exists.
        0   the mechanism was excluded by a stronger, exclusive variant call.

    Decay-triggered loss of function IS exclusive, and outranks even a
    curated gain-of-function allele: a destroyed transcript makes no protein
    to gain a function. A LOFTEE-rescued score of 2 does not outrank curation.

    The grade decides applicability directly: score 2 lands in
    variant_mechanism_applicable, score 1 in variant_mechanism_uncertain.
    Neither of those two columns is read here. PVS1 in particular reads no
    loss-of-function score at all: it gates on whether loss of function is a
    disease mechanism OF THE GENE, from five independent gene-level sources,
    and grades a null variant's strength from the consequence itself.

    Inheritance and penetrance are intentionally outside this function. They
    are parsed from ``variant_inheritance`` and ``variant_penetrance`` by
    _variant_condition_masks; keeping the two adapters separate preserves the
    distinction between what the variant does and how the selected conditions
    are inherited.
    """
    required = {
        "var_plausible_patho_mechs",
        "variant_effect",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise KeyError(
            "variant-level mechanism annotations are required; missing columns: "
            + ", ".join(missing)
        )

    def text_series(column: str) -> pd.Series:
        return df[column].fillna("").astype(str)

    profile = text_series("var_plausible_patho_mechs")
    accepted = profile
    has_profile = profile.str.strip().ne("")

    rec_accepted = accepted.str.contains(r"(?:^|;)recessive(?:_[^;]+)?(?:;|$)", regex=True)
    dom_accepted = accepted.str.contains(r"(?:^|;)dominant(?:_[^;]+)?(?:;|$)", regex=True)
    x_recessive_accepted = accepted.str.contains(
        r"(?:^|;)x_linked_recessive(?:_[^;]+)?(?:;|$)", regex=True
    )
    x_dominant_accepted = accepted.str.contains(
        r"(?:^|;)x_linked_dominant(?:_[^;]+)?(?:;|$)", regex=True
    )
    x_unspecified_accepted = accepted.str.contains(
        r"(?:^|;)x_linked_unspecified(?:_[^;]+)?(?:;|$)", regex=True
    )

    def profile_mechanism(inheritance: str, mechanism: str) -> pd.Series:
        return profile.str.contains(
            rf"(?:^|;){inheritance}_{mechanism}(?:;|$)",
            regex=True,
        )

    has_rec_lof_history = profile_mechanism("recessive", "LOF")
    has_rec_gof_history = profile_mechanism("recessive", "GOF")
    has_rec_dn_history = profile_mechanism("recessive", "DN")
    has_rec_unresolved_history = profile.str.contains(r"(?:^|;)recessive(?:;|$)", regex=True)
    has_dom_lof_history = profile_mechanism("dominant", "LOF")
    has_dom_gof_history = profile_mechanism("dominant", "GOF")
    has_dom_dn_history = profile_mechanism("dominant", "DN")
    has_dom_unresolved_history = profile.str.contains(r"(?:^|;)dominant(?:;|$)", regex=True)
    has_x_recessive_lof_history = profile_mechanism("x_linked_recessive", "LOF")
    has_x_recessive_gof_history = profile_mechanism("x_linked_recessive", "GOF")
    has_x_recessive_dn_history = profile_mechanism("x_linked_recessive", "DN")
    has_x_recessive_unresolved_history = profile.str.contains(
        r"(?:^|;)x_linked_recessive(?:;|$)", regex=True
    )
    has_x_dominant_lof_history = profile_mechanism(
        "x_linked_dominant", "LOF"
    )
    has_x_dominant_gof_history = profile_mechanism(
        "x_linked_dominant", "GOF"
    )
    has_x_dominant_dn_history = profile_mechanism(
        "x_linked_dominant", "DN"
    )
    has_x_dominant_unresolved_history = profile.str.contains(
        r"(?:^|;)x_linked_dominant(?:;|$)", regex=True
    )

    effect = text_series("variant_effect")
    lof_score = pd.to_numeric(
        df.get("variant_lof_score", pd.Series(0, index=df.index)),
        errors="coerce",
    ).fillna(0)
    gof_score = pd.to_numeric(
        df.get("variant_gof_score", pd.Series(0, index=df.index)),
        errors="coerce",
    ).fillna(0)
    dn_score = pd.to_numeric(
        df.get("variant_dn_score", pd.Series(0, index=df.index)),
        errors="coerce",
    ).fillna(0)
    is_exact_gof = (
        gof_score.eq(2)
        | effect.str.contains(r"(?:^exact_known_|\+)GOF(?:\+|$)", regex=True)
    )
    is_exact_dn = (
        dn_score.eq(2)
        | effect.str.contains(
            r"(?:^exact_known_|\+)DOMINANT_NEGATIVE(?:\+|$)",
            regex=True,
        )
    )
    is_uncertain = effect.eq("uncertain")
    # Score 1 has two meanings that variant_effect keeps distinct: actual but
    # non-exclusive LOF evidence, or one of three parallel hypotheses when no
    # mechanism can be ruled out. Only the former is a predicted LOF variant.
    is_predicted_lof = ~is_uncertain & (
        lof_score.ge(1) | effect.eq("predicted_LOF_high_confidence")
    )

    return {
        "modern_profile": has_profile,
        "is_exact_gof": is_exact_gof.astype(bool),
        "is_exact_dn": is_exact_dn.astype(bool),
        "is_exact_nonlof": (is_exact_gof | is_exact_dn).astype(bool),
        "is_predicted_lof": is_predicted_lof.astype(bool),
        "is_uncertain": is_uncertain.astype(bool),
        "has_recessive_compatible": rec_accepted.astype(bool),
        "has_dominant_compatible": dom_accepted.astype(bool),
        "has_x_linked_recessive_compatible": x_recessive_accepted.astype(bool),
        "has_x_linked_dominant_compatible": x_dominant_accepted.astype(bool),
        "has_x_linked_unspecified_compatible": x_unspecified_accepted.astype(bool),
        "has_rec_lof_history": has_rec_lof_history.astype(bool),
        "has_rec_gof_history": has_rec_gof_history.astype(bool),
        "has_rec_dn_history": has_rec_dn_history.astype(bool),
        "has_rec_unresolved_history": has_rec_unresolved_history.astype(bool),
        "has_dom_lof_history": has_dom_lof_history.astype(bool),
        "has_dom_gof_history": has_dom_gof_history.astype(bool),
        "has_dom_dn_history": has_dom_dn_history.astype(bool),
        "has_dom_unresolved_history": has_dom_unresolved_history.astype(bool),
        "has_x_recessive_lof_history": has_x_recessive_lof_history.astype(bool),
        "has_x_recessive_gof_history": has_x_recessive_gof_history.astype(bool),
        "has_x_recessive_dn_history": has_x_recessive_dn_history.astype(bool),
        "has_x_recessive_unresolved_history": has_x_recessive_unresolved_history.astype(bool),
        "has_x_dominant_lof_history": has_x_dominant_lof_history.astype(bool),
        "has_x_dominant_gof_history": has_x_dominant_gof_history.astype(bool),
        "has_x_dominant_dn_history": has_x_dominant_dn_history.astype(bool),
        "has_x_dominant_unresolved_history": has_x_dominant_unresolved_history.astype(bool),
    }
