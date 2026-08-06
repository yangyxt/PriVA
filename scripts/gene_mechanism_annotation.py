"""Variant-to-condition selection and DataFrame annotation.

This is the joining stage: it combines a query variant effect with compatible
condition histories and attaches the stable output contract to table rows.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from gene_mechanism_common import (
    DEFAULT_CLINGEN_DOSAGE,
    DEFAULT_HGNC_TABLE,
    DEFAULT_HPO_COLLAPSED,
    DEFAULT_HPO_CONDITION_MECHANISM_CACHE,
    DEFAULT_LOEUF_TABLE,
    UNRESOLVED_MECHANISM,
    VARIANT_MECHANISM_OUTPUT_COLUMNS,
    VARIANT_MECHANISM_SCORE_KEYS,
    _clean,
    _deduplicate_by,
    _exact_mechanisms_from_effect,
)
from gene_mechanism_conditions import (
    NON_MENDELIAN_INHERITANCE,
    condition_cache_mechanism_assertions,
    gene_has_lof_mechanism_history,
    gene_inheritance_from_constraint,
    normalize_inheritance,
    normalize_inheritance_histories,
    normalize_penetrance,
)
from gene_mechanism_resources import GeneMechanismHub
from gene_mechanism_variants import infer_query_variant_effect


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
    if inheritance in NON_MENDELIAN_INHERITANCE:
        # These histories exist to preserve a conservative criterion gate, not
        # to create molecular-mechanism profile tags such as polygenic_LOF.
        return ""
    mechanism = _clean(history.get("mechanism")).upper() or "UNRESOLVED"
    if mechanism == "UNRESOLVED":
        return inheritance or "uncertain"
    mechanism = "DN" if mechanism == "DOMINANT_NEGATIVE" else mechanism
    return f"{inheritance}_{mechanism}" if inheritance else mechanism


def _compact_profile_tags(assertions: list[dict[str, Any]]) -> list[str]:
    tags = list(
        dict.fromkeys(
            tag
            for assertion in assertions
            if (tag := _mechanism_profile_tag(assertion))
        )
    )
    inheritance_values = (
        "x_linked_recessive",
        "x_linked_dominant",
        "x_linked_unspecified",
        "mitochondrial",
        "recessive",
        "dominant",
        "y_linked",
    )

    def split_profile_tag(tag: str) -> tuple[str, str]:
        for inheritance in inheritance_values:
            prefix = inheritance + "_"
            if tag.startswith(prefix):
                return inheritance, tag[len(prefix):]
        return "", ""

    split_tags = [split_profile_tag(tag) for tag in tags]
    inheritance_with_mechanism = {
        inheritance for inheritance, mechanism in split_tags if mechanism
    }
    qualified_mechanisms = {
        mechanism for _inheritance, mechanism in split_tags if mechanism
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

    * an exact known GOF or DN allele selects that known mechanism;
    * a high-confidence predicted LOF allele selects curated or deduced LOF
      histories;
    * an allele with no mechanism-specific evidence keeps LOF, GOF and DN as
      parallel plausible hypotheses; and
    * every case also keeps ``UNRESOLVED`` histories, because an absent
      condition-mechanism assertion cannot establish incompatibility.

    Exact curated mechanisms are exclusive among KNOWN mechanisms: a
    consequence-based loss-of-function prediction cannot reintroduce a known
    incompatible history after an exact match. It does not exclude a condition
    whose mechanism was never asserted.

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
    elif effect == "uncertain":
        allowed = set(VARIANT_MECHANISM_SCORE_KEYS)
    else:
        allowed = set()

    histories: list[dict[str, Any]] = []
    for assertion in assertions:
        mechanism = _clean(assertion.get("mechanism")).upper()
        if mechanism != UNRESOLVED_MECHANISM and mechanism not in allowed:
            continue
        for inheritance, x_linked in normalize_inheritance_histories(
            assertion.get("allelic_requirement"),
            assertion.get("hpo_inheritance_modes"),
        ):
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
                        assertion.get("penetrance_statuses")
                        or assertion.get("penetrance_hpo_ids")
                    ),
                }
            )
    return histories


HISTORY_IDENTITY_FIELDS = ("mechanism", "condition", "inheritance", "x_linked", "penetrance")


def _classify_variant_applicability(
    histories: list[dict[str, Any]],
    scores: dict[str, int],
) -> dict[str, Any]:
    """Split the selected histories into established and merely possible.

    Step 3 has already removed the histories this variant's mechanism cannot
    reach, so nothing here is incompatible. What remains is the one distinction
    the scores already draw, read straight off them rather than re-derived:

        score 2  established   -> applicable
        score 1  plausible     -> uncertain
        UNRESOLVED mechanism   -> uncertain

    This deliberately does not consider what any downstream criterion will do
    with the answer. Whether a particular criterion fires, and how strongly, is
    that criterion's own business; this function only reports how well the
    variant's mechanism is established for each history.
    """
    applicable: list[str] = []
    uncertain: list[str] = []
    for history in histories:
        mechanism = _clean(history.get("mechanism")).upper()
        established = (
            mechanism != UNRESOLVED_MECHANISM
            and int(scores.get(mechanism, 0) or 0) == 2
        )
        tag = _mechanism_profile_tag(history)
        if not tag:
            continue
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

    Condition histories come only from the variant's own annotation and the
    HPO condition cache. Gene-wide signals -- a ClinVar pathogenic history, a
    constrained LOEUF, a high average AlphaMissense score, a ClinGen dosage
    call -- cannot create a condition or molecular mechanism. Constraint is
    used only as a clearly labelled inheritance fallback when no surviving
    condition history supplies inheritance.
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
    condition_records = hub._load_condition_cache()
    assertion_cache: dict[str, list[dict[str, Any]]] = {}
    # A gene-level answer, so it is computed once per gene rather than once per
    # row, the same way assertion_cache already works.
    lof_history_cache: dict[str, bool] = {}
    for _, row in out.iterrows():
        gene = row[symbol_col]
        symbol = hub.resolve_symbol(gene)
        gene_record = condition_records.get(symbol, {})
        if symbol not in assertion_cache:
            assertion_cache[symbol] = condition_cache_mechanism_assertions(gene_record)
        assertions = list(assertion_cache[symbol])
        if symbol not in lof_history_cache:
            lof_history_cache[symbol] = gene_has_lof_mechanism_history(gene_record)

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
            {
                "LOF": effect_call["variant_lof_score"],
                "GOF": effect_call["variant_gof_score"],
                "DOMINANT_NEGATIVE": effect_call["variant_dn_score"],
            },
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
        row_chromosome = _clean(row.get("chrom", "")).lower().removeprefix("chr")
        matched_x_linked = (
            row_chromosome == "x" or any(h["x_linked"] for h in histories)
        )
        if matched_inheritance:
            selected_mechanisms = {
                _clean(history.get("mechanism")).upper() for history in histories
            }
            selected_known = bool(selected_mechanisms - {UNRESOLVED_MECHANISM})
            selected_unresolved = UNRESOLVED_MECHANISM in selected_mechanisms
            if selected_known and selected_unresolved:
                basis = "matched_and_unresolved_history"
            elif selected_unresolved:
                basis = "unresolved_condition_history"
            else:
                basis = "matched_history"
        else:
            # Scope and mechanism gates apply to individual conditions only.
            # Once those conditions have been discarded, they cannot veto an
            # independent gene-level inheritance fallback.
            inferred_inheritance = gene_inheritance_from_constraint(
                symbol,
                clingen=hub._load_clingen(),
                loeuf=hub._load_loeuf(),
                chromosome=row.get("chrom", ""),
            )
            matched_inheritance = (
                [inferred_inheritance] if inferred_inheritance else []
            )
            matched_x_linked = inferred_inheritance.startswith("x_")
            basis = "gene_constraint"

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
        variant_outputs["gene_lof_mechanism_history"].append(
            "true" if lof_history_cache[symbol] else "false"
        )
        variant_outputs["variant_penetrance"].append(
            ";".join(
                dict.fromkeys(
                    h["penetrance"] for h in histories if h["penetrance"]
                )
            )
        )

    out[output_col] = plausible_mechanism_values
    for column, values in variant_outputs.items():
        out[column] = values
    return out
