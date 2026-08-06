"""Shared binary penetrance vocabulary.

PriVA needs one conservative answer for observation-based ACMG criteria:
whether an apparently unaffected carrier can be informative.  The normalized
condition-level values are therefore deliberately binary:

``complete``
    Explicit complete or high penetrance.

``incomplete``
    Explicit incomplete, moderate, or low penetrance, plus late onset.  This
    value takes precedence when a condition also has a complete assertion.

An empty string means that no usable penetrance assertion was supplied.  Age-
dependent penetrance, variable expressivity, and other onset terms are not
converted to a third ``unknown`` category; they simply do not populate the
penetrance axis.
"""

from __future__ import annotations

import re
from typing import Any


HPO_COMPLETE_PENETRANCE_TERMS = {
    "HP:0034950",  # Typified by complete penetrance
}

HPO_HIGH_PENETRANCE_TERMS = {
    "HP:4000158",  # Typified by high penetrance (80% to <100%)
}

HPO_EXPLICIT_INCOMPLETE_PENETRANCE_TERMS = {
    "HP:0003829",  # Typified by incomplete penetrance
    "HP:4000159",  # Typified by moderate penetrance
    "HP:4000160",  # Typified by low penetrance
}

# Of the onset terms, only explicit late onset blocks an unaffected observation.
HPO_OBSERVATION_CONFOUNDING_TERMS = {
    "HP:0003584",  # Late onset
}

HPO_INCOMPLETE_PENETRANCE_EQUIVALENT_TERMS = (
    HPO_EXPLICIT_INCOMPLETE_PENETRANCE_TERMS
    | HPO_OBSERVATION_CONFOUNDING_TERMS
)

# Only terms that produce one of the two downstream values belong here.  In
# particular HP:0003831 (age-dependent penetrance) and HP:0003828 (variable
# expressivity) are intentionally absent.
HPO_PENETRANCE_STATUS_BY_TERM = {
    **{term: "complete" for term in HPO_COMPLETE_PENETRANCE_TERMS},
    **{term: "complete" for term in HPO_HIGH_PENETRANCE_TERMS},
    **{
        term: "incomplete"
        for term in HPO_INCOMPLETE_PENETRANCE_EQUIVALENT_TERMS
    },
}


def _as_values(values: Any) -> list[str]:
    if isinstance(values, str):
        return [part.strip() for part in re.split(r"[;,]", values) if part.strip()]
    return [str(value).strip() for value in (values or []) if str(value).strip()]


def _normalized_hpo_ids(values: Any) -> set[str]:
    return {value.upper() for value in _as_values(values)}


def recognized_penetrance_hpo_ids(values: Any = ()) -> list[str]:
    """Return recognized binary-penetrance HPO IDs in deterministic order."""
    return sorted(_normalized_hpo_ids(values) & set(HPO_PENETRANCE_STATUS_BY_TERM))


def _resolve_statuses(statuses: Any = ()) -> str:
    normalized = {str(value).strip().lower() for value in (statuses or []) if str(value).strip()}
    if "incomplete" in normalized:
        return "incomplete"
    if "complete" in normalized or "high" in normalized:
        return "complete"
    return ""


def normalize_penetrance_hpo_ids(values: Any = ()) -> str:
    """Normalize condition-linked HPO IDs to ``complete``, ``incomplete``, or empty."""
    statuses = [
        HPO_PENETRANCE_STATUS_BY_TERM[hpo_id]
        for hpo_id in recognized_penetrance_hpo_ids(values)
    ]
    return _resolve_statuses(statuses)


def normalize_penetrance_text(values: Any = ()) -> str:
    """Normalize controlled source text such as G2P cross-cutting modifiers."""
    text = " ; ".join(_as_values(values)).lower().replace("_", "-")
    text = re.sub(r"\s+", " ", text)
    incomplete_phrases = (
        "incomplete penetrance",
        "moderate penetrance",
        "low penetrance",
        "late onset",
    )
    if any(phrase in text for phrase in incomplete_phrases):
        return "incomplete"
    if "complete penetrance" in text or "high penetrance" in text:
        return "complete"
    return ""


def normalize_penetrance_evidence(
    *,
    hpo_ids: Any = (),
    raw_values: Any = (),
) -> str:
    """Combine source assertions with incomplete-over-complete precedence."""
    return _resolve_statuses(
        [
            normalize_penetrance_hpo_ids(hpo_ids),
            normalize_penetrance_text(raw_values),
        ]
    )


def normalize_penetrance_values(values: Any = ()) -> str:
    """Normalize cache statuses, HPO IDs, or controlled raw values.

    This compatibility adapter lets runtime readers consume both the modern
    condition cache's normalized status list and older HPO-ID-only assertions.
    Literal ``unknown`` is treated as no assertion and is never emitted.
    """
    raw_values = _as_values(values)
    statuses = [value.lower() for value in raw_values]
    direct = _resolve_statuses(statuses)
    return _resolve_statuses(
        [
            direct,
            normalize_penetrance_hpo_ids(raw_values),
            normalize_penetrance_text(raw_values),
        ]
    )
