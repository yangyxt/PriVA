"""Shared HPO penetrance and carrier-observability vocabulary.

PriVA's benign and segregation criteria ask a narrower question than the HPO
ontology does: can an apparently unaffected carrier be treated as evidence
against pathogenicity?  A condition can make that observation uninformative
through literal incomplete penetrance, through age-dependent or delayed onset,
or through variable expressivity.  Those source concepts must remain distinct
in the condition cache, but they deliberately collapse to the same downstream
``incomplete`` category.

The boundary is conservative but not unlimited:

* high penetrance remains ``high``; it is not treated as incomplete;
* sex-limited expression is not treated as a penetrance equivalent;
* slowly progressive disease describes change after onset, not whether disease
  is present, so it is not a penetrance equivalent either.

Keeping the vocabulary here prevents the cache builder, mechanism hub, and
legacy inheritance bridge from silently using different HPO identifiers.
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

# These terms directly state that penetrance is incomplete or age-dependent.
HPO_EXPLICIT_INCOMPLETE_PENETRANCE_TERMS = {
    "HP:0003829",  # Typified by incomplete penetrance
    "HP:0003831",  # Typified by age-related disease onset / age-dependent penetrance
    "HP:4000159",  # Typified by moderate penetrance
    "HP:4000160",  # Typified by low penetrance
}

# These remain onset or expressivity assertions in HPO.  PriVA treats them as
# incomplete only for the operational question of whether an unaffected
# observation can contradict the condition.
HPO_OBSERVATION_CONFOUNDING_TERMS = {
    "HP:0034857",  # Typified by highly variable age of onset
    "HP:0003581",  # Adult onset
    "HP:0011462",  # Young adult onset
    "HP:0003596",  # Middle age onset
    "HP:0003584",  # Late onset
    "HP:0003587",  # Insidious onset (exact synonym: gradual onset)
    "HP:0003828",  # Variable expressivity
}

HPO_INCOMPLETE_PENETRANCE_EQUIVALENT_TERMS = (
    HPO_EXPLICIT_INCOMPLETE_PENETRANCE_TERMS
    | HPO_OBSERVATION_CONFOUNDING_TERMS
)

# Values stored in the condition cache.  The original HPO identifier and its
# provenance remain beside the value, so "incomplete" never erases why the
# condition was placed in that downstream category.
HPO_PENETRANCE_STATUS_BY_TERM = {
    "HP:0003829": "incomplete",
    "HP:0003831": "incomplete",
    "HP:0034950": "complete",
    "HP:4000158": "high",
    "HP:4000159": "moderate",
    "HP:4000160": "low",
    **{term: "incomplete" for term in HPO_OBSERVATION_CONFOUNDING_TERMS},
}


def _normalized_hpo_ids(values: Any) -> set[str]:
    if isinstance(values, str):
        raw_values = re.split(r"[;,]", values)
    else:
        raw_values = list(values or [])
    return {
        str(value).strip().upper()
        for value in raw_values
        if str(value).strip()
    }


def normalize_penetrance_hpo_ids(values: Any = ()) -> str:
    """Return ``complete``, ``high``, ``incomplete``, or ``unknown``.

    ``incomplete`` takes precedence because any linked assertion that makes an
    unaffected carrier non-informative must protect downstream logic.  ``high``
    remains separate from both complete and incomplete; when high and complete
    assertions conflict, high is the more cautious category because it does not
    activate rules reserved for explicitly complete penetrance.
    """
    hpo_ids = _normalized_hpo_ids(values)
    if hpo_ids & HPO_INCOMPLETE_PENETRANCE_EQUIVALENT_TERMS:
        return "incomplete"
    if hpo_ids & HPO_HIGH_PENETRANCE_TERMS:
        return "high"
    if hpo_ids & HPO_COMPLETE_PENETRANCE_TERMS:
        return "complete"
    return "unknown"
