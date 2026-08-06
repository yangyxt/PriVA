"""Shared constants and side-effect-free helpers for gene mechanism modules.

This module owns the contracts used by more than one stage. It performs no
resource loading and is never executed as a pipeline entry point.
"""

from __future__ import annotations

import logging
import math
import re
from pathlib import Path
from typing import Any

import numpy as np

from hpo_penetrance import HPO_PENETRANCE_STATUS_BY_TERM


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
    DATA_DIR / "patho_mechanism" / "hpo_condition_mechanism_cache.json.gz"
)
HPO_CONDITION_MECHANISM_SCHEMA_VERSION = "1.2"
# The non-LOF cache is the mechanism input. It used to be chosen over an older
# combined cache only when it existed on disk, but it is the file the repository
# ships, so that choice had one possible outcome and the other name is gone.
DEFAULT_MECHANISM_JSON = (
    DATA_DIR / "patho_mechanism" / "gene_nonlof_mechanism_curated_assertions.json.gz"
)
DEFAULT_DDG2P_MECHANISM_EVIDENCE = (
    DATA_DIR / "patho_mechanism" / "gene_pathogenic_mechanism_evidence.tsv"
)
DEFAULT_GOFCARDS_EXACT_GOF_HGVSP = DATA_DIR / "gofcards" / "gofcards_exact_gof.json.gz"

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
UNRESOLVED_MECHANISM = "UNRESOLVED"
# Mechanisms that can be established exclusively for one allele, meaning the
# other two are then held to be absent. Loss of function belongs here only by
# the decay route: if nonsense-mediated decay destroys the transcript there is
# no protein left to gain a function or to poison a complex.
EXACT_SEQUENCE_MECHANISMS = {"GOF", "DOMINANT_NEGATIVE", "LOF"}
# The subset a curated allele match can establish. Loss of function is absent
# because PriVA has no curated loss-of-function allele database; its exclusive
# claim comes from the decay prediction, not from curation.
CURATED_EXACT_MECHANISMS = {"GOF", "DOMINANT_NEGATIVE"}
VARIANT_MECHANISM_SCORE_KEYS = ("LOF", "GOF", "DOMINANT_NEGATIVE")
# Every source that may state a condition's mechanism. GoFCards is included:
# a curated gain-of-function allele in a gene IS that gene's curated history
# for the condition it was curated against, and 97 mechanism blocks have no
# other source at all. The old worry -- that one allele would give every
# variant in the gene a gain-of-function history -- does not arise, because
# select_condition_histories_for_variant keeps only the histories the query
# variant's own mechanism reaches. A predicted loss-of-function allele never
# sees the gain-of-function history. When no mechanism-specific variant evidence
# exists, LOF, GOF and DN remain parallel score-1 hypotheses, so all three kinds
# of asserted history remain visible beside mechanism-unresolved conditions.
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
    "penetrance_raw",
    "penetrance_hpo_ids",
    "normalized_penetrance",
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
HPO_PENETRANCE_TERMS = set(HPO_PENETRANCE_STATUS_BY_TERM)
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
    "HP:0034857",  # Highly variable age of onset
    "HP:0003587",  # Insidious / gradual onset
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
}
DDG2P_DOMINANT_LOF_INHERITANCE = {
    "monoallelic_autosomal",
}
DDG2P_X_LINKED_DOMINANT_LOF_INHERITANCE = {
    "monoallelic_X_heterozygous"
}
DDG2P_X_LINKED_UNSPECIFIED_LOF_INHERITANCE = {
    "monoallelic_X",
    "monoallelic_X_hemizygous",
}
DDG2P_Y_LINKED_LOF_INHERITANCE = {"monoallelic_Y_hemizygous"}
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
    # A property of the GENE, not of this variant: does the condition cache
    # record loss of function as a disease mechanism for it at all. PVS1 reads
    # this and decides strength for the particular null variant itself.
    "gene_lof_mechanism_history",
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
