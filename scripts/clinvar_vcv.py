#!/usr/bin/env python3
"""Stream ClinVar VCV XML and retain high-review exact GoFCards matches."""

from __future__ import annotations

import contextlib
import gzip
import io
import json
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, BinaryIO, Callable, Iterable, Iterator

import pandas as pd


REVIEW_STARS = {
    "no classification provided": 0,
    "no assertion criteria provided": 0,
    "no classifications from unflagged records": 0,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, multiple submitters, no conflicts": 2,
    "reviewed by expert panel": 3,
    "practice guideline": 4,
}

ASSEMBLY_TO_BUILD = {
    "GRCh37": "hg19",
    "GRCh38": "hg38",
}

# The allele-record shape and the flattener live here, in the lower-level
# module, so both consumers share one definition without an import cycle.



def iter_gofcards_variants(cache: dict[str, Any]) -> Iterator[tuple[str, str, dict[str, Any]]]:
    """Walk the GoFCards cache, yielding (gene symbol, variant id, variant block).

    This is the single reader for the nested cache.  It does not reshape
    anything: consumers receive the variant block exactly as stored, with its
    ``record`` and its per-assembly ``genomic`` and ``transcripts`` intact, and
    read the fields they need by their real paths.

    Structure walked here (see docs/GOFCARDS_CACHE_SCHEMA.md):

        genes -> <SYMBOL> -> variants -> <loc_...._grch37>
                                          |- record.eligibility / source / evidence
                                          '- assemblies -> hg19|hg38
                                                            |- genomic {chrom,pos,ref,alt}
                                                            '- transcripts -> <tx>
                                                                              |- by_hgvsc
                                                                              '- by_hgvsp
    """
    for symbol, gene in (cache.get("genes") or {}).items():
        for variant_id, variant in (gene.get("variants") or {}).items():
            yield symbol, variant_id, variant


def clean_text(value: Any) -> str:
    """Normalize a field to text, mapping every placeholder for "empty" to "".

    VEP writes `-` for an absent field and GoFCards writes `_`, `.` or `-`.
    Treating those as text would produce match keys like `hgvsp_key="-"`, which
    would collide across every variant that has no protein change.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "_", ".", "-"} else text


def normalize_chrom(value: Any) -> str:
    return re.sub(r"^chr", "", clean_text(value), flags=re.IGNORECASE)


def normalize_allele(value: Any) -> str:
    text = clean_text(value)
    return "" if text in {"-", "."} else text.upper()


def normalize_int(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    text = str(value).strip()
    return text[:-2] if text.endswith(".0") else text


def variant_id_of(chrom: Any, start: Any, ref: Any, alt: Any) -> str:
    """Mint the stable variant identifier, `loc_<chrom>:<start>:<ref>-><alt>_grch37`.

    This lives beside the cache readers because it is part of the cache
    contract: it is the second-level key the cache is organized under, and both
    the normalizer that writes the cache and every consumer that looks a variant
    up must mint it the same way. Two copies of this rule would be two ways to
    name the same variant.

    GoFCards publishes no identifier that is both complete and unique -- its
    ``Accession`` field is a ClinVar cross-reference covering only 61% of the
    catalogue -- so identity has to be constructed.  Chromosome, start,
    reference and alternate are the minimum that separates every record:
    chromosome and position alone collide for 104 variants, because hotspot
    codons carry several different substitutions at one base (TP53 chr17:7577120
    is C>T, C>A and C>G).

    The ``grch37`` suffix is deliberate and honest.  GoFCards states coordinates
    on GRCh37 only, so the identifier is minted from that build.  It is an opaque
    label: never parse it, and never read a GRCh38 position from it.

    Missing alleles are written ``-`` exactly as GoFCards writes them, so an
    insertion reads ``loc_6:91261902:-->TACTAC_grch37``.  Compound variants need
    no special handling: the alleles are written verbatim, so an in-cis pair of
    substitutions reads ``loc_7:140453136:AC->CT_grch37``.
    """
    def allele(value: Any) -> str:
        text = normalize_allele(value)
        return text if text else "-"

    return (f"loc_{normalize_chrom(chrom)}:{normalize_int(start)}:"
            f"{allele(ref)}->{allele(alt)}_grch37")


@contextlib.contextmanager
def open_text(path: Path | str, mode: str = "rt") -> Iterator[Any]:
    """Open a text file, transparently gzipped when the name ends in ``.gz``.

    Every cache in this pipeline is read and written through here, so the
    decision to compress one is a change of filename and nothing else.

    Writes keep neither a modification time nor a stored filename in the gzip
    header, so identical content always compresses to identical bytes. These
    caches are tracked in git: with either field present, every rebuild would
    look like a change and add another copy to the repository's history, and a
    file written to a temporary name would differ from the same content written
    directly.
    """
    path = Path(path)
    if str(path).endswith(".gz"):
        if "r" in mode:
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                yield handle
        else:
            with path.open("wb") as raw, \
                    gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed, \
                    io.TextIOWrapper(compressed, encoding="utf-8") as handle:
                yield handle
    else:
        with path.open(mode, encoding="utf-8") as handle:
            yield handle


def load_gofcards_cache(path: Path) -> dict[str, Any]:
    """Read the nested GoFCards cache and check it is the right shape."""
    with open_text(path) as handle:
        cache = json.load(handle)
    if "genes" not in cache:
        raise ValueError(f"GoFCards cache has no 'genes' block: {path}")
    meta = cache.get("metadata") or {}
    if meta.get("source") != "GoFCards" or meta.get("mechanism") != "GOF":
        raise ValueError(f"GoFCards cache does not declare GoFCards/GOF: {path}")
    return cache


def gofcards_variant_is_eligible(variant: dict[str, Any]) -> bool:
    """Return whether a variant may be used as gain-of-function evidence.

    Only the explicit ``eligible`` state qualifies.  Every quarantine state --
    gene discordance and reviewed non-GOF mechanisms alike -- is retained in the
    cache for audit and must never match.
    """
    return ((variant.get("record") or {}).get("eligibility") or {}).get("status") == "eligible"


def gofcards_variant_genomic_keys(variant: dict[str, Any]) -> dict[str, str]:
    """Return {assembly: 'chrom|pos|ref|alt'} for every build the variant has.

    A variant that failed liftover simply has no hg38 entry, so the caller sees
    one key instead of two rather than an empty string it has to test for.
    """
    keys: dict[str, str] = {}
    for assembly, block in (variant.get("assemblies") or {}).items():
        coords = block.get("genomic") or {}
        if coords.get("chrom") and coords.get("pos") is not None:
            keys[assembly] = (f"{coords['chrom']}|{coords['pos']}"
                              f"|{coords.get('ref','')}|{coords.get('alt','')}")
    return keys


def _text(element: ET.Element | None) -> str:
    return "" if element is None or element.text is None else element.text.strip()


def _child(element: ET.Element | None, name: str) -> ET.Element | None:
    if element is None:
        return None
    for item in element:
        if item.tag.rsplit("}", 1)[-1] == name:
            return item
    return None


def _children(element: ET.Element | None, name: str) -> list[ET.Element]:
    if element is None:
        return []
    return [item for item in element if item.tag.rsplit("}", 1)[-1] == name]


def _path(element: ET.Element | None, *names: str) -> ET.Element | None:
    current = element
    for name in names:
        current = _child(current, name)
        if current is None:
            break
    return current


def _descendants(element: ET.Element | None, name: str) -> Iterable[ET.Element]:
    if element is None:
        return []
    return (
        item
        for item in element.iter()
        if item.tag.rsplit("}", 1)[-1] == name
    )


def _nonempty(data: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in data.items()
        if value is not None and value != "" and value != [] and value != {}
    }


def _int_or_text(value: str | None) -> int | str | None:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except ValueError:
        return value


def review_stars(review_status: str) -> int:
    """Map ClinVar's exact aggregate germline review status to stars."""
    normalized = re.sub(r"\s+", " ", review_status.strip().lower())
    return REVIEW_STARS.get(normalized, 0)


def normalize_vcf_key(chromosome: str, position: str, ref: str, alt: str) -> str:
    chrom = chromosome.strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    if chrom == "MT":
        chrom = "M"
    return "|".join([chrom, position.strip(), ref.strip().upper(), alt.strip().upper()])


def _gofcards_assertion_key(row: dict[str, Any]) -> tuple[str, str]:
    """Return the curated source gene plus stable compact allele identity."""
    symbol = str(
        row.get("GoFCards_HGNC_Symbol") or row.get("HGNC_Symbol") or ""
    ).strip().upper()
    allele = str(
        row.get("gofcards_variant_id") or row.get("allele_key") or ""
    ).strip()
    if not symbol or not allele:
        raise ValueError("GoFCards row lacks a curated gene or stable allele ID")
    return symbol, allele


def partition_gofcards_variants(
    cache: dict[str, Any],
) -> tuple[list[tuple[str, str, dict[str, Any]]], list[tuple[str, str, dict[str, Any]]]]:
    """Split the cache into variants that may match and variants that may not.

    Eligibility is a property of the variant and is already decided upstream, so
    this is a straight partition on that verdict rather than the row-level
    reconciliation the flat table used to need.
    """
    eligible: list[tuple[str, str, dict[str, Any]]] = []
    quarantined: list[tuple[str, str, dict[str, Any]]] = []
    for symbol, variant_id, variant in iter_gofcards_variants(cache):
        target = eligible if gofcards_variant_is_eligible(variant) else quarantined
        target.append((symbol, variant_id, variant))
    return eligible, quarantined


def gofcards_variation_id(variant: dict[str, Any]) -> str:
    """Return GoFCards' own ClinVar VariationID for a variant, if it has one.

    Present for roughly 61% of the catalogue -- only the variants that are in
    ClinVar at all. Where it exists it is an exact identifier and beats matching
    on coordinates, which is why the injection step tries it first.
    """
    annotations = ((variant.get("record") or {}).get("annotations") or {})
    return str(annotations.get("clinvar_variation_id") or "").strip()


GofcardsVariant = tuple[str, str, dict[str, Any]]


def index_gofcards_variants(
    cache: dict[str, Any],
    key: Callable[[str, str, dict[str, Any]], Iterable[Any]],
    *,
    eligible_only: bool = True,
) -> dict[Any, list[GofcardsVariant]]:
    """Index the cache under whatever key the caller declares.

    The cache has exactly one identity -- gene, then variant identifier -- and
    every other key is a secondary route to that same identity. So there is one
    traversal and one index builder, and each consumer supplies only the key it
    needs. That is why this replaces the several near-identical index builders
    that each walked the cache and rebuilt a flat row of their own.

    The stored value is always ``(symbol, variant_id, variant)`` with the
    variant left nested and uncopied. A projection here would be a second
    description of the same data, free to drift from the first -- which is
    exactly how the matching path went silently dead once already.

    ``key`` returns an iterable because one variant legitimately lands under
    several keys: a genomic index holds it under both GRCh37 and GRCh38.
    """
    index: dict[Any, list[GofcardsVariant]] = defaultdict(list)
    for symbol, variant_id, variant in iter_gofcards_variants(cache):
        if eligible_only and not gofcards_variant_is_eligible(variant):
            continue
        for entry_key in key(symbol, variant_id, variant):
            index[entry_key].append((symbol, variant_id, variant))
    return dict(index)


def gofcards_genomic_index_key(
    _symbol: str, _variant_id: str, variant: dict[str, Any]
) -> Iterable[tuple[str, str]]:
    """Key a variant by (build, normalized VCF allele) on every build it has."""
    for build, genomic_key in gofcards_variant_genomic_keys(variant).items():
        yield build, normalize_vcf_key(*genomic_key.split("|"))


def gofcards_variation_id_index_key(
    _symbol: str, _variant_id: str, variant: dict[str, Any]
) -> Iterable[str]:
    """Key a variant by GoFCards' own ClinVar VariationID, when it has one."""
    variation_id = gofcards_variation_id(variant)
    if variation_id:
        yield variation_id


def _open_xml(path: Path) -> BinaryIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rb")
    return path.open("rb")


def _allele_locations(allele: ET.Element) -> dict[str, dict[str, Any]]:
    locations: dict[str, dict[str, Any]] = {}
    for sequence in _descendants(_child(allele, "Location"), "SequenceLocation"):
        assembly = sequence.get("Assembly", "").strip()
        build = ASSEMBLY_TO_BUILD.get(assembly)
        chrom = sequence.get("Chr", "").strip()
        pos = sequence.get("positionVCF", "").strip()
        ref = sequence.get("referenceAlleleVCF", "").strip()
        alt = sequence.get("alternateAlleleVCF", "").strip()
        if not build or not all((chrom, pos, ref, alt)):
            continue
        locations[assembly] = _nonempty(
            {
                "assembly_accession_version": sequence.get(
                    "AssemblyAccessionVersion", ""
                ),
                "chromosome": chrom,
                "position_vcf": _int_or_text(pos),
                "reference_allele_vcf": ref,
                "alternate_allele_vcf": alt,
                "vcf_key": normalize_vcf_key(chrom, pos, ref, alt),
            }
        )
    return locations


def _allele_genes(allele: ET.Element) -> list[dict[str, str]]:
    genes = []
    for gene in _children(_child(allele, "GeneList"), "Gene"):
        genes.append(
            _nonempty(
                {
                    "symbol": gene.get("Symbol", ""),
                    "hgnc_id": gene.get("HGNC_ID", ""),
                    "entrez_id": gene.get("GeneID", ""),
                    "relationship_type": gene.get("RelationshipType", ""),
                }
            )
        )
    return genes


def _allele_hgvs(allele: ET.Element) -> list[dict[str, str]]:
    rows = []
    for hgvs in _children(_child(allele, "HGVSlist"), "HGVS"):
        nucleotide = _child(hgvs, "NucleotideExpression")
        protein = _child(hgvs, "ProteinExpression")
        expression = _text(_child(nucleotide, "Expression")) or _text(
            _child(protein, "Expression")
        )
        if expression:
            rows.append(
                _nonempty(
                    {
                        "type": hgvs.get("Type", ""),
                        "expression": expression,
                        "sequence_accession_version": (
                            nucleotide.get("sequenceAccessionVersion", "")
                            if nucleotide is not None
                            else protein.get("sequenceAccessionVersion", "")
                            if protein is not None
                            else ""
                        ),
                    }
                )
            )
    return rows


def _variant_alleles(classified: ET.Element) -> list[tuple[ET.Element, str, str]]:
    """Return aggregate variant alleles with explicit classification scope."""
    direct = _child(classified, "SimpleAllele")
    if direct is not None:
        return [(direct, "simple_allele", "allele")]

    haplotype = _child(classified, "Haplotype")
    if haplotype is not None:
        return [
            (allele, "haplotype", "haplotype_component")
            for allele in _children(haplotype, "SimpleAllele")
        ]

    genotype = _child(classified, "Genotype")
    if genotype is None:
        return []
    rows: list[tuple[ET.Element, str, str]] = []
    for child in genotype:
        tag = child.tag.rsplit("}", 1)[-1]
        if tag == "SimpleAllele":
            rows.append((child, "genotype", "genotype_component"))
        elif tag == "Haplotype":
            rows.extend(
                (allele, "genotype", "genotype_haplotype_component")
                for allele in _children(child, "SimpleAllele")
            )
    return rows


def _condition_values(rcv: ET.Element) -> list[dict[str, str]]:
    values = []
    condition_list = _child(rcv, "ClassifiedConditionList")
    for condition in _children(condition_list, "ClassifiedCondition"):
        values.append(
            _nonempty(
                {
                    "name": _text(condition),
                    "database": condition.get("DB", ""),
                    "id": condition.get("ID", ""),
                }
            )
        )
    return values


def _condition_tokens(conditions: list[dict[str, str]]) -> set[str]:
    tokens: set[str] = set()
    for condition in conditions:
        db = condition.get("database", "").strip().lower()
        ident = condition.get("id", "").strip().lower()
        name = re.sub(r"\s+", " ", condition.get("name", "").strip().lower())
        if ident:
            tokens.add(ident)
            if db:
                tokens.add(f"{db}:{ident}")
        if name:
            tokens.add(f"name:{name}")
    return tokens


def _trait_tokens(assertion: ET.Element) -> set[str]:
    tokens: set[str] = set()
    for trait in _descendants(_child(assertion, "TraitSet"), "Trait"):
        for xref in _children(trait, "XRef"):
            db = xref.get("DB", "").strip().lower()
            ident = xref.get("ID", "").strip().lower()
            if ident:
                tokens.add(ident)
                if db:
                    tokens.add(f"{db}:{ident}")
        for preferred in _descendants(_child(trait, "Name"), "ElementValue"):
            value = re.sub(r"\s+", " ", _text(preferred).lower())
            if value:
                tokens.add(f"name:{value}")
    return tokens


def _trait_mapping_by_assertion(classified: ET.Element) -> dict[str, list[dict[str, str]]]:
    mapping: dict[str, list[dict[str, str]]] = defaultdict(list)
    for item in _children(_child(classified, "TraitMappingList"), "TraitMapping"):
        assertion_id = item.get("ClinicalAssertionID", "")
        medgen = _child(item, "MedGen")
        row = _nonempty(
            {
                "clinical_assertion_id": assertion_id,
                "trait_type": item.get("TraitType", ""),
                "mapping_type": item.get("MappingType", ""),
                "mapping_value": item.get("MappingValue", ""),
                "mapping_ref": item.get("MappingRef", ""),
                "medgen_name": medgen.get("Name", "") if medgen is not None else "",
                "medgen_cui": medgen.get("CUI", "") if medgen is not None else "",
            }
        )
        if assertion_id:
            mapping[assertion_id].append(row)
    return dict(mapping)


def _mapping_tokens(rows: list[dict[str, str]]) -> set[str]:
    tokens: set[str] = set()
    for row in rows:
        cui = row.get("medgen_cui", "").strip().lower()
        name = re.sub(r"\s+", " ", row.get("medgen_name", "").strip().lower())
        if cui:
            tokens.update((cui, f"medgen:{cui}"))
        if name:
            tokens.add(f"name:{name}")
    return tokens


def _attribute_values(assertion: ET.Element, wanted_type: str) -> list[str]:
    values = []
    for attribute_set in _children(assertion, "AttributeSet"):
        attribute = _child(attribute_set, "Attribute")
        if attribute is not None and attribute.get("Type") == wanted_type:
            value = _text(attribute)
            if value:
                values.append(value)
    return sorted(set(values))


def _cooccurrence(item: ET.Element) -> dict[str, Any]:
    alleles = []
    for allele in _children(item, "AlleleDescSet"):
        significance = _child(allele, "ClinicalSignificance")
        alleles.append(
            _nonempty(
                {
                    "name": _text(_child(allele, "Name")),
                    "relative_orientation": _text(
                        _child(allele, "RelativeOrientation")
                    ),
                    "zygosity": _text(_child(allele, "Zygosity")),
                    "clinical_significance": _text(
                        _child(significance, "Description")
                    )
                    or _text(significance),
                }
            )
        )
    return _nonempty(
        {
            "asserted_variant_zygosity": _text(_child(item, "Zygosity")),
            "count": _int_or_text(_text(_child(item, "Count"))),
            "cooccurring_alleles": alleles,
        }
    )


def _observation(item: ET.Element) -> dict[str, Any]:
    sample = _child(item, "Sample")
    data = []
    for observed_data in _children(item, "ObservedData"):
        attribute = _child(observed_data, "Attribute")
        if attribute is not None:
            data.append(
                _nonempty(
                    {
                        "type": attribute.get("Type", ""),
                        "value": _text(attribute),
                    }
                )
            )
    return _nonempty(
        {
            "origin": _text(_child(sample, "Origin")),
            "affected_status": _text(_child(sample, "AffectedStatus")),
            "sex": _text(_child(sample, "Sex")),
            "number_tested": _int_or_text(_text(_child(sample, "NumberTested"))),
            "number_males": _int_or_text(_text(_child(sample, "NumberMales"))),
            "number_females": _int_or_text(_text(_child(sample, "NumberFemales"))),
            "number_chromosomes_tested": _int_or_text(
                _text(_child(sample, "NumberChrTested"))
            ),
            "observed_data": data,
            "cooccurrence_sets": [
                _cooccurrence(value)
                for value in _children(item, "Co-occurrenceSet")
            ],
        }
    )


def _clinical_assertion(
    assertion: ET.Element,
    trait_mappings: list[dict[str, str]],
    linkage_method: str,
) -> dict[str, Any]:
    accession = _child(assertion, "ClinVarAccession")
    classification = _child(assertion, "Classification")
    germline = _child(classification, "GermlineClassification")
    observations = [
        _observation(item)
        for item in _children(_child(assertion, "ObservedInList"), "ObservedIn")
    ]
    zygosity_counts: Counter[str] = Counter()
    for observation in observations:
        for cooccurrence in observation.get("cooccurrence_sets", []):
            value = cooccurrence.get("asserted_variant_zygosity", "")
            if value:
                count = cooccurrence.get("count", 1)
                zygosity_counts[value] += count if isinstance(count, int) else 1
    return _nonempty(
        {
            "clinical_assertion_id": assertion.get("ID", ""),
            "scv_accession": accession.get("Accession", "") if accession is not None else "",
            "scv_version": _int_or_text(accession.get("Version")) if accession is not None else None,
            "record_status": _text(_child(assertion, "RecordStatus")),
            "contributes_to_aggregate_classification": (
                assertion.get("ContributesToAggregateClassification", "").lower()
                == "true"
            ),
            "submitter": _nonempty(
                {
                    "name": accession.get("SubmitterName", "") if accession is not None else "",
                    "organization_id": accession.get("OrgID", "") if accession is not None else "",
                    "organization_category": accession.get("OrganizationCategory", "") if accession is not None else "",
                }
            ),
            "classification": _nonempty(
                {
                    "germline_significance": _text(germline),
                    "review_status": _text(_child(classification, "ReviewStatus")),
                    "date_last_evaluated": classification.get(
                        "DateLastEvaluated", ""
                    )
                    if classification is not None
                    else "",
                }
            ),
            "submitted_mode_of_inheritance": _attribute_values(
                assertion, "ModeOfInheritance"
            ),
            "penetrance": _attribute_values(assertion, "Penetrance"),
            "trait_linkage_method": linkage_method,
            "trait_mappings": trait_mappings,
            "observed_zygosity_counts": dict(sorted(zygosity_counts.items())),
            "observations": observations,
        }
    )


def _eligible_conditions(classified: ET.Element, min_review_stars: int) -> list[dict[str, Any]]:
    eligible = []
    for rcv in _children(_child(classified, "RCVList"), "RCVAccession"):
        germline = _path(rcv, "RCVClassifications", "GermlineClassification")
        if germline is None:
            continue
        status = _text(_child(germline, "ReviewStatus"))
        stars = review_stars(status)
        if stars < min_review_stars:
            continue
        description = _child(germline, "Description")
        conditions = _condition_values(rcv)
        eligible.append(
            {
                "rcv_accession": rcv.get("Accession", ""),
                "rcv_version": _int_or_text(rcv.get("Version")),
                "title": rcv.get("Title", ""),
                "trait_set_id": _child(rcv, "ClassifiedConditionList").get(
                    "TraitSetID", ""
                )
                if _child(rcv, "ClassifiedConditionList") is not None
                else "",
                "conditions": conditions,
                "condition_tokens": _condition_tokens(conditions),
                "germline_classification": _nonempty(
                    {
                        "clinical_significance": _text(description),
                        "review_status": status,
                        "review_stars": stars,
                        "date_last_evaluated": description.get(
                            "DateLastEvaluated", ""
                        )
                        if description is not None
                        else "",
                        "submission_count": _int_or_text(
                            description.get("SubmissionCount")
                            if description is not None
                            else None
                        ),
                    }
                ),
            }
        )
    return eligible


def _attach_scvs(classified: ET.Element, conditions: list[dict[str, Any]]) -> None:
    mapping_by_id = _trait_mapping_by_assertion(classified)
    assertions = _children(_child(classified, "ClinicalAssertionList"), "ClinicalAssertion")
    for condition in conditions:
        condition_tokens = condition.pop("condition_tokens")
        matches = []
        unlinked = 0
        for assertion in assertions:
            if _text(_child(assertion, "RecordStatus")) != "current":
                continue
            if assertion.get("ContributesToAggregateClassification", "").lower() != "true":
                continue
            classification = _child(assertion, "Classification")
            if _child(classification, "GermlineClassification") is None:
                continue
            assertion_id = assertion.get("ID", "")
            mapping_rows = mapping_by_id.get(assertion_id, [])
            mapping_tokens = _mapping_tokens(mapping_rows)
            trait_tokens = _trait_tokens(assertion)
            if condition_tokens & mapping_tokens:
                linkage_method = "trait_mapping_clinical_assertion_id"
            elif condition_tokens & trait_tokens:
                linkage_method = "scv_trait_identifier_or_name"
            else:
                unlinked += 1
                continue
            matches.append(
                _clinical_assertion(
                    assertion,
                    trait_mappings=mapping_rows,
                    linkage_method=linkage_method,
                )
            )
        condition["matched_scvs"] = matches
        condition["scv_linkage"] = {
            "method": "condition identifiers/names plus TraitMapping ClinicalAssertionID",
            "matched_current_contributing_germline_scvs": len(matches),
            "unlinked_current_contributing_germline_scvs": unlinked,
        }

    scv_to_rcvs: dict[str, set[str]] = defaultdict(set)
    for condition in conditions:
        rcv_accession = str(condition.get("rcv_accession", ""))
        for scv in condition.get("matched_scvs", []):
            key = str(
                scv.get("clinical_assertion_id") or scv.get("scv_accession") or ""
            )
            if key and rcv_accession:
                scv_to_rcvs[key].add(rcv_accession)
    for condition in conditions:
        for scv in condition.get("matched_scvs", []):
            key = str(
                scv.get("clinical_assertion_id") or scv.get("scv_accession") or ""
            )
            linked_rcvs = sorted(scv_to_rcvs.get(key, set()))
            if len(linked_rcvs) > 1:
                scv["trait_linkage_ambiguous_across_eligible_rcvs"] = True
                scv["linked_eligible_rcv_accessions"] = linked_rcvs


def _match_rows_for_allele(
    allele: ET.Element,
    exact_lookup: dict[tuple[str, str], list[GofcardsVariant]],
) -> tuple[dict[str, dict[str, Any]], list[tuple[str, str, dict[str, Any], str]]]:
    """Find the GoFCards variants sharing this ClinVar allele's coordinates.

    The route is attached here rather than stored in the index, because how a
    link was found is a property of the match, not of the variant. The same
    variant can be reached by coordinates on one archive and by GoFCards' own
    ClinVar identifier on another.
    """
    locations = _allele_locations(allele)
    rows: list[tuple[str, str, dict[str, Any], str]] = []
    seen: set[tuple[str, str]] = set()
    for assembly, location in locations.items():
        build = ASSEMBLY_TO_BUILD[assembly]
        key = location["vcf_key"]
        for symbol, variant_id, variant in exact_lookup.get((build, key), []):
            # The identifier is coordinate-derived, so two genes curated at the
            # same allele share it; the gene has to be part of the identity or
            # one of the two claims would be dropped before it can be judged.
            marker = (symbol, variant_id)
            if marker in seen:
                continue
            seen.add(marker)
            rows.append((symbol, variant_id, variant, "genomic_coordinates"))
    return locations, rows


def _make_match(
    archive: ET.Element,
    allele: ET.Element,
    classification_scope: str,
    component_context: str,
    locations: dict[str, dict[str, Any]],
    matched_rows: list[tuple[str, str, dict[str, Any], str]],
    conditions: list[dict[str, Any]],
) -> dict[str, Any]:
    genes = _allele_genes(allele)
    clinvar_symbols = {gene.get("symbol", "") for gene in genes if gene.get("symbol")}
    matched_symbols = {symbol for symbol, _vid, _v, _route in matched_rows if symbol}
    concordant = sorted(clinvar_symbols & matched_symbols)
    assignment_symbols = (
        concordant if clinvar_symbols else sorted(matched_symbols)
    )
    kept = (
        [row for row in matched_rows if row[0] in concordant]
        if clinvar_symbols
        else matched_rows
    )
    # Only what identifies the match: which variant, in which gene, found how.
    # Not a copy of the variant -- the injection step reads this to nest the
    # ClinVar record back onto the variant the identifier already names.
    retained_rows = [
        {"HGNC_Symbol": symbol, "gofcards_variant_id": variant_id, "link_route": route}
        for symbol, variant_id, _variant, route in kept
    ]
    discarded_symbols = sorted(
        {
            symbol
            for symbol, _vid, _v, _route in matched_rows
            if symbol and symbol not in assignment_symbols
        }
    )
    return {
        "gene_symbols": assignment_symbols,
        "ClinVar_VCV": {
            "mechanism_context": ["GOF"],
            "match": {
                "method": "exact_normalized_vcf_allele",
                "gene_concordance": bool(concordant),
                "clinvar_gene_symbols": sorted(clinvar_symbols),
                "exact_cache_gene_symbols": sorted(matched_symbols),
                "discarded_gene_discordant_symbols": discarded_symbols,
                "matched_gofcards_records": retained_rows,
            },
            "variation": _nonempty(
                {
                    "variation_id": archive.get("VariationID", ""),
                    "vcv_accession": archive.get("Accession", ""),
                    "vcv_version": _int_or_text(archive.get("Version")),
                    "variation_name": archive.get("VariationName", ""),
                    "variation_type": archive.get("VariationType", ""),
                    "record_type": archive.get("RecordType", ""),
                    "record_status": _text(_child(archive, "RecordStatus")),
                    "date_created": archive.get("DateCreated", ""),
                    "date_last_updated": archive.get("DateLastUpdated", ""),
                    "classification_scope": classification_scope,
                    "matched_component_context": component_context,
                    "component_allele_id": allele.get("AlleleID", ""),
                    "component_variation_id": allele.get("VariationID", ""),
                    "component_name": _text(_child(allele, "Name")),
                    "component_variant_type": _text(_child(allele, "VariantType")),
                    "canonical_spdi": _text(_child(allele, "CanonicalSPDI")),
                    "genes": genes,
                    "locations": locations,
                    "hgvs": _allele_hgvs(allele),
                }
            ),
            "condition_assertions": conditions,
            "allelic_requirement": {
                "value": "unresolved",
                "reason": (
                    "ClinVar observed zygosity and submitted mode of inheritance are "
                    "evidence, not a variant-condition causal allelic requirement."
                ),
            },
        },
    }


def stream_parse_clinvar_vcv(
    xml_path: Path,
    exact_lookup: dict[tuple[str, str], list[GofcardsVariant]],
    min_review_stars: int = 2,
    variation_id_lookup: dict[str, list[GofcardsVariant]] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """Stream VCV XML; retain archives that match a GoFCards variant.

    Two routes decide whether an archive is kept, tried in this order:

      1. ``variation_id_lookup`` -- GoFCards publishes its own ClinVar
         VariationID for roughly 61% of its catalogue. Where it exists it is an
         exact identifier and needs no coordinate comparison at all.
      2. ``exact_lookup`` -- normalized genomic allele on either build, for the
         remaining variants and for any ClinVar record GoFCards has not
         cross-referenced.

    Each retained record carries ``link_route`` saying which one found it.
    """
    matches: list[dict[str, Any]] = []
    stats: Counter[str] = Counter()
    with _open_xml(xml_path) as handle:
        context = ET.iterparse(handle, events=("start", "end"))
        _, root = next(context)
        for event, archive in context:
            if event != "end" or archive.tag.rsplit("}", 1)[-1] != "VariationArchive":
                continue
            stats["variation_archives_seen"] += 1
            if stats["variation_archives_seen"] % 100000 == 0:
                print(
                    "clinvar_vcv_progress: "
                    f"variation_archives={stats['variation_archives_seen']}, "
                    f"exact_archives={stats['exact_variant_archives']}, "
                    f"retained_matches={stats['retained_component_matches']}",
                    file=sys.stderr,
                )
            status = _text(_child(archive, "RecordStatus"))
            classified = _child(archive, "ClassifiedRecord")
            if status != "current":
                stats["skipped_noncurrent"] += 1
            elif classified is None:
                stats["skipped_included_record"] += 1
            else:
                by_variation_id = (variation_id_lookup or {}).get(
                    str(archive.get("VariationID", "")).strip(), []
                )
                candidate_alleles = []
                for allele, scope, component_context in _variant_alleles(classified):
                    locations, rows = _match_rows_for_allele(allele, exact_lookup)
                    if by_variation_id:
                        # GoFCards named this exact ClinVar record, so keep it
                        # whether or not the coordinates also line up.
                        known = {(symbol, vid) for symbol, vid, _v, _r in rows}
                        rows = list(rows) + [
                            (symbol, vid, variant, "gofcards_variation_id")
                            for symbol, vid, variant in by_variation_id
                            if (symbol, vid) not in known
                        ]
                    if rows:
                        candidate_alleles.append(
                            (allele, scope, component_context, locations, rows)
                        )
                if by_variation_id:
                    stats["archives_matched_by_gofcards_variation_id"] += 1
                if candidate_alleles:
                    stats["exact_variant_archives"] += 1
                    conditions = _eligible_conditions(classified, min_review_stars)
                    if conditions:
                        _attach_scvs(classified, conditions)
                        for allele, scope, component_context, locations, rows in candidate_alleles:
                            match = _make_match(
                                archive,
                                allele,
                                scope,
                                component_context,
                                locations,
                                rows,
                                json.loads(json.dumps(conditions)),
                            )
                            if match.get("gene_symbols"):
                                matches.append(match)
                                stats["retained_component_matches"] += 1
                                stats[f"retained_{scope}"] += 1
                            else:
                                stats[
                                    "skipped_gene_discordant_component_matches"
                                ] += 1
                        stats["eligible_rcv_conditions"] += len(conditions)
                    else:
                        stats["exact_archives_without_eligible_rcv"] += 1
            archive.clear()
            root.clear()
    stats["retained_gene_assignments"] = sum(
        len(match.get("gene_symbols", [])) for match in matches
    )
    return matches, dict(sorted(stats.items()))


def flatten_match_audit(matches: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for match in matches:
        payload = match["ClinVar_VCV"]
        variation = payload["variation"]
        for symbol in match.get("gene_symbols", []):
            for condition in payload.get("condition_assertions", []):
                classification = condition.get("germline_classification", {})
                rows.append(
                    {
                        "gene_symbol": symbol,
                        "vcv_accession": variation.get("vcv_accession", ""),
                        "variation_id": variation.get("variation_id", ""),
                        "variation_name": variation.get("variation_name", ""),
                        "classification_scope": variation.get("classification_scope", ""),
                        "matched_component_context": variation.get("matched_component_context", ""),
                        "rcv_accession": condition.get("rcv_accession", ""),
                        "condition": "|".join(
                            item.get("name", "")
                            for item in condition.get("conditions", [])
                            if item.get("name")
                        ),
                        "clinical_significance": classification.get(
                            "clinical_significance", ""
                        ),
                        "review_status": classification.get("review_status", ""),
                        "review_stars": classification.get("review_stars", ""),
                        "matched_scv_count": len(condition.get("matched_scvs", [])),
                    }
                )
    return rows
