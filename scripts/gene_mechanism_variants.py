r"""Query-variant mechanism scoring.

Definitions here describe what the current allele plausibly does. Resource
loading and pipeline orchestration remain outside this module.

The mechanism decision tree is deliberately kept beside its implementation:

    upstream exact-allele annotation
      canonical non-LOF cache (currently eligible curated GoFCards alleles;
      ClinVar supplies condition links but pathogenicity alone is not GOF/DN)
      match route: HGNC symbol + HGVSp, then HGVSc, then normalized allele
                                |
                                v
    transcript-destroying LOF?  stop/frameshift/splice-frameshift + NMD
      yes -> LOF=2, GOF=0, DN=0                         [exclusive]
      no
       |
       v
    exact curated GOF or DN score=2?
      yes -> asserted mechanism(s)=2, all others=0      [exclusive]
      no
       |
       v
    predicted LOF evidence?
      LOFTEE HC/OS; VEP LOF consequence; PriVA splice/5UTR LOF;
      NMD-escaping truncation
      yes -> LOF=1 or 2; retain only LOF-compatible histories
      no  -> LOF=1, GOF=1, DN=1; retain every explicit mechanism history
             plus condition-level UNRESOLVED histories [parallel hypotheses]

Only a score of 2 is exclusive. Three parallel score-1 values mean that the
variant evidence cannot rule out any of PriVA's supported mechanisms; they are
not three positive molecular findings.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any, Callable

import pandas as pd

from gene_mechanism_common import (
    CURATED_EXACT_MECHANISMS,
    DEFAULT_GOFCARDS_EXACT_GOF_HGVSP,
    EXACT_SEQUENCE_MECHANISMS,
    VARIANT_MECHANISM_SCORE_KEYS,
    _bounded_mechanism_score,
    _bool_value,
    _clean,
    _exact_mechanisms_from_effect,
    _genomic_match_key,
    _norm_hgvsp,
    _normalize_assembly,
    _normalize_sequence_mechanism,
    _requested_genomic_key_types,
    _split_annotation_tokens,
)


def infer_query_variant_effect(row: dict[str, Any] | pd.Series) -> dict[str, Any]:
    """Infer one exclusive variant effect while retaining suppressed predictions.

    Exact curated GoF or dominant-negative evidence has score ``2`` and defines
    the allowed mechanism set for this allele. Generic consequence, NMD, or
    LOFTEE signals remain visible as audit evidence but cannot create a competing
    LoF call. Without an exact assertion, prediction-based LoF evidence receives
    score ``1``. LOFTEE ``OS`` is its other-splice category. If none of those
    observations identifies or excludes a mechanism, all three mechanisms score
    ``1`` in parallel. That is an uncertainty representation, not positive
    evidence for three simultaneous effects.
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

    curated_mechanisms = {
        mechanism for mechanism, score in scores.items() if score == 2
    } & CURATED_EXACT_MECHANISMS
    predicted_lof_evidence = list(dict.fromkeys(predicted_lof_evidence))
    suppressed_evidence: list[str] = []
    evidence: list[str] = []
    decay_lof = "NMD_PREDICTED_LOF" in predicted_lof_evidence

    if decay_lof:
        # Highest precedence, above even a curated allele. Nonsense-mediated
        # decay destroys the transcript, so there is no protein left to gain a
        # function or to poison a complex. A curated gain-of-function claim for
        # this position cannot survive the message being degraded, so it is
        # recorded as suppressed rather than allowed to stand.
        suppressed_evidence = [
            f"CANONICAL_EXACT_{mechanism}" for mechanism in sorted(curated_mechanisms)
        ]
        scores = {mechanism: 0 for mechanism in scores}
        scores["LOF"] = 2
        evidence.extend(predicted_lof_evidence)
        effect = "exact_known_LOF"
    elif curated_mechanisms:
        # A curated gain-of-function or dominant-negative allele outranks every
        # remaining loss-of-function signal, including a LOFTEE rescue: those
        # are predictions, this is a curator's call on this exact allele.
        for mechanism in scores:
            scores[mechanism] = 2 if mechanism in curated_mechanisms else 0
        evidence.extend(
            f"CANONICAL_EXACT_{mechanism}" for mechanism in sorted(curated_mechanisms)
        )
        suppressed_evidence = predicted_lof_evidence
        effect = "exact_known_" + "+".join(sorted(curated_mechanisms))
    elif predicted_lof_evidence:
        # No decay and no curated allele. LOFTEE's high-confidence call still
        # establishes the loss; anything weaker leaves the transcript possibly
        # intact, so the loss is plausible rather than established.
        scores["LOF"] = max(
            scores["LOF"], 2 if "LOFTEE_HC" in predicted_lof_evidence else 1
        )
        evidence.extend(predicted_lof_evidence)
        effect = "predicted_LOF_high_confidence"
    else:
        # Absence of a mechanism-specific observation cannot make all known
        # condition mechanisms incompatible. Parallel score-1 hypotheses retain
        # LOF, GOF and dominant-negative histories until stronger variant evidence
        # selects one; condition-level UNRESOLVED histories are retained beside
        # them by select_condition_histories_for_variant().
        scores = {mechanism: max(score, 1) for mechanism, score in scores.items()}
        effect = "uncertain"
    exact_mechanisms = _exact_mechanisms_from_effect(effect)
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

class ExactVariantMatcher:
    """Lazy exact-allele indexes used by :class:`GeneMechanismHub`.

    The matcher owns variant normalization, route precedence, and match
    reporting. Resource resolution is supplied as callbacks, so this module
    never imports the resource service or the top-level hub facade.
    """

    def __init__(
        self,
        *,
        mechanism_json: str | Path,
        load_mechanisms: Callable[[], dict[str, dict[str, Any]]],
        try_resolve_symbol: Callable[[Any], str],
        resolved_symbol_key: Callable[[Any], str],
    ) -> None:
        self.mechanism_json = Path(mechanism_json)
        self._load_mechanisms_callback = load_mechanisms
        self._try_resolve_symbol_callback = try_resolve_symbol
        self._resolved_symbol_key_callback = resolved_symbol_key
        self._canonical_exact_nonlof_rows: list[dict[str, Any]] | None = None
        self._gofcards_by_symbol_hgvsp: (
            dict[tuple[str, str], list[dict[str, str]]] | None
        ) = None
        self._gofcards_by_symbol_hgvsc: (
            dict[tuple[str, str], list[dict[str, str]]] | None
        ) = None
        self._gofcards_by_symbol_genomic: (
            dict[tuple[str, str, str, str], list[dict[str, str]]] | None
        ) = None

    def _load_mechanisms(self) -> dict[str, dict[str, Any]]:
        return self._load_mechanisms_callback()

    def _try_resolve_symbol(self, gene_symbol: Any) -> str:
        return self._try_resolve_symbol_callback(gene_symbol)

    def _resolved_symbol_key(self, gene_symbol: Any) -> str:
        return self._resolved_symbol_key_callback(gene_symbol)

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
