#!/usr/bin/env python3
"""Check a GoFCards cache and the provenance it records.

The arguments are the cache path and any KEY=VALUE pairs. `priva
gofcards-cache-validate` is the installer's check; this script serves the main
preflight, which runs before priva_binary is recorded.
"""
import gzip
import json
import re
import sys

path = sys.argv[1]
expected_provenance = sys.argv[2:]
opener = gzip.open if path.endswith(".gz") else open
with opener(path, "rt", encoding="utf-8") as handle:
    cache = json.load(handle)

for key in ("metadata", "genes"):
    if key not in cache:
        raise SystemExit(f"{path} has no top-level {key!r} block")
meta = cache["metadata"]
if meta.get("source") != "GoFCards":
    raise SystemExit(f"{path} metadata does not declare GoFCards as the source")
if meta.get("mechanism") != "GOF":
    raise SystemExit(f"{path} metadata does not declare GOF as the mechanism")

for pair in expected_provenance:
    key, _, want = pair.partition("=")
    if not want:
        continue
    # A dotted key walks into a nested block, so a caller can ask about
    # "clinvar.injector_sha256" as easily as about a value recorded at the top.
    have = meta
    for part in key.split("."):
        have = have.get(part) if isinstance(have, dict) else None
    if str(have) != want:
        raise SystemExit(
            f"{path} was built with {key}={have!r}, but {want!r} is configured now"
        )

# loc_<chrom>:<start>:<ref>-><alt>_grch37
VARIANT_ID = re.compile(r"^loc_[0-9XYMT]+:\d+:[ACGTN-]+->[ACGTN-]+_grch37$")
ELIGIBLE = "eligible"
ALLOWED = {
    ELIGIBLE,
    "quarantined_gene_discordance",
    "quarantined_allele_gene_discordance",
    "quarantined_missing_source_gene",
    "quarantined_reviewed_lof",
    "quarantined_reviewed_mixed",
    "quarantined_reviewed_dominant_negative",
    "quarantined_reviewed_uncertain",
    "quarantined_reviewed_exclude",
    "quarantined_mechanism_review_required",
}
genes = variants = eligible = quarantined = with_protein = with_hg38 = evidence = 0
hgvsc_keys = hgvsp_keys = 0
for symbol, gene in cache["genes"].items():
    genes += 1
    if not symbol.strip():
        raise SystemExit(f"{path} has an empty gene symbol key")
    for vid, variant in (gene.get("variants") or {}).items():
        variants += 1
        if not VARIANT_ID.match(vid):
            raise SystemExit(f"{path}: {symbol}/{vid} is not a well-formed variant identifier")
        record = variant.get("record") or {}
        state = (record.get("eligibility") or {}).get("status", "")
        if state not in ALLOWED:
            raise SystemExit(f"{path}: {symbol}/{vid} has unknown eligibility {state!r}")
        assemblies = variant.get("assemblies") or {}
        if not assemblies:
            raise SystemExit(f"{path}: {symbol}/{vid} is present on no assembly")
        protein = coding = False
        for assembly, block in assemblies.items():
            if assembly not in ("hg19", "hg38"):
                raise SystemExit(f"{path}: {symbol}/{vid} has unknown assembly {assembly!r}")
            coords = block.get("genomic") or {}
            if not coords.get("chrom") or not isinstance(coords.get("pos"), int):
                raise SystemExit(f"{path}: {symbol}/{vid}/{assembly} has no usable coordinates")
            if not coords.get("ref") or not coords.get("alt"):
                raise SystemExit(f"{path}: {symbol}/{vid}/{assembly} has no reference or alternate")
            for transcript, view in (block.get("transcripts") or {}).items():
                by_c = view.get("by_hgvsc") or {}
                by_p = view.get("by_hgvsp") or {}
                hgvsc_keys += len(by_c)
                hgvsp_keys += len(by_p)
                if by_c: coding = True
                if by_p: protein = True
                for hgvsp, coding_keys in by_p.items():
                    if not isinstance(coding_keys, list) or not coding_keys:
                        raise SystemExit(f"{path}: {symbol}/{vid} protein key {hgvsp!r} maps to no coding change")
                    for ck in coding_keys:
                        if ck not in by_c:
                            raise SystemExit(f"{path}: {symbol}/{vid} protein key {hgvsp!r} points at unknown {ck!r}")
        if "hg38" in assemblies: with_hg38 += 1
        if protein: with_protein += 1
        if state == ELIGIBLE:
            eligible += 1
            if not (protein or coding or assemblies):
                raise SystemExit(f"{path}: {symbol}/{vid} is eligible but offers no match route")
        else:
            quarantined += 1
            if not (record.get("eligibility") or {}).get("reason"):
                raise SystemExit(f"{path}: {symbol}/{vid} is quarantined without a reason")
        evidence += len(record.get("evidence") or [])

if variants == 0:
    raise SystemExit(f"{path} contains no variants")
if eligible == 0:
    raise SystemExit(f"{path} contains no runtime-eligible variants")
if with_protein == 0:
    raise SystemExit(f"{path} contains no variant with a protein-change key")
print(
    f"gofcards_genes={genes}; variants={variants}; eligible={eligible}; "
    f"quarantined={quarantined}; with_protein_key={with_protein}; on_hg38={with_hg38}; "
    f"hgvsc_keys={hgvsc_keys}; hgvsp_keys={hgvsp_keys}; evidence_entries={evidence}"
)
