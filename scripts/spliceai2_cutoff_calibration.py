#!/usr/bin/env python3
"""Calibrate SpliceAI2 cutoffs against SpliceAI v1 on public ClinVar SNVs (aggregate counts only).

spliceai2_cutoff_calibration.py CLINVAR_VCF V1_SNV_VCF SPLICEAI2_SNV_VCF OUT_PREFIX

The values in spliceai_cutoffs.py come from the first five rows of OUT_PREFIX.cutoffs.tsv
(no_protein_change, 1+ stars), run on the 2026-10-19 GRCh38 ClinVar release with
bcftools on PATH. The run reads about 75 GB of score files and needs about 30 minutes
with 12 CPUs.

Steps, in time order:
 1. Known variants: ClinVar single-nucleotide variants classified pathogenic or
    likely pathogenic ("pathogenic") or benign or likely benign ("benign"), with
    at least one review star. Each goes into one group by ClinVar's molecular
    consequence:
      no_protein_change  only intron_variant and/or synonymous_variant, so a
                         pathogenic one most likely acts through splicing;
      canonical_splice   splice_donor_variant or splice_acceptor_variant.
 2. Scores: read with bcftools from each score file at those positions. For each
    allele the record of the ClinVar gene is used when present, otherwise the
    record with the largest delta score.
 3. Only variants scored by both files are compared, so both versions are judged
    on the same variants.
 4. For each PriVA decision, the SpliceAI2 cutoff that makes the mistake PriVA
    must avoid as often as SpliceAI v1 does at its current cutoff:
      benign_flagged               share of benign variants that pass the cutoff
                                   (for decisions that add pathogenic weight);
      pathogenic_called_no_impact  share of pathogenic variants called "no
                                   splicing impact" (for the benign decision).
    A 95% interval comes from 1,000 bootstrap resamples of the variants. The
    cutoff that best reproduces SpliceAI v1's own yes/no calls (highest F1) is
    reported beside it.
Writes OUT_PREFIX.coverage.tsv, OUT_PREFIX.curves.tsv and OUT_PREFIX.cutoffs.tsv;
prints counts and the cutoff table. Never prints variant rows.
"""
from concurrent.futures import ThreadPoolExecutor
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

PATHOGENIC = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
BENIGN = {"Benign", "Likely_benign", "Benign/Likely_benign"}
STARS = {"criteria_provided,_single_submitter": 1,
         "criteria_provided,_multiple_submitters,_no_conflicts": 2,
         "reviewed_by_expert_panel": 3, "practice_guideline": 4}
NO_PROTEIN_CHANGE = {"intron_variant", "synonymous_variant"}
CANONICAL = {"splice_donor_variant", "splice_acceptor_variant"}
SCORES = {"largest": [2, 3, 4, 5], "acceptor_gain": [2], "donor_loss": [5]}
# Every PriVA comparison of a score from the selected SpliceAI files with a fixed cutoff:
# (decision, score, comparison, SpliceAI v1 cutoff, mistake held at the v1 rate, where)
DECISIONS = [
    ("splice_altering", "largest", ">", 0.8, "benign_flagged",
     "splicing_var_analysis.py SpliceAI_interpretation cutoff"),
    ("nearby_acceptor_rescues_lost_acceptor", "acceptor_gain", ">", 0.5, "agreement_with_v1",
     "splicing_var_analysis.py SpliceAI_interpretation DS_AG > 0.5"),
    ("strong_enough_to_compare_with_known_pathogenic_splice_variant", "largest", ">=", 0.5, "benign_flagged",
     "acmg_ps1_pm5_known_variant.py PS1"),
    ("no_splicing_impact", "largest", "<", 0.1, "pathogenic_called_no_impact",
     "acmg_pp3_bp4_bp7_insilico.py BP4"),
    ("donor_loss_counted_as_putative_nmd", "donor_loss", ">", 0.5, "benign_flagged",
     "clinvar_stat_variants.py gene statistics (installation time)"),
]
GRID = np.arange(101)
BOOTSTRAPS = 1000


def known_variants(clinvar):
    command = ["bcftools", "query", "-i", 'INFO/CLNVC="single_nucleotide_variant"', "-f",
               "%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/MC\t%INFO/GENEINFO\n", clinvar]
    table = pd.read_table(subprocess.Popen(command, stdout=subprocess.PIPE).stdout, header=None, dtype=str,
                          names=["chrom", "pos", "ref", "alt", "clnsig", "revstat", "mc", "geneinfo"])
    print("ClinVar single-nucleotide variants: %d" % len(table))
    table["label"] = np.where(table["clnsig"].isin(PATHOGENIC), "pathogenic",
                              np.where(table["clnsig"].isin(BENIGN), "benign", ""))
    table["stars"] = table["revstat"].map(STARS).fillna(0).astype(int)

    def group(mc):
        terms = {item.split("|")[-1] for item in mc.split(",")} if mc not in (".", "") else set()
        if terms & CANONICAL:
            return "canonical_splice"
        if terms and terms <= NO_PROTEIN_CHANGE:
            return "no_protein_change"
        return ""
    table["group"] = table["mc"].fillna("").map(group)
    table = table[(table["label"] != "") & (table["stars"] >= 1) & (table["group"] != "")
                  & (table["chrom"] != "chrM")].copy()
    table["genes"] = table["geneinfo"].fillna("").map(
        lambda text: frozenset(item.split(":")[0] for item in text.split("|") if item and item != "."))
    table["pos"] = table["pos"].astype(int)
    return table.drop(columns=["clnsig", "revstat", "mc", "geneinfo"]).reset_index(drop=True)


def query_chromosome(path, chrom, positions, strip_chr, folder, label):
    region_file = os.path.join(folder, "%s_%s.tsv" % (label, chrom))
    name = chrom[3:] if strip_chr and chrom.startswith("chr") else chrom
    with open(region_file, "w") as stream:
        for pos in positions:
            stream.write("%s\t%d\n" % (name, pos))
    output = subprocess.run(["bcftools", "query", "-R", region_file, "-f",
                             "%POS\t%REF\t%ALT\t%INFO/SpliceAI\n", path],
                            stdout=subprocess.PIPE, text=True, check=True).stdout
    rows = []
    for line in output.splitlines():
        pos, ref, alt, info = line.split("\t")
        for entry in info.split(","):
            parts = entry.split("|")
            if parts[0] != alt or len(parts) < 6:
                continue
            try:
                deltas = [int(round(float(parts[i]) * 100)) for i in range(2, 6)]
            except ValueError:
                continue
            rows.append((chrom, int(pos), ref, alt, parts[1], max(deltas), deltas[0], deltas[3]))
    return rows


def read_scores(path, known, strip_chr, label, threads, folder):
    positions = known.groupby("chrom")["pos"].apply(lambda values: sorted(set(values)))
    with ThreadPoolExecutor(threads) as pool:
        parts = pool.map(lambda item: query_chromosome(path, item[0], item[1], strip_chr, folder, label),
                         positions.items())
        rows = [row for part in parts for row in part]
    scores = pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "symbol",
                                         "largest", "acceptor_gain", "donor_loss"])
    # One record per allele: the ClinVar gene's record if present, else the largest delta.
    merged = known[["chrom", "pos", "ref", "alt", "genes"]].merge(scores, on=["chrom", "pos", "ref", "alt"])
    merged["gene_match"] = [symbol in genes for symbol, genes in zip(merged["symbol"], merged["genes"])]
    merged = merged.sort_values(["gene_match", "largest"], ascending=False)
    chosen = merged.drop_duplicates(["chrom", "pos", "ref", "alt"])
    print("%s: %d score records at known positions; %d known variants scored, %d by the ClinVar gene's record"
          % (label, len(scores), len(chosen), int(chosen["gene_match"].sum())))
    return chosen.set_index(["chrom", "pos", "ref", "alt"])[list(SCORES)].add_prefix(label + "_")


def passing(counts, comparison):
    """Number of variants passing each cutoff 0..100 (scores are in hundredths)."""
    at_or_above = counts[::-1].cumsum()[::-1]
    if comparison == ">=":
        return at_or_above
    if comparison == ">":
        return np.append(at_or_above[1:], 0)
    return counts.sum() - at_or_above


def rates(counts, comparison):
    total = counts.sum()
    return passing(counts, comparison) / total if total else np.full(len(GRID), np.nan)


def rate(counts, comparison, cutoff):
    return rates(counts, comparison)[cutoff]


def choose(v1_counts, v2_counts, comparison, v1_cutoff, mistake):
    """SpliceAI2 cutoff holding the mistake rate at or below SpliceAI v1's at its cutoff."""
    side = "benign" if mistake == "benign_flagged" else "pathogenic"
    target = rate(v1_counts[side], comparison, v1_cutoff)
    curve = rates(v2_counts[side], comparison)
    allowed = GRID[curve <= target + 1e-12]
    if not len(allowed):
        return np.nan, target
    # Decisions that add pathogenic weight pass fewer benign variants as the cutoff rises:
    # take the lowest allowed cutoff. "No impact" (score below cutoff) calls fewer
    # pathogenic variants as the cutoff falls: take the highest allowed cutoff.
    return (allowed.min() if mistake == "benign_flagged" else allowed.max()), target


def agreement_cutoff(v1_scores, v2_scores, comparison, v1_cutoff):
    """SpliceAI2 cutoff whose yes/no calls best match SpliceAI v1's at its cutoff (highest F1)."""
    reference = compare(v1_scores, comparison, v1_cutoff)
    calls = passing(histogram(v2_scores), comparison)
    both = passing(histogram(v2_scores[reference]), comparison)
    denominator = reference.sum() + calls
    f1 = np.where(denominator > 0, 2 * both / np.maximum(denominator, 1), 0.0)
    best = int(np.argmax(f1))
    return best, f1[best]


def compare(scores, comparison, cutoff):
    return scores > cutoff if comparison == ">" else scores >= cutoff if comparison == ">=" else scores < cutoff


def histogram(values):
    return np.bincount(values, minlength=101)


def main():
    clinvar, v1_path, v2_path, prefix = sys.argv[1:5]
    threads = int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))
    known = known_variants(clinvar)
    print("known variants kept: %d" % len(known))
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(prefix))) as folder:
        v1 = read_scores(v1_path, known, True, "v1", threads, folder)
        v2 = read_scores(v2_path, known, False, "spliceai2", threads, folder)
    table = known.set_index(["chrom", "pos", "ref", "alt"]).join(v1).join(v2)
    table["scored_v1"] = table["v1_largest"].notna()
    table["scored_spliceai2"] = table["spliceai2_largest"].notna()

    table["scored_both"] = table["scored_v1"] & table["scored_spliceai2"]
    coverage = table.groupby(["group", "label"]).agg(
        known=("stars", "size"), scored_v1=("scored_v1", "sum"), scored_spliceai2=("scored_spliceai2", "sum"),
        scored_both=("scored_both", "sum"))
    coverage.to_csv(prefix + ".coverage.tsv", sep="\t")
    print(coverage.to_string())

    common = table[table["scored_both"]].copy()
    for column in [c for c in common if c.startswith(("v1_", "spliceai2_"))]:
        common[column] = common[column].astype(int)
    rng = np.random.default_rng(20261029)
    curves, cutoffs = [], []
    selections = [("no_protein_change", 1), ("no_protein_change", 2), ("canonical_splice", 1), ("both_groups", 1)]
    for group, min_stars in selections:
        subset = common[(common["stars"] >= min_stars)
                        & ((common["group"] == group) if group != "both_groups" else True)]
        by_label = {label: subset[subset["label"] == label] for label in ("benign", "pathogenic")}
        for score in SCORES:
            for comparison in (">", ">=", "<"):
                for source in ("v1", "spliceai2"):
                    for label, rows in by_label.items():
                        curve = rates(histogram(rows[source + "_" + score].to_numpy()), comparison)
                        curves.extend((group, min_stars, score, comparison, source, label, cutoff / 100, value)
                                      for cutoff, value in zip(GRID, curve))
        for decision, score, comparison, v1_cutoff, mistake, where in DECISIONS:
            c1 = int(round(v1_cutoff * 100))
            counts = {source: {label: histogram(rows[source + "_" + score].to_numpy())
                               for label, rows in by_label.items()} for source in ("v1", "spliceai2")}
            agree, f1 = agreement_cutoff(subset["v1_" + score].to_numpy(), subset["spliceai2_" + score].to_numpy(),
                                         comparison, c1)
            if mistake == "agreement_with_v1":
                chosen, target, interval = agree, np.nan, []
                for _ in range(BOOTSTRAPS):
                    pick = rng.integers(0, len(subset), len(subset))
                    interval.append(agreement_cutoff(subset["v1_" + score].to_numpy()[pick],
                                                     subset["spliceai2_" + score].to_numpy()[pick], comparison, c1)[0])
            else:
                chosen, target = choose(counts["v1"], counts["spliceai2"], comparison, c1, mistake)
                interval = []
                for _ in range(BOOTSTRAPS):
                    resampled = {}
                    for label, rows in by_label.items():
                        pick = rng.integers(0, len(rows), len(rows))
                        resampled[label] = (rows["v1_" + score].to_numpy()[pick],
                                            rows["spliceai2_" + score].to_numpy()[pick])
                    interval.append(choose({label: histogram(pair[0]) for label, pair in resampled.items()},
                                           {label: histogram(pair[1]) for label, pair in resampled.items()},
                                           comparison, c1, mistake)[0])
            low, high = np.nanpercentile(np.array(interval, dtype=float), [2.5, 97.5])
            row = {"decision": decision, "where": where, "group": group, "min_stars": min_stars,
                   "benign": len(by_label["benign"]), "pathogenic": len(by_label["pathogenic"]),
                   "comparison": comparison, "v1_cutoff": v1_cutoff, "mistake_held": mistake,
                   "v1_mistake_rate": target, "spliceai2_cutoff": chosen / 100,
                   "spliceai2_cutoff_95pct_low": low / 100, "spliceai2_cutoff_95pct_high": high / 100,
                   "spliceai2_cutoff_matching_v1_calls": agree / 100, "matching_f1": round(f1, 4)}
            for source, cutoff in (("v1", c1), ("spliceai2", chosen)):
                for label in ("benign", "pathogenic"):
                    row["%s_%s_passing" % (source, label)] = (
                        rate(counts[source][label], comparison, int(cutoff)) if cutoff == cutoff else np.nan)
            cutoffs.append(row)
    pd.DataFrame(curves, columns=["group", "min_stars", "score", "comparison", "source", "label", "cutoff",
                                  "share_passing"]).to_csv(prefix + ".curves.tsv", sep="\t", index=False)
    result = pd.DataFrame(cutoffs)
    result.to_csv(prefix + ".cutoffs.tsv", sep="\t", index=False, float_format="%.4g")
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(result.drop(columns=["where"]).to_string(index=False, float_format=lambda x: "%.4g" % x))


if __name__ == "__main__":
    main()
