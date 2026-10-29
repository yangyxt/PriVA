"""SpliceAI cutoffs for PriVA's splicing and ACMG rules, by score version.

PriVA's rules were written for SpliceAI v1 scores. SpliceAI2 scores the same
variants lower, so each rule has its own SpliceAI2 cutoff. Each SpliceAI2 cutoff
makes the mistake its rule must avoid as often as SpliceAI v1 does at the
original cutoff. This was measured on the 988,870 benign and 3,988 pathogenic
ClinVar single-nucleotide variants (2026-10-19 release, 1+ review stars) that
change no amino acid (intronic or synonymous) and are scored by both versions.

  rule (where)                     mistake held at the v1 rate     v1    SpliceAI2
  splice_altering                  benign flagged: 0.059%          0.8   0.54
    (SpliceAI_interpretation, feeds PVS1/PP3)
  nearby_acceptor_rescue           no ClinVar truth exists; the    0.5   0.45
    (SpliceAI_interpretation)      cutoff best reproducing v1's
                                   own calls is used
  known_pathogenic_comparison      benign passing: 0.30%           0.5   0.26
    (PS1 splice comparison)
  no_splicing_impact               pathogenic called "no splicing  0.1   0.07
    (BP4)                          impact": 7.0%
  donor_loss_putative_nmd          benign counted: 0.058%          0.5   0.34
    (ClinVar gene statistics built at installation)

At these cutoffs SpliceAI2 also catches more pathogenic variants than SpliceAI
v1 (splice-altering: 63% versus 43%). 95% bootstrap intervals are within
0.02 of each value except nearby_acceptor_rescue (0.38-0.52). Reproduce with
spliceai2_cutoff_calibration.py, which reads only the public ClinVar VCF and the
two score files.

The version comes from the headers of the configured spliceai_snv_prescore and
spliceai_indel_prescore files, which must be the same version and match the
configured assembly (SpliceAI2 is GRCh38 only). Indels use the cutoffs measured
on single-nucleotide variants.

Command line: spliceai_cutoffs.py CONFIG [SNV_VCF INDEL_VCF]
prints the version and its cutoffs, or exits with an error when the files differ.
"""
from pathlib import Path
import sys

import yaml

from spliceai_cache_input import canonical_assembly, inspect_header

CUTOFFS = {
    "vcf_v1": {
        "splice_altering": 0.8,
        "nearby_acceptor_rescue": 0.5,
        "known_pathogenic_comparison": 0.5,
        "no_splicing_impact": 0.1,
        "donor_loss_putative_nmd": 0.5,
    },
    "vcf_spliceai2": {
        "splice_altering": 0.54,
        "nearby_acceptor_rescue": 0.45,
        "known_pathogenic_comparison": 0.26,
        "no_splicing_impact": 0.07,
        "donor_loss_putative_nmd": 0.34,
    },
}
VERSION_NAMES = {"vcf_v1": "SpliceAI v1", "vcf_spliceai2": "SpliceAI2"}


def score_version(snv=None, indel=None, assembly=None):
    """Return 'vcf_v1' or 'vcf_spliceai2' for the given score files.

    Missing or unset files are skipped so the annotation step can report them;
    with neither file present the original SpliceAI v1 cutoffs apply.
    """
    expected = canonical_assembly(assembly) if assembly not in (None, "", "null") else None
    found = {}
    for kind, value in (("SNV", snv), ("indel", indel)):
        if value in (None, "", "null") or not Path(value).is_file():
            continue
        header = inspect_header(value)
        if expected and header["assembly"] not in (None, expected):
            raise ValueError(f"The SpliceAI {kind} score file is {header['assembly']}, "
                             f"but the configured assembly is {expected}: {value}")
        found[kind] = header["format"]
    versions = set(found.values())
    if len(versions) > 1:
        raise ValueError("The SpliceAI SNV and indel score files are different versions: " +
                         ", ".join(f"{kind} {VERSION_NAMES.get(name, name)}" for kind, name in found.items()))
    version = next(iter(versions), "vcf_v1")
    if version not in CUTOFFS:
        raise ValueError("SpliceAI2 publisher tables must be converted to indexed VCF by the installer first")
    return version


def read_config(config_path):
    settings = {}
    if config_path and Path(config_path).is_file():
        settings = yaml.safe_load(Path(config_path).read_text()) or {}
    if not isinstance(settings, dict):
        raise ValueError("PriVA config must be a mapping")
    return settings


def cutoffs_for_config(config_path=None):
    """Cutoffs for the SpliceAI score version configured in config_path."""
    settings = read_config(config_path)
    return CUTOFFS[score_version(settings.get("spliceai_snv_prescore"), settings.get("spliceai_indel_prescore"),
                                 settings.get("assembly"))]


def main():
    if len(sys.argv) not in (2, 4):
        sys.exit(__doc__)
    try:
        settings = read_config(sys.argv[1])
        snv, indel = (sys.argv[2:4] if len(sys.argv) == 4 else
                      (settings.get("spliceai_snv_prescore"), settings.get("spliceai_indel_prescore")))
        version = score_version(snv, indel, settings.get("assembly"))
    except (OSError, EOFError, ValueError) as error:
        sys.exit("Cannot use the SpliceAI score files: " + str(error))
    print(VERSION_NAMES[version] + " cutoffs: " +
          ", ".join(f"{name} {value}" for name, value in CUTOFFS[version].items()))


if __name__ == "__main__":
    main()
