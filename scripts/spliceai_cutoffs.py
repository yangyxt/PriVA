"""SpliceAI cutoffs for PriVA's splicing and ACMG rules, by score version.

PriVA's rules were written for SpliceAI v1 scores. SpliceAI2 (Jaganathan et al.,
"A unified framework for quantitative splicing and transcript prediction", 2026)
reports the same four delta scores, and its authors recommend their own
thresholds on the largest of them (their summary score). Their README maps each
threshold to the SpliceAI v1 tier it replaces:

  tier                      SpliceAI v1   SpliceAI2
  high recall                  0.2          0.1
  balance of both              0.5          0.25
  high precision               0.8          0.5

The manuscript calls variants scoring above 0.1 predicted splice-altering and
prioritises disease variants at 0.25. Each PriVA rule takes the SpliceAI2
threshold of the tier its SpliceAI v1 cutoff belongs to:

  rule (where)                                       v1    SpliceAI2
  splice_altering (SpliceAI_interpretation, PVS1/PP3) 0.8   0.5   high precision
  nearby_acceptor_rescue (SpliceAI_interpretation)    0.5   0.25  balance
  known_pathogenic_comparison (PS1)                   0.5   0.25  balance
  no_splicing_impact (BP4; score below the cutoff)    0.1   0.1   not predicted
                                                                  splice-altering
  donor_loss_putative_nmd (ClinVar gene statistics)   0.5   0.25  balance

The version comes from the headers of the configured spliceai_snv_prescore and
spliceai_indel_prescore files, which must be the same version and match the
configured assembly (SpliceAI2 is GRCh38 only). Indels use the same cutoffs.

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
        "splice_altering": 0.5,
        "nearby_acceptor_rescue": 0.25,
        "known_pathogenic_comparison": 0.25,
        "no_splicing_impact": 0.1,
        "donor_loss_putative_nmd": 0.25,
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
