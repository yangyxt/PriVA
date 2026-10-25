# PriVA scripts

Installation and execution follow this relationship:

`install_utils.sh main_install` → configured software/resources →
`annotation_vcf.sh` → `prioritization_vcf_per_fam.sh` → family outputs.

| Entry point or group | Role |
| --- | --- |
| `install_utils.sh` | Dispatch `main_install` or one named installation function with a configuration file |
| `installer_config_paths.py`, `installer_prerequisites.py`, `cache_resources.py` | Resolve storage paths, check prerequisites and describe resource lifecycles |
| Reference, VEP, BioMart, annotation-cache and CADD helpers | Acquire/build resources called by the corresponding installer functions |
| `annotation_vcf.sh` | Shared official VEP annotation and the configured CADD route |
| `filtration_vcf_per_fam.sh`, `prioritization_vcf_per_fam.sh` | Family selection, evidence assignment and ranking |
| `acmg_*.py`, `gene_mechanism_*.py` | Scientific runtime calculations and mechanism interpretation |

Resource destinations come from `config.yaml`; generated caches, logs and run
outputs belong in configured external storage. Source tests and their fixed
inputs are under `tests/`. See the [user manual](../README.md) for prerequisites,
commands, scientific-profile requirements and output meanings.
