# PriVA scripts

Installation and execution run as Nextflow workflows; the scripts here are the
helpers their steps call and the Python originals of the `priva` commands:

`nextflow run main.nf --action install` (one task per step of
`conf/install_steps.json`) → configured software/resources →
`nextflow run main.nf` (annotation, filtration and prioritization tasks, which
run in `priva`) → family outputs.

| Entry point or group | Role |
| --- | --- |
| `bootstrap_nextflow_runtime.sh`, `bootstrap_priva_rust.py` | Install the pinned Nextflow driver and build the `priva` program before the installation starts |
| `installer_config_paths.py`, `installer_prerequisites.py` | Resolve storage paths and check prerequisites (the Python originals of `priva installer-config-paths` and `priva installer-prerequisites`) |
| Reference, VEP, BioMart, annotation-cache and CADD helpers | Acquire/build resources; install steps run them or their `priva` ports |
| `cadd_conda.py`, `cadd_scripts/` | Find or prepare the Conda that upstream `CADD.sh` needs, and PriVA's faster copies of two CADD scripts |
| `das_install_state.py`, `das_input_continuity.py` | Track completed DAS outputs; reuse after a ClinVar input change requires a separate exact consumer-equivalence proof |
| `acmg_*.py`, `gene_mechanism_*.py` | Scientific calculations and mechanism interpretation, the originals the Rust ACMG step was compared with |

The shell installer (`install_utils.sh`), the cluster scripts and the analysis
shell scripts (`annotation_vcf.sh`, `filtration_vcf_per_fam.sh`,
`prioritization_vcf_per_fam.sh`, `common_bash_utils.sh`) were removed; commit
841e200 is the last tree that has them.

Resource destinations come from `config.yaml`; generated caches, logs and run
outputs belong in configured external storage. Source tests and their fixed
inputs are under `tests/`. See the [user manual](../README.md) for prerequisites,
commands, scientific-profile requirements and output meanings.
