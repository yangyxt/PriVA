"""Call `priva task manifest` for tests that inspect the Nextflow resource set.

The program is `PRIVA_BINARY`, else `priva` on PATH. A refusal of the task
contract (exit 2) is raised as ContractError carrying the message.
"""
import json
import os
from pathlib import Path
import subprocess
import tempfile

import yaml

PRIORITY_KEYS = (
    "clinvar_aa_stat", "clinvar_splice_stat", "clinvar_pm5_index", "all_intolerant_domains",
    "alphamissense_intolerant_domains", "pm1_regions_pkl", "alphamissense_vcf",
    "alphamissense_tranx_domain_map", "clinvar_patho_af_stat", "clinvar_patho_exon_af_stat",
    "interpro_mapping_pickle", "clingen_map", "gene_dosage_sensitivity",
    "hpo_condition_mechanism_json", "gene_nonlof_mechanism_json",
)


class ContractError(Exception):
    pass


def program():
    return os.environ.get("PRIVA_BINARY") or "priva"


def manifest(config, source, mode="annotated", checksum=False):
    """The `resources` records of a manifest made from this configuration."""
    with tempfile.TemporaryDirectory(prefix="priva_task_cli_",
                                     dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR"))) as folder:
        config_path = Path(folder) / "legacy.yaml"
        config_path.write_text(yaml.safe_dump(config, width=float("inf")))
        output = Path(folder) / "resources.json"
        command = [program(), "task", "manifest", "--config", str(config_path),
                   "--legacy-root", str(source), "--mode", "raw" if mode == "raw" else "annotated",
                   "--output", str(output)] + (["--checksum"] if checksum else [])
        result = subprocess.run(command, capture_output=True, text=True)
        marker = "PriVA compatibility contract failed: "
        if result.returncode == 2 and marker in result.stderr:
            raise ContractError(result.stderr.split(marker, 1)[1].strip())
        if result.returncode != 0:
            raise RuntimeError(f"{command} exit {result.returncode}: {result.stderr}")
        return json.loads(output.read_text())["resources"]


def runtime_resources(config, source, mode="annotated"):
    """Each resource key and the path the shell steps read."""
    return {key: Path(record["path"]) for key, record in manifest(config, source, mode).items()}
