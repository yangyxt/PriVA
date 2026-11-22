import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_mondo_disease_scope import (  # noqa: E402
    annotate_hpo_assertions,
    build_disease_scope_registry,
)


def _write_mondo(tmp_path: Path) -> Path:
    mondo = tmp_path / "mondo.obo"
    mondo.write_text(
        """format-version: 1.2
data-version: mondo/releases/test/mondo-simple.owl

[Term]
id: MONDO:0000001
name: disease

[Term]
id: MONDO:0003847
name: hereditary disease
is_a: MONDO:0000001 ! disease

[Term]
id: MONDO:0005070
name: neoplasm
is_a: MONDO:0000001 ! disease

[Term]
id: MONDO:1000001
name: inherited test disease
xref: OMIM:1 {source="MONDO:equivalentTo"}
is_a: MONDO:0003847 ! hereditary disease

[Term]
id: MONDO:1000002
name: hereditary cancer test disease
xref: OMIM:2 {source="MONDO:equivalentTo"}
is_a: MONDO:0003847 ! hereditary disease
is_a: MONDO:0005070 ! neoplasm

[Term]
id: MONDO:1000003
name: neoplasm with unresolved etiology
xref: Orphanet:3 {source="MONDO:equivalentTo"}
is_a: MONDO:0005070 ! neoplasm

[Term]
id: MONDO:1000004
name: susceptibility candidate
subset: omim_susceptibility
xref: OMIM:4 {source="MONDO:equivalentTo"}
is_a: MONDO:0000001 ! disease

[Term]
id: MONDO:1000006
name: inherited susceptibility without neoplasm ancestry
subset: omim_susceptibility
xref: OMIM:6 {source="MONDO:equivalentTo"}
is_a: MONDO:0003847 ! hereditary disease
""",
        encoding="utf-8",
    )
    return mondo


def _write_assertions(tmp_path: Path) -> Path:
    assertions = pd.DataFrame(
        [
            ["GENE1", "OMIM:1", "HP:0000006", "-", "TAS", "OMIM:1"],
            ["GENE2", "OMIM:2", "HP:0003829", "-", "PCS", "PMID:2"],
            ["GENE3", "ORPHA:3", "HP:0003584", "1/2", "PCS", "PMID:3"],
            ["GENE4", "OMIM:4", "HP:0010982", "-", "TAS", "OMIM:4"],
            ["GENE5", "OMIM:5", "HP:0001250", "1/1", "PCS", "PMID:5"],
            ["GENE6", "OMIM:6", "HP:0001250", "1/1", "PCS", "PMID:6"],
        ],
        columns=[
            "gene_symbol",
            "disease_id",
            "hpo_id",
            "frequency",
            "evidence",
            "reference",
        ],
    )
    path = tmp_path / "assertions.tsv"
    assertions.to_csv(path, sep="\t", index=False)
    return path


def test_builds_conservative_mondo_hpo_scope_registry(tmp_path: Path) -> None:
    registry = build_disease_scope_registry(
        _write_mondo(tmp_path),
        _write_assertions(tmp_path),
    ).set_index("disease_id")

    assert registry.loc["OMIM:1", "disease_scope"] == "mendelian_non_neoplastic"
    assert registry.loc["OMIM:1", "priva_scope"] == "include"
    assert registry.loc["OMIM:2", "disease_scope"] == "mendelian_cancer_predisposition"
    assert registry.loc["OMIM:2", "priva_scope"] == "include"
    assert registry.loc["ORPHA:3", "mondo_id"] == "MONDO:1000003"
    assert registry.loc["ORPHA:3", "disease_scope"] == "neoplastic_uncertain"
    assert registry.loc["ORPHA:3", "priva_scope"] == "review"
    assert registry.loc["OMIM:4", "disease_scope"] == "complex_or_non_monogenic"
    assert registry.loc["OMIM:4", "priva_scope"] == "exclude"
    assert registry.loc["OMIM:5", "disease_scope"] == "uncertain"
    assert (
        registry.loc["OMIM:6", "disease_scope"]
        == "mendelian_susceptibility_uncertain"
    )
    assert registry.loc["OMIM:6", "priva_scope"] == "review"


def test_manual_override_can_confirm_somatic_only_disease(tmp_path: Path) -> None:
    overrides = tmp_path / "overrides.tsv"
    overrides.write_text(
        "disease_id\tdisease_scope\tpriva_scope\tscope_review_status\t"
        "scope_evidence\tscope_reference\tscope_note\n"
        "ORPHA:3\tsomatic_oncogenic\texclude\tmanually_confirmed\t"
        "manual_review\tPMID:99\tSomatic-only condition\n",
        encoding="utf-8",
    )
    registry = build_disease_scope_registry(
        _write_mondo(tmp_path),
        _write_assertions(tmp_path),
        manual_overrides=overrides,
    ).set_index("disease_id")

    assert registry.loc["ORPHA:3", "disease_scope"] == "somatic_oncogenic"
    assert registry.loc["ORPHA:3", "priva_scope"] == "exclude"
    assert registry.loc["ORPHA:3", "scope_review_status"] == "manually_confirmed"


def test_scope_is_attached_to_every_assertion(tmp_path: Path) -> None:
    assertion_path = _write_assertions(tmp_path)
    assertions = pd.read_csv(assertion_path, sep="\t", dtype=str)
    registry = build_disease_scope_registry(_write_mondo(tmp_path), assertion_path)

    annotated = annotate_hpo_assertions(assertions, registry)

    assert len(annotated) == len(assertions)
    assert annotated["priva_scope"].isna().sum() == 0
    assert set(annotated["priva_scope"]) == {"include", "exclude", "review"}


def test_failed_second_publication_is_rebuilt_on_warm_retry(tmp_path):
    import json
    import os
    import subprocess
    from unittest import TestCase, mock
    import build_mondo_disease_scope as subject
    mondo = _write_mondo(tmp_path)
    raw = tmp_path / "raw"
    raw.mkdir()
    assertions = _write_assertions(raw)
    frame = pd.read_csv(assertions, sep="\t")
    base = raw / "genes_to_phenotype.assertions.base.tsv.gz"
    frame.to_csv(base, sep="\t", index=False)
    registry, scoped = tmp_path / "registry.tsv.gz", tmp_path / "scoped.tsv.gz"
    subject.build_and_write_resources(mondo, base, registry, annotated_assertions_output=scoped)
    old_scoped = scoped.read_bytes()
    frame.loc[len(frame)] = ["GENE7", "OMIM:7", "HP:0000006", "-", "TAS", "OMIM:7"]
    frame.to_csv(base, sep="\t", index=False)
    actual_replace = os.replace
    def fail_second(source, destination):
        if Path(destination) == scoped:
            raise OSError("injected second publication failure")
        return actual_replace(source, destination)
    with mock.patch.object(subject.os, "replace", fail_second):
        with TestCase().assertRaisesRegex(OSError, "second publication"):
            subject.build_and_write_resources(mondo, base, registry, annotated_assertions_output=scoped)
    assert scoped.read_bytes() == old_scoped
    assert len(pd.read_csv(registry, sep="\t")) == 7
    assert len(pd.read_csv(scoped, sep="\t")) == 6
    for name in ("genes_to_phenotype.txt", "phenotype.hpoa"):
        path = raw / name
        path.write_text("local raw fixture unused by warm route\n")
        os.utime(path, ns=(1, 1))
    config = tmp_path / "config.json"
    selected = {"hpo_raw_dir": str(raw), "hpo_assertions": str(scoped), "mondo_obo": str(mondo),
                "mondo_disease_scope_registry": str(registry),
                "mondo_disease_scope_overrides": str(tmp_path / "absent-overrides.tsv")}
    config.write_text(json.dumps(selected))
    script = r'''
source "$1"
read_yaml() { "$TEST_PYTHON" -c 'import json,sys;print(json.load(open(sys.argv[1])).get(sys.argv[2],""))' "$@"; }
update_yaml() { "$TEST_PYTHON" -c 'import json,sys;p,k,v=sys.argv[1:];d=json.load(open(p));d[k]=v;open(p,"w").write(json.dumps(d))' "$@"; }
download_resource_atomic() { echo UNEXPECTED_DOWNLOAD >&2; return 98; }
wget() { echo UNEXPECTED_DOWNLOAD >&2; return 98; }
python() {
    if [[ "${FORBID_BUILD:-0}" == 1 && "$1" == */build_mondo_disease_scope.py ]]; then
        echo UNEXPECTED_BUILD >&2; return 97
    fi
    "$TEST_PYTHON" "$@"
}
mondo_hpo_scope_install "$2"
'''
    env = {key: value for key, value in os.environ.items() if not key.startswith("PRIVA_FORCE_")}
    env.update(TEST_PYTHON=sys.executable, PYTHONDONTWRITEBYTECODE="1")
    command = ["bash", "-c", script, "fixture", str(ROOT / "scripts/install_utils.sh"), str(config)]
    retry = subprocess.run(command, env=env, text=True, capture_output=True, timeout=45)
    assert retry.returncode == 0, retry.stderr
    assert "UNEXPECTED" not in retry.stdout + retry.stderr
    assert len(pd.read_csv(registry, sep="\t")) == 7
    assert len(pd.read_csv(scoped, sep="\t")) == 7
    expected = subject.annotate_hpo_assertions(frame.astype(str), subject.build_disease_scope_registry(mondo, base))
    actual = pd.read_csv(scoped, sep="\t", dtype=str, keep_default_na=False)
    pd.testing.assert_frame_equal(actual, expected.astype(str))
    def identity(path):
        value = path.stat()
        # Validation reads may update atime; publication identity must not change.
        return value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns, value.st_mode
    before = {path: (path.read_bytes(), identity(path)) for path in (registry, scoped, config)}
    warm = subprocess.run(command, env=dict(env, FORBID_BUILD="1"), text=True, capture_output=True, timeout=45)
    assert warm.returncode == 0, warm.stderr
    assert "UNEXPECTED" not in warm.stdout + warm.stderr
    for path in (registry, scoped):
        assert (path.read_bytes(), identity(path)) == before[path]
    assert config.read_bytes() == before[config][0]
