"""Final installer routing/status regressions; acquisition and scientific calls are stubs."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
import unittest


ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "scripts/install_utils.sh"
PRIORITIZATION = ROOT / "scripts/prioritization_vcf_per_fam.sh"
PRELUDE = r'''
source "$1"
DATA_DIR="$TEST_DATA"
log() { printf '%s\n' "$*" >&2; }
event() {
    "$TEST_PYTHON" -c 'import json,os,sys; open(os.environ["TEST_EVENTS"],"a").write(json.dumps(sys.argv[1:])+"\n")' "$@"
}
read_yaml() {
    "$TEST_PYTHON" -c 'import json,sys; v=json.load(open(sys.argv[1])).get(sys.argv[2]); print("null" if v is None else str(v).lower() if isinstance(v,bool) else v)' "$1" "$2"
}
update_yaml() {
    event update "$2" "$3"
    [[ "$TEST_FAIL" != "$2" ]] || return 71
    "$TEST_PYTHON" -c 'import json,sys; p,k,v=sys.argv[1:]; d=json.load(open(p)); d[k]=v; open(p,"w").write(json.dumps(d))' "$@"
}
update_or_append_yaml() { update_yaml "$@"; }
wget() {
    event download "$@"
    local destination
    while [[ $# -gt 0 ]]; do
        if [[ "$1" == -O ]]; then destination=$2; break; fi
        shift
    done
    cp "$TEST_FASTA" "$destination"
}
git() { event git "$@"; return 72; }
python() { "$TEST_PYTHON" "$@"; }
'''


class FinalInstallerRoutes(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="priva_final_routes_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.cache = self.folder / "cache"
        self.cache.mkdir()
        self.data = self.folder / "source_data"
        self.data.mkdir()
        self.config_file = self.folder / "config.json"
        self.events = self.folder / "events.jsonl"
        self.events.touch()
        self.fasta = self.folder / "public_fixture.fa.gz"
        self.write_fasta(self.fasta)
        self.config = {"assembly": "hg19", "vep_plugins_cachedir": str(self.cache), "vep_installed_version": "116"}
        self.env = dict(os.environ, TEST_PYTHON=sys.executable, TEST_EVENTS=str(self.events),
                        TEST_CONFIG=str(self.config_file), TEST_CACHE=str(self.cache), TEST_DATA=str(self.data),
                        TEST_FASTA=str(self.fasta), TEST_FAIL="none", TMPDIR=str(self.folder))

    @staticmethod
    def write_fasta(path, assembly="GRCh37", compressed=True):
        text = f">ENST00000000233.1 cds chromosome:{assembly}:1:1:6:1\nATGTAA\n"
        if compressed:
            with path.open("wb") as output:
                with gzip.GzipFile(fileobj=output, mode="wb", filename="", mtime=0) as handle:
                    handle.write(text.encode())
        else:
            path.write_text(text)

    def run_shell(self, command, setup="", fail="none"):
        self.config_file.write_text(json.dumps(self.config))
        overrides = ""
        previous = os.environ.get("PRIVA_TEST_BEFORE_SOURCE_ROOT")
        if previous:
            text = (Path(previous) / "scripts/install_utils.sh").read_text()
            for name in ("UTRAnnotator_install", "CDS_FASTA_install", "AlphaMissense_install", "MaveDB_install",
                         "pext_install", "gene_pathogenic_mechanism_cache_install", "gofcards_clinvar_injection_install",
                         "gene_nonlof_mechanism_cache_install", "hpo_condition_mechanism_cache_install", "mondo_hpo_scope_install"):
                overrides += re.search(r"^function " + name + r"\s*\(\)\s*\{.*?(?=^function |\Z)",
                                       text, re.M | re.S).group() + "\n"
        result = subprocess.run(["bash", "--noprofile", "--norc", "-c",
                                 PRELUDE + overrides + setup + "\n" + command, "test", str(INSTALLER)],
                                env=dict(self.env, TEST_FAIL=fail), cwd=self.folder,
                                capture_output=True, text=True, timeout=45)
        return result, [json.loads(line) for line in self.events.read_text().splitlines()]

    def test_utr_initial_clone_preserves_caller_cwd_and_relative_input(self):
        self.config["utr_annotator_file"] = ""
        selected = self.folder / "relative_input.vcf"
        selected.write_text("retained relative input\n")
        setup = r'''
git() {
    event clone "$PWD" "$@"
    mkdir -p UTRannotator
    printf '37 resource\n' > UTRannotator/uORF_5UTR_GRCh37_PUBLIC.txt
    printf '38 resource\n' > UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt
}
'''
        result, events = self.run_shell(r'''
original_pwd=$PWD
UTRAnnotator_install "$TEST_CONFIG" "$TEST_CACHE" || exit $?
[[ $PWD == "$original_pwd" && -s relative_input.vcf ]]
''', setup)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events[0][0], "clone")
        self.assertEqual(Path(events[0][1]).parent, self.cache)
        self.assertTrue(Path(events[0][1]).name.startswith(".UTRannotator.clone."))
        self.assertEqual(json.loads(self.config_file.read_text())["utr_annotator_file"],
                         str(self.cache / "UTRannotator/uORF_5UTR_GRCh37_PUBLIC.txt"))

    def test_cds_explicit_gzip_and_plain_inputs_are_preserved_without_managed_acquisition(self):
        for compressed in (True, False):
            with self.subTest(compressed=compressed):
                selected = self.folder / ("selected input's.fa.gz" if compressed else "selected input's.fa")
                self.write_fasta(selected, compressed=compressed)
                before = selected.read_bytes()
                self.config["cds_fasta_file"] = str(selected)
                self.events.write_text("")
                result, events = self.run_shell('CDS_FASTA_install "$TEST_CONFIG"')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(events, [])
                self.assertEqual(json.loads(self.config_file.read_text()), self.config)
                self.assertEqual(selected.read_bytes(), before)
                self.assertFalse((self.cache / "Ensembl_CDS").exists())

    def test_cds_missing_custom_fails_without_download_or_configuration_change(self):
        self.config["cds_fasta_file"] = str(self.folder / "missing.fa.gz")
        result, events = self.run_shell('CDS_FASTA_install "$TEST_CONFIG"')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, [])
        self.assertEqual(json.loads(self.config_file.read_text()), self.config)

    def test_cds_configured_symlink_to_managed_input_keeps_exact_spelling(self):
        target = self.cache / "Ensembl_CDS/release-116/Homo_sapiens.GRCh37.cds.all.fa.gz"
        target.parent.mkdir(parents=True)
        self.write_fasta(target)
        selected = self.folder / "selected alias.fa.gz"
        selected.symlink_to(target)
        self.config["cds_fasta_file"] = str(selected)
        before = target.read_bytes()
        result, events = self.run_shell('CDS_FASTA_install "$TEST_CONFIG"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, [])
        self.assertEqual(json.loads(self.config_file.read_text()), self.config)
        self.assertEqual(target.read_bytes(), before)

    def test_cds_invalid_custom_retains_bytes_and_path(self):
        path = self.folder / "selected.fa.gz"
        for kind in ("truncated", "text", "empty_sequence", "wrong_assembly"):
            with self.subTest(kind=kind):
                if kind == "truncated":
                    path.write_bytes(self.fasta.read_bytes()[:-5])
                elif kind == "wrong_assembly":
                    self.write_fasta(path, "GRCh38")
                else:
                    with gzip.open(path, "wt") as handle:
                        handle.write("not FASTA\n" if kind == "text" else ">ENST00000000233.1\n")
                self.config["cds_fasta_file"] = str(path)
                before = path.read_bytes()
                self.events.write_text("")
                result, events = self.run_shell('CDS_FASTA_install "$TEST_CONFIG"')
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(events, [])
                self.assertEqual(path.read_bytes(), before)
                self.assertEqual(json.loads(self.config_file.read_text()), self.config)

    def test_cds_default_acquisition_and_reuse_with_aliases(self):
        for assembly in ("hg19", "GRCh37", "hg38", "GRCh38"):
            with self.subTest(assembly=assembly):
                self.config["assembly"] = assembly
                native = "GRCh37" if assembly in ("hg19", "GRCh37") else "GRCh38"
                self.write_fasta(self.fasta, native)
                self.config["cds_fasta_file"] = ""
                self.events.write_text("")
                # The optional cache warm-up alone is a stub; validation is real.
                setup = r'''
python() {
    if [[ "$1" == */alternative_start_codon.py ]]; then event warmup "$@"; return 0; fi
    "$TEST_PYTHON" "$@"
}
'''
                result, events = self.run_shell('CDS_FASTA_install "$TEST_CONFIG"', setup)
                self.assertEqual(result.returncode, 0, result.stderr)
                target = self.cache / f"Ensembl_CDS/release-116/Homo_sapiens.{native}.cds.all.fa.gz"
                self.assertEqual(target.read_bytes(), self.fasta.read_bytes())
                self.assertEqual(json.loads(self.config_file.read_text())["cds_fasta_file"], str(target))
                self.assertLessEqual(sum(event[0] == "download" for event in events), 1)
                for event in events:
                    if event[0] == "download":
                        self.assertIn("/release-116/", " ".join(event))
                        self.assertNotIn("/current", " ".join(event))

    def test_cds_missing_release_fails_before_acquisition(self):
        self.config.pop("vep_installed_version")
        result, events = self.run_shell('CDS_FASTA_install "$TEST_CONFIG"')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("numeric vep_installed_version", result.stderr)
        self.assertEqual(events, [])

    def test_cds_download_failure_retains_staging_without_publishing_or_registering(self):
        result, events = self.run_shell('CDS_FASTA_install "$TEST_CONFIG"', r'''
wget() { event download "$@"; printf 'partial bytes\n' > "${@: -1}"; return 72; }
''')
        self.assertNotEqual(result.returncode, 0)
        target = self.cache / "Ensembl_CDS/release-116/Homo_sapiens.GRCh37.cds.all.fa.gz"
        self.assertFalse(target.exists())
        self.assertEqual(Path(str(target) + ".partial.gz").read_bytes(), b"partial bytes\n")
        self.assertEqual([event[0] for event in events], ["download"])

    def test_cds_bad_download_retained_and_update_failure_propagates(self):
        result, _ = self.run_shell('CDS_FASTA_install "$TEST_CONFIG"', r'''
wget() { printf 'corrupt text\n' > "${@: -1}"; }
''')
        self.assertNotEqual(result.returncode, 0)
        target = self.cache / "Ensembl_CDS/release-116/Homo_sapiens.GRCh37.cds.all.fa.gz"
        self.assertFalse(target.exists())
        self.assertTrue(Path(str(target) + ".partial.gz").exists())
        self.events.write_text("")
        result, _ = self.run_shell('CDS_FASTA_install "$TEST_CONFIG"', fail="cds_fasta_file")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(json.loads(self.config_file.read_text()), self.config)

    def test_alpha_invalid_supplied_source_stops_without_overwrite_or_registration(self):
        supplied = self.folder / "invalid.tsv.gz"
        supplied.write_bytes(b"retain invalid explicitly supplied bytes")
        self.config["alphamissense_prescore"] = str(supplied)
        result, events = self.run_shell('AlphaMissense_install "$TEST_CONFIG" "$TEST_CACHE" hg19 /plugins')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(supplied.read_bytes(), b"retain invalid explicitly supplied bytes")
        self.assertIn("publisher size mismatch", result.stderr)
        self.assertFalse(any(event[0] == "download" for event in events))
        self.assertFalse(any(event[0] == "update" for event in events))

    def test_alpha_hg38_helper_failure_propagates_without_other_assembly_work(self):
        directory = self.cache / "AlphaMissense"
        directory.mkdir()
        existing = directory / "AlphaMissense_hg19.tsv.gz"
        existing.write_text("fixture\n")
        index = Path(str(existing) + ".tbi")
        index.write_text("fixture index\n")
        os.utime(index, (time.time() + 5, time.time() + 5))
        result, events = self.run_shell('AlphaMissense_install "$TEST_CONFIG" "$TEST_CACHE" hg38 /plugins', r'''
python() { event alpha-helper "$@"; return 73; }
''')
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(sum(event[0] == "alpha-helper" for event in events), 1)
        self.assertIn("--assembly", events[0])
        self.assertIn("hg38", events[0])
        self.assertEqual(existing.read_text(), "fixture\n")
        self.assertFalse(any(event[0] == "update" for event in events))

    def test_pext_configuration_failures_propagate_before_success(self):
        self.config.update(pext_enabled=True, pext_tissues="")
        for key in ("pext_mean_bw", "pext_tissue_bw_dir"):
            with self.subTest(key=key):
                self.events.write_text("")
                result, events = self.run_shell('pext_install "$TEST_CONFIG"', fail=key)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(events[-1][0:2], ["update", key])

    def mechanism_inputs(self):
        raw = self.cache / "raw"
        raw.mkdir(exist_ok=True)
        schema = self.folder / "schema.json"
        schema.write_text("{}\n")
        hgnc = self.folder / "hgnc.tsv"
        hgnc.write_text("fixture\n")
        gof = self.folder / "gof.json.gz"
        gof.write_text("fixture\n")
        self.config.update(gene_mechanism_cache_dir=str(self.cache), gene_mechanism_raw_dir=str(raw),
                           gene_nonlof_mechanism_schema=str(schema), hgnc_table=str(hgnc),
                           gene_nonlof_mechanism_json=str(self.cache / "new_nonlof.json.gz"),
                           hpo_condition_mechanism_json=str(self.cache / "new_hpo_mechanism.json.gz"),
                           gofcards_exact_gof_cache=str(gof))
        return raw

    @staticmethod
    def builder_stub():
        return r'''
python() {
    event builder "$@"
    if [[ "$*" == *--fetch-sources-only* && "$*" != *--clinvar-vcv* ]]; then
        mkdir -p "$TEST_CACHE/raw/clinvar_vcv"
        printf 'managed fixture\n' > "$TEST_CACHE/raw/clinvar_vcv/ClinVarVCVRelease_00-latest_weekly.xml.gz"
    fi
}
validate_gene_pathogenic_mechanism_evidence() { return 0; }
validate_gene_nonlof_mechanism_cache() { return 0; }
validate_hpo_condition_mechanism_cache() { return 0; }
'''

    def test_custom_vcv_selection_is_passed_to_both_builder_phases_and_not_replaced(self):
        raw = self.mechanism_inputs()
        selected = self.folder / "pinned source's.xml.gz"
        selected.write_bytes(b"XML structure is separately validated by operational helper tests\n")
        self.config["clinvar_vcv_xml"] = str(selected)
        result, events = self.run_shell(r'''
gene_nonlof_mechanism_cache_install "$TEST_CONFIG" sources-only || exit $?
gene_nonlof_mechanism_cache_install "$TEST_CONFIG"
''', self.builder_stub())
        self.assertEqual(result.returncode, 0, result.stderr)
        commands = [event for event in events if event[0] == "builder"]
        self.assertEqual(len(commands), 2)
        for command in commands:
            self.assertIn("--clinvar-vcv", command)
            self.assertEqual(command[command.index("--clinvar-vcv") + 1], str(selected))
        self.assertIn("--fetch-sources-only", commands[0])
        self.assertNotIn("--fetch-sources-only", commands[1])
        self.assertEqual(json.loads(self.config_file.read_text())["clinvar_vcv_xml"], str(selected))
        self.assertFalse((raw / "clinvar_vcv").exists())

    def test_missing_custom_vcv_fails_before_builder_acquisition(self):
        self.mechanism_inputs()
        self.config["clinvar_vcv_xml"] = str(self.folder / "missing.xml.gz")
        result, events = self.run_shell('gene_nonlof_mechanism_cache_install "$TEST_CONFIG" sources-only', self.builder_stub())
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, [])

    def test_managed_vcv_default_still_uses_original_weekly_builder_and_registers_path(self):
        raw = self.mechanism_inputs()
        result, events = self.run_shell('gene_nonlof_mechanism_cache_install "$TEST_CONFIG" sources-only', self.builder_stub())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("--clinvar-vcv", events[0])
        self.assertEqual(json.loads(self.config_file.read_text())["clinvar_vcv_xml"],
                         str(raw / "clinvar_vcv/ClinVarVCVRelease_00-latest_weekly.xml.gz"))

    def test_explicit_vcv_symlink_to_managed_source_is_preserved_and_digest_bound(self):
        raw = self.mechanism_inputs()
        canonical = raw / "clinvar_vcv/ClinVarVCVRelease_00-latest_weekly.xml.gz"
        canonical.parent.mkdir()
        canonical.write_text("explicit alias identity fixture\n")
        selected = self.folder / "selected_alias.xml.gz"
        selected.symlink_to(canonical)
        self.config["clinvar_vcv_xml"] = str(selected)
        before = canonical.read_bytes()
        result, events = self.run_shell('gene_nonlof_mechanism_cache_install "$TEST_CONFIG" sources-only', self.builder_stub())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--clinvar-vcv", events[0])
        self.assertEqual(events[0][events[0].index("--clinvar-vcv") + 1], str(selected))
        self.assertEqual(canonical.read_bytes(), before)
        self.assertEqual(json.loads(self.config_file.read_text())["clinvar_vcv_xml"], str(selected))
        workdir = self.folder / "work"
        workdir.mkdir()
        (workdir / "gofcards_exact_gof.normalized.json.gz").write_text("normalized fixture\n")
        self.config["gofcards_workdir"] = str(workdir)
        self.events.write_text("")
        result, events = self.run_shell('gofcards_clinvar_injection_install "$TEST_CONFIG"', self.injection_stub())
        self.assertEqual(result.returncode, 0, result.stderr)
        command = next(event for event in events if event[0] == "injection")
        self.assertEqual(command[command.index("--clinvar-vcv-sha256") + 1], hashlib.sha256(before).hexdigest())
        self.assertEqual(command[command.index("--clinvar-vcv") + 1], str(selected))
        self.assertEqual(canonical.read_bytes(), before)
        # Selecting the literal managed pathname retains its weekly acquisition contract.
        self.config["clinvar_vcv_xml"] = str(canonical)
        self.events.write_text("")
        result, events = self.run_shell('gene_nonlof_mechanism_cache_install "$TEST_CONFIG" sources-only', self.builder_stub())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("--clinvar-vcv", events[0])

    def test_sources_only_registration_failure_propagates(self):
        self.mechanism_inputs()
        result, _ = self.run_shell('gene_nonlof_mechanism_cache_install "$TEST_CONFIG" sources-only',
                                   self.builder_stub(), fail="clinvar_vcv_xml")
        self.assertNotEqual(result.returncode, 0)

    def test_custom_injection_cannot_reuse_cache_citing_different_xml_and_updates_are_guarded(self):
        self.mechanism_inputs()
        workdir = self.folder / "work"
        workdir.mkdir()
        (workdir / "gofcards_exact_gof.normalized.json.gz").write_text("normalized fixture\n")
        selected = self.folder / "pinned.xml.gz"
        selected.write_text("XML input is separately validated by builder\n")
        self.config.update(gofcards_workdir=str(workdir), clinvar_vcv_xml=str(selected))
        setup = self.injection_stub()
        result, events = self.run_shell('gofcards_clinvar_injection_install "$TEST_CONFIG"', setup)
        self.assertEqual(result.returncode, 0, result.stderr)
        injection = next(event for event in events if event[0] == "injection")
        self.assertEqual(injection[injection.index("--clinvar-vcv") + 1], str(selected))
        self.assertEqual(injection[injection.index("--clinvar-vcv-sha256") + 1],
                         hashlib.sha256(selected.read_bytes()).hexdigest())
        for key in ("gofcards_exact_gof_cache", "gofcards_clinvar_injector_script"):
            with self.subTest(key=key):
                self.events.write_text("")
                result, events = self.run_shell('gofcards_clinvar_injection_install "$TEST_CONFIG"', setup, fail=key)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(events[-1][0:2], ["update", key])

    @staticmethod
    def injection_stub():
        return r'''
validate_gofcards_exact_gof_cache() {
    "$TEST_PYTHON" - "$@" <<'CACHE_IDENTITY_STUB'
import json, sys
try:
    cache = json.load(open(sys.argv[1]))
except (ValueError, FileNotFoundError):
    cache = {}
for argument in sys.argv[2:]:
    key, value = argument.split("=", 1)
    if key not in ("clinvar.source", "clinvar.source_sha256"):
        continue
    if cache.get("metadata", {}).get("clinvar", {}).get(key.split(".", 1)[1]) != value:
        raise SystemExit(73)
CACHE_IDENTITY_STUB
}
python() {
    event injection "$@"
    "$TEST_PYTHON" - "$@" <<'INJECTOR_IDENTITY_STUB'
import json, sys
arguments = sys.argv[1:]
path = arguments[arguments.index("--clinvar-vcv") + 1]
sha = arguments[arguments.index("--clinvar-vcv-sha256") + 1] if "--clinvar-vcv-sha256" in arguments else None
out = arguments[arguments.index("--out-json") + 1]
json.dump({"metadata": {"clinvar": {"source": path, "source_sha256": sha}}}, open(out, "w"))
INJECTOR_IDENTITY_STUB
}
'''

    def test_custom_vcv_changed_bytes_at_same_path_reinjects_and_unchanged_reuses(self):
        self.mechanism_inputs()
        workdir = self.folder / "work"
        workdir.mkdir()
        (workdir / "gofcards_exact_gof.normalized.json.gz").write_text("normalized fixture\n")
        selected = self.folder / "pinned.xml.gz"
        selected.write_text("first XML identity fixture\n")
        self.config.update(gofcards_workdir=str(workdir), clinvar_vcv_xml=str(selected))
        setup = self.injection_stub()
        for content, expect_injection in ((None, True), (None, False), ("changed at same path\n", True)):
            with self.subTest(content=content, expect_injection=expect_injection):
                if content is not None:
                    selected.write_text(content)
                self.events.write_text("")
                result, events = self.run_shell('gofcards_clinvar_injection_install "$TEST_CONFIG"', setup)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(any(event[0] == "injection" for event in events), expect_injection)
                metadata = json.loads(Path(self.config["gofcards_exact_gof_cache"]).read_text())["metadata"]["clinvar"]
                self.assertEqual(metadata["source_sha256"], hashlib.sha256(selected.read_bytes()).hexdigest())

    def test_cds_path_with_spaces_and_quote_remains_one_runtime_argument(self):
        text = PRIORITIZATION.read_text()
        function = re.search(r"^function assign_acmg_criteria\s*\([^\n]*\)\s*\{.*?^\}", text, re.M | re.S).group()
        setup = function + r'''
interpret_splicing_annotations() { return 0; }
interpret_utr_annotations() { return 0; }
check_path() { return 0; }
check_vcf_validity() { return 1; }
display_table() { return 0; }
python() { event acmg "$@"; }
'''
        for compressed in (True, False):
            with self.subTest(compressed=compressed):
                selected = self.folder / ("selected input's.fa.gz" if compressed else "selected input's.fa")
                self.write_fasta(selected, compressed=compressed)
                self.config["cds_fasta_file"] = str(selected)
                self.events.write_text("")
                result, events = self.run_shell('assign_acmg_criteria -i absent.tsv -c "$TEST_CONFIG" -t 1', setup)
                self.assertEqual(result.returncode, 0, result.stderr)
                command = next(event for event in events if event[0] == "acmg")
                index = command.index("--cds_fasta_path")
                self.assertEqual(command[index:index + 2], ["--cds_fasta_path", str(selected)])

    def test_mechanism_configuration_update_failures_propagate(self):
        self.mechanism_inputs()
        assertions = self.folder / "assertions.tsv.gz"
        evidence = self.folder / "evidence.tsv"
        dosage = self.folder / "dosage.tsv"
        mondo = self.folder / "mondo.obo"
        registry = self.folder / "registry.tsv.gz"
        for path in (assertions, evidence, dosage, registry):
            path.write_text("fixture\n")
        mondo.write_text("data-version: mondo/fixture\n")
        self.config.update(hpo_assertions=str(assertions), ddg2p_mechanism_evidence=str(evidence),
                           gene_dosage_sensitivity=str(dosage), mondo_obo=str(mondo), mondo_disease_scope_registry=str(registry))
        hpo_raw = self.folder / "hpo_raw"
        hpo_raw.mkdir()
        self.config["hpo_raw_dir"] = str(hpo_raw)
        for name in ("genes_to_phenotype.txt", "phenotype.hpoa", "genes_to_phenotype.assertions.base.tsv.gz"):
            path = hpo_raw / name
            path.write_text("fixture\n")
            os.utime(path, (time.time() + 1000, time.time() + 1000))
        for path in (registry, assertions):
            os.utime(path, (time.time() + 2000, time.time() + 2000))
        stages = {"gene_pathogenic_mechanism_cache_install": ("ddg2p_mechanism_evidence",),
                  "gene_nonlof_mechanism_cache_install": ("gene_nonlof_mechanism_builder_script", "gene_nonlof_mechanism_schema", "gene_nonlof_mechanism_json"),
                  "hpo_condition_mechanism_cache_install": ("hpo_condition_mechanism_builder_script", "hpo_condition_mechanism_json"),
                  "mondo_hpo_scope_install": ("hpo_assertions", "mondo_obo", "mondo_disease_scope_registry")}
        for stage, keys in stages.items():
            for key in keys:
                with self.subTest(stage=stage, key=key):
                    self.events.write_text("")
                    result, events = self.run_shell(f'{stage} "$TEST_CONFIG"',
                                                    self.builder_stub() + '\nvalidate_mondo_hpo_scope_cache() { return 0; }\n', fail=key)
                    self.assertNotEqual(result.returncode, 0, result.stderr)
                    self.assertEqual(events[-1][0:2], ["update", key])

    def test_mavedb_hg19_update_failure_is_not_reported_as_success(self):
        result, _ = self.run_shell('MaveDB_install "$TEST_CONFIG" "$TEST_CACHE"', fail="mavedb_file")
        self.assertNotEqual(result.returncode, 0)

    def test_grch38_and_hg38_pass_identical_mavedb_metadata_to_acmg(self):
        source_root = Path(os.environ.get("PRIVA_TEST_BEFORE_SOURCE_ROOT", ROOT))
        text = (source_root / "scripts/prioritization_vcf_per_fam.sh").read_text()
        function = re.search(r"^function assign_acmg_criteria\s*\([^\n]*\)\s*\{.*?^\}", text, re.M | re.S).group()
        self.config.update({key: str(self.folder / "resource") for key in (
            "clinvar_gene_stat", "repeat_region_file_name", "clinvar_aa_stat", "clinvar_splice_stat",
            "all_intolerant_domains", "alphamissense_intolerant_domains", "alphamissense_tranx_domain_map",
            "pm1_regions_pkl", "clinvar_patho_af_stat", "clinvar_patho_exon_af_stat", "interpro_mapping_pickle",
            "clingen_map", "gene_dosage_sensitivity", "hpo_condition_mechanism_json", "gene_nonlof_mechanism_json")})
        commands = []
        setup = function + r'''
interpret_splicing_annotations() { return 0; }
interpret_utr_annotations() { return 0; }
check_path() { return 0; }
check_vcf_validity() { return 1; }
display_table() { return 0; }
python() { event acmg "$@"; }
'''
        for assembly in ("hg38", "GRCh38"):
            self.events.write_text("")
            self.config["assembly"] = assembly
            result, events = self.run_shell('assign_acmg_criteria -i absent.tsv -c "$TEST_CONFIG" -t 1', setup)
            self.assertEqual(result.returncode, 0, result.stderr)
            command = next(event for event in events if event[0] == "acmg")
            self.assertIn("--mavedb_metadata_tsv", command)
            commands.append(command)
        self.assertEqual(commands[0], commands[1])


if __name__ == "__main__":
    unittest.main()
