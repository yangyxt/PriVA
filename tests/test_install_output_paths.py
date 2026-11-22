"""Installer output routing with genuine gzip bytes; scientific parser calls are stubbed."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest


ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "scripts/install_utils.sh"
CLINGEN_EREPO = "ClinVar Variation Id\tApplied Evidence Codes (Met)\n123\tPS3\n"
CLINGEN_DOSAGE = (
    "# comment\n#Gene Symbol\tGenomic Location\tHaploinsufficiency Score\t"
    "Haploinsufficiency Description\tHaploinsufficiency Disease ID\tTriplosensitivity Score\n"
    "GENE\tchr1:1-2\t3\tSufficient evidence\tMONDO:0000001\t0\n"
)
PRELUDE = r'''
source "$1"
DATA_DIR="$TEST_SOURCE_DATA"
log() { printf '%s\n' "$*" >&2; }
event() {
    "$TEST_PYTHON" -c 'import json,os,sys; open(os.environ["TEST_EVENTS"],"a").write(json.dumps(sys.argv[1:])+"\n")' "$@"
}
read_yaml() {
    "$TEST_PYTHON" -c 'import json,sys; v=json.load(open(sys.argv[1])).get(sys.argv[2]); print("null" if v is None else v)' "$1" "$2"
}
update_yaml() {
    event update "$2" "$3"
    [[ "$TEST_FAIL" != "update_$2" ]] || return 71
    "$TEST_PYTHON" -c 'import json,sys; p,k,v=sys.argv[1:]; d=json.load(open(p)); d[k]=v; open(p,"w").write(json.dumps(d))' "$@"
}
wget() {
    event download "$@"
    local destination=${@: -1} stage
    case "$*" in
        *interpro.xml.gz*) stage=xml ;;
        *ClinGen_gene_curation_list*) stage=dosage ;;
        *) stage=erepo ;;
    esac
    if [[ "$TEST_FAIL" == "download_$stage" ]]; then
        printf 'retained partial bytes\n' > "$destination"
        return 72
    fi
    case "$stage" in
        xml) cp "$TEST_XML" "$destination" ;;
        dosage) cp "$TEST_CLINGEN_DOSAGE" "$destination" ;;
        erepo) cp "$TEST_CLINGEN_EREPO" "$destination" ;;
    esac
}
python() {
    # Keep scientific parsers stubbed; execute source validation/receipts for real.
    if [[ "$1" == *mechanism_source_input.py ]]; then "$TEST_PYTHON" "$@"; return; fi
    event parser "$@"
    [[ "$TEST_FAIL" != parser ]] || return 73
    "$TEST_PYTHON" - "$@" <<'PARSER_STUB'
import gzip, pathlib, pickle, sys
arguments = sys.argv[1:]
if pathlib.Path(arguments[0]).name == "protein_domain_mapping.py":
    output = arguments[arguments.index("--mapping_output") + 1]
    pathlib.Path(arguments[arguments.index("--report") + 1]).write_text("manual curation report fixture\n")
else:
    output = arguments[arguments.index("-o") + 1]
with gzip.open(output, "wb") as stream:
    pickle.dump({"routing_fixture": True}, stream)
PARSER_STUB
}
'''


class InstallerOutputPaths(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="priva_output_paths_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.source_data = self.folder / "source_data"
        for name in ("InterPro", "clingen"):
            directory = self.source_data / name
            directory.mkdir(parents=True)
            (directory / "curated_source.tsv").write_text("curated source remains unchanged\n")
        self.cache = self.folder / "data side's cache"
        self.mapping = self.cache / "maps" / "selected_mapping.pkl.gz"
        self.dosage = self.cache / "dosage" / "selected_dosage.tsv"
        self.report = self.source_data / "InterPro/Interpro_format_report.txt"
        self.xml_fixture = self.folder / "xml_fixture.gz"
        with gzip.open(self.xml_fixture, "wt") as stream:
            stream.write("<interprodb/>\n")
        self.config = {"assembly": "hg19", "interpro_mapping_pickle": "Interpro_entry_mapping.pkl.gz",
                       "interpro_xml_gz_url": "https://fixtures.invalid/interpro.xml.gz",
                       "alphamissense_vep_vcf": str(self.folder / "alpha.vcf.gz"),
                       "clinvar_vcf": str(self.folder / "clinvar.vcf.gz"),
                       "clingen_map": "clingen_map.hg19.pkl.gz",
                       "gene_dosage_sensitivity": "gene_dosage_sensitivity.hg19.tsv"}
        self.config_file = self.folder / "config's path.json"
        self.events_file = self.folder / "events.jsonl"
        self.events_file.touch()
        self.clingen_erepo_fixture = self.folder / "erepo_fixture.tsv"
        self.clingen_dosage_fixture = self.folder / "dosage_fixture.tsv"
        self.clingen_erepo_fixture.write_text(CLINGEN_EREPO)
        self.clingen_dosage_fixture.write_text(CLINGEN_DOSAGE)
        self.env = dict(os.environ, TEST_SOURCE_DATA=str(self.source_data), TEST_EVENTS=str(self.events_file),
                        TEST_XML=str(self.xml_fixture), TEST_CONFIG=str(self.config_file),
                        TEST_PYTHON=sys.executable, TEST_FAIL="none", TMPDIR=str(self.folder),
                        TEST_CLINGEN_EREPO=str(self.clingen_erepo_fixture),
                        TEST_CLINGEN_DOSAGE=str(self.clingen_dosage_fixture))

    def run_function(self, name, fail="none", setup=""):
        self.config_file.write_text(json.dumps(self.config))
        result = subprocess.run(["bash", "--noprofile", "--norc", "-c",
                                 PRELUDE + setup + f'\n{name} "$TEST_CONFIG"', "test", str(INSTALLER)],
                                env=dict(self.env, TEST_FAIL=fail), text=True, capture_output=True, timeout=40)
        events = [json.loads(line) for line in self.events_file.read_text().splitlines()]
        return result, events

    def source_snapshot(self):
        return {str(path.relative_to(self.source_data)): path.read_bytes()
                for path in self.source_data.rglob("*") if path.is_file() and path != self.report}

    def explicit_interpro(self):
        self.config["interpro_mapping_pickle"] = str(self.mapping)
        self.report = self.mapping.parent / "Interpro_format_report.txt"

    def explicit_clingen(self):
        self.config["clingen_map"] = str(self.mapping)
        self.config["gene_dosage_sensitivity"] = str(self.dosage)

    def current_interpro(self, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        xml = path.parent / "interpro.xml.gz"
        xml.write_bytes(self.xml_fixture.read_bytes())
        Path(str(xml) + ".sha256").write_text(hashlib.sha256(xml.read_bytes()).hexdigest() + "\n")
        path.write_bytes(b"existing mapping retained\n")
        newest = max(time.time(), (ROOT / "scripts/protein_domain_mapping.py").stat().st_mtime) + 10
        os.utime(path, (newest, newest))
        return xml

    def test_interpro_exact_absolute_output_routes_cache_metadata_and_keeps_manual_report(self):
        self.explicit_interpro()
        before = self.source_snapshot()
        result, events = self.run_function("InterPro_parsing")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([event[0] for event in events], ["download", "parser", "update"])
        self.assertEqual(events[1], ["parser", str(ROOT / "scripts/protein_domain_mapping.py"),
                                     "--xml", str(self.mapping.parent / "interpro.xml.gz"),
                                     "--mapping_output", str(self.mapping),
                                     "--vcf", self.config["alphamissense_vep_vcf"], "--report", str(self.report)])
        self.assertTrue(self.mapping.is_file())
        xml = self.mapping.parent / "interpro.xml.gz"
        self.assertEqual(xml.read_bytes(), self.xml_fixture.read_bytes())
        self.assertEqual(Path(str(xml) + ".sha256").read_text().strip(), hashlib.sha256(xml.read_bytes()).hexdigest())
        self.assertEqual(self.report.read_text(), "manual curation report fixture\n")
        self.assertEqual(self.source_snapshot(), before)

    def test_interpro_managed_generation_does_not_require_source_report_directory(self):
        self.explicit_interpro()
        for path in (self.source_data / "InterPro").iterdir():
            path.unlink()
        (self.source_data / "InterPro").rmdir()
        result, _ = self.run_function("InterPro_parsing")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(self.report.is_file())
        self.assertFalse((self.source_data / "InterPro").exists())

    def test_interpro_absolute_current_cache_reuses_exact_input_without_source_writes(self):
        self.explicit_interpro()
        xml = self.current_interpro(self.mapping)
        before = self.source_snapshot()
        original = (self.mapping.read_bytes(), xml.read_bytes())
        result, events = self.run_function("InterPro_parsing")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, [])
        self.assertEqual((self.mapping.read_bytes(), xml.read_bytes()), original)
        self.assertEqual(self.source_snapshot(), before)

    def test_interpro_current_bare_name_is_registered_at_legacy_path(self):
        selected = self.source_data / "InterPro/Interpro_entry_mapping.pkl.gz"
        self.current_interpro(selected)
        result, events = self.run_function("InterPro_parsing")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, [["update", "interpro_mapping_pickle", str(selected)]])
        self.assertEqual(json.loads(self.config_file.read_text())["interpro_mapping_pickle"], str(selected))

    def test_interpro_bare_default_generation_keeps_original_output_and_report_paths(self):
        result, events = self.run_function("InterPro_parsing")
        self.assertEqual(result.returncode, 0, result.stderr)
        command = next(event for event in events if event[0] == "parser")
        self.assertEqual(command[command.index("--mapping_output") + 1],
                         str(self.source_data / "InterPro/Interpro_entry_mapping.pkl.gz"))
        self.assertEqual(command[command.index("--report") + 1], str(self.report))

    def test_interpro_absolute_download_failure_retains_selected_partial_and_configuration(self):
        self.explicit_interpro()
        before = self.source_snapshot()
        result, events = self.run_function("InterPro_parsing", fail="download_xml")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual([event[0] for event in events], ["download"])
        self.assertEqual((self.mapping.parent / "interpro.xml.gz.partial.gz").read_bytes(), b"retained partial bytes\n")
        self.assertFalse((self.mapping.parent / "interpro.xml.gz").exists())
        self.assertEqual(json.loads(self.config_file.read_text()), self.config)
        self.assertEqual(self.source_snapshot(), before)

    def test_interpro_absolute_parser_and_configuration_failures_propagate(self):
        self.explicit_interpro()
        for failure in ("parser", "update_interpro_mapping_pickle"):
            with self.subTest(failure=failure):
                self.events_file.write_text("")
                before = self.source_snapshot()
                result, events = self.run_function("InterPro_parsing", fail=failure)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(json.loads(self.config_file.read_text()), self.config)
                self.assertEqual(self.source_snapshot(), before)

    def test_interpro_absolute_parent_failure_stops_before_acquisition(self):
        self.explicit_interpro()
        self.cache.mkdir(parents=True)
        self.mapping.parent.write_text("blocks cache directory\n")
        result, events = self.run_function("InterPro_parsing")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, [])

    def test_interpro_absolute_identity_mismatch_is_retained_without_rebuild(self):
        self.explicit_interpro()
        xml = self.current_interpro(self.mapping)
        Path(str(xml) + ".sha256").write_text("0" * 64 + "\n")
        result, events = self.run_function("InterPro_parsing")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, [])
        self.assertEqual(self.mapping.read_bytes(), b"existing mapping retained\n")

    def test_clingen_absolute_outputs_are_independent_and_preserve_source_tree(self):
        self.explicit_clingen()
        before = self.source_snapshot()
        result, events = self.run_function("ClinGen_deploy")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([event[0] for event in events], ["download", "parser", "update", "download", "update"])
        self.assertEqual(events[1], ["parser", str(ROOT / "scripts/prepare_clingen_map.py"),
                                     "-v", self.config["clinvar_vcf"], "-t", str(self.mapping.parent / "clingen_erepo.tsv"),
                                     "-o", str(self.mapping)])
        self.assertEqual(self.dosage.read_text(), CLINGEN_DOSAGE.split("\n", 1)[1])
        downloads = [event[-1] for event in events if event[0] == "download"]
        for staged, target in zip(downloads, (self.mapping.parent / "clingen_erepo.tsv", self.dosage)):
            self.assertEqual(Path(staged).parent, target.parent)
            self.assertTrue(Path(staged).name.startswith(target.name + ".download."))
            self.assertNotEqual(Path(staged), target)
        self.assertEqual((self.mapping.parent / "clingen_erepo.tsv").read_text(), CLINGEN_EREPO)
        receipt = json.loads(Path(str(self.dosage) + ".source.json").read_text())
        self.assertEqual(receipt["sha256"], hashlib.sha256(self.dosage.read_bytes()).hexdigest())
        self.assertTrue(self.mapping.is_file())
        self.assertEqual(self.source_snapshot(), before)
        configured = json.loads(self.config_file.read_text())
        self.assertEqual(configured["clingen_map"], str(self.mapping))
        self.assertEqual(configured["gene_dosage_sensitivity"], str(self.dosage))

    def test_clingen_current_selected_sources_reuse_without_download_in_original_parser_order(self):
        self.explicit_clingen()
        self.mapping.parent.mkdir(parents=True)
        self.dosage.parent.mkdir(parents=True)
        raw = self.mapping.parent / "clingen_erepo.tsv"
        raw.write_text("retained exact eRepo source\n")
        self.dosage.write_text("#Gene\tScore\nGENE\t3\n")
        result, events = self.run_function("ClinGen_deploy")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([event[0] for event in events], ["parser", "update", "update"])
        self.assertEqual(raw.read_text(), "retained exact eRepo source\n")

    def test_clingen_only_absolute_map_keeps_bare_dosage_legacy_destination(self):
        self.config["clingen_map"] = str(self.mapping)
        result, events = self.run_function("ClinGen_deploy")
        self.assertEqual(result.returncode, 0, result.stderr)
        configured = json.loads(self.config_file.read_text())
        self.assertEqual(configured["clingen_map"], str(self.mapping))
        self.assertEqual(configured["gene_dosage_sensitivity"],
                         str(self.source_data / "clingen/gene_dosage_sensitivity.hg19.tsv"))

    def test_clingen_only_absolute_dosage_keeps_bare_map_legacy_destination(self):
        self.config["gene_dosage_sensitivity"] = str(self.dosage)
        result, events = self.run_function("ClinGen_deploy")
        self.assertEqual(result.returncode, 0, result.stderr)
        configured = json.loads(self.config_file.read_text())
        self.assertEqual(configured["clingen_map"], str(self.source_data / "clingen/clingen_map.hg19.pkl.gz"))
        self.assertEqual(configured["gene_dosage_sensitivity"], str(self.dosage))

    def test_clingen_absolute_erepo_failure_retains_partial_and_never_uses_source_tree(self):
        self.explicit_clingen()
        before = self.source_snapshot()
        result, events = self.run_function("ClinGen_deploy", fail="download_erepo")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual([event[0] for event in events], ["download"])
        canonical = self.mapping.parent / "clingen_erepo.tsv"
        stages = list(canonical.parent.glob(canonical.name + ".download.*"))
        self.assertEqual(len(stages), 1)
        self.assertEqual(stages[0].read_bytes(), b"retained partial bytes\n")
        self.assertEqual(events[0][-1], str(stages[0]))
        self.assertFalse(canonical.exists())
        self.assertFalse(self.mapping.exists())
        self.assertEqual(self.source_snapshot(), before)
        self.assertEqual(json.loads(self.config_file.read_text()), self.config)

    def test_clingen_parser_or_map_configuration_failure_stops_before_dosage_mutation(self):
        self.explicit_clingen()
        self.dosage.parent.mkdir(parents=True)
        self.dosage.write_bytes(b"existing dosage retained\n")
        for failure in ("parser", "update_clingen_map"):
            with self.subTest(failure=failure):
                self.events_file.write_text("")
                result, events = self.run_function("ClinGen_deploy", fail=failure)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(self.dosage.read_bytes(), b"existing dosage retained\n")
                self.assertFalse(any(event[:2] == ["update", "gene_dosage_sensitivity"] for event in events))

    def test_clingen_dosage_download_failure_stops_before_filter_and_registration(self):
        self.explicit_clingen()
        before = self.source_snapshot()
        result, events = self.run_function("ClinGen_deploy", fail="download_dosage")
        self.assertNotEqual(result.returncode, 0)
        stages = list(self.dosage.parent.glob(self.dosage.name + ".download.*"))
        self.assertEqual(len(stages), 1)
        self.assertEqual(stages[0].read_bytes(), b"retained partial bytes\n")
        self.assertEqual(events[-1][-1], str(stages[0]))
        self.assertFalse(self.dosage.exists())
        self.assertFalse(Path(str(self.dosage) + ".source.json").exists())
        self.assertEqual(self.source_snapshot(), before)
        self.assertEqual(events[-1][0], "download")
        self.assertFalse(any(event[:2] == ["update", "gene_dosage_sensitivity"] for event in events))

    def test_clingen_dosage_filter_failure_preserves_downloaded_bytes(self):
        self.explicit_clingen()
        before = self.source_snapshot()
        result, events = self.run_function("ClinGen_deploy", setup="awk() { return 75; }")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.dosage.read_text(), CLINGEN_DOSAGE)
        receipt = json.loads(Path(str(self.dosage) + ".source.json").read_text())
        self.assertEqual(receipt["sha256"], hashlib.sha256(self.dosage.read_bytes()).hexdigest())
        self.assertEqual(self.source_snapshot(), before)
        self.assertFalse(any(event[:2] == ["update", "gene_dosage_sensitivity"] for event in events))


if __name__ == "__main__":
    unittest.main()
