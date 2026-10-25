"""AlphaMissense stage routing and failure barriers; scientific procedures stubbed."""
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
import test_alphamissense_annotation_base as BASE_FIXTURE


ROOT = Path(__file__).resolve().parents[1]
PRELUDE = r'''
source "$1"
log() { printf '%s\n' "$*" >&2; }
read_yaml() { "$TEST_PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1])).get(sys.argv[2], ""))' "$1" "$2"; }
update_yaml() { printf 'UPDATE:%s:%s\n' "$2" "$3" >> "$TEST_EVENTS"; "$TEST_PYTHON" -c 'import json,sys;p,k,v=sys.argv[1:];d=json.load(open(p));d[k]=v;open(p,"w").write(json.dumps(d))' "$@"; }
check_vcf_validity() { [[ -s $1 && -s $1.tbi ]]; }
check_vcf_infotags() { [[ -s $1 ]]; }
fixture_vcf() {
    "$TEST_PYTHON" - "$1" "$TEST_CONFIG" "$TEST_FAILURE" "$2" <<'ALPHA_TEST_VCF_PY'
import json,subprocess,sys
path,config,failure,kind=sys.argv[1:]
settings=json.load(open(config)); major=settings["vep_installed_version"]
assembly="GRCh37.p13" if settings["assembly"] in ("hg19","GRCh37") else "GRCh38.p14"
if failure == "header_major": major="115"
if failure == "header_assembly": assembly="GRCh37.p13"
header=f'##VEP="v{major}.3" assembly="{assembly}"\n' if kind == "annotated" else ""
if failure == "header_missing": header=""
text=("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n"+header+
      '##INFO=<ID=CSQ,Number=.,Type=String,Description="fixture">\n'+
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"+
      "chr1\t1\t.\tA\tC\t.\tPASS\tCSQ=C|missense_variant\n")
with open(path,"wb") as output: subprocess.run(["bgzip","-c"],input=text.encode(),stdout=output,check=True)
subprocess.run(["tabix","-f","-p","vcf",path],check=True)
ALPHA_TEST_VCF_PY
}
python() {
    if [[ $1 == */alphmis_tsv2vcf.py ]]; then
        printf 'CONVERT:%s:%s\n' "$2" "$3" >> "$TEST_EVENTS"
        [[ $TEST_FAILURE != convert ]] || return 72
        fixture_vcf "$3" converted
    elif [[ $1 == */stat_protein_domain_amscores.py ]]; then
        printf 'STAT:%s:%s\n' "$2" "$3" >> "$TEST_EVENTS"
        [[ $TEST_FAILURE != stat ]] || return 73
        printf 'statistics fixture\n' > "$3"
    elif [[ $1 == */alphamissense_cache.py && ${2:-} != check-annotation ]]; then
        printf 'RAW:%s\n' "$*" >> "$TEST_EVENTS"
        printf '%s/AlphaMissense/AlphaMissense_hg38.tsv.gz\n' "$TEST_CACHE"
    elif [[ $1 == */am_pick_intolerant_domains.py ]]; then
        printf 'DOMAINS:%s\n' "$*" >> "$TEST_EVENTS"
        [[ $TEST_FAILURE != domains ]] || return 75
        "$TEST_PYTHON" -c '
from pathlib import Path
import sys
args=sys.argv[1:]
folder=Path(args[args.index("--output_dir")+1])
assembly=args[args.index("--assembly")+1]
for stem, suffix in (("domain_tolerance_analysis", "tsv"), ("transcript_pp_domain_mapping", "pkl"), ("all_intolerant_domains", "pkl")):
    (folder / f"{stem}.{assembly}.{suffix}").write_bytes((stem + " fixture bytes").encode())
' "$@"
    else
        "$TEST_PYTHON" "$@"
    fi
}
basic_vep_annotation() {
    printf 'VEP:%s\n' "$*" >> "$TEST_EVENTS"
    [[ $TEST_FAILURE != vep ]] || return 74
    local input=$2 output=${2/.vcf*/.vep.vcf.gz}
    fixture_vcf "$output" annotated
    printf 'stats fixture\n' > "${input/.vcf*/.vep.stats.html}"
}
'''


class AlphaConfiguredOutputs(unittest.TestCase):
    def setUp(self):
        self.base_fixture = BASE_FIXTURE.AnnotationBaseFixtures()
        self.base_fixture.setUp()
        self.addCleanup(self.base_fixture.doCleanups)
        temporary = tempfile.TemporaryDirectory(prefix="alpha_outputs_", dir=os.environ.get("PRIVA_TEST_TMPDIR"))
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.source = self.folder / "external original.tsv.gz"
        self.source.write_bytes(b"source identity retained")
        Path(str(self.source) + ".tbi").write_bytes(b"source index fixture")
        self.managed = self.folder / "managed root"
        self.events = self.folder / "events"
        self.events.touch()
        self.config = self.folder / "config.json"
        self.settings = dict(self.base_fixture.cfg, cache_root=str(self.folder),
                             alphamissense_prescore=str(self.source),
                             alphamissense_vcf=str(self.managed / "custom raw.vcf.gz"),
                             alphamissense_vep_vcf=str(self.managed / "custom annotated.vcf.gz"),
                             alphamissense_vep_base_vcf=str(self.managed / "complete base.vcf.gz"),
                             alphamissense_pd_stat=str(self.managed / "custom statistic.pkl"))
        self.config.write_text(json.dumps(self.settings))
        self.env = dict(os.environ, TEST_PYTHON=sys.executable, TEST_EVENTS=str(self.events),
                        TEST_CACHE=str(self.managed), TEST_CONFIG=str(self.config), TEST_FAILURE="none")

    def invoke(self, command, failure="none"):
        return subprocess.run(["bash", "--noprofile", "--norc", "-c", PRELUDE + command,
                               "test", str(ROOT / "scripts/install_utils.sh")],
                              env=dict(self.env, TEST_FAILURE=failure), capture_output=True, text=True)

    def test_annotation_honors_both_configured_paths_and_preserves_original_source(self):
        result = self.invoke('AlphaMissense_anno "$TEST_CONFIG"')
        self.assertEqual(result.returncode, 0, result.stderr)
        with gzip.open(self.settings["alphamissense_vcf"], "rt") as converted:
            self.assertIn("chr1\t1\t.\tA\tC", converted.read())
        with gzip.open(self.settings["alphamissense_vep_vcf"], "rt") as output:
            self.assertIn('##VEP="v116.3" assembly="GRCh38.p14"', output.read())
        self.assertEqual(self.source.read_bytes(), b"source identity retained")
        self.assertEqual(list(self.folder.glob("external*.vcf*")), [])
        lines = self.events.read_text().splitlines()
        self.assertTrue(lines[0].startswith("CONVERT:"))
        self.assertTrue(lines[1].startswith("VEP:"))
        self.assertEqual(sum(line.startswith("UPDATE:") for line in lines), 3)

    def completed_annotation(self, metadata='##VEP="v116.3" assembly="GRCh38.p14"'):
        result = self.invoke('AlphaMissense_anno "$TEST_CONFIG"')
        self.assertEqual(result.returncode, 0, result.stderr)
        converted = Path(self.settings["alphamissense_vcf"])
        annotated = Path(self.settings["alphamissense_vep_vcf"])
        if metadata != '##VEP="v116.3" assembly="GRCh38.p14"':
            with gzip.open(annotated, "wt") as output:
                output.write("##fileformat=VCFv4.2\n" + metadata + "\n" +
                         "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
                         "chr1\t1\t.\tA\tC\t.\tPASS\tCSQ=C|missense_variant\n")
        self.events.write_text("")
        return converted, annotated

    def test_matching_major_and_assembly_reuse_preserves_outputs_and_source(self):
        converted, annotated = self.completed_annotation()
        before = [(path.read_bytes(), path.stat().st_mtime_ns, path.stat().st_ctime_ns)
                  for path in (self.source, converted, annotated)]
        result = self.invoke('AlphaMissense_anno "$TEST_CONFIG"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.events.read_text(), "")
        self.assertEqual([(path.read_bytes(), path.stat().st_mtime_ns, path.stat().st_ctime_ns)
                          for path in (self.source, converted, annotated)], before)

    def test_major_change_rebuilds_in_original_order_and_selects_new_annotation(self):
        self.completed_annotation()
        self.settings["vep_installed_version"] = "117"
        self.base_fixture.vep.write_text("#!/bin/sh\nprintf 'ensembl-vep : 117.3\\n'\n")
        cache = Path(self.settings["vep_cache_dir"]) / "homo_sapiens_merged/117_GRCh38"
        cache.mkdir()
        (cache / "info.txt").write_text("new selected release fixture\n")
        self.config.write_text(json.dumps(self.settings))
        result = self.invoke('AlphaMissense_anno "$TEST_CONFIG"')
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = self.events.read_text().splitlines()
        self.assertTrue(lines[0].startswith("CONVERT:"))
        self.assertTrue(lines[1].startswith("VEP:"))
        with gzip.open(self.settings["alphamissense_vep_vcf"], "rt") as output:
            self.assertIn('##VEP="v117.3"', output.read())
        self.assertEqual(self.source.read_bytes(), b"source identity retained")

    def test_matching_assembly_alias_reuses_but_other_assembly_rebuilds(self):
        self.completed_annotation()
        self.settings["assembly"] = "GRCh38"
        self.config.write_text(json.dumps(self.settings))
        result = self.invoke('AlphaMissense_anno "$TEST_CONFIG"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.events.read_text(), "")
        self.settings["assembly"] = "hg19"
        cache = Path(self.settings["vep_cache_dir"]) / "homo_sapiens_merged/116_GRCh37"
        cache.mkdir()
        (cache / "info.txt").write_text("selected other assembly fixture\n")
        self.config.write_text(json.dumps(self.settings))
        result = self.invoke('AlphaMissense_anno "$TEST_CONFIG"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("CONVERT:", self.events.read_text())
        with gzip.open(self.settings["alphamissense_vep_vcf"], "rt") as output:
            self.assertIn('assembly="GRCh37.p13"', output.read())

    def test_invalid_runtime_header_is_restored_from_verified_base_without_vep(self):
        for metadata in ('', '##VEP="v116.3"',
                         '##VEP="v116.3" assembly="GRCh38"\n##VEP="v117" assembly="GRCh38"',
                         '##VEP="v116.3" assembly="GRCh38" assembly="GRCh37"'):
            with self.subTest(metadata=metadata):
                self.events.write_text("")
                self.completed_annotation(metadata)
                result = self.invoke('AlphaMissense_anno "$TEST_CONFIG"')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertNotIn("CONVERT:", self.events.read_text())
                self.assertNotIn("VEP:", self.events.read_text())
                with gzip.open(self.settings["alphamissense_vep_vcf"], "rt") as restored:
                    self.assertIn('##VEP="v116.3" assembly="GRCh38.p14"', restored.read())

    def test_incompatible_new_annotation_keeps_previous_outputs_and_configuration(self):
        for failure in ("header_major", "header_assembly", "header_missing"):
            with self.subTest(failure=failure):
                self.events.write_text("")
                converted, annotated = self.completed_annotation('##VEP="v115" assembly="GRCh38"')
                receipt = Path(self.settings["alphamissense_vep_base_vcf"] + ".annotation-base.json")
                old = json.loads(receipt.read_text())
                old["science"]["converter_sha256"] = "0" * 64
                receipt.write_text(json.dumps(old))
                paths = (self.source, converted, annotated, Path(str(converted) + ".tbi"), Path(str(annotated) + ".tbi"))
                before = [path.read_bytes() for path in paths]
                result = self.invoke('AlphaMissense_anno "$TEST_CONFIG"', failure)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual([path.read_bytes() for path in paths], before)
                self.assertNotIn("UPDATE:", self.events.read_text())

    def test_conversion_and_annotation_failure_never_publish_configured_outputs(self):
        for failure in ("convert", "vep"):
            with self.subTest(failure=failure):
                self.events.write_text("")
                result = self.invoke('AlphaMissense_anno "$TEST_CONFIG"', failure)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(Path(self.settings["alphamissense_vep_vcf"]).exists())
                self.assertNotIn("UPDATE:", self.events.read_text())

    def test_statistic_honors_configured_destination(self):
        result = self.invoke('AlphaMissense_stat "$TEST_CONFIG"')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(Path(self.settings["alphamissense_pd_stat"]).read_text(), "statistics fixture\n")
        self.assertFalse(Path(self.settings["alphamissense_vep_vcf"].replace(".vcf.gz", ".prot.domain.stats.pkl")).exists())

    def test_failed_statistic_does_not_replace_previous_output_or_register_it(self):
        target = Path(self.settings["alphamissense_pd_stat"])
        target.parent.mkdir()
        target.write_bytes(b"retain previous statistic")
        os.utime(target, (1, 1))
        result = self.invoke('AlphaMissense_stat "$TEST_CONFIG"', "stat")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(target.read_bytes(), b"retain previous statistic")
        self.assertNotIn("UPDATE:", self.events.read_text())

    def test_named_raw_install_dispatches_selected_assembly(self):
        result = self.invoke('AlphaMissense_install --config "$TEST_CONFIG" --cache-dir "$TEST_CACHE" --assembly hg38')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--assembly hg38", self.events.read_text())
        self.assertIn("UPDATE:alphamissense_prescore:", self.events.read_text())

    def test_domain_outputs_honor_three_configured_managed_destinations(self):
        self.domain_inputs()
        keys = ("alphamissense_intolerant_domains", "alphamissense_tranx_domain_map", "all_intolerant_domains")
        for key in keys:
            self.settings[key] = str(self.managed / key / "selected output")
        self.config.write_text(json.dumps(self.settings))
        result = self.invoke('AlphaMissense_pick_intolerant_domains "$TEST_CONFIG"')
        self.assertEqual(result.returncode, 0, result.stderr)
        for key in keys:
            self.assertTrue(Path(self.settings[key]).is_file())
            self.assertIn(f"UPDATE:{key}:{self.settings[key]}", self.events.read_text())
        self.assertEqual(sum(line.startswith("UPDATE:") for line in self.events.read_text().splitlines()), 3)

    def test_domain_failure_preserves_previous_selected_outputs(self):
        self.domain_inputs()
        keys = ("alphamissense_intolerant_domains", "alphamissense_tranx_domain_map", "all_intolerant_domains")
        for key in keys:
            path = self.managed / key / "selected output"
            path.parent.mkdir(parents=True)
            path.write_bytes(b"retain selected domain output")
            os.utime(path, (1, 1))
            self.settings[key] = str(path)
        self.config.write_text(json.dumps(self.settings))
        result = self.invoke('AlphaMissense_pick_intolerant_domains "$TEST_CONFIG"', "domains")
        self.assertNotEqual(result.returncode, 0)
        for key in keys:
            self.assertEqual(Path(self.settings[key]).read_bytes(), b"retain selected domain output")
        self.assertNotIn("UPDATE:", self.events.read_text())

    def domain_inputs(self):
        # The real completion guard runs before the labelled scientific stub.
        # Supply its required inputs so these cases reach the intended writer.
        stats = Path(self.settings["alphamissense_pd_stat"])
        stats.parent.mkdir(parents=True, exist_ok=True)
        stats.write_bytes(b"completed domain-statistics fixture")
        interpro = self.folder / "selected-interpro.pkl"
        interpro.write_bytes(b"selected InterPro fixture")
        self.settings["interpro_mapping_pickle"] = str(interpro)


if __name__ == "__main__":
    unittest.main()
