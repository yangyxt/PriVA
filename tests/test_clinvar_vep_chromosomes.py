"""ClinVar VEP per chromosome: the install steps D1_split and D1_chrN, and the
gather in _ClinVar_VCF_deploy, with real bgzip, tabix and bcftools.

VEP is replaced by a stand-in that adds a CSQ field and VEP's own header
lines, so the test covers splitting, the per-chromosome calls and gathering,
not the annotation itself.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ("bgzip", "tabix", "bcftools", "jq")

PRELUDE = r'''
set -uo pipefail
event() { printf '%s\n' "$*" >> "$TEST_EVENTS"; }
log() { printf '%s\n' "$*"; }
read_yaml() { awk -v k="$2" 'index($0, k ": ") == 1 { value = substr($0, length(k) + 3); found = 1 } END { print (found ? value : "null") }' "$1"; }
update_or_append_yaml() { event "set $2"; printf '%s: %s\n' "$2" "$3" >> "$1"; }
check_vcf_validity() { [[ -s "$1" && -s "$1.tbi" ]]; }
liftover_from_ucsc_to_GRCh() { event rename; cp "$TEST_CONVERTED" "$3"; cp "$TEST_CONVERTED.tbi" "$3.tbi"; }
bash() {
    case "$2" in
        anno_VEP_cache_worker)
            local input= output= forks=
            while [[ $# -gt 0 ]]; do
                case "$1" in --input-vcf) input=$2 ;; --output-vcf) output=$2 ;; --forks) forks=$2 ;; esac
                shift
            done
            event "vep $(basename "$input") forks $forks"
            awk -v name="$(basename "$input")" '
                /^##/ { print; next }
                /^#CHROM/ {
                    print "##VEP=\"v116.3\" time=\"" name "\" assembly=\"GRCh38\""
                    print "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Gene\">"
                    print "##VEP-command-line=vep -i " name
                    print; next
                }
                $5 == "." { next }
                { OFS = "\t"; $8 = $8 ";CSQ=" $5 "|G1"; print }' "$input" > "$output" ;;
        anno_VEP_cache_candidate_phase) event "candidate ${*: -1}" ;;
        *) return 8 ;;
    esac
}
'''

HEADER = """##fileformat=VCFv4.2
##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="ClinVar allele">
##contig=<ID=chr1>
##contig=<ID=chr2>
##contig=<ID=chrM>
##contig=<ID=NT_1.1>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""
RECORDS = [("chr1", 100, "A", "G"), ("chr1", 200, "C", "T"), ("chr2", 50, "G", "A"),
           ("chrM", 10, "T", "C"), ("chrM", 20, "A", "."), ("NT_1.1", 5, "C", "G")]


@unittest.skipUnless(all(shutil.which(tool) for tool in TOOLS), "needs " + ", ".join(TOOLS))
class ClinVarVepChromosomeTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name).resolve()
        self.cache = self.directory / "cache"
        self.cache.mkdir()
        self.converted = self.write_vcf(self.directory / "converted.vcf", RECORDS)
        raw = self.write_vcf(self.directory / "raw.vcf", RECORDS)
        release = self.directory / "inputs.ready.json"
        release.write_text(json.dumps({"sources": {"vcf": {"path": str(raw)}}}))
        (self.directory / "data/liftover").mkdir(parents=True)
        (self.directory / "data/liftover/GRC_to_ucsc.contig.map.txt").write_text("1\tchr1\n")
        self.events = self.directory / "events.txt"
        self.events.touch()
        self.config = self.directory / "config.yaml"
        self.config.write_text(f"assembly: hg38\nclinvar_vcf_dir: {self.cache}\n"
                               f"clinvar_release_inputs: {release}\nthreads: 2\n")
        source = (ROOT / "scripts/install_utils.sh").read_text()
        self.functions = source[source.index("# ClinVar VEP per chromosome (install workflow)."):
                                source.index("function ClinVar_Gene_stat () {")]
        self.env = dict(os.environ, TEST_EVENTS=str(self.events), TEST_CONVERTED=str(self.converted),
                        DATA_DIR=str(self.directory / "data"), SCRIPT_DIR=str(ROOT / "scripts"))

    def write_vcf(self, plain: Path, records) -> Path:
        body = "".join(f"{c}\t{p}\t.\t{r}\t{a}\t.\t.\tALLELEID={i}\n" for i, (c, p, r, a) in enumerate(records, 1))
        plain.write_text(HEADER + body)
        subprocess.run(["bgzip", "-f", str(plain)], check=True)
        path = Path(str(plain) + ".gz")
        subprocess.run(["tabix", "-f", "-p", "vcf", str(path)], check=True)
        return path

    def run_shell(self, command):
        result = subprocess.run(["bash", "-c", PRELUDE + self.functions + "\n" + command],
                                env=self.env, text=True, capture_output=True)
        events = self.events.read_text().splitlines()
        self.events.write_text("")
        return result, events

    def config_value(self, key):
        values = [line.split(": ", 1)[1] for line in self.config.read_text().splitlines()
                  if line.startswith(key + ": ")]
        return values[-1]

    def split_and_annotate(self, pieces=("chr1", "chr2", "chr3", "chrOther")):
        result, events = self.run_shell(f'ClinVar_vep_split "{self.config}"')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(events, ["rename", "set clinvar_vep_split"])
        for name in pieces:
            result, _ = self.run_shell(f'ClinVar_vep_chrom "{self.config}" {name}')
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return Path(self.config_value("clinvar_vep_split"))

    def test_split_writes_every_piece_with_the_whole_header(self):
        manifest = json.loads(self.split_and_annotate(pieces=()).read_text())
        self.assertEqual(len(manifest["pieces"]), 25)
        counts = {piece["name"]: piece["records"] for piece in manifest["pieces"]}
        self.assertEqual((counts["chr1"], counts["chr2"], counts["chrOther"], counts["chrX"]), (2, 1, 3, 0))
        self.assertEqual(manifest["records"], len(RECORDS))
        with gzip.open(self.converted, "rt") as handle:
            body = "".join(line for line in handle if not line.startswith("#"))
        self.assertEqual(manifest["body_sha256"], hashlib.sha256(body.encode()).hexdigest())
        pieces = {piece["name"]: Path(piece["path"]).read_text() for piece in manifest["pieces"]}
        for text in pieces.values():
            self.assertTrue(text.startswith(HEADER))
        self.assertEqual([line.split("\t")[0] for line in pieces["chrOther"][len(HEADER):].splitlines()],
                         ["chrM", "chrM", "NT_1.1"])

    def test_each_piece_runs_vep_with_22_forks_and_empty_pieces_skip_vep(self):
        self.run_shell(f'ClinVar_vep_split "{self.config}"')
        result, events = self.run_shell(f'ClinVar_vep_chrom "{self.config}" chr1')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(events, ["vep chr1.vcf forks 22", "set clinvar_vep_chr1"])
        result, events = self.run_shell(f'ClinVar_vep_chrom "{self.config}" chrX')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(events, ["set clinvar_vep_chrX"])
        self.assertTrue(self.config_value("clinvar_vep_chrX").endswith("vep_chrX/output.vcf.gz"))

    def gather(self, target):
        manifest = self.config_value("clinvar_vep_split")
        return self.run_shell(f'_ClinVar_vep_gather "{self.config}" "{manifest}" "{self.converted}" "{target}"')

    def test_gather_sorts_the_pieces_between_the_candidate_records(self):
        self.split_and_annotate()
        target = self.cache / "clinvar.hg38.vep.vcf.gz"
        shutil.copy(self.converted, target)
        result, events = self.gather(target)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(events, ["candidate begin", "candidate complete"])
        with gzip.open(target, "rt") as handle:
            lines = handle.read().splitlines()
        body = [line.split("\t") for line in lines if not line.startswith("#")]
        # chrM:20 has ALT '.', which VEP leaves out.
        self.assertEqual([(row[0], row[1]) for row in body],
                         [("chr1", "100"), ("chr1", "200"), ("chr2", "50"), ("chrM", "10"), ("NT_1.1", "5")])
        self.assertTrue(all(";CSQ=" in row[7] for row in body))
        self.assertEqual(sum(line.startswith("##VEP-command-line=") for line in lines), 1)
        self.assertTrue(Path(str(target) + ".tbi").exists())

    def test_gather_refuses_records_other_than_the_split_ones(self):
        self.split_and_annotate()
        self.converted = self.write_vcf(self.directory / "other.vcf", RECORDS[:-1])
        target = self.cache / "clinvar.hg38.vep.vcf.gz"
        shutil.copy(self.converted, target)
        before = target.read_bytes()
        result, events = self.gather(target)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("made from other ClinVar records", result.stdout)
        self.assertEqual(events, [])
        self.assertEqual(target.read_bytes(), before)

    def test_gather_refuses_a_piece_with_another_vep_header(self):
        self.split_and_annotate()
        output = Path(self.config_value("clinvar_vep_chr2"))
        with gzip.open(output, "rt") as handle:
            text = handle.read().replace('Format: Allele|Gene', 'Format: Allele|Gene|Extra')
        plain = output.with_suffix("")
        plain.write_text(text)
        subprocess.run(["bgzip", "-f", str(plain)], check=True)
        target = self.cache / "clinvar.hg38.vep.vcf.gz"
        shutil.copy(self.converted, target)
        result, events = self.gather(target)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("VEP header of chr2 differs", result.stdout)
        self.assertEqual(events, [])

    def test_gather_refuses_a_missing_piece_output(self):
        self.split_and_annotate(pieces=("chr1", "chr2"))
        target = self.cache / "clinvar.hg38.vep.vcf.gz"
        shutil.copy(self.converted, target)
        result, events = self.gather(target)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("clinvar_vep_chrOther", result.stdout)
        self.assertEqual(events, [])


if __name__ == "__main__":
    unittest.main()
