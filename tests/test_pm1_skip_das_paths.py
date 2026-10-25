"""Run real tiny PM1 builds for empty/default and explicitly selected DAS paths.

python -m unittest discover -s tests -p test_pm1_skip_das_paths.py -v
No DAS model is rebuilt, no publisher cache is downloaded, and no family is run.
"""
from contextlib import redirect_stdout
import importlib.util
import io
from pathlib import Path
import pickle
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

import yaml


SOURCE = Path(__file__).resolve().parents[1] / "scripts/prepare_pm1_regions.py"
spec = importlib.util.spec_from_file_location("pm1_skip_das_paths", SOURCE)
pm1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pm1)


class ExistingDasPathTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory(prefix="pm1-paths-")
        self.root = Path(self.temporary.name)
        self.previous_root = pm1.PRIVA_DIR
        pm1.PRIVA_DIR = str(self.root)
        self.addCleanup(setattr, pm1, "PRIVA_DIR", self.previous_root)
        self.addCleanup(self.temporary.cleanup)
        for relative in ("data/HCSeeker", "data/gnomAD_RMC", "cache"):
            (self.root / relative).mkdir(parents=True)
        (self.root / "data/HCSeeker/HC_spots.hg38.tsv").write_text(
            "type\tgene\tunused1\tunused2\tisoform\tstart\tend\n"
            "hotspot\tGENE1\t.\t.\tNM_000001.2\t3\t4\n"
            "coldspot\tGENE1\t.\t.\tNM_000001.2\t5\t8\n")
        (self.root / "data/gnomAD_RMC/gnomad_v2.1.1_rmc_flat.hg38.tsv").write_text(
            "transcript\tstart_aa_raw\tend_aa_raw\toe\tp\n"
            "ENST000001\tArg10\tGly12\t0.1\t0.001\n"
            "ENST000002\tArg2\tGly3\t0.9\t0.7\n")
        self.mapping = self.root / "mapping.tsv"
        self.mapping.write_text(
            "Gene stable ID\tTranscript stable ID\tHGNC symbol\tUniProtKB/Swiss-Prot ID\tRefSeq mRNA ID\n"
            "ENSG1\tENST000001\tGENE1\tP1\tNM_000001\n"
            "ENSG1\tENST000003\tGENE1\tP1\tNM_000001\n")
        self.default_domains = self.root / "cache/all_intolerant_domains.hg38.pkl"
        self.default_domains.write_bytes(pickle.dumps({"ENSG1:Pfam:PF00001"}))
        self.config = {"alphamissense_pd_stat": str(self.root / "cache/model.pkl")}
        self.args = SimpleNamespace(config=str(self.root / "config.yaml"),
            output=str(self.root / "pm1.pkl"), assembly="hg38", threads=1,
            nm_to_enst_tsv=str(self.mapping), skip_das=True,
            rmc_oe_threshold=0.4, rmc_p_bh_threshold=0.05)

    def build(self, expected):
        Path(self.args.config).write_text(yaml.safe_dump(self.config))
        with redirect_stdout(io.StringIO()):
            pm1.prepare_pm1_regions(self.args)
        with Path(self.args.output).open("rb") as stream:
            observed = pickle.load(stream)
        self.assertEqual(observed, {
            "das_intolerant_domains": expected,
            "hcseeker_hotspots": {"ENST000001": {3, 4}, "ENST000003": {3, 4}},
            "rmc_constrained": {"ENST000001": {10, 11, 12}},
            "metadata": {"assembly": "hg38", "rmc_oe_threshold": 0.4,
                "rmc_p_bh_threshold": 0.05, "n_das_intolerant_domains": len(expected),
                "n_hcseeker_hotspot_transcripts": 2, "n_hcseeker_hotspot_residues": 4,
                "n_rmc_constrained_transcripts": 1, "n_rmc_constrained_residues": 3}})

    def test_absent_key_uses_original_default(self):
        self.build({"ENSG1:Pfam:PF00001"})

    def test_null_key_uses_default_instead_of_type_error(self):
        self.config["all_intolerant_domains"] = None
        self.build({"ENSG1:Pfam:PF00001"})

    def test_empty_key_uses_default_instead_of_missing_file_exit(self):
        self.config["all_intolerant_domains"] = ""
        self.build({"ENSG1:Pfam:PF00001"})

    def test_explicit_file_is_used_without_default_substitution(self):
        selected = self.root / "custom domains.pkl"
        selected.write_bytes(pickle.dumps({"ENSG2:Pfam:PF00002"}))
        self.config["all_intolerant_domains"] = str(selected)
        self.build({"ENSG2:Pfam:PF00002"})

    def test_missing_explicit_file_fails_even_when_default_exists(self):
        self.config["all_intolerant_domains"] = str(self.root / "missing.pkl")
        Path(self.args.config).write_text(yaml.safe_dump(self.config))
        with self.assertLogs(pm1.logger, level="ERROR"), self.assertRaises(SystemExit) as raised:
            pm1.prepare_pm1_regions(self.args)
        self.assertEqual(raised.exception.code, 1)
        self.assertFalse(Path(self.args.output).exists())


if __name__ == "__main__":
    unittest.main()
