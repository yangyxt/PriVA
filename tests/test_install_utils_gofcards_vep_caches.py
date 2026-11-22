"""GoFCards prepares exact-release official caches only when its build needs them."""
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import fake_priva_binary

ROOT = Path(__file__).resolve().parents[1]
DATA_TMP = Path(os.environ.get("PRIVA_TEST_TMPDIR", tempfile.gettempdir())) / "priva-gofcards-cache-tests"
SPEC = importlib.util.spec_from_file_location("gofcards_prerequisites", ROOT / "scripts/installer_prerequisites.py")
PREREQS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PREREQS)


class GoFCardsMergedCaches(unittest.TestCase):
    def setUp(self):
        DATA_TMP.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(prefix=ROOT.name + "_", dir=DATA_TMP)
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.cache = self.folder / "cache"
        self.plugins = self.folder / "plugins"
        self.plugin_cache = self.folder / "plugin-cache"
        for path in (self.cache, self.plugins, self.plugin_cache):
            path.mkdir()
        self.references = {}
        for assembly, length in (("hg19", 249250621), ("hg38", 248956422)):
            path = self.folder / (assembly + ".fa")
            path.write_text(">chr1\nA\n")
            Path(str(path) + ".fai").write_text(f"chr1\t{length}\t6\t1\t2\n")
            self.references[assembly] = path
        self.chain = self.folder / "hg19ToHg38.over.chain"
        self.chain.write_text("chain 1 chr1 249250621 + 0 1 chr1 248956422 + 0 1 1\n1\n")
        self.hgnc = self.folder / "hgnc.tsv"
        self.hgnc.write_text("fixture\n")
        self.config = {"assembly": "hg38", "ref_genome": str(self.references["hg38"]),
                       "vep_installed_version": "116", "vep_cache_dir": str(self.cache),
                       "vep_plugins_dir": str(self.plugins), "vep_plugins_cachedir": str(self.plugin_cache),
                       "gofcards_hg19_fasta": str(self.references["hg19"]),
                       "gofcards_hg38_fasta": str(self.references["hg38"]),
                       "gofcards_hg19_to_hg38_chain": str(self.chain),
                       "gofcards_workdir": str(self.folder / "work"), "gofcards_exact_gof_cache": str(self.folder / "fresh_deployed.json.gz"), "hgnc_table": str(self.hgnc)}
        self.config_path = self.folder / "config.json"
        self.events = self.folder / "events.jsonl"
        self.api = self.folder / "fake-api.py"
        self.api.write_text('''import json,os,pathlib,sys
args=sys.argv[1:]
def option(name): return args[args.index(name)+1]
assembly=option('--VEP_ASSEMBLY')
assembly={'hg19':'GRCh37','hg38':'GRCh38'}.get(assembly,assembly)
release=option('--VEP_VERSION')
with open(os.environ['TEST_EVENTS'],'a') as f: f.write(json.dumps({'event':'CACHE','assembly':assembly,'release':release})+'\\n')
if os.environ.get('TEST_FAIL_ASSEMBLY')==assembly: raise SystemExit(73)
destination=pathlib.Path(option('--VEP_CACHEDIR'))/'homo_sapiens_merged'/(os.environ.get('TEST_WRONG_RELEASE',release)+'_'+assembly)
(destination/'1').mkdir(parents=True,exist_ok=True)
(destination/'info.txt').write_text('species\\thomo_sapiens\\nassembly\\t'+assembly+'\\nvar_type\\ttabix\\n')
for name in ('chr_synonyms.txt','1/all_vars.gz','1/all_vars.gz.csi'): (destination/name).write_text('fixture\\n')
''')
        self.env = dict(os.environ, TEST_PYTHON=sys.executable, TEST_CONFIG=str(self.config_path),
                        TEST_EVENTS=str(self.events), TEST_API=str(self.api), PYTHONDONTWRITEBYTECODE="1",
                        PRIVA_INSTALL_INSPECT_ROOTS=os.pathsep.join((str(ROOT), str(self.folder))),
                        # The build route checks its inputs with priva installer-prerequisites.
                        **fake_priva_binary.environment(self.folder))
        for key in ("PRIVA_SMOKE_SKIP_VEP", "PRIVA_FORCE_GOFCARDS_CACHE", "PRIVA_FORCE_ALL_CACHES"):
            self.env.pop(key, None)


    def test_standalone_validator_cannot_silently_defer_missing_caches(self):
        with self.assertRaisesRegex(ValueError, "merged VEP 116 cache"):
            PREREQS.validate_gofcards(self.config, main_preflight=True)


if __name__ == "__main__":
    unittest.main()
