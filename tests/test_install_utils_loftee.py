"""LoFtee installer failures/reuse with actual small gzip, SQLite and BigWig resources."""
import gzip
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import unittest

from test_conservation_resource_check import make_bigwig, make_reference


ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "scripts/install_utils.sh"
PRELUDE = r'''
source "$1"
log() { printf '%s\n' "$*" >&2; }
event() { printf '%s\n' "$*" >> "$TEST_EVENTS"; }
read_yaml() {
    "$TEST_PYTHON" -c 'import json,sys; print(json.load(open(sys.argv[1])).get(sys.argv[2], ""))' "$1" "$2"
}
update_yaml() {
    event "update $2"
    [[ "$TEST_FAIL" != update ]] || return 71
    "$TEST_PYTHON" -c 'import json,sys; p,k,v=sys.argv[1:]; d=json.load(open(p)); d[k]=v; open(p,"w").write(json.dumps(d))' "$@"
}
update_or_append_yaml() { update_yaml "$@"; }
priva_env_vars_set() { export "$@"; }
python() { "$TEST_PYTHON" "$@"; }
git() { event "git $*"; return 72; }
wget() {
    local destination=${@: -1} name fixture
    name=${destination##*/}
    case "$name" in
        human_ancestor.fa.gz) fixture=ancestor.fa.gz ;;
        human_ancestor.fa.gz.fai) fixture=ancestor.fa.gz.fai ;;
        phylocsf_gerp.sql.gz) fixture=conservation.sql.gz ;;
        gerp_conservation_scores.homo_sapiens.bw) fixture=gerp.bw ;;
        *) return 73 ;;
    esac
    event "download $name"
    if [[ "$TEST_FAIL" == "$name" ]]; then
        printf 'retained partial bytes\n' > "$destination"
        return 74
    fi
    cp "$TEST_ASSETS/$fixture" "$destination"
}
unset PERL5LIB
'''


class LoFteeInstallerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="priva_loftee_", dir=os.environ.get("PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.cache = self.folder / "cache"
        self.plugins = self.folder / "plugins"
        self.cache.mkdir()
        self.plugins.mkdir()
        self.events = self.folder / "events.txt"
        self.events.touch()
        self.config = self.folder / "config.json"
        self.env = dict(os.environ, TEST_EVENTS=str(self.events), TEST_PYTHON=sys.executable,
                        TEST_CONFIG=str(self.config), TEST_CACHE=str(self.cache), TEST_FAIL="none",
                        TMPDIR=str(self.folder))
        self.configure("hg19")

    def configure(self, assembly):
        self.assembly = assembly
        native = "hg19" if assembly in ("hg19", "GRCh37") else "hg38"
        length = 249250621 if native == "hg19" else 248956422
        self.assets = self.folder / (assembly + "_assets")
        self.assets.mkdir()
        self.ancestor = self.assets / "ancestor.fa.gz"
        with gzip.open(self.ancestor, "wb") as stream:
            stream.write(b">1\nA\n")
        Path(str(self.ancestor) + ".fai").write_text(f"1\t{length}\t3\t1\t2\n")
        self.sql = self.assets / "conservation.sql"
        with sqlite3.connect(self.sql) as database:
            database.execute("CREATE TABLE conservation (score REAL)")
            database.execute("INSERT INTO conservation VALUES (1.0)")
        with gzip.open(str(self.sql) + ".gz", "wb") as stream:
            stream.write(self.sql.read_bytes())
        self.reference = make_reference(self.assets, native)
        self.gerp = make_bigwig(self.assets, "gerp", {"1": length}) if native == "hg38" else None
        self.repository = self.plugins / "LoFtee_repos" / ("loftee-" + assembly)
        self.repository.mkdir(parents=True)
        (self.repository / "LoF.pm").write_text("# fixture plugin bytes\n")
        self.settings = {"assembly": assembly, "vep_plugins_dir": str(self.plugins),
                         "human_ancestor_fasta": str(self.ancestor), "ref_genome": str(self.reference),
                         "loftee_conservation_file": str(self.sql), "gerp_bigwig": str(self.gerp or "")}
        self.save()
        self.env["TEST_ASSETS"] = str(self.assets)

    def save(self):
        self.config.write_text(json.dumps(self.settings))

    def fresh(self):
        for name in ("human_ancestor_fasta", "loftee_conservation_file", "gerp_bigwig"):
            self.settings[name] = ""
        self.save()

    def run_shell(self, setup="", suffix=""):
        body = PRELUDE + setup + '\nLoFtee_install "$TEST_CACHE" "$TEST_CONFIG"' + suffix
        result = subprocess.run(["bash", "--noprofile", "--norc", "-u", "-c", body, "test", str(INSTALLER)],
                                text=True, capture_output=True, env=self.env, cwd=self.folder, timeout=60)
        return result, self.events.read_text().splitlines()

    def test_actual_resources_reuse_for_all_assembly_aliases_without_updates(self):
        for assembly in ("hg19", "GRCh37", "hg38", "GRCh38"):
            with self.subTest(assembly=assembly):
                if assembly != "hg19":
                    self.configure(assembly)
                before = self.config.read_bytes()
                result, events = self.run_shell()
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(events, [])
                self.assertEqual(self.config.read_bytes(), before)
                self.assertIn("already installed", result.stderr)

    def test_fresh_hg19_under_nounset_uses_original_three_sources_without_gerp_url(self):
        self.fresh()
        result, events = self.run_shell()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, ["download human_ancestor.fa.gz", "download human_ancestor.fa.gz.fai",
                                  "download phylocsf_gerp.sql.gz", "update human_ancestor_fasta",
                                  "update loftee_conservation_file", "update gerp_bigwig"])
        registered = json.loads(self.config.read_text())
        self.assertEqual(registered["human_ancestor_fasta"], str(self.cache / "LoFtee/hg19/human_ancestor.fa.gz"))
        self.assertEqual(registered["gerp_bigwig"], "")

    def test_fresh_grch38_validates_bigwig_before_registration(self):
        self.configure("GRCh38")
        self.fresh()
        result, events = self.run_shell()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("download gerp_conservation_scores.homo_sapiens.bw", events)
        self.assertEqual(events[-3:], ["update human_ancestor_fasta", "update loftee_conservation_file", "update gerp_bigwig"])
        self.assertIn('"scientific_qualified": false', result.stdout)

    def test_every_download_failure_propagates_and_retains_partial_without_registering(self):
        self.configure("hg38")
        self.fresh()
        for name in ("human_ancestor.fa.gz", "human_ancestor.fa.gz.fai", "phylocsf_gerp.sql.gz",
                     "gerp_conservation_scores.homo_sapiens.bw"):
            with self.subTest(name=name):
                shutil.rmtree(self.cache / "LoFtee", ignore_errors=True)
                self.events.write_text("")
                self.env["TEST_FAIL"] = name
                before = self.config.read_bytes()
                result, events = self.run_shell()
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(any(event.startswith("update ") for event in events))
                self.assertEqual(self.config.read_bytes(), before)
                stages = list((self.cache / "LoFtee/hg38").glob(".acquire.*/" + name))
                self.assertEqual(len(stages), 1)
                self.assertEqual(stages[0].read_bytes(), b"retained partial bytes\n")
                self.assertFalse((self.cache / "LoFtee/hg38" / name).exists())
                self.assertIn("explicit recovery", result.stderr)

    def test_sql_decompression_failure_is_not_masked_and_compressed_input_is_retained(self):
        self.fresh()
        result, events = self.run_shell('gunzip() { return 75; }')
        self.assertNotEqual(result.returncode, 0)
        retained = list((self.cache / "LoFtee/hg19").glob(".acquire.*/phylocsf_gerp.sql.gz"))
        self.assertEqual(len(retained), 1)
        self.assertEqual(retained[0].read_bytes(),
                         Path(str(self.sql) + ".gz").read_bytes())
        self.assertFalse(any(event.startswith("update ") for event in events))

    def test_successful_download_with_invalid_sql_body_is_retained_without_registration(self):
        with gzip.open(str(self.sql) + ".gz", "wb") as stream:
            stream.write(b"successful transfer is not SQLite\n")
        self.fresh()
        result, events = self.run_shell()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(any(event.startswith("update ") for event in events))
        retained = list((self.cache / "LoFtee/hg19").glob(".acquire.*/phylocsf_gerp.sql"))
        self.assertEqual(len(retained), 1)
        self.assertFalse((self.cache / "LoFtee/hg19/phylocsf_gerp.sql").exists())
        self.assertEqual(retained[0].read_bytes(),
                         b"successful transfer is not SQLite\n")

    def test_interrupted_acquisitions_leave_no_canonical_payload_and_retry_succeeds(self):
        self.configure("hg38")
        self.fresh()
        self.env["TEST_FAIL"] = "human_ancestor.fa.gz.fai"
        result, events = self.run_shell()
        self.assertNotEqual(result.returncode, 0)
        target = self.cache / "LoFtee/hg38"
        self.assertFalse((target / "human_ancestor.fa.gz").exists())
        self.assertFalse((target / "human_ancestor.fa.gz.fai").exists())
        self.assertFalse(any(event.startswith("update ") for event in events))
        failed_stages = list(target.glob(".acquire.*"))
        self.assertEqual(len(failed_stages), 1)
        self.env["TEST_FAIL"] = "none"
        self.events.write_text("")
        result, events = self.run_shell()
        self.assertEqual(result.returncode, 0, result.stderr)
        registered = json.loads(self.config.read_text())
        self.assertEqual(Path(registered["human_ancestor_fasta"]).read_bytes(), self.ancestor.read_bytes())
        self.assertEqual(Path(registered["loftee_conservation_file"]).read_bytes(), self.sql.read_bytes())
        self.assertEqual(list(target.glob(".acquire.*")), failed_stages)

    def test_sqlite_without_tables_is_rejected_without_acquisition(self):
        self.sql.unlink()
        with sqlite3.connect(self.sql) as database:
            database.execute("PRAGMA user_version=1")
        result, events = self.run_shell()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("SQLite has no tables", result.stderr)
        self.assertEqual(events, [])

    def test_git_clone_failure_stops_before_acquisition_and_registration(self):
        self.fresh()
        shutil.rmtree(self.repository)
        result, events = self.run_shell()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(len(events), 1)
        self.assertTrue(events[0].startswith("git clone https://github.com/konradjk/loftee.git "))

    def test_interrupted_repository_clone_keeps_canonical_path_free_and_retry_succeeds(self):
        self.configure("GRCh38")
        self.fresh()
        shutil.rmtree(self.repository)
        partial = r'''
        git() {
            event "git $*"
            [[ "$1" == clone ]] || return 73
            mkdir -p "${@: -1}"
            printf 'partial repository bytes\n' > "${@: -1}/LoF.pm"
            return 72
        }
        '''
        result, events = self.run_shell(partial)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.repository.exists())
        self.assertFalse(any(event.startswith("download ") or event.startswith("update ") for event in events))
        failed = list(self.repository.parent.glob(".loftee-GRCh38.clone.*"))
        self.assertEqual(len(failed), 1)
        self.assertEqual((failed[0] / "loftee/LoF.pm").read_text(), "partial repository bytes\n")
        self.events.write_text("")
        complete = r'''
        git() {
            event "git $*"
            if [[ "$1" == clone ]]; then
                mkdir -p "${@: -1}"
                printf '# complete selected plugin\n' > "${@: -1}/LoF.pm"
            fi
        }
        '''
        result, events = self.run_shell(complete)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.repository / "LoF.pm").read_text(), "# complete selected plugin\n")
        self.assertEqual(list(self.repository.parent.glob(".loftee-GRCh38.clone.*")), failed)
        self.assertEqual(json.loads(self.config.read_text())["loftee_repo"], str(self.repository))

    def test_clone_branch_setup_and_worktree_reuse_preserve_caller_directory(self):
        self.configure("GRCh38")
        self.fresh()
        shutil.rmtree(self.repository)
        self.env["TEST_REPOSITORY"] = str(self.repository)
        setup = r'''
        git() {
            event "git $*"
            if [[ "$1" == clone ]]; then
                mkdir -p "${@: -1}"
                printf '# fixture plugin bytes\n' > "${@: -1}/LoF.pm"
            fi
        }
        caller_directory=$PWD
        '''
        suffix = '\nstatus=$?\n[[ "$status" -eq 0 && "$PWD" == "$caller_directory" ]]'
        result, events = self.run_shell(setup, suffix)
        self.assertEqual(result.returncode, 0, result.stderr)
        git_events = [event for event in events if event.startswith("git ")]
        staged = git_events[0].removeprefix("git clone https://github.com/konradjk/loftee.git ")
        self.assertEqual(git_events[1:], [f"git -C {staged} checkout grch38", f"git -C {staged} pull"])
        self.assertFalse(Path(staged).exists())
        # Ancillary cache recovery must not mutate a selected worktree revision.
        self.events.write_text("")
        (self.repository / ".git").write_text("gitdir: fixture\n")
        self.fresh()
        result, events = self.run_shell(setup, suffix)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(any(event.startswith("git ") for event in events))

    def test_conventional_git_directory_repair_preserves_existing_plugin_bytes_without_pull(self):
        (self.repository / ".git").mkdir()
        before = (self.repository / "LoF.pm").read_bytes()
        self.fresh()
        result, events = self.run_shell()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(any(event.startswith("git ") for event in events))
        self.assertEqual((self.repository / "LoF.pm").read_bytes(), before)

    def test_worktree_revision_is_preserved_during_ancillary_recovery(self):
        (self.repository / ".git").write_text("gitdir: fixture\n")
        self.fresh()
        result, events = self.run_shell()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(any(event.startswith("git ") for event in events))

    def test_explicit_repository_is_used_without_creating_default_replacement(self):
        selected = self.folder / "selected LoFtee checkout"
        self.repository.rename(selected)
        self.settings["loftee_repo"] = str(selected)
        self.save()
        before = self.config.read_bytes(), (selected / "LoF.pm").read_bytes()
        result, events = self.run_shell()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(events, [])
        self.assertEqual((self.config.read_bytes(), (selected / "LoF.pm").read_bytes()), before)
        self.assertFalse(self.repository.exists())

    def test_missing_custom_repository_or_resource_stops_before_network(self):
        for key in ("loftee_repo", "human_ancestor_fasta", "loftee_conservation_file"):
            with self.subTest(key=key):
                saved = self.settings.get(key)
                self.settings[key] = str(self.folder / "missing custom input")
                self.save()
                before = self.config.read_bytes()
                result, events = self.run_shell()
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(events, [])
                self.assertEqual(self.config.read_bytes(), before)
                if saved is None:
                    self.settings.pop(key)
                else:
                    self.settings[key] = saved

    def test_existing_invalid_gzip_is_rejected_and_retained_without_network(self):
        self.ancestor.write_bytes(b"not gzip")
        result, events = self.run_shell()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, [])
        self.assertEqual(self.ancestor.read_bytes(), b"not gzip")

    def test_wrong_assembly_or_malformed_ancestor_index_is_rejected(self):
        index = Path(str(self.ancestor) + ".fai")
        for content in ("1\t248956422\t3\t1\t2\n", "bad index\n"):
            with self.subTest(content=content):
                index.write_text(content)
                result, events = self.run_shell()
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(events, [])
                self.assertEqual(index.read_text(), content)

    def test_explicit_ancestor_missing_index_does_not_get_guessed_default_index(self):
        Path(str(self.ancestor) + ".fai").unlink()
        result, events = self.run_shell()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, [])
        self.assertIn("Supply the matching index", result.stderr)

    def test_invalid_sqlite_or_empty_sqlite_is_rejected_without_replacement(self):
        for content in (b"not sqlite", b"SQLite format 3\x00truncated"):
            with self.subTest(content=content):
                self.sql.write_bytes(content)
                result, events = self.run_shell()
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(events, [])
                self.assertEqual(self.sql.read_bytes(), content)

    def test_invalid_or_wrong_assembly_grch38_bigwig_is_rejected_without_registration(self):
        self.configure("hg38")
        wrong = make_bigwig(self.assets, "wrong", {"1": 249250621})
        for content in (b"not BigWig", wrong.read_bytes()):
            with self.subTest(bytes=len(content)):
                self.gerp.write_bytes(content)
                result, events = self.run_shell()
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(events, [])
                self.assertEqual(self.gerp.read_bytes(), content)

    def test_configuration_write_failure_propagates_after_validating_resources(self):
        self.fresh()
        self.env["TEST_FAIL"] = "update"
        before = self.config.read_bytes()
        result, events = self.run_shell()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events[-1], "update human_ancestor_fasta")
        self.assertEqual(self.config.read_bytes(), before)

    def test_plugin_cache_directory_failure_stops_before_downloads(self):
        (self.cache / "LoFtee").write_text("blocks directory creation\n")
        result, events = self.run_shell()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(events, [])


if __name__ == "__main__":
    unittest.main()
