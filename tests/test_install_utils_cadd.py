"""CADD_install: upstream CADD-scripts and data installed the original way.

Fixtures use real files and archives; only the network (wget) and the
position-index builder are stubbed. No CADD data is downloaded and no CADD
model is run.
"""
import hashlib
import io
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile
import unittest
import zipfile

ROOT = Path(__file__).resolve().parents[1]
INSTALLER = ROOT / "scripts/install_utils.sh"
UPSTREAM = {
    "src/scripts/extract_scored.py": b"# fixture upstream lookup\n",
    "src/scripts/lib/tools/MMSplice.py": b"# fixture upstream MMSplice\n",
}

SHELL = r'''
source "$1" >/dev/null || exit
log() { printf '%s\n' "$*" >&2; }
priva_configure_cache_paths() { return 0; }
update_or_append_yaml() { printf '%s=%s\n' "$2" "$3" >> "$FIXTURE/config-updates.txt"; }
read_yaml() {
    local value
    value=$(awk -F '\t' -v key="$2" '$1 == key {print $2}' "$FIXTURE/settings.tsv")
    printf '%s\n' "${value:-null}"
}
wget() {
    # wget -c URL -O TARGET: serve the fixture named after the URL, appending
    # to TARGET like a resumed download.
    printf '%s\n' "$2" >> "$FIXTURE/downloads.txt"
    local served="$FIXTURE/served/$(basename "$2")"
    [[ -f $served ]] || return 8
    cat "$served" >> "$4"
}
python() {
    printf '%s\n' "$*" >> "$FIXTURE/python-calls.txt"
    case ${1##*/} in
        cadd_scripts_install.py) command python3 "$@" ;;
        cadd_prescore_index.py)
            if [[ $2 == check ]]; then
                [[ -f $3.posidx ]]
            else
                [[ -z ${INDEX_FAIL:-} ]] && printf 'fixture index\n' > "$3.posidx"
            fi ;;
        *) return 97 ;;
    esac
}
CADD_install "$FIXTURE/config.yaml" hg38
'''


def md5(payload):
    return hashlib.md5(payload).hexdigest()


class CaddInstall(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="cadd_install_", dir=os.environ.get("PRIVA_TEST_TMPDIR"))
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "served").mkdir()
        self.parent = self.root / "CADD"
        self.base = self.parent / "CADD-scripts-1.7.2"
        self.prescored = self.base / "data/prescored/GRCh38_v1.7/incl_anno"
        self.settings = {"cadd_version": "v1.7", "cadd_parent_dir": str(self.parent), "threads": "4",
                         "cadd_zip_download_url": "https://example.org/v1.7.2.zip",
                         "cadd_GRCh38_anno_url": "https://example.org/GRCh38_v1.7.tar.gz"}
        for kind, name in (("snv", "whole_genome_SNVs_inclAnno.tsv.gz"),
                           ("indel", "gnomad.genomes.r4.0.indel_inclAnno.tsv.gz")):
            data, index = ("%s rows\n" % kind).encode(), ("%s index\n" % kind).encode()
            self.serve(name, data)
            self.serve(name + ".tbi", index)
            self.settings.update({"cadd_GRCh38_%s_anno_url" % kind: "https://example.org/" + name,
                                  "cadd_GRCh38_%s_anno_md5" % kind: md5(data),
                                  "cadd_GRCh38_%s_anno_tbi_md5" % kind: md5(index)})
        annotations = io.BytesIO()
        with tarfile.open(fileobj=annotations, mode="w:gz") as archive:
            for name, payload in (("GRCh38_v1.7/vep/info.txt", b"cache\n"), ("GRCh38_v1.7/mmsplice/genes.gtf", b"gtf\n")):
                member = tarfile.TarInfo(name)
                member.size = len(payload)
                archive.addfile(member, io.BytesIO(payload))
        self.serve("GRCh38_v1.7.tar.gz", annotations.getvalue())
        self.settings["cadd_GRCh38_anno_md5"] = md5(annotations.getvalue())
        release = io.BytesIO()
        with zipfile.ZipFile(release, "w") as archive:
            script = zipfile.ZipInfo("CADD-scripts-1.7.2/CADD.sh")
            script.external_attr = 0o755 << 16
            archive.writestr(script, "#!/bin/sh\n")
            for relative, payload in UPSTREAM.items():
                archive.writestr("CADD-scripts-1.7.2/" + relative, payload)
        self.serve("v1.7.2.zip", release.getvalue())

    def serve(self, name, payload):
        (self.root / "served" / name).write_bytes(payload)

    def install(self, **environment):
        (self.root / "settings.tsv").write_text("".join("%s\t%s\n" % item for item in self.settings.items()))
        (self.root / "config.yaml").write_text("fixture: true\n")
        # Fixture upstream bytes stand in for the reviewed CADD files.
        blobs = {relative: hashlib.sha1(b"blob %d\0" % len(payload) + payload).hexdigest()
                 for relative, payload in UPSTREAM.items()}
        helper = self.root / "cadd_scripts_install.py"
        source = (ROOT / "scripts/cadd_scripts_install.py").read_text()
        source = source.replace("e88f103891da518a150d704f50d6bd3752970613", blobs["src/scripts/extract_scored.py"])
        source = source.replace("d9cf8353b7a23874101bbee8f978e18210d705f0", blobs["src/scripts/lib/tools/MMSplice.py"])
        source = source.replace('Path(__file__).resolve().parent / "cadd_scripts"',
                                "Path(%r)" % str(ROOT / "scripts/cadd_scripts"))
        helper.write_text(source)
        script = SHELL.replace('CADD_install "$FIXTURE', 'SCRIPT_DIR="$FIXTURE"\nCADD_install "$FIXTURE')
        result = subprocess.run(["bash", "--noprofile", "--norc", "-c", script, "test", str(INSTALLER)],
                                env=dict(os.environ, FIXTURE=str(self.root), **environment),
                                capture_output=True, text=True, timeout=60)
        return result

    def lines(self, name):
        path = self.root / name
        return path.read_text().splitlines() if path.exists() else []

    def index_builds(self):
        return [line for line in self.lines("python-calls.txt") if "cadd_prescore_index.py build" in line]

    def test_fresh_install_downloads_unpacks_installs_copies_and_indexes(self):
        result = self.install()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self.base / "CADD.sh").stat().st_mode & 0o100)
        self.assertEqual((self.base / "data/annotations/GRCh38_v1.7/vep/info.txt").read_bytes(), b"cache\n")
        self.assertFalse((self.base / "data/annotations/GRCh38_v1.7.tar.gz").exists())
        self.assertEqual(sorted(p.name for p in self.prescored.iterdir()), [
            "gnomad.genomes.r4.0.indel_inclAnno.tsv.gz", "gnomad.genomes.r4.0.indel_inclAnno.tsv.gz.posidx",
            "gnomad.genomes.r4.0.indel_inclAnno.tsv.gz.tbi", "whole_genome_SNVs_inclAnno.tsv.gz",
            "whole_genome_SNVs_inclAnno.tsv.gz.posidx", "whole_genome_SNVs_inclAnno.tsv.gz.tbi"])
        for relative, (copy, _) in {"src/scripts/extract_scored.py": ("extract_scored.py", None),
                                    "src/scripts/lib/tools/MMSplice.py": ("MMSplice.py", None)}.items():
            self.assertEqual((self.base / relative).read_bytes(), (ROOT / "scripts/cadd_scripts" / copy).read_bytes())
        self.assertEqual(self.lines("config-updates.txt"), ["cadd_base_dir=" + str(self.base)])
        self.assertEqual(len(self.index_builds()), 2)
        for call in self.index_builds():
            self.assertTrue(call.endswith(".tsv.gz --threads 4"), call)
        self.assertEqual([p.name for p in self.parent.iterdir()], [self.base.name])

    def test_complete_install_downloads_nothing_and_never_rehashes(self):
        self.assertEqual(self.install().returncode, 0)
        (self.root / "downloads.txt").unlink()
        (self.root / "python-calls.txt").unlink()
        for name in ("whole_genome_SNVs_inclAnno.tsv.gz", "gnomad.genomes.r4.0.indel_inclAnno.tsv.gz"):
            # Different bytes under the same name: a rerun must not hash or replace them.
            (self.prescored / name).write_bytes(b"changed after install\n")
        result = self.install()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.lines("downloads.txt"), [])
        self.assertEqual(self.index_builds(), [])

    def test_resumed_download_is_published_only_after_its_md5_matches(self):
        self.assertEqual(self.install().returncode, 0)
        target = self.prescored / "whole_genome_SNVs_inclAnno.tsv.gz"
        target.unlink()
        Path(str(target) + ".posidx").unlink()
        # A stopped download left "snv "; the resumed request serves the rest.
        self.serve(target.name, b"rows\n")
        Path(str(target) + ".part").write_bytes(b"snv ")
        result = self.install()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(target.read_bytes(), b"snv rows\n")
        self.assertFalse(Path(str(target) + ".part").exists())

    def test_md5_mismatch_publishes_nothing_and_keeps_the_partial_file(self):
        self.settings["cadd_GRCh38_snv_anno_md5"] = "0" * 32
        result = self.install()
        self.assertNotEqual(result.returncode, 0)
        target = self.prescored / "whole_genome_SNVs_inclAnno.tsv.gz"
        self.assertFalse(target.exists())
        self.assertTrue(Path(str(target) + ".part").exists())
        self.assertIn("differs from", result.stderr)

    def test_bad_annotation_archive_publishes_no_annotation_folder(self):
        self.serve("GRCh38_v1.7.tar.gz", b"not a tar archive")
        self.settings["cadd_GRCh38_anno_md5"] = md5(b"not a tar archive")
        result = self.install()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self.base / "data/annotations/GRCh38_v1.7").exists())

    def test_existing_folder_without_cadd_sh_is_retained(self):
        self.base.mkdir(parents=True)
        (self.base / "notes.txt").write_text("user files")
        self.settings["cadd_base_dir"] = str(self.base)
        result = self.install()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(sorted(p.name for p in self.base.iterdir()), ["notes.txt"])
        self.assertEqual(self.lines("downloads.txt"), [])

    def test_failed_index_stops_the_install(self):
        result = self.install(INDEX_FAIL="1")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(sorted(p.name for p in self.prescored.glob("*.posidx")), [])
        self.assertIn("Failed to index", result.stderr)


if __name__ == "__main__":
    unittest.main()
