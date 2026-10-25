"""Tiny real InterPro parser/publication checks and labelled fault injection.

No main installation or publisher acquisition occurs. Synthetic XML/VCF fixtures
exercise real helpers. Only explicit I/O/per-entry failures and a shell output
postcondition boundary are injected; exact successful bytes use a frozen parser.
"""
import contextlib
import gzip
import importlib.util
import io
import os
from pathlib import Path
import pickle
import runpy
import stat
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
PARSER = ROOT / "scripts/protein_domain_mapping.py"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


DOMAIN = module("interpro_atomic_subject", PARSER)
LEGACY = module("interpro_atomic_frozen", ROOT / "tests/fixtures/interpro_parser_before_atomic_20261018.py")
PATHS = module("interpro_atomic_paths", ROOT / "scripts/installer_config_paths.py")
XML = b'''<interprodb><release version="synthetic"/>
<interpro id="IPR000001" type="Domain" short_name="First"><name>First name</name>
<class_list><classification class_type="GO" id="GO:0000001"><category>category</category><description>description</description></classification></class_list>
<abstract>First <p>abstract</p> text.</abstract><member_list>
<db_xref db="Pfam" dbkey="PF00001"/><db_xref db="Pfam" dbkey="PF00001"/>
<db_xref db="SMART" dbkey="SM00001"/><db_xref db="Pfam"/>
</member_list></interpro>
<interpro id="IPR000002"><name>Second name</name><member_list><db_xref db="PFAM" dbkey="PF00001"/></member_list></interpro>
<interpro><member_list><db_xref db="PFAM" dbkey="NORMAL_OMISSION"/></member_list></interpro>
</interprodb>'''


class InterProAtomicInstall(unittest.TestCase):
    def setUp(self):
        retained = os.environ.get("PRIVA_INTERPRO_FIXTURE_ROOT")
        if retained:
            self.folder = Path(retained) / self._testMethodName
            self.folder.mkdir(parents=True, exist_ok=False)
        else:
            temporary = tempfile.TemporaryDirectory(prefix="interpro_atomic_fixture_", dir=os.environ.get(
                "PRIVA_TEST_TMPDIR", os.environ.get("TMPDIR")))
            self.addCleanup(temporary.cleanup)
            self.folder = Path(temporary.name)
        self.xml = self.folder / "interpro.xml.gz"
        self.write_xml(XML)
        self.vcf = self.folder / "synthetic.vcf"
        self.vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n'
                            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Synthetic fixture. Format: Allele|DOMAINS">\n'
                            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                            'chr1\t1\t.\tA\tC\t.\t.\tCSQ=C|Pfam:PF00001\n')
        self.output = self.folder / "mapping.pkl.gz"
        self.report = self.folder / "comparison.txt"

    def write_xml(self, content):
        with gzip.open(self.xml, "wb") as handle:
            handle.write(content)

    def snapshot(self, path):
        if not path.exists():
            return None
        value = path.stat()
        return path.read_bytes(), value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns, stat.S_IMODE(value.st_mode)

    def command(self, output=None):
        return [sys.executable, str(PARSER), "--xml", str(self.xml), "--mapping_output", str(output or self.output),
                "--vcf", str(self.vcf), "--report", str(self.report)]

    def cli(self):
        result = subprocess.run(self.command(), text=True, capture_output=True, timeout=30)
        (self.folder / "cli.stdout.log").write_text(result.stdout)
        (self.folder / "cli.stderr.log").write_text(result.stderr)
        return result

    def run_injected_cli(self, output=None):
        captured = io.StringIO()
        with mock.patch.object(sys, "argv", self.command(output)[1:]), \
                contextlib.redirect_stderr(captured), contextlib.redirect_stdout(captured):
            with self.assertRaises(SystemExit) as error:
                runpy.run_path(str(PARSER), run_name="__main__")
        self.assertNotEqual(error.exception.code, 0)
        (self.folder / "injected-cli.log").write_text(captured.getvalue())

    def assert_no_temporary(self):
        self.assertEqual(list(self.folder.glob(".*.write-*")), [])

    def shell(self, parser_boundary=""):
        config = self.folder / "config.yaml"
        text = PATHS.render_config((ROOT / "config.yaml").read_text(), {
            "cache_root": str(self.folder), "tmp_dir": str(self.folder), "assembly": "hg38",
            "vep_installed_version": "116", "interpro_mapping_pickle": str(self.output),
            "interpro_xml_gz_url": "https://fixtures.invalid/interpro.xml.gz", "alphamissense_vep_vcf": str(self.vcf)})
        config.write_text(text)
        before = config.read_bytes()
        script = r'''
source "$1"
wget() { printf 'UNEXPECTED_NETWORK\n' >&2; return 99; }
curl() { printf 'UNEXPECTED_NETWORK\n' >&2; return 99; }
''' + parser_boundary + '\nInterPro_parsing "$2"\n'
        result = subprocess.run(["bash", "--noprofile", "--norc", "-c", script, "fixture",
                                 str(ROOT / "scripts/install_utils.sh"), str(config)],
                                text=True, capture_output=True, timeout=30)
        (self.folder / "shell.stdout.log").write_text(result.stdout)
        (self.folder / "shell.stderr.log").write_text(result.stderr)
        self.assertNotIn("UNEXPECTED_NETWORK", result.stdout + result.stderr)
        return result, config, before

    def test_success_has_exact_frozen_bytes_protocol_order_and_multiplicity(self):
        original = LEGACY.parse_interpro_xml_for_mapping(str(self.xml))
        current = DOMAIN.parse_interpro_xml_for_mapping(str(self.xml))
        self.assertEqual(list(current), [("PFAM", "PF00001"), ("SMART", "SM00001")])
        self.assertEqual(len(current[("PFAM", "PF00001")]), 3)
        self.assertEqual([value[0] for value in current[("PFAM", "PF00001")]],
                         ["IPR000001", "IPR000001", "IPR000002"])
        self.assertEqual(current, original)
        expected = pickle.dumps(original, protocol=pickle.HIGHEST_PROTOCOL)
        (self.folder / "frozen.expected.pkl").write_bytes(expected)
        DOMAIN._write_interpro_mapping_atomic(current, self.output)
        self.assertEqual(self.output.read_bytes(), expected)
        self.assertEqual(self.output.read_bytes()[:2], bytes((0x80, pickle.HIGHEST_PROTOCOL)))
        result = self.cli()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.output.read_bytes(), expected)
        self.assert_no_temporary()

    def test_empty_full_mapping_fails_cli_and_preserves_old_final(self):
        self.write_xml(b"<interprodb><empty/></interprodb>")
        self.output.write_bytes(b"previous final fixture")
        before = self.snapshot(self.output)
        result = self.cli()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("No usable InterPro member mappings", result.stderr)
        self.assertEqual(self.snapshot(self.output), before)
        self.assert_no_temporary()

    def test_empty_full_mapping_without_old_final_stays_absent(self):
        self.write_xml(b"<interprodb/>")
        result = self.cli()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.output.exists())
        self.assert_no_temporary()

    def test_malformed_full_mapping_fails_and_keeps_old_final(self):
        self.write_xml(b"<interprodb><broken>")
        self.output.write_bytes(b"previous final fixture")
        before = self.snapshot(self.output)
        with self.assertRaises(DOMAIN.ET.ParseError):
            DOMAIN.parse_interpro_xml_for_mapping(str(self.xml))
        result = self.cli()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.snapshot(self.output), before)
        self.assert_no_temporary()

    def test_missing_xml_raises_instead_of_returning_empty_mapping(self):
        with self.assertRaises(FileNotFoundError):
            DOMAIN.parse_interpro_xml_for_mapping(str(self.folder / "absent.xml.gz"))

    def test_unexpected_entry_exception_cannot_return_partial_mapping(self):
        class BrokenEntry:
            tag = "interpro"
            def get(self, key):
                if key == "id":
                    return "BROKEN_FIXTURE"
                raise RuntimeError("injected unexpected entry failure")
            def clear(self):
                pass

        good = DOMAIN.ET.fromstring('<interpro id="IPR1"><member_list><db_xref db="PFAM" dbkey="PF1"/></member_list></interpro>')
        events = iter([("end", DOMAIN.ET.Element("release")), ("end", good), ("end", BrokenEntry())])
        with mock.patch.object(DOMAIN.ET, "iterparse", return_value=events):
            with self.assertRaisesRegex(RuntimeError, "injected unexpected entry"):
                DOMAIN.parse_interpro_xml_for_mapping(str(self.xml))

    def test_fatal_second_parse_failure_propagates_cli_before_publication(self):
        self.output.write_bytes(b"previous final fixture")
        before = self.snapshot(self.output)
        actual = DOMAIN.ET.iterparse
        calls = []
        def fail_second(*args, **kwargs):
            calls.append(1)
            if len(calls) == 2:
                raise OSError("injected full-map read failure")
            return actual(*args, **kwargs)
        with mock.patch.object(DOMAIN.ET, "iterparse", side_effect=fail_second):
            self.run_injected_cli()
        self.assertEqual(len(calls), 2)
        self.assertEqual(self.snapshot(self.output), before)
        self.assert_no_temporary()

    def test_dump_flush_close_replace_failures_keep_final_or_absence_and_owned_cleanup(self):
        actual_open = open
        for failure in ("dump", "flush", "close", "replace"):
            for old in (False, True):
                with self.subTest(failure=failure, old=old):
                    destination = self.folder / f"{failure}-{old}.pkl"
                    if old:
                        destination.write_bytes(b"previous final fixture")
                        destination.chmod(0o640)
                    before = self.snapshot(destination)
                    retained = self.folder / ("." + destination.name + ".write-retained")
                    retained.mkdir()
                    marker = retained / "previous-attempt"
                    marker.write_bytes(b"retained older staging fixture")

                    def partial_dump(data, handle, protocol):
                        self.assertEqual(protocol, pickle.HIGHEST_PROTOCOL)
                        handle.write(b"partial pickle bytes")
                        raise OSError("injected partial dump failure")

                    class FaultyFile:
                        def __init__(self, *args, **kwargs):
                            self.handle = actual_open(*args, **kwargs)
                        def __enter__(self):
                            return self
                        def write(self, value):
                            return self.handle.write(value)
                        def flush(self):
                            if failure == "flush":
                                raise OSError("injected flush failure")
                            return self.handle.flush()
                        def __exit__(self, kind, value, traceback):
                            self.handle.close()
                            if failure == "close" and kind is None:
                                raise OSError("injected close failure")

                    with contextlib.ExitStack() as stack:
                        if failure == "dump":
                            stack.enter_context(mock.patch.object(DOMAIN.pickle, "dump", side_effect=partial_dump))
                        elif failure == "replace":
                            stack.enter_context(mock.patch.object(DOMAIN.os, "replace", side_effect=OSError("injected replace failure")))
                        else:
                            stack.enter_context(mock.patch.object(DOMAIN, "open", FaultyFile, create=True))
                        with self.assertRaises(OSError):
                            DOMAIN._write_interpro_mapping_atomic({("PFAM", "PF1"): [("fixture",)]}, destination)
                    self.assertEqual(self.snapshot(destination), before)
                    self.assertEqual(list(self.folder.glob("." + destination.name + ".write-*")), [retained])
                    self.assertEqual(marker.read_bytes(), b"retained older staging fixture")

    def test_cli_dump_and_replace_failures_are_nonzero_with_and_without_old_final(self):
        for failure in ("dump", "replace"):
            for old in (False, True):
                with self.subTest(failure=failure, old=old):
                    destination = self.folder / f"cli-{failure}-{old}.pkl"
                    if old:
                        destination.write_bytes(b"previous final fixture")
                    before = self.snapshot(destination)
                    def partial_dump(data, handle, protocol):
                        handle.write(b"partial pickle bytes")
                        raise OSError("injected CLI dump failure")
                    target = mock.patch.object(DOMAIN.pickle, "dump", side_effect=partial_dump) if failure == "dump" else \
                        mock.patch.object(DOMAIN.os, "replace", side_effect=OSError("injected CLI replace failure"))
                    with target:
                        self.run_injected_cli(destination)
                    self.assertEqual(self.snapshot(destination), before)
                    self.assert_no_temporary()

    def test_old_mode_and_new_umask_are_preserved_on_success(self):
        self.output.write_bytes(b"previous final fixture")
        self.output.chmod(0o640)
        DOMAIN._write_interpro_mapping_atomic({"fixture": [1, 2]}, self.output)
        self.assertEqual(stat.S_IMODE(self.output.stat().st_mode), 0o640)
        new = self.folder / "new.pkl"
        previous = os.umask(0o027)
        try:
            DOMAIN._write_interpro_mapping_atomic({"fixture": [1, 2]}, new)
        finally:
            os.umask(previous)
        self.assertEqual(stat.S_IMODE(new.stat().st_mode), 0o640)
        self.assert_no_temporary()

    def test_report_write_failure_is_nonzero_and_does_not_replace_mapping(self):
        self.output.write_bytes(b"previous final fixture")
        before = self.snapshot(self.output)
        self.report.mkdir()
        result = self.cli()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.snapshot(self.output), before)

    def test_newer_empty_external_mapping_is_rebuilt_by_real_parser(self):
        self.output.touch()
        timestamp = max(time.time_ns(), PARSER.stat().st_mtime_ns, self.xml.stat().st_mtime_ns) + 1000000
        os.utime(self.output, ns=(timestamp, timestamp))
        result, config, before = self.shell()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("already generated. Skip", result.stderr)
        self.assertEqual(self.output.read_bytes(), pickle.dumps(LEGACY.parse_interpro_xml_for_mapping(str(self.xml)),
                                                                protocol=pickle.HIGHEST_PROTOCOL))

    def test_stale_previous_final_and_config_stay_exact_on_real_parse_failure(self):
        self.write_xml(b"<interprodb><empty/></interprodb>")
        self.output.write_bytes(b"previous final fixture")
        os.utime(self.output, ns=(1, 1))
        before = self.snapshot(self.output)
        result, config, config_before = self.shell()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.snapshot(self.output), before)
        self.assertEqual(config.read_bytes(), config_before)

    def test_genuine_nonempty_bundle_is_reused_without_parser_or_acquisition(self):
        data = self.folder / "source_data"
        bundle = data / "InterPro/Interpro_entry_mapping.pkl.gz"
        bundle.parent.mkdir(parents=True)
        # Opaque on purpose: genuine bundles require existence/nonzero only.
        bundle.write_bytes(b"opaque original bundled snapshot fixture")
        config = self.folder / "config.yaml"
        config.write_text(PATHS.render_config((ROOT / "config.yaml").read_text(),
                                              {"interpro_mapping_pickle": str(bundle)}))
        before = self.snapshot(bundle), config.read_bytes()
        script = r'''
source "$1"
DATA_DIR=$3
python() { printf 'UNEXPECTED_PARSER\n'; return 99; }
wget() { printf 'UNEXPECTED_ACQUISITION\n'; return 99; }
InterPro_parsing "$2"
'''
        result = subprocess.run(["bash", "--noprofile", "--norc", "-c", script, "fixture",
                                 str(ROOT / "scripts/install_utils.sh"), str(config), str(data)],
                                text=True, capture_output=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("UNEXPECTED", result.stdout + result.stderr)
        self.assertEqual((self.snapshot(bundle), config.read_bytes()), before)

    def test_shell_rejects_missing_or_empty_output_from_zero_exit_boundary(self):
        for produced in ("missing", "empty"):
            with self.subTest(produced=produced):
                self.output.unlink(missing_ok=True)
                boundary = r'''
# Explicit output-postcondition fault injection; the real parser is covered above.
python() {
    if [[ "$1" != */protein_domain_mapping.py ]]; then command python "$@"; return; fi
''' + (r'''while [[ $# -gt 0 ]]; do
    if [[ "$1" == --mapping_output ]]; then : > "$2"; break; fi
    shift
done
''' if produced == "empty" else "") + '\nreturn 0\n}\n'
                result, config, before = self.shell(boundary)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("Missing or empty generated InterPro mapping", result.stderr)
                self.assertEqual(config.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
