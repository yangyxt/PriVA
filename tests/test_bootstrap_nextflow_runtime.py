"""Driver installer contracts using tiny labelled fake downloads and JVM assets.

Only a temporary copy's URL/checksum constants are changed. No production pin
override, real network acquisition, Java execution or scientific workflow is used.
Run with a data-side TMPDIR/cwd; fixtures are removed after each completed test.
"""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tarfile
import tempfile
import unittest


SOURCE = Path(__file__).resolve().parents[1] / "scripts/bootstrap_nextflow_runtime.sh"
JAVA_DIRECTORY = "jdk-21.0.12.1+1-jre"
NF_FILENAME = "nextflow-25.10.4-dist"


class BootstrapRuntimeTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="fake-nextflow-bootstrap-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.prefix = self.root / "runtime with spaces"
        self.stage = Path(str(self.prefix) + ".install")
        self.requests = self.root / "fake-curl-requests.jsonl"
        self.fake_bin = self.root / "fake-bin"
        self.fake_bin.mkdir()
        self.payloads = self.root / "tiny-fake-downloads"
        self.payloads.mkdir()
        self.scratch = self.root / "configured temporary directory"
        self.scratch.mkdir()
        self.helper = self.root / "bootstrap-copy.sh"
        self.environment = dict(os.environ, PATH=str(self.fake_bin) + os.pathsep + os.environ["PATH"],
                                FAKE_REQUESTS=str(self.requests), FAKE_PREFIX=str(self.prefix),
                                FAKE_INTERRUPTED=str(self.root / "interrupted"),
                                TMPDIR=str(self.scratch), FAKE_EXPECTED_TMPDIR=str(self.scratch),
                                PYTHONDONTWRITEBYTECODE="1")
        self.write_curl()
        self.write_payloads()

    def write_curl(self):
        # This executable is a labelled transport stub. It never opens a socket.
        script = self.fake_bin / "curl"
        script.write_text("""#!/usr/bin/env python3
import json, os
from pathlib import Path
import sys
from urllib.parse import unquote, urlparse
args = sys.argv[1:]
assert '--continue-at' in args and args[args.index('--continue-at') + 1] == '-'
assert args[args.index('--connect-timeout') + 1] == '30'
assert args[args.index('--speed-limit') + 1] == '1024'
assert args[args.index('--speed-time') + 1] == '120'
assert not os.path.lexists(os.environ['FAKE_PREFIX']), 'Published before validation'
url, destination = args[-3], Path(args[-1])
offset = destination.stat().st_size if destination.exists() else 0
with open(os.environ['FAKE_REQUESTS'], 'a') as requests:
    requests.write(json.dumps({'url': url, 'offset': offset}) + '\\n')
if os.environ.get('FAKE_NETWORK_DENIED'):
    sys.exit(97)
assert urlparse(url).scheme == 'file', 'Only tiny fixture files are permitted'
data = Path(unquote(urlparse(url).path)).read_bytes()
mode = os.environ.get('FAKE_CURL_MODE', '')
interrupted = Path(os.environ['FAKE_INTERRUPTED'])
if mode == 'corrupt':
    data = b'deliberately corrupt fake download\\n'
if mode in ('interrupt-once', 'valid-then-interrupt') and not interrupted.exists():
    destination.write_bytes(data if mode == 'valid-then-interrupt' else data[:len(data)//2])
    interrupted.write_text('The fake transfer was interrupted, exit18.\\n')
    sys.exit(18)
if offset > len(data):
    sys.exit(33)
with destination.open('ab') as output:
    output.write(data[offset:])
""")
        script.chmod(0o755)

    def write_payloads(self, *, java_version="21.0.12.1", nf_version="25.10.4",
                       hang_java=False, hang_nextflow=False):
        distribution = self.payloads / NF_FILENAME
        distribution.write_text("""#!/usr/bin/env bash
# Tiny fake Nextflow: assert the actual production probe stays local and offline.
set -euo pipefail
runtime=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
[[ $# == 1 && $1 == -version && $NXF_OFFLINE == true ]]
[[ $NXF_HOME == "$runtime/nxf-home" && $PWD == "$NXF_HOME" ]]
[[ $TMPDIR == "${FAKE_EXPECTED_TMPDIR:-$NXF_HOME}" ]]
[[ $NXF_TEMP == "$TMPDIR" && $XDG_CACHE_HOME == "$NXF_HOME" ]]
[[ $JAVA_HOME == "$runtime/jdk-21.0.12.1+1-jre" && $NXF_JAVA_HOME == "$JAVA_HOME" ]]
[[ $NXF_VER == 25.10.4 && $NXF_SYNTAX_PARSER == v2 ]]
# FAKE_DELAY
printf '%s\\n' 'N E X T F L O W  ~  version VERSION build 11177 (tiny fake asset)'
""".replace("VERSION", nf_version).replace("# FAKE_DELAY", "sleep 60" if hang_nextflow else ":"))
        tree = self.payloads / JAVA_DIRECTORY
        (tree / "bin").mkdir(parents=True, exist_ok=True)
        (tree / "lib").mkdir(exist_ok=True)
        java = tree / "bin/java"
        java.write_text("""#!/usr/bin/env bash
# Tiny fake JVM: version output only, no JVM is launched.
set -euo pipefail
[[ ${!#} == -version ]]
runtime=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
expected_tmpdir=${FAKE_EXPECTED_TMPDIR:-$runtime/nxf-home}
[[ $TMPDIR == "$expected_tmpdir" ]]
found_tmpdir=false
for argument in "$@"; do
    [[ $argument != "-Djava.io.tmpdir=$expected_tmpdir" ]] || found_tmpdir=true
done
[[ $found_tmpdir == true ]]
# FAKE_DELAY
printf '%s\\n' 'openjdk version "VERSION" (tiny fake asset)' 'Temurin-21.0.12.1+1' >&2
""".replace("VERSION", java_version).replace("# FAKE_DELAY", "sleep 60" if hang_java else ":"))
        java.chmod(0o755)
        (tree / "release").write_text("Clearly labelled tiny fake JRE release metadata.\n")
        (tree / "lib/fake-library.bin").write_bytes(b"tiny fake JVM library\n")
        link = tree / "lib/release-link"
        if not link.is_symlink():
            link.symlink_to("../release")
        archive = self.payloads / "temurin-jre.tar.gz"
        with tarfile.open(archive, "w:gz") as stream:
            stream.add(tree, arcname=JAVA_DIRECTORY)
        text = SOURCE.read_text()
        for key, value in {"nf_url": distribution.as_uri(), "nf_sha": self.digest(distribution),
                           "java_url": archive.as_uri(), "java_sha": self.digest(archive)}.items():
            text, count = re.subn(r"^" + key + r"='[^']*'$", lambda unused: key + "='" + value + "'",
                                 text, flags=re.MULTILINE)
            self.assertEqual(count, 1)
        self.helper.write_text(text)

    @staticmethod
    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def invoke(self, procedure="install", *, extra_env=None, legacy=False):
        environment = dict(self.environment, **(extra_env or {}))
        arguments = [str(self.prefix)] if legacy else [procedure, "--prefix", str(self.prefix)]
        return subprocess.run(["bash", str(self.helper), *arguments], cwd=self.root,
                              env=environment, capture_output=True, text=True, timeout=20)

    def require_success(self, result):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def recorded_requests(self):
        return [json.loads(line) for line in self.requests.read_text().splitlines()] if self.requests.exists() else []

    def prefix_snapshot(self):
        return {str(path.relative_to(self.prefix)): (path.lstat().st_ino, path.lstat().st_mtime_ns,
                os.readlink(path) if path.is_symlink() else self.digest(path))
                for path in self.prefix.rglob("*") if path.is_file() or path.is_symlink()}

    def test_cold_install_is_published_only_after_validation(self):
        self.require_success(self.invoke())
        self.assertEqual(len(self.recorded_requests()), 2)
        self.assertTrue((self.prefix / "bin/nextflow").is_file())
        self.assertFalse(list(self.root.glob(self.prefix.name + ".partial.*")))
        self.assertEqual((self.stage / NF_FILENAME).stat().st_ino,
                         (self.prefix / "bin" / NF_FILENAME).stat().st_ino)
        self.require_success(self.invoke("check", extra_env={"FAKE_NETWORK_DENIED": "1"}))

    def test_owned_and_old_layout_warm_reuse_need_no_network_or_rewrite(self):
        self.require_success(self.invoke())
        before = self.prefix_snapshot()
        self.require_success(self.invoke(extra_env={"FAKE_NETWORK_DENIED": "1"}))
        self.assertEqual(before, self.prefix_snapshot())
        # The exact old production layout predates the new ownership marker.
        (self.prefix / ".priva-nextflow-runtime").unlink()
        before = self.prefix_snapshot()
        self.require_success(self.invoke(extra_env={"FAKE_NETWORK_DENIED": "1"}, legacy=True))
        self.require_success(self.invoke("check", extra_env={"FAKE_NETWORK_DENIED": "1"}))
        self.assertEqual(before, self.prefix_snapshot())
        self.assertEqual(len(self.recorded_requests()), 2)

    def test_supplied_tmpdir_controls_probes_without_moving_persistent_state(self):
        # A stale caller NXF_TEMP must not override main_install's selected TMPDIR.
        self.environment["NXF_TEMP"] = str(self.root / "unselected scratch")
        self.require_success(self.invoke())
        self.require_success(self.invoke("check", extra_env={"FAKE_NETWORK_DENIED": "1"}))
        self.assertTrue((self.prefix / "nxf-home").is_dir())
        self.assertFalse((self.root / "unselected scratch").exists())
        self.assertEqual(list(self.scratch.iterdir()), [])

    def test_absent_tmpdir_keeps_standalone_runtime_fallback(self):
        self.environment.pop("TMPDIR", None)
        self.environment.pop("FAKE_EXPECTED_TMPDIR", None)
        self.environment["NXF_TEMP"] = str(self.root / "unselected scratch")
        self.require_success(self.invoke(legacy=True))
        self.require_success(self.invoke("check", extra_env={"FAKE_NETWORK_DENIED": "1"}))
        self.assertTrue((self.prefix / "nxf-home").is_dir())
        self.assertFalse((self.root / "unselected scratch").exists())

    def test_interrupted_download_resumes_existing_bytes(self):
        result = self.invoke(extra_env={"FAKE_CURL_MODE": "interrupt-once"})
        self.assertEqual(result.returncode, 18, result.stderr)
        self.assertFalse(self.prefix.exists())
        partial = (self.stage / NF_FILENAME).read_bytes()
        self.assertGreater(len(partial), 0)
        self.require_success(self.invoke())
        self.assertEqual(self.recorded_requests()[1]["offset"], len(partial))
        self.assertTrue((self.prefix / "bin" / NF_FILENAME).read_bytes().startswith(partial))

    def test_checksum_valid_bytes_survive_interruption_without_redownload(self):
        result = self.invoke(extra_env={"FAKE_CURL_MODE": "valid-then-interrupt"})
        self.assertEqual(result.returncode, 18, result.stderr)
        self.assertFalse((self.stage / (NF_FILENAME + ".complete")).exists())
        before = (self.stage / NF_FILENAME).stat().st_ino
        self.require_success(self.invoke())
        self.assertEqual(len(self.recorded_requests()), 2)  # One call per asset, never reacquired NF.
        self.assertEqual(before, (self.prefix / "bin" / NF_FILENAME).stat().st_ino)

    def test_corrupt_completed_download_is_retained_and_fails_again(self):
        result = self.invoke(extra_env={"FAKE_CURL_MODE": "corrupt"})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("checksum differs", result.stderr)
        bad = self.stage / NF_FILENAME
        before = bad.read_bytes(), bad.stat().st_mtime_ns
        result = self.invoke(extra_env={"FAKE_NETWORK_DENIED": "1"})
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(before, (bad.read_bytes(), bad.stat().st_mtime_ns))
        self.assertEqual(len(self.recorded_requests()), 1)
        self.assertFalse(self.prefix.exists())

    def test_unrelated_prefix_and_unowned_staging_directory_are_retained(self):
        self.prefix.mkdir()
        original = self.prefix / "unrelated.txt"
        original.write_text("Keep this unrelated directory.\n")
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(original.read_text(), "Keep this unrelated directory.\n")
        original.unlink()
        self.prefix.rmdir()
        self.stage.mkdir()
        original = self.stage / "unrelated.txt"
        original.write_text("Keep this unowned staging directory.\n")
        self.assertNotEqual(self.invoke().returncode, 0)
        self.assertEqual(original.read_text(), "Keep this unowned staging directory.\n")
        self.assertFalse(self.recorded_requests())

    def test_symlinked_prefix_is_rejected(self):
        untouched = self.root / "unrelated-target"
        untouched.mkdir()
        self.prefix.symlink_to(untouched, target_is_directory=True)
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("symlinked prefix", result.stderr)
        self.assertEqual(list(untouched.iterdir()), [])
        self.assertFalse(self.recorded_requests())

    def test_symlinked_staging_child_is_rejected_without_target_changes(self):
        self.assertEqual(self.invoke(extra_env={"FAKE_CURL_MODE": "interrupt-once"}).returncode, 18)
        target = self.root / "unrelated-target"
        target.write_bytes(b"Do not modify this unrelated file.\n")
        (self.stage / NF_FILENAME).unlink()
        (self.stage / NF_FILENAME).symlink_to(target)
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(target.read_bytes(), b"Do not modify this unrelated file.\n")
        self.assertEqual(len(self.recorded_requests()), 1)

    def test_wrong_java_version_is_not_published(self):
        self.write_payloads(java_version="21.0.12.2")
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Java version differs", result.stderr)
        self.assertFalse(self.prefix.exists())
        self.assertTrue(list(self.root.glob(self.prefix.name + ".partial.*")))

    def test_wrong_nextflow_version_is_not_published(self):
        self.write_payloads(nf_version="25.10.5")
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Nextflow version differs", result.stderr)
        self.assertFalse(self.prefix.exists())

    def test_changed_launcher_is_rejected_and_retained(self):
        self.require_success(self.invoke())
        launcher = self.prefix / "bin/nextflow"
        launcher.write_text(launcher.read_text() + "# Unrelated launcher edit.\n")
        before = launcher.read_bytes()
        result = self.invoke(extra_env={"FAKE_NETWORK_DENIED": "1"})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Launcher differs", result.stderr)
        self.assertEqual(launcher.read_bytes(), before)
        self.assertEqual(len(self.recorded_requests()), 2)

    def test_changed_extracted_jre_library_is_rejected(self):
        self.require_success(self.invoke())
        library = self.prefix / JAVA_DIRECTORY / "lib/fake-library.bin"
        library.write_bytes(b"Altered tiny fake JVM library.\n")
        result = self.invoke("check")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("JRE bytes differ", result.stderr)
        self.assertEqual(library.read_bytes(), b"Altered tiny fake JVM library.\n")

    def test_concurrent_installer_lock_fails_before_any_acquisition(self):
        lock = Path(str(self.prefix) + ".install.lock")
        with lock.open("a") as stream:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
            result = self.invoke()
        self.assertEqual(result.returncode, 75, result.stderr)
        self.assertFalse(self.prefix.exists())
        self.assertFalse(self.stage.exists())
        self.assertFalse(self.recorded_requests())

    def shorten_probe_timeout(self):
        # Check the real production wall/kill flags, but shorten only the tiny
        # fixture's timeout executable so a hung fake JVM does not cost90s.
        timeout = shutil.which("timeout")
        self.assertIsNotNone(timeout)
        wrapper = self.fake_bin / "timeout"
        wrapper.write_text("#!/usr/bin/env bash\nset -euo pipefail\n"
                           "[[ $1 == --kill-after=10s && $2 == 90s ]]\nshift 2\nexec "
                           + shlex.quote(timeout) + " --kill-after=0.2s 0.2s \"$@\"\n")
        wrapper.chmod(0o755)

    def test_hung_java_probe_times_out_without_publication(self):
        self.write_payloads(hang_java=True)
        self.shorten_probe_timeout()
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Java version check failed", result.stderr)
        self.assertFalse(self.prefix.exists())

    def test_hung_nextflow_probe_times_out_without_publication(self):
        self.write_payloads(hang_nextflow=True)
        self.shorten_probe_timeout()
        result = self.invoke()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Nextflow version check failed", result.stderr)
        self.assertFalse(self.prefix.exists())

    def test_normalized_root_prefix_is_rejected_before_installation(self):
        for prefix in ("/", "//", "/./"):
            with self.subTest(prefix=prefix):
                result = subprocess.run(["bash", str(self.helper), "install", "--prefix", prefix],
                                        cwd=self.root, env=self.environment, capture_output=True,
                                        text=True, timeout=10)
                self.assertNotEqual(result.returncode, 0)
                if prefix != "/":
                    self.assertIn("Refusing the filesystem root", result.stderr)
        self.assertFalse(self.recorded_requests())


if __name__ == "__main__":
    unittest.main()
