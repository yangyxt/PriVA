# PriVA tests

Installer diagnosis, 2026-10-15: all 28 checks in
`test_install_utils_bootstrap.py` pass, including the regression that a zero
VEP-installer exit without cache files must not record success. The test also
checks unattended, release-pinned cache/plugin flags. The actual Slurm rerun
reached the correct VEP113 GRCh37 archive download and stopped at a deliberate
timeout; that checkpoint does not establish complete installation.

`test_install_utils_resources.py` adds nine small command/fixture checks for
configured-reference propagation, all-chromosome scheduling, failure barriers,
Conservation ordering/reuse, and the configured AlphaMissense motif output.
External tools are replaced for these orchestration tests; they do not claim
scientific liftover or annotation equivalence.

`test_cluster_licensed_import.py` covers explicit resource subsets and the
unchanged combined default. Tiny synthetic BGZF/index fixtures run the frozen
installer functions and verify resource selection, receipt scope, unchanged
unselected configuration, and preserved source bytes. Preflight cases reject
incomplete/unconfirmed acquisitions, missing or extra roles, wrong assemblies,
wrong cache roots or variant types, partial filenames, missing/escaped indexes,
size/checksum mismatches, and existing destinations. Source/index mutations during checksum
stop before config creation or installer invocation; mutations during the
installer prevent a completion receipt. These tests require the installed
`priva_acmg` tools (`bgzip`, `tabix`, `yq`) and PyYAML; they never read the
MBP-managed downloads or perform real licensed-cache imports.

Cluster deployment integrity: `test_cluster_deployment.py` checks atomic
publication, corruption rejection, unchanged-cache reuse, and rejection of
unidentified or modified FASTA/container resources without network or solver
operations. Real installation qualification is recorded separately by Slurm.

Updated: 2026-09-14. Existing test modules cover the workflows and data contracts
identified by their module docstrings. This inventory records the startup tests
added alongside the installer change; other upstream tests are unchanged.

| File | Purpose and dependencies | Status |
| --- | --- | --- |
| test_install_utils_bootstrap.py | Standard-library tests for installer startup and cache paths. The existing InstallerBootstrapTests use temporary fake managers and require fixture approval. New CachePathTests exercise defaults, custom names, old-path relocation, validation, comment-preserving writes, repeat runs, check mode, startup failure propagation and CADD selection using only memory and mocked commands/filesystem calls. Neither class downloads resources. Requires Python 3, Bash and common shell utilities. | Updated 2026-09-14; 14 focused checks passed on macOS and RHEL 8.10; based on PriVA commit 810e60165a663bd2e1223809af9f32d664a9d543. |
| README.md | Test change inventory and execution instructions. | Added 2026-09-13. |

Run the startup suite from the repository root:

```bash
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -p test_install_utils_bootstrap.py -v
```

Run only the cache-path tests, without creating temporary fixtures:

```bash
PYTHONDONTWRITEBYTECODE=1 python tests/test_install_utils_bootstrap.py CachePathTests -v
```

The full discovery command above also runs InstallerBootstrapTests, which creates
temporary files/directories; obtain fixture approval before running it.
