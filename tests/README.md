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

`test_cluster_alpha.py` checks mapping preflight before annotation, source
identity rejection, the original checksum-pinned motif command's option order,
updated VCF propagation, timestamp reuse without overwrite, and failure/status
propagation. Mapping hash or metadata changes during work prevent completion.
These are command/fixture tests; annotation and motif computation are replaced,
and no mapping resource or scientific result is qualified by them.

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

The licensed-import fixtures do not require Git history: each current checkout
script is accepted only when its SHA256 equals the recorded bytes of installer
`054b4944055054148005016192647a8ed0384eef`. The four hashes are in
`REFERENCE_SHA256` in the test module. A changed source fails before creating the
fixture snapshot; there is no fallback to a newer installer. If these scripts
later change, explicitly retain/review the frozen oracle rather than replacing
its identity silently. CI installs Ubuntu's `tabix` package (including `bgzip`)
and `yq`, records their package versions, and retains real installer fixture
execution. This test tooling does not qualify production scientific software.

`test_cluster_primate_validation.py` compares the validator with the literal
original formatting pipeline, including duplicated/tied rows, header/blank/CRLF
handling and unterminated records. Dropped/replaced rows, changed ordering,
leading/trailing byte differences, partial-sort errors, source/derived gzip
errors and changed file identities must fail; owned children are reaped and
scratch removed on tested failure paths. The licensed-import suite also runs
an actual frozen formatter with an injected partial-output sort failure and
requires exact validation to prevent its otherwise false-success receipt.
Neither suite formats a production cache. PID-only SIGTERM cleanup is not
claimed by these tests.

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

## SpliceAI2 preparation checks

Current dispatch requirements: assembly chooses fresh GRCh37 v1 / GRCh38 v2;
decompressed headers choose compatible VCF passthrough versus native84 TSV
adaptation. Only spliceai_snv_prescore and spliceai_indel_prescore remain in
YAML. Sidecars and actual cache headers retain provenance/scientific gates.
The agreed helper API is:

```bash
python scripts/spliceai_cache_input.py install --config CONFIG --cache-dir CACHE_DIR --snv SNV --indel INDEL
```

SpliceAI_install invokes scripts/spliceai_cache_input.py install automatically, retaining the existing prompts and two resource keys. The helper detects gzip by magic bytes, then classifies a bounded decompressed header as vcf_v1, vcf_spliceai2 or the exact native 84-column SpliceAI2 TSV schema; filenames, extensions and score magnitudes do not select the predictor. An unmarked VCF receives the legacy-compatible vcf_v1 routing label; that label alone does not prove predictor identity. Existing compatible VCF pairs pass through, including existing GRCh38 v1; fresh GRCh38 selects v2. Wrong assembly, mixed types or failing v2 preparation fail without fallback. SpliceAI2 VCF passthrough additionally requires a completed adjacent adapter.ready.json with exact SNV/indel/index membership and complete published-catalogue status; it performs no conversion. Reject .partial or .building path components, including symlink targets, and refuse in-progress preparation or staged outputs. A final VCF renamed before receipt publication is not ready. V1 VCF passthrough requires no additional adapter receipt.

For native TSV, the helper discovers the fixed adjacent SpliceAI2_GRCh38.downloads-complete.json sidecar. The selected SNV/indel inputs must match their manifest roles and the full 47-shard catalogue. It calls the existing low-level SpliceAI2_install to prepare or reuse a fully validated receipt, then updates only spliceai_snv_prescore and spliceai_indel_prescore. The pinned priva converter is resolved through PATH; no YAML converter or manifest option is added. Provenance and scientific qualification remain in sidecars/cache headers. The fixed adjacent sidecar is published at yangyxt/caches/SpliceAI2_GRCh38.downloads-complete.json, SHA256 fbc6c9fcd1859d6f1605f9343d5659b8420d7397c4152da19daab34c6194e49d. The PATH entry yangyxt/softwares/bin/priva resolves to the verified optimized 225aa4a converter. This changes no YAML keys beyond the original two paths.

Two-key interface qualification: code 35779c4578509ce01218a052a8e7e1d8a785a75b passed 239 Python tests, zero failures/errors/skips, in 190.848 seconds: 39 helper, 25 installer, 14 runtime-guard, 21 Nextflow-adapter and 140 existing deployment/driver/BGZ tests. The installer-shell cases mock the helper boundary; they are not real conversion evidence. Log: yangyxt_dt/projects/PriVA_upgrade/logs/spliceai-two-key-integration-20261016.log. This supersedes the earlier 197-test interface checkpoint without adding its count; unchanged Rust evidence remains 41 tests and the earlier real-projection evidence remains 512 unique rows.

Actual resource checks used final helper SHA256 f69d6ea8d6179d68ae44cad849b3c6c7912a879f05e305a36acffba09c81bb88: all 47 native TSV and four v1 VCF headers classified correctly, including assembly metadata. Real SpliceAI_install reused the existing v1 pair for hg19/GRCh37 and hg38/GRCh38 with configuration bytes, inode and mtime unchanged and resource stat identities unchanged. Receipt: yangyxt_dt/projects/PriVA_upgrade/validation/spliceai-two-key-real-install-20261016/receipt.json. These are bounded header/reuse checks, not new body checksums, full VEP execution or model/scientific qualification.

Nextflow stub verification passed all four scenarios: raw input, resume with all four tasks cached, preannotated input and stale-environment rejection. Log: yangyxt_dt/projects/PriVA_upgrade/logs/spliceai-two-key-nextflow-20261016.log; retained validation/log artifacts: validation/spliceai-two-key-nextflow-20261016/ under that data project. The run used code 35779c4 plus the documentation state at execution; its recorded source_digest identifies that exact state. Subsequent documentation edits are not attributed to the run.

Frozen-run boundary: job 1514732 continues unchanged with installer f8de4b3 and the previously verified converter. Its prior explicit-manifest/six-key interface is historical execution provenance, not the new installation contract. After successful completion, validated output artifacts may be imported through the two original keys; do not attribute those artifacts or the 197 Python tests to the new dispatcher. The new two-key interface is qualified separately at 35779c4; full scientific acceptance remains pending.

The older suite results below retain their frozen source identities; the current
239-test qualification above covers the completed two-key interface:

`test_spliceai2_deployment.py` exercises manifest completeness/identity checks,
GENCODE context agreement, source preservation, conversion/index/report failures
and preparation receipts with synthetic fixtures. `test_install_utils_spliceai2.py`
checks assembly/profile routing, unchanged v1 behavior, full-receipt registration
and refusal to activate the new scientific profile. These are preparation tests;
they do not convert the licensed 47-shard catalogue or qualify v2 ACMG evidence.

`rust/tests/spliceai2_projection.rs` covers the native 84-column contract, exact
top0 string projection, POS/model-strand gene joins, duplicate/ambiguity rejection
and CLI behavior. The revised policy retains all context IDs/intervals and
coalesces only same-strand symbol labels within each source prediction; opposite-
strand symbol conflicts and duplicate native predictions still fail. Final integrated test totals belong to their producing source
revision; the tests do not establish full VEP lookup coverage.

`spliceai2/test_official_plugin.pl` passed 11 subtests / 38 assertions against the
actual pinned release-113 plugin and Ensembl allele matcher. It also rejects a
wrong plugin identity before loading. Only tabix candidates/file registration and
light VEP objects are fixtures; real reference/cache/index retrieval is not
exercised. See [the audit](../docs/spliceai2-plugin-audit.md) for hashes and runtime
instructions. Synthetic records contain no licensed cache rows.

```bash
PRIVA_SPLICEAI2_BIN=/absolute/built/priva \
  python -m unittest discover -s tests -p test_spliceai2_deployment.py -v
python -m unittest discover -s tests -p test_install_utils_spliceai2.py -v
python -m unittest discover -s tests -p test_spliceai2_runtime_guard.py -v
cargo test --manifest-path rust/Cargo.toml --locked --test spliceai2_projection
# Requires installed pinned official VEP/plugin modules and priva_acmg Perl:
perl tests/spliceai2/test_official_plugin.pl \
  --plugin-file /absolute/VEP_plugin/SpliceAI.pm \
  --vep-lib /absolute/priva_acmg/share/ensembl-vep-113.2-0
```

The driver integration cases skip when `PRIVA_SPLICEAI2_BIN` is absent; set it
to the reviewed real converter to reproduce the no-skip result. Use a Slurm
allocation on this cluster. The [preparation reference](../docs/nextflow_migration.md#spliceai2-preparation)
distinguishes complete published coverage from bounded subset fixtures and
scientific activation. No GPU/model inference or new resource downloads are
required for the synthetic checks.

Historical combined run before the stable-ID/same-strand-symbol policy update:
**165 tests passed, zero failures/skips, 85.702 seconds**
(119 existing/route + 14 runtime guards + 21 driver + 11 deployment). The driver
suite used real Rust/bgzip/tabix with synthetic scores/context/plugin identity,
including a complete 47-shard fixture and CRC/mutation failures. Its retained log
`yangyxt_dt/projects/PriVA_upgrade/logs/spliceai2-install-integration-20261015.log`
records integration `ad9af5d` and tested driver/test content hashes. Separately,
Rust source `428ce4c` passed 30 projection + 9 existing tests, and the pinned
Perl harness passed 11 subtests / 38 assertions. These totals have different
boundaries and do not establish real resource conversion or scientific parity.
Current Rust 225aa4a (integration 32212d7) passed 32 projection + 9 existing tests.
The final Python suite at f8de4b32beaf3bf5185ef53f75be162833055153 passed 197 tests,
zero failures/skips, in 190.847 seconds: 169 installer/deployment (including 25 driver)
+ 20 Nextflow adapter + 8 BGZ-validator. Log:
`yangyxt_dt/projects/PriVA_upgrade/logs/spliceai2-final-integration-20261015.log`.

Bounded real-cache check: 256 SNV plus 256 indel prefix rows from chr1, all positive-strand, passed an independent column-name-based full-record/order oracle and real bgzip/tabix queries with zero mismatches. Verifier f8de4b3 first checked the debug converter from 225aa4a; receipt: yangyxt_dt/projects/PriVA_upgrade/validation/spliceai2-real-projection-225aa4a-debug/receipt.json. The same 512 source rows were independently rechecked against the installed optimized release with zero mismatches; this remains 512 unique rows, not 1,024. This is real precomputed-row conversion, not model inference, whole-cache coverage, full VEP execution, public-sample or supplied-family acceptance.

Installed release: source 225aa4aa774b4be7f53576e1cde56af71e67bd2c at
`yangyxt/softwares/priva-tools-225aa4a/bin/priva`, SHA256
`bbc3b7288ed2ea63f1f8455ed87b8a9d08c7bc1fe0a562c1d6134525e2a12dc6`.
The same 25 driver cases passed again against this release. Unique suite totals
remain tied to their checkpoint: 197 Python then, 239 Python now, with unchanged
41 Rust tests. A stale shared Cargo artifact was detected
before installation; forced rebuilding and an installed smoke check corrected
it, while the receipt retains the failed attempt. The complete execution manifest/handoff were saved in DataLad commit 5f5946b
before launch. Full preparation job 1514732 is running on c81o-03 with frozen installer f8de4b3 and the
installed release (4 CPUs, 16 GB, 7 days); these reruns and the submission do not
establish full conversion or VEP qualification. Separately, the frozen job has
passed all 47 full source SHA256 checks and started projecting SNV chr1; no
complete conversion receipt exists yet.

Separately, the
complete real metadata check passed 20,083 unique stable-ID/context mappings,
zero ID/bounds/strand/biotype mismatches and ten same-strand symbol groups; all
versioned GTF IDs remain retained. Closed evidence is under
`yangyxt_dt/projects/PriVA_upgrade/validation/spliceai2-real-context-20261015/`.
This checks gene context, not score conversion or indexed-VEP execution.

Scientific boundary for the simplified installer: inspect the actual selected cache headers and sidecar provenance to prevent unqualified v2 scientific use, including preannotated import and resume; no new enable/disable YAML toggle is introduced. Successful installation is separate from v2 ACMG/ClinVar, full VEP, public-sample and supplied-family acceptance. Unmarked external annotations still require pinned scientific provenance because predictor version cannot be inferred from numeric values. The runtime guards now inspect actual cache headers at the two original keys, including preannotated launch/task paths; the current 21 Nextflow-adapter and 14 runtime-guard tests are included once in the 239-test total. Earlier f8de4b3 checks remain historical; scientific acceptance is still pending.
