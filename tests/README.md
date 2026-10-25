# PriVA tests

Task read-only profile,24October2026 CT:35 tests in `tests/nextflow/test_*.py`
passed in closed1524526, followed by actual Bubblewrap write-denial/exit-code
checks and default-versus-boundary Nextflow stub output equality. Raw, all-four
cached resume, annotated and changed-profile rejection cases passed. These are
synthetic engineering fixtures, not model/family/public or whole-installation
acceptance. Failed1524525 and its effective-config diagnosis remain preserved.
Sources are `bin/priva_task_boundary.py`, `conf/task_readonly.config` and the
run/collect controller in the source-side PriVA_upgrade/scripts directory;
runtime/logs are in its equivalent yangyxt_dt `validation/nextflow-readonly-check-*`
and `logs/nextflow-readonly-check-*` paths.

RefSeq acquisition, 2026-10-19: `test_cadd_refseq_resources.py` passes 23 focused
synthetic checks per checkout with mocked HTTP and temporary files. Run
`python3 -B -m unittest discover -s tests -p test_cadd_refseq_resources.py -v`.
It checks release/publisher/identity/resume/publication boundaries and local
read-only verification. Closed evidence is
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/deployment/cadd-refseq-candidate-helper-tests-20261019-v1.log`;
source/tests are in
`/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/{PriVA,PriVA-nextflow}/{scripts,tests}`.
These checks establish acquisition behavior, not real source coverage or scoring.

Current installer recovery checkpoint, 2026-10-18: InterPro consumers honor the
selected raw/gzip mapping; unexpected chromosome failures abort publication.
HPO/MONDO and DAS multi-output/intermediate caches now have checked recovery.
DAS CLI and PM1 regeneration require explicit completion after the original
writes; early exits cannot register stale results. Twelve labelled boundary
cases and 24 existing reader tests passed per checkout; full-model checks remain separate.
DAS metadata binds to the selected VCF, assembly and extraction source; only
managed derived metadata regenerates when that association is missing/stale.
Ordinary valid pickle values/order and bundled existence/nonzero policy remain.
Config-file arguments are quoted through installer dispatch. The root README
has exact focused-test receipts and separate pending full/scientific acceptance.
Source/tests live here; all logs and synthetic runtime fixtures are under
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/{logs,deployment,validation,tmp}`;
maintained harnesses are under
`/stornext/snfs170/next-gen/scratch/yangyxt/projects/PriVA_upgrade/scripts/`.


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

The licensed-import fixtures do not require Git history:
`fixtures/installer-054b494/scripts/*` retains the exact oracle, and each script
must match the recorded bytes of installer
`054b4944055054148005016192647a8ed0384eef`. A current-checkout fallback is accepted
only when its SHA256 still matches that historical oracle. The four hashes are in
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

Cluster deployment integrity: `test_cluster_deployment.py` checks the pinned
Washington GRCh37 v1.7 data/index pair and matching checksum-source host/release,
atomic publication, full-MD5 reuse of a completed partial without HTTP,
incomplete/corrupt partials, size gates and full rechecking of existing finals.
It also rejects unidentified or modified FASTA/container resources without
network or solver operations. At source `20774bd56a706d8b82a5640bf23b758216544136`,
all **20 deployment + 31 licensed-import tests passed** in Slurm step
1514363.254. The initial run retained 21 licensed-fixture failures caused by a
missing historical helper; the exact `054b494` helper was restored with SHA256
`c5e345f6f3af31aff1e2a724e248f0f32a1085afe7b347183022196b359ea6cf`,
without changing expected hashes. Both closed logs follow
`yangyxt_dt/projects/PriVA_upgrade/logs/cadd-prescore-recovery-tests-20261016{,-final}.log`
and include before/after tested source hashes. These are installation fixtures;
real resource, VEP, model and family qualification are recorded separately.

Real recovery 1514744 subsequently COMPLETED (0:0) at
2026-10-16T08:02:48Z with all four full MD5 checks and a deployed receipt.
It reused the existing SNV/index and indel bytes; only the matching indel index
was downloaded. Bounded real indexed-pair validation passed in Slurm step
1514363.262 at 08:03:30Z: 24 contigs per pair, matching 135-column schemas,
every contig's first row/exact-coordinate refetch, stable file stat identities,
expected receipt/pins and no remaining partials. This audit did not repeat
full-body MD5/CRC; the completed recovery supplies full MD5 evidence.
DataLad `9e70725` records the closed receipt/log/report evidence. These checks
remain separate from the 51 fixture tests and from real VEP/model/family qualification.

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

SpliceAI_install invokes scripts/spliceai_cache_input.py install automatically, retaining the existing prompts and two resource keys. The helper detects gzip by magic bytes, then classifies a bounded decompressed header as vcf_v1, vcf_spliceai2 or the exact native 84-column SpliceAI2 TSV schema; filenames, extensions and score magnitudes do not select the predictor. An unmarked VCF receives the legacy-compatible vcf_v1 routing label; that label alone does not prove predictor identity. Existing compatible VCF pairs pass through, including existing GRCh38 v1; fresh GRCh38 selects v2. Wrong assembly, mixed types or failing v2 preparation fail without fallback. SpliceAI2 VCF passthrough additionally requires a completed adjacent adapter.ready.json with exact SNV/indel/index membership and complete published-catalogue status. This checks completed-publication metadata and recorded membership, not a new full VCF/index-body checksum; it performs no conversion. Reject .partial or .building path components, including symlink targets, and refuse in-progress preparation or staged outputs. A final VCF renamed before receipt publication is not ready. V1 VCF passthrough requires no additional adapter receipt.

For native TSV, the helper discovers the fixed adjacent SpliceAI2_GRCh38.downloads-complete.json sidecar. The selected SNV/indel inputs must match their manifest roles and the full 47-shard catalogue. The existing low-level SpliceAI2_install preparation/reuse path performs full checksum verification and prepares or reuses a validated receipt, then updates only spliceai_snv_prescore and spliceai_indel_prescore. Fresh conversion builds/caches the bundled standalone installation converter; no YAML converter or manifest option is added. Provenance and scientific qualification remain in sidecars/cache headers. The current adjacent sidecar at yangyxt/caches/SpliceAI2_GRCh38.downloads-complete.json has SHA256 bf546817c03689359138d9230e3587fe225764a06546b3914c8819afbbd2bdde: all 47 original source identities remain unchanged; the qualified plugin pin is116. The historical manifest remains preserved. The historical complete preparation used yangyxt/softwares/bin/priva selecting priva-tools-spliceai2-mnv-20261016/bin/priva, SHA25670eefad75c6ad5a8001fcd87ca6970c3b0b1495f21b55f8c7775db67122c567a; that production executable is retained as a comparison reference. This changes no YAML keys beyond the original two paths.

Two-key interface qualification: code 35779c4578509ce01218a052a8e7e1d8a785a75b passed 239 Python tests, zero failures/errors/skips, in 190.848 seconds: 39 helper, 25 installer, 14 runtime-guard, 21 Nextflow-adapter and 140 existing deployment/driver/BGZ tests. The installer-shell cases mock the helper boundary; they are not real conversion evidence. Log: yangyxt_dt/projects/PriVA_upgrade/logs/spliceai-two-key-integration-20261016.log. This supersedes the earlier 197-test interface checkpoint without adding its count; unchanged Rust evidence remains 41 tests and the earlier real-projection evidence remains 512 unique rows.

Actual resource checks used final helper SHA256 f69d6ea8d6179d68ae44cad849b3c6c7912a879f05e305a36acffba09c81bb88: all 47 native TSV and four v1 VCF headers classified correctly, including assembly metadata. Real SpliceAI_install reused the existing v1 pair for hg19/GRCh37 and hg38/GRCh38 with configuration bytes, inode and mtime unchanged and resource stat identities unchanged. Receipt: yangyxt_dt/projects/PriVA_upgrade/validation/spliceai-two-key-real-install-20261016/receipt.json. These are bounded header/reuse checks, not new body checksums, full VEP execution or model/scientific qualification.

Nextflow stub verification passed all four scenarios: raw input, resume with all four tasks cached, preannotated input and stale-environment rejection. Log: yangyxt_dt/projects/PriVA_upgrade/logs/spliceai-two-key-nextflow-20261016.log; retained validation/log artifacts: validation/spliceai-two-key-nextflow-20261016/ under that data project. The run used code 35779c4 plus the documentation state at execution; its recorded source_digest identifies that exact state. Subsequent documentation edits are not attributed to the run.

Historical frozen-run boundary: job 1514732 used installer f8de4b3 and converter225aa4a, then failed at the publisher MNV on indel chr1 line11929. It produced no completed adapter receipt. Its explicit-manifest/six-key interface and197 Python tests remain historical provenance. Recovery retains all47 original download identities and the two-key public installer interface; only validated complete outputs can later be imported. Full scientific acceptance remains pending.

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
PRIVA_SPLICEAI2_INSTALL_BIN=/absolute/built/spliceai2-converter \
  python -m unittest discover -s tests -p test_spliceai2_deployment.py -v
python -m unittest discover -s tests -p test_install_utils_spliceai2.py -v
python -m unittest discover -s tests -p test_spliceai2_runtime_guard.py -v
cargo test --manifest-path rust/Cargo.toml --locked --test spliceai2_projection
# Requires installed pinned official VEP/plugin modules and priva_acmg Perl:
perl tests/spliceai2/test_official_plugin.pl \
  --plugin-file /absolute/VEP_plugin/SpliceAI.pm \
  --vep-lib /absolute/priva_acmg/share/ensembl-vep-113.2-0
```

The driver integration cases skip when `PRIVA_SPLICEAI2_INSTALL_BIN` is absent; set it
to the reviewed standalone installation converter to reproduce the no-skip result. Use a Slurm
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
Historical Rust 225aa4a (integration 32212d7) passed 32 projection + 9 existing tests.
The final Python suite at f8de4b32beaf3bf5185ef53f75be162833055153 passed 197 tests,
zero failures/skips, in 190.847 seconds: 169 installer/deployment (including 25 driver)
+ 20 Nextflow adapter + 8 BGZ-validator. Log:
`yangyxt_dt/projects/PriVA_upgrade/logs/spliceai2-final-integration-20261015.log`.

Bounded real-cache check: 256 SNV plus 256 indel prefix rows from chr1, all positive-strand, passed an independent column-name-based full-record/order oracle and real bgzip/tabix queries with zero mismatches. Verifier f8de4b3 first checked the debug converter from 225aa4a; receipt: yangyxt_dt/projects/PriVA_upgrade/validation/spliceai2-real-projection-225aa4a-debug/receipt.json. The same 512 source rows were independently rechecked against the installed optimized release with zero mismatches; this remains 512 unique rows, not 1,024. This is real precomputed-row conversion, not model inference, whole-cache coverage, full VEP execution, public-sample or supplied-family acceptance.

Historical installed release: source 225aa4aa774b4be7f53576e1cde56af71e67bd2c at
`yangyxt/softwares/priva-tools-225aa4a/bin/priva`, SHA256
`bbc3b7288ed2ea63f1f8455ed87b8a9d08c7bc1fe0a562c1d6134525e2a12dc6`.
The same 25 driver cases passed again against this release. Unique suite totals
remain tied to their checkpoint: 197 Python then, 239 Python now, with unchanged
41 Rust tests. A stale shared Cargo artifact was detected
before installation; forced rebuilding and an installed smoke check corrected
it, while the receipt retains the failed attempt. The complete execution manifest/handoff were saved in DataLad commit 5f5946b
before launch. Historical preparation job 1514732 ran on c81o-03 and failed on the MNV with frozen installer f8de4b3 and the
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

Scientific boundary for the simplified installer: inspect the actual selected cache headers and sidecar provenance to prevent unqualified v2 scientific use, including preannotated import and resume; no new enable/disable YAML toggle is introduced. Successful installation is separate from v2 ACMG/ClinVar, full VEP, public-sample and supplied-family acceptance. Unmarked external annotations still require pinned scientific provenance because predictor version cannot be inferred from numeric values. The runtime guards now inspect actual cache headers at the two original keys, including preannotated launch/task paths; the earlier `35779c4` checkpoint included 21 Nextflow-adapter and 14 runtime-guard tests once in its 239-test total; the focused PED checks below are a later, separate test run. Earlier f8de4b3 checks remain historical; scientific acceptance is still pending.

PED compatibility, 2026-10-16: the Nextflow adapter accepts phenotype `0/-9`
as missing affection and preserves sex, parent IDs, PED bytes and family row
order. Added generic adapter regressions cover both missing parental states,
full-cohort sample membership, invalid codes and the affected-first guard.
Scientific filtering and sex-dependent interpretation remain unchanged; these
contract checks do not establish scientific family acceptance.

Installer repair coverage, 2026-10-16: `test_install_utils_*.py`,
`test_installer_{prerequisites,boundaries}.py`, `test_spliceai_cache_input.py`,
`test_spliceai2_runtime_guard.py` and `test_vcf_bgz_validation.py` cover shell
routing/failure propagation, manual prerequisites, CADD checksum/resume/extraction,
and SpliceAI header/index/registration gates. Tests use temporary synthetic inputs,
fake external tools where marked, local HTTP recovery fixtures, and small real
BGZF/tabix operations. Startup-only tests stub the separately tested prerequisite
stage. `fixtures/installer-054b494/scripts/{install_utils,annotation_vcf,common_bash_utils}.sh`
are unchanged historical oracles, copied from the existing Nextflow test suite;
they are never production entry points. No real-model, full-cache, family or
public-sample scientific acceptance is established by these tests.

Installer-repair checkpoint before BioMart integration: 210 distinct focused tests passed for this checkout in Slurm allocation 1514363. Bash syntax and whitespace checks passed; the actual example config preflight failed as intended on unresolved original host paths. Closed logs and before-edit backups are retained under `yangyxt_dt/projects/PriVA_upgrade/validation/installer-repair-20261016/`. Earlier failed fixture attempts remain recorded, followed by their successful corrected reruns.

BioMart integration verification, 2026-10-16: `test_biomart_nm_to_enst.py` and
`test_installer_prerequisites.py` now cover completed export/schema, exact snapshot
pinning, row multiplicity/order, missing or stale receipts, interrupted publication,
concurrent writers, explicit manual inputs and assembly gating. Both repositories
passed 231 focused installer regression tests in Slurm 1514363. The subsequent
snapshot-lock change passed all 41 affected tests per repository (15 helper plus
26 prerequisites, including three new tests); these overlap the earlier suite
and must not be summed as distinct tests. Acquisition/error-path tests use small
synthetic or mocked responses, not a repeat public download.

Separately, the complete real BioMart export from job 1514785 passed schema,
completion and checksum validation; both unchanged scientific readers passed
independent multiplicity/order and HCSeeker coverage checks. Final real-bundle
installer checks passed two offline reuses per checkout, preserved bundle/config
identities and rejected implicit GRCh38 fallback without modifying configuration.
Bash syntax, whitespace and unchanged scientific/other-script hashes passed.
Closed evidence is in `yangyxt_dt/projects/PriVA_upgrade/validation/biomart-mapping-20261016/`;
`*-installer-tests.log` records the broad checkpoint, `*-snapshot-lock-tests-attempt1.log`
records the final affected suite, and `final-installation-verification.{json,log}`
records the final source identities and real offline checks. Coverage is partial:
1,469/1,549 hg19 HCSeeker rows map. See the root README's BioMart record for the
80 omitted rows' downstream impact, source pin and full size/SHA256/provenance.
No original-table equivalence, full AlphaMissense calibration, public-sample or
family scientific acceptance is claimed. Source changes remain uncommitted.

Reference-route verification, 2026-10-16: the final corrected-state suite passed
**82 tests per checkout** in Slurm1514363, including real samtools tiny-FASTA
indexing, mocked HTTP resume/corruption, assembly/provenance failures, quoted
atomic config updates, and a validator-backed GoFCards reuse case proving no
opposite-assembly acquisition. These are fixtures, not production FASTA downloads.
Live publisher manifests/ranged FASTA responses and the full small chain checksum
were checked separately. Final logs are
`yangyxt_dt/projects/PriVA_upgrade/logs/reference_reuse_PriVA{,-nextflow}_20261016.log`.
Earlier `reference_{resources,integration}_*.log` attempts remain retained, including
the corrected error-message failure; their counts must not be added to the final
82. Scientific Python remained unchanged; full scientific acceptance is pending.

CADD runtime verification, 2026-10-16: 15 focused runtime tests plus 20 existing
deployment tests passed (35 total) in Slurm1514363. Tests cover wrapper conflicts,
model/resource readiness, image-source pin checks, offline environment probes,
reference-only VEP110 help-banner validation and receipt publication barriers.
Real recheck1514901 passed separately; first attempt1514899 failed on unsupported
VEP `--version` before publication. `cadd-runtime-complete-20261016*-tests.log` and
`cadd-runtime-vep-option-probe.json` under the data project's `validation/` retain
fixture/probe evidence. No scoring/model inference, baseline profiling rerun,
shared-VEP, public-sample or family scientific acceptance is implied.

VEP release-selection verification, 2026-10-16:

Validation on HGSC, 2026-10-16: **165 focused fixture tests passed per checkout**
in Slurm1514363.346 (330 executions; the earlier 51-test checkpoint overlaps).
This covers requested versions, manager failures/conflicts, API/helper mismatches,
plugin receipts/custom edits/interrupted upgrades, cache/index failures and affected
installer routes. Pytest was unavailable: six unittest files and all 17 GoFCards
functions ran through a retained in-memory unittest adapter supplying only `tmp_path`;
no test dependency was installed. The earlier missing-pytest attempt is retained.
Read-only real VEP113 executable/API/installer checks and the existing GRCh37 cache
checks passed separately, including the final CSI-only cache guard. No real VEP116
installation/download, model run, public-sample or family acceptance was performed.
Closed commands, harness and producer hashes are retained under
`yangyxt_dt/projects/PriVA_upgrade/logs/vep-release-*.log`.

Real VEP116 deployment follow-up, 2026-10-16: Slurm1514974 passed56 focused
bootstrap/cache/plugin/failure tests per checkout after removing the redundant
`--freeze-installed` flag; exact other-package pins remain. These cases overlap the
prior165-case suite. Actual solver probes1514972 and the completed116.3upgrade1514975
separately proved a VEP-only transaction with527non-VEP identities unchanged;
software import/runtime checks passed. These checks do not qualify scientific
annotation, models or family results. The Nextflow checkout additionally passed
12new release/profile routing cases and20existing deployment tests, including
rejection of a113software receipt for a116config. Closed evidence uses
`yangyxt_dt/projects/PriVA_upgrade/logs/vep116-*`.

SpliceAI2 MNV recovery, 2026-10-16: `spliceai2_to_vep_vcf.py` now explicitly
qualifies official113/116 plugin identities; `test_spliceai2_plugin_identity.py`
passes 5 cases per checkout. The shared optimized Rust converter accepts intact
multi-base substitutions in publisher indel shards. Full deployment/test limits
and recovery job1515009 are recorded in the repository README. No scientific
Python, YAML resource keys or annotation sequence was changed by this correction.

ClinVar regeneration verification, 2026-10-17: `test_clinvar_regeneration.py` adds
seven native unittest groups for explicit exact-version transcript mapping,
missingness, strict donor-loss threshold, category precedence/review filtering,
retained NMD label behavior, the strict CLNCSQ default and staged annotation
contracts. Slurm 1515126 passed all seven per checkout, 1,024 exact synthetic
gene-group comparisons with frozen scientific source `09223ef`, and two exact
8-row/43-column reader comparisons. Seven reproduced frozen-source missingness
exceptions are explicitly repaired to absent evidence; no other mismatches were
observed. This verifies reader behavior, not full new-cache or family acceptance.
Closed regression evidence is in
`yangyxt_dt/projects/PriVA_upgrade/deployment/clinvar_candidate_20261011_vep116/regressions.ready.json`.

Unified CADD installer verification, 2026-10-17: **53 focused tests passed in the
original checkout and 88 in the Nextflow checkout**, in Slurm1514363. These include
17 new stage/receipt/pin/active-writer fixtures, 36 cache integrity/resume cases
(including two no-duplicate-archive guards), and the Nextflow checkout's 15 runtime
plus 20 deployment regressions. Exact final logs use
`yangyxt_dt/projects/PriVA_upgrade/logs/cadd-runtime-installer-tests-{primary,nextflow}-final-20261017.log`.
The separate real calculation-only GRCh37 and GRCh38 checks verified the existing
image, Python2 route/imports, original source/model/PHRED assets and standalone
873e9c1e executable; receipts use `deployment/cadd-install.GRCh*.calculate.*.ready.json`.
No scores were calculated. Complete annotation verification remains queued in
1515129 after 1515113. Installation evidence does not qualify shared VEP features,
public data or family results; the NaN CLI and experimental cancellation gates
remain separate.

2026-10-17: `test_install_utils_cadd_runtime.py` adds7 shell-only cases for stage/config forwarding and main-install success/failure propagation; no downloads or model inference. Root verification log is `yangyxt_dt/projects/PriVA_upgrade/logs/installer-root-checks-20261017.log`, saved in DataLad243dca6.

ClinVar installer follow-up, 2026-10-17: Slurm 1515135 passed nine reader/contract
groups and ten stage/failure-publication groups per checkout. The corrected
native VEP116 pilot 1515131 produced 1,508 AA changes and 3,790 splice entries.
Full-reference join audit 1515132 found zero differences across 4,615,346
GRCh38 AlleleID mappings, with 4,554,749 exact-version MANE targets. Coverage
remains incomplete: the 4,594,726-record raw VCF has 278 AlleleIDs absent from
the GRCh38 summary and 40,067 summary matches without an exact MANE target;
64 of 369 high-review pathogenic pilot targets are absent from VEP. The original
exclusions are retained, so these rows can omit evidence in new gene statistics.
No transcript substitution or complete scientific qualification is established.
The isolated full build is queued in original order; completed status requires
its final receipt, separately from these tests and pilot checks.

CADD source-acquisition follow-up, 2026-10-17: the focused suites now pass58
primary and93 Nextflow cases, including20 runtime-install cases and38 cache/shell
cases. The original configured-parent bootstrap fixture now asserts the exact
pinned source-only helper/config/destination, one bounded parent search, no
prompt/environment search and failure propagation. That precise fixture passes
in both checkouts with unchanged before/after source hashes; logs use
`cadd-runtime-installer-source-tests-{primary,nextflow}-20261017.log` and
`cadd-runtime-source-bootstrap-fixture-correction-20261017.log` under the existing
data-side project log directory.

Installer shell follow-up, 2026-10-17: all **56 focused tests per checkout**
(37 prerequisite and 19 resource groups) passed in existing Slurm allocation
1514363, step424. New behavior checks compare complete alias/canonical motif,
PM1 and domain-filter calls and destinations, exercise LoFtee reuse with native
`bash -u` for all four assembly labels, and check the explicit ClinVar candidate
diagnostic. Bash syntax and whitespace checks passed. External scientific
builders are stubbed in these routing tests; existing mapping/receipt checks run
against real small fixtures and bundled anchors. Closed merged evidence is
`yangyxt_dt/projects/PriVA_upgrade/logs/clinvar_installer_edgecases_1514363_20261017.log`.

Conservation integrity follow-up, 2026-10-17: **47 tests per checkout** passed in
Slurm allocation1514363, step468: 16 actual small-BigWig helper tests, 12 shell
barrier tests with a stubbed network writer, and 19 existing resource groups.
Coverage includes corrupt/truncated/wrong-assembly inputs, shared reference
lengths, extra contigs, exact default size/BSD pins, preserved custom identity,
publication races/failures and YAML failure barriers. Source-routing stubs remain
explicitly separate from the native reader tests. Bash syntax/whitespace passed.
Fixtures used `PRIVA_TEST_TMPDIR`/`TMPDIR` under the authorized data-side project
`tmp/`; closed logs are `yangyxt_dt/projects/PriVA_upgrade/logs/conservation_*_1514363_20261017.log`.
No production download, profile activation or scientific acceptance was tested.

Installer edge-failure checks, 2026-10-17: `test_install_utils_loftee.py` uses actual small gzip/FAI, SQLite and BigWig fixtures to cover reuse, acquisition/validation/configuration failures, hg19 nounset handling and caller-directory preservation. `test_install_utils_clinvar.py` additionally covers canonical cache registration and untouched explicit missing/incomplete/corrupt custom inputs. `test_install_utils_download_failures.py` checks InterPro exact-source retention, reuse, rebuilding and identity failures. Existing LoFtee routing stubs remain distinct from native format checks. Run generated fixtures/logs under the authorized data-side project root via `PRIVA_TEST_TMPDIR`/`TMPDIR`; these tests perform no network transfers or scientific qualification.

Configured HPO path checks, 2026-10-17: `test_hpo_runtime_paths.py` covers explicit/default routing, missing inputs, quoted paths, freshness and main-call order; `test_installer_prerequisites.py` checks data-side builder targets and existing path guards. Counterpart tests: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/tests/test_hpo_runtime_paths.py and /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/tests/test_installer_prerequisites.py. The Nextflow manifest identity/default/safety checks are in /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests/test_hpo_manifest_paths.py. Conversion stubs test shell routing separately; the scientific pandas implementation is unchanged. All new fixtures/logs use the authorized data-side project root.

Explicit output-route checks, 2026-10-17: `test_install_output_paths.py` separates stubbed scientific parsers from actual gzip/SHA256, path, timestamp and native awk behavior. It covers exact data-side paths, independent ClinGen destinations, retained InterPro source/cache reuse, failure barriers, legacy compatibility and unchanged curated-source payloads. The manual InterPro report is checked separately from cache metadata. Counterpart tests: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/tests/test_install_output_paths.py. Together with existing download-failure tests, 33 checks passed per checkout under Slurm1514363; all generated fixtures/logs use the authorized data-side project root.

Managed-reference cache checks, 2026-10-17: `test_reference_resources.py` now covers installer opt-in versus direct `--data-dir`, all four assembly aliases, new managed defaults, exact legacy-default redirection, preserved original/old-managed identities, missing custom inputs, orphan/staging retention, chain direction and matching preflight/install plans. Tiny real FASTAs use native samtools indexing; HTTP bytes are explicit fixtures. Existing merged-cache metadata fixtures now include the required official-format species/assembly/tabix fields and index. All 42 focused checks passed in each checkout under Slurm1514363, with new fixtures/logs under the authorized data-side project. Counterpart tests: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/tests/test_reference_resources.py. No live reference/index, production configuration or frozen job was changed.

Final installer routing/I/O verification, 2026-10-17: `test_install_utils_{final_routes,cds_io}.py` and `test_selected_clinvar_vcv.py` cover UTR cwd, exact CDS/XML paths, missing/invalid inputs, retry/configuration failures, GRCh38 metadata aliasing, custom XML content identity and immutable injection/full-build gates. Small genuine gzip/plain FASTAs exercise the actual CDS reader and identical alternative-start outputs while proving source bytes unchanged; XML structure/digest checks read genuine small XML streams. Acquisition, injection matching and full scientific builders are explicitly stubbed where marked. Nextflow adds `test_cds_manifest_paths.py` for actual private copies/inventory and old parent-directory `.fa.gz` alias hazards.

The closed Slurm1514363 v3 run passed38 cases in PriVA and46 in Nextflow (the latter includes existing HPO manifest5), with scoped source hashes unchanged. All8 retained-before method checks deliberately reproduced the old failures, with no fixture errors; those expected failures are labelled separately and the runner exits0. The final two-line lexical XML-alias change then passed all7 affected VCV cases per checkout, again with unchanged source hashes. Counts overlap. Earlier fixture attempts remain retained; they are not qualification runs. Closed logs: `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/logs/installer_final_routes*_1514363_20261017.log` and `installer_final_vcv_alias_1514363_20261017.log`. Scientific parser/filter bodies and CDS calculations are unchanged; only acquisition/provenance/caller I/O changed. These checks do not establish production scientific acceptance. Counterpart tests: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/tests. Nextflow adapter tests: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests/test_cds_manifest_paths.py.

Portable configuration/bundled-data fixtures, 2026-10-17: Slurm allocation1514363 passed38 targeted cases per checkout in `installer-portable-config-{PriVA,PriVA-nextflow}-20261017-v2.log`, then all3 LOEUF size/selection cases per checkout in `installer-portable-loeuf-{PriVA,PriVA-nextflow}-20261017.log`, under `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/logs/`. New modules are `test_installer_config_paths.py`, `test_install_utils_bundled_data.py` and `test_install_utils_bundled_loeuf.py`; they cover basic-path failures before environment changes, immutable effective preflight, old managed-root migration, external selections, opaque nonzero canonical bundles, force protection, empty/missing fallback, assembly aliases, competing external/bundled eRepo input and independent external GoFCards source fetches. Earlier v1/pilot failures are retained as superseded fixture/policy diagnostics. Existing fresh-build fixtures explicitly select external outputs so they continue exercising external validation rather than bundled shortcuts. No downloads, model inference or scientific reader changes occurred. Counterpart tests: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/tests/test_installer_config_paths.py, /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/tests/test_install_utils_bundled_data.py and /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/tests/test_install_utils_bundled_loeuf.py.

Standalone installation regression, 2026-10-17: `test_spliceai2_converter.py` covers private compiler acquisition (mocked), host compiler selection, compiler-free completed-build reuse, tamper refusal and actual streaming allele/gene/strand/order/score/error behavior. `PRIVA_SPLICEAI2_INSTALL_BIN` selects the standalone installation tool; `PRIVA_SPLICEAI2_REFERENCE_BIN` enables exact row/diagnostic comparison against the existing pinned production converter, without changing it. `test_spliceai2_deployment.py` runs real gzip/BGZF/tabix subset publication with the standalone tool and records its independent source/build identity. The installer/header tests cover automatic standalone dispatch and retain completed-adapter reuse. These fixtures do not rerun the full licensed catalogue or qualify ACMG/CADD use.

Closed standalone-converter verification, 2026-10-18 UTC: allocation 1514363 step 651 passed 128 cases in each checkout (256 total; zero failures/errors/skips): converter 13, deployment 26, installer 26, cache-input 42, plugin-identity 5, runtime-guard 14 and existing MaveDB 2. Both registry-route and shell-syntax checks also pass. Actual publisher chr1 prefixes converted 64 SNV and 65 indel-shard rows, including the first MNV found within 11,928 indel source rows; projected rows and JSON diagnostics exactly match frozen converter SHA256 `70eefad75c6ad5a8001fcd87ca6970c3b0b1495f21b55f8c7775db67122c567a`. The source files and reference binary stayed unchanged; no complete source-body hash reread or full catalogue conversion was performed.

The completed receipt, exact source hashes and public-score samples are in `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/standalone-spliceai2-installer-20261017-v4/`; logs match `standalone-spliceai2-installer-*-1514363-20261017-v4.log` under `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/logs/`. The separately completed binary/build receipt remain under `standalone-spliceai2-installer-20261017-v2/tools/5fb2cbf5de613641e0df2fbed78a177c97ae8c5555cbfafe854922545cdd4bc5/` in that validation root. Earlier attempts are retained as diagnostics: v1 was stopped during test-owned compiler acquisition after PATH selected system Rust 1.66; v2/v3 ran inside the agent filename-filtered view, where VCF creation returned ENOENT (v2 also lacked yq). The successful new srun step uses canonical priva_acmg Python 3.12.3/yq and existing Rust 1.90, without compiler acquisition. Tests/source remain in /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests; generated evidence stays in the data-side paths above.

ClinVar candidate boundaries, 2026-10-18: `test_clinvar_cache_candidate.py` exercises small genuine BGZF/tabix files with a deterministic VEP executable fixture: candidate preparation, identical runtime v1 arguments, unchanged plugin order/cutoff, retained production guard, managed output ownership and failure publication. `test_clinvar_gene_refresh.py` adds absent/zero baseline reconstruction, managed reuse, explicit-input errors and failed/concurrent-build preservation. `test_spliceai2_runtime_guard.py` retains its frozen-source command comparison and puts disposable fixtures in `PRIVA_TEST_TMPDIR`/`TMPDIR`. Sources are in `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests`; closed Slurm1514363 logs/receipts are under `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/clinvar-cache-candidate-20261018-v*`. No real model inference, cache downloads or scientific qualification occur in these tests.

ClinVar profile workspace verification, 2026-10-18: final `clinvar-cache-candidate-20261018-v3` passed61 focused cases per checkout (122 total) under Slurm1514363. `test_clinvar_gene_refresh.py` now checks distinct same-release destinations after actual source/reference changes, preserved prior failed/completed workspaces, unchanged-content mtime reuse and matching legacy workspace reuse. Earlier v1/v2 receipts remain separate historical fixture evidence.

Physical SpliceAI2 adoption fixtures, 2026-10-18: `test_spliceai2_relocation.py` uses real small BGZF/tabix files to check independent copied inodes, unchanged original bytes/mtime/ctime, exact scientific metadata preservation, original-receipt retention, root aliases, owned-copy resumption and rejection of linked staging, corrupted bytes/receipts, missing original locks and occupied destinations. Relocated copies validate without opening historical source paths and can be imported again. Sources are in `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests`; generated fixture logs and exact-source receipts are under `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/spliceai2-physical-adoption-20261018-v*`. These tests neither copy the real49GB cache nor rerun its converter or activate scientific use.

Selected plugin input/recovery fixtures, 2026-10-18: `test_indexed_plugin_cache.py` checks real BGZF/TBI repair, preserved zero-based/header settings and exact PrimateAI parity against the original shell pipeline. `test_install_utils_explicit_plugin_inputs.py` checks selected paths, source preservation, managed publication and interruption/retry. Existing LoFtee/SpliceVault/UTR fixtures cover staged resource/repository recovery, unchanged existing checkout identity and configured scratch paths. The final changed-route run passed 44 cases per checkout (88 total) under Slurm 1514363 step 718 with unchanged source; earlier complete helper runs and failed UTR fixture diagnostics remain separate. Mirrored tests are in /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/tests; commands, logs and exact-source receipts are in /stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/plugin-selected-inputs-20261018-v*. These fixtures do not perform full publisher-cache conversion or scientific qualification.

Mapping-bound cache upgrade fixtures, 2026-10-18: `test_mapping_artifact_routes.py` covers current reuse, changed mapping/anchors, old byte preservation, invalid/orphan receipts, occupied outputs, completed-target registration, root aliases and link/root escapes. Scientific builders are bounded argument-recording fixtures; receipt validation and named CLI routing run normally. Existing BioMart and installer prerequisite checks also ran: final v2 passed 82 cases per checkout (164 total), no failures/errors/skips, under Slurm 1514363 step 737 with unchanged selected source. Failed v1 preserves the old alias fixture orphan-receipt diagnostic; fixture resets now remove artifact/receipt pairs together. Mirrored tests are in /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/tests; exact commands, logs and source identities are in /stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/mapping-artifact-routes-20261018-v{1,2}. No model runs or real cache regeneration occurred.

gnomAD annotation failure/retry fixtures, 2026-10-18: `test_gnomad_annotation_failures.py` contains 32 tests with tiny real bcftools/htslib/GNU Parallel data. It checks producer/status propagation, strict joblogs, hidden shard staging, interrupted VCF/TBI publication, failed index repair, occupied index/marker paths, valid index links, expected-shard concat order and exact original rows/INFO/missingness/contig behavior. Closed Slurm1514363 focused versions passed 27, 30 and finally 6 affected tests per checkout (54, 60 and 12 total); these overlapping counts describe separate source versions. Full integrated verification is recorded separately. Exact pre-edit/fixed sources, synthetic inputs/results and closed-save manifests are under `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/gnomad-annotation-propagation-20261018-v1`; matching `gnomad-annotation-propagation-1514363.*-20261018-v*.{json,log}` receipts/logs are in its data-side project's `deployment`/`logs` folders. Mirrored maintained tests are in `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/tests`. No publisher downloads, large cache conversion or scientific activation occurs.

<!-- nextflow-main-driver-20261019-start -->
Nextflow driver automation (2026-10-19): `main_install` calls
`Nextflow_runtime_install CONFIG`, which derives `cache_root/runtime/priva-nextflow-25.10.4`.
The helper exposes `install|check --prefix`; existing pins/launcher are unchanged.
152 affected/dispatch tests passed, followed by six dispatch tests with the added
snapshot-failure case;16 helper fixtures use explicitly fake download/JVM assets.
Real cold/warm/check passed in1514363.1291. These checks do not run a workflow or
qualify models. Maintained source is this checkout; closed logs/receipts are under
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/nextflow-driver-{integration,real-install}-20261019-v1/`,
`logs/nextflow-bootstrap-unit-20261019-v{1,2,3}.log`, and
`deployment/nextflow-driver-main-integration-real-20261019-v1.json`.
The later root-alias correction passed nine focused cases plus real canonical
and alias warm reuse; it preserves lexical normalization and config bytes.
Closed evidence: `deployment/nextflow-driver-root-alias-20261019-v1.json` and
`validation/nextflow-driver-root-alias-20261019-v1/*.log` under the same data project.
<!-- nextflow-main-driver-20261019-end -->

Original producer publication fixtures, 2026-10-24 CT: `test_alphamissense_original_cache.py --fixture-root PATH [--receipt PATH]` retains tiny synthetic inputs and interrupted states in the selected runtime tree. It checks source/input/index/output changes and missingness, native evidence integrity, distinct original/builder receipt schemas, incomplete/short copies, nonowned artifacts, aliases/root escapes, three publication recovery points, identical-content source paths and actual installer routing. Final v3 passed 44 tests per checkout (88 total) under Slurm1514363, one CPU/3GiB; v1 retains one overly narrow diagnostic-regex test failure and v2 passed its then-current 41 tests per checkout. Source syntax checks also passed. No real cache load, original scientific calculation, canonical handoff or full installation acceptance is represented by these fixtures. Mirrored suite: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/tests/test_alphamissense_original_cache.py. Closed fixtures/logs/exact-source receipts: /stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/alphamissense-original-handoff-20261024-v{1,2,3}.

DAS source/receipt portability fixtures, 24 October 2026 CT: `test_das_install_completion.py` now checks identical-content source relocation, same-size/mtime content changes, changed/missing sources, strict selected scientific data, ordered filter parameters, bundled existence/nonzero-only handling, old incompatible receipts, actual other-checkout warm reuse and exact scientific argv order. It also distinguishes a genuine bundled InterPro mapping from an external same-basename file or external symlink target. Each checkout passed22 active tests; one existing optional historical-predecessor reproduction was skipped (23 discovered cases per checkout). Both installer shell syntax checks passed. Tests use a labelled synthetic writer and do not execute domain science or publish a canonical cache. Mirrored tests: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/tests/test_das_install_completion.py. Closed logs, results and exact source snapshots: /stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/das-state-portability-20261024-v1.

Optional Rust condition masks, 24 October 2026 CT: `priva condition-masks --help`
documents the typed CSV boundary for the unchanged `_variant_condition_masks`
function (reference source `09223ef`, file SHA256 `c7e8da80…`). It emits all13
boolean masks in Python dictionary order. `rust/tests/condition_masks.rs` checks
column-phase errors, exact semicolon/final-newline matching, case conversion,
nullable missingness, duplicate labels, generic index labels/names and Polars
input preservation. `condition_masks_oracle.py capture --source-file PATH
--output PATH` executes the exact source function AST with actual pandas;
`condition_masks_compare.py compare --binary PATH --binary-sha256 SHA --oracle
PATH --oracle-manifest-sha256 SHA --output PATH` checks every value, source
Series name, boolean dtype, opaque index token/dtype and single string/None
index name, including exact source exceptions. Per-column dtype declarations
preserve heterogeneous duplicate-label failures; input dtype inference is not
performed. Arbitrary Python objects, categorical/Arrow/time dtypes and CLI
MultiIndex reconstruction remain outside this boundary's contract.

`condition_masks_validate.py prepare|execute` freezes sources/toolchain and runs
the new comparison plus existing mechanism/scoring cases inside Slurm. Generated
fixtures, snapshots, executables, build targets and logs belong under
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/{validation,build}/rust-condition-masks-*`;
the maintained sources/tests are in this checkout. Production still uses Python.
Main-install1524365 uses its separate frozen5:09PM CT source and neither builds
nor selects this optional Rust path. These checks are separate from full
criterion-assignment, ranking, family and installer acceptance.

Frozen validation `rust-condition-masks-20261024-v1` completed in Slurm1524393
at5:40PM CDT,24October2026 (compute,2CPUs/8GiB,2minutes, native exit0).
All70 Rust tests passed:16 condition,9 mechanism,9 scoring and36 existing
SpliceAI2. Actual-Python condition comparison passed81 cases/1,071rows,
including24 expected exceptions, plus4 transport guards; values and metadata
matched. The retained mechanism comparison passed133 cases/2,349rows and
scoring matched19,754rows plus16 failure/overwrite checks. No mismatch occurred.
Source manifest SHA256 `2542f02f3a9a88a58ddbab18a31d38a18d9ae7e8734a1240e1c97e6d9fb1f8e9`;
the separate candidate binary SHA256 is
`2f8d7945610d80a99151769aa5c2622cf37ac1818210f6b1f6c336e9007f3b27`.
The validation controller retains exact source/toolchain/oracle pins, individual
comparison results and logs in the data-side validation subtree above. Build
targets/private Cargo state stay in its data-side build counterpart.

DAS pool memory repair, 24 October 2026 CT: `test_das_spawn_pools.py` exposes
named `prepare`/`run` commands (usage in its docstring) and freezes the actual
pre-fix source and candidate dependencies. Eleven bounded real-worker fixtures
passed in Slurm1514363.3040: visualization, KS, both dip rescues and chromosome
metadata retain exact task payloads, seeds, per-task numerical returns, PNG bytes,
metadata order, input nonmutation and error propagation. All five independent
pools now select spawn explicitly; the global start method remains unchanged.
Existing metadata-cache and DAS-completion suites also passed46 active checks
per checkout in Slurm1514363.3041 (47 discovered, one optional historical skip).
The five constructors and one explanatory comment are the complete source diff
from frozen SHA256 `57cc7be79168bcbc10aebc12efe94ea019fd33c148f35b5ba65360fc29b648c3`.
The two checkouts and all tested dependencies are byte-identical. Native unordered
completion sequences and table files are retained; no sorting or seed-policy
change was introduced. These small fixtures do not establish a full-cache RAM
bound, deterministic full-run row order or whole-installation acceptance.
Maintained mirrored source/tests are under
`/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/{PriVA,PriVA-nextflow}/{scripts,tests}`;
all exact-source snapshots, synthetic outputs and logs are under
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/das-spawn-pools-20261024-v1`.

Optional Rust pext conversion, 24 October 2026 CT: `priva coerce-pext --input
INPUT.csv --output NEW.csv` exposes only the pinned
`acmg_consequence.py::coerce_pext_value` scalar boundary. The primitive Polars
adapter returns Float64, retains the column name and row order, and emits IEEE
NaN for missing results. Callers explicitly select Python or NumPy boolean
boxing: the original treats Python `bool` as numeric and NumPy `bool_` as text.
The typed CSV retains opaque primitive index tokens, including duplicates;
it neither infers a pandas dtype nor reconstructs an arbitrary pandas index.

The actual source checks `pd.isna` first, converts numeric scalars directly,
then parses comma/ampersand text and takes Python's ordered `max`. Its docstring
says finite, but its code retains infinities and some text NaNs. Tests compare
float bits, signed zero, NaN sign/order, exact skip tokens, decimal/underscore
grammar and the frozen CPython3.12.3 Unicode15 decimal digits. Big integer
conversion preserves supported source OverflowErrors. Malformed transport is a
separate guard failure. Objects, arrays, bytes, complex, float16/longdouble,
categorical/time/decimal/list/struct dtypes and surrogate strings require the
unchanged Python implementation; no automatic fallback or production dispatch
was added.

`rust/tests/coerce_pext_oracle.py capture --source-file PATH --output PATH`
captures the unchanged original AST with real pandas/NumPy;
`coerce_pext_validate.py prepare|execute` freezes sources/toolchain and reuses
that previously captured oracle for exact native comparison, plus the existing
condition/mechanism/scoring cases. Maintained Rust, oracle and test sources stay
in this checkout. Generated oracles, source snapshots, candidate executable and
logs are under
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/rust-coerce-pext-20261024-v1`;
private Cargo/target/temp products stay in its data-side build counterpart
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/build/rust-coerce-pext-20261024-v1`.
Main1524420 uses separate frozen6:39PM CT sources and does not build or select
this later optional module. This boundary does not establish installation,
criterion-assignment, ranking, family or public scientific acceptance.

Frozen coerce-pext validation1524423 completed0:0 at7:04:04PM CDT,
24October2026 (compute,2CPUs/8GiB,2m12s; batch MaxRSS628,044KiB).
All84 Rust tests passed, including14 new scalar/Polars tests. The oracle was
captured before implementation in1514363.3099; all5,786 supported scalar
outputs matched exact f64 bits, and18 source OverflowErrors matched category
and message. Eight malformed/unsupported transport guards plus overwrite,
repeated-option, empty and duplicate-label checks passed. Nine arbitrary-object
diagnostics remain explicitly outside parity. Existing85 condition cases,
133 mechanism cases and19,754 scoring rows plus16 failure checks also passed.
No mismatch occurred. Source manifest SHA256
`9f11a3bb3527e555321d07e472557508f80bd601b3c7649ba534b661f44e6dd4`;
candidate binary SHA256
`a3b3bb42d69f5b98f9a074f7c5e3626fd72f43d080401a9b6f4ce7a0443214f8`.
The original validation JSON retains inherited environment step label0;
scheduler closure identifies the actual1524423 batch wrapper, with no numeric
step0 accounting record. Closed intake records this distinction explicitly.

Full-ranking/control oracle groundwork, 24 October 2026 CT:
`rust/tests/ranking_control_oracle.py prepare|capture|validate` freezes actual
`09223ef` modules and captures unchanged Python execution with real SciPy.
Final `ranking-control-oracle-20261024-v3` passed 341 checks across 165 paired
traced/untraced cases (117 values, 48 original exceptions), preserving typed
tables, input mutation/serialization, scientific arguments, source checkpoints
and first failure phases. The preceding v1 observer-cache mismatch and failed
diagnostic logs, plus preparation-only v2, remain retained. The matrix covers
all 16 control-column patterns, true count/AF route priority, CDF/Beta calls,
threshold neighbors and R1–R8 ranking/caller edge cases. The current-progress
addendum in [the source ledger](../docs/ranking-source-ledger-09223ef.md) records
the exact scope, identities and the original mixed-AF object-mask failure.
Maintained harness/docs are in this checkout; all snapshots, synthetic inputs,
pickles, typed states, traces and logs reside under
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/{validation/ranking-control-oracle-20261024-v*,logs/ranking-control-oracle-*}`.
No Rust integration, scientific Python/installer change or production selection
was performed; this fixture oracle does not qualify a family or full pipeline.

Opt-in R8 terminal ordering, 24 October 2026 CT:
[`terminal_order_oracle.py`](../rust/tests/terminal_order_oracle.py) exposes
`capture --full-oracle CLOSED_ROOT --output NEW_ROOT`; the separate integer
supplement exposes `capture --base-oracle CLOSED_R8_ROOT --output NEW_ROOT`.
[`terminal_order_validate.py`](../rust/tests/terminal_order_validate.py)
provides named `prepare`/`execute` procedures, source/toolchain/oracle pins and
an explicit reuse option that requires identical Rust/scientific source. It
retains whole-run failures even when reusing completed compilation/tests.

`priva terminal-order --input INPUT.csv --output NEW.csv --metadata-output NEW.csv
--column-dtype POSITION=DTYPE` is a separate typed boundary; zero-based positions
refer to data columns after the index token pair. It emits a source-row
permutation, temporary maxima and ranks plus a column/axis plan. The adapter
gathers original payload buffers and preserves actual empty RangeIndex metadata.
The existing `final-order` body and interface remain unchanged. Supported dtypes,
missingness, integer promotion and explicit Python fallbacks are detailed in
[the R8 ledger addendum](../docs/ranking-source-ledger-09223ef.md).

Final closure 1514363.3253 passed **84 exact primary outputs / six original
failures / one Python fallback**, and separately **15 exact integer outputs /
zero failures / zero fallbacks**. The actual Int64 checkpoint appears in both
sets. All 97 Rust tests, prior PEXT/condition/mechanism/scoring comparisons and
60 old `final-order` fixtures/five rejections passed. Initial extractor failures,
the Int64 deferral failure 1524435, narrower pass 1524440 and empty-axis adapter
failure 1524441 remain preserved. Sources/build provenance and binary SHA are in
the ledger; later Rust/test sources are separate from main installation 1524420's frozen
6:39 PM CT installer snapshot. All fixtures/source snapshots/binaries/logs are
under data-side `validation/rust-terminal-order-20261024-v*` and
`logs/rust-terminal-order-*`; Cargo/targets are excluded under `build/`.
No C1–R7 numerics, scientific Python, installer or production dispatch changed;
these tests do not qualify a family or complete pipeline.


Opt-in R1 copy/drop/column-copy/filter plan, 24 October 2026 CT:
[`ranking_setup_oracle.py`](../rust/tests/ranking_setup_oracle.py) captures the
four unchanged original statements before Rust using actual full-ranker R2-entry
checkpoints and R1 failures. [`ranking_setup_validate.py`](../rust/tests/ranking_setup_validate.py)
exposes named `prepare`/`execute` procedures with source, toolchain and oracle
pins; guarded binary reuse permits only the two R1 test-harness files to change.
`priva rank-setup --input INPUT.csv --output NEW.csv --metadata-output NEW.csv
--column-dtype POSITION=DTYPE` returns source row positions and an ordered
source-column copy/drop plan. The comparator applies every emitted position to
original payload buffers and retains opaque axes/metadata. R1 performs no score
arithmetic, PED lookup or C1/R2–R7 numerical work.

Capture 1514363.3276 retains 129 cases (118 original outputs/11 failures), with
unchanged immediate source input pickle bytes. Final 1514363.3308 completed 0:0:
**117 exact outputs/eight original failures/four Python fallbacks**, with zero
mismatches. Duplicate required/assigned labels are explicit fallbacks, separate
from supported original failures; unrelated duplicate labels retain positions.
Build 1524493.0 passed all 105 Rust tests and prior R8/PEXT/condition/mechanism/
scoring/old-final-order regressions. The final check reused its unchanged binary
after strengthening column-plan application in the Python adapter. Both checks
remain closed evidence. Source/binary pins and the exact primitive/fallback
contract are in [the R1 ledger addendum](../docs/ranking-source-ledger-09223ef.md).
Fixtures, source snapshots, binary and logs stay under the data-side
`validation/rust-ranking-setup-20261024-v1` and `logs/rust-ranking-setup-*`;
Cargo/targets stay under `build/rust-ranking-setup-*`. The R8 algorithm/tests,
scientific Python, installer and production/family dispatch are unchanged. This
later optional module is separate from main 1524420's frozen installer snapshot
and does not establish complete ranking, family or public scientific acceptance.

Opt-in R2 pedigree context/mask, 24 October 2026 CT:
[`pedigree_context_oracle.py`](../rust/tests/pedigree_context_oracle.py) captures
the actual source statements L501–524 and full-ranker R3-entry locals/failures.
[`pedigree_context_validate.py`](../rust/tests/pedigree_context_validate.py)
exposes named `prepare`/`execute` procedures with source/toolchain/oracle pins.
`priva pedigree-context --input INPUT.csv --family-input SCALAR.csv
--pedigree-input PED.csv --output NEW.csv --metadata-output NEW.csv
--column-dtype POSITION=DTYPE --pedigree-column-dtype POSITION=DTYPE` returns
source pedigree position, primitive context, raw Boolean mask/null bits and
structured warning metadata. Omitting pedigree input preserves the original
no-PED branch. The mask counts `1` across the complete sample string.

The comparator reconstructs each selected row through a one-row DataFrame take
and checks full source row/Series metadata, raw nullable-Boolean data beneath
null bits, warning arguments/order, original failure phase/type/message and
immediate input pickle bytes. Final 1524517/.0 completed 0:0: **197 exact
outputs/53 original failures/28 deliberate Python fallbacks**, zero mismatches;
all preceding procedure differentials passed again. All 118 Rust tests passed
in 1524516.0 before its later adapter failure. That whole failure and the prior
empty-selection failure 1524515 remain preserved. The final check reused only
the unchanged binary/completed Rust tests, then reran every differential phase.

Numeric nullable sample accessor/cache effects, NumPy family scalar/object and
native numeric family-column semantics, ambiguous PED labels, nontext proband
IDs and arbitrary objects explicitly retain Python. Two actual untraced source
dtype-cache serialization changes remain recorded, separate from value equality.
Source/binary pins, supported types and the sex-predicate-only boundary are in
[the R2 ledger addendum](../docs/ranking-source-ledger-09223ef.md). Maintained
sources stay here; all fixtures/snapshots/binaries/logs are under data-side
`validation/rust-pedigree-context-20261024-v*` and `logs/rust-pedigree-context-*`,
with build/Cargo caches separately under `build/rust-pedigree-context-*`.
R1/R8 algorithms/tests, numerical/scientific Python, installer and production
dispatch remain unchanged; this module is separate from main 1524420's frozen
snapshot and provides no complete-ranking/family/public qualification.

Opt-in R3 control multiplication/assignment, 24 October 2026 CT:
[`control_multiply_oracle.py`](../rust/tests/control_multiply_oracle.py) captures
actual L528 after unchanged Python C1–C5/SciPy, bound to full-ranker checkpoints.
[`control_multiply_validate.py`](../rust/tests/control_multiply_validate.py)
provides named `prepare`/`execute` procedures and guarded binary/test reuse.
`priva control-multiply --input INPUT.csv --output NEW.csv --metadata-output NEW.csv
--row-index-kind flat --column-dtype POSITION=DTYPE --axis-metadata TOKEN`
emits exact score bits, source positions/index tokens and the existing target
column position. The public Rust function mutates that score column in memory;
the CSV boundary is for optional procedure validation, not a production process
per source statement.

Qualified inputs use native `float64` operands, flat opaque axes and at most
1,000,000 rows. Other dtypes/promotion, ambiguous required labels, MultiIndex
warnings, larger numexpr routes, unsupported objects and two-NaN operands defer
before mutation. The actual Int64/nullable-Float64 checkpoints remain Python.
Fallback runs only original L528 on current buffers; a future full-Python
dispatcher must decide before any prior stage executes. Earlier warnings,
C1–C5 and mutations must never be replayed after a partial Rust workflow.

Final 1524524/.0 completed 0:0 at 10:54 PM CDT: **88 exact outputs/three original
failures/25 Python fallbacks**, zero mismatches. Six deferred source failures
remain separate from the three supported failures. The 116-case oracle and the
separate preceding C1–C5 failure are retained. All 127 Rust tests passed in
1524519/.0 before its later observer-rehydration failure; final validation reused
the unchanged binary and reran every earlier procedure differential successfully.
Both failed attempts remain FAILED. The first exposed exact two-NaN payload
differences; the second exposed dtype-cache key-string aliases in serialized
observer state. Corrected tests reconstruct genuine prepared entry state and
match immediate before/after bytes, float bits, metadata, warnings and failure
phase without resetting caches or waiving mismatches.

Source/binary pins and the exact fallback contract are in
[the R3 ledger addendum](../docs/ranking-source-ledger-09223ef.md). Fixtures,
diagnoses, snapshots, binary and logs stay under data-side
`validation/rust-control-multiply-20261024-v1` and `logs/rust-control-multiply-*`;
Cargo/targets stay separately under `build/rust-control-multiply-*`. R1/R2/R8
algorithms/tests, scientific Python, installer and production dispatch are
unchanged. This later opt-in source is separate from main 1524420's frozen
snapshot and establishes no full-ranking/family/public acceptance.


Opt-in R4 compatibility masks/column assignment, 24 October 2026 CT:
[`ranking_compatibility_oracle.py`](../rust/tests/ranking_compatibility_oracle.py)
captures actual original L537–597 and genuine full-function L537/L599 states
before Rust. Its `run_current` observer accepts the actual current buffers;
fallback never rehydrates an entry snapshot. The named `prepare`/`execute`
procedures in [`ranking_compatibility_validate.py`](../rust/tests/ranking_compatibility_validate.py)
freeze source/toolchain/oracles and run all earlier procedure regressions.

`priva ranking-compatibility --input INPUT.csv --heterozygosity-input MASK.csv
--heterozygosity-dtype bool --proband-is-male true --proband-is-female false
--output NEW.csv --metadata-output NEW.csv --column-dtype POSITION=DTYPE`
returns ordered mechanism/compatibility masks, raw Boolean/null bits, opaque
metadata origins, alias roles, exact INFO counts and existing target position.
The public Rust API directly consumes the existing R2 Context in memory. This
CSV boundary is for optional validation, not a production process per source line.
Only compatibility initialization/assignment occurs here; multiplication stays
at R6 after original gene-list operations.

Final 1524530/.0 completed 0:0 at 11:41 PM CDT: **103 exact outputs/nine original
failures/35 Python fallbacks**, zero mismatches; 14 deferred source failures stay
in the fallback count. All 136 Rust tests (nine R4 tests), R3/R2/R1/R8/integer,
PEXT/condition/mechanism/scoring and old-final-order regressions passed. The oracle
retains 147 cases and eight separate preceding-stage failures. The first observer
phase-label failure and failed Rust attempt 1524529 remain preserved. The latter
exposed three numeric-conversion attrs fields, fixed with an explicit fresh
metadata plan and regression; no semantic or byte tolerance was introduced.

All source masks, raw storage, attrs/flags/axes/aliases, INFO/warnings, column
positions and immediate input mutation bytes are compared. Unsupported nullable
helper/raw-score cache, chromosome/object and index-alignment representations
retain the original R4 stage on current buffers, without prior-stage replay.
The exact domain, source/binary pins and fallback rule are in
[the R4 ledger addendum](../docs/ranking-source-ledger-09223ef.md). Generated
states/snapshots/binaries/diagnostics and logs stay under data-side
`validation/rust-ranking-compatibility-20261024-v*` and
`logs/rust-ranking-compatibility-*`; build/Cargo/targets stay under the separate
`build/rust-ranking-compatibility-*`. Scientific Python, installers, earlier
procedure bodies, C1–C5/SciPy, R6 and production dispatch remain unchanged. This
later source is separate from main 1524420's frozen snapshot and does not qualify
full ranking, a family or a public benchmark.


### Optional R5→R6 gene-weight procedure — 25 October 2026 CT

`priva ranking-gene-weights --help` documents the typed fixture boundary;
ordered Rust workflows use `ranking_gene_weights` and
`Penalty::from_compatibility` in memory. It preserves dispensable 0.9, relevant
1.2 and selected compatibility multiplication as three separate operations.
`None`, selected empty files, exact symbol membership, column positions, opaque
axes and raw score/null bits are distinct. Only a literal plain-text parser
subset and native float64 route are qualified. Unsupported pandas parsing,
inference, objects/dtypes/errors, alignment and NaN-payload routes retain the
original block once on current R5-entry buffers, with no prior-stage replay.

The source-only oracle has 203 cases, including 77 actual full-function L599/L621
boundaries. Final native1514363.3624 completed0:0 at12:24:56AMCDT: **120 exact
outputs and83 Python fallbacks (28 original failures)**, zero mismatches.
All147 Cargo tests passed, including11 new tests and direct R2→R4→R5/R6 API
composition. The binary/tests are reused under unchanged code/toolchain pins
from failed1524546; only comparator/validator changed. Its original fallback
pickle failure and the earlier capture/format failures remain preserved.

The fallback mismatch was a genuine process-global pandas Float64 dtype cache
context: actual original R8 sorting creates `index_class`. The successor records
a separate original full-source fixture before any paired R5 buffers/Rust case;
cold/warm original proofs, exact mutation bytes/aliases and all source warnings
remain visible. No cache reset or tolerance was used. Earlier Python
qualifications are bound by unchanged body/test hashes rather than rerun.

`ranking_gene_weights_{oracle,compare,validate}.py` expose named capture/inspect,
compare and format/prepare/execute procedures. The [R5→R6 ledger
addendum](../docs/ranking-source-ledger-09223ef.md) records the exact contract,
binary/source/oracle pins and exclusions. Generated evidence/logs stay under
data-side `validation/rust-ranking-gene-weights-20261024-v{1,2}` and
`logs/rust-ranking-gene-weights-*`; build/Cargo/targets stay separately under
`build/rust-ranking-gene-weights-*`. Scientific Python, C1–C5/SciPy, R7, installers
and production dispatch remain unchanged; no full-ranking/family qualification.


### Optional R7 PEXT procedure — 25 October 2026 CT

`priva ranking-pext --help` describes the fixture boundary. Ordered native
workflows use `ranking_pext` and `Modulation::from_compatibility` in memory on
current R6 buffers. It retains first-tissue/mean selection, the existing scalar
parser, parameter/log/assignment order, original axes/column positions and
separate floating operations. Explicit boxed scalar origins preserve pandas map
behavior. Qualified applied inputs have 1–1M rows, flat opaque axes, native f64
scores/valid finite float parameters and an aligned Boolean mask. Unsupported
representations/warnings/numerics defer atomically to original R7 once on current
buffers, without repeating prior stages. No selected column skips unused checks.

**1524555/.0 completed0:0 at1:17:21AMCDT:** 162 Rust tests (15 R7), 154 exact native
outputs and47 Python fallbacks, including33 original failures; zero mismatches
across201 cases. Compare checks include local/input pickle bytes and aliases,
raw intermediate/score bits, nulls, dtype/axes/attrs/flags, logs/warnings and
failure phase/partial mutation. The root170-case original oracle and31 added
edges retain actual pandas-boxed parser inputs.

Failed1524551,1524552 and1514363.3687 remain preserved. Their distinct fixture,
metadata-name alias and actual Series-power dispatch findings were not waived.
R7's Series `np.power` takes the original scalar exponent1/0.5/2 branches;
the ndarray-only probe is separately scoped. Actual-Series252-case/807012-value
proof matches source scalar dispatch exactly at opt0/opt2, while generic powf
has32 mismatches. Source cutoff-neighbor and signed-zero witnesses are regressions.

`ranking_pext_{oracle,compare,validate}.py` expose capture/inspect, compare and
format/prepare/execute procedures. The [R7 ledger addendum](../docs/ranking-source-ledger-09223ef.md)
records exact domain, native/source/binary pins, dispatch provenance and fallback.
Generated evidence/logs live on the data side under
`validation/{rust-ranking-pext,ranking-pext-*}-20261025-v1` and matching logs;
Cargo/targets remain separately under `caches/runtime/ranking-pext-build-*`.
Earlier bodies/tests, scientific Python, C1–C5/SciPy, installers and production
selection remain unchanged. This later optional source does not qualify a
complete ranking workflow, family or public benchmark.
