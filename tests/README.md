# PriVA tests

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

Configured HPO path checks, 2026-10-17: `test_hpo_runtime_paths.py` covers explicit/default routing, missing inputs, quoted paths, freshness and main-call order; `test_installer_prerequisites.py` checks data-side builder targets and existing path guards. Counterpart tests: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests/test_hpo_runtime_paths.py and /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests/test_installer_prerequisites.py. The Nextflow manifest identity/default/safety checks are in /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests/test_hpo_manifest_paths.py. Conversion stubs test shell routing separately; the scientific pandas implementation is unchanged. All new fixtures/logs use the authorized data-side project root.

Explicit output-route checks, 2026-10-17: `test_install_output_paths.py` separates stubbed scientific parsers from actual gzip/SHA256, path, timestamp and native awk behavior. It covers exact data-side paths, independent ClinGen destinations, retained InterPro source/cache reuse, failure barriers, legacy compatibility and unchanged curated-source payloads. The manual InterPro report is checked separately from cache metadata. Counterpart tests: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests/test_install_output_paths.py. Together with existing download-failure tests, 33 checks passed per checkout under Slurm1514363; all generated fixtures/logs use the authorized data-side project root.

Managed-reference cache checks, 2026-10-17: `test_reference_resources.py` now covers installer opt-in versus direct `--data-dir`, all four assembly aliases, new managed defaults, exact legacy-default redirection, preserved original/old-managed identities, missing custom inputs, orphan/staging retention, chain direction and matching preflight/install plans. Tiny real FASTAs use native samtools indexing; HTTP bytes are explicit fixtures. Existing merged-cache metadata fixtures now include the required official-format species/assembly/tabix fields and index. All 42 focused checks passed in each checkout under Slurm1514363, with new fixtures/logs under the authorized data-side project. Counterpart tests: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests/test_reference_resources.py. No live reference/index, production configuration or frozen job was changed.

Final installer routing/I/O verification, 2026-10-17: `test_install_utils_{final_routes,cds_io}.py` and `test_selected_clinvar_vcv.py` cover UTR cwd, exact CDS/XML paths, missing/invalid inputs, retry/configuration failures, GRCh38 metadata aliasing, custom XML content identity and immutable injection/full-build gates. Small genuine gzip/plain FASTAs exercise the actual CDS reader and identical alternative-start outputs while proving source bytes unchanged; XML structure/digest checks read genuine small XML streams. Acquisition, injection matching and full scientific builders are explicitly stubbed where marked. Nextflow adds `test_cds_manifest_paths.py` for actual private copies/inventory and old parent-directory `.fa.gz` alias hazards.

The closed Slurm1514363 v3 run passed38 cases in PriVA and46 in Nextflow (the latter includes existing HPO manifest5), with scoped source hashes unchanged. All8 retained-before method checks deliberately reproduced the old failures, with no fixture errors; those expected failures are labelled separately and the runner exits0. The final two-line lexical XML-alias change then passed all7 affected VCV cases per checkout, again with unchanged source hashes. Counts overlap. Earlier fixture attempts remain retained; they are not qualification runs. Closed logs: `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/logs/installer_final_routes*_1514363_20261017.log` and `installer_final_vcv_alias_1514363_20261017.log`. Scientific parser/filter bodies and CDS calculations are unchanged; only acquisition/provenance/caller I/O changed. These checks do not establish production scientific acceptance. Counterpart tests: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests. Nextflow adapter tests: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests/test_cds_manifest_paths.py.

Portable configuration/bundled-data fixtures, 2026-10-17: Slurm allocation1514363 passed38 targeted cases per checkout in `installer-portable-config-{PriVA,PriVA-nextflow}-20261017-v2.log`, then all3 LOEUF size/selection cases per checkout in `installer-portable-loeuf-{PriVA,PriVA-nextflow}-20261017.log`, under `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/logs/`. New modules are `test_installer_config_paths.py`, `test_install_utils_bundled_data.py` and `test_install_utils_bundled_loeuf.py`; they cover basic-path failures before environment changes, immutable effective preflight, old managed-root migration, external selections, opaque nonzero canonical bundles, force protection, empty/missing fallback, assembly aliases, competing external/bundled eRepo input and independent external GoFCards source fetches. Earlier v1/pilot failures are retained as superseded fixture/policy diagnostics. Existing fresh-build fixtures explicitly select external outputs so they continue exercising external validation rather than bundled shortcuts. No downloads, model inference or scientific reader changes occurred. Counterpart tests: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests/test_installer_config_paths.py, /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests/test_install_utils_bundled_data.py and /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests/test_install_utils_bundled_loeuf.py.

Standalone installation regression, 2026-10-17: `test_spliceai2_converter.py` covers private compiler acquisition (mocked), host compiler selection, compiler-free completed-build reuse, tamper refusal and actual streaming allele/gene/strand/order/score/error behavior. `PRIVA_SPLICEAI2_INSTALL_BIN` selects the standalone installation tool; `PRIVA_SPLICEAI2_REFERENCE_BIN` enables exact row/diagnostic comparison against the existing pinned production converter, without changing it. `test_spliceai2_deployment.py` runs real gzip/BGZF/tabix subset publication with the standalone tool and records its independent source/build identity. The installer/header tests cover automatic standalone dispatch and retain completed-adapter reuse. These fixtures do not rerun the full licensed catalogue or qualify ACMG/CADD use.

Closed standalone-converter verification, 2026-10-18 UTC: allocation 1514363 step 651 passed 128 cases in each checkout (256 total; zero failures/errors/skips): converter 13, deployment 26, installer 26, cache-input 42, plugin-identity 5, runtime-guard 14 and existing MaveDB 2. Both registry-route and shell-syntax checks also pass. Actual publisher chr1 prefixes converted 64 SNV and 65 indel-shard rows, including the first MNV found within 11,928 indel source rows; projected rows and JSON diagnostics exactly match frozen converter SHA256 `70eefad75c6ad5a8001fcd87ca6970c3b0b1495f21b55f8c7775db67122c567a`. The source files and reference binary stayed unchanged; no complete source-body hash reread or full catalogue conversion was performed.

The completed receipt, exact source hashes and public-score samples are in `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/standalone-spliceai2-installer-20261017-v4/`; logs match `standalone-spliceai2-installer-*-1514363-20261017-v4.log` under `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/logs/`. The separately completed binary/build receipt remain under `standalone-spliceai2-installer-20261017-v2/tools/5fb2cbf5de613641e0df2fbed78a177c97ae8c5555cbfafe854922545cdd4bc5/` in that validation root. Earlier attempts are retained as diagnostics: v1 was stopped during test-owned compiler acquisition after PATH selected system Rust 1.66; v2/v3 ran inside the agent filename-filtered view, where VCF creation returned ENOENT (v2 also lacked yq). The successful new srun step uses canonical priva_acmg Python 3.12.3/yq and existing Rust 1.90, without compiler acquisition. Tests/source remain in /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/tests; generated evidence stays in the data-side paths above.

ClinVar candidate boundaries, 2026-10-18: `test_clinvar_cache_candidate.py` exercises small genuine BGZF/tabix files with a deterministic VEP executable fixture: candidate preparation, identical runtime v1 arguments, unchanged plugin order/cutoff, retained production guard, managed output ownership and failure publication. `test_clinvar_gene_refresh.py` adds absent/zero baseline reconstruction, managed reuse, explicit-input errors and failed/concurrent-build preservation. `test_spliceai2_runtime_guard.py` retains its frozen-source command comparison and puts disposable fixtures in `PRIVA_TEST_TMPDIR`/`TMPDIR`. Sources are in `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/tests`; closed Slurm1514363 logs/receipts are under `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/clinvar-cache-candidate-20261018-v*`. No real model inference, cache downloads or scientific qualification occur in these tests.

ClinVar profile workspace verification, 2026-10-18: final `clinvar-cache-candidate-20261018-v3` passed61 focused cases per checkout (122 total) under Slurm1514363. `test_clinvar_gene_refresh.py` now checks distinct same-release destinations after actual source/reference changes, preserved prior failed/completed workspaces, unchanged-content mtime reuse and matching legacy workspace reuse. Earlier v1/v2 receipts remain separate historical fixture evidence.

Physical SpliceAI2 adoption fixtures, 2026-10-18: `test_spliceai2_relocation.py` uses real small BGZF/tabix files to check independent copied inodes, unchanged original bytes/mtime/ctime, exact scientific metadata preservation, original-receipt retention, root aliases, owned-copy resumption and rejection of linked staging, corrupted bytes/receipts, missing original locks and occupied destinations. Relocated copies validate without opening historical source paths and can be imported again. Sources are in `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA/tests`; generated fixture logs and exact-source receipts are under `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/spliceai2-physical-adoption-20261018-v*`. These tests neither copy the real49GB cache nor rerun its converter or activate scientific use.

Selected plugin input/recovery fixtures, 2026-10-18: `test_indexed_plugin_cache.py` checks real BGZF/TBI repair, preserved zero-based/header settings and exact PrimateAI parity against the original shell pipeline. `test_install_utils_explicit_plugin_inputs.py` checks selected paths, source preservation, managed publication and interruption/retry. Existing LoFtee/SpliceVault/UTR fixtures cover staged resource/repository recovery, unchanged existing checkout identity and configured scratch paths. The final changed-route run passed 44 cases per checkout (88 total) under Slurm 1514363 step 718 with unchanged source; earlier complete helper runs and failed UTR fixture diagnostics remain separate. Mirrored tests are in /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests; commands, logs and exact-source receipts are in /stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/plugin-selected-inputs-20261018-v*. These fixtures do not perform full publisher-cache conversion or scientific qualification.

Mapping-bound cache upgrade fixtures, 2026-10-18: `test_mapping_artifact_routes.py` covers current reuse, changed mapping/anchors, old byte preservation, invalid/orphan receipts, occupied outputs, completed-target registration, root aliases and link/root escapes. Scientific builders are bounded argument-recording fixtures; receipt validation and named CLI routing run normally. Existing BioMart and installer prerequisite checks also ran: final v2 passed 82 cases per checkout (164 total), no failures/errors/skips, under Slurm 1514363 step 737 with unchanged selected source. Failed v1 preserves the old alias fixture orphan-receipt diagnostic; fixture resets now remove artifact/receipt pairs together. Mirrored tests are in /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests; exact commands, logs and source identities are in /stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/mapping-artifact-routes-20261018-v{1,2}. No model runs or real cache regeneration occurred.

gnomAD annotation failure/retry fixtures, 2026-10-18: `test_gnomad_annotation_failures.py` contains 32 tests with tiny real bcftools/htslib/GNU Parallel data. It checks producer/status propagation, strict joblogs, hidden shard staging, interrupted VCF/TBI publication, failed index repair, occupied index/marker paths, valid index links, expected-shard concat order and exact original rows/INFO/missingness/contig behavior. Closed Slurm1514363 focused versions passed 27, 30 and finally 6 affected tests per checkout (54, 60 and 12 total); these overlapping counts describe separate source versions. Full integrated verification is recorded separately. Exact pre-edit/fixed sources, synthetic inputs/results and closed-save manifests are under `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/gnomad-annotation-propagation-20261018-v1`; matching `gnomad-annotation-propagation-1514363.*-20261018-v*.{json,log}` receipts/logs are in its data-side project's `deployment`/`logs` folders. Mirrored maintained tests are in `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests`. No publisher downloads, large cache conversion or scientific activation occurs.

Original producer publication fixtures, 2026-10-24 CT: `test_alphamissense_original_cache.py --fixture-root PATH [--receipt PATH]` retains tiny synthetic inputs and interrupted states in the selected runtime tree. It checks source/input/index/output changes and missingness, native evidence integrity, distinct original/builder receipt schemas, incomplete/short copies, nonowned artifacts, aliases/root escapes, three publication recovery points, identical-content source paths and actual installer routing. Final v3 passed 44 tests per checkout (88 total) under Slurm1514363, one CPU/3GiB; v1 retains one overly narrow diagnostic-regex test failure and v2 passed its then-current 41 tests per checkout. Source syntax checks also passed. No real cache load, original scientific calculation, canonical handoff or full installation acceptance is represented by these fixtures. Mirrored suite: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests/test_alphamissense_original_cache.py. Closed fixtures/logs/exact-source receipts: /stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/alphamissense-original-handoff-20261024-v{1,2,3}.

DAS source/receipt portability fixtures, 24 October 2026 CT: `test_das_install_completion.py` now checks identical-content source relocation, same-size/mtime content changes, changed/missing sources, strict selected scientific data, ordered filter parameters, bundled existence/nonzero-only handling, old incompatible receipts, actual other-checkout warm reuse and exact scientific argv order. It also distinguishes a genuine bundled InterPro mapping from an external same-basename file or external symlink target. Each checkout passed22 active tests; one existing optional historical-predecessor reproduction was skipped (23 discovered cases per checkout). Both installer shell syntax checks passed. Tests use a labelled synthetic writer and do not execute domain science or publish a canonical cache. Mirrored tests: /stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow/tests/test_das_install_completion.py. Closed logs, results and exact source snapshots: /stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/das-state-portability-20261024-v1.

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
