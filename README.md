# PriVA

<!-- priva-resumed-current-20261024 -->
**Installation passed; migration continues — 12:19 AM CDT, 25 October 2026.** Both configured
GRCh38/VEP116 literal main_install runs passed, in order, with each stage's own
native CADD component check and final VEP/PEXT/resource checks. Main1524420
completed0:0 at12:08:54 AM CDT; terminal1524421 completed0:0 at12:09:56 AM,
25 October. Its closed receipt records installation_accepted=true; root verified
the closed native rows and current installer-script identities. This reused
completed caches; no fresh-host rebuild was required.

Original AlphaMissense statistics, DAS/PM1 and ClinVar AF are deployed. The
optional exact merger stays isolated. No installer job remains live or required
input is missing in this run. Scientific and family acceptance remain separate:
full UTR qualification1524532/1524533, complete CADD feature/profile qualification,
remaining Rust/Nextflow integration, family sorting/ACMG and public benchmark are
pending. The three-allele real-model CADD boundary comparison1524545 passed with
exact features, raw/PHRED scores and ordered outputs; this is a bounded control.
Keep original CADD Python; cadd-rs stays suspended. HGNC_ID-only metadata
mismatches are accepted; other scientific rows, joins, order and scores remain
protected.

The maintained [project status](/stornext/snfs170/next-gen/scratch/yangyxt/projects/PriVA_upgrade/README.md)
and [migration plan](/stornext/snfs170/next-gen/scratch/yangyxt/projects/PriVA_upgrade/migration_plan.html)
carry current evidence and the remaining goal. Earlier dated job states and
pause statements below are historical and superseded by the explicit resume.
<!-- priva-resumed-current-20261024-end -->

<!-- utrannotator-local-reporting-20261024 -->
**UTRAnnotator local reporting fix — 24 October 2026, 11:33 AM CDT.** Installed profile
`utrannotator-116-all-effects-20261024-v1` preserves the five calculations and
reports every active effect with paired details in fixed order. Both checkouts
also fix per-effect dispatch. Slurm validation passed 24 Python fixture tests,
56 Perl checks and a real AGRN case reporting both start-gain and stop-gain;
the installed-module replay exactly matched the candidate. Both installers
recognize the separately recorded local patch. Existing annotations retain their
original producing profile. Full-cache/fork/chromosome and scientific
qualification remain pending. The user's later resume is active, and the
investigated HGNC_ID-only metadata differences are accepted. Earlier dated job
states below are historical observations.

Source, original backup, exact diff and restore command:
`/stornext/snfs170/next-gen/scratch/yangyxt/projects/PriVA_upgrade/scripts/utrannotator_all_effects_20261024_v1/README.md`.
Closed runtime/test evidence and logs:
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/caches/runtime/utrannotator-all-effects-20261024-v1/`.
<!-- utrannotator-local-reporting-20261024-end -->


<!-- installer-current-20261022-start -->
**Historical installer checkpoint — 24 October 2026, 9:56 AM CDT.**

`AlphaMissense annotation/motif/AF + ClinVar preserved → domain-statistic OOM repair → exact validation → new original/Nextflow literal acceptance`

The **71,697,556-row, 25-contig VEP116 AlphaMissense annotation, full-input motif
and AF stages are complete** and preserved. Original job **1522950 failed at
6:30 AM CDT**: AlphaMissense domain statistics exhausted its 96-GiB allocation
(step.0 OUT_OF_MEMORY0:125; parent FAILED1:0). ClinVar VEP, sort/index, gnomAD
and AA/splice-statistic generation had finished; latest-code reuse remains
subject to existing producer/profile gates. Domain partials are not ready.
Nextflow and both latest-source literal acceptance runs have not started.
The chromosome-checkpoint/ordered-merge helper is implemented in both
installers; its scientific collector is unchanged. All **42 fixture cases** pass
(native1514363.2759, zero skips). **Array1524246 chr21 is complete** (4m59s, child peak0.65GiB); independent
comparison1524249 exactly matches all208 genes against the unchanged collector.
All25 chromosome collectors and verification1524268 completed0:0.
Ordered merge1524277 exhausted128GiB after32m34s. All25 completed
checkpoints are preserved and reverified; source-identical merge1524279 was
submitted on ordinary medium1CPU/176GiB, without repeating chromosome collection. Full reference1524250 failed
at176GiB/OOM. Ordinary-maximum retry1524276 ended NODE_FAIL when its node
became unresponsive; no reference pickle exists and OOM is not established.
Fresh1524278 retains the unmodified22-worker code/profile on bigmem24CPU/256GiB
for memory headroom beyond ordinary-node capacity; last observed RSS was184.54GiB. Full-model comparison and both latest
literal installers remain pending; no VEP annotation is repeated.
The first pilot array1524243 failed before collection because `/usr/bin/time`
was absent; the corrected wrapper measures its child with Python instead, retaining
the failed launch evidence. The frozen failed run and prior source-only success-adapter/v8
fixtures remain history and cannot certify this failed installation.

Target: **10 CPUs per chromosome, up to 12 simultaneous chromosome jobs**, subject
to ordinary-queue availability and aggregate user limits. Eight configured VEP116.3
forks account for its possible nine children plus parent. ClinVar chr1–8 have
verified preserved serial outputs; the dated inventory does not certify later
chromosomes. The fork8/chromosome candidate remains **disabled** because real
exact comparisons found HGNC_ID and UTRAnnotator differences. Existing bytes,
failed comparisons and the original producing profile are preserved.

Both source checkouts now check isolated `diptest`/`pyarrow` imports before
annotation, reject aliased ClinVar AF destinations before writes, and require
configured assembly/VEP-major provenance before ClinVar warm reuse or repair.
Incompatible managed annotations and their overwrite members are retained under
`<cache_root>/.installer-retained/clinvar-annotations/` before the existing rebuild.
Incompatible or unverifiable explicitly supplied annotations are preserved and
reported. Recovery retries accept node-dependent device-number changes only;
all other recorded identity fields and immediate header/move checks stay strict.
Historical identities are retained separately from retry observations. Genuine
bundled snapshots still use existence/nonzero checks only.

Validation: 116 bootstrap fixtures plus real isolated environment imports, 64 AF
alias fixtures, 224 header/recovery fixtures, and **126 affected retry/gene-refresh
fixtures** passed under Slurm. These are overlapping engineering suites, not a
single independent-test total. The final retry step1514363.2705 completed0:0 at
2:51:29 PM CDT. Bounded real AlphaMissense and ClinVar pilot header checks passed;
no full-model, public-data or family ranking/ACMG acceptance is claimed. CADD
Rust remains suspended; scientific-profile qualification is a separate gate.

Matching reports and detailed source/test evidence:
`/stornext/snfs170/next-gen/scratch/yangyxt/projects/PriVA_upgrade/{README.md,migration_plan.md,migration_plan.html}`
and `scripts/main-install-handoff-20261017/README.md` there. Cache payloads,
runtime products and all logs are under
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/`.
Older dated checkpoints below are historical and superseded here.
<!-- installer-current-20261022-end -->
<!-- mapped-artifact-publication-20261019-start -->
**Historical installer checkpoint — 19 October 2026, 19:05 UTC:** both latest full
`main_install` acceptance runs remain pending. AlphaMissense1515341 is RUNNING;
checker1515343 waits on its dependency. The final AlphaMissense publication event
has not been observed. The running frozen producer is unchanged.

**Completed in both current installers:** motif and PM1 builds now write into
private staging, create the existing mapping receipt, then publish the pair with
interruption recovery. Scientific scripts, input selection and calculation order
are unchanged. Unknown supplied outputs without valid mapping receipts remain
refused; prior outputs and failed staging files are preserved.

`Staged output → existing mapping receipt → recoverable pair publication → reuse`

Slurm1514363.1478/1479 passed **25 focused/routing tests per checkout**;
1514363.1482/1483 passed **9 affected legacy tests per checkout**. These are tiny
fixture producers, not full-cache or scientific acceptance. Initial setgid-stage
failures are retained. The separately reviewed frozen handoff v2 passed **22
recovery tests and61 controller tests** in1514363.1480/1481, both COMPLETED0:0.
It binds recovery of a legacy interrupted motif pair to a pre-publication absence
receipt and preserves the interrupted bytes while selecting one fresh managed
output path. No stop, recovery capture or replacement main submission has run.

Sources and handoff instructions:
`/stornext/snfs170/next-gen/scratch/yangyxt/projects/PriVA_upgrade/scripts/main-install-handoff-20261017/README.md`.
Closed receipts/logs under
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/` match
`deployment/mapped-artifact-publication-20261019-v1.json` and
`{deployment,logs}/motif-recovery-{guards,controller}-20261019-v2.{json,log}`.
CADD shared-feature, public-data and family qualification remain pending.
<!-- mapped-artifact-publication-20261019-end -->

<!-- cadd-native-runtime-20261019-start -->
**Native CADD integration — 19 October 2026:** both working trees replace the
legacy `CADD.sh` caller with one shared official VEP run, the original ordered
publisher-prescore lookup, complete feature enrichment for true misses, and the
unchanged four-stage Python raw/PHRED calculation. Exact duplicate rows, header
selection and final CADD sorting are retained. The private model view contains
CSQ only; the full PriVA VEP artifact is preserved. Whole-genome wrapper joins
use streams and temporary SQLite storage under `tmp_dir`.

Native both-stage installation now publishes a resource manifest and an
**unqualified** candidate at
`<cache_root>/runtime/cadd-1.7.2/profiles/<assembly>.vep<release>.<digest>.installed.json`.
Installation does not activate it. Analysis requires an explicitly qualified
`cadd_profile` (or `--cadd-profile`); missing selection, prescore-policy or
resource evidence stops execution. CADD Rust remains suspended.

The native/core changes passed **132 affected fixture checks** in
**1514363.1338**, plus **2 source-drift checks** in **1514363.1342**.
These used tiny fake assets/backends. Shared workflow fixtures passed **17/17
per checkout** in **1514363.1341**, then the final explicit-profile guard passed
**2/2 per checkout** in **1514363.1343**. Independent final source review found
no remaining integration blocker. The initial shell fixture run retained two
invalid synthetic-header failures; corrected fixtures passed. Logs
`cadd-native-deployment-20261019-*.log` and
`cadd-shared-workflow-tests-20261019-*.log` are under
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/logs/`.
Source and matching plans are under
`/stornext/snfs170/next-gen/scratch/yangyxt/projects/PriVA_upgrade/`.

Both latest full `main_install` runs remain pending. AlphaMissense **1515341**
is still running; checker **1515343** waits on it. CADD feature diagnostic
**1515599** closed COMPLETED (0:0) at 17:39:41 UTC on 19 October, but its
shared-feature compatibility gate failed. Both scenarios produced complete
151-column tables with exact row/FeatureID order. The merged scenario selected
`NM_007294.4`, which is missing from the pinned peptide FASTA and GTF; 96 ESM
entries retain the original missing-peptide fallback. Execution success does
not qualify these features. No transcript replacement or gate relaxation was
applied; scientific/public/family acceptance remains pending. The completed
**80/80 profiling comparisons** remain retained without rerunning.

**CLI evidence scope clarified, 19 October:** the direct original-Python PHRED
CLI already passed **10/10 edge cases** on 18 October, including NaN for both
assemblies: exact exit codes, stdout and stderr. See the [closed CLI receipt](/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/validation/cadd-native-runtime-20261018/cli-parity.ready.json).
The NaN exit 0 to exit 1 mismatch remains a failed **suspended `4194be3` artifact wrapper**
gate, separate from the experimental **PID-only cancellation** failure. This
does not qualify the complete artifact interface or accept malformed score tables;
strict artifact checks and scientific/public/family gates remain unchanged.
<!-- cadd-native-runtime-20261019-end -->

<!-- cadd-refseq-acquisition-20261019-start -->
**RefSeq resource candidate — 19 October:** the separate
`scripts/cadd_refseq_resources.py plan|acquire|check --config CONFIG` helper
acquires four pinned NCBI originals matching merged VEP116 source
`GCF_000001405.40-RS_2025_08`. It passed23 fixture tests per checkout.
Real acquisition/check1515736 completed0:0 at18:10:18UTC:242,277,346bytes, with
publisher MD5/length and local SHA256 checks. Files are under
`cache_root/references/RefSeq/GRCh38/GCF_000001405.40-RS_2025_08/`;
`cache_resources.json` records this separate experimental snapshot.
Exact one-transcript validation1514363.1443 passed10 named checks for
`NM_007294.4`→`NP_009225.1`, full sequences,23 exons/22 CDS segments and
all96 diagnostic positions, with zero mismatches. The validator passed16 focused
tests, including installed-Perl extraction. Whole-RefSeq and model compatibility
remain unqualified; the complete proof is recorded in the matching project plan.
This does not change `main_install`, the old288-row diagnostic, its missingness,
the chosen transcript policy or the disabled CADD qualification gates.
<!-- cadd-refseq-acquisition-20261019-end -->


<!-- nextflow-main-driver-20261019-start -->
**Nextflow driver automation — 19 October 2026, 06:57 UTC:** the Nextflow
checkout's `main_install` now installs/reuses pinned Nextflow25.10.4 and
Temurin21.0.12.1+1 at `cache_root/runtime/priva-nextflow-25.10.4`, after
manual/configuration/bundled-snapshot preflight and before reference acquisition.
Original PriVA has no new Nextflow requirement. The modular bootstrap retains
its upstream pins and launcher, resumes owned downloads, retains failed bytes,
and publishes only a verified completed runtime. No extra YAML key is added.

Slurm1514363.1286 passed152 integration/affected tests in94.265s. After adding
snapshot-failure coverage, the six dispatch tests passed again in1514363.1288;
16 bootstrap cases passed with tiny fake download/JVM assets. These fixtures are
separate from the real test in1514363.1291: historical-runtime verification,
managed cold installation (11.907s), warm reuse (7.339s), and managed verification
all passed. Payload bytes/file identities were unchanged across warm checks;
source/configuration were unchanged. No workflow or scientific model was launched.
Both new suites are included in CI, but remote CI has not been rerun for this edit.

The real driver receipt is
`deployment/nextflow-driver-main-integration-real-20261019-v1.json` under
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade`.
The running AlphaMissense build is untouched. Full latest-source `main_install`
acceptance for both checkouts, shared-CADD qualification and family acceptance
remain pending; driver installation does not open those gates.
Root-alias correction (19 October 2026): the driver now normalizes the configured
root lexically, then resolves that existing root, matching the shared cache
resolver. YAML stays unchanged; linked runtime/staging children are still refused.
Nine focused dispatch tests and real canonical/alias warm reuse passed in1514363.1296
(7.384s /7.339s for the real calls), with source/configuration preserved. This is
reused driver verification, not a repeated cold install or full main acceptance.
Receipt: `deployment/nextflow-driver-root-alias-20261019-v1.json` under the same
data-side project. Earlier cold/helper evidence remains tied to its source.

<!-- nextflow-main-driver-20261019-end -->

<!-- installer-tmpdir-20261019-start -->
**Temporary-directory maintenance — 19 October 2026:** both CADD runtime
helpers select explicit `--tmp-root` / `PRIVA_CADD_TMP_ROOT`, then YAML `tmp_dir`,
then the legacy receipt-directory `tmp` fallback. Nextflow's offline driver
probes honor inherited `TMPDIR` over stale `NXF_TEMP`; persistent runtime state
stays under the managed runtime root.

CADD fixtures passed **48/48 per checkout** in Slurm 1514363.1301; bootstrap
fixtures passed **18/18** in Slurm 1514363.1302. Actual managed-driver warm reuse
passed in **7.202s** in Slurm 1514363.1303, preserving 262 payload file/link
identities and the configuration/source. Both prepared-config CADD resolvers selected YAML `tmp_dir`;
no model ran. Receipt:
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/deployment/installer-tmpdir-native-20261019-v1.json`.
Latest `main_install` acceptance for both checkouts and scientific/family
qualification remain pending.
<!-- installer-tmpdir-20261019-end -->

<!-- installation-current-closure-start -->
**Current checkpoint — 19 October 2026, 02:19 UTC:** installer1515341 remains
RUNNING; checker1515343 waits on its dependency. Fresh-host installation is deferred.
The full SpliceAI2/ClinVar qualification1515186 has **FAILED1:0** after VEP completed,
at its all-source allele/ID/ClinVar-field preservation audit. The first missing
record is a nonvariant (ClinVar836156, chr1:942451 T>., no_sequence_alteration);
a raw-only inventory found1,341 source-proven nonvariant skips matching the count
gap. Full audit1515478 now confirms exactly those1,341 omissions:4,593,423 retained
rows match, with zero unexpected additions/losses or retained-row order changes.
The original all-source preservation contract remains false and is retained. Postcheck1515187 separately
**FAILED1:0** at an environment-state hash guard before annotation review. Both
failures are retained; neither permits scientific activation or a relaxed gate.
Full identity/order audit1515478 COMPLETED0:0; no VEP/scoring was repeated. AF retry
1515465 closed COMPLETED0:0 with all 200,000 row/order/INFO/index checks passed.
The latest AF warm-reuse fix passed34 focused tests per checkout, then both latest
wrappers reused the completed200,000-row artifact with unchanged bytes/index/file
identity in1514363.1225 (5.501s /5.467s). It corrects header-ID matching without
changing annotation commands. Focused evidence is in DataLad9a99a67; the real warm
receipt is validation/alphamissense-af-warm-check-20261019-v1/receipt.json under
the data-side project. Full latest-code installer acceptance remains pending. A separately reviewed
recovery procedure passed13 new plus20 affected cases: only after genuine completed
AlphaMissense publication, close the obsolete checker/writer and verify the closed
resource before running both latest entrypoints. Old cancellation is not acceptance.
No job has been stopped or successor submitted.
The latest raw-input chromosome-map correction passed two portable real-tool tests
per checkout in1514363.1252; all34 gnomAD tests per checkout also passed on the same
source in1514363.1248. AF/cache-producing functions are unchanged, so no large
cache rebuild is needed for this correction. The conservative ClinVar source
fingerprint remains enforced. Latest annotation source is064a134f; full installer
acceptance is still pending.
Separate CADD enrichment diagnostic1515463 FAILED1:0 at01:59:37UTC on19October
during ESM missense at its fixed deadline; later feature stages were not reached.
This does not block installer acceptance or establish a numerical mismatch.
At that earlier checkpoint, native CADD components/assets were installed but the
caller still used legacy CADD.sh/Snakemake and its second-VEP110 route. The native
shared-VEP replacement is now implemented and undergoing integration validation,
as recorded above; full main_install and scientific acceptance remain pending.
Original CADD specifies Ensembl-only
per-gene selection; new oracle checks will test that reference behavior.
The current project README/MD/HTML record recovery evidence. Earlier dated running
checkpoints below are historical, not the current status.
<!-- installation-current-closure-end -->


<!-- independent-install-checks-start -->
**Parallel installer fixes — 18 October 2026:** fresh-install and interrupted-build fixes are implemented
and tested in both working checkouts. Full installation acceptance remains pending.

| Independent thread | Completed evidence | Scope |
| --- | --- | --- |
| Manual input readiness | Early gate after environment setup, before resource acquisition. Prompts fill only unset keys together. **131 distinct focused cases per checkout passed** before the later cache-error fixes. | Synthetic helper/PTY/local-HTTP cases. Actual prepared files separately passed metadata-only readiness without changing the config. |
| Fresh pext acquisition | Fixed Bash variable scope that could omit selected tissues. Real cold downloads obtained **four tracks per checkout**; all eight passed native-reader queries and offline reuse. | Actual resource acquisition; transfer interruption/retry was simulated. |
| InterPro failure handling | Empty/error results now fail; mappings publish atomically after successful serialization. Existing final mappings survive failed builds. **51 focused cases per checkout passed**, including exact frozen-parser pickle bytes/order/multiplicity and injected I/O failures. | Real helpers on synthetic XML/VCF; bundled reuse remains existence/nonzero only. No full publisher rebuild was used as acceptance. |
| GRCh37 SpliceVault readiness | Cold installation now uses the same usable-row/index check as warm reuse before updating config. **17 focused cases per checkout passed**, including real conversion and cold/warm equivalence. | Synthetic mapped/unmapped variants, reference and chain. Header-only failed candidates are retained but not registered; liftover/filter/order rules are unchanged. |
| InterPro consumer contract | Configured mapping reaches DAS and standalone PM1; DAS, PVS1 and splicing readers select raw/gzip by bytes. **24 focused cases per checkout passed**, plus real imports and CLI help. | Tiny real producer/reader fixtures; ordinary pickle reads, no extra bundled integrity scan. Mapping values/order/multiplicity and writer bytes are unchanged. |
| Chromosome error propagation | Unexpected worker errors now abort the parent before publication; missing required exon CSQ schema fails nonzero. **6 test methods / 29 subcases per checkout passed**. | Real indexed tiny VCFs and multiprocessing, differential against frozen v7. Normal record skips and successful values/order are unchanged. |
| HPO/MONDO recovery | HPO writes publish atomically. A newer registry with older scoped assertions triggers a rebuild instead of mixed-generation reuse. **8 tests per checkout passed**. | Tiny real helpers, failed writes/publication, restored 7/7 pairing and unchanged warm reuse. |
| Domain-cache completion | A small metadata-only receipt is written after the scientific command, all three selected outputs and all YAML updates succeed. Failed/unreceipted builds retry; warm reuse restores missing output keys and receipts exclude node-local device IDs. **11 tests per checkout passed**. | Real shell/receipt helper with a synthetic scientific writer; failures after each output, config-write failures, changed inputs, custom paths and warm reuse. Not a full DAS model test. Legacy unreceipted derived domains rebuild once. |
| Domain early-exit handling | DAS explicitly reports completion only after all original writes; its CLI and standalone PM1 regeneration caller reject incomplete execution. **12 tests per checkout passed**, plus **24 existing reader tests** and independent source review. | Actual producer/callers/shell with labelled scientific helper stubs; exact successful output bytes and operation order match the frozen function. Early exits preserve older results but cannot register them as a new completed build. Standalone PM1's separate warm-mtime reuse remains unchanged; main_install uses the earlier gated DAS stage. |
| Domain intermediate association | Selected VCF, assembly, extraction source and completed output must match before reuse. Atomic data/sidecar publication preserves prior data on failed writes. **24 tests per checkout passed**, with exact original pickle bytes/order on tiny real extraction. | Old unassociated or changed-source derived metadata regenerates from the selected VCF. This corrects stale-source reuse; no filtering thresholds or extraction rules change. Warm-read races and partial-write recovery are covered. |
| Config filename dispatch | Quoted config-file arguments through main/plugin calls and affected reads. **5 new tests per checkout passed**, plus **17 existing plugin-input tests** on the same production code. | Stubbed stage calls and a real selected UTR input read; no acquisition. Existing whitespace restrictions on managed storage roots are unchanged. |
| PEXT before first annotation | Moved the sole pext installation call before ClinVar and recorded its map dependency. **7 routing tests per checkout passed**, plus shell syntax and both 39-resource maps. | Acquisition/VEP stubs verify ready-before-consumer and failure abort. Existing HGSC tracks were already present; this does not certify older warm VCFs for every configured annotation column. ClinVar statistics do not consume PEXT. |
| GRCh37 reference/chain filenames | Shell quoting and CrossMap argv lists preserve spaces/apostrophes in accepted filenames. **18 tests per checkout passed**, including actual cold/warm indexed conversion. | Includes the earlier 17 plugin-input cases. Liftover/filter/row order is unchanged; managed storage-root whitespace remains unsupported. The initial PATH preflight failed before any test and is retained. |
| Latest-code acceptance | Successor controller passed **48 fixtures and independent review**. The fixes are mirrored across original and Nextflow checkouts. | Active **1515341** uses earlier frozen v7/v4 and excludes these fixes. The reviewed recovery route can reuse finalized AlphaMissense after controlled closure of the obsolete pair; both latest literal `main_install` entrypoints and their final resource/terminal checks remain required. No successor has been submitted. |

These test counts describe separate checks and must not be summed as one suite.
Closed evidence under `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade`:
`deployment/installer_fresh_recovery_20261018_v3.json`,
`deployment/manual-cache-real-readiness-20261018.json`,
`deployment/pext-cold-acquisition-20261018-v4.json`,
`deployment/interpro-atomic-tests-20261018-v2.json`,
`deployment/splicevault-cold-readiness-tests-20261018-v1.json`, and
`deployment/accepted-predecessor-tests-20261018-v3.json`,
`deployment/interpro-consumer-tests-20261018-v2.json`,
`deployment/chromosome-worker-failures-20261018-v3.json`,
`deployment/mondo-hpo-recovery-tests-20261018-v3.json`, and
`deployment/das-install-completion-20261018-v3.json`,
`deployment/domain-metadata-cache-tests-20261018-v1.json`, and
`deployment/installer-config-space-tests-20261018-v2.json`,
`deployment/das-early-completion-probe-20261018-v1.json`, and
`deployment/das-completion-contract-tests-20261018-v1.json`, and
`deployment/interpro-consumer-tests-20261018-v3.json`, with corresponding logs.
Initial observed defects and earlier fixture failures remain retained.
Additional closed order/path evidence: `deployment/pext-main-order-tests-20261018-v1.json`
and `deployment/splicevault-selected-path-tests-20261018-v2.json`, with matching
logs; DataLad853ad4c /117cb2c. These are scoped implementation checks, not full
installation or scientific acceptance.


Source remains in both checkouts' `scripts/` and `tests/`; maintained check harnesses
are under `/stornext/snfs170/next-gen/scratch/yangyxt/projects/PriVA_upgrade/scripts/`.
The data-side `validation/interpro-{output-probe,atomic-tests}-20261018-v*/**` and
`tmp/splicevault-cold-reuse-20261018-v1/**` contain synthetic diagnostic fixtures,
not production caches. Scientific, public-data and family gates remain unchanged.

Live checkpoint **18 October 2026, 17:42 UTC**: **1515341 RUNNING** (lowspace,
24 CPUs / 96 GiB) continues the frozen original installer’s VEP AlphaMissense stage;
its intermediate output was 11,717,059,342 bytes at 17:11:09 UTC on the owning node.
This is progress metadata, not completed-cache evidence. **1515186 RUNNING** is the
independent SpliceAI2 qualification. **1515343 / 1515187 PENDING** on dependencies.
The running snapshot excludes these new fixes; latest-source full acceptance remains pending.
The separate DAS early-return false-success defect is fixed and covered by the
12-case completion checks above; full-model and latest-source installation acceptance remain pending.
<!-- independent-install-checks-end -->

<!-- cpu-scaling-status-start -->
**Sequential installer acceptance — 18 October 2026, 15:29 UTC:** **1515341 RUNNING** on
`lowspace/c81n-40`, 24 CPUs / 96 GiB, explicit VEP fork22. Literal updated original
`main_install` (frozen v7) began at **15:28:08 UTC**; Nextflow (frozen v4) follows
only on success in the same allocation, then bounded VEP/pext resource checks.
**1515343 PENDING** on `afterany:1515341` requires exact closed allocation, batch,
results and configs before accepting installation. Full software acceptance is
pending; scientific/public/family qualification remains separate.

```text
1515341: original v7 [RUNNING] → Nextflow v4 [WAITING FOR SUCCESS]
  → bounded VEP/pext checks → 1515343 [PENDING; afterany:1515341]
  → exact terminal job/batch/results/config checks → acceptance [PENDING]
```

**Measured CPU choice:** lowspace benchmark **1515332 COMPLETED 0:0** on 20,000
identical rows: 8 CPUs/fork8 **393.477 s**, 24 CPUs/fork22 **190.790 s** (**2.062×**).
Pre-sort and post-sort comparisons had zero differences; headers matched except
execution time/fork metadata. This is sample-only evidence, not full-table or
family acceptance. Bigmem trial1515330 was cancelled without a completed scaling comparison;
use bigmem only for genuine high-memory need when other suitable queues are
unavailable. Trial1515331 failed before VEP because GNU time was missing; the
portable timer correction is retained. Independent1515186 running and1515187
pending remain protected and unchanged.

**Controlled replacement history:** old1515311 →1515309 →1515293 →1515291 →1515289
were cancelled leaf-first at **15:20:24–15:20:38 UTC**; all four dependents never
started, and original partial outputs remain retained. Evidence:
`deployment/installer_cpu_scaling_stop_20261018.json`. Held1515339 was cancelled
before starting after its resource-parser guard failed; evidence:
`deployment/sequential_cpu_guard_held_stop_20261018.json`. The isolated v2 driver
passed **53 fixtures** and independent review; its initial fixture-only failure
and v1's 51 passing fixtures remain historical evidence. No scientific gate changed.

**Current controls:** `main-install-handoff-20261017/attempt_cpu_v2/` holds the
plan, submission, runtime configs and logs; future `acceptance.closed.json` is
pending. Plan SHA256 `7184284feebb5e934c76b5a573b99730c2da5dfc0965dc93b67098f5dab4fac5`;
submission `22ad5703e50871012f52a84f9f7bf63e72db848d8dd9aa4f071db6926ecf627b`.
Checker receipt: `deployment/sequential_cpu_terminal_submission_20261018_v2.json`;
benchmark receipt: `deployment/vep_scaling_lowspace_20261018.json`. These paths are
relative to `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade`;
maintained helper source is under
`/stornext/snfs170/next-gen/scratch/yangyxt/projects/PriVA_upgrade/scripts/main-install-handoff-20261017/cpu_scaling_candidate_20261018_v2/`.
<!-- cpu-scaling-status-end -->

<a id="installer-retry-20261018"></a>

**Historical 12:30 UTC chain — superseded by the 15:20 UTC controlled stop:**
original1515289 and never-started Nextflow1515291 / closure1515293 are CANCELLED.
Full software installation acceptance awaits successful jobs, results and final
config/resource checks. Scientific/public/family qualification remains separate.
Frozen v5/v2 source passed v8: **875 PriVA / 905 Nextflow fixtures**, zero failures/errors/skips,
with source and harness unchanged; the original/Nextflow/closure retry helpers
passed **37 / 53 / 43 fixtures**. These checks are not completed installation.

**Historical source boundary:** the cancelled1515289/1515291 chain selected frozen v5/v2;
later ClinGen acquisition/retry and four configured mechanism source URL fixes are
complete, tested and independently reviewed; they were outside that cancelled run and are included in current v7/v4 sources.
Focused v3 passed **57/57 unittest + 14/14 existing pytest cases per checkout**,
with zero failures/errors/skips and source unchanged. A separate check passed on a
**250-line, 41,334-byte nonbundled publisher dosage prefix** and a tiny ClinGen
parity fixture using the original and current source. Full cache and family
acceptance remain pending. URL settings apply only to actual acquisition or
regeneration; genuine bundles retain existence and nonzero size checks. No
scientific policy changed. The corrected warm-cache rerun now runs as1515341 using v7/v4 sources. Closed evidence: [focused v3](/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/deployment/clingen-source-urls-1514363.856-20261018-v3.json)
and [bounded format/parity](/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/deployment/clingen-source-parity-1514363.857-20261018.json).

**Earlier broad checks:** combined v11 passed **923 PriVA / 953 Nextflow tests**,
zero failures/errors/skips, with source and harness unchanged. The subsequent,
isolated CADD receipt portability fix passed **70/70 focused tests per checkout**
and independent review; that final fix is outside v11. Only the mount-local device
ID is excluded from persisted reference identities; inode, size and timestamps
remain strict. Retained v9 fixture failures and their 41/41 focused correction
remain recorded. A bounded live eRepo fetch passed the actual validator and original
pandas conversion for 20 complete rows; this is not full-download validation.
The sequential follow-up driver passed 38 fixtures and independent review.
**Final revision replacement:** the 39-resource audit found that CADD_install could
ignore a supplied custom `cadd_base_dir`. Both installers now preserve that selection,
retain invalid supplied paths, and keep the legacy alternate-directory route.
Final CADD tests passed **42/42 per checkout**, plus **150 related cases per checkout**
and shell syntax checks on unchanged production source. Independent review passed.
The initial fixture omitted two required config fields; its failed result is retained.
Evidence: `deployment/cadd-selected-source-20261018-v[12].json` and matching logs.
Jobs **1515303 and1515305 were cancelled before start**
so they cannot accept the older code. No main_install attempt or result was created.
The old five-job chain, including queued1515309/1515311, was cancelled at15:20 UTC. Its retained controls remain under `main-install-handoff-20261017/attempt_revision_v2/`. Current1515341 reuses the same frozen v7/v4 snapshots; checker1515343 enforces closed success. Full installation acceptance remains pending. Cancellation evidence:
`deployment/cadd-source-selection-queued-replacement-20261018.json`.

**CADD candidate resources installed:** both installers now acquire or reuse the
selected-release peptide/GTF files under `cache_root`; the cache map covers 39
resource groups. On this host, step **1514363.884** adopted the verified release116
raw/plain files from closed job1515153 into
`<cache_root>/references/Ensembl116/GRCh38/`, then passed offline verification.
Original files and source receipts remained unchanged. Closed evidence is
`deployment/cadd-selected-transcripts-{tests-20261018-v3,adoption-20261018}.json`
and matching data-side logs. Model/package pins are unchanged. This completes the
resource installation route; shared-VEP feature/scoring integration remains disabled
until its separate scientific qualification. Final warm-cache main_install
acceptance awaits current1515341 and terminal checker1515343.



**Installation tools and CADD audit — 18 October 2026 UTC:** SpliceAI2
TSV conversion now has its own bundled, standard-library-only Rust source,
`scripts/spliceai2_converter.rs`. The Python installation driver builds/caches
that small tool independently; it does not build or require the production
`priva` binary, Cargo or Polars. Existing validated cache receipts remain reusable
without a compiler. The obsolete `MaveDB_deployment` function and registry route
are removed; active `MaveDB_install`, bundled data and the approved GRCh37 omission
remain. Focused converter verification is recorded in the matching project plan.

**Converter verification:** Slurm 1514363.651 passed **128 focused tests per
checkout**, zero failures/errors/skips, with tested source hashes unchanged.
The standalone tool also produced byte-identical VCF output and JSON diagnostics
for **129 actual publisher rows** (64 SNVs, 65 indels including the negative-strand
GG→AA MNV) against the frozen converter. This is bounded conversion evidence,
not SpliceAI model, full-cache, public-genome or family acceptance. Existing cache
receipts are reused; no complete cache was reconverted. Fresh compiler acquisition
is fixture-tested only. Earlier interrupted/environment-failed attempts are retained.

The full AlphaMissense **variant prescores are external Google DeepMind releases**,
not Git bundles. Original `config.yaml` selected
`/paedyl01/disk1/yangyxt/public_data/VEP_plugins_caches/AlphaMissense/AlphaMissense_hg19.tsv.gz`
and matching derived files there. The original installer already had public
hg19/hg38 download URLs. HG19 acquisition job1514583 retained the 622,293,310-byte
publisher file (MD5 `3f14fba08c60b09a90a2ae14b7802e4c`); GRCh38 acquisition1515167,
verified by1515172, retained 642,961,469 bytes (MD5 `9fd167735f16a1b87da6eb3e4c25fcb5`).
The former `inputs/AlphaMissense/...` directory is historical acquisition storage,
not the installer destination. Both installers now use selected-assembly acquisition,
validated external reuse/index repair and configured derived-output paths beneath
`<cache_root>/VEP_plugin_cache/AlphaMissense/`. Invalid supplied originals are preserved.
The actual Git bundle is the different 1,500,483-byte per-gene summary
`data/alphamissense/alphamissense_mean_score.tsv`; reuse requires only existence/nonzero.
[Publisher hg19](https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg19.tsv.gz)
· [publisher hg38](https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz).

**Current-host acceptance:** fresh-host installation is deferred. Verify dependency
order and both latest installer entrypoints using prepared manual files and reusable
completed caches. Code changes require affected tests to rerun and affected derived
caches to rebuild; earlier results certify only their recorded source revision.

Fresh CADD runtime installation now defaults to native original-Python environments;
Apptainer is optional for an explicitly selected existing reference installation.
Job **1515248 completed 0:0**: 46 core packages installed, and both 18-row completed
annotation fixtures matched every original stage exactly (features, raw, hierarchy,
PHRED), against both the pinned image and archived original outputs. Ten additional
native/reference CLI edge comparisons matched stdout/stderr/exit exactly, including
the original NaN malformed-output behavior; this does not make that artifact valid.
Native feature dependency job **1515250 failed** on a seqiolib/TensorFlow package
constraint. Recovery **1515253 completed0:0**: all four native environments are
installed and passed offline imports, including MMSplice. Both assembly receipts
and publisher/member asset verification are complete (1h33m41s, peak RSS634536KiB).
Model-runtime smoke **1515447 completed0:0** at21:42:27UTC: all five ESM models,
RegSeq and MMSplice ran on synthetic sequences with finite correctly shaped outputs,
exact within-process repeats and unchanged source/model identities. This is model
loading/CPU execution, not shared-VEP feature generation or independent scientific
reference parity. Closed evidence is saved in DataLad **b19a477** under the data-side
project's `validation/cadd-model-smoke-20261018-v1/` and matching logs.
The reference distribution versions remain pinned.
Existing 80/80 profiling measurements were not rerun; the historical experimental NaN-validation
and PID-only cancellation gates remain separate. The old `Calculate_CADD` caller
used `CADD.sh -m`, separate Snakemake Conda environments and VEP110. The current
working-tree replacement uses the guarded native/shared-VEP route documented
above. Its installation and scientific acceptance remain separate; installing
the native environments alone does not establish a qualified raw-input run.
Public/family acceptance remains pending. CADD Rust remains suspended. See the
[shared-VEP/native-Python design and remaining gates](../../projects/PriVA_upgrade/migration_plan.html#lean-cadd-20261018).

**Resource lifecycle map:** [`cache_resources.json`](cache_resources.json) inventories
40 resource groups: the 39 shared with original PriVA plus the pinned Nextflow/Java
runtime installed by `Nextflow_runtime_install` under `cache_root/runtime`. Entries record
source URLs/version templates, selected pins, bundled paths, dependencies and
installation/refresh entrypoints. `maintained` means retain a refresh route, not a
promise of constant upstream updates; unverified cadence is recorded explicitly.
Snapshots reuse bundled data when supplied, while licensed/large external
snapshots stay external. Derived caches name the inputs that can require a rebuild.
Only the new ClinVar gene-stat age rule is automatically enabled by this map;
existing external refresh behavior is unchanged. Version templates are not proof
that arbitrary releases exist or are scientifically compatible.

```bash
python scripts/cache_resources.py list --refreshable
python scripts/cache_resources.py show --resource clinvar_gene_statistics
python scripts/cache_resources.py refresh-plan --resource alphamissense_prescores
python scripts/cache_resources.py validate --check-routes
```

These inventory commands never download, modify or validate cache payloads. The
refresh plan lists dependencies for review, not a scheduler or an execution order.
The ClinVar date/source policy is read by `scripts/clinvar_gene_refresh.py`.
[NCBI's download documentation](https://www.ncbi.nlm.nih.gov/clinvar/docs/downloads/)
and [Ensembl's release documentation](https://static.ensembl.org/info/docs/tools/vep/script/vep_download.html)
provide the recorded upstream feed/release evidence; other entries retain their
configured publisher endpoints and explicitly unverified update cadence.

**Earlier verification — 17 October 2026:** the final resource-map/installer
regression passed **562 PriVA / 592 Nextflow cases** under Slurm1514363, with zero
failures/errors/skips and unchanged source/harness hashes. At that checkpoint the maps
matched and validated all38 groups (15 snapshots,17 maintained distributions,6 derived).
These are software/fixture checks. A separate bounded live NCBI check verified
ClinVar `fileDate=2026-10-11` and stable publisher MD5s for both assemblies; it did
not acquire a full payload or regenerate gene statistics. The detector accepts
both ISO and compact dates. The initial parser failure and the broader run's one
obsolete GoFCards fixture per checkout are retained; the fixture now mocks the
new wrapper while keeping its intended failure assertions. No tolerances or
scientific gates were relaxed. Final controller fixtures:36 passed; the actual
SpliceAI2 guard probe correctly refused the unqualified GRCh38 route.

**Closed mechanism candidate:** job1515158 completed0:0 at
2026-10-17T22:20:45Z. All7stage receipts,9 frozen scientific sources from09223ef,
and39 receipted output hashes were verified. HPO/MONDO, mechanism, both-assembly
GoFCards/ClinVar and integrated HPO caches are built as candidates. The original
baseline is unchanged; production activation and scientific acceptance remain
false. These completed candidates do not remove the SpliceAI2 guard or establish
public/BH14233 acceptance.

**Earlier integration check — 18 October 2026:** full v5 passed
**843 PriVA / 873 Nextflow cases**, zero failures/errors/skips, with tested source
and harness unchanged. Receipt: data-side
`deployment/installer-policy-20261018-autonomy-v5.json`.
Previous v4 also passed805 /835 cases; it predates the final upgrade/path fixes.
Retained v1 collected 735 / 765 cases with four failures and nine skips per
checkout; v2 was interrupted in a compiler fixture; v3 collected 804 / 834 with
one outdated UTR fixture failure per checkout and zero errors/skips.
Focused suites passed 128 AlphaMissense/path cases, 44 native CADD contracts,
60 plugin v2 cases and 44 LoFtee/UTR v4 cases per checkout. These fixture checks
remain separate from scientific/public/family acceptance. No CI run is claimed.

**Single-command installation acceptance — GRCh38/BH14233, 18 October 2026:**
**Historical cancelled retry controls (superseded at 15:20 UTC):** original job **1515289** started at **10:24:08 UTC** from the 340-file `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-install-snapshots/main_install_20261018_v5` snapshot, with retained controls and closed log under data-side `main-install-handoff-20261017/attempt_v5/`. Nextflow job **1515291** uses the 430-file `main_install_nextflow_20261018_v2` snapshot in the same external snapshot directory and controls under `attempt_nextflow_v2/`; it cannot enter `main_install` before exact original allocation/batch/result/final-config success. Closure job **1515293** then requires both exact installers to close successfully before the unchanged bounded VEP116/pext checks. Its future result is `main-install-handoff-20261017/closure_v2/readiness.json`; no actual readiness result exists yet. Final integration evidence is `deployment/installer-policy-20261018-autonomy-v8.json` under `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade`. Prior sources, configs, logs and incomplete outputs are retained.

<a id="installer-stop-20261018"></a>

**Historical controlled stop:** original job **1515274 was CANCELLED at 09:43:53 UTC, 18 October 2026** (batch closed at **09:43:54 UTC**), during AlphaMissense VEP116 annotation and before the gnomAD stage. The actual `main_install` had begun at **08:25:20 UTC**, after staging, and reused the existing VEP116 environment. The stop followed a source audit of gnomAD error propagation; no partial-publication failure was observed. Its 339-file source remains frozen at `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-install-snapshots/main_install_20261018_v4`; the plan, retained runtime config and closed log remain under data-side `main-install-handoff-20261017/attempt_v4/`. The corrected guards subsequently passed final v8 fixtures; the replacement chain is recorded above. This cancelled attempt supplies no completed installation or scientific/public/family acceptance.

Closed cancellation evidence: [`deployment/installer-gnomad-failure-guard-stop-20261018.json`](/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/deployment/installer-gnomad-failure-guard-stop-20261018.json), SHA256 `53ed262f503f47b59458beac91a3f44453c297f420ac2961b5d118e45f5e2138`.

**Historical Nextflow submission — cancelled before start:** job **1515277**, submitted at **09:02:15 UTC, 18 October 2026**, was **CANCELLED at 09:43:53 UTC** while pending on `afterok:1515274`; its `main_install` never started. Its 429-file source remains frozen at `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-install-snapshots/main_install_nextflow_20261018`; retained controls remain under data-side `main-install-handoff-20261017/attempt_nextflow_v1/`. The wrapper's **37 focused fixtures** and independent review remain closed evidence. Its design preserves scientific selections, rebases genuine bundled/code defaults and uses isolated temporary/output paths, with no wrapper payload staging or cancellation. This cancelled attempt supplies no real Nextflow installation acceptance; the later1515291 replacement was also cancelled before start; current1515341 is recorded above.

**Staging recovery:** the preceding job1515270 failed before main installation on a 600-second copy timeout. Recovery byte-verified 151 completed copies, retained the incomplete CADD copy and broken provenance link separately, and used verified immutable prescore hardlinks with independent indexes. Six frozen README sync amendments were preserved and the recorded bodies restored; new snapshots avoid automatic README mutation. The final controller passes **56 fixtures**; no source-hash or scientific gate was relaxed.

**Retained failed submission:** job1515267 failed1:0 before entering the controller: Slurm executes a copied batch script, so its relative lookup could not find `handoff.py`. No `main_install`, cache staging or superseded-job cancellation occurred. The failed source, log and submission are retained in `attempt_v2`; the explicit engine-path fix now passes46 controller/spool fixtures before attempt_v3. Core installer v5 remains843 /873 passed.

Slurm **1515262 was cancelled before launch** after its deliberate hold; its
06:59:37 UTC submission and source snapshot remain preserved. No installer or
cache staging ran in that attempt.
The AlphaMissense release/assembly reuse guard now passes 76 focused cases per
checkout, including a read-only check of the completed real VEP116 pilot header.
Mapping-transition routing passes 82 focused cases per checkout: intact older
motif/PM1 artifacts remain preserved while the selected mapping gets a distinct
managed path. Custom-path fixes pass 115 cases per checkout: an absolute reverse
chain is adopted under cache_root and honored by both consumers; the existing
repeat filename setting rejects path syntax clearly. Bundled defaults retain
existence/nonzero reuse. Full integration v5 passed against unchanged sources.
The previous v4 evidence predates these corrections. The prior frozen
source is `source_snapshots/main_install_20261018` (336 files). At the earlier held checkpoint, the literal installer
had not started. Nine closed producer receipts and
200 staging operations passed read-only preflight: 169 copies, 28 hardlinks,
three hardlink trees and no symlinks. The unproven unversioned GRCh38 CDS copy was
removed; installation will fetch the official release116 source. The controller
passes 46 unique fixtures while retaining producer and file-identity guards.
The installation-only ClinVar route preserves the production scientific guard.
CADD prerequisite1515253 closed successfully (0:0). The retained handoff design
verifies its receipts before staging, exact superseded-installer cancellation
and actual execution. The later1515289 chain was cancelled at15:20 UTC; current1515341 runs the updated original/Nextflow stages. Independent scientific
qualification and unrelated jobs stay protected. No one-command installation
success is claimed.

```text
1515341 original v7 [RUNNING] → Nextflow v4 [WAITING FOR SUCCESS]
  → bounded VEP/pext resource checks
  → 1515343 terminal checker [PENDING; afterany:1515341]
  → full software installation acceptance [INCOMPLETE]
  → scientific/public/family qualification remains separate

Latest ClinVar fileDate ≥ bundled statistic UTC mtime date + 3 calendar months?
  no  → reuse the nonempty bundled statistic
  yes → pin inputs → original VCF/AA/gene build order → select managed result
```

The ClinVar rule is an explicitly requested **metadata-age exception** to bundled
reuse. It reads the file's modification time, without content/schema/checksum
validation of the bundled statistic, and preserves the original file. Month ends
are clamped (for example, 31 January + three months is 30 April). The latest data
date comes from the publisher VCF header, not the local download time. New source
identities and the pinned MANE input are recorded separately; an explicit
incompatible VCF must not be silently replaced. A successful managed result is
selected by its full path so the actual workflow consumes it. Network/rebuild
failures must be reported rather than treated as successful freshness checks.

**ClinVar AA/splice producer compatibility — 18 October 2026:** direct reuse now
requires a completed pair association with the current producer and VCF/index/
output metadata. Interrupted generation, publication or YAML writes cannot mark
the pair ready. Existing outer refresh receipts already bind the scientific
producer and outputs and keep their stronger reuse check. This does not validate
bundled snapshots or change scientific calculations. Both checkouts passed
33 focused tests each under Slurm1514363.1108, including a two-record indexed-VCF
byte/order comparison with frozen v7 and failure/retry cases. Independent review
passed; DataLad caa025f retains the closed evidence. Full latest-source installer
acceptance remains pending.

AlphaMissense's early-return issue concerned an **external downloaded prescore**.
The current fix validates external reuse and repairs a missing/broken index only
when the data are valid. Genuine bundled snapshots remain existence/nonzero-only.
Conda/Mamba/Micromamba is an existing documented prerequisite, already satisfied
on HGSC; it is not an outstanding installation gap.

**Current correction — 17 October 2026:** CADD Rust work and the `cadd-rs`
integration are suspended. Native runtime installation provides the original
upstream Python stages through raw and PHRED scores; production integration into
the legacy workflow remains gated. The old repository,
binaries, receipts and benchmarks remain historical evidence; their earlier
integration requirement is superseded. Shared official VEP annotation plus complete
CADD feature enrichment remains the design. Portable first-install paths and release116 mapping selection are implemented;
verification is recorded below. Scientific and public/family acceptance remain pending.

**Earlier installer checkpoint — 17 October 2026:** Slurm **1515211 completed0:0**:
**526 PriVA / 549 PriVA-nextflow tests passed**, zero failures/errors/skips, with
source and harness hashes unchanged. Both installers reuse canonical bundled data
using existence and nonzero size only; external inputs keep their existing routes.
The reviewed ClinGen eRepo edge case is fixed, including direct calls, without
changing existing external-source precedence. AlphaMissense external reuse was unchanged
at that checkpoint; its current acquisition/reuse route is described above.

The retained first combined run1515209 had five failures and one error per checkout:
two path expectations and four GoFCards/basic-configuration fixtures still described
the previous install setup. Their route assertions were retained with the correct
selected paths and inputs. Separate negative/positive tests demonstrate the actual
ClinGen fix. These are installer contracts/fixtures; the independent18-row CADD
comparisons use real models, but public-data and family acceptance remain pending.

**Historical installation checkpoint — 17 October 2026:** installation and migration are **not complete**.

```text
VEP116 + both core caches [installed] → assembly-matched annotation
SpliceAI2 full cache [prepared]       → bounded real VEP check [passed]
                                     → derived caches [building] → scientific gates [pending]
CADD original-Python runtime [checked] → annotation resources [deploying]
                                     → shared-VEP feature qualification [pending]
```

- **Completed:** VEP116.3 and official release116 plugins; both merged caches; four explicitly approved obsolete113 payloads removed by successful recovery **1515116**. The earlier synthetic core checks establish cache loading/CSQ rows only. SpliceAI2 **1515009 completed 0:0** after all47 source checks and full conversion/index publication. The real SNV/indel/MNV check passes **96 scored transcript rows with zero exact score/gene mismatches**; hg19 remains on v1.
- **Building:** gnomAD GRCh37 **1514729**; CADD GRCh38 annotation assets **1515113**, followed by Python-only runtime/resource verification **1515208** (supersedes never-started1515170/1515129/1515160); AlphaMissense116 annotation/motif **1515130** after a96-record BRCA1 pilot passed. Corrected ClinVar pilot1515131 and join audit1515132 passed. Its full source-order chain is active: **1515136 VEP →1515137 population →1515138 AA/gene →1515139 Alpha domains/PM1 →1515140 AF**, with gnomAD/Alpha producer dependencies. The earlier wrong-LoFtee pilot remains failed qualification.
- **Installer fixes:** both checkouts expose `CADD_runtime_install CONFIG [annotate|calculate|both]`, called with `both` by `main_install`. The direct Python wrapper passed both 18-row real-model fixtures, exactly matching features, raw scores, hierarchy and final PHRED output against fresh pinned upstream Python and the archived pure-Python baseline. No cadd-rs binary is required. Installation does not enable the experimental scorer or qualify shared-VEP features.
- **Remaining:** finish and validate derived caches, matching GRCh38 plugin resources and mechanism inputs; qualify complete shared-VEP CADD features; complete PriVA Rust modules and public/family acceptance. hg19 MaveDB omission is approved. New BioMart/ClinVar/VEP116 candidates are not identical to unavailable original resources.

The SpliceAI2 adapter preserves the scored allele (including GG→AA), scores and
threshold provenance. Earlier bounded evidence remains:44 Rust tests,25 driver
fixtures,5 identity cases per checkout,12,000 actual indel rows and the overlapping
512-row regression; mocked candidate retrieval is distinct from the new real VEP
check. The first real test added `--check_ref` and returned intergenic annotations;
removing that extra test-only option restored expected features. The passing test
follows original PriVA options and independently verifies REF with `samtools faidx`.
Both attempts are retained. V2 scientific thresholds, public-data and family gates
remain pending. Complete preparation evidence is saved in DataLad `55b9249`. The17 October closed pilots, installer checks and VEP recovery are saved in `243dca6` (70 explicit files; active outputs excluded).

ClinVar coverage is explicit:64/369 pilot target transcripts are absent from VEP;
278of4,594,726 source AlleleIDs lack a summary match and40,067 joined rows lack an
exact MANE target. Original exclusions remain and can omit evidence in regenerated
gene statistics; candidate coverage/scientific activation stays false. Pinned
GoFCards UCSC reference/chain and ClinVar VCV acquisition completed in1515133/1515134
(DataLad `ac93655`/`b4d9a26`). Closed ClinVar pilots and installer follow-up are saved
in `2341411`. GRCh38 preparation1515142 closed **FAILED1:0** only at the frozen Conservation
URL404, after seven resource stages completed, including MaveDB. Separate
Conservation recovery1515153 remains active; mechanism source-only1515143 completed. GoFCards summary acquisition1515156 completed0:0 after the retained1515148
publisher502 attempt: all2034 unique alleles have summaries and GRCh38 positions;
1245 have ClinVar accessions. The closed receipt and both immutable snapshot hashes
were verified. Full mechanism1515158 is now running; active consumer inputs remain
unannexed until a safe boundary. Actual GoFCards pilot1515145 passed three
genes and239 transcript annotations across both assemblies, with no reference
rejects or liftover misses. Closed1515146 acquired G2P/Orphadata and preserved the
pinned HPO/MONDO bytes (DataLad `b858c73`); full mechanism qualification is pending.

The live full ClinVar116 run1515136 has retained LOFTEE exceptions (`Slice start
cannot be greater than slice end`) and transcript-assembly diagnostics. Official
VEP can catch a plugin exception and still exit successfully; the existing schema
check alone cannot qualify affected evidence. Read-only all-warning/affected-record audit **1515197 is queued** after the
producer; eight synthetic fixtures passed. It will inspect the retained VEP-only
artifact and closed log, with PID prefixes kept distinct from variant identifiers.
Later population processing copies that artifact before mutation, so the original
annotation remains available. Derived builds may proceed, but full scientific
qualification awaits explicit warning and coverage assessment; no plugin behavior
or tolerance is changed.

The old Conservation URL returns404. Both installers now explicitly select the
archived release114 **91-mammal** GRCh38 candidate, preserving the configured
species set. Release113 and114 files differ; neither is claimed identical to the
unavailable old-current object or LoFtee's separate GERP track. Acquisition runs in1515153 alongside release116 peptide/GTF preparation;
bounded original Conservation/MAX plugin check1515163 waits for it. Profile
qualification remains pending. No92-mammal substitution is made. Independent
GRCh38 AlphaMissense raw bytes are verified by1515172 after1515167 failed
an additive chromosome-prefix check; no redownload or score change was needed.
Fixed20261011 ClinVar raw1515168 completed: all4,594,764 rows survive the original
chromosome-alias transform. AlphaMissense96-row pilot1515174 waits for Conservation
acquisition1515153 and actual plugin check1515163. Full GRCh38 ClinVar annotation
remains gated on SpliceAI2 scientific qualification; it needs its own complete
join/coverage audit and cannot reuse hg19 coverage counts.

An isolated release116 BioMart attempt1515149 failed closed: the registry returned
mart116 metadata, but the dataset request returned Ensembl's service-unavailable
HTML. No table was acquired and no mapping changed. The explicit mart113 GRCh38
candidate remains separate from VEP116 compatibility/coverage qualification.
Ensembl's [23 September service report](https://www.ensembl.info/2026/09/23/ensembl-biomart-sep-2026-status-report/)
confirms that the old mirrors are retired. The official
[release116 mart dump](https://ftp.ensembl.org/pub/release-116/mysql/ensembl_mart_116/)
supplied the completed release116 candidate in **1515184 (0:0)**. The five-column
export has691,498 rows /670,670 transcripts and28,785,905 bytes; SHA256
`dfae01f608b6a9093c2b351ca3fe162af0aacc363fa4a09ba2e4dd5de4850243`.
Publisher checksums, schema, full-table coverage and join multiplicity checks
passed, including288,242 untranslated transcripts and66 retained duplicate rows.
Its deterministic dump order is new; neither original BioMart response order nor
identity to the unavailable original table is claimed. Consumer/HCSeeker/VEP116
coverage is audited below; compatibility qualification remains pending and the
existing explicit mart113 choices are retained; new VEP116/GRCh38 installation targets the separate FTP116 candidate, pending scientific qualification.

CADD profiling1514598 is closed with **80/80 exact table matches**; no rerun is
needed. Its **NaN CLI exit gate** and the separate experimental **PID-only
cancellation gate** remain unresolved. Experimental calculation stays disabled.

Final narrow installer fixes preserve legacy cache paths, canonicalize assembly
aliases only inside the motif/domain/PM1 wrappers, and initialize the LoFtee cache
directory before use. The corrected installed-environment harness passed179
installer tests per checkout; its earlier missing-PATH failure is retained. Later
Conservation integrity checks passed47 cases per checkout, including actual BigWig
reads, assembly/reference matching and invalid-download barriers. These overlap
existing suites; they are not scientific acceptance. The14 SpliceAI runtime guards
also pass per checkout, retaining the explicit GRCh38 no-fallback policy.
Independent review caught a receipt-label error: both upstream CADD assembly
configs and ESM peptide labels specify110; GRCh37 coordinate-reference tracks
have lineage75. These are now separate corrected provenance fields, with26
CADD runtime tests passing per checkout. The historical v4 verifier1515170 was cancelled before starting; Python-only
replacement1515208 retains the corrected metadata; underlying model/scoring assets and historical receipts are unchanged.

The earlier combined installer regression **1515201 completed 0:0**: **446 tests
for PriVA and 454 for PriVA-nextflow**, zero failures/errors/skips, with source
and harness hashes unchanged. This supersedes the earlier 283-case checkpoint
for the tested operational routes; receipts retain each source identity. These
are installer fixtures and contracts, separate from real-model/public/family
acceptance.

The final fixes preserve the caller directory during UTRAnnotator installation,
honor explicit CDS/ClinVar VCV paths (including symlink aliases), and bind custom
VCV reuse/injection/builds to the selected XML path and SHA256. CDS cache names
cannot overwrite plain FASTA inputs or Nextflow private copies, and quoted paths
reach consumers intact. GRCh38 aliases include MaveDB; hg19 still omits it as
approved. Confirmed download/index/config-write failures now propagate. Scientific
parsing, filters and calculation order are unchanged; no replacement input is guessed.

Both installers now prepare missing opposite-assembly VEP caches for the original
GoFCards builder, preserve explicit ClinVar inputs, propagate LoFtee failures, and
retain InterPro XML with observed SHA256 for reuse. Focused checks passed69 cases
per checkout; GoFCards checks passed49 each, with12 focused portable checks each.
Actual hg19 and GRCh38 LoFtee validation passed, including all55 shared GRCh38
reference chromosomes. These are installation checks, not family acceptance.

GRCh38 CADD prescores were separately absent: **1515185 is running**, using all
four matching Washington data/index MD5 pins. Its683,150,170,257-byte acquisition
writes to `/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade/resources/CADD/prescored/GRCh38_v1.7/incl_anno/`.
No partial is published; exact MD5 and bounded index queries precede compatibility
links. This is independent of annotation archive1515113 and Python-only verifier1515208.

Alpha37 coverage audit1515181 now precedes the original in-place derived stage
1515139. SpliceAI2/ClinVar38 pilot1515183 **completed0:0**:191 original alleles,
101 exactly scored alleles,2,138 scored transcript rows and zero score mismatches;
90 coverage-negative alleles remain. Ten alternative-contig REF bases are unavailable
in the selected FASTA, without adding check_ref or dropping rows. Eighteen official
plugin diagnostics concern an unrecognized version token in CSQ descriptions;
explicit v2 provenance and all numeric fields are intact. These diagnostics are
retained, not suppressed. Full annotation-only qualification **1515186 failed**
its all-source preservation audit after VEP completed; independent postcheck
**1515187 failed** an environment-state hash guard before annotation review.
Recovery is under investigation; scientific activation stays gated. Original consumer thresholds and publisher v2 guidance
remain distinct. Closed pilot evidence: DataLad8c3bb86.

Configured HPO paths now reach both shell consumers and the Nextflow resource
manifest. Fifteen focused path/default/freshness tests passed per checkout, plus
five Nextflow identity tests. The earlier overrestrictive inspection-fixture
allowlist failure is retained; its corrected case passed. Scientific Python
execution/mutation order is unchanged. InterPro/ClinGen explicit output paths now
pass33 focused output/failure checks per checkout; manual-curation InterPro reports
stay with source reports while binary/XML identity bundles follow the selected
data path. Managed-reference placement under existing `cache_root` passed **42
focused checks per checkout**. Fresh managed FASTA/index/receipts and liftOver
chains use that cache root; existing selected paths and their validation remain
unchanged. Missing explicit custom inputs fail without a replacement.

BioMart116 audit **1515188 failed1:0** at its final source-identity gate after
scanning all5,066 regions /1,237,995 physical model records. Concurrent installer
reference-routing changes altered four helper pins; scientific readers/data and
scanner pins stayed unchanged. The failed run and completed scan are retained.
Recovery 1515195 failed final JSON serialization on a nullable source-cache tag;
that failed run and incomplete receipt remain retained. Corrected read-only audit
**1515199 completed 0:0**, preserving unknown tags explicitly. Registration
**1515200 completed 0:0**: both installers recognize the exact release-116 FTP
provenance for GRCh38, reject hg19, and leave both configurations unchanged.
The FTP route passed 54 focused tests per checkout; actual offline reuse also
passed. Existing explicit mart113 selections and frozen jobs remain unchanged. New VEP116/GRCh38 installation targets the separate complete FTP116 candidate; scientific qualification is pending.

HCSeeker mapped spots improve from 1,664 to 1,669 of 1,670; `NM_153447` remains
unmapped. Across 1,836 mapped spot/ENST pairs, 7 lack a native cache model, 174
lack a model on the spot chromosome, 6 lack a peptide, and 27 exceed peptide
bounds. These categories overlap; including exact RefSeq-version uncertainty,
202 pairs across 137 spots are flagged, plus the wholly unmapped spot. The full
mapping has 3,218 ENST IDs absent under every native tag and 47 with an unknown
source tag. No affected row is discarded and no protein position is remapped.
Motif retains ordered multiplicity; PM1 retains the first occurrence of each
distinct NM-to-ENST pair, including all distinct ENST mappings. New-table hotspot
anchors gain 162 and lose 82 transcript/residue pairs; coldspots gain 339 and
lose 558. Compatibility and derived-resource qualification remain pending.

The future SpliceAI2 header fix passed **45 Rust tests** and actual three-allele
VEP116 comparison: complete body, gene/FeatureID, raw scores and order are exact.
Read-only postcheck **1515191 completed0:0**; official header diagnostics fall from
18 to zero while v2.0 provenance stays explicit. The original1515189 harness
failure (expected9 diagnostics instead of18) remains recorded. This qualifies the
source fix only; current deployed binary, large caches and active1515186 inputs
are unchanged. Full scientific qualification remains pending.

Deployment data and all logs for this checkout are under
`/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade`
(`resources/*`, `derived_caches/*`, `validation/*`, `deployment/*`, `logs/*`).
Established implementation stays in
`/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow`;
plans/reports remain under `/stornext/snfs170/next-gen/scratch/yangyxt/projects/PriVA_upgrade`.
The data-side README points back to both repositories and the plan folder.

Closed BioMart audit/registration and final installer evidence are recorded in
DataLad `143ef3d` (64 explicit data-side paths, including retained failed attempts).
All selected hashes were checked before intake; active outputs/logs were excluded.
The read-only final checkpoint independently confirmed the 446/454 test results
and unchanged source, plus successful BioMart receipts with qualification false.

The updated17 October storage boundary applies to new work: source/configuration
and plans remain under `yangyxt`; new generated payloads, caches, runtimes and all
logs belong under `yangyxt_dt`. Existing frozen writers retain their recorded
paths until safe relocation with compatibility links; no active cache is moved.


A comprehensive pipeline for ACMG/AMP-based variant interpretation and variant prioritization in genetic analysis.

## Nextflow / Rust migration branch

This branch develops a Nextflow workflow and one opt-in `priva` Rust/Polars
binary, with faithful function order and exact comparison to source `09223ef`.
See [execution and deployment](docs/nextflow_migration.md) and
[Rust scope and exact comparison](docs/rust-migration.md). The
[frozen ranking ledger](docs/ranking-source-ledger-09223ef.md) maps remaining
control/ranking mutations, SciPy dependencies and parser/failure contracts.

```text
One official VEP run → PriVA annotation view → family interpretation / Rust ports
                    ↘ CADD feature enrichment → complete annotation artifact
                       → direct original Python raw + PHRED
```

`main.nf` isolates existing PriVA tasks; `cadd.nf` develops independent annotation
and calculation stages. Shared-VEP transcript selection and complete feature
mapping must be qualified before production activation. CADD retains its original
Python calculation, scientific dependencies and models. Profile first; enable
Python optimizations only after exact score/row/order/missingness/failure checks
and measured timing/memory comparisons against the pinned unmodified reference.

The supplied **GRCh38 trio is delivered for input/implementation testing**;
input staging and packet metadata checks are complete. Twenty-five focused
adapter tests and eight packet metadata checks pass; multi-sample import remains pending.
Its separate single-sample VCFs are
not jointly genotyped. The workflow needs one multi-sample family VCF; import
remains to be implemented and validated. No merge or missing-call imputation
has been performed. A tested genotype-preserving merge may support implementation
parity when labeled nonjoint, with missing calls unknown. Unknown PED affection
does not supply phenotype-1 controls.
Input/control contracts, documented sex and qualified GRCh38 resources remain
acceptance gates; see the [current acceptance gate](docs/nextflow_migration.md#current-acceptance-gate).
No family annotation/ranking acceptance or public 1000 Genomes acceptance is claimed.

Rust currently covers final ordering and a bounded criteria-summary/posterior
interface, not all criterion assignment. Cluster software installation has passed;
cache deployment and calibration are still in progress. Synthetic checks, real
CADD model tests, public benchmark acceptance and family acceptance are distinct.
The public benchmark and supplied-family exact ACMG/ranking gates remain pending.

The actual installer reproduced and fixed VEP's unattended false success. The
separate reference/VEP deployment job 1514579 has now completed (0:0) at
2026-10-16T08:21:30Z; deployment evidence is recorded below. Full installation
and scientific qualification remain pending.
Both 18-row CADD assembly fixtures now match original Python exactly through
Nextflow. One `NaN` edge still changes the wrapper exit status when original
Python emits a malformed table; this mismatch remains reported and unaccepted.

Baseline profiling **1514598 is closed: FAILED (1:0)**, 2 h 21 m 24 s,
ended **2026-10-16T02:51:05Z** (scheduler 2026-10-15T21:51:05-05:00).
Frozen harness `e3b7818ed4362fda089a6928a184e3c4cf6ea894` completed all ten
reports: 20 assembly evaluations and **80/80 exact model-table comparisons**,
zero table mismatches. The only recorded mismatch is the NaN PHRED CLI exit:
four of five synthetic edges passed; NaN stdout bytes match, but exit 0 becomes
1 without a completed artifact. The bundled report keeps acceptance=false;
nine scaled reports pass only their bounded model inputs. The shell retained
that first failure, completed all nine scaled invocations, then exited 1.
See the [closed-run provenance](docs/nextflow_migration.md#cadd-baseline-profile-1514598).

The separate experimental CADD batching branch has passed 44 exact original
transformation comparisons and 3 helper cancellation cases. Completed Slurm
benchmark 1514607 passed all 24 paired synthetic trials, including 96 exact
intermediate/final-table comparisons. At 100,008 rows, median calculation time
fell from 556–578 seconds to about 279 seconds, with similar sampled summed
process-tree RSS (this estimate can count shared pages repeatedly). These results
do not qualify shared annotation, failure behavior, public data or families.
The historical tested original-Python integration used standalone
`4194be32320292cb8acd290ee0210c19456e2573`; that integration is now suspended.
The later wrapper audit failed PID-only SIGTERM cleanup for the experimental
shell, profiler and standalone CLI; process-group termination passed and the
audit reaped all remaining children. The candidate remains disabled.

SpliceAI v1 originals/indexes and both PrimateAI v0.2 originals are acquired.
The hg19 SpliceAI import 1514615 and PrimateAI import 1514620 completed with
exit 0 and valid receipts. The PrimateAI receipt includes exact external
source-to-formatted validation; cancelled job 1514618 never ran. Original hg38
files remain retained and unselected. Separate main-config promotion updated only
the three selected hg19 resource-path keys with a receipt; scientific qualification
remains pending. gnomAD job 1514601 failed its `.vcf.bgz` extension validation;
replacement 1514729 is running with the minimal fix (eight focused tests passed).
At that earlier checkpoint, AlphaMissense 1514589 was held with reason
`JobHeldUser` for the unavailable original NM-to-ENST mapping. The separately
authorized [BioMart table is now acquired](#biomart-mapping); it does not establish
original-profile identity or release that frozen job. No evidence is silently omitted.
The [installation dependency ledger](docs/installation-dependency-ledger-09223ef.md)
records the ordered calibration recipe and separates bundled ClinVar baseline
statistics from its currently unqualified fresh-build path.

Resource update, **2026-10-16T08:46:04Z** (SpliceAI2/held-job check around 08:45Z):

| Installation work | Verified state | Remaining gate |
| --- | --- | --- |
| GRCh37 reference + VEP113 merged cache 1514579 | COMPLETED (0:0) at 2026-10-16T08:21:30Z after 8 h 39 m 49 s. Both stage receipts are deployed. Full archive checks, FASTA/index samples and bounded canonical-cache layout/header checks passed at 08:46:04Z, using frozen `fc21869`. The existing config already selects the reference and cache. | No real VEP/plugin annotation was run by this audit. Shared-VEP/CADD, public-sample and supplied-family acceptance remain pending. AlphaMissense remains held for the original NM-to-ENST mapping. |
| CADD prescore recovery 1514744 | COMPLETED (0:0), 13 m 15 s, with a deployed receipt for the four Washington GRCh37 v1.7 files. Existing SNV/index and indel bytes passed full MD5 rechecks; only the missing matching index was downloaded. Real indexed-pair checks passed all 24 contigs in each SNV/indel pair. Source `20774bd` has 20 deployment and 31 licensed-import fixture tests passed. | Deployment and bounded indexed access are verified. Real VEP/model, shared-annotation, public-sample and supplied-family acceptance remain pending. |
| SpliceVault 1514603 | COMPLETED (0:0) at 2026-10-16T07:03:13Z; the original release113-to-GRCh37 derived resource is built, verified and already configured. | Real VEP use, scientific qualification and supplied-family acceptance remain pending. |

The [deployment checkpoint](docs/nextflow_migration.md#current-cluster-build-checkpoint)
records the reference/VEP audit, recovery guards, closed test evidence and SpliceVault verification.
Reference/VEP job 1514579 used frozen source
`fc2186901029acc88cd41bb004828778c1fd1cb8`, distinct from the CADD recovery source.
DataLad `80bb1f0` records its closed receipts/log and
`validation/vep-reference-receipt-verified-1514579.json`, with descriptive metadata
recorded and all five selected paths clean.
Recovery source `20774bd` passed [CI 38869590348](https://github.com/yangyxt/PriVA/actions/runs/38869590348).
DataLad `9e70725` records the closed recovery receipt/log and indexed-pair report/log.
The documentation base `b171466c31d9e695af22c943e824f024b238b482` passed
[CI 38871801476](https://github.com/yangyxt/PriVA/actions/runs/38871801476), checked
at 08:46:31Z; these documentation edits are subsequent.
Installation progress leaves the original-Python pin and experimental-backend
acceptance boundaries above unchanged.

## Overview

This repository contains scripts and configuration files for an automated variant prioritization pipeline based on the ACMG/AMP guidelines. PriVA combines transcript-level annotation, population and functional data, gene-disease mechanism and inheritance context, and family information to evaluate ACMG/AMP criteria and produce ranked variant lists for each family. Annotation resources are deployed once and can be refreshed; each cohort run then performs annotation, family-level filtering, criteria matching, and prioritization.

![PriVA workflow (Figure 1)](PriVA_figure1.png)

## Requirements

- Micromamba, Mamba, or Conda for environment management
- Sufficient disk space for annotation resources (>1TB)
- Reference genome file

## Pedigree File Requirements

If a valid pedigree table is provided, please note:
1. The family sample ID should be consistent between the pedigree table and the VCF file.
2. The proband of the family should stay as the first record of the family in the pedigree table
3. Every family should have a unique identifier in the pedigree table

## SpliceAI2 GRCh38 preparation

Current installation direction: keep only the original spliceai_snv_prescore and spliceai_indel_prescore YAML keys. Choose assembly once at installation: fresh GRCh37 uses v1 and fresh GRCh38 uses v2. Inspect decompressed headers automatically: an existing compatible VCF passes through; native 84-column SpliceAI2 TSV invokes the narrow adapter. A prepared pair is registered under the two original keys and installation succeeds. Resource/model identity, hashes, schema and qualification remain in sidecars and cache headers, without extra selector or converter YAML switches. Existing compatible GRCh38 v1 VCF pairs pass through; fresh GRCh38 installations select v2. Wrong assembly, mixed input types and failed v2 preparation never trigger a fallback.

The [preparation route](docs/nextflow_migration.md#spliceai2-preparation) keeps the
original 84-column tables and creates a separate top0 VCF for the official plugin,
using exact scored-region/GENCODE48 gene context. The Rust converter streams rows;
the Python deployment driver verifies the 47-shard manifest and controls
`gzip`/`bgzip`/`tabix`. It downloads nothing. GENCODE48 has passed its resource
checks. Download job 1514724 completed (47 final shards, zero partials), and the
complete 20,083-region stable-ID join passed with exact bounds/strand/biotype.
Full GTF IDs/versions and `_PAR_Y` remain in context; ten same-strand symbol
groups coalesce labels only per prediction, while opposite-strand symbol contexts
fail. The execution manifest is complete and saved in DataLad commit 5f5946b.
All 47 full source SHA256 checks passed in frozen job 1514732. At the bounded
2026-10-16 ~08:45Z scheduler/log check it was RUNNING on c81o-03 and its log
had reached SNV chr4; complete conversion and indexed-VEP checks remain pending.

```text
Assembly → decompressed headers → VCF passthrough / native84 adapter
                                → indexed VCF pair / two original YAML keys
                                → sidecars + cache-header scientific gate
```

The simplified installer succeeds after preparing/registering the indexed pair
under the two original resource keys; provenance stays in sidecars.

SpliceAI_install invokes scripts/spliceai_cache_input.py install automatically, retaining the existing prompts and two resource keys. The helper detects gzip by magic bytes, then classifies a bounded decompressed header as vcf_v1, vcf_spliceai2 or the exact native 84-column SpliceAI2 TSV schema; filenames, extensions and score magnitudes do not select the predictor. An unmarked VCF receives the legacy-compatible vcf_v1 routing label; that label alone does not prove predictor identity. Existing compatible VCF pairs pass through, including existing GRCh38 v1; fresh GRCh38 selects v2. Wrong assembly, mixed types or failing v2 preparation fail without fallback. SpliceAI2 VCF passthrough additionally requires a completed adjacent adapter.ready.json with exact SNV/indel/index membership and complete published-catalogue status. This checks completed-publication metadata and recorded membership, not a new full VCF/index-body checksum; it performs no conversion. Reject .partial or .building path components, including symlink targets, and refuse in-progress preparation or staged outputs. A final VCF renamed before receipt publication is not ready. V1 VCF passthrough requires no additional adapter receipt.

```bash
python scripts/spliceai_cache_input.py install --config CONFIG --cache-dir CACHE_DIR --snv SNV --indel INDEL
```

For native TSV, the helper discovers the fixed adjacent SpliceAI2_GRCh38.downloads-complete.json sidecar. The selected SNV/indel inputs must match their manifest roles and the full 47-shard catalogue. The existing low-level SpliceAI2_install preparation/reuse path performs full checksum verification and prepares or reuses a validated receipt, then updates only spliceai_snv_prescore and spliceai_indel_prescore. The driver builds or reuses the standalone installation converter from scripts/spliceai2_converter.rs; no production priva executable or extra YAML option is required. Provenance and scientific qualification remain in sidecars/cache headers. The current adjacent sidecar at yangyxt/caches/SpliceAI2_GRCh38.downloads-complete.json has SHA256 bf546817c03689359138d9230e3587fe225764a06546b3914c8819afbbd2bdde: all 47 original source identities remain unchanged; the qualified plugin pin is116. The historical manifest remains preserved. Historical full-cache conversion used the PATH entry yangyxt/softwares/bin/priva selecting priva-tools-spliceai2-mnv-20261016/bin/priva, SHA25670eefad75c6ad5a8001fcd87ca6970c3b0b1495f21b55f8c7775db67122c567a. This changes no YAML keys beyond the original two paths.

Two-key interface qualification: code 35779c4578509ce01218a052a8e7e1d8a785a75b passed 239 Python tests, zero failures/errors/skips, in 190.848 seconds: 39 helper, 25 installer, 14 runtime-guard, 21 Nextflow-adapter and 140 existing deployment/driver/BGZ tests. The installer-shell cases mock the helper boundary; they are not real conversion evidence. Log: yangyxt_dt/projects/PriVA_upgrade/logs/spliceai-two-key-integration-20261016.log. This supersedes the earlier 197-test interface checkpoint without adding its count; unchanged Rust evidence remains 41 tests and the earlier real-projection evidence remains 512 unique rows.

Remote CI for integrated source `2cc8ecb5326f19a35bb86fc4fd98cc1cf4d61963` passed [Actions 38854687155](https://github.com/yangyxt/PriVA/actions/runs/38854687155). The earlier eb16493 run also passed; the f5fc377 run was cancelled after being superseded. Local two-key interface evidence remains tied to its tested 35779c4 source.

Actual resource checks used final helper SHA256 f69d6ea8d6179d68ae44cad849b3c6c7912a879f05e305a36acffba09c81bb88: all 47 native TSV and four v1 VCF headers classified correctly, including assembly metadata. Real SpliceAI_install reused the existing v1 pair for hg19/GRCh37 and hg38/GRCh38 with configuration bytes, inode and mtime unchanged and resource stat identities unchanged. Receipt: yangyxt_dt/projects/PriVA_upgrade/validation/spliceai-two-key-real-install-20261016/receipt.json. These are bounded header/reuse checks, not new body checksums, full VEP execution or model/scientific qualification.

Nextflow stub verification passed all four scenarios: raw input, resume with all four tasks cached, preannotated input and stale-environment rejection. Log: yangyxt_dt/projects/PriVA_upgrade/logs/spliceai-two-key-nextflow-20261016.log; retained validation/log artifacts: validation/spliceai-two-key-nextflow-20261016/ under that data project. The run used code 35779c4 plus the documentation state at execution; its recorded source_digest identifies that exact state. Subsequent documentation edits are not attributed to the run.

Historical frozen-run boundary: job 1514732 used installer f8de4b3 and converter225aa4a, then failed at the publisher MNV on indel chr1 line11929. It produced no completed adapter receipt. Its explicit-manifest/six-key interface and197 Python tests remain historical provenance. Recovery retains all47 original download identities and the two-key public installer interface; only validated complete outputs can later be imported. Full scientific acceptance remains pending.
GRCh37 v1 and legacy CADD masked-v1 inputs remain unchanged. The real-plugin
boundary harness passed 11 subtests / 38 assertions with synthetic tabix/VEP
objects. Historical checks before the stable-ID/same-strand-symbol update: Rust
passed 30 new plus 9 existing tests; the integrated installer/driver run passed 165 tests with zero failures/skips, including real Rust/bgzip/tabix
on synthetic inputs. These historical totals remain tied to their producing source.
Current Rust 225aa4a passed 32 projection + 9 existing tests; Python f8de4b3 passed
197 tests with zero failures/skips (169 installer/deployment + 20 adapter + 8 BGZ).
A bounded real-cache check passed 256 SNV + 256 indel chr1 positive-strand rows
with exact full-record/order and bgzip/tabix checks. It does not establish model
inference, whole-cache, full VEP, public-sample or family acceptance.

All 47 publisher SHA256 values and VERIFIED acquisition receipts are now supplied
in [caches/README.md](../../caches/README.md). The complete execution manifest and handoff were saved in DataLad commit 5f5946b before launch. Full preparation job
1514732 is running on c81o-03 (4 CPUs, 16 GB, 7 days), using frozen installer f8de4b3 and
the isolated GRCh38 preparation config; all 47 full source SHA256 checks passed.
The bounded 2026-10-16 ~08:45Z check confirmed RUNNING and a log through SNV
chr4; complete conversion is pending. Optimized
source 225aa4aa774b4be7f53576e1cde56af71e67bd2c is installed at
`yangyxt/softwares/priva-tools-225aa4a/bin/priva`, SHA256
`bbc3b7288ed2ea63f1f8455ed87b8a9d08c7bc1fe0a562c1d6134525e2a12dc6`.
The installed release passed the same 25 driver cases and same 512 real-prefix
rows again with zero mismatches; the frozen checkpoint had 197 Python, 41 Rust and 512 rows; the new interface
has 239 Python tests, with unchanged 41 Rust tests and 512 unique rows. A stale shared Cargo artifact was caught before deployment, rebuilt and
checked at the installed path; the receipt retains that failed attempt. Full
conversion and real VEP/scientific acceptance remain pending. See the linked route
for the exact manifest, GTF/plugin identities, commands and pending gates.

Scientific boundary for the simplified installer: inspect the actual selected cache headers and sidecar provenance to prevent unqualified v2 scientific use, including preannotated import and resume; no new enable/disable YAML toggle is introduced. Successful installation is separate from v2 ACMG/ClinVar, full VEP, public-sample and supplied-family acceptance. Unmarked external annotations still require pinned scientific provenance because predictor version cannot be inferred from numeric values. The runtime guards now inspect actual cache headers at the two original keys, including preannotated launch/task paths; the earlier 35779c4 checkpoint included 21 Nextflow-adapter and 14 runtime-guard tests once in its 239-test total. The later 25 focused PED-adapter tests are separate evidence; the full suite has not been rerun. Earlier f8de4b3 checks remain historical; scientific acceptance is still pending.

## Installation

VEP release selection follows the existing `vep_installed_version` key in both
installers. At the user's request, the live HGSC environment was upgraded in place
to **VEP116.3** in Slurm1514975, completed2026-10-16 19:36:07 UTC. Full software
checks passed; all 527 non-VEP package identities were unchanged. All 12 required
official plugins loaded with API 116 in job 1514978; LoFtee was unchanged. Both
assembly merged116 caches are deployed; recovery1515116 completed on17 October.
The approved four obsolete113 core/cache/package payloads were removed after
verification. CADD110 reference assets remain required. Historical
AlphaMissense1514589 stays held because its frozen wrapper requires113; the new
isolated116 candidate runs as1515130. Derived and scientific qualification remain
pending. Successful module imports alone did not prove the runtime plugin search
path: the ClinVar pilot exposed a comma-separated directory bug, now corrected
with a representative rerun1515131.

```text
Requested VEP major → matching executable/API → official release plugins
                                             → same-release assembly cache
```

For example, `vep_installed_version: "116"` requests release 116 when running
`bash scripts/install_utils.sh main_install CONFIG`. This selects an available
release rather than a hard-coded latest version. The existing `priva_acmg`
environment is changed in place when a different major is requested; no second
VEP environment is created. Every other installed Conda package is constrained to
its current name/version/build. Unavailable or conflicting releases fail without
falling back or replacing the YAML request with the old installed version.

`scripts/vep_plugins_release.py` verifies staged module names and records hashes
before successful reuse. A durable transition record supports retries after a
first upgrade; unknown later edits and same-major customized plugins are preserved.
Plugin files publish individually, so an interrupted transition is incomplete until
a verified retry succeeds. Missing/incomplete caches, wrong assembly metadata or
missing converted indexes stop installation. No cache is deleted automatically.
The installer does not establish equal annotations across releases: a changed VEP
release requires qualification of derived caches, transcript mappings, CADD shared
features and family results before adopting the new scientific profile. Manual
annotation-cache inputs remain required as documented below.

The current Nextflow cluster recipe now pins116.3 for the authorized replacement;
historical113source snapshots, before-change configs and scientific references
remain preserved. The suspended standalone4194be3 evidence remains retained. Active CADD calculation
uses original Python directly; reference-only VEP110 assets retain their pinned model role. The experimental CADD path remains disabled.

The route follows the official [VEP installer protocol](https://github.com/Ensembl/ensembl-vep/blob/release/116.3/INSTALL.pl)
and [same-release cache recommendation](https://grch37.ensembl.org/info/docs/tools/vep/script/vep_cache.html).
The common installation protocol permits major-release selection; availability
and scientific compatibility still require validation for each selected release.

Earlier validation on HGSC, 2026-10-16: **165 focused fixture tests passed per checkout**
in Slurm1514363.346 (330 executions; the earlier 51-test checkpoint overlaps).
This covers requested versions, manager failures/conflicts, API/helper mismatches,
plugin receipts/custom edits/interrupted upgrades, cache/index failures and affected
installer routes. Pytest was unavailable: six unittest files and all 17 GoFCards
functions ran through a retained in-memory unittest adapter supplying only `tmp_path`;
no test dependency was installed. The earlier missing-pytest attempt is retained.
Read-only real VEP113 executable/API/installer checks and the existing GRCh37 cache
checks passed separately, including the final CSI-only cache guard. At that earlier checkpoint, no real VEP116 installation/download, model run,
public-sample or family acceptance had been performed.
Closed commands, harness and producer hashes are retained under
`yangyxt_dt/projects/PriVA_upgrade/logs/vep-release-*.log`.

Actual116installation exposed a manager behavior absent from the fixtures:
`--freeze-installed` froze VEP itself. Both installers now omit that redundant flag
while retaining exact version/build constraints for every other installed Conda
package. Real install/update solver plans changed only VEP113.2→116.3; **56 focused
bootstrap tests per checkout** passed in Slurm1514974. This overlaps the earlier
suite and is not an additional distinct-test total. Original diagnostic113software
and its failed cache partial were removed; production113corecaches will be removed
only after verified116replacement. Keep licensed caches and scientific model assets.



<a id="biomart-mapping"></a>

BioMart mapping acquired on HGSC, 2026-10-16: the new table is bundled as
`data/hgnc/ensembl_GRCh37_mart116_nm_to_enst.tsv`, with its adjacent
`.provenance.json` receipt. It contains **109,530 rows / 5,084,679 bytes**;
SHA256 **`f730640e2411807e7bc38dad141acafc16c3ffca5e4d8ef500976ebd586ddbc7`**.
Its pinned source is [Ensembl GRCh37 BioMart](https://grch37.ensembl.org/biomart/martservice),
`ENSEMBL_MART_ENSEMBL` / `ensembl_mart_116`, dataset `hsapiens_gene_ensembl`,
assembly **GRCh37.p13**. The mart release is **116**, while the underlying
[GRCh37 gene annotation is based on release 75](https://grch37.ensembl.org/index.html).
These identifiers are distinct from the installed VEP release; VEP was not changed.
The failed/timed-out release-113 archive export was not relabeled as this source.

Acquisition job **1514785 completed (0:0)** at **2026-10-16 14:25:24 UTC**.
The retained raw response includes BioMart's terminal `[success]`; registry/dataset
metadata matched before and after export. The receipt records the exact query,
metadata, raw/table hashes and frozen producer SHA256
`944e3881be3f550c2cbcba554e3851f74839e78b83b0a1d569ddcfcdabd793f3`.
The five attributes remain in this order: `ensembl_gene_id`,
`ensembl_transcript_id`, `hgnc_symbol`, `uniprotswissprot`, `refseq_mrna`.
No table sorting, ID rewriting, row deduplication or first-transcript selection
was performed. There are **35,555 NM accessions / 35,812 distinct NM–ENST pairs**;
177 NM accessions have multiple ENST mappings (maximum 10). All 179 exact duplicate
rows and 72 repeated NM–ENST rows remain. A legitimate `NC_000016` row is retained;
both unchanged readers ignore non-NM accessions. The initial job 1514784 failed
because the first validator incorrectly rejected that row; its closed raw response
and failure log remain retained.

For `hg19`/`GRCh37`, an unset `nm_to_enst_tsv` selects or retrieves this pinned
bundle; explicit supplied paths remain supported. Both motif and PM1 calls receive
the resolved path. Existing derived outputs require matching mapping/anchor-input receipts and
output size/SHA256 before reuse; a changed or unrecorded mapping requires an explicit new output path
or recovery of the old artifact. GRCh38 selects its separately pinned release113 bundle described below;
this GRCh37 bundle is never an implicit fallback.

GRCh38 acquisition recovered on 2026-10-16: the original release113 BioMart
endpoint responded from HGSC at 18:35 UTC. After a completed three-gene probe,
Slurm **1514907 completed (0:0)** at **18:37:02 UTC** in 28 seconds. The exact
[October 2024 BioMart](https://oct2024.archive.ensembl.org/biomart/martservice)
source is `ensembl_mart_113`, dataset `hsapiens_gene_ensembl`, **GRCh38.p14**;
its gene annotation was updated July 2024. The table is bundled as
`data/hgnc/ensembl_GRCh38_mart113_nm_to_enst.tsv*`: **150,650 rows / 7,324,483 bytes**,
SHA256 `bb5f64c3c5eef06da92c95e50068d716c7e288629733efea6eb34dbab7bc1540`.
The five-column query, terminal success, before/after metadata and all row order
and multiplicity were verified. It retains 123,845 ENSTs, 66,718 NM IDs, 74,240 distinct
NM–ENST pairs, 56 duplicate full rows and 85 repeated NM–ENST pairs; 3,918 NM IDs
have multiple ENSTs (maximum 42). The initial 1514890 timeout is retained as failure
evidence, not evidence of a wrong URL or retirement. No 115/116 source was selected.

Both unchanged scientific readers passed independent full-table and actual-anchor
checks in the existing PriVA environment. GRCh38 HCSeeker coverage is
**1,664/1,670 rows / 846/848 NM IDs**. The unmapped `NM_001320643` contributes two
coldspot rows and `NM_001321073` four hotspot rows; these contribute no HCSeeker
anchors under the unchanged readers. This limited coverage may affect motif/PM1
calibration; RMC ENST anchors are separate. Identity to the unavailable original
mapping, derived calibration, public-data and family acceptance remain unclaimed.
No held AlphaMissense job was released. Closed evidence is retained under
`yangyxt_dt/projects/PriVA_upgrade/validation/biomart-grch38-113-*`.

Both repositories passed **28 mapping-helper tests and 91 installer prerequisite,
reference, bootstrap and boundary checks** each. Both assemblies' real bundles
passed offline verification/reuse. Missing explicit inputs, wrong-assembly bundles
and stale or absent provenance are rejected. These installer/fixture checks do not
establish calibration or family equivalence. Closed integration logs use
`yangyxt_dt/projects/PriVA_upgrade/logs/mapping_grch38_*.log`.

The mapping entry point is `bash scripts/install_utils.sh nm_to_enst_mapping_install CONFIG [RAW_OUTPUT]`;
`main_install` also calls it after prerequisite checks. Existing genuine bundled
tables require only existence and nonzero size for installer reuse. Fresh retrieval
and external-cache reuse retain their pinned source/provenance checks, including
assembly and row order. A compatibility symlink to an external table is not a bundle.
Scientific mapping compatibility remains a separate qualification step.

Coverage against bundled hg19 HCSeeker is **1,469/1,549 rows** and **754/805 NM IDs**:
882/911 hotspot rows and 587/638 coldspot rows map. The remaining **80 rows / 51 NM IDs**
are listed in the independent coverage receipt. The unchanged readers omit those
29 hotspot and 51 coldspot rows from the HCSeeker contribution, which can affect
motif calibration and PM1 evidence; RMC's direct ENST mappings are separate.
This newly acquired table is **not claimed identical to the unavailable original**,
and cross-references alone do not prove protein-sequence identity. Original-profile
parity, scientific qualification, public-sample and family acceptance remain pending.
No held AlphaMissense job was released by this mapping task.
Both repositories passed 231 installer regression tests, followed by 41 affected
mapping tests per checkout after exact snapshot pinning. Real-bundle offline reuse,
assembly gating and unchanged scientific-source checks also passed in Slurm 1514363.
See `tests/README.md` for fixture-versus-real evidence and retained log paths.

Installer repairs, 2026-10-16: `main_install CONFIG` bootstraps the environment,
then resolves the explicitly pinned or supplied mapping and checks remaining
inputs. Configure `nm_to_enst_tsv` as described above; resolve reported
old-host paths explicitly. A fresh GoFCards build needs both reference assemblies and the forward
hg19-to-hg38 chain. The new `assembly_references_install CONFIG` shell function
prepares both FASTAs and samtools `.fai` indexes; `gofcards_reference_install CONFIG`
also prepares the forward chain. Main installation prepares the selected reference
and, for fresh GoFCards builds, both assemblies/chain. Both merged VEP caches must
match the installed release; missing opposite-assembly merged caches are now prepared automatically.

Missing standard reference destinations use the exact UCSC original bigZips
[hg19 FASTA](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz)
(MD5 `806c02398f5ac5da8ffd6da2d1d5d1a9`) and
[hg38 FASTA](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
(MD5 `1c9dcaddfa41027f17cd8f7a82c7293b`). The
[forward chain](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz)
is pinned to MD5 `35887f73fe5e2231656504d1f6430900`. Publisher manifests and
bounded live FASTA ranges were verified; the full small chain matched its MD5.
The helper resumes staged downloads, checks source size/hash, validates assembly,
builds indexes and records source receipts before updating only reference keys.
Explicit original references retain their sequence/masking/contigs; a missing
selected DRAGEN/other reference fails instead of silently substituting UCSC.
The final reference/bootstrap/prerequisite suite passed **82 tests per checkout**
in Slurm1514363, including a completed-GoFCards reuse case with no opposite-assembly
acquisition. No full production FASTA download was repeated for these tests.
Run substantial acquisition/indexing under Slurm. Direct reference functions use
the prepared PriVA environment (Python/PyYAML and samtools); `main_install`
bootstraps that environment before invoking resource functions. HPO assertions
must use this checkout's `data/hpo/genes_to_phenotype.assertions.tsv.gz`, which the
unchanged runtime reads. Valid bundled ClinVar gene statistics are reused; missing
or corrupt statistics block the incomplete fresh-build route instead of running
it with guessed original inputs.
The user approved the existing hg19/GRCh37 MaveDB omission on 2026-10-16.
The hg19 config now explicitly leaves `mavedb_file` empty. MaveDB-only functional
PS3/BS3 support is unavailable; independent ClinVar evidence and the GRCh38
MaveDB route are retained.

Licensed caches remain user-supplied. Existing SpliceAI v1 VCFs are unchanged and
require matching `.tbi` indexes; bounded header/index/first-record checks precede
registration, without claiming full-cache integrity. GRCh38 native SpliceAI2 TSVs
use the adapter with the complete source manifest and gene context. Its standalone
installation converter is built automatically from the bundled Rust source; a
compatible host compiler is reused, or an isolated pinned compiler is acquired
through the existing environment manager. GRCh37 stays on v1. V2 preparation does not
qualify ACMG/ClinVar or CADD use. The original annotation/calibration entry points
reject unqualified v2 inputs.

CADD cache downloads now retain/reuse verified bytes, stage incomplete files and
validate data and indexes separately. GRCh37 uses the approved historical
Washington pair. GRCh38 data/index pins match the same-source
[Washington v1.7 MD5SUMs](https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/MD5SUMs),
verified 2026-10-16; both index pins are now set in `config.yaml`. This verifies
manifest pairing, not downloaded cache integrity or scientific qualification;
the existing acceptance gates remain in place. Live-manifest configuration checks
and five focused CADD fixture tests per checkout passed in Slurm 1514363.
Existing mismatched files are
retained and reported. Annotation archives
are retained for verified extraction/reuse. Scientific Python, official VEP and
plugins, model selection and score thresholds are unchanged. See
[scripts/README.md](scripts/README.md) and [tests/README.md](tests/README.md) for
implementation and test scope. Full scientific installation/acceptance remains
pending; successful installer fixtures do not establish it.


CADD HGSC runtime verification, 2026-10-16: the expanded Nextflow cluster runtime
stage passed **35 fixture tests** and real Slurm **1514901 (0:0)**. It verifies the
unchanged official image SHA256, five wrapper routes, Python2/ESM/RegSeq/MMSplice
imports, packaged MMSplice assets, and both assemblies' model/PHRED files. That historical check pinned source
`1fead0c1e02d3f0a8599bdb04e6b3a2dd70b4866`, the image and standalone `4194be3`.
The standalone dependency is now suspended; the source/model assets remain pinned. The first runtime attempt1514899 failed on unsupported VEP `--version`;
its evidence is retained, and the supported `--help` check now confirms110.1.
The packaged VEP is for reference validation; production still shares PriVA's one
official VEP run. No model inference or profiling was repeated.

The new receipt is
`yangyxt_dt/projects/PriVA_upgrade/deployment/cadd-runtime.verified-20261016.ready.json`.
GRCh37's checked runtime-asset group is 11/11 readable; the corresponding14 GRCh38
annotation assets remain missing. The five-link zero-copy hub at
`yangyxt/softwares/cadd-1.7.2/esm-hub` preserves original weights; a future qualified
profile must select it via `CADD_ESM_MODEL_DIR`/`--esm-model-dir`. No profile was
rewired. Complete unified installation, shared-VEP compatibility, the separate
NaN CLI and experimental PID-only cancellation gates, public-data and family
acceptance remain pending; the experimental backend is disabled. Closed profiling
1514598 remains 80/80 exact table matches and was not rerun.

The repository includes `scripts/install_utils.sh` to deploy supported resources from a YAML configuration. Original manual inputs and the prerequisites listed above must be supplied before a complete run.

**Carefully read comments in `config.yaml` to prepare for execution of install_utils.sh and PriVA**

`main_install` detects an available Micromamba, Mamba, or Conda executable,
then checks for the fixed environment `priva_acmg`. If absent, it creates that
environment from the repository's `acmg_conda_rough.yml`; if present, it reuses
it. The installer activates and checks the environment before reading the
project configuration with its `yq` and continuing resource installation.
The bundled environment includes `yq`, so it is not needed beforehand.

```bash
bash scripts/install_utils.sh main_install path/to/config.yaml
```

The environment file and name are built into the installer; the only argument
is the project configuration, which specifies resource paths and settings.
The legacy `conda_env_yaml` configuration value is ignored. Existing environments
are checked for required tools, VEP packages and usable `diptest`/`pyarrow` isolated Python
imports before resource annotation. Both modules are explicit dependencies of the
bundled environment. Missing modules or binary import failures stop installation
early; existing scientific packages are not silently upgraded. Miniconda installations are detected through their
`conda` executable.

Executable discovery checks `MAMBA_EXE`, then `micromamba`/`mamba` on `PATH`,
then `CONDA_EXE`/`conda`. It retains the selected executable or wrapper path.
Conda-backed Mamba installations use their sibling Conda for activation and saved
environment variables. Keep the manager's root, package cache, and environment
directories configured to the intended storage location. The installer does not
initialize shell profiles or choose a new manager root.

To prepare only the environment, the existing entry point remains available:

```bash
bash scripts/install_utils.sh conda_install_vep
```

Activation inside a child installer does not activate the caller's shell. Enter
the environment yourself before invoking individual resource functions directly.
The first resource installation remains assembly-specific; configure and prepare
its reference and other settings before running the full installer. With
`cache_root`, the six selected resource directories are created automatically.

Startup regression tests use isolated fake managers and do not install packages
or download resources:

```bash
python -m unittest discover -s tests -p test_install_utils_bootstrap.py -v
```

Installer startup updated 2026-09-13. Implementation and helper metadata are in
`scripts/README.md`; regression-test metadata are in `tests/README.md`.

`main_install` also deploys the gene-mechanism resources used by ACMG step 3 in
dependency order: pinned HPO/MONDO data, the broader LOF/non-LOF mechanism
cache, exact normalized GoFCards variants, the canonical non-LOF/ClinVar VCV
cache, and the integrated HPO-framed condition cache. All builders, parsers,
schemas, and outputs are inside the PriVA tree; `llm_gene_reranker` is not an
installation dependency.

The complete canonical JSON and ClinVar VCV/RCV/SCV field contract is documented
in `GENE_NONLOF_MECHANISM_JSON_STRUCTURE.md` and enforced by the bundled JSON
Schema.

The supported full manual refresh is the same entry point as installation:

```bash
bash scripts/install_utils.sh main_install path/to/config.yaml
```

Users who want automatic refreshes should schedule that command with their own
cron or cluster scheduler. PriVA does not install a separate weekly job. Each
run checks source-specific refresh intervals and rebuilds downstream caches only
when an input checksum or implementation file changes. A forced mechanism
refresh can be requested with `PRIVA_FORCE_ALL_CACHES=1`.

The mechanism stages can also be run independently for debugging:

```bash
bash scripts/install_utils.sh gene_pathogenic_mechanism_cache_install path/to/config.yaml
bash scripts/install_utils.sh gofcards_exact_gof_cache_install path/to/config.yaml
bash scripts/install_utils.sh gene_nonlof_mechanism_cache_install path/to/config.yaml sources-only
bash scripts/install_utils.sh gofcards_clinvar_injection_install path/to/config.yaml
bash scripts/install_utils.sh gene_nonlof_mechanism_cache_install path/to/config.yaml
bash scripts/install_utils.sh mondo_hpo_scope_install path/to/config.yaml
bash scripts/install_utils.sh hpo_condition_mechanism_cache_install path/to/config.yaml
bash scripts/install_utils.sh mechanism_resource_install path/to/config.yaml
```

The GoFCards steps must run in the order shown. `gofcards_exact_gof_cache_install`
writes only a build intermediate; `gofcards_clinvar_injection_install` turns that
into the one deployed cache every consumer reads. The `sources-only` call in
between downloads the weekly ClinVar VCV XML that injection needs, because the
component that owns that download is the non-LOF builder, and the non-LOF
builder's ordinary run cannot start until injection has finished. See
`docs/GOFCARDS_CACHE_SCHEMA.md` for the full ordering.

The integrated JSON is published by writing and validating a complete temporary
snapshot before atomically replacing the live cache. It nests external
mechanism histories and variants under an HPO gene/condition only through exact
OMIM, ORPHA, or MONDO identifiers; unresolved records remain under the gene's
`unmapped_evidence` section.

The MONDO/HPO installer downloads pinned releases, builds one disease-scope
registry, and adds the resulting scope and provenance fields to each HPO
assertion. Confirmed non-monogenic diseases are excluded from PriVA's flattened
inheritance/penetrance context. Neoplastic diseases without positive hereditary
evidence remain review-required rather than being assumed somatic-only. Manual
decisions are stored in `data/mondo/disease_scope_overrides.tsv` and take
precedence when the registry is rebuilt.

## Smoke test

The smoke test is a small end-to-end demo for checking that PriVA can run on the current machine. It uses curated demo assets in `example/smoke_data/` and does not run the full production annotation-resource installation.

From a fresh clone, first create and activate the PriVA conda environment:

```bash
git clone <PriVA_repo_url>
cd PriVA
mamba env create -f acmg_conda.yml
mamba activate priva_acmg
```

If `mamba` is unavailable, use `conda env create -f acmg_conda.yml` instead.

Then run the smoke test for one or both supported assemblies:

```bash
./example/run_smoke_test.sh hg19
./example/run_smoke_test.sh hg38
```

The wrapper checks required command-line tools and bundled smoke resources, downloads the configured UCSC reference FASTA if needed, creates the FASTA index with `samtools faidx`, and runs the PriVA Snakemake workflow on the 20-variant demo set. Expected outputs are written under `example/smoke_runs/`.

The smoke test does not require Git LFS on current `main`; all current smoke assets are stored as normal Git files.

## Configuration File
Before the first installation, set four values in `config.yaml`:

```yaml
cache_root: /absolute/data/storage/PriVA/caches
tmp_dir: /absolute/data/storage/PriVA/tmp
assembly: hg38                  # hg19/GRCh37 or hg38/GRCh38
vep_installed_version: "116"     # requested VEP major
```

The portable template contains no old-host paths. Create these two storage
directories on the host before running the installer. Their paths must be absolute
and dedicated, without spaces or shell characters. `main_install` stops before
environment setup, downloads or config changes if either directory is missing,
is a regular file, or is a broken directory link. It does not create these two
root directories; managed subdirectories are created later as needed.
The installer supports selecting
recent VEP majors; availability of matching plugins, caches and scientifically
qualified resources is a separate requirement. BioMart defaults are explicitly
pinned for supported profiles; an arbitrary VEP major does not authorize a mapping
release substitution.

```text
Choose cache_root + tmp_dir + assembly + VEP major
  → check settings → prepare environment → resolve paths and check inputs
  → acquire/reuse resources in the existing scientific build order
  → set family/cohort paths separately before a pipeline run
```

Run the existing entry point after editing those settings:

```bash
bash scripts/install_utils.sh main_install config.yaml
```

`ref_genome` is optional at installation: leave it empty for the managed UCSC
assembly/reference-index route, or supply an existing assembly-correct FASTA.
An explicitly selected missing reference fails; the installer does not guess a
replacement. The original GoFCards build prepares both assemblies, matching VEP
caches and the forward liftOver chain when needed.

Six managed directory groups use the optional `cache_folders` mapping:

```yaml
cache_folders:
  vep_cache: VEP_cache
  vep_plugin: VEP_plugin
  vep_plugin_cache: VEP_plugin_cache
  gnomad: gnomAD
  clinvar: ClinVar
  cadd: CADD
```

Each value is one folder name beneath `cache_root`, using letters, digits,
underscores, dots or hyphens. Existing managed paths follow changes to these six
roots. Explicitly supplied external cache inputs are adopted without changing their
bytes, with originals retained and the configuration updated only after successful
adoption. Genuine checkout snapshots stay in the repository; source scripts stay
source files. Completed SpliceAI2 projections retain original receipt provenance
when copied to their managed location. Additional unset managed paths derive from `cache_root` for
references, HPO/MONDO, runtime products and resources without a reusable canonical
bundle. Existing canonical InterPro, ClinGen, GoFCards and mechanism tables are
reused from the checkout; new generated outputs use the data storage root. Bundled
curated tables/scripts resolve against this checkout. CADD defaults to `runtime/cadd-1.7.2` beneath
`cache_root`, using original upstream Python. The suspended `cadd-rs` repository
and executable are not installation requirements.

For **data genuinely bundled in the repository**, installer reuse checks only
that the file exists and has nonzero size. It adds no content, schema or checksum
validation. An external downloaded cache does not become bundled data because a
compatibility symlink is under the checkout. New external acquisitions retain their
resource-specific validation; scientific readers and profile qualification remain
separate. The newly authorized AlphaMissense route validates only external publisher
prescores; bundled gene means remain existence/nonzero-only.

To preview portable defaults without acquiring resources or changing files, use
the installed environment:

```bash
python scripts/installer_config_paths.py resolve --config config.yaml --source-root "$PWD" --check
```

The existing six-group resolver also supports a preview; it additionally shows
old managed-root remapping:

```bash
bash scripts/install_utils.sh priva_configure_cache_paths config.yaml --check
```

Manual acquisition remains required for the licensed SpliceAI and PrimateAI
resources. Supply their completed matching-assembly files through
`spliceai_snv_prescore`, `spliceai_indel_prescore`, and `primateai_prescore`.
After environment setup, `main_install` checks all three selections before
resource acquisition or whole-table annotation. In a terminal, it prompts only
for unset keys and saves them together after all paths pass readiness checks.
Batch runs stop with the missing keys; invalid explicit selections are preserved
and reported. Prepared selections do not cause a config write. This early check
covers file existence, nonzero size, readability and unfinished filenames; the
resource installers retain their format, assembly and index checks.
hg19/GRCh37 uses SpliceAI v1; hg38/GRCh38 uses SpliceAI2. The
installer detects supported input format from the decompressed header and adapts
SpliceAI2 TSV input for the official VEP plugin. Never supply `.partial` files.

Before a **pipeline run**, separately fill `input_vcf`, `ped_file`, `output_dir`
and hub paths, plus any selected cohort/PP1/control inputs. These are not first-install
requirements. Preserve the configured scientific thresholds and family metadata;
missing optional data are not replaced with guessed values.

## Running PriVA

The recommended way to run the pipeline is using Snakemake, which manages dependencies and execution order.

**Prerequisites:**

1. Activate Conda Environment: Ensure the `priva_acmg` conda environment is activated:
   ```bash
   conda activate priva_acmg
   ```
2. Configuration File: Have your `config.yaml` file ready with all necessary paths and parameters correctly specified.

**Running the Pipeline:**

Navigate to your working directory (or specify absolute paths in the command).
The config file can be duplicated and customized for batch-specific runs.

**Basic Command:**

```bash
snakemake --snakefile /path/to/your/PriVA/Snakefile \
          --configfile /path/to/your/config.yaml \
          --cores <number_of_cores>
```

* `--snakefile`: Path to the pipeline's `Snakefile`.
* `--configfile`: Path to your specific configuration file.
* `--cores` (or `-j`): Maximum number of CPU cores Snakemake can use.

**Bash Function Wrapper Example:**

You can wrap the Snakemake command in a bash function for easier execution. Add this function to your `.bashrc` or `.bash_profile`, or save it in a script file and source it.

```bash
run_priva_pipeline() {
    # --- Configuration ---
    local snakefile_path="/path/to/your/PriVA/Snakefile" # CHANGE THIS to the actual path
    local default_config="/path/to/your/default_config.yaml" # Optional: Set a default config
    local default_cores=8
    local default_log_dir="snakemake_logs"

    # --- Argument Parsing ---
    local config_file="${1:-$default_config}"
    local num_cores="${2:-$default_cores}"
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    local log_dir="${3:-$default_log_dir}"
    local log_file="${log_dir}/priva_run_${timestamp}.log"

    # --- Input Validation ---
    if [[ ! -f "$snakefile_path" ]]; then
        echo "ERROR: Snakefile not found at $snakefile_path"
        return 1
    fi
    if [[ -z "$config_file" ]] || [[ ! -f "$config_file" ]]; then
        echo "ERROR: Config file not found or not specified: '$config_file'"
        echo "Usage: run_priva_pipeline <config_file_path> [num_cores] [log_directory]"
        return 1
    fi
    if ! [[ "$num_cores" =~ ^[0-9]+$ ]] || [[ "$num_cores" -lt 1 ]]; then
        echo "ERROR: Number of cores must be a positive integer: '$num_cores'"
        return 1
    fi

    # --- Environment Check ---
    if [[ -z "$CONDA_DEFAULT_ENV" ]] || [[ "$CONDA_DEFAULT_ENV" != "priva_acmg" ]]; then
        echo "WARNING: Conda environment 'priva_acmg' might not be active. Attempting to activate..."
        conda activate priva_acmg
        if [[ $? -ne 0 ]]; then
            echo "ERROR: Failed to activate conda environment 'priva_acmg'. Please activate it manually."
            return 1
        fi
        echo "Activated conda environment 'priva_acmg'."
    fi

    # --- Execution ---
    echo "Starting PriVA pipeline..."
    echo "Snakefile: $snakefile_path"
    echo "Config File: $config_file"
    echo "Cores: $num_cores"
    mkdir -p "$log_dir"
    echo "Log File: $log_file"
    echo "-----------------------------------------"

    nohup snakemake --snakefile "$snakefile_path" \
                    --configfile "$config_file" \
                    --cores "$num_cores" \
                    --printshellcmds \
                    --rerun-incomplete \
                    --verbose \
                    > "$log_file" 2>&1 &

    local job_pid=$!
    echo "Pipeline submitted in the background. PID: $job_pid"
    echo "Monitor progress with: tail -f $log_file"
    echo "-----------------------------------------"
}

# --- How to use the function ---
# 1. Source the script or add the function to your .bashrc/.bash_profile
# 2. Run it:
#    run_priva_pipeline /path/to/my_experiment_config.yaml 16 /path/to/my/logs
#    # Or using defaults:
#    run_priva_pipeline /path/to/my_experiment_config.yaml
```


**Explanation of the Function:**

1. **Configuration:** Set the path to your `Snakefile` and optionally a default config file path.
2. **Arguments:** Takes the config file path (required), number of cores (optional, defaults to 8), and log directory (optional, defaults to `snakemake_logs`) as arguments.
3. **Validation:** Checks if the Snakefile and config file exist and if the core count is valid.
4. **Environment Check:** Verifies if the `priva_acmg` conda environment is active and tries to activate it if not.
5. **Execution:** Creates a log directory, constructs the `snakemake` command with useful options (`--printshellcmds`, `--verbose`, `--rerun-incomplete`), and runs it in the background using `nohup`. It prints the PID and log file location for monitoring.

**Comprehensive Example Command (from original README):**

This shows a full command often used for running on a cluster or server:

```bash
# Activate environment first!
# conda activate priva_acmg

nohup snakemake --snakefile /path/to/PriVA/Snakefile \
                --cores 50 \
                --configfile /path/to/config.yaml \
                --printshellcmds \
                --verbose \
                > /path/to/snakemake_pipeline.log 2>&1 &
```

* `nohup ... &`: Runs the command in the background, detached from the terminal, logging output to the specified file.
* `--printshellcmds`: Shows the actual shell commands being executed by Snakemake.
* `--verbose`: Provides more detailed output from Snakemake.

For family-specific analysis, ensure your pedigree file (`ped_file` in `config.yaml`) is properly formatted and specified in the config file.


## Annotation Resources

The pipeline integrates multiple annotation resources to comprehensively evaluate variant pathogenicity:

1. **Splicing disruption effect**: 
   - SpliceAI (VEP_plugin)
   - SpliceVault (VEP_Plugin)

2. **5-UTR uORF disruption effect**: 
   - UTRannotator (VEP_plugin)

3. **General Deleterious effect prediction**: 
   - CADD (original Python calculation runtime)
   - PrimateAI (VEP_plugin)
   - AlphaMissense (VEP_plugin)
   - LOFTEE (VEP_plugin)

4. **Haplo-insufficiency**: 
   - LOEUF (VEP_plugin)
   - AlphaMissense_mean_score_per_gene

5. **Transcript disruption effect prediction**: 
   - VEP

6. **Clinical variants**: 
   - ClinVar VCF (bcftools annotate)

7. **Population-wise allele frequency + number of homozygous carriers**: 
   - gnomADv4 VCFs (bcftools annotate)

8. **Conservation**: 
   - Conservation (VEP_plugin)


## Pipeline Structure

The pipeline consists of three main steps:

1. **Annotation**: Annotates variants with multiple resources (Rule `annotate_variants`)
   - VEP annotation with plugins
   - CADD scoring
   - gnomAD frequency annotation
   - ClinVar annotation

2. **Filtration**: Filters variants based on family information (if available) (Rules `filter_variants_per_family` / `filter_variants_no_family`)
   - Family-specific filtering when pedigree is provided
   - General filtering without family information

3. **Prioritization**: Prioritizes variants according to ACMG guidelines (Rules `prioritize_variants_per_family` / `prioritize_variants_no_family`)
   - Generates TSV output with prioritized variants
   - Produces ACMG-specific reports


## Output Files

The pipeline produces several output files in the directory specified by `output_dir` in your `config.yaml`:

1. Annotated VCF files (`<input_base>.anno.vcf.gz`)
2. CADD score files (`<input_base>.anno.cadd.tsv`)
3. Filtered VCF files (`<input_base>.anno.filtered.vcf.gz` or `<input_base>.anno.<family>.filtered.vcf.gz`)
4. Prioritized variant tables (`<input_base>.anno.filtered.tsv` or `<input_base>.anno.<family>.filtered.tsv`)
5. ACMG-classified variant tables (`<input_base>.anno.filtered.acmg.tsv` or `<input_base>.anno.<family>.filtered.acmg.tsv`)

Where `<input_base>` is derived from your `input_vcf` filename, and `<family>` corresponds to family IDs found in your `ped_file`.

## Workflow Modes

The pipeline automatically adapts based on the presence of the `ped_file` specified in the `config.yaml`:

- **With family information**: If `ped_file` exists, uses pedigree data to perform family-specific filtering and prioritization for each family listed.
- **Without family information**: If `ped_file` does not exist or is not specified, performs general filtering and prioritization on the entire VCF.

## License

PriVA is licensed under the Apache License 2.0 with Commons Clause. This means:
- You can use, modify, and share it freely for academic and clinical research.
- Commercial use is prohibited without my permission, including selling PriVA as a product or involving PriVA as a part of charging service.

## Contact

Maintainer: yangyxt@hku.hk, yangyxt@gmail.com
