# Nextflow and Rust migration

The target combines one official VEP annotation pass, complete CADD feature
enrichment and the original CADD Python scoring path. PriVA-owned functions move
to one Rust/Polars binary under exact source-order contracts. Shell adapters and
mature scientific tools remain in use.

```mermaid
flowchart LR
  I[VCF + fixed scientific profile] --> V[One official VEP run]
  V --> P[PriVA transcript view]
  V --> C[CADD transcript view + mature feature tools]
  C --> A[Complete checked annotation artifact]
  A --> S[Pinned cadd-rs interface / original Python raw + PHRED]
  P --> R[PriVA interpretation / staged Rust modules]
  S --> R
  R --> G[Exact ACMG + ordering acceptance]
```

The current compatibility workflow remains gated for raw annotation. A tested
stub route is not a complete production shared-VEP adapter. Do not activate it
until resource, transcript-selection, model and acceptance gates pass.

<a id="spliceai2-preparation"></a>

## SpliceAI2 GRCh38 preparation and manifest

This route prepares an annotation-only view for the retained official SpliceAI
VEP plugin. It preserves the native 84-column tables and emits top0 scores with
adapter-derived gene context. It does not activate PriVA ACMG interpretation or
replace CADD’s legacy masked-v1 features. GRCh37/hg19 keeps the original v1 route.

```text
Completed native score shards ────────────────┐
Official scored regions + exact GENCODE48 ── verify ── Rust POS/strand join
                                                     ↓
                        BGZF VCF + tabix + preparation receipt
                                                     ↓
                          separate indexed-VEP / scientific activation gates
```

Run preparation inside a suitable Slurm allocation with the pinned `priva`
executable and `gzip`, `bgzip`, `tabix` on PATH. The driver uses Python's standard
library; installer config registration additionally uses the existing PyYAML
runtime. It performs no download or model inference. The output location must be
new: existing prepared outputs are not overwritten.

```bash
# Registration writes only namespaced spliceai2_* preparation metadata.
bash scripts/install_utils.sh SpliceAI2_install \
  /absolute/grch38-preparation.config.yaml \
  /absolute/spliceai2.downloads-complete.json \
  /absolute/new/SpliceAI2/GRCh38 \
  /absolute/pinned/priva

# Standalone preparation; does not edit a PriVA config.
python scripts/spliceai2_to_vep_vcf.py \
  --manifest /absolute/spliceai2.downloads-complete.json \
  --assembly GRCh38 --priva-bin /absolute/pinned/priva \
  --output-dir /absolute/new/SpliceAI2/GRCh38 \
  --receipt /absolute/new/SpliceAI2/GRCh38/adapter.ready.json \
  --threads 2
```

The separate preparation config can select:

```yaml
assembly: GRCh38
spliceai_profile: v2
spliceai2_manifest: /absolute/spliceai2.downloads-complete.json
spliceai2_converter_bin: /absolute/pinned/priva
```

Do not change an active hg19 config for this operation. Fresh GRCh38 resolves to
v2; existing generic GRCh38 SpliceAI paths require explicit `spliceai_profile`
selection. Explicit `v1` preserves the reference route. Requesting v2 on hg19
fails without fallback. After successful v2 preparation, the assembly-profile
installer deliberately returns status **2**, stopping scientific activation.
Calling `SpliceAI2_install` directly can succeed in preparation/registration;
that success still does not qualify downstream analysis.

The acquisition owner supplies a JSON manifest with these exact fields:

| Field | Required type/value |
| --- | --- |
| `schema_version` | Integer `1` |
| `resource`, `version`, `assembly` | Strings `SpliceAI2`, `2.0`, `GRCh38` |
| `source_revision` | `cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9` |
| `status`, `confirmed_by` | `downloads_complete`; nonempty owner confirmation string |
| `files` | Array of 47 unique shard objects; each has `kind` (`snv` or `indel`), `chrom`, absolute `path`, positive integer `bytes`, lowercase 64-hex `sha256` |
| `scored_regions` | Object with absolute `path`, positive integer `bytes`, lowercase 64-hex `sha256` |
| `gencode_gtf` | Same file identity object, plus `release: GENCODE48` |
| `vep_spliceai_plugin` | Same file identity object, with pinned plugin SHA256 below |

All 24 SNV chromosomes (`chr1`–`chr22`, `chrX`, `chrY`) and 23 indel chromosomes
(the same set without `chrY`) are mandatory for full preparation. Source basenames
are exactly `SpliceAI2_GRCh38_{snv,indel}_chr*.tsv.gz`. Partial names/targets,
downloading inventories, missing/extra/duplicate shards and wrong hashes fail.
The scored-region parser requires exactly
`gene_id,chrom,gene_start,gene_end,strand` as tab-separated columns; real publisher
schema/interval qualification remains a required acquisition-time check.

For each region, exact GENCODE gene ID, contig, bounds, strand and protein-coding
biotype must agree. No version stripping or symbol aliases are used. One global
0-based-half-open or 1-based-closed source convention must explain every region;
the resulting context is 1-based closed. Same-chromosome symbol-to-multiple-ID
ambiguities fail, including nonoverlapping/opposite-strand cases; cross-chromosome
symbols do not collide in the plugin lookup.

Rust validates every native column and streams ascending positions. It attaches
each source POS/model strand to all containing verified genes, emits symbols in
lexical order, and preserves source POS/REF/ALT and selected numerical strings.
The projected INFO contains only one `SpliceAI=ALT|SYMBOL|eight-values` entry per
VCF line, ordered AG/AL/DG/DL scores then AG/AL/DG/DL distances. It does not apply
masking, normalization, strand flipping, cutoffs or score recalculation. The
native top3/junction/reference/alternate/summary fields remain in the originals.
Ref/alt probabilities are checked in [0,1], directional deltas in [-1,1]; this
retains the original signed top-k calculation and literal `-0.0` values.

The pinned model takes sequence, strand and allele lengths, not gene identity.
The attached gene label is derived context, not gene-conditioned inference; a
predicted splice site may lie outside the containing gene. V2 metadata thresholds
0.1/0.25/0.5 are not ACMG evidence strengths. Missing coverage is not zero and does
not trigger a v1 fallback.

Verified GTF acquisition: job **1514725** completed successfully. File
`/stornext/snfs170/next-gen/scratch/yangyxt/caches/gene_annotations/gencode.v48.annotation.gtf.gz`
is **58,971,717 bytes**, SHA256
`37a298c2b57e41988cfe12667134106523c69c24a1807347679d19e88aa1f4a9`.
Publisher MD5 `3eb94f7f200b42eb85d56b3afae5b0e7`, gzip integrity and release/header
checks passed. Resource receipt:
`yangyxt_dt/projects/PriVA_upgrade/deployment/spliceai2-gencode48.ready.json`.
This does not establish scored-region agreement or complete cache conversion.
Sources: [GENCODE v48](https://www.gencodegenes.org/human/release_48.html),
[publisher MD5 inventory](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/MD5SUMS).

The installed official plugin is release113 commit
`10054085053745b62ab20129f13cd02f06f7e183`, SHA256
`500a534c403bcbcc15753718768f0a1f352cc6b93f9119e5eba25745f7205b9a`.
The separately inspected SpliceAI2 source is commit
`1dcd2a1882a8223095e7529982ae70bcb658340a`; it is not established as the code
revision that produced the pinned dataset. Original and derived score resources
retain [Illumina's terms](https://github.com/Illumina/SpliceAI2/blob/1dcd2a1882a8223095e7529982ae70bcb658340a/LICENSE),
including redistribution/public-tool restrictions; this route does not publish
resources or enable a public lookup service.

Successful full preparation publishes
`SpliceAI2_GRCh38.{snv,indel}.vep.vcf.gz{,.tbi}` and
`SpliceAI2_GRCh38.gene_context.tsv`, followed by `adapter.ready.json` last.
The receipt records source/driver/converter/tool identities, row counts, mapping
convention, hashes and `status: prepared_annotation_only`,
`coverage: complete_published`, `scientific_qualified: false`. The latter coverage
means the published catalogue, not every possible indel or genome region.
Registration verifies output/index hashes and writes only
`spliceai2_{snv_prescore,indel_prescore,adapter_receipt,status,version,assembly}`;
active generic SpliceAI paths and unrelated config values are preserved.

For bounded fixtures only, the direct driver accepts `--validation-subset` with
at least one SNV and one indel shard. It emits
`status: validated_subset_annotation_only`, `coverage: validation_subset`;
the full installer rejects this as an installation receipt. It is not a way to
activate an incomplete cache. Converter/source/compression/index failures prevent
a completion receipt; PID-only cancellation behavior is not qualified here.

The [official-plugin harness](spliceai2-plugin-audit.md) passed 11 subtests / 38
assertions using real plugin/Ensembl allele matching and synthetic candidates.
Its tabix and light VEP objects are fixtures; no reference, real index or full VEP
execution is exercised. Rust source `428ce4c` (integrated as `c432668`) passed
30 projection plus 9 existing tests. The combined Slurm run passed **165 tests**
with zero failures/skips in 85.702 seconds: 119 existing/route, 14 runtime guards,
21 driver and 11 deployment tests. The driver suite uses real Rust, bgzip and
tabix with synthetic scores/context/plugin identity, including full 47-shard
concatenation and CRC/mutation failures. Log
`yangyxt_dt/projects/PriVA_upgrade/logs/spliceai2-install-integration-20261015.log`
records integration `ad9af5d` and exact driver/test content hashes; those
uncommitted-file hashes, not a later checkout, identify the tested driver.
Bash/YAML/diff checks also passed. No completed real 47-shard conversion,
v2 ACMG/CADD activation or supplied-family parity is claimed.

Nextflow activation guard: shared bin/priva_task.py validation rejects explicit v2 scientific configuration in the launch/task path, including preannotated import/resume, and rejects SpliceAI2-marked input headers. Namespaced prepared-resource metadata remains allowed. A separate adapter suite passed 20 tests (16 existing + 4 new) under Slurm, covering launch and all three task phases before output; this is separate from the 165 installer/deployment tests. Unmarked external preannotated artifacts still require pinned scientific-profile provenance: predictor version cannot be inferred from numeric scores.

## Frozen source and branch

- Branch: `feature/nextflow-rust-migration`.
- Reference commit: `09223ef5a8b0ceb4c6d61055c6eefe62b51ba903`.
- Reference incorporates the eight previously reviewed local changes on top of
  `810e60165a663bd2e1223809af9f32d664a9d543`.
- Development worktree: `/stornext/snfs170/next-gen/scratch/yangyxt/git_repos/PriVA-nextflow`.
- The original `PriVA` checkout keeps its original uncommitted changes.

The branch and worktree serve complementary purposes: the branch records the
migration; the worktree keeps its files separate from the current checkout.
The HTML design is in `projects/PriVA_upgrade/migration_plan.html` under the
approved `yangyxt` root. See [Rust scope and parity](rust-migration.md).

## Current acceptance gate

The user will supply the family VCF/PED, assembly, working configuration and
existing result or its producing software environment. Until then, no real-family
analysis or identical ACMG/ranking result is claimed. A historic result qualifies
as a reference only when its source, software, parameters and resource identities
are established; otherwise regenerate with the frozen source.

Use the same input bytes, phenotype/sex/proband context, original cohort/control
context, resource releases, numerical settings and ClinGen override on both sides.
Select a family after cohort annotation when the original cohort defines control
frequencies. Do not silently reduce the cohort before annotation. An invalid PED
is an error. Transcript multiplicity and tied-row order are part of acceptance.

Compare both final tables using `bin/priva-compare --help`. Required equality
includes `sort_index`, `pext_sort_index` when present, `variant_rank`, variant and
transcript row order, every ACMG criterion/strength, combined judgement and score.
The comparator also checks all other serialized fields; it does not hide numeric
differences behind a tolerance or silently deduplicate records.

## Storage

Code, installation scripts and narrative remain under `yangyxt`. Inputs, outputs,
work files, build products and retained logs are under the project-scoped dataset:

```
/stornext/snfs170/next-gen/scratch/yangyxt_dt/projects/PriVA_upgrade
  inputs/       selected family and provenance
  results/      reference/candidate results and comparison reports
  validation/   synthetic tests, distinct from scientific acceptance
  logs/         completed setup and validation logs
  work/         ignored Nextflow task work, retained for resume
  build/        ignored Rust build products
  tmp/          ignored temporary files
```

Save named completed results with DataLad and attach project/purpose/status
git-annex metadata after all writers finish. Keep Nextflow controller state and
task work together for resume. Neither an output directory nor a Git commit
contains all state needed to resume an interrupted pipeline.

Save/materialize analysis inputs with DataLad before the first production
launch and keep their storage identity stable. Annexing inputs afterwards can
change resolved paths or modification times and invalidate the recorded run
identity. The synthetic resume assertion was made before evidence archival.

Existing reference/cache directories can remain under `yangyxt/caches`; analysis
must not rebuild or alter them. A host-side metadata check through Slurm job
1514363 on 2026-10-15 found only README files
in `caches/{CADD,ClinVar,gnomAD,VEP}` and no `yangyxt_dt/caches`. The repository
config additionally refers to nonexistent `VEP_cache`, `VEP_plugin` and
`VEP_plugin_cache` folders and old `/paedyl01` paths. Those values are not a
deployed resource bundle. Supply the actual cache set before a real run. This
finding applies to these locations, not every possible PriVA installation.
The agent's filename-filtered directory view alone cannot establish missing
genomic data. Run data-facing commands in authorized Slurm jobs with the existing
filter policy intact; the host check is recorded in `logs/cache-path-check.log`.

## Driver installation

`scripts/bootstrap_nextflow_runtime.sh /absolute/new/prefix` installs the official
Nextflow **25.10.4** standalone distribution and Temurin **21.0.12.1+1** JRE for
Linux x86_64, verifying the pinned upstream SHA-256 digests before extraction.
It refuses an existing prefix and makes no shell-startup changes. An interrupted
prefix is incomplete; inspect it and choose a new prefix for a retry.

The installed local launcher is:

```
/stornext/snfs170/next-gen/scratch/yangyxt/softwares/priva-nextflow-25.10.4/bin/nextflow
```

It selects strict parser v2, a private `NXF_HOME`, and a bounded driver JVM.
This installs the workflow driver only. The legacy scientific environment,
VEP/plugins, CADD and resource deployment require their own verified identity.
The old Dockerfile and unpinned `acmg_conda_rough.yml` are not a reproducibility
lock and have not been used to claim production readiness. Freeze the actual
reference environment before solving or upgrading its scientific dependencies.

See the official [installation documentation](https://docs.seqera.io/nextflow/install).
An offline run must also have all required plugins, environments/images and
resources available. `NXF_OFFLINE=true` alone cannot stop subprocess downloads.

## Execution contract

Run `nextflow run /absolute/PriVA-nextflow --help` for current parameters.
The compatibility workflow accepts BGZF VCF plus index, a proband-first PED,
one selected family, assembly, legacy YAML, resource manifest, software identity
and a new output directory. Preannotated mode also requires the matching CADD
annotation TSV; an arbitrary score-only table is insufficient.

**Current deployment gate:** real raw-VCF annotation is deliberately refused
until CADD is qualified with enforced read-only shared resources and fully
prepared software. The unchanged official launcher may otherwise create Conda
environments under the shared CADD tree. Raw stubs test wiring only. The first
real compatibility run can start from the original annotated cohort VCF and
paired CADD TSV when its complete prioritization caches and software are ready.
This gate is an outstanding implementation/deployment step, not evidence that
raw annotation has passed.

Each task owns mutable VCF/index/config/output/temp files. Reference resources
are read from their existing locations. Paths unsafe for the unchanged shell
scripts are rejected. Source content, config and resource identity are explicit
task inputs. Resource metadata inventory is distinct from content verification;
the final family gate needs a frozen, identified cache set.

VCF hub policy must be identical for reference and candidate. The legacy VCF hub
can carry cohort-specific control frequencies, so a hub from another cohort is
not a valid shortcut. Shared hub writes are isolated from the original caches.
Do not switch only one side to a different hub policy and call it orchestration
parity. All output publication is downstream of task success.

The `local` profile is for an existing allocation or a workstation. `slurm` is
for a verified driver that can submit child jobs, with explicit site resources.
The `apptainer` profile needs an already prepared, identified image; it does not
make the current container recipe deployable. `test` is for synthetic wiring.
Use `-stub-run` only for wiring validation, never scientific acceptance.

Launch-time resource checks run even on a fully cached `-resume`. The output
directory records the input/source/resource/settings identity and rejects a
different analysis, avoiding stale published tables. Per-stage provenance and
the validation contract are published alongside uniquely named run reports.
The software ID identifies a caller-supplied frozen environment/image; it is not
an automatic environment lock. Full production software freezing remains part
of the deployment gate.

## Initial validation evidence

On 2026-10-15, validation passed 15 adapter tests, 14 exact-comparator tests,
60 Rust/Python terminal-order comparisons plus five rejected-input cases, and
the 27 retained environment-bootstrap tests. The Rust oracle executed the
terminal sorting block from the frozen scientific source using pandas 3.0.6;
it did not execute or reimplement all upstream scoring/ACMG decisions.

The pinned Nextflow runtime completed raw and preannotated synthetic workflows,
reused all four tasks with `-resume`, preserved an additional cohort sample,
and rejected output reuse with a changed software identity. Successful evidence
is in the data project's `validation/nextflow-integration-v3/{validation,logs}`;
its `work/` tree remains available for resume. These artificial fixtures are
not a real-family test. GitHub Actions later passed the expanded contract, Rust and both Nextflow
stub suites at `56b7762` ([run 38828310089](https://github.com/yangyxt/PriVA/actions/runs/38828310089)).
The integrated source-digest ordering regression was fixed and 16 adapter checks
passed. These CI checks contain no scientific annotation or real-family analysis.

The driver must run in a valid allocation/approved service. On HGSC analysis
nodes, work exceeding any of 2 CPUs, 4 GB or 2 minutes belongs in Slurm. Local
tasks and the JVM share the parent allocation; set process resources to leave
capacity for the driver. A 12 GB or greater CADD task does not fit the current
10 GB interactive allocation and needs its own adequately sized allocation.

## Shared VEP and original CADD calculation

CADD scripts **1.7.2** are pinned at
`1fead0c1e02d3f0a8599bdb04e6b3a2dd70b4866`; the score model is **1.7**.
Production annotation reuses PriVA's official VEP run, then mature upstream
ESM/RegSeq/MMSplice and non-VEP track enrichment produce the complete feature
artifact. Validate assembly, schema, missingness, all required resources, exact
FeatureID joins and order. Qualify CADD's transcript selection without removing
PriVA transcript rows. The separate original-VEP route is a reference-only
qualification tool. VEP113 peptide/GTF compatibility requires evidence; directory
names alone cannot establish it.

`cadd.nf -c conf/cadd.config --cadd_mode full|annotate|calculate` exposes separate
stages. See its `--help` and pinned profile templates under `conf/cadd/` for the
implemented contract. Profiles must be configured and qualified; templates do
not constitute installed resources. Calculate-only consumes a completed artifact
and requires no annotation cache installation. A validated bridge to the pinned
standalone `cadd-rs` interface must preserve both receipt identities and the
identical table; their distinct receipt schemas are not interchangeable.

The standalone scorer must run all four original stages: transformation,
prediction, consequence hierarchy and PHRED conversion. The previous standalone
pin `79d767a` included native hierarchy/PHRED and is a historical experiment;
candidate `4194be3` runs all four original Python stages and passes contract
tests. Both 18-row upstream assembly fixtures match direct original Python
exactly through features, raw scores, hierarchy and PHRED, including full rows
and order. Both real Nextflow calculation routes and one cached resume pass.
Shared annotation and broad scientific qualification remain pending.
Keep the original numerical libraries, other scientific packages and model
assets pinned. An optimized Python backend must be explicitly selected.

| Qualification | Required evidence |
| --- | --- |
| Profiling | Per-stage and total elapsed/CPU time, peak RSS, inputs/source/environment/resource/profile hashes, fixed CPU/thread limits. Repeat reference and candidate fairly; report setup and warm/cold effects. |
| Measured changes | NumPy/pandas vectorization, batching or independent parallel work only for an observed bottleneck; source operation/mutation order and failure behavior remain fixed. |
| Exactness | Raw and PHRED scores, full FeatureID rows, cell values/order, duplicates, missingness, intermediate numerics and failures; every mismatch reported with no tolerance relaxation. |
| Edges | Empty/header-only, missing features, duplicates, transcript tiers/ties, PHRED boundaries, nonfinite/malformed values, failures and whole-variant chunk boundaries. |
| Public benchmark — pending | Later pin a 1000 Genomes sample or agreed dataset; record assembly, release/accession, hashes, selection, sample/variant/transcript coverage and exclusions. Validate both annotation and scoring. |
| Family — pending | User-supplied family with identical fixed resources/profile; exact ACMG criteria, scores and ranking. |

Upstream model fixtures qualify calculation only. They do not satisfy public-data
annotation acceptance or family acceptance. Keep failed optimization experiments
and their mismatch reports; retain the unmodified executable reference.

The nonfinite edge suite found one explicit CLI mismatch: original PHRED code
prints a diagnostic into stdout for raw `nan` and exits zero; the standalone
preserves the output bytes but rejects the malformed table and exits one.
This failure-behavior gate remains failed. Compare original and optimized Python
at the same scorer boundary, retaining stdout/stderr, exit/signal and malformed
outputs. Record artifact validation separately; do not call malformed output a
valid score artifact. Independent profiling of passing model inputs continues.

`CaddSelection.pm` delegates selection to official VEP methods and emits native
ordinals while keeping every PriVA transcript row. `cadd_select_view.py` validates
coverage/context before recovering the CADD view. Its 25 Python and 18 native
method fixture checks establish mechanics only. The policy fragment leaves the
CADD transcript source pool unqualified: Ensembl-only needs RefSeq-only overlap
handling; Ensembl plus RefSeq needs compatible peptide/GTF/model coverage. Neither
is activated by default. PriVA's merged cache and `use_transcript_ref` remain
fixed, and their CADD-derived profile requires separate qualification.

## Current cluster build checkpoint

Installed `priva_acmg` software passed import and executable checks after adding
missing VEP Perl helpers. Its readiness receipt and restorable `environment.conda-explicit.txt` lock are
stored in the data project's `deployment/`. The earlier `environment.explicit.txt`
is a historical URL inventory, not a valid restore file. The pinned resource installers under
`scripts/deploy_cluster_*.py` verify content identity, keep partial downloads
inactive and write stage receipts. `scripts/check_cluster_readiness.py` reports
remaining gaps without equating file presence with scientific qualification.

VEP/cache archives, CADD annotations/prescores and gnomAD inputs are separate
Slurm jobs; completed shared-plugin and HPO/pext stages do not imply full readiness.
The original SpliceAI v1 and PrimateAI collections are acquired under
`yangyxt/caches/{SpliceAI_v1,PrimateAI_v0_2}`. Hg19 SpliceAI import 1514615 and
PrimateAI import 1514620 both completed with exit 0 and valid receipts. The latter
retains the original formatter and exact external source-to-formatted validation,
including every child exit and full decompressed output comparison. The earlier
queued-import checkpoint is superseded; cancelled job 1514618 never ran. Original
hg38 files stay retained and unselected. Main-config promotion is separately
tracked and does not follow automatically from an import receipt.

GENCODE48 is verified for the separate SpliceAI2 preparation route; the final
scored-region TSV remains absent at the latest metadata check, and full v2
conversion is pending. gnomAD derived job 1514601 failed extension validation for
`.vcf.bgz`; its repair is separate. AlphaMissense 1514589 remains held for the
original NM-to-ENST mapping. Hg19 MaveDB policy, derived calibration caches and
shared-VEP/CADD qualification remain explicit gates. Do not duplicate the
owner-managed acquisition or consume partial files. No scientific evidence has
been silently excluded.

The actual original `main_install` was reproduced in Slurm 1514592 with an
isolated copied environment. VEP's update prompt exited zero on batch EOF,
leaving its cache/plugins absent while PriVA continued. The fix disables that
prompt, pins the cache release, retains official Conda VEP/API and checks cache
and plugin outputs. All 28 bootstrap tests pass. Corrected real rerun 1514593
reached the intended VEP113 GRCh37 download, then deliberately timed out (124);
that bounded checkpoint is not a complete installation.

Nine additional installer orchestration regressions pass after correcting the
gnomAD reference argument and all-chromosome scheduling, failure barriers,
Conservation-before-AlphaMissense ordering, Conservation reuse and configured
motif output. Scientific Python calculations remain unchanged.
