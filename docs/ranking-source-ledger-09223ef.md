**Remaining PriVA ranking migration: preserve the source before extending Rust**

The remaining work produces `control_common_index` and `sort_index`, then feeds
the final ordering step. Existing Rust scoring/ordering tests cover bounded
interfaces; they do not establish fidelity for the full Python ranking path.
This is an audit and test plan. No scientific Python or Rust was changed.

```mermaid
flowchart TD
    I[Production or re-sort parser: distinct input contracts] --> R1[R1–R2: copy, protein-coding filter, pedigree]
    R1 --> C[C1–C5: control routes, SciPy CDF and Beta quantile]
    C --> R3[R3–R4: control multiplication and mechanism masks]
    R3 --> R5[R5–R6: dispensable, relevant, then compatibility multiplication]
    R5 --> R7[R7: select PEXT column, parse, power, multiply]
    R7 --> R8[R8: coordinate maximum, six-key sort, first-appearance rank]
    R8 --> O[Ranked table; existing ACMG score and criteria retained]
    R8 -. bounded subset already implemented .-> B[Rust final-order: finite scores and nonmissing identifiers]
```

Solid arrows are execution order. The dotted arrow marks existing limited Rust
coverage, not completion of the preceding steps.

**Pinned evidence**

| Item | Identity / source boundary |
| --- | --- |
| Scientific source | `09223ef5a8b0ceb4c6d61055c6eefe62b51ba903`; line numbers below refer to this commit. |
| Ranking module | [acmg_scoring.py](../scripts/acmg_scoring.py), SHA256 `9111a76e3fc94a064e14d97b535bb1dd2132f7c61836f8e8892477d4c5e69646`; `compute_control_common` L390–471; `sort_and_rank_variants` L474–697. |
| Called helpers | [acmg_variant_mechanism.py](../scripts/acmg_variant_mechanism.py) `_variant_mechanism_masks`, implementation L468–591; [acmg_consequence.py](../scripts/acmg_consequence.py) `coerce_pext_value`, L28–49. |
| Probe runtime | Installed `priva_acmg`: Python 3.12.3, NumPy 1.26.4, pandas 2.2.3, SciPy 1.14.1. Small in-memory probes imported the frozen `softwares/priva-reference-09223ef/scripts` source under the approved `yangyxt` root inside Slurm allocation 1514363. No family data or licensed downloads were read. |
| Rust boundary inspected | `6a7761ce19462f927bf95f0901936736154be93d`: `rust/src/acmg_scoring/*` implements four summary/posterior helpers; `rust/src/cli.rs::final_order` implements terminal grouping/order/rank. |
| Existing qualification | Rust worker reports 19,754 exact scoring rows plus 16 failure cases, and 60 ordering fixtures plus 5 error cases. These are existing bounded claims; this audit does not broaden them. |

**C1–C5: control-common mutation ledger**

This helper mutates and returns the supplied DataFrame; it does not copy it.
Defaults are `cutoff=0.01`, `conf_level=0.999`. The docstring's fixed 0.5 penalty
description is outdated: the implementation grades flagged rows using Beta.

| Step / source | Reads and decisions, in execution order | Mutation / observable behavior |
| --- | --- | --- |
| C1, L400–411 | If none of `control_AC`, `control_AN`, `control_AF`, `control_nhomalt` exists, return early. Otherwise apply `pd.to_numeric(errors='coerce')` separately to all four; valid AN means `AN.fillna(0)>0`. | Early path writes `control_common=False`, then `control_common_index=1.0`. A partially absent schema is **not** converted into four null Series: missing `.get()` values become scalar NaN and later Series methods can raise. |
| C2, L413–425 | If the named proband exists, stringify its sample value, split at `:`, and recognize only GT `1/1`, `1\|1`, or `1` as homozygous. Priority: valid AN + homozygous + nonnull nhomalt; otherwise valid AN + nonnull AC; otherwise valid AN + nonnull AF. | Three mutually exclusive masks retain the original index. AF never substitutes for missing/nonpositive AN. |
| C3, L427–447 | Start false. Homozygous route: `k=int(max(nhomalt,0))`, `n=int(max(round(AN/2),1))`. AC route: `k=int(max(AC,0))`, `n=int(max(round(AN),1))`. Call `binom.cdf(k,n,cutoff)>conf_level`, strictly. AF route compares `AF>cutoff`, strictly. | Boolean assignments are index-selected; SciPy receives positional NumPy arrays. Fractional counts truncate after lower clipping; trial counts use pandas ties-to-even rounding. No upper clipping of counts. |
| C4, L449–462 | Write `control_common`. Independently fill original AC/AN nulls with 0 and cast to float. For common rows with positive AN, compute `q=beta.ppf(0.01, AC+1, AN-AC+1)`; other rows keep q=0. | Beta always uses **original AC and AN**, including homozygous/AF routes; it does not reuse rounded trials or nhomalt. No extra shape validation or repair occurs. |
| C5, L463–471 | For common rows use `clip(1-q,0.01,0.5)`; otherwise 1.0. | Write `control_common_index` using positional `np.where`. NaN Beta results remain NaN. Log flagged count and return the same object. |

The CDF asks whether the observed count is already above almost all mass under
the configured frequency. Beta then supplies a lower frequency bound used to
grade the penalty. These are separate calculations: for AN 100, AC missing,
AF 0.9 and nhomalt missing, the AF route flags the row but Beta uses AC 0, and
the resulting penalty is 0.5. Substituting AF×AN or homozygote counts changes
the current behavior.

**R1–R8: ranking mutation ledger**

| Step / source | Exact order and result |
| --- | --- |
| R1, L492–498 | Copy input. Drop only `haplo_insuf_index` and `zygosity_inheritance_mechanism_index` if present. Initialize `sort_index` from `ACMG_quant_score`, then retain exact `BIOTYPE == 'protein_coding'`. Keep surviving row indices; do not reset or deduplicate here. |
| R2, L501–524 | When both PED and family are non-None, select the first matching affected row (`Phenotype` in string/integer 2), obtain its sample column, normalize Sex, and compute heterozygosity by `.str.count('1')==1` on the **entire sample string**. This differs from C2's GT-only parser. With no PED/family, use an index-aligned false Series. Missing affected row/sample can raise; no alternate proband policy. |
| R3, L527–528 | Run C1–C5, then multiply `sort_index *= control_common_index`. Preserve this first floating multiplication and its NaNs. |
| R4, L537–597 | Read required hub mechanism fields through `_variant_mechanism_masks`. Detect X by uppercase chromosome text and removal of leading `CHR`. Dominant LOF/bare dominant history permits a single alternate allele generally; males at X additionally accept X-recessive or X-dominant LOF/bare histories, females at X accept X-dominant LOF/bare histories. Penalty mask: raw `ACMG_quant_score>0.89` AND R2 heterozygosity AND predicted LOF AND not exact GOF AND no compatible history. Write compatibility 1.0/0.9, but defer multiplication. DN-without-LOF is logged, not an independent penalty. |
| R5, L599–611 | If dispensable list is not None: read with `header=None`, `names=['SYMBOL']`, `comment='#'`; exact set membership gives 0.9/1.0, then multiply. If relevant list is not None: same parsing gives 1.2/1.0, then multiply. No whitespace/case normalization. Explicit None disables that branch; existing unused index columns can remain. |
| R6, L613–616 | Multiply by compatibility 0.9 only on the previously computed R4 mask. Do not combine/reorder the three gene/compatibility multiplications. |
| R7, L621–665 | Use only the first nonempty configured tissue, replacing hyphens with underscores. If that column is absent, try `PEXT_MEAN`; do not try subsequent tissues. Validate parameters only if a column is selected. Map `coerce_pext_value`, fill null with 1, clip 0..1, divide by cutoff, clip 0..1, then apply the exact NumPy power expression. Write `pext_sort_index` for all rows; multiply only rows whose **raw** ACMG score exceeds0.89, regardless of predicted-LOF status. |
| R8, L668–697 | Group by typed `(chrom,pos,ref,alt)` and broadcast maximum `sort_index`. Sort by maximum score↓, **raw ACMG score↓**, chrom↑, pos↑, ref↑, alt↑. Factorize existing `variant_id` by first appearance and add 1; if column absent, factorize colon-joined `astype(str)` coordinates. Drop `max_variant_score`; return sorted rows without resetting indices. |

R7's local calculation is `floor + (1-floor) * power(clip(pext/cutoff,0,1),shape)`
below cutoff, otherwise 1. At defaults, pext 0.02 gives 0.8894427190999916 and a
raw score 0.95 becomes 0.844970583144992 when earlier penalties are 1. The exact
operation order matters; the source does not collapse all multipliers into one
product or round the final sorting index.

**Input contracts and pandas behavior that Rust must retain**

| Boundary | Contract / fidelity consequence |
| --- | --- |
| Ranking library | Required: `ACMG_quant_score`, `BIOTYPE`, coordinate fields, `var_plausible_patho_mechs`, `variant_effect`; `SYMBOL` when a gene list is used; named sample when PED selects a proband. Numeric ACMG scores are assumed, not coerced. `variant_id` is optional. Control schema must exercise both no-columns and partial-column behavior. |
| Production CLI | `acmg_scoring.py` has no direct CLI. [acmg_criteria_assign.py](../scripts/acmg_criteria_assign.py) L196 loads annotation TSV with inferred dtypes/default NA handling and drops duplicate complete rows. PED parsing L217–237 explicitly uses strings, preserves NA-like text, and supports headerless six-column PED. The final scorer receives the **post-assignment** frame, not the initial TSV frame. Ranking options include tissues, cutoff 0.1, floor 0.8, shape 0.5, optional relevant list and bundled dispensable list by default. |
| Production outputs | L533–564 writes the criteria matrix before ranking; then sets posterior/class, ranks, and overwrites the annotation TSV with `index=False`. Protein-coding filtering does not retroactively filter that saved matrix. Compare criteria separately and preserve row association. |
| Re-sort CLI | [resort_step3_only.py](../scripts/resort_step3_only.py) reads existing TSV/PED with inferred dtypes, does not support the production PED fallback, removes 10 existing ranking columns first, and defaults dispensable list to None. It skips missing files, successfully parsed empty frames, or frames without `ACMG_quant_score`, and overwrites the supplied TSV after ranking. Parser errors propagate: zero-byte or whitespace-only inputs raise `pandas.errors.EmptyDataError`. These caller differences must not be silently harmonized. |
| Dtypes / missingness | Preserve pandas `to_numeric` coercion, inferred integer/float/object columns, default NA tokens, native typed coordinate sorting and gene-list membership. A numeric-only chromosome column can sort numerically while a mixed X/autosome column sorts as text. Nullable extension dtypes need their own oracle cases; they are not interchangeable with float NaN or Polars null. |
| Alignment / column order | Preserve index labels and a separate row-position identity, including gaps and duplicates. Series `.loc` arithmetic aligns labels; `.values` and `np.where` are positional. Overwriting an existing column keeps its position; newly created columns append. Optional skipped branches do not generally remove stale columns. |
| Grouping / ranking | Default pandas groupby drops null coordinate keys; max skips null scores and all-null groups retain null maxima. Sorting places NaNs last. Equal six-key rows need source-order evidence. Coordinate grouping and `variant_id` factorization are different identities: missing IDs receive code −1 and therefore rank 0; distinct coordinates may share an ID. Absent-ID stringification must match pandas, including null spellings. |
| Helper parsing | Mechanism tags are case-sensitive, semicolon-delimited exact matches without trimming each tag; missing required fields raise. Optional LOF/GOF/DN scores use numeric coercion and null→0. `coerce_pext_value` parses comma/ampersand segments and Python `max`; despite its “finite” comment it accepts infinities and some NaN spellings. Do not replace it with a finite-only or null-ignoring reduction. |

**Numerical boundary: keep the actual SciPy calls**

PriVA calls only `scipy.stats.binom.cdf` and `scipy.stats.beta.ppf` in C1–C5;
it does not fit a distribution or use a beta-binomial test. SciPy 1.14.1's
binomial implementation floors k and calls compiled `scipy.special` `_binom_cdf`;
its beta quantile calls `_beta_ppf`. The latter uses SciPy's Boost inverse-beta
wrapper and error policy. Public distribution wrappers also apply domain and
support handling. Preserve the exact installed build, not only function names.
See [SciPy 1.14.1 binomial source](https://github.com/scipy/scipy/blob/v1.14.1/scipy/stats/_discrete_distns.py),
[beta source](https://github.com/scipy/scipy/blob/v1.14.1/scipy/stats/_continuous_distns.py),
and [compiled beta wrapper](https://github.com/scipy/scipy/blob/v1.14.1/scipy/special/boost_special_functions.h).

An alternative CDF, survival-tail reformulation, inverse solver, count clipping,
or numerical tolerance is not proposed. Preserve comparison thresholds,
ties-to-even rounding, integer conversion, binary floating values, NaN masks,
and NumPy power/multiplication order. A one-bit difference can alter a tie or a
strict threshold. Keep this numerical stage on the pinned Python/SciPy path
until a source-faithful Rust implementation has exact differential evidence.

**Representative oracle observations, already checked**

These are small probes, not a full ranking qualification. All use the pinned
source/runtime above. Full-ranker examples have one protein-coding variant,
raw score 0.95, empty mechanism history, effect `uncertain`, no PED, and explicit
`dispensable_gene_list=None`.

| Input | Observed Python result | Existing Rust `final-order` implication |
| --- | --- | --- |
| No control columns | `False`, index 1.0; helper returns the same frame object. | No upstream computation exists in Rust yet. |
| Only `control_AF=0.9` | `AttributeError: 'numpy.float64' object has no attribute 'fillna'`. | Capture failure parity; do not silently synthesize missing columns. |
| Only `control_AN=100` | `AttributeError: 'numpy.float64' object has no attribute 'notna'`. | Same unsupported upstream scope. |
| AN 2, AC 4, AF 0.9, nhomalt NaN | Common=True, index NaN; full ranker returns `sort_index=NaN`, `variant_rank=1`. | Rust explicitly rejects nonfinite `sort_index`. This is a full-path gap requiring extension/reference tests, **not legitimate full-ranker equivalence**. |
| AN 100, AC NaN, AF 0.9, nhomalt NaN | Common=True, index 0.5. | Needs the original AF-route/Beta combination. |
| PEXT `NAN&0.02` | Parser NaN→fill 1; sort 0.95, hex `0x1.e666666666666p-1`. | Terminal Rust accepts the finite prepared score but does not compute it. |
| PEXT `0.02&NAN` | Parser 0.02; sort 0.844970583144992, hex `0x1.b09ffbf9612c1p-1`. | Order-sensitive parser is an unported upstream dependency. |
| PEXT `inf&0.2` / `-inf` | Parser returns +∞ / −∞; ranker later clips to 1 / 0. | Parser's nonfinite value is not itself grounds to reject the full Python path. |

Current Rust additionally requires nonmissing `variant_id`/coordinates, positive
integer positions, finite raw/final scores, unique headers, and explicit
`--chrom-type string|integer`. Python's absent/missing-ID branch, null grouping,
nonfinite scores and broader inferred dtypes are outside that validated subset.
Existing bounded results remain useful; they must not be advertised as full
ranking or family parity.

**Differential-test plan and the next independent module**

| Work | Required capture / acceptance |
| --- | --- |
| Freeze the oracle | Record scientific commit and module hashes, exact Conda lock/build identities, SciPy/NumPy native-library identities, parameters, and gene-list bytes/hashes. Capture the actual ranking-entry DataFrame with ordered columns, index labels, dtypes, null masks, row-position IDs and float bit patterns. Retain caller-specific parsing evidence; TSV alone loses that state. |
| Capture C and R checkpoints | Use a separate harness against unchanged frozen source: capture route masks, k/n arrays, CDF outputs, Beta inputs/outputs, every newly assigned column and each sequential `sort_index` state. Use independent deep copies for isolated helper cases; also test helper in-place mutation and ranker input preservation. Record exceptions and their first failing stage. |
| Controls matrix | All 16 present/absent column combinations; bad numeric text; zero/negative/missing AN; fractional/half-integer counts; AC>AN; infinities; all three route collisions/priorities; GT suffixes; values around strict cutoffs; large counts; binary-float neighbors of resulting thresholds. Compare booleans, integer arrays, float bits/NaN masks and failure behavior exactly. |
| Ranking matrix | Filtered rows with gapped/duplicate indices; empty frame; PED first-affected selection and caller dtype differences; full sample strings containing additional 1s; male/female/unknown Sex at X/autosomes; every exact history/effect tag; missing mechanism fields; list None versus empty/missing files and NA-like symbols; overlap of both gene lists; raw-score neighbors of 0.89; tissue fallback; PEXT strings/nulls/infinities/NaN token order; invalid parameters with/without a selected column. |
| Terminal extension | Include group ties/all-null groups, null coordinate keys, string versus numeric chromosomes/positions, missing or absent IDs, repeated IDs across coordinates, stale columns and output column positions. Preserve full row order and every rank; compare against the full Python function, not only already normalized Rust inputs. |
| First independent Rust piece | Port `_variant_mechanism_masks` as its own source-mapped module **after** capturing its exact masks/coercions/errors. It is a Boolean/string/numeric-coercion boundary without CDF/PPF or PEXT power. Qualify it independently, then call it at R4; moving it before C1 changes which error happens first. The scalar `coerce_pext_value` parser is a second separable candidate once order-sensitive NaN cases are frozen. |
| Integration gate | No reassociation, parallel reduction, null-policy improvement or approximate acceptance. Require exact intermediate/output equality on the declared supported domain before extending that domain. Preserve Python for unsupported cases; later compare the supplied family against the same complete cache snapshot, including exact `sort_index`, ordering/rank, criteria matrix and criteria judgments. |

Only the representative in-memory probes above were executed during this audit.
The comprehensive oracle capture, Rust extensions, complete family run and
scientific equivalence checks remain to be performed.
