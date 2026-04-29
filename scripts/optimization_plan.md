# PriVA Performance Optimization Plan: VCF→TSV Skip + PS1/PM5 Speedup

## Problem

The ClinVar hg38 re-run took ~3h47m. Two stages dominate:

| Stage | Duration | % of total |
|-------|----------|-----------|
| VCF→TSV conversion (`combine_annotations.py`) | 80 min | 35% |
| PS1/PM5 ClinVar lookup (`check_aa_pathogenic` via `mp.Pool.starmap`) | 105 min | 46% |
| Everything else | 42 min | 19% |

## Fix 1: Skip VCF→TSV when TSV is solid (shell script)

### Current behavior (line 31-35 of `prioritization_vcf_per_fam.sh`)

```
[[ -f ${output_tab} ]] &&
[[ ${output_tab} -nt ${input_vcf} ]] &&
[[ ${output_tab} -nt ${SCRIPT_DIR}/combine_annotations.py ]] &&
{ log "up to date, skip"; return 0; }
```

This checks: TSV exists, newer than VCF, newer than the Python script.

### Problem

When the ACMG assignment step (`acmg_criteria_assign.py`) crashes mid-run, it overwrites the
`.filtered.tsv` with a partial file (104 columns instead of 123+). On the next run, the shell
script sees the TSV is newer than the VCF and skips conversion — but the TSV is incomplete.

Conversely, when only `acmg_criteria_assign.py` changed (not `combine_annotations.py`), the
VCF→TSV conversion runs unnecessarily because the TSV timestamp was updated by the ACMG step.

### Fix

Add a column-count check: the TSV must have the expected annotation columns from
`combine_annotations.py` (specifically `vep_consq_lof` which is the last column added by
the splicing/UTR interpretation steps). If the column exists and the file is newer than
both the VCF and the script, skip.

```bash
# In prepare_combined_tab():
[[ -f ${output_tab} ]] &&
[[ ${output_tab} -nt ${input_vcf} ]] &&
[[ ${output_tab} -nt ${SCRIPT_DIR}/combine_annotations.py ]] &&
check_table_column ${output_tab} "vep_consq_lof" &&
{ log "up to date, skip"; return 0; }
```

`check_table_column` already exists in the codebase (used by `interpret_splicing_annotations`).

### Location

`/paedyl01/disk1/yangyxt/PriVA/scripts/prioritization_vcf_per_fam.sh` line 31-35

---

## Fix 2: PS1/PM5 ClinVar lookup optimization

### Current workflow (105 min for 3.4M rows)

```
PS1_PM5_criteria()
  │
  ├─ Load clinvar_aa_dict.pkl (transcript → position → {HGVSp: {CLNSIG, CLNREVSTAT}})
  ├─ Load clinvar_splice_dict.pkl
  │
  ├─ df.to_dict('records')          ← 3.4M rows → 3.4M dicts (slow, high memory)
  │
  ├─ Build args list:               ← 3.4M tuples, each containing a full row dict
  │   for each row:
  │     (row_dict, clinvar_aa_dict[transcript], clinvar_splice_dict[transcript], pvs1, patho, status)
  │
  └─ mp.Pool.starmap(check_aa_pathogenic, args, chunksize=...)
       │
       └─ Per row: dict lookup by Feature → Protein_position → HGVSp comparison
          + check_splice_pathogenic fallback
```

### Bottlenecks

1. **`df.to_dict('records')`**: Converts 3.4M rows × 100+ columns into Python dicts.
   Each dict copies every column value. For 3.4M rows this is ~10 GB of Python objects.

2. **`mp.Pool.starmap` with row dicts**: Each row dict must be pickled to send to worker
   processes. 3.4M pickle round-trips dominate IPC overhead.

3. **Per-row function call overhead**: `check_aa_pathogenic` does simple dict lookups but
   is called 3.4M times through multiprocessing.

### Proposed optimization: vectorized pandas approach

Replace the per-row multiprocessing with vectorized pandas operations:

```python
def PS1_PM5_criteria_vectorized(df, clinvar_aa_dict_pkl, clinvar_splice_dict_pkl,
                                 ps3_clinvar_patho, pvs1_criteria, threads=10):
    clinvar_aa_dict = pickle.load(...)
    clinvar_splice_dict = pickle.load(...)

    # Step 1: Build a flat lookup table from clinvar_aa_dict
    # {(transcript, protein_position, HGVSp): (is_pathogenic, review_stars)}
    # This is a one-time O(ClinVar_size) operation
    clinvar_flat = {}
    for transcript, positions in clinvar_aa_dict.items():
        for pos, variants in positions.items():
            for hgvsp, entry in variants.items():
                is_patho = any("Pathogenic" in sig and status.get(rev, 0) >= 2
                              for sig, rev in zip(entry['CLNSIG'], entry['CLNREVSTAT']))
                if is_patho:
                    clinvar_flat[(transcript, pos, hgvsp)] = True

    # Step 2: Vectorized PS1 check — exact HGVSp match
    df['_lookup_key'] = list(zip(df['Feature'], df['Protein_position'], df['HGVSp']))
    ps1_mask = df['_lookup_key'].isin(clinvar_flat)

    # Step 3: Vectorized PM5 check — same (transcript, position) but different HGVSp
    clinvar_positions = {(t, p) for t, p, _ in clinvar_flat}
    df['_pos_key'] = list(zip(df['Feature'], df['Protein_position']))
    same_position = df['_pos_key'].isin(clinvar_positions)
    pm5_mask = same_position & ~ps1_mask & ~df['Consequence'].str.contains('synonymous')

    # Step 4: Splice fallback for rows without HGVSp
    # (keep multiprocessing only for the small splice subset)

    df.drop(columns=['_lookup_key', '_pos_key'], inplace=True)
    ...
```

### Expected speedup

- Eliminates `to_dict('records')` (saves ~2 min + memory)
- Eliminates 3.4M pickle round-trips (saves ~90 min)
- Flat dict lookup via `isin()` is O(n) with hash set, not O(n × workers × pickle)
- Splice fallback only runs on the small subset without HGVSp (~5% of rows)

**Estimated: 105 min → ~5 min**

### Location

`/paedyl01/disk1/yangyxt/PriVA/scripts/acmg_criteria_assign.py` lines 983-1049

### Risk

- Must verify that the vectorized PS1/PM5 logic produces identical results to the per-row version
- The `PS1_PM5` coexistence case (both PS1 and PM5 from different ClinVar entries at same position)
  needs careful handling in the vectorized version
- Splice pathogenic fallback still needs per-row logic but only for the small non-HGVSp subset
