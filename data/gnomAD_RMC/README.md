# gnomAD Regional Missense Constraint (RMC) Data

## Source

gnomAD v2.1.1 Regional Missense Constraint (RMC) identifies sub-genic regions
with significant depletion (or enrichment) of missense variants relative to a
neutral expectation model.

- **Source**: [gnomAD v2.1.1](https://gnomad.broadinstitute.org/)
- **Publication**: [Karczewski et al., Nature 2020](https://doi.org/10.1038/s41586-020-2308-7)
- **Original assembly**: GRCh37 (hg19) — gnomAD v2.1.1 is based on GRCh37

## Files

| File | Assembly | Origin | Regions | Transcripts | Genes |
|------|----------|--------|---------|-------------|-------|
| `gnomad_v2.1.1_rmc_flat.hg19.tsv` | GRCh37 | Original (from gnomAD) | 17,842 | 5,127 | 5,118 |
| `gnomad_v2.1.1_rmc_flat.hg38.tsv` | GRCh38 | Derived via protein-coordinate transfer | 16,411 | 4,745 | 4,741 |

## Columns

### hg19 (original, all columns)

| Column | Description |
|--------|-------------|
| `transcript` | Ensembl transcript ID (stable, no version; e.g., `ENST00000000233`) |
| `gene_name` | Gene symbol |
| `gene_id` | Ensembl gene ID |
| `contig` | Chromosome (bare integer; hg19-specific) |
| `genomic_start` | Genomic start position (hg19-specific) |
| `genomic_end` | Genomic end position (hg19-specific) |
| `start_aa_raw` | Start amino acid with position (e.g., `Arg97`) |
| `end_aa_raw` | End amino acid with position (e.g., `Asp96`) |
| `obs` | Observed missense variant count |
| `exp` | Expected missense variant count (under neutral model) |
| `oe` | Observed/Expected ratio |
| `chisq` | Chi-squared test statistic |
| `p` | P-value from chi-squared test |

### hg38 (derived, genomic columns dropped)

Same as hg19 but **without** `contig`, `genomic_start`, `genomic_end` (these are
hg19-specific and cannot be directly transferred to hg38).

## Key Metrics for Anchor Classification

The raw `p` column tests whether oe differs from 1, but inflates with region size and is
not directional. **Benjamini-Hochberg FDR correction** must be applied across all RMC
intervals to produce `p_bh` before classification.

```
Constrained (positive anchor):   oe < 0.4  AND  p_bh < 0.05
Indeterminate (exclude):         0.4 <= oe <= 0.8  AND  p_bh < 0.05
Unconstrained (negative anchor): oe > 0.8  OR  p_bh > 0.05
```

## Cross-Assembly Mapping Method (hg19 → hg38)

Since RMC amino acid positions (`start_aa_raw`, `end_aa_raw`) are in **protein space**,
genomic liftover is inappropriate. Instead, protein-coordinate transfer is used:

1. For each RMC region, check if its ENST transcript has an identical protein
   sequence in both GRCh37 and GRCh38 (using `transcripts_identical_between_GRCh38_and_GRCh37.tsv`,
   verified by `aa_md5` match)
2. If identical → transfer amino acid positions directly (protein space is assembly-invariant)
3. If not identical → **drop** the region
4. Genomic columns (`contig`, `genomic_start`, `genomic_end`) are dropped since they
   are hg19-specific and cannot be directly transferred

**Script**: `test_acmg_auto/generate_rmc_hg38.py`

## Limitations

- **1,431 regions (8.0%) dropped** during hg19 → hg38 transfer, affecting 382 transcripts
  in 378 genes. These are transcripts whose protein sequence differs between GRCh37 and
  GRCh38 annotation (updated gene models, alternative canonical isoforms, or coding
  sequence corrections).
- **No genomic coordinates in hg38 file**: The `contig`, `genomic_start`, `genomic_end`
  columns are hg19-specific. For analyses requiring hg38 genomic coordinates, a separate
  CDS→genome mapping step would be needed using the hg38 transcript annotation.
- **gnomAD v4.x does NOT provide sub-genic RMC**: The gnomAD v4.1 constraint release
  only contains gene-level oe metrics. There is no region-level RMC for GRCh38 natively.
  The derived hg38 file here is the only available sub-genic RMC data for hg38.
- **Potential for stale RMC boundaries**: RMC region boundaries were determined using
  gnomAD v2 variant data (mostly exome data from ~125k individuals). The statistical
  power to detect constraint varies by region length and variant density. Short regions
  with few expected variants may have unreliable oe estimates (check `obs` and `exp`
  columns).
