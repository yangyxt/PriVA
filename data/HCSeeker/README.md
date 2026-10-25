# HCSeeker Hotspot/Coldspot Data

## Source

HCSeeker identifies statistically significant somatic mutation hotspots and coldspots
in protein-coding genes using ClinVar pathogenic/benign variant distributions.

- **Publication**: [HCSeeker paper](https://doi.org/10.1002/humu.24647)
- **Original data**: Supplementary Table 2 from the HCSeeker publication
- **Original assembly**: GRCh38 (hg38)

## Files

| File | Assembly | Origin | Records |
|------|----------|--------|---------|
| `HC_spots.hg38.tsv` | GRCh38 | Original (from publication) | 1,669 (hotspots + coldspots) |
| `HC_spots.hg19.tsv` | GRCh37 | Derived via protein-coordinate transfer | 1,549 |

## Columns

| Column | Description |
|--------|-------------|
| `type` | `hotspot` or `coldspot` |
| `chr` | Chromosome number (bare integer, not `chr` prefix) |
| `gene` | Gene symbol |
| `HGNC_ID` | HGNC ID (e.g., `HGNC:7882`) |
| `iso` | RefSeq transcript accession with version (e.g., `NM_024408.4`) |
| `aa_start_pos` | Start amino acid position (1-based) |
| `aa_end_pos` | End amino acid position (1-based) |
| `PLP_number` | Number of pathogenic/likely pathogenic ClinVar variants in the region |
| `BLB_number` | Number of benign/likely benign ClinVar variants in the region |
| `profile-coefficient` | HCSeeker profile coefficient |
| `odd_benign` | Odds score for benign enrichment |
| `odd_path` | Odds score for pathogenic enrichment |

## Cross-Assembly Mapping Method (hg38 → hg19)

Since HCSeeker coordinates are in **protein space** (amino acid positions on NM_ transcripts),
genomic liftover is inappropriate. Instead, protein-coordinate transfer is used:

1. Map each record's NM_ transcript to ENST via MANE v1.5
2. Check if the ENST has an identical protein sequence in both GRCh38 and GRCh37
   (using `transcripts_identical_between_GRCh38_and_GRCh37.tsv`, verified by `aa_md5` match)
3. If identical → transfer amino acid positions directly (protein space is assembly-invariant)
4. If not identical → **drop** the record (protein-space coordinates may not correspond)

**Script**: `test_acmg_auto/generate_hcseeker_hg19.py`

## Limitations

- 120 spots (7.2%) dropped during hg38 → hg19 transfer due to non-identical protein sequences
  between assemblies. These correspond to transcripts where the coding sequence differs
  between GRCh37 and GRCh38 (e.g., due to updated gene models, alternative splice sites,
  or reference sequence corrections).
- The `iso` (NM_ accession) in the hg19 file retains the hg38 RefSeq version number.
  The transcript version may differ in hg19 annotation databases, but the protein sequence
  is verified identical so the amino acid coordinates remain valid.
- HCSeeker was developed using ClinVar data; benchmark analyses using HCSeeker anchors
  for ClinVar-derived validation could introduce circularity. Use for bandwidth calibration
  only, not as independent validation ground truth.

BioMart coverage audit, 2026-10-16: the separately acquired GRCh37.p13 mart116
mapping covers 1,469 of the 1,549 hg19 rows (882/911 hotspots, 587/638 coldspots),
representing 754/805 bare NM accessions. The unchanged consumers skip 80 unmapped
rows (29 hotspot, 51 coldspot), affecting the available HCSeeker calibration/PM1
contribution; no replacement or transcript choice was invented. Mapping/table
multiplicity and row order are retained. This audit does not establish new
protein-sequence identity or equivalence to the unavailable original mapping.
See [the mapping provenance](../../README.md#biomart-mapping) for the exact source,
5,084,679-byte table SHA256 and scientific gates. The current hg38 TSV contains
1,670 data rows (a measured correction to the earlier narrative count); the new
GRCh37 mapping is not authorized as a GRCh38 fallback.
