# HCSeeker resources

Bundled hotspot/coldspot tables from [HCSeeker](https://doi.org/10.1002/humu.24647).
`HC_spots.hg38.tsv` contains GRCh38-source records;
`HC_spots.hg19.tsv` is the GRCh37 protein-coordinate transfer snapshot.
Amino-acid positions are one-based on the stated RefSeq transcript, not genomic coordinates.

| Install step | Consumer | Purpose |
| --- | --- | --- |
| C10 (`modules/local/install_steps/das_domains.nf`) | `priva das-domains` (from `scripts/am_pick_intolerant_domains.py`) | Coldspot-overlap filtering |
| C11 (`modules/local/install_steps/pm1_regions.nf`) | `priva pm1-regions` (from `scripts/prepare_pm1_regions.py`) | Hotspot contribution to the PM1 cache |

These steps select the assembly-matched bundled table. RefSeq-to-Ensembl
mapping is supplied by step A5; unmatched transcripts do
not acquire substitute coordinates. Cross-assembly transfer is restricted to
identical protein sequences. HCSeeker is ClinVar-derived and is not independent
ClinVar validation truth. See the [user manual](../../README.md).
