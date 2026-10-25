# HCSeeker resources

Bundled hotspot/coldspot tables from [HCSeeker](https://doi.org/10.1002/humu.24647).
`HC_spots.hg38.tsv` contains GRCh38-source records;
`HC_spots.hg19.tsv` is the GRCh37 protein-coordinate transfer snapshot.
Amino-acid positions are one-based on the stated RefSeq transcript, not genomic coordinates.

| Installer function in `scripts/install_utils.sh` | Consumer | Purpose |
| --- | --- | --- |
| `AlphaMissense_pick_missense_constrained_segments` | `scripts/am_pick_missense_constrained_segments.py` | Hotspot/coldspot calibration anchors |
| `AlphaMissense_pick_intolerant_domains` | `scripts/am_pick_intolerant_domains.py` | Coldspot-overlap filtering |
| `prepare_pm1_regions` | `scripts/prepare_pm1_regions.py` | Hotspot contribution to the PM1 cache |

These functions select the assembly-matched bundled table. RefSeq-to-Ensembl
mapping is supplied by `nm_to_enst_mapping_install`; unmatched transcripts do
not acquire substitute coordinates. Cross-assembly transfer is restricted to
identical protein sequences. HCSeeker is ClinVar-derived and is not independent
ClinVar validation truth. See the [user manual](../../README.md).
