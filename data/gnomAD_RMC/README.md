# gnomAD regional missense constraint

Bundled [gnomAD v2.1.1](https://gnomad.broadinstitute.org/) regional missense
constraint snapshots. `gnomad_v2.1.1_rmc_flat.hg19.tsv` is the GRCh37 source;
`gnomad_v2.1.1_rmc_flat.hg38.tsv` transfers protein intervals only where the
Ensembl transcript protein sequence is identical between assemblies.

| Installer function in `scripts/install_utils.sh` | Consumer | Purpose |
| --- | --- | --- |
| `AlphaMissense_pick_missense_constrained_segments` | `scripts/am_pick_missense_constrained_segments.py` | RMC calibration anchors |
| `prepare_pm1_regions` | `scripts/prepare_pm1_regions.py` | Regional constraint contribution to PM1 |

The tables contain transcript identifiers, amino-acid bounds, observed/expected
counts, `oe` and statistical fields. GRCh37 genomic coordinates are not reused
as GRCh38 coordinates; the transferred table omits those genomic columns.
Installer functions choose the assembly-matched snapshot. These are calibration
resources, not ready-made PM1 assignments. See the [user manual](../../README.md).
