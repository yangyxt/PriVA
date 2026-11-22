# gnomAD regional missense constraint

Bundled [gnomAD v2.1.1](https://gnomad.broadinstitute.org/) regional missense
constraint snapshots. `gnomad_v2.1.1_rmc_flat.hg19.tsv` is the GRCh37 source;
`gnomad_v2.1.1_rmc_flat.hg38.tsv` transfers protein intervals only where the
Ensembl transcript protein sequence is identical between assemblies.

| Install step | Consumer | Purpose |
| --- | --- | --- |
| C11 (`modules/local/install_steps/pm1_regions.nf`) | `priva pm1-regions` (from `scripts/prepare_pm1_regions.py`) | Regional constraint contribution to PM1 |

The tables contain transcript identifiers, amino-acid bounds, observed/expected
counts, `oe` and statistical fields. GRCh37 genomic coordinates are not reused
as GRCh38 coordinates; the transferred table omits those genomic columns.
The install step chooses the assembly-matched snapshot. These are calibration
resources, not ready-made PM1 assignments. See the [user manual](../../README.md).
