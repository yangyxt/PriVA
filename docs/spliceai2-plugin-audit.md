# Official SpliceAI plugin boundary test

`tests/spliceai2/test_official_plugin.pl` tests synthetic projected records through
the installed, unchanged release-113 SpliceAI plugin and the actual Ensembl allele
matcher. It does not contain licensed SpliceAI2 cache rows or reproduce inference.
The initial run passed **11 subtests / 38 assertions** with `priva_acmg` Perl
5.32.1 and VEP 113.2.

```text
Invented VCF records → real parse_data / run / Ensembl allele matching
                              ↑
                    fixture tabix candidates + light VEP objects
```

Run inside an authorized Slurm allocation, with the environment's `bin` on PATH:

```bash
perl tests/spliceai2/test_official_plugin.pl \
  --plugin-file /path/to/VEP_plugin/SpliceAI.pm \
  --vep-lib /path/to/priva_acmg/share/ensembl-vep-113.2-0
```

The harness rejects different plugin or matcher bytes. Audited identities:

| Implementation | Identity |
| --- | --- |
| SpliceAI.pm | VEP_plugins commit `10054085053745b62ab20129f13cd02f06f7e183`; SHA256 `500a534c403bcbcc15753718768f0a1f352cc6b93f9119e5eba25745f7205b9a` |
| Ensembl Variation Utils Sequence.pm | Installed VEP 113.2; SHA256 `a85aa34bb337798dd45599270340aa56740f2d53496be1590e610480ff8f6d2e` |

Cases cover exact eight-field projection, multiple alleles/genes, anchored and
minimal indels, padded substitutions, negative VariationFeature orientation,
exact gene symbols/HGNC fallback, absent records and the inclusive 0.25 cutoff
flag. Low scores and `-0.0` remain present. The real plugin's duplicate-gene
overwrite, comma-packed gene loss and extra-INFO contamination are demonstrated;
the producer must prevent those formats.

The plugin does **not** check transcript/model strand: the region/GTF adapter
must establish that association. Its VariationFeature strand controls allele
orientation and is a different concept. `add_file` and `get_data` are fixtures;
the latter supplies candidates without testing real tabix retrieval. Lightweight
objects expose transcript/VariationFeature APIs; the minimal-insertion branch
still calls the real Ensembl `VariationFeature::to_VCF_record` implementation.
No real reference sequence, chromosome-alias resolution, BGZF index, VEP cache,
actual converter output or full VEP execution is qualified by these checks.
SpliceAI2 scientific interpretation, CADD compatibility and family parity remain
separate gates.
