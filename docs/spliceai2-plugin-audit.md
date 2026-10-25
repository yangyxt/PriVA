# Official SpliceAI plugin boundary test

## Current release 116 qualification — 2026-10-16

The adapter accepts two explicit plugin SHA256 identities, retains the 113 route,
and records the selected plugin's actual SHA256, release and matching behavior
in each receipt. The current official [release-116 plugin](https://github.com/Ensembl/VEP_plugins/blob/release/116/SpliceAI.pm)
is pinned to `d08ce7d8c2e6b3229cef638beff145e465e224f0e3fcd43fd49d501ebc2740b1`.

```text
Publisher chr1:957127 GG>AA, model strand −
    → verified NOC2L context (ENSG00000188976.11, strand −)
    → unchanged GG>AA VCF record
    → official 116: exact position + REF + ALT + gene symbol
```

| Boundary | Historical 113 | Current 116 |
| --- | --- | --- |
| Allele matching | Ensembl normalization matcher | Exact position/REF/ALT; VF indels containing `-` first convert to VCF |
| Padded versus minimal substitution | Matcher can normalize | Different representations do not match |
| Negative model/transcript strand | Does not flip source alleles | Does not flip source alleles |
| Negative VariationFeature (VF) strand | Variant orientation | Variant orientation; reverse-complement before exact comparison |

`GG>AA` is already minimal and remains one intact MNV. The new
`tests/spliceai2/test_official_plugin_116.pl` pins the current plugin and actual
Ensembl 116 API bytes, constructs real VF/transcript/allele objects, and replaces
only plugin file registration and candidate retrieval. Its synthetic run passed
**4 subtests / 24 assertions** in Slurm step **1514363.355**, covering all eight
lexical values, allele/gene/locus negatives, allele number, VF orientation and
the changed representation matching. The optional `--projected-vcf` argument
reads a bounded local artifact for the exact NOC2L MNV; licensed score strings
are not embedded in the repository. Separate source-to-projection validation
checks publisher values and provenance. This boundary test does not qualify
real tabix retrieval, full VEP execution or scientific interpretation.

## Historical release 113 evidence

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
