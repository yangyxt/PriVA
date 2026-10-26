# HGNC master and transcript mappings

These are separate inputs:

```text
HGNC master snapshot       → gene symbols, aliases and identifiers
Assembly/release mapping   → NM accessions ↔ Ensembl transcripts → DAS / PM1
```

| Bundled file | Source and version |
| --- | --- |
| `hgnc_non_alt_loci_set_sha256_890b18743cc8.tsv` | HGNC master snapshot, identified by content SHA256. Publisher release, Ensembl release and assembly binding are unknown. |
| `non_alt_loci_set.tsv` | Relative compatibility link to the same HGNC master bytes. |
| `ensembl_GRCh37_mart116_nm_to_enst.tsv` | Ensembl GRCh37.p13, mart116 with annotation-base75 lineage. |
| `ensembl_GRCh38_mart113_nm_to_enst.tsv` | Ensembl GRCh38.p14, October 2024 archive / mart113. |
| `ensembl_GRCh38_mart116_ftp_candidate.tsv` | Ensembl GRCh38.p14 / mart116, deterministic projection from pinned official MySQL FTP tables. Scientific coordinate qualification remains separate. |

The Ensembl tables retain the original five-column layout: gene ID, transcript ID, HGNC symbol, UniProt ID and RefSeq mRNA ID. Row order, duplicate rows and blank RefSeq fields are retained. The FTP projection records its own deterministic order; it does not claim the BioMart service's response order or original export identity. Adjacent provenance files describe each selected snapshot.

[`nm_to_enst_mapping_install`](../../scripts/install_utils.sh) selects these filenames from assembly and VEP major when `nm_to_enst_tsv` is unset; an explicit mapping stays selected. VEP116 on GRCh38 selects FTP116; VEP113 and legacy calls without a release retain mart113. Bundled reuse checks existence and nonzero size only. External or symlinked Ensembl mappings retain source/assembly validation. Raw acquisition and projection work belongs outside the checkout, under `cache_root/BioMart` or a selected external workspace.

The HGNC master is read by the [gene mechanism resources](../../scripts/gene_mechanism_resources.py). Ensembl mappings are read by [domain selection](../../scripts/am_pick_intolerant_domains.py) and [PM1 preparation](../../scripts/prepare_pm1_regions.py). Changing an input does not qualify a new scientific profile or permit reuse of derived results with a different mapping identity.
