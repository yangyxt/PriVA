# InterPro resources

`Interpro_entry_mapping.pkl.gz` is the Python mapping snapshot. Its Rust counterpart,
`Interpro_entry_mapping.v1.rkyv`, uses `PRIVARKV` schema 9, version 1. Both retain the complete
mapping from `(database, accession)` to an ordered list of InterPro entries, including GO
three-tuples and the full abstract. The Rust archive pools repeated texts and is read with
the checked `cache_convert::interpro_mapping::InterproCache` reader.

`scripts/install_utils.sh` provides `InterPro_parsing` for the `interpro_mapping_pickle`
configuration resource. `scripts/protein_domain_mapping.py` produces the Python mapping;
`priva cache-convert interpro_mapping_pickle --pickle SOURCE --output DESTINATION` converts
that mapping to the Rust format. The original pickle remains available for Python users.

| Rust resource | Bytes | SHA256 |
| --- | ---: | --- |
| `Interpro_entry_mapping.v1.rkyv` | 39787928 | `29eacbb8f832b4ede0fa44e46a05d83120c02966ca1e081b2db3a9fd577bd4d1` |
