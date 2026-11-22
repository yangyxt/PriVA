# ClinGen resources

`clingen_map.hg19.pkl.gz` and `clingen_map.hg38.pkl.gz` are assembly-specific Python mapping
snapshots. The corresponding `.v1.rkyv` files use `PRIVARKV` schema 10, version 1 and are read
with the checked `cache_convert::clingen_map::ClingenCache` reader. Evidence codes remain
text; a present float NaN code is represented by `None`, separately from an absent variant.

`scripts/install_utils.sh` provides `ClinGen_deploy` for the `clingen_map` configuration
resource, selecting the matching assembly. `scripts/prepare_clingen_map.py` produces the
Python mappings; `priva cache-convert clingen_map --pickle SOURCE --output DESTINATION`
converts either assembly to the Rust format. The original pickles and the gene dosage
sensitivity tables remain available.

| Rust resource | Bytes | SHA256 |
| --- | ---: | --- |
| `clingen_map.hg19.v1.rkyv` | 632556 | `8adb180d3ac36c2ea1454dd45d7c7bd98dfcb519f58514b64184d63386fad8eb` |
| `clingen_map.hg38.v1.rkyv` | 632632 | `2dfd8104a1c899e7c404215bce77d5969ebf37eb63fddafd93e2277ea13ff765` |
