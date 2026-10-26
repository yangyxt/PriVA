# VEP support files

`csq_interpretations.tsv` supplies consequence metadata to PriVA's annotation readers.

`UTRAnnotator116_all_effects.pm` retains Ensembl's Apache-2.0 license and the
original biological methods, while reporting all calculated UTR effects in a
deterministic order. `vep_plugins_release.py ensure-utr-reporting`, called by
`install_utils.sh vep_cache_api_install`, deploys it only over the matching
official VEP116 source. Other releases remain unchanged. The installer retains
the original module and an exact patch beside the installed plugins.
