# PriVA's copies of two CADD-scripts files

PriVA runs upstream CADD-scripts v1.7.2 (`CADD.sh` and its Snakemake workflow)
with all of CADD's own data. `main_install` (`CADD_install`, through
`scripts/cadd_scripts_install.py`) replaces two of its scripts with the copies
here, which write the same output faster. The copies are derived from
CADD-scripts and remain under CADD's license (`CADD-LICENSE`), not PriVA's.

| File | Replaces in `cadd_base_dir` | Change |
| --- | --- | --- |
| `extract_scored.py` | `src/scripts/extract_scored.py` | Prescore lookup: starts each lookup from the position index (`TABLE.posidx`, written by `scripts/cadd_prescore_index.py` during `CADD_install`) a few BGZF blocks before the variant, and reads nearby variants' rows in one pass. Same rows, same order. Runs under CADD's Python 2.7. When `TABLE.priva_compact` and its index match the table and `PRIVA_CADD_PRESCORE_BIN` names the `priva` program (set by `Calculate_CADD`), it hands the lookup to `priva cadd-prescore` instead (see below). |
| `MMSplice.py` | `src/scripts/lib/tools/MMSplice.py` | Saves the exon table read from CADD's Ensembl 110 GTF in the folder named by `PRIVA_CADD_MMSPLICE_EXON_CACHE` (set by PriVA's CADD calls) and reloads it on later runs. Same output. |
| `CADD-LICENSE` | — | CADD-scripts license and disclaimer |

Only files holding exactly the reviewed upstream bytes (identical in v1.7.2 and
v1.7.3) or an earlier PriVA copy are replaced. When a copy here changes, add
the Git blob ID of the version it replaces to `PREVIOUS` in
`scripts/cadd_scripts_install.py`.

## Compact copies for `priva cadd-prescore` (Nextflow branch, 2026-11-14)

`CADD_install` also runs `priva cadd-prescore-compact build TABLE` for each
prescored table. It writes `TABLE.priva_compact`, the six columns PriVA reads
from CADD's output (`#Chrom`, `Pos`, `Ref`, `Alt`, `FeatureID`, `PHRED`) in the
table's row order, in BGZF blocks of about 4 KB of text, and
`TABLE.priva_compact.idx`, the first position and file offset of every block.
Neither name ends in `.tsv.gz`, so CADD's prescore step does not take them for
prescored tables. `priva cadd-prescore` decompresses only the blocks holding a
wanted position, on several threads (`PRIVA_CADD_PRESCORE_THREADS`, default
the available CPUs up to 8).

Its output is the Python lookup's, with one planned difference: a found row
keeps all of the table's columns, but only the six kept columns carry values;
the others (including `RawScore`) are empty. CADD's `join` rule sorts these rows
together with the rows CADD scores itself, so the joined table keeps its
column layout, and PriVA reads only the six columns.
