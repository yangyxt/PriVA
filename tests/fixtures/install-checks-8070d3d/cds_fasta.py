import gzip
from pathlib import Path
import re
import sys

path = Path(sys.argv[1])
assembly = {"hg19": "GRCh37", "GRCh37": "GRCh37", "hg38": "GRCh38", "GRCh38": "GRCh38"}[sys.argv[2]]
opener = gzip.open if str(path).endswith(".gz") else open
records = 0
has_sequence = False
with opener(path, "rt", encoding="utf-8") as handle:
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if records and not has_sequence:
                raise ValueError(f"CDS FASTA contains an empty sequence: {path}")
            if not line[1:].strip():
                raise ValueError(f"CDS FASTA contains an empty header: {path}")
            declared = set(re.findall(r"\bGRCh(?:37|38)\b", line))
            if declared and declared != {assembly}:
                raise ValueError(f"CDS FASTA header conflicts with {assembly}: {path}")
            records += 1
            has_sequence = False
        else:
            if not records:
                raise ValueError(f"Not a CDS FASTA: {path}")
            has_sequence = True
if not records or not has_sequence:
    raise ValueError(f"CDS FASTA has no complete sequence: {path}")
print(f"Validated CDS FASTA structure: {path}; records={records}")
