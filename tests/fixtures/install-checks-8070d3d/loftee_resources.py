import gzip
from pathlib import Path
import sqlite3
import sys
try:
    ancestor, conservation, assembly = sys.argv[1:]
    with gzip.open(ancestor, "rt") as stream:
        if not stream.readline().startswith(">"):
            raise ValueError("LoFtee ancestor is not FASTA")
    lengths = {}
    with open(ancestor + ".fai", encoding="utf-8") as stream:
        for line in stream:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 5 or not all(value.isdigit() for value in fields[1:5]) or int(fields[1]) < 1:
                raise ValueError("Malformed LoFtee ancestor FASTA index")
            name = fields[0].removeprefix("chr")
            if name in lengths and lengths[name] != int(fields[1]):
                raise ValueError("Conflicting LoFtee ancestor chromosome lengths")
            lengths[name] = int(fields[1])
    expected = 249250621 if assembly in ("hg19", "GRCh37") else 248956422
    if lengths.get("1") != expected:
        raise ValueError("LoFtee ancestor FASTA index has the wrong assembly")
    with sqlite3.connect(Path(conservation).resolve().as_uri() + "?mode=ro", uri=True) as database:
        if database.execute("PRAGMA quick_check").fetchall() != [("ok",)]:
            raise ValueError("LoFtee conservation SQLite integrity check failed")
        if not database.execute("SELECT name FROM sqlite_master WHERE type='table' LIMIT 1").fetchone():
            raise ValueError("LoFtee conservation SQLite has no tables")
except (OSError, ValueError, sqlite3.Error) as error:
    print(f"ERROR: Invalid LoFtee resource: {error}", file=sys.stderr)
    sys.exit(1)
