import json
from pathlib import Path
import re
import sys

path = Path(sys.argv[1])
try:
    records = json.loads(path.read_text())
    if isinstance(records, dict):
        records = records.get("packages")
    if not isinstance(records, list) or not records:
        raise ValueError("package list must be a nonempty JSON array")
    pins = {}
    for record in records:
        if not isinstance(record, dict):
            raise ValueError("installed package records must be JSON objects")
        name = record["name"]
        build = record.get("build_string") or record.get("build")
        if (record.get("channel") == "pypi" or record.get("platform") == "pypi"
                or build == "pypi_0" or name == "ensembl-vep"):
            continue
        version = record["version"]
        if not all(isinstance(value, str) and re.fullmatch(r"[A-Za-z0-9_.+!-]+", value)
                   for value in (name, version, build)):
            raise ValueError(f"cannot pin installed package {name}")
        if name in pins:
            raise ValueError(f"duplicate installed package {name}")
        pins[name] = f"{name}={version}={build}"
    if not pins:
        raise ValueError("no installed Conda dependencies could be pinned")
    path.write_text("\n".join(pins[name] for name in sorted(pins)) + "\n")
except (OSError, ValueError, KeyError, TypeError) as exc:
    print(f"ERROR: Cannot preserve installed Conda packages while changing VEP: {exc}", file=sys.stderr)
    sys.exit(1)
