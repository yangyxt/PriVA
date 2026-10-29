#!/usr/bin/env python3
"""Install PriVA's faster copies of two CADD-scripts files (scripts/cadd_scripts).

cadd_scripts_install.py --source CADD_DIR [--state]

extract_scored.py (the prescore lookup) and MMSplice.py write the same output
as upstream CADD's own files, faster. A file holding the reviewed upstream
bytes (identical in CADD-scripts v1.7.2 and v1.7.3) or an earlier PriVA copy
is replaced atomically, keeping its mode. Any other content is reported and
left unchanged, so CADD keeps running its own code there.
"""
import argparse
import hashlib
import os
from pathlib import Path
import sys
import tempfile


COPIES = Path(__file__).resolve().parent / "cadd_scripts"
# CADD path: (PriVA copy, Git blob ID of the reviewed upstream file)
FILES = {
    "src/scripts/extract_scored.py": ("extract_scored.py", "e88f103891da518a150d704f50d6bd3752970613"),
    "src/scripts/lib/tools/MMSplice.py": ("MMSplice.py", "d9cf8353b7a23874101bbee8f978e18210d705f0"),
}
# Git blob IDs of earlier PriVA copies. When a copy changes, add the ID of the
# version it replaces, so installed copies of that version are updated too.
PREVIOUS = set()


def blob_id(payload):
    return hashlib.sha1(b"blob %d\0" % len(payload) + payload).hexdigest()


def state(source, relative):
    path = Path(source) / relative
    if path.is_symlink() or not path.is_file():
        return "missing"
    payload = path.read_bytes()
    copy, upstream = FILES[relative]
    if payload == (COPIES / copy).read_bytes():
        return "priva"
    if blob_id(payload) == upstream:
        return "upstream"
    return "previous priva" if blob_id(payload) in PREVIOUS else "other"


def replace(target, payload, expected):
    """Write payload over target atomically if target still holds expected."""
    handle, temporary = tempfile.mkstemp(prefix="." + target.name + ".", dir=target.parent)
    try:
        with os.fdopen(handle, "wb") as stream:
            stream.write(payload)
        os.chmod(temporary, target.stat().st_mode & 0o777)
        if target.read_bytes() != expected:
            raise RuntimeError(f"{target} changed during replacement; retained")
        os.replace(temporary, target)
    finally:
        if os.path.lexists(temporary):
            os.unlink(temporary)


def install(source):
    for relative, (copy, _) in FILES.items():
        target = Path(source) / relative
        current = state(source, relative)
        if current in ("upstream", "previous priva"):
            replace(target, (COPIES / copy).read_bytes(), target.read_bytes())
            print(f"Installed PriVA's {relative} over the {current} file")
        elif current == "priva":
            print(f"PriVA's {relative} is already installed")
        else:
            print(f"Left {target} unchanged ({current}: not a reviewed upstream or PriVA file); "
                  "CADD keeps its own code there")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", type=Path, required=True, help="CADD-scripts directory (cadd_base_dir)")
    parser.add_argument("--state", action="store_true", help="Only report what each file holds")
    args = parser.parse_args()
    try:
        if not (args.source / "CADD.sh").is_file():
            raise RuntimeError(f"CADD.sh was not found in {args.source}")
        if args.state:
            for relative in FILES:
                print(relative, state(args.source, relative), sep="\t")
        else:
            install(args.source)
    except (OSError, RuntimeError) as error:
        sys.exit("ERROR: " + str(error))


if __name__ == "__main__":
    main()
