import hashlib
import json
import os
from pathlib import Path
import sys
import yaml

config, manifest, output = map(Path, sys.argv[1:])
content = config.read_bytes()
settings = yaml.safe_load(content)
if not isinstance(settings, dict) or settings.get('assembly') not in ('hg38', 'GRCh38'):
    raise SystemExit('SpliceAI2 preparation requires hg38/GRCh38; no fallback is performed.')
if not manifest.is_absolute() or not output.is_absolute():
    raise SystemExit('SpliceAI2 manifest and output directory must be absolute paths.')

def partial(path):
    return any('.partial' in part.lower() for part in Path(path).parts)

if (partial(manifest) or partial(output) or manifest.resolve() != manifest
        or output.resolve() != output or Path(os.path.abspath(output)) != output or not manifest.is_file()):
    raise SystemExit('SpliceAI2 requires a complete regular manifest and a non-partial output path.')
raw = manifest.read_bytes()
data = json.loads(raw)
if (not isinstance(data, dict) or data.get('status') != 'downloads_complete'
        or not isinstance(data.get('confirmed_by'), str) or not data['confirmed_by'].strip()
        or data.get('assembly') != 'GRCh38' or data.get('version') != '2.0'
        or data.get('source_revision') != 'cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9'):
    raise SystemExit('SpliceAI2 requires an owner-confirmed complete GRCh38 v2.0 source manifest at the pinned revision.')

def inspect_paths(value):
    if isinstance(value, dict):
        for key, item in value.items():
            if key == 'path' and (not isinstance(item, str) or partial(item)):
                raise SystemExit('SpliceAI2 does not consume .partial source paths.')
            inspect_paths(item)
    elif isinstance(value, list):
        for item in value:
            inspect_paths(item)

inspect_paths(data)
if output.exists():
    if not output.is_dir() or any(path.name.startswith('.building-') for path in output.iterdir()):
        raise SystemExit('SpliceAI2 output is incomplete or still being prepared; refusing reuse.')
    if not (output / 'adapter.ready.json').exists() and any(output.iterdir()):
        raise SystemExit('SpliceAI2 output has no completed receipt; refusing duplicate conversion.')
print(hashlib.sha256(content).hexdigest(), hashlib.sha256(raw).hexdigest())
