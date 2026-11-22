import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import sys
import tempfile
import yaml

config, manifest, output, receipt = map(Path, sys.argv[1:5])
config_sha, manifest_sha = sys.argv[5:7]

def identity(path):
    value = path.stat()
    return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)

def sha256(path):
    print('SpliceAI2: verify SHA256 ' + str(path), flush=True)
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

def no_links_or_partial(path):
    if not path.is_absolute() or Path(os.path.abspath(path)) != path:
        raise SystemExit('SpliceAI2 receipt paths must be absolute and canonical.')
    for part in (path, *path.parents):
        if '.partial' in part.name.lower() or part.is_symlink():
            raise SystemExit('SpliceAI2 refuses symlinks or .partial paths: ' + str(path))

no_links_or_partial(output)
no_links_or_partial(receipt)
if receipt.parent != output or not receipt.is_file() or receipt.stat().st_size == 0:
    raise SystemExit('SpliceAI2 requires a nonempty regular adapter receipt under OUTPUT_DIR.')
lock_path = output / '.install.lock'
no_links_or_partial(lock_path)
lock = lock_path.open('a')
try:
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
except BlockingIOError:
    raise SystemExit('SpliceAI2 preparation is active; refusing receipt reuse.')
if any(path.name.startswith('.building-') for path in output.iterdir()):
    raise SystemExit('SpliceAI2 output is incomplete or still being prepared; refusing reuse.')
tracked = {receipt: identity(receipt), manifest: identity(manifest), config: identity(config)}
data = json.loads(receipt.read_bytes())
if not isinstance(data, dict):
    raise SystemExit('SpliceAI2 adapter receipt must be an object.')
expected = dict(schema_version=1, resource='SpliceAI2', version='2.0', assembly='GRCh38',
                status='prepared_annotation_only', scientific_qualified=False,
                coverage='complete_published',
                source_revision='cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9',
                source_manifest_sha256=manifest_sha)
for key, value in expected.items():
    if data.get(key) != value or type(data.get(key)) is not type(value):
        raise SystemExit('SpliceAI2 receipt does not qualify full cache preparation: ' + key)
if sha256(manifest) != manifest_sha:
    raise SystemExit('SpliceAI2 source manifest changed during preparation.')
original = config.read_bytes()
if hashlib.sha256(original).hexdigest() != config_sha:
    raise SystemExit('PriVA config changed during SpliceAI2 preparation; refusing to overwrite it.')
outputs = data.get('outputs')
if not isinstance(outputs, dict) or set(outputs) != {'snv', 'indel'}:
    raise SystemExit('SpliceAI2 receipt requires both SNV and indel outputs.')
paths = {}
seen = set()
for role in ('snv', 'indel'):
    item = outputs[role]
    if not isinstance(item, dict):
        raise SystemExit('SpliceAI2 output entry must be an object: ' + role)
    for key, checksum_key in (('path', 'sha256'), ('index_path', 'index_sha256')):
        value, checksum = item.get(key), item.get(checksum_key)
        if not isinstance(value, str) or not isinstance(checksum, str) or not re.fullmatch('[0-9a-f]{64}', checksum):
            raise SystemExit('SpliceAI2 output path/checksum is missing or invalid: ' + role)
        path = Path(value)
        no_links_or_partial(path)
        if output not in path.parents or not path.is_file() or path.stat().st_size == 0 or path in seen:
            raise SystemExit('SpliceAI2 output must be a distinct nonempty regular file under OUTPUT_DIR: ' + str(path))
        seen.add(path)
        tracked[path] = identity(path)
        if key == 'path':
            if type(item.get('bytes')) is not int or item['bytes'] != path.stat().st_size:
                raise SystemExit('SpliceAI2 output size mismatch: ' + role)
            if not path.name.endswith(('.vcf.gz', '.vcf.bgz')):
                raise SystemExit('SpliceAI2 output must have a BGZF VCF filename: ' + role)
            paths[role] = str(path)
        elif value != item['path'] + '.tbi':
            raise SystemExit('SpliceAI2 tabix index must match its VCF path: ' + role)
        if sha256(path) != checksum:
            raise SystemExit('SpliceAI2 output SHA256 mismatch: ' + str(path))

updates = {
    'spliceai_snv_prescore': paths['snv'],
    'spliceai_indel_prescore': paths['indel'],
}
settings = yaml.safe_load(original)
text = original.decode()
for key, value in updates.items():
    pattern = re.compile(r'^' + re.escape(key) + r':[^\r\n]*', re.MULTILINE)
    if len(pattern.findall(text)) > 1:
        raise SystemExit('Duplicate SpliceAI2 config key: ' + key)
    replacement = key + ': ' + json.dumps(value)
    if pattern.search(text):
        text = pattern.sub(lambda match: replacement, text)
    else:
        text = text + ('' if not text or text.endswith('\n') else '\n') + replacement + '\n'
modified = yaml.safe_load(text)
if modified != dict(settings, **updates):
    raise SystemExit('SpliceAI2 config update would alter unrelated settings; refusing it.')
for path, before in tracked.items():
    if identity(path) != before:
        raise SystemExit('SpliceAI2 verification input changed: ' + str(path))
# Preserve all existing keys/comments except the original two cache paths.
# Preparation provenance and scientific status remain in the sidecar receipt.
name = None
try:
    with tempfile.NamedTemporaryFile(mode='w', dir=config.parent, prefix=config.name + '.', delete=False) as stream:
        name = stream.name
        stream.write(text)
        stream.flush()
        os.fsync(stream.fileno())
    os.chmod(name, stat.S_IMODE(config.stat().st_mode))
    os.replace(name, config)
finally:
    if name and os.path.exists(name):
        os.unlink(name)
print('SpliceAI2 preparation recorded; scientific activation remains pending.', flush=True)
