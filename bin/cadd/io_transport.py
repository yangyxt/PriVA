#!/usr/bin/env python3
"""Preserve UTF-8 text around unchanged upstream CADD predictors.

regseq --bgzip PROGRAM --work-dir DIRECTORY -- SCRIPT ORIGINAL_ARGUMENTS...
mmsplice --input EXACT_INPUT -- SCRIPT ORIGINAL_ARGUMENTS...

RegSeq writes a fresh plain VCF, which mature bgzip compresses and this utility
verifies before exclusive publication. MMSplice changes only the implicit text
encoding of its exact input. Scientific arguments, packages and models remain
upstream. This module uses Python 3.6-compatible standard-library interfaces.
"""
import argparse
import codecs
import gzip
import hashlib
import os
from pathlib import Path
import runpy
import stat
import subprocess
import sys
import tempfile


REGSEQ_OPTIONS = ['--variants', '--model', '--weights', '--reference', '--genome', '--output']
MMSPLICE_OPTIONS = ['-i', '-g', '-f']
CHUNK = 1024 * 1024


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def identity(path):
    """Observe only an ordinary script/input, without following a final link."""
    value = os.lstat(str(path))
    require(stat.S_ISREG(value.st_mode), 'Expected an ordinary input or script file')
    return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)


def source_pin(path):
    before = identity(path)
    digest = hashlib.sha256()
    with open(str(path), 'rb') as stream:
        for block in iter(lambda: stream.read(CHUNK), b''):
            digest.update(block)
    require(identity(path) == before, 'Predictor source changed while reading')
    return before, digest.hexdigest()


def original_argv(argv, options):
    require(len(argv) == 1 + 2 * len(options) and argv[1::2] == options
            and all(isinstance(value, str) and value for value in argv),
            'Expected the unchanged original predictor option sequence')
    require(Path(argv[0]).is_absolute(), 'Predictor script must use its original absolute path')
    return list(argv)


def invoke_predictor(argv):
    """Match direct-script import/argv context, restoring both on every exit."""
    previous_argv, previous_path = sys.argv, sys.path
    require(isinstance(previous_path, list) and previous_path, 'Missing interpreter import context')
    try:
        sys.argv = list(argv)
        sys.path = [os.path.dirname(os.path.abspath(argv[0]))] + previous_path[1:]
        try:
            runpy.run_path(argv[0], run_name='__main__')
        except SystemExit as error:
            if error.code not in (None, 0):
                raise
    finally:
        sys.argv, sys.path = previous_argv, previous_path


def plain_digest(path):
    before = identity(path)
    digest, size = hashlib.sha256(), 0
    decoder = codecs.getincrementaldecoder('utf-8')('strict')
    with open(str(path), 'rb') as stream:
        for block in iter(lambda: stream.read(CHUNK), b''):
            decoder.decode(block)
            digest.update(block)
            size += len(block)
        decoder.decode(b'', final=True)
    require(identity(path) == before and size == before[2], 'Private plain output changed while verifying')
    require(size > 0, 'Private plain output is empty')
    return size, digest.hexdigest(), before


def compress(program, plain, partial, log):
    with open(str(partial), 'xb', buffering=0) as target, open(str(log), 'xb') as stderr:
        subprocess.run([program, '-c', str(plain)], stdout=target, stderr=stderr, check=True)
        os.fsync(target.fileno())


def compressed_digest(path):
    digest, size = hashlib.sha256(), 0
    # Reading to EOF also checks the CRC of every gzip/BGZF member.
    with gzip.open(str(path), 'rb') as stream:
        for block in iter(lambda: stream.read(CHUNK), b''):
            digest.update(block)
            size += len(block)
    return size, digest.hexdigest()


def publish(partial, output):
    # link(2) refuses existing regular, empty and dangling destination paths.
    os.link(str(partial), str(output), follow_symlinks=False)
    os.unlink(str(partial))


def regseq(argv, bgzip, work_dir):
    argv = original_argv(argv, REGSEQ_OPTIONS)
    require(os.environ.get('LC_ALL') == 'C' and os.environ.get('PYTHONUTF8') == '1'
            and getattr(sys.flags, 'utf8_mode', 0) == 1,
            'RegSeq requires command-local PYTHONUTF8=1 with the existing LC_ALL=C')
    source_before = source_pin(argv[0])
    input_path = os.path.abspath(argv[2])
    input_before = identity(input_path)
    output = os.path.abspath(argv[-1])
    root = os.path.abspath(work_dir)
    require(os.path.isdir(root) and not os.path.islink(root)
            and os.path.commonpath([root, output]) == root and output != root
            and output.endswith('.gz') and os.path.isdir(os.path.dirname(output)),
            'RegSeq output must be a .gz file within its existing task work directory')
    require(not os.path.lexists(output), 'Refusing an existing RegSeq output')
    work = Path(tempfile.mkdtemp(prefix='regseq-utf8-', dir=root))
    plain, partial = work / 'regseq.vcf', work / 'regseq.bgzf.partial'
    actual = list(argv)
    actual[-1] = str(plain)
    invoke_predictor(actual)
    require(identity(input_path) == input_before and source_pin(argv[0]) == source_before,
            'RegSeq input or predictor source changed during execution')
    expected = plain_digest(plain)
    compress(bgzip, plain, partial, work / 'compression.log')
    require(compressed_digest(partial) == expected[:2],
            'Compressed output differs from the complete UTF-8 plain output')
    require(identity(plain) == expected[2] and identity(input_path) == input_before
            and source_pin(argv[0]) == source_before, 'RegSeq input/source/plain output changed before publication')
    publish(partial, output)
    return expected[:2]


def mmsplice(argv, input_path):
    argv = original_argv(argv, MMSPLICE_OPTIONS)
    require(os.environ.get('LC_ALL') == 'C' and os.environ.get('PYTHONIOENCODING') == 'utf-8'
            and codecs.lookup(sys.stdout.encoding).name == 'utf-8' and sys.stdout.errors == 'strict',
            'MMSplice requires command-local PYTHONIOENCODING=utf-8 with the existing LC_ALL=C')
    exact_input = os.path.abspath(input_path)
    require(os.path.abspath(argv[2]) == exact_input, 'MMSplice input differs from its original -i operand')
    input_before, source_before = identity(exact_input), source_pin(argv[0])
    original_open = gzip.open

    def exact_open(filename, mode='rb', compresslevel=9, encoding=None, errors=None, newline=None):
        try:
            actual = os.path.abspath(os.fsdecode(filename))
        except TypeError:
            actual = None
        if actual == exact_input and mode == 'rt' and encoding is None:
            encoding = 'utf-8'
        return original_open(filename, mode, compresslevel=compresslevel,
                             encoding=encoding, errors=errors, newline=newline)

    try:
        gzip.open = exact_open
        invoke_predictor(argv)
    finally:
        gzip.open = original_open
    require(identity(exact_input) == input_before and source_pin(argv[0]) == source_before,
            'MMSplice input or predictor source changed during execution')


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    sub = result.add_subparsers(dest='procedure')
    command = sub.add_parser('regseq')
    command.add_argument('--bgzip', required=True)
    command.add_argument('--work-dir', required=True)
    command.add_argument('predictor_argv', nargs=argparse.REMAINDER)
    command = sub.add_parser('mmsplice')
    command.add_argument('--input', required=True)
    command.add_argument('predictor_argv', nargs=argparse.REMAINDER)
    return result


def main():
    command = parser()
    args = command.parse_args()
    if args.procedure is None:
        command.error('a procedure is required')
    argv = args.predictor_argv
    if argv[:1] == ['--']:
        argv = argv[1:]
    if args.procedure == 'regseq':
        regseq(argv, args.bgzip, args.work_dir)
    else:
        mmsplice(argv, args.input)


if __name__ == '__main__':
    main()
