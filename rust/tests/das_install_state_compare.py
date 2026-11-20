#!/usr/bin/env python3
"""Compare DAS completion receipts against the unchanged Python helper.

check --producer SCRIPT --binary PRIVA --work-root NEW_DIRECTORY
oracle --producer SCRIPT --exception JSON -- ORIGINAL_ARGUMENTS

Fixtures use only generated metadata and small files. No cache is regenerated.
The comparison retains stdout, exact exception class/message, receipt bytes,
mode, input preservation and both cross-language check directions.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
import time

sys.dont_write_bytecode = True


def save(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')


def identity(path):
    path = Path(path)
    return {'path': str(path), 'bytes': path.stat().st_size,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def oracle(args):
    spec = importlib.util.spec_from_file_location('das_state_reference', args.producer)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    last = {}
    def observe(frame, event, arg):
        if event == 'exception' and frame.f_code.co_filename == str(args.producer):
            category, value, _ = arg
            if category is not SystemExit:
                last.update(category=category.__name__, message=str(value))
        return observe
    sys.argv = [str(args.producer), *args.arguments]
    sys.settrace(observe)
    try:
        result = module.main()
    finally:
        sys.settrace(None)
        if last:
            save(args.exception, last)
    raise SystemExit(result)


def check(args):
    root = args.work_root.absolute()
    root.mkdir()
    producer, binary = args.producer.absolute(), args.binary.absolute()
    source_pin, binary_pin = identity(producer), identity(binary)
    result = {'producer': source_pin, 'binary': binary_pin, 'cases': [], 'unqualified': []}
    index = 0

    def run(arguments, base, label, language):
        folder = base / label
        folder.mkdir()
        if language == 'python':
            command = [sys.executable, str(Path(__file__).absolute()), 'oracle',
                       '--producer', str(producer), '--exception', str(folder/'exception.json'), '--', *arguments]
        else:
            command = [str(binary), 'das-install-state', *arguments]
        started = time.monotonic()
        completed = subprocess.run(command, cwd=base, capture_output=True, timeout=30,
            env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1', POLARS_MAX_THREADS='1', OMP_NUM_THREADS='1'))
        (folder/'stdout.log').write_bytes(completed.stdout)
        (folder/'stderr.log').write_bytes(completed.stderr)
        observation = {'exit': completed.returncode,
                       'stdout_sha256': hashlib.sha256(completed.stdout).hexdigest(),
                       'seconds': time.monotonic()-started}
        if completed.returncode == 1:
            if language == 'python' and (folder/'exception.json').exists():
                observation.update(json.loads((folder/'exception.json').read_text()))
            else:
                text = completed.stderr.decode()
                match = re.search(r'(?:^|\n)(?:priva: )?([A-Za-z]+Error): (.*)', text, re.S)
                if match:
                    observation.update(category=match[1], message=match[2].rstrip('\n'))
                else:
                    observation.update(category=None, message=text)
        return observation, completed.stdout

    def pair(arguments, base, name, expected_exit=None):
        nonlocal index
        index += 1
        py, pstdout = run(arguments, base, f'{index:03d}_python', 'python')
        native, nstdout = run(arguments, base, f'{index:03d}_native', 'native')
        keys = ['exit', 'stdout_sha256'] + (['category', 'message'] if py['exit'] == 1 else [])
        exact = all(py.get(key) == native.get(key) for key in keys)
        if expected_exit is not None:
            exact &= py['exit'] == expected_exit
        result['cases'].append({'name': name, 'python': py, 'native': native, 'exact': exact})
        return py, native, pstdout, nstdout

    def fixture(name):
        base = root/name
        base.mkdir()
        for name, content in [('input.bin', b'input\n'), ('output.bin', b'output\n'),
                              ('source.py', b'# source\n'), ('bundle.bin', b'bundled\n')]:
            (base/name).write_bytes(content)
        (base/'empty').write_bytes(b'')
        (base/'directory').mkdir()
        (base/'link.bin').symlink_to('input.bin')
        (base/'broken.bin').symlink_to('absent')
        (base/'loop').symlink_to('loop')
        common = ['--profile', 'fixture-v1', '--input', 'input.bin', '--source', 'analysis=source.py',
                  '--bundled', 'mapping=bundle.bin', '--parameter', 'iterations=1000']
        return base, common

    capture_cases = [
        ('ordinary', []), ('duplicate_inputs', ['--input', 'input.bin', '--input', 'link.bin']),
        ('parameter_order', ['--parameter', 'z=2', '--parameter', 'a=1=2', '--parameter', 'z=3']),
        ('duplicate_sources', ['--source', 'analysis=source.py']),
        ('duplicate_bundles', ['--bundled', 'mapping=bundle.bin']),
        ('unicode', ['--profile', 'éλ', '--parameter', '名=值']),
        ('blank_profile', ['--profile', '']), ('linked_input', ['--input', 'link.bin']),
        ('lexical_dotdot', ['--input', 'directory/../input.bin']),
        ('redundant_slashes', ['--input', './directory//.././input.bin']),
        ('bad_input', ['--input', 'absent']), ('empty_input', ['--input', 'empty']),
        ('directory_input', ['--input', 'directory']), ('loop_input', ['--input', 'loop']),
        ('broken_input', ['--input', 'broken.bin']),
        ('bad_source', ['--source', 'other=absent']), ('empty_source', ['--source', 'other=empty']),
        ('directory_source', ['--source', 'other=directory']),
        ('bad_bundle', ['--bundled', 'other=absent']), ('empty_bundle', ['--bundled', 'other=empty']),
        ('directory_bundle', ['--bundled', 'other=directory']), ('loop_bundle', ['--bundled', 'other=loop']),
        ('bad_named_source', ['--source', 'missing']), ('empty_role', ['--source', '=source.py']),
        ('empty_source_value', ['--source', 'analysis=']), ('bad_named_bundle', ['--bundled', '=bundle.bin']),
        ('bad_parameter', ['--parameter', 'missing']), ('empty_parameter', ['--parameter', 'value=']),
        ('space_parameter', ['--parameter', ' = ']),
        ('input_failure_first', ['--input', 'absent', '--source', 'missing']),
        ('source_failure_first', ['--source', 'other=absent', '--bundled', 'missing']),
        ('bundle_failure_first', ['--bundled', 'other=absent', '--parameter', 'missing']),
    ]
    for name, extra in capture_cases:
        base, common = fixture(name)
        pair(['capture', *common, *extra], base, name)
    for name, chmod_target in [('unread_bundle', 'bundle.bin'), ('unread_source', 'source.py')]:
        base, common = fixture(name)
        (base/chmod_target).chmod(0)
        try:
            pair(['capture', *common], base, name)
        finally:
            (base/chmod_target).chmod(0o600)

    for name, change in [
        ('roundtrip', None), ('changed_input_mtime', 'input_mtime'),
        ('changed_input_size', 'input_size'), ('changed_source', 'source'),
        ('same_source_new_inode', 'source_copy'), ('same_source_new_path', 'source_path'),
        ('changed_bundle_bytes', 'bundle'), ('empty_bundle_after_record', 'bundle_empty'),
        ('changed_output', 'output'), ('changed_input_mode_only', 'input_mode'),
        ('reordered_parameters', 'parameters'), ('reordered_duplicate_inputs', 'inputs'),
    ]:
        base, common = fixture(name)
        if change == 'inputs': common += ['--input', 'link.bin']
        _, _, captured, _ = pair(['capture', *common], base, name+'_capture', 0)
        expected = captured.decode().strip()
        receipt = base/'receipt.json'
        arguments = ['record', *common, '--output', 'output.bin', '--receipt', 'receipt.json', '--expected-inputs', expected]
        py, _ = run(arguments, base, 'record_python', 'python')
        pybytes = receipt.read_bytes() if receipt.exists() else None
        pymode = stat.S_IMODE(receipt.stat().st_mode) if receipt.exists() else None
        cross_n, _ = run(['check', *common, '--output', 'output.bin', '--receipt', 'receipt.json'], base, 'check_python_by_native', 'native')
        native, _ = run(arguments, base, 'record_native', 'native')
        nbytes = receipt.read_bytes() if receipt.exists() else None
        nmode = stat.S_IMODE(receipt.stat().st_mode) if receipt.exists() else None
        cross_p, _ = run(['check', *common, '--output', 'output.bin', '--receipt', 'receipt.json'], base, 'check_native_by_python', 'python')
        result['cases'].append({'name': name+'_publication', 'python': py, 'native': native,
            'receipt_bytes_equal': pybytes == nbytes, 'receipt_modes': [pymode,nmode],
            'cross_checks': [cross_n,cross_p], 'exact': pybytes == nbytes and pymode == nmode == 0o600
                and all(x['exit'] == 0 for x in (py,native,cross_n,cross_p))})
        if change == 'input_mtime':
            s=(base/'input.bin').stat();os.utime(base/'input.bin', ns=(s.st_atime_ns,s.st_mtime_ns+1))
        elif change == 'input_size': (base/'input.bin').write_bytes(b'input different\n')
        elif change == 'source': (base/'source.py').write_bytes(b'# changed\n')
        elif change in ('source_copy','source_path'):
            (base/'new_source.py').write_bytes((base/'source.py').read_bytes())
            if change == 'source_copy': (base/'new_source.py').replace(base/'source.py')
            else: common[common.index('analysis=source.py')] = 'analysis=new_source.py'
        elif change == 'bundle': (base/'bundle.bin').write_bytes(b'different bundled bytes\n')
        elif change == 'bundle_empty': (base/'bundle.bin').write_bytes(b'')
        elif change == 'output': (base/'output.bin').write_bytes(b'different output\n')
        elif change == 'input_mode': (base/'input.bin').chmod(0o400)
        elif change == 'parameters': common += ['--parameter', 'before=1']
        elif change == 'inputs':
            at=common.index('input.bin');common[at]='link.bin';common[-1]='input.bin'
        before = receipt.read_bytes()
        pair(['check', *common, '--output', 'output.bin', '--receipt', 'receipt.json'], base, name+'_check')
        result['cases'].append({'name': name+'_preservation', 'exact': receipt.read_bytes() == before,
            'receipt_preserved': receipt.read_bytes() == before})

    for label, raw in [('empty',''),('not_json','bad'),('wrong','{}'),('null','null'),
                       ('truncated','{"a":'),('nan','NaN'),('array','[]')]:
        base, common = fixture('expected_'+label)
        receipt=base/'receipt.json';receipt.write_bytes(b'preserve\n')
        pair(['record',*common,'--output','output.bin','--receipt','receipt.json','--expected-inputs',raw],
             base,'expected_'+label,1)
        result['cases'].append({'name':'expected_'+label+'_preserve','exact':receipt.read_bytes()==b'preserve\n'})
    for label, stamp in [('exact_float', 1<<60), ('rounded_float', (1<<60)+3)]:
        base,common=fixture('timestamp_'+label)
        os.utime(base/'input.bin',ns=(stamp,stamp))
        _,_,captured,_=pair(['capture',*common],base,label+'_capture',0)
        expected=json.loads(captured);expected['inputs'][0]['mtime_ns']=float(stamp)
        (base/'receipt.json').write_bytes(b'preserve\n')
        pair(['record',*common,'--output','output.bin','--receipt','receipt.json',
              '--expected-inputs',json.dumps(expected)],base,label+'_record',0 if label=='exact_float' else 1)
        # A forged float in an otherwise valid receipt must also be rejected
        # whenever it lost nanoseconds, in either comparison direction.
        plain=json.loads(captured)
        outputs=json.loads(run(['capture','--profile','x','--input','output.bin'],base,'output_identity','python')[1])['inputs']
        receipt={'version':3,**plain,'outputs':outputs}
        receipt['inputs'][0]['mtime_ns']=float(stamp)
        (base/'receipt.json').write_text(json.dumps(receipt))
        pair(['check',*common,'--output','output.bin','--receipt','receipt.json'],base,label+'_check',
             0 if label=='exact_float' else 1)
    for label, content in [('empty',b''),('not_json',b'bad'),('wrong',b'{}'),('utf8',b'\xff'),('crlf',b'{\r\n}')]:
        base, common = fixture('receipt_'+label)
        (base/'receipt.json').write_bytes(content)
        pair(['check',*common,'--output','output.bin','--receipt','receipt.json'],base,'receipt_'+label,1)
    for label, options in [('missing_output',['--output','absent']), ('empty_output',['--output','empty']),
                           ('directory_output',['--output','directory'])]:
        base,common=fixture(label)
        _,_,capture,_=pair(['capture',*common],base,label+'_capture',0)
        (base/'receipt.json').write_bytes(b'preserve\n')
        pair(['record',*common,*options,'--receipt','receipt.json','--expected-inputs',capture.decode()],base,label,1)
        result['cases'].append({'name':label+'_preserve','exact':(base/'receipt.json').read_bytes()==b'preserve\n'})
    for label, tail in [('missing_required',[]),('missing_value',['--input']),('unknown',['--other','x']),
                        ('ambiguous',['--p','x']),('unexpected_output',['--output','output.bin'])]:
        base, common = fixture('cli_'+label)
        pair(['capture',*([] if label=='missing_required' else common),*tail],base,'cli_'+label,2)
    result['source_unchanged'] = identity(producer) == source_pin
    result['binary_unchanged'] = identity(binary) == binary_pin
    result['failed'] = [v['name'] for v in result['cases'] if not v['exact']]
    result['passed'] = not result['failed'] and result['source_unchanged'] and result['binary_unchanged']
    save(root/'comparison.json', result)
    print(json.dumps({'passed':result['passed'],'cases':len(result['cases']),'failed':result['failed']}))
    if not result['passed']: raise SystemExit(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    p=commands.add_parser('check')
    for name in ('producer','binary','work-root'):p.add_argument('--'+name,type=Path,required=True)
    p.set_defaults(function=check)
    p=commands.add_parser('oracle');p.add_argument('--producer',type=Path,required=True)
    p.add_argument('--exception',type=Path,required=True);p.add_argument('arguments',nargs=argparse.REMAINDER)
    p.set_defaults(function=oracle)
    args=parser.parse_args()
    if getattr(args,'arguments',None) and args.arguments[0]=='--':args.arguments=args.arguments[1:]
    args.function(args)


if __name__ == '__main__':main()
