#!/usr/bin/env python3
"""Capture actual pinned coerce_pext_value before the independent Rust port.

capture --source-file PYTHON_SOURCE --output DATA_ORACLE
Compile the original AST verbatim with actual pandas/NumPy. Typed primitive
cases cross the CSV boundary; unsupported object/array diagnostics are retained
separately and cannot count as Rust parity. No production input is read.
"""
import argparse
import ast
import csv
import hashlib
import json
import os
import pickle
import platform
import random
import struct
import unicodedata
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


SOURCE_SHA = 'a9efd0964f25e84dc01ad276d4f0905a10492a5d1f94c2b277573ac94244d354'
FUNCTION_SHA = '67d26c68d7149f47719c3b8db2289895c0c4f77fc5227cc4ef288fe972a55dbe'
COMMIT = '09223ef5a8b0ceb4c6d61055c6eefe62b51ba903'


def bits(value):
    return struct.pack('>d', value).hex()


def double(raw):
    return struct.unpack('>d', bytes.fromhex(raw))[0]


def primitives(zeroes):
    for name, value in [('none', None), ('pd_na', pd.NA), ('pd_nat', pd.NaT),
                        ('numpy_nat', np.datetime64('NaT'))]:
        yield name, value, ['null', '']
    for value in (False, True):
        yield f'python_bool_{value}', value, ['bool', str(int(value))]
        yield f'numpy_bool_{value}', np.bool_(value), ['numpy_bool', str(int(value))]
    for value in (-2**63, -2**53 - 1, -1, 0, 1, 2**53 + 1, 2**63 - 1):
        yield f'python_int_{value}', value, ['int', str(value)]
    for value in (2**63, 2**64 - 1):
        yield f'python_uint_{value}', value, ['uint', str(value)]
    for dtype in (np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32, np.uint64):
        for value in (np.iinfo(dtype).min, np.iinfo(dtype).max):
            scalar = dtype(value)
            yield f'{dtype.__name__}_{value}', scalar, ['uint' if value > 2**63 - 1 else 'int', str(value)]
    for power in (65, 100, 511, 1023, 1024, 2000):
        for offset in (-1, 0, 1):
            for sign in (-1, 1):
                value = sign * (2**power + offset)
                yield f'bigint_{power}_{offset}_{sign}', value, ['bigint', str(value)]
    # Integers straddling f64 overflow rounding, including the exact midpoint.
    threshold = 2**1024 - 2**970
    for offset in (-2**918, -1, 0, 1, 2**918):
        for sign in (-1, 1):
            value = sign * (threshold + offset)
            yield f'bigint_overflow_midpoint_{offset}_{sign}', value, ['bigint', str(value)]
    for raw in ('0000000000000000', '8000000000000000', '0000000000000001', '8000000000000001',
                '0010000000000000', '7fefffffffffffff', 'ffefffffffffffff', '7ff0000000000000',
                'fff0000000000000', '7ff8000000000000', 'fff8000000000000', '7ff8000000000123',
                'fff8000000000456', '7ff0000000000001', 'fff0000000000001', '3fb999999999999a'):
        for family in ('python', 'numpy'):
            value = double(raw)
            if family == 'numpy': value = np.float64(value)
            yield f'{family}_float_{raw}', value, ['float', raw]
    for raw in ('00000000', '80000000', '00000001', '00800000', '7f7fffff', '7f800000',
                'ff800000', '7fc00123', 'ffc00456', '3dcccccd'):
        value = np.frombuffer(bytes.fromhex(raw), dtype='>f4')[0]
        yield f'numpy_float32_{raw}', value, ['float32', raw]
    strings = ['', '.', 'NA', 'NaN', 'nan', 'NAN', 'Nan', '+nan', '-NaN', '-nan', '+NaN',
               'inf', '-inf', '+INF', 'Infinity', '-INFINITY', 'infinite', 'nan(1)', 'NaN(payload)',
               'None', 'NULL', '<NA>', 'True', 'False', 'N/A', '1.#INF', '0x1.8p2', '0X10',
               '0', '-0', '+0', '-0.0', '+0.0', '.5', '1.', '+.5', '-.0', '1e3', '1e-3',
               '1_000', '1_2.3_4e+0_1', '_1', '1_', '1__2', '1_.0', '1._0', '1_e2', '1e_2',
               '1e+_2', '1e2_', '1e', '.', '+', '--1', '++1', '1 2', '1\x002', '1\x00',
               '−1', '²', 'Ⅷ', '½', '①', '１.２３', '١.٢_٣', '۱۲۳e-۲', '𝟙.𝟚e𝟛',
               '1,2', '1;2', '1|2', '1&&2', 'NA&.&bad', '1.2&-2,3e-1',
               '1e10000', '-1e10000', '1e-10000', '-1e-10000',
               '1.7976931348623157e308', '1.7976931348623158e308', '1.7976931348623159e308',
               '4.9406564584124654e-324', '2.4703282292062327e-324', '2.4703282292062328e-324']
    whitespace = [' ', '\t', '\n', '\r', '\v', '\f', '\x1c', '\x1d', '\x1e', '\x1f',
                  '\x85', '\xa0', '\u1680', '\u2000', '\u2007', '\u200a', '\u2028', '\u2029',
                  '\u202f', '\u205f', '\u3000', '\u200b', '\ufeff', '\x00']
    for space in whitespace:
        strings += [space + '.125' + space, space + '-NAN' + space, '1' + space + '2', space + 'NA' + space]
    atoms = ['0', '-0', '0.0', '-0.0', '.5', '1', '-1', 'nan', 'NaN', 'NAN', '+nan', '-nan', 'inf', '-inf', '.', 'bad', '']
    strings += [a + separator + b for a in atoms for b in atoms for separator in ('&', ',')]
    strings += ['NAN&0&1', '0&NAN&1', '-nan&inf&1', '1&-nan&inf', '-0&0&-0', '0&-0&0',
                'nan&NAN&0', 'NaN&-nan&0', 'NA&NaN&nan&NAN', '.&NA&bad']
    for zero in zeroes:
        strings += [chr(zero + digit) for digit in range(10)]
        strings += [chr(zero + 1) + '.' + chr(zero + 2) + '_' + chr(zero + 3),
                    chr(zero + 1) + 'e-' + chr(zero + 2)]
    rng = random.Random(2849)
    for _ in range(1200):
        value = double(f'{rng.getrandbits(64):016x}')
        if np.isfinite(value):
            strings += [repr(value), format(value, '.17g'), format(value, '.24g')]
    for _ in range(500):
        digits = ''.join(str(rng.randrange(10)) for _ in range(rng.randrange(1, 80)))
        point = rng.randrange(len(digits) + 1)
        strings.append(rng.choice(['', '+', '-']) + digits[:point] + '.' + digits[point:] +
                       rng.choice(['e', 'E']) + str(rng.randrange(-450, 451)))
    for i, text in enumerate(strings):
        yield f'text_{i}', text, ['str', text.encode().hex()]


def capture(args):
    raw = args.source_file.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == SOURCE_SHA
    text = raw.decode()
    function = next(node for node in ast.parse(text).body if isinstance(node, ast.FunctionDef) and node.name == 'coerce_pext_value')
    assert hashlib.sha256(ast.get_source_segment(text, function).encode()).hexdigest() == FUNCTION_SHA
    namespace = {'pd': pd, 'np': np}
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(args.source_file), 'exec'), namespace)
    original = namespace['coerce_pext_value']
    zeroes = [code for code in range(0x110000) if unicodedata.decimal(chr(code), None) == 0]
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / 'original_acmg_consequence.py').write_bytes(raw)
    (args.output / 'oracle_harness.py').write_bytes(Path(__file__).read_bytes())
    successful, errors = [], []
    for name, value, wire in primitives(zeroes):
        before = pickle.dumps(value, protocol=5)
        with warnings.catch_warnings(record=True) as warning:
            warnings.simplefilter('always')
            try:
                outcome = {'status': 'value', 'bits': bits(original(value))}
            except Exception as error:
                outcome = {'status': 'error', 'error_type': type(error).__name__, 'error_message': str(error)}
        assert pickle.dumps(value, protocol=5) == before, name
        entry = {'name': name, 'python_type': type(value).__module__ + '.' + type(value).__name__, 'wire': wire,
                 'input_unchanged': True, 'warnings': [type(item.message).__name__ + ': ' + str(item.message) for item in warning], **outcome}
        (successful if outcome['status'] == 'value' else errors).append(entry)
    with (args.output / 'values.input.csv').open('x', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['index_kind', 'index_value', 'value_kind', 'value'])
        for i, entry in enumerate(successful): writer.writerow(['int', '-2' if i % 3 == 0 else '7', *entry['wire']])
    with (args.output / 'values.expected.csv').open('x', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['index_kind', 'index_value', 'value_kind', 'value'])
        for i, entry in enumerate(successful): writer.writerow(['int', '-2' if i % 3 == 0 else '7', 'float', entry['bits']])
    for i, entry in enumerate(errors):
        entry['input'] = f'failure_{i}.input.csv'
        with (args.output / entry['input']).open('x', newline='') as stream:
            writer = csv.writer(stream); writer.writerow(['index_kind', 'index_value', 'value_kind', 'value'])
            writer.writerow(['str', '6475706c6963617465', *entry['wire']])
    diagnostics = []
    for name, value in [('bytes', b'.5'), ('list', [1, 2]), ('empty_list', []), ('array', np.array([1., 2.])),
                        ('zero_dim_array', np.array(.5)), ('dict', {'x': 1}), ('complex', 1 + 2j),
                        ('numpy_float16', np.float16(.5)), ('numpy_longdouble', np.longdouble('.1'))]:
        before = pickle.dumps(value, protocol=5)
        with warnings.catch_warnings(record=True) as warning:
            warnings.simplefilter('always')
            try: outcome = {'status': 'value', 'bits': bits(original(value))}
            except Exception as error: outcome = {'status': 'error', 'error_type': type(error).__name__, 'error_message': str(error)}
        assert pickle.dumps(value, protocol=5) == before
        diagnostics.append(dict(name=name, unsupported=True, input_unchanged=True,
                                warnings=[type(item.message).__name__ for item in warning], **outcome))
    manifest = dict(schema='priva.coerce-pext-python-oracle.v1', source_commit=COMMIT, source_sha256=SOURCE_SHA,
        native=dict(job_id=os.environ.get('SLURM_JOB_ID'), step_id=os.environ.get('SLURM_STEP_ID')),
        function_source_sha256=FUNCTION_SHA, function_lines=[function.lineno, function.end_lineno],
        runtime=dict(python=platform.python_version(), pandas=pd.__version__, numpy=np.__version__, unicode=unicodedata.unidata_version),
        harness_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), decimal_zeroes=zeroes,
        successful_cases=successful, source_failures=errors, unsupported_diagnostics=diagnostics,
        execution='unchanged actual source AST with real pandas/NumPy; scalar input pickle bytes unchanged',
        production_backend_selected=False, scientific_acceptance=False)
    manifest['files'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(args.output.iterdir()) if p.is_file()}
    (args.output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(dict(values=len(successful), source_errors=len(errors), unsupported_diagnostics=len(diagnostics),
                         runtime=manifest['runtime'], decimal_zeroes=len(zeroes),
                         manifest_sha256=hashlib.sha256((args.output / 'manifest.json').read_bytes()).hexdigest())))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    command = parser.add_subparsers(dest='command', required=True).add_parser('capture')
    command.add_argument('--source-file', type=Path, required=True)
    command.add_argument('--output', type=Path, required=True)
    capture(parser.parse_args())


if __name__ == '__main__':
    main()
