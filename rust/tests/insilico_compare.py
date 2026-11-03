"""Compare the complete optional insilico module with SAVED genuine originals.

No original procedure is re-executed. Unknown whole profiles are reported as
original-only admission cases, not replacement parity. Native and explicitly
mixed operation profiles are counted separately, with every raw mismatch saved.
This driver is only a procedure; run it through a separately reviewed frozen
Slurm controller/source manifest. It never builds or activates a binary.
"""
import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import numpy as np
import pandas as pd
import acmg_pp3_bp4_bp7_insilico as reference
import priva_rust_insilico as native
import priva_rust_insilico_source as guards

CAPTURE_SHA = '3f622a4199d7b5132824fae56b43ca570c683aa8fcbd626a6f2e942ef767613e'
CASES_SHA = '1980e1121e92822da760e3aee9b7b32634665f17a6abade16307b21815bd36c3'
PLAN_SHA = 'e744dc760d56c8b9c6354043110ccdb76dd9986629a608cfe65e58207083bd57'
RESULT_SHAS = {'0': 'e407a84214f5efe03e4fb5fa8b6ca5ac56092a3ff7bd28a38ed5828117f27af6',
               '1': '76f07c3dbd12cfa13dc6b46a2145093e6a27308f8c04fb5f95c2b6599cf2228e'}


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path, digest, name):
    if sha(path) != digest: raise RuntimeError('saved authority source changed')
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def save(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write('\n')


def protected(outcome):
    error = outcome['error']
    # Physical adapter/native transport frames remain retained in raw evidence.
    # They are not asserted identical to the unmodified Python call stack.
    return dict(value=outcome['value'], error=None if error is None else
                {name: error[name] for name in ('type', 'message', 'args')},
                events=outcome['events'], warnings=outcome['warnings'],
                input_before=outcome['input_before'], input_after=outcome['input_after'])


def compare(args):
    saved = Path(args.saved_root)
    output = Path(args.output)
    if output.exists(): raise RuntimeError('preserve every prior comparison')
    if sha(args.binary) != args.binary_sha256: raise RuntimeError('exact qualified native ELF required')
    if sha(saved / 'plan.json') != PLAN_SHA: raise RuntimeError('exact saved original plan required')
    plan = json.loads((saved / 'plan.json').read_bytes())
    old_result_path = saved / ('seed-' + args.seed) / 'result.json'
    if sha(old_result_path) != RESULT_SHAS[args.seed]: raise RuntimeError('exact genuine seed result required')
    old_result = json.loads(old_result_path.read_bytes())
    expected = plan['specifications']
    if (not old_result['completed'] or old_result['seed'] != args.seed
            or old_result['plan_sha256'] != PLAN_SHA
            or len(old_result['cases']) != len(expected)
            or [(row['id'], row['group'], row['procedure']) for row in old_result['cases']]
                != [(row['id'], row['group'], row['procedure']) for row in expected]):
        raise RuntimeError('complete ordered saved original inventory required')
    capture = load(saved / 'controller.py', CAPTURE_SHA, 'saved_insilico_capture')
    fixtures = load(saved / 'cases.py', CASES_SHA, 'saved_insilico_fixtures')
    # Retain the actual canonical namespace/path used by the saved producer,
    # including warning filenames. Loading definitions is not an original call;
    # the entry guard below forbids executing either original procedure.
    reference = load(plan['reference']['path'], guards.REFERENCE_SHA,
                     'acmg_pp3_bp4_bp7_insilico')
    output.mkdir()
    bridge = native.Bridge(binary=args.binary, binary_sha256=args.binary_sha256,
                           reference=reference, allow_original_fallback=False)
    # The observer sees the same real logger/default objects. Only its callable
    # entrypoints are a proxy; the canonical original module/code is never
    # patched, executed, re-captured or substituted in source admission.
    def pp3(frame, *values, **kwargs):
        return bridge.call('PP3_BP4_criteria', frame, *values, **kwargs)
    pp3.__defaults__ = reference.PP3_BP4_criteria.__defaults__
    def bp7(frame): return bridge.call('BP7_criteria', frame)
    proxy = SimpleNamespace(logger=reference.logger, PP3_BP4_criteria=pp3,
                            BP7_criteria=bp7, __file__=reference.__file__)
    summary = dict(completed=False, seed=args.seed, cases=[], mismatches=0,
        admission_only=0, exact_by_mode=Counter(), actual_processes=0,
        native_operations=Counter(), mature_operations=Counter(), new_original_calls=0,
        limits='Finite saved-oracle scope. Admission-only cases are NOT replacement parity. Native/mixed execution and unchanged Python fallback remain distinct. Physical stack/line-trace identities are retained, not asserted identical.')
    original_codes = {reference.PP3_BP4_criteria.__code__, reference.BP7_criteria.__code__}
    # Keep its actual getter profile untouched. A local sys proxy wraps only
    # the driver's trace installation and forwards exactly its returned trace
    # decision. A forbidden original entry is stopped BEFORE its body, retained
    # as a blocking observer failure, never allowed as a hidden baseline rerun.
    original_entries = []
    class TraceProxy:
        def __getattr__(self, name): return getattr(sys, name)
        @staticmethod
        def settrace(callback):
            if callback is None:
                sys.settrace(None)
                return
            def guarded(frame, event, value):
                if event == 'call' and frame.f_code in original_codes:
                    original_entries.append(frame.f_code.co_name)
                    raise RuntimeError('forbidden original procedure reexecution')
                following = callback(frame, event, value)
                return guarded if following is not None else None
            sys.settrace(guarded)
    old_sys = capture.sys
    capture.sys = TraceProxy()
    try:
        for spec, pin in zip(expected, old_result['cases']):
            old_path = saved / ('seed-' + args.seed) / (spec['id'] + '.json')
            if pin['pin']['path'] != str(old_path) or sha(old_path) != pin['pin']['sha256']:
                raise RuntimeError('saved outcome path/content changed')
            old = json.loads(old_path.read_bytes())
            if (old['id'], old['group'], old['procedure']) != (spec['id'], spec['group'], spec['procedure']):
                raise RuntimeError('saved outcome identity changed')
            actual = capture.observe_case(spec, fixtures, proxy, np, pd)
            if original_entries:
                raise RuntimeError('original entry was attempted; no acceptance')
            metrics = bridge.last_metrics
            if metrics['mode'] == 'original':
                status = 'unsupported-original-only; no original rerun'
                summary['admission_only'] += 1
                differences = []
            else:
                before, after = protected(old), protected(actual)
                differences = [name for name in before if before[name] != after[name]]
                status = 'mismatch' if differences else 'exact'
                summary['mismatches'] += bool(differences)
                if not differences: summary['exact_by_mode'][metrics['mode']] += 1
                summary['actual_processes'] += metrics['processes']
                summary['native_operations'].update(metrics['native_operations'])
                summary['mature_operations'].update(metrics['mature_fallbacks'])
            save(output / (spec['id'] + '.json'), dict(id=spec['id'], status=status,
                differences=differences, saved_original_sha256=pin['pin']['sha256'],
                original=old if differences else None, actual=actual, metrics=metrics))
            summary['cases'].append(dict(id=spec['id'], group=spec['group'], procedure=spec['procedure'],
                                         status=status, differences=differences))
        if sha(reference.__file__) != guards.REFERENCE_SHA or sha(args.binary) != args.binary_sha256:
            raise RuntimeError('reference/native source changed during comparison')
        summary['completed'] = True
    finally:
        capture.sys = old_sys
        summary['new_original_entry_attempts'] = original_entries
        save(output / 'result.json', summary)
    print(json.dumps({key: value for key, value in summary.items() if key != 'cases'}, sort_keys=True))
    return 1 if summary['mismatches'] else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    run = commands.add_parser('compare')
    for name in ('saved-root', 'binary', 'binary-sha256', 'output'): run.add_argument('--' + name, required=True)
    run.add_argument('--seed', choices=('0', '1'), required=True)
    args = parser.parse_args()
    raise SystemExit(compare(args))


if __name__ == '__main__': main()
