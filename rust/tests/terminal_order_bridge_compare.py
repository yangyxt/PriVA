#!/usr/bin/env python3
"""Compare complete original ranking with the optional native R8 bridge.

compare --reference SOURCE_ROOT --candidate SOURCE_ROOT --binary FILE
        --output NEW_DIRECTORY
worker --source SOURCE_ROOT --reference SOURCE_ROOT --cases FILE
       --mode python|native --profile fresh|source-bootstrap --output FILE
       [--binary FILE]

The Python reference is frozen before the adapter exists. C1-C5 and R1-R7 stay
source-executed. Outputs retain ordered axes, dtype/null storage, IEEE bits,
attributes and label flags; original failures and logs are compared as well.
"""
import argparse
import ast
import copy
import hashlib
import importlib.util
import io
import json
import logging
import os
from pathlib import Path
import pickle
import subprocess
import sys
import traceback
import warnings

import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def source_proof(reference, candidate):
    original = ast.parse((reference / 'scripts/acmg_scoring.py').read_text())
    updated = ast.parse((candidate / 'scripts/acmg_scoring.py').read_text())
    ref = next(node for node in original.body if isinstance(node, ast.FunctionDef)
               and node.name == 'sort_and_rank_variants')
    new = next(node for node in updated.body if isinstance(node, ast.FunctionDef)
               and node.name == ref.name)
    ref, new = copy.deepcopy(ref), copy.deepcopy(new)
    assert [node.arg for node in new.args.args[-2:]] == ['terminal_order_binary', 'ranking_binary']
    for _ in range(2):
        assert isinstance(new.args.defaults[-1], ast.Constant) and new.args.defaults[-1].value is None
        new.args.args.pop(); new.args.defaults.pop()
    ref.body.pop(0); new.body.pop(0)  # Only the API docstring gains an argument.
    # The R1 opt-in keeps every original setup statement in its fallback branch.
    # Prove that exact branch and the small dispatch, rather than accepting an
    # arbitrary new prefix or relabeling the entire ranker as an unchanged file.
    setup_dispatch = ast.parse('''native_setup = None
if ranking_binary is not None:
    from priva_rust_rank_setup import try_rank_setup
    native_setup = try_rank_setup(df, ranking_binary)
''').body
    assert ast.dump(ast.Module(body=new.body[:2], type_ignores=[]), include_attributes=False) == ast.dump(
        ast.Module(body=setup_dispatch, type_ignores=[]), include_attributes=False)
    setup = new.body[2]
    assert isinstance(setup, ast.If) and ast.unparse(setup.test) == 'native_setup is None'
    assert ast.dump(ast.Module(body=setup.orelse, type_ignores=[]), include_attributes=False) == ast.dump(
        ast.parse('df = native_setup'), include_attributes=False)
    assert ast.dump(ast.Module(body=setup.body, type_ignores=[]), include_attributes=False) == ast.dump(
        ast.Module(body=ref.body[:4], type_ignores=[]), include_attributes=False)
    new.body[:3] = setup.body
    hooks = [node for node in new.body if isinstance(node, ast.If)
             and ast.unparse(node.test) == 'terminal_order_binary is not None or ranking_binary is not None']
    assert len(hooks) == 1
    hook = hooks[0]
    position = new.body.index(hook)
    assert isinstance(new.body[position + 1], ast.Assign)
    assert ast.unparse(new.body[position + 1].value) == "df.groupby(['chrom', 'pos', 'ref', 'alt'])"
    new.body.remove(hook)
    assert ast.dump(ref, include_attributes=False) == ast.dump(new, include_attributes=False)
    # The original function and its full numerical dependencies are otherwise
    # byte-identical; no stripped/reordered helper or NumPy/SciPy substitute.
    for name in ['acmg_consequence.py', 'acmg_variant_mechanism.py']:
        assert sha(reference / 'scripts' / name) == sha(candidate / 'scripts' / name)
    original.body = [node for node in original.body if not (isinstance(node, ast.FunctionDef) and node.name == ref.name)]
    updated.body = [node for node in updated.body if not (isinstance(node, ast.FunctionDef) and node.name == ref.name)]
    assert ast.dump(original, include_attributes=False) == ast.dump(updated, include_attributes=False)
    return dict(reference_scoring_sha256=sha(reference / 'scripts/acmg_scoring.py'),
                candidate_scoring_sha256=sha(candidate / 'scripts/acmg_scoring.py'),
                all_original_statements_and_prior_source_order_unchanged=True)


def extra_cases(helper):
    def item(name, frame, pedigree=None, **kwargs):
        options = dict(fam_name='F' if pedigree is not None else None, dispensable_gene_list=None)
        options.update(kwargs)
        return dict(name='bridge_' + name, kind='ranking', frame=frame, ped=pedigree,
                    kwargs=options, tags=['complete-ranker-integration'])
    for name, index in [
        ('multi_index', pd.MultiIndex.from_tuples([('b', 2), ('a', 9), ('b', 2), ('z', -3)], names=['first', 'second'])),
        ('reverse_range', pd.RangeIndex(9, 1, -2, name='reverse')),
        ('duplicate_index', pd.Index(['same', None, 'same', ''], dtype=object, name='labels'))]:
        data = helper.base_frame(4, ACMG_quant_score=[.7, .9, .9, .8])
        data.index = index
        data.columns.name = 'original-column-axis'
        data['payload'] = [[1, 2], {'a': 1}, None, ['raw']]
        yield item(name, data)
    data = helper.base_frame(4)
    data.flags.allows_duplicate_labels = False
    yield item('label_flags', data)
    data = helper.base_frame(3)
    data['chrom'] = pd.Categorical(['2', '1', '2'], categories=['2', '1'], ordered=True)
    yield item('categorical_python_fallback', data)
    data = helper.base_frame(3)
    data['pos'] = pd.Series([100, pd.NA, 102], dtype=object)
    yield item('object_pd_na_coordinate', data)
    data = helper.base_frame(3, chrom=['a:b', 'a', 'a'], pos=[1., np.nan, 2.]).drop(columns='variant_id')
    yield item('absent_id_float_colon_fallback', data)
    data = helper.base_frame(4, ACMG_quant_score=[-0., 0., np.inf, -np.inf], pos=[1, 1, 2, 2])
    data.insert(1, 'max_variant_score', ['old'] * 4)
    data.insert(2, 'variant_rank', ['old'] * 4)
    yield item('overwrite_columns_and_nonfinite', data)


def call_case(scoring, case, mode, binary, helper, native_stages='r8'):
    # A fresh load retains the actual prepared input, without normalizing dtypes
    # or reading metadata before the source function starts.
    data, pedigree, kwargs = pickle.loads(pickle.dumps(
        (case['frame'], case['ped'], case['kwargs']), protocol=5))
    source_calls = []
    native_calls = []
    setup_calls = []
    def profile(pyframe, event, arg):
        if event == 'call' and pyframe.f_code.co_name in (
                'compute_control_common', '_variant_mechanism_masks', 'coerce_pext_value'):
            source_calls.append(pyframe.f_code.co_name)
        if event == 'return' and pyframe.f_code.co_name == 'try_terminal_order' and (
                Path(pyframe.f_code.co_filename).name == 'priva_rust_terminal_order.py'):
            native_calls.append(arg is not None)
        if event == 'return' and pyframe.f_code.co_name == 'try_rank_setup' and (
                Path(pyframe.f_code.co_filename).name == 'priva_rust_rank_setup.py'):
            setup_calls.append(arg is not None)
    if mode == 'native':
        # Observe returns through the existing profiler. Do not import the
        # adapter early: the real ranker first imports it at its R8 entry after
        # original C1-C5/R1-R7, including in the first cold-context native call.
        kwargs['ranking_binary' if native_stages == 'r1-r8' else 'terminal_order_binary'] = str(binary)
    log = io.StringIO()
    handler = logging.StreamHandler(log)
    handler.setFormatter(logging.Formatter('%(levelname)s %(name)s %(message)s'))
    root_logger = logging.getLogger()
    old_level = root_logger.level
    root_logger.setLevel(logging.INFO)
    root_logger.addHandler(handler)
    previous = sys.getprofile()
    output, error = None, None
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            sys.setprofile(profile)
            try:
                output = scoring.sort_and_rank_variants(data, pedigree, **kwargs)
            except Exception as failure:
                error = dict(type=helper.typename(failure), message=str(failure), args=helper.encode(failure.args),
                    source_trace=[dict(function=row.name, statement=row.line) for row in
                        traceback.extract_tb(failure.__traceback__) if Path(row.filename).name in
                        ['acmg_scoring.py', 'acmg_consequence.py', 'acmg_variant_mechanism.py']])
            finally:
                sys.setprofile(previous)
    finally:
        root_logger.removeHandler(handler)
        root_logger.setLevel(old_level)
    return dict(name=case['name'], output=helper.encode(output), exception=error,
        input_after=helper.encode(data), pedigree_after=helper.encode(pedigree), source_calls=source_calls,
        warnings=[dict(type=helper.typename(value.message), message=str(value.message)) for value in caught],
        log=log.getvalue(), native_calls=native_calls, setup_calls=setup_calls)


def worker(args):
    sys.path.insert(0, str(args.source / 'scripts'))
    helper = load('complete_ranking_fixture_helper', args.reference / 'rust/tests/ranking_control_oracle.py')
    import acmg_scoring as scoring
    assert Path(scoring.__file__).parent == args.source / 'scripts'
    cases = pickle.loads(args.cases.read_bytes())
    if args.profile == 'source-bootstrap':
        bootstrap = load('unchanged_source_bootstrap', args.reference / 'scripts/acmg_scoring.py')
        # This executes separate original source on separate input buffers; the
        # measured current-stage input is never rehydrated or precomputed.
        for case in cases:
            call_case(bootstrap, case, 'python', None, helper)
    reports = [call_case(scoring, case, args.mode, args.binary, helper, args.native_stages) for case in cases]
    with args.output.open('x') as stream:
        json.dump(dict(cases=reports, versions=dict(python=sys.version, numpy=np.__version__, pandas=pd.__version__),
            native={key: os.environ.get(key) for key in ['SLURM_JOB_ID', 'SLURM_STEP_ID', 'SLURMD_NODENAME']}), stream, indent=2)
        stream.write('\n')


def compare(args):
    args.output.mkdir(parents=True, exist_ok=False)
    proof = source_proof(args.reference, args.candidate)
    helper = load('complete_ranking_fixture_builder', args.reference / 'rust/tests/ranking_control_oracle.py')
    matrix, _ = helper.case_matrix(args.output)
    cases = [case for case in matrix if case['kind'] == 'ranking'] + list(extra_cases(helper))
    if args.native_stages == 'r1-r8':
        setup = load('unchanged_R1_fixture_builder', args.candidate / 'rust/tests/ranking_setup_oracle.py')
        cases.extend(dict(name='r1_complete_' + name, kind='ranking', frame=frame, ped=None,
                          kwargs=dict(fam_name=None, dispensable_gene_list=None),
                          tags=['complete-ranker-r1-integration']) for name, frame in setup.extras(helper))
    prepared = args.output / 'cases.pkl'
    prepared.write_bytes(pickle.dumps(cases, protocol=5))
    reports, comparisons = {}, []
    failure = None
    try:
        for profile in ['fresh', 'source-bootstrap']:
            for label, source, mode in [('reference', args.reference, 'python'),
                                         ('default', args.candidate, 'python'),
                                         ('native', args.candidate, 'native')]:
                name = label + '-' + profile
                output = args.output / (name + '.json')
                command = [sys.executable, '-B', __file__, 'worker', '--source', str(source),
                    '--reference', str(args.reference), '--cases', str(prepared), '--mode', mode,
                    '--profile', profile, '--output', str(output), '--binary', str(args.binary),
                    '--native-stages', args.native_stages]
                with (args.output / (name + '.log')).open('x') as stream:
                    result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
                assert result.returncode == 0, (name, result.returncode)
                reports[name] = json.loads(output.read_text())
            reference = reports['reference-' + profile]['cases']
            for label in ['default', 'native']:
                actual = reports[label + '-' + profile]['cases']
                assert len(actual) == len(reference) == len(cases)
                for expected, observed in zip(reference, actual):
                    native_calls = observed.pop('native_calls')
                    setup_calls = observed.pop('setup_calls')
                    wanted = {key: value for key, value in expected.items() if key not in ('native_calls', 'setup_calls')}
                    assert observed == wanted, (profile, label, expected['name'],
                        [key for key in wanted if observed[key] != wanted[key]])
                    comparisons.append(dict(name=expected['name'], profile=profile, mode=label,
                        exact=True, native_used=bool(native_calls and native_calls[0]),
                        fallback=bool(native_calls and not native_calls[0]),
                        setup_native_used=bool(setup_calls and setup_calls[0]),
                        setup_fallback=bool(setup_calls and not setup_calls[0]),
                        original_exception=expected['exception'] is not None,
                        earlier_failure_did_not_reach_R8=expected['exception'] is not None and not native_calls))
        native = [item for item in comparisons if item['mode'] == 'native']
        assert any(item['native_used'] for item in native), 'no genuine native complete-function result'
        assert any(item['fallback'] for item in native), 'no genuine current-buffer Python fallback'
        assert any(item['earlier_failure_did_not_reach_R8'] for item in native)
        if args.native_stages == 'r1-r8':
            assert any(item['setup_native_used'] for item in native), 'no genuine native R1 complete-function result'
            assert any(item['setup_fallback'] for item in native), 'no genuine current-buffer R1 fallback'
    except Exception as error:
        failure = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        native = [item for item in comparisons if item['mode'] == 'native']
        summary = dict(successful=failure is None, failure=failure, source_proof=proof,
            cases=len(cases), profiles=['fresh', 'source-bootstrap'], exact_comparisons=len(comparisons),
            native_outputs=sum(item['native_used'] for item in native),
            Python_fallbacks=sum(item['fallback'] for item in native),
            native_stages=args.native_stages,
            setup_native_outputs=sum(item['setup_native_used'] for item in native),
            setup_Python_fallbacks=sum(item['setup_fallback'] for item in native),
            original_failures=sum(item['original_exception'] for item in native),
            earlier_failures_without_native_call=sum(item['earlier_failure_did_not_reach_R8'] for item in native),
            binary_sha256=sha(args.binary), fixture_sha256=sha(prepared),
            family_acceptance=False, public_acceptance=False, production_default_changed=False)
        (args.output / 'comparison.json').write_text(json.dumps(summary, indent=2) + '\n')
        (args.output / 'case-comparisons.json').write_text(json.dumps(comparisons, indent=2) + '\n')
        print(json.dumps(summary), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='procedure', required=True)
    command = sub.add_parser('compare')
    for name in ['reference', 'candidate', 'binary', 'output']:
        command.add_argument('--' + name, required=True, type=Path)
    command.add_argument('--native-stages', choices=['r8', 'r1-r8'], default='r8')
    command = sub.add_parser('worker')
    for name in ['source', 'reference', 'cases', 'output']:
        command.add_argument('--' + name, required=True, type=Path)
    command.add_argument('--binary', type=Path)
    command.add_argument('--mode', required=True, choices=['python', 'native'])
    command.add_argument('--profile', required=True, choices=['fresh', 'source-bootstrap'])
    command.add_argument('--native-stages', choices=['r8', 'r1-r8'], default='r8')
    args = parser.parse_args()
    globals()[args.procedure](args)


if __name__ == '__main__':
    main()
