#!/usr/bin/env python3
"""Portable complete-module PM4/BP3 differential comparison.

Use --binary PATH [--output-root NEW_DIRECTORY]. Existing saved original
outcomes can be supplied with --reference-outcomes PATH; worker execution stays
source-bound and unsupported/mature-tool failures remain exact. Generated inputs,
logs and outcomes belong in the supplied output directory, never the repository.
"""
import argparse
import base64
import contextlib
import hashlib
import importlib.util
import io
import json
import logging
import os
from pathlib import Path
import resource
import shutil
import stat
import subprocess
import sys
import tempfile
import time
import traceback
import warnings

REFERENCE_SHA = '08ef5af382801363727c6aec4915ca5751267067b43ef7996e3715e3fd1bc515'
REPOSITORY = Path(__file__).resolve().parents[2]
REFERENCE_NAME = 'acmg_pm4_bp3_protein_length'
RAW_FIELDS = ('st_dev', 'st_ino', 'st_size', 'st_mtime_ns', 'st_ctime_ns', 'st_mode', 'st_nlink')


def pin(path):
    path = Path(path)
    before = path.lstat()
    assert stat.S_ISREG(before.st_mode), str(path)
    raw = path.read_bytes()
    after = path.lstat()
    assert [getattr(before, key) for key in RAW_FIELDS] == [getattr(after, key) for key in RAW_FIELDS]
    return dict(path=str(path), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest(),
                raw7=[getattr(after, key) for key in RAW_FIELDS])


def write_json(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, ensure_ascii=True, allow_nan=False)
        handle.write('\n')


def native():
    return {key: os.environ[key] for key in (
        'SLURM_JOB_ID', 'SLURM_STEP_ID', 'SLURM_CPUS_PER_TASK',
        'TZ', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'POLARS_MAX_THREADS') if key in os.environ}


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def worker(args):
    source = Path(args.source_root)
    work = Path(args.work_root)
    assert work.is_dir() and not any(work.iterdir())
    reference_path = source
    assert pin(reference_path)['sha256'] == REFERENCE_SHA
    os.environ['PATH'] = str(Path(sys.executable).parent) + os.pathsep + os.environ['PATH']
    import numpy as np
    import pandas as pd
    import pybedtools
    from pybedtools import helpers, settings
    cases = load_module('pm4_bp3_cases', REPOSITORY / 'tests/fixtures/pm4_bp3_cases.py')
    state = load_module('pm4_bp3_state', REPOSITORY / 'tests/fixtures/pm4_bp3_state.py')
    module = load_module(REFERENCE_NAME, reference_path)
    fixture = cases.build(args.case)
    df, pvs, loc = fixture['frame'], fixture['pvs'], fixture['loc']
    setup_callbacks = state.atom(list(cases.CALLBACKS))
    before = dict(frame=state.frame(df), pvs=state.output(pvs), loc=state.output(loc))
    cache_before = {key: id(value) for key, value in df._item_cache.items()}
    backing_before = [df._mgr.iget_values(i) for i in range(len(df.columns))]
    cases.CALLBACKS.clear()
    regions = work / 'regions.bed'
    if fixture['regions'] is not None:
        regions.write_bytes(fixture['regions'])

    logs, mature_calls, source_calls, tempnames = [], [], [], []
    logger = module.logger
    # Observe genuine records at Logger.handle rather than installing a custom
    # emit callback. The adapter admits only ordinary mature logging handlers.
    logger.handlers = [logging.StreamHandler(io.StringIO())]
    if fixture['log_failure'] is not None:
        logger.handlers.append(cases.RaiseAtLog(fixture['log_failure']))
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    prototype = None
    if args.prototype_root is not None:
        sys.path.insert(0, str(args.prototype_root))
        import priva_rust_pm4_bp3 as adapter
        binding, digest = adapter._package_profile()
        assert digest == args.package_binding_sha256
        prototype = adapter._Bridge(module, binary=args.binary, binary_sha256=args.binary_sha256,
            scratch=args.scratch, package_binding=binding,
            package_binding_sha256=args.package_binding_sha256)
    codes = {pybedtools.BedTool.from_dataframe.__code__: 'BedTool.from_dataframe',
             pybedtools.BedTool._tmp.__code__: 'BedTool._tmp',
             helpers.call_bedtools.__code__: 'call_bedtools'}
    reference_codes = {module.find_overlaps_bedtools_efficient.__code__: 'find_overlaps_bedtools_efficient',
                       module.PM4_BP3_criteria.__code__: 'PM4_BP3_criteria'}
    diagnostics = {}
    observer_failures = []

    def profile(frame, event, arg):
        try:
            if frame.f_code is logging.Logger.handle.__code__ and event == 'call':
                record = frame.f_locals['record']
                logs.append(dict(name=record.name, level=record.levelname, levelno=record.levelno,
                                 message=record.getMessage()))
            if frame.f_code in reference_codes and event in ('call', 'return'):
                source_calls.append(dict(event=event, function=reference_codes[frame.f_code]))
                if event == 'return' and frame.f_code is module.PM4_BP3_criteria.__code__:
                    diagnostics.update({key: state.output(value) for key, value in frame.f_locals.items()
                        if key not in ('df', 'pvs1_criteria', 'repeat_regions_file', 'loc_intol_domain')})
            label = codes.get(frame.f_code)
            if label is not None and event in ('call', 'return'):
                item = dict(event=event, function=label)
                if event == 'call' and label == 'call_bedtools':
                    item['cmds'] = list(frame.f_locals['cmds'])
                    item['tmpfn'] = frame.f_locals['tmpfn']
                if event == 'call' and label == 'BedTool.from_dataframe':
                    item['dataframe'] = state.frame(frame.f_locals['df'])
                if event == 'return' and label == 'BedTool._tmp' and type(arg) is str:
                    tempnames.append(arg)
                    item['path'] = arg
                if event == 'return' and label == 'BedTool.from_dataframe' and arg is not None:
                    item['path'] = arg.fn
                mature_calls.append(item)
        except BaseException:
            observer_failures.append(traceback.format_exc())
        return profile

    sequence = tempfile._get_candidate_names()
    rng = sequence.rng
    rng_before = rng.getstate()
    tempfile_before, keep_before = tempfile.tempdir, settings.KEEP_TEMPFILES
    profile_before = sys.getprofile()
    assert profile_before is None
    rng.seed(39271801)
    tempfile.tempdir = str(work)
    settings.KEEP_TEMPFILES = True
    np_before = np.geterr()
    np_call_before = np.geterrcall()
    np.seterr(**({key: fixture['np_errors'] for key in np_before} if fixture['np_errors'] else np_before))
    result, error, raw_traceback = None, None, None
    started = time.perf_counter()
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            sys.setprofile(profile)
            try:
                if args.case.startswith('overlap-'):
                    function = module.find_overlaps_bedtools_efficient if prototype is None else prototype.find_overlaps_bedtools_efficient
                    result = function(df, str(regions), method=fixture['method'])
                else:
                    function = module.PM4_BP3_criteria if prototype is None else prototype.PM4_BP3_criteria
                    result = function(df, pvs, str(regions), loc)
            except BaseException as exc:
                error = dict(type=state.kind(exc), message=str(exc), args=state.atom(exc.args))
                raw_traceback = traceback.format_exc()
            finally:
                sys.setprofile(profile_before)
            recorded_warnings = [dict(category=state.kind(w.message), message=str(w.message),
                                      filename=w.filename, lineno=w.lineno) for w in caught]
    finally:
        elapsed = time.perf_counter() - started
        sys.setprofile(profile_before)
        np_after, np_call_after = np.geterr(), np.geterrcall()
        np.seterr(**np_before)
        np.seterrcall(np_call_before)
        rng.setstate(rng_before)
        tempfile.tempdir = tempfile_before
        # The child owns this retention-only setting until it exits. Restoring
        # False here would let pybedtools' atexit cleanup delete captured files.
        settings.KEEP_TEMPFILES = True
    assert not observer_failures, observer_failures
    after = dict(frame=state.frame(df), pvs=state.output(pvs), loc=state.output(loc))
    alias_state = dict(cache_previous_objects_preserved={key: df._item_cache.get(key) is not None and
        id(df._item_cache[key]) == identity for key, identity in cache_before.items()},
        original_backing_buffers_still_shared=[bool(np.shares_memory(value, df._mgr.iget_values(i)))
            if type(value) is np.ndarray and type(df._mgr.iget_values(i)) is np.ndarray else
            value is df._mgr.iget_values(i) for i, value in enumerate(backing_before)])
    value = dict(captured=True, case=args.case, native=native(), source=pin(reference_path),
        versions=dict(python=sys.version, numpy=np.__version__, pandas=pd.__version__, pybedtools=pybedtools.__version__),
        function_status='raised' if error else 'returned', error=error, raw_traceback=raw_traceback,
        before=before, after=after, aliases=alias_state, output=state.output(result),
        logs=logs, warnings=recorded_warnings, callbacks=state.atom(list(cases.CALLBACKS)),
        fixture_setup_callbacks=setup_callbacks, mature_calls=mature_calls, source_calls=source_calls,
        prototype_events=None if prototype is None else prototype.events,
        prototype_masks=None if prototype is None or prototype.last_masks is None else
            {name: state.array(value) for name, value in prototype.last_masks.items()},
        diagnostic_locals_not_cross_language_contract=diagnostics, temporary_names=tempnames,
        pybedtools_tempfiles=list(pybedtools.BedTool.TEMPFILES),
        files=[pin(path) for path in sorted(work.iterdir()) if path.is_file()],
        numpy_state=dict(before=np_before, during=np_after,
            callback_unchanged=np_call_before is np_call_after),
        measurement=dict(function_wall_seconds=elapsed, process_maxrss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
        fixture_context=dict(tempfile_rng_state_restored=rng.getstate() == rng_before,
            tempfile_method_unchanged=True, mature_methods_unchanged=True,
            private_retention_only=True, keep_tempfiles_initial=keep_before,
            keep_tempfiles_until_isolated_process_exit=True, reference_callable_unmodified=True))
    write_json(args.output, value)



def differences(expected, actual, prefix=''):
    if type(expected) is not type(actual): return [dict(path=prefix, category='type', expected_type=type(expected).__name__, actual_type=type(actual).__name__)]
    if type(expected) is dict:
        if set(expected) != set(actual): return [dict(path=prefix, category='keys', expected=sorted(expected), actual=sorted(actual))]
        return [item for key in expected for item in differences(expected[key], actual[key], prefix + '/' + key)]
    if type(expected) is list:
        if len(expected) != len(actual): return [dict(path=prefix, category='length', expected=len(expected), actual=len(actual))]
        return [item for i, (left, right) in enumerate(zip(expected, actual)) for item in differences(left, right, prefix + '/' + str(i))]
    return [] if expected == actual else [dict(path=prefix, category='value', expected=expected, actual=actual)]


def logical_bool(value):
    if value['type'].endswith('.Series'): value = value['values']
    if value['type'].endswith('.BooleanArray'):
        data, mask = logical_bool(value['data']), logical_bool(value['mask'])
        assert not any(mask), 'unexpected missing private mask in supported source domain'
        return data
    assert value['type'] == 'numpy.ndarray' and value['dtype'] == '|b1'
    return list(base64.b64decode(value['raw_C_order_base64']))


def protected(expected, actual):
    result = [item for key in CONTRACT for item in differences(expected[key], actual[key], '/' + key)]
    source_files = [(Path(item['path']).name, item['bytes'], item['sha256']) for item in expected['files']]
    actual_files = [(Path(item['path']).name, item['bytes'], item['sha256']) for item in actual['files']]
    result += differences(source_files, actual_files, '/mature_file_names_sizes_sha256')
    if actual.get('prototype_masks') is not None:
        for name, actual_value in actual['prototype_masks'].items():
            expected_value = expected['diagnostic_locals_not_cross_language_contract'][name]
            result += differences(logical_bool(expected_value), logical_bool(actual_value), '/private_logical_mask/' + name)
    return result


def wrapper(target, binary, mode):
    path = target / (mode + '.py')
    program = '#!' + str(Path(sys.executable).resolve()) + '\nimport os,sys\n'
    if mode == 'prepare-exit': program += "if sys.argv[1:3]==['pm4-bp3','prepare-bed']: sys.exit(23)\n"
    elif mode == 'mask-exit': program += "if sys.argv[1:3]==['pm4-bp3','assign']: sys.exit(29)\n"
    elif mode == 'mask-corrupt':
        program += "if sys.argv[1:3]==['pm4-bp3','assign']:\n p=sys.argv[sys.argv.index('--output')+1]\n with open(p,'x') as h: h.write('invalid\\n0\\n')\n sys.exit(0)\n"
    else: raise ValueError(mode)
    program += 'os.execv(' + repr(str(binary)) + ', [' + repr(str(binary)) + '] + sys.argv[1:])\n'
    with path.open('x') as handle: handle.write(program)
    path.chmod(0o700)
    return path


CONTRACT = ('function_status', 'error', 'before', 'after', 'aliases', 'output',
    'logs', 'warnings', 'callbacks', 'fixture_setup_callbacks', 'mature_calls',
    'temporary_names', 'pybedtools_tempfiles')
FAULTS = (('mask-exit','whole-base'), ('mask-exit','whole-nullable-booleans'),
    ('mask-exit','whole-stale-id-cache'), ('prepare-exit','overlap-all'),
    ('prepare-exit','whole-base'), ('mask-corrupt','whole-base'),
    ('default','whole-base'), ('default','overlap-all'))
CI_CASES = ('whole-base', 'whole-string-python', 'whole-nullable-booleans',
    'whole-warm-cache', 'whole-stale-ref-cache', 'whole-stale-id-cache',
    'whole-stale-consequence-cache', 'overlap-any', 'overlap-all',
    'overlap-full-ref-boundary', 'overlap-invalid-regions')


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def compare(args):
    # These are mandatory test dependencies: a missing tool is a failed setup,
    # never a skip or a successful always-fallback native comparison.
    import pybedtools
    command = subprocess.run(['bedtools', '--version'], capture_output=True, check=False)
    if (command.returncode or command.stdout.strip() != b'bedtools v2.31.1'
            or pybedtools.__version__ != '0.10.0'):
        raise RuntimeError('PM4/BP3 comparison requires bedtools2.31.1 and pybedtools0.10.0')
    source = Path(args.reference_source).resolve()
    assert sha(source) == REFERENCE_SHA
    binary = Path(args.binary).resolve(strict=True); binary_sha = sha(binary)
    if args.output_root is None:
        target = Path(tempfile.mkdtemp(prefix='priva-pm4-bp3-compare-'))
    else:
        target = Path(args.output_root); target.mkdir()
    print(json.dumps(dict(comparison_output_root=str(target))), flush=True)
    cases = load_module('pm4_bp3_inventory_cases', REPOSITORY/'tests/fixtures/pm4_bp3_cases.py')
    inventory = cases.descriptors()
    if args.case:
        names=set(args.case); inventory=[row for row in inventory if row['name'] in names]
        assert {row['name'] for row in inventory}==names
    elif not args.all_cases:
        inventory=[row for row in inventory if row['name'] in CI_CASES]
        assert len(inventory)==len(CI_CASES)
    work = Path(args.fixture_work_root) if args.fixture_work_root else target/'active-private-fixture'
    assert not work.exists()
    scratch=target/'native-transport'; scratch.mkdir()
    scripts=REPOSITORY/'scripts'
    sys.path.insert(0, str(scripts))
    import priva_rust_pm4_bp3 as adapter
    binding,digest=adapter._package_profile()
    assert sha(binding)==digest
    wrappers={mode:wrapper(target,binary,mode) for mode in ('prepare-exit','mask-exit','mask-corrupt')}
    mismatches=[]; rows=[]
    result=dict(passed=False,comparisons=rows,mismatches=mismatches,binary_sha256=binary_sha,
        reference_source_sha256=sha(source),original_saved_outcomes_reused=args.reference_outcomes is not None)
    def execute(case,mode,out):
        out.mkdir(); work.mkdir()
        argv=[sys.executable,str(Path(__file__).resolve()),'worker','--case',case,
            '--source-root',str(source),'--work-root',str(work),'--output',str(out/'outcome.json')]
        if mode!='original':
            argv+=['--prototype-root',str(scripts),'--scratch',str(scratch),'--package-binding-sha256',digest]
            if mode!='default':
                selected=binary if mode=='native' else wrappers[mode]
                argv+=['--binary',str(selected),'--binary-sha256',sha(selected)]
        with (out/'stdout.log').open('x') as stdout,(out/'stderr.log').open('x') as stderr:
            done=subprocess.run(argv,stdout=stdout,stderr=stderr,timeout=120,check=False)
        work.rename(out/'retained-private-fixture')
        write_json(out/'relocation.json',dict(command=argv,returncode=done.returncode,child_reaped=True,
            original=str(work),retained=str(out/'retained-private-fixture')))
        assert done.returncode==0,(case,mode)
        return json.loads((out/'outcome.json').read_text())
    try:
        references={}
        for row in inventory:
            case=row['name']
            if args.reference_outcomes:
                path=Path(args.reference_outcomes)/'cases'/case/'outcome.json'
                references[case]=json.loads(path.read_text())
            else: references[case]=execute(case,'original',target/('original-'+case))
        selected=[('native',row['name']) for row in inventory]
        selected+=[(mode,case) for mode,case in FAULTS if case in references]
        for ordinal,(mode,case) in enumerate(selected):
            actual=execute(case,mode,target/('case-%03d-%s-%s'%(ordinal,mode,case)))
            differences=protected(references[case],actual)
            events=actual['prototype_events']
            native=any(value['event']=='native_whole_published' for value in events)
            if mode=='native' and case in ('whole-base','whole-string-python','whole-nullable-booleans',
                    'whole-warm-cache','whole-stale-ref-cache','whole-stale-id-cache','whole-stale-consequence-cache'):
                assert native,('supported fixture did not execute native procedure',case)
            if mode=='native' and case in ('overlap-any','overlap-all','overlap-unicode-ref',
                    'overlap-empty-ref','overlap-pos-zero','overlap-pos-negative',
                    'overlap-missing-regions','overlap-invalid-regions','overlap-full-ref-boundary'):
                assert any(value['event']=='native_result' and value.get('procedure')=='prepare-bed'
                    for value in events),('supported helper did not execute native procedure',case)
            intersections=sum(value['event']=='call' and value['function']=='call_bedtools' for value in actual['mature_calls'])
            if mode in ('mask-exit','mask-corrupt'):
                assert intersections==1
                assert any(value['event']=='original_tail_after_one_intersection' for value in events)
                assert not any(value['event']=='original_overlap_before_bedtools' for value in events)
            rows.append(dict(case=case,mode=mode,exact=not differences,native_whole_published=native,
                bedtools_intersections=intersections))
            mismatches.extend(dict(case=case,mode=mode,**value) for value in differences)
        assert not mismatches,'protected whole-module mismatch'
        result['passed']=True
    finally:
        write_json(target/'comparison.json',result)
    print(json.dumps(dict(passed=True,comparisons=len(rows),native_whole=sum(row['native_whole_published'] for row in rows))))


def main():
    parser=argparse.ArgumentParser(description=__doc__); commands=parser.add_subparsers(dest='procedure',required=True)
    command=commands.add_parser('worker')
    for name in ('case','source-root','work-root','output'): command.add_argument('--'+name,required=True)
    command.add_argument('--prototype-root',type=Path); command.add_argument('--binary',type=Path)
    command.add_argument('--binary-sha256'); command.add_argument('--package-binding-sha256'); command.add_argument('--scratch',type=Path)
    command=commands.add_parser('compare'); command.add_argument('--binary',type=Path,required=True)
    command.add_argument('--output-root',type=Path)
    command.add_argument('--reference-source',type=Path,default=REPOSITORY/'scripts/acmg_pm4_bp3_protein_length.py')
    command.add_argument('--reference-outcomes',type=Path); command.add_argument('--fixture-work-root',type=Path)
    command.add_argument('--case',action='append')
    command.add_argument('--all-cases',action='store_true')
    arguments=sys.argv[1:]
    if not arguments or arguments[0].startswith('--'):
        arguments=['compare',*arguments]
    args=parser.parse_args(arguments); globals()[args.procedure](args)


if __name__=='__main__': main()
