#!/usr/bin/env python3
"""Compare complete original ranking with optional R1/R2/R3/R4/R5-R6/R8 bridges.

compare --reference SOURCE_ROOT --candidate SOURCE_ROOT --binary FILE
        --output NEW_DIRECTORY
worker --source SOURCE_ROOT --reference SOURCE_ROOT --cases FILE
       --mode python|native --profile fresh|source-bootstrap --output FILE
       [--binary FILE]

The Python reference is frozen before the adapter exists. C1-C5 and R7 stay
source-executed; each preceding ranking stage retains its exact fallback. Outputs retain ordered axes, dtype/null storage, IEEE bits,
attributes and label flags; original failures and logs are compared as well.
R4 replaces only the original mechanism-helper call at its source position;
raw call traces remain recorded and its one explicit replacement is checked.
"""
import argparse
import ast
from contextlib import nullcontext
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
from unittest.mock import patch

import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _ranker_definition(tree):
    definitions = [node for node in tree.body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                   and node.name == 'sort_and_rank_variants']
    assert len(definitions) == 1
    definition = definitions[0]
    assert isinstance(definition, ast.FunctionDef) and not definition.decorator_list
    return definition


def source_proof(reference, candidate):
    original = ast.parse((reference / 'scripts/acmg_scoring.py').read_text())
    updated = ast.parse((candidate / 'scripts/acmg_scoring.py').read_text())
    ref = _ranker_definition(original)
    new = _ranker_definition(updated)
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
    context_dispatch = ast.parse('''native_context = None
if ranking_binary is not None:
    from priva_rust_pedigree_context import try_pedigree_context
    native_context = try_pedigree_context(df, ped_df, fam_name, ranking_binary)
''').body
    assert ast.dump(ast.Module(body=new.body[4:6], type_ignores=[]), include_attributes=False) == ast.dump(
        ast.Module(body=context_dispatch, type_ignores=[]), include_attributes=False)
    context = new.body[6]
    assert isinstance(context, ast.If) and ast.unparse(context.test) == 'native_context is None'
    assert ast.dump(ast.Module(body=context.body, type_ignores=[]), include_attributes=False) == ast.dump(
        ast.Module(body=[ref.body[4]], type_ignores=[]), include_attributes=False)
    accepted = ast.parse('''proband = native_context['proband']
proband_row = native_context['proband_row']
proband_is_male = native_context['proband_is_male']
proband_is_female = native_context['proband_is_female']
proband_het = native_context['proband_het']
if native_context['source_warning'] == 'unknown_sex':
    pass
elif native_context['source_warning'] == 'no_pedigree':
    pass
''').body
    # Warning calls are taken from the frozen original statement, including
    # their real argument expressions and order; no reconstructed log text.
    accepted[-1].body = [copy.deepcopy(ref.body[4].body[5].body[0])]
    accepted[-1].orelse[0].body = [copy.deepcopy(ref.body[4].orelse[0])]
    assert ast.dump(ast.Module(body=context.orelse, type_ignores=[]), include_attributes=False) == ast.dump(
        ast.Module(body=accepted, type_ignores=[]), include_attributes=False)
    new.body[4:7] = context.body
    product_dispatch = ast.parse('''native_product = None
if ranking_binary is not None:
    from priva_rust_control_multiply import try_control_multiply
    native_product = try_control_multiply(df, ranking_binary)
''').body
    assert ast.dump(ast.Module(body=new.body[6:8], type_ignores=[]), include_attributes=False) == ast.dump(
        ast.Module(body=product_dispatch, type_ignores=[]), include_attributes=False)
    product = new.body[8]
    assert isinstance(product, ast.If) and ast.unparse(product.test) == 'native_product is None'
    assert ast.dump(ast.Module(body=product.body, type_ignores=[]), include_attributes=False) == ast.dump(
        ast.Module(body=[ref.body[6]], type_ignores=[]), include_attributes=False)
    assert ast.dump(ast.Module(body=product.orelse, type_ignores=[]), include_attributes=False) == ast.dump(
        ast.parse('df.loc[:, "sort_index"] = native_product'), include_attributes=False)
    new.body[6:9] = product.body
    compatibility_dispatch = ast.parse('''native_compatibility = None
if ranking_binary is not None:
    from priva_rust_ranking_compatibility import try_ranking_compatibility
    native_compatibility = try_ranking_compatibility(
        df, proband_is_male, proband_is_female, proband_het, ranking_binary)
''').body
    assert ast.dump(ast.Module(body=new.body[7:9], type_ignores=[]), include_attributes=False) == ast.dump(
        ast.Module(body=compatibility_dispatch, type_ignores=[]), include_attributes=False)
    compatibility = new.body[9]
    assert isinstance(compatibility, ast.If) and ast.unparse(compatibility.test) == 'native_compatibility is None'
    assert ast.dump(ast.Module(body=compatibility.body, type_ignores=[]), include_attributes=False) == ast.dump(
        ast.Module(body=ref.body[7:18], type_ignores=[]), include_attributes=False)
    roles = ('mechanism_masks', 'chrom', 'x_locus', 'has_single_alt_lof_compatible', 'exact_gof',
             'asserted_lof_effect', 'dn_without_hi_or_ambiguous', 'heterozygous_lof_incompatible',
             'high_acmg_lof', 'compatibility_penalty_rows', 'pext_modulation_rows')
    accepted_compatibility = ast.parse('\n'.join(
        name + " = native_compatibility['" + name + "']" for name in roles)).body
    assert ast.dump(ast.Module(body=compatibility.orelse, type_ignores=[]), include_attributes=False) == ast.dump(
        ast.Module(body=accepted_compatibility, type_ignores=[]), include_attributes=False)
    new.body[7:10] = compatibility.body
    weights_dispatch = ast.parse('''native_weights = None
if ranking_binary is not None:
    from priva_rust_gene_weights import try_gene_weights, gene_membership
    native_weights = try_gene_weights(
        df, dispensable_gene_list, relevant_gene_list, compatibility_penalty_rows, ranking_binary)
''').body
    assert ast.dump(ast.Module(body=new.body[21:23], type_ignores=[]), include_attributes=False) == ast.dump(
        ast.Module(body=weights_dispatch, type_ignores=[]), include_attributes=False)
    weights = new.body[23]
    assert isinstance(weights, ast.If) and ast.unparse(weights.test) == 'native_weights is None'
    assert ast.dump(ast.Module(body=weights.body, type_ignores=[]), include_attributes=False) == ast.dump(
        ast.Module(body=ref.body[21:24], type_ignores=[]), include_attributes=False)
    accepted_weights = ast.parse('''if dispensable_gene_list is not None:
    dispensable_genes = native_weights['dispensable_genes']
    df["dispensable_gene_index"] = 1.0
    df.loc[gene_membership(df, native_weights['dispensable_selected']), "dispensable_gene_index"] = 0.9
    df.loc[:, "sort_index"] = native_weights['dispensable_product']
if relevant_gene_list is not None:
    relevant_genes = native_weights['relevant_genes']
    df["relevant_gene_index"] = 1.0
    df.loc[gene_membership(df, native_weights['relevant_selected']), "relevant_gene_index"] = 1.2
    df.loc[:, "sort_index"] = native_weights['relevant_product']
df.loc[compatibility_penalty_rows, "sort_index"] = native_weights['compatibility_product']
''').body
    assert ast.dump(ast.Module(body=weights.orelse, type_ignores=[]), include_attributes=False) == ast.dump(
        ast.Module(body=accepted_weights, type_ignores=[]), include_attributes=False)
    new.body[21:24] = weights.body
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
                all_original_statements_and_prior_source_order_unchanged=True,
                original_R2_statement_and_warning_calls_exact=True,
                original_R3_statement_and_control_call_order_exact=True,
                original_R4_block_INFO_and_two_assignment_order_exact=True,
                original_R5_R6_statements_and_three_product_order_exact=True)


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


def _source_call_projection(raw, replacements, compatibility_calls):
    """Permit only the recorded native replacement at the original R4 slot."""
    assert len(compatibility_calls) <= 1
    if not replacements:
        assert not any(compatibility_calls)
        return raw
    assert compatibility_calls == [True]
    assert replacements == [dict(position=1, name='_variant_mechanism_masks')]
    assert raw and raw[0] == 'compute_control_common' and '_variant_mechanism_masks' not in raw
    assert all(name == 'coerce_pext_value' for name in raw[1:])
    return [raw[0], '_variant_mechanism_masks', *raw[1:]]


def _gene_parser_projection(raw, replacements, calls, enabled):
    """Keep raw parser calls; permit only explicit enabled native replacements."""
    assert len(calls) <= 1
    assert enabled in ([], ['dispensable'], ['relevant'], ['dispensable', 'relevant'])
    if calls == [True]:
        assert not raw and replacements == enabled
        return list(enabled)
    assert not replacements
    return raw


def call_case(scoring, case, mode, binary, helper, native_stages='r8'):
    # A fresh load retains the actual prepared input, without normalizing dtypes
    # or reading metadata before the source function starts.
    data, pedigree, kwargs = pickle.loads(pickle.dumps(
        (case['frame'], case['ped'], case['kwargs']), protocol=5))
    source_calls = []
    native_calls = []
    setup_calls = []
    pedigree_calls = []
    product_calls = []
    compatibility_calls = []
    compatibility_replacements = []
    r2_contexts = []
    r3_outputs = []
    r4_outputs = []
    r4_entry_seen = False
    gene_calls, gene_parser_calls, gene_parser_replacements = [], [], []
    r5_snapshots, gene_enabled = [], []
    r5_entry_seen, r6_finished = False, False
    source = ast.parse(Path(scoring.__file__).read_text())
    ranker = _ranker_definition(source)
    info_lines = [node.lineno for node in ranker.body if isinstance(node, ast.Expr)
        and isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Attribute)
        and node.value.func.attr == 'info' and node.value.args
        and isinstance(node.value.args[0], ast.JoinedStr) and any(
            isinstance(part, ast.Constant) and isinstance(part.value, str) and
            part.value.startswith('Zygosity/inheritance/mechanism compatibility signals for ranking: ')
            for part in node.value.args[0].values)]
    assert len(info_lines) == 1
    ranker_code = scoring.sort_and_rank_variants.__code__
    weight_dispatch = [node for node in ranker.body if isinstance(node, ast.If)
                       and ast.unparse(node.test) == 'native_weights is None']
    original_weights = weight_dispatch[0].body if weight_dispatch else [
        node for node in ranker.body if isinstance(node, ast.If)
        and ast.unparse(node.test) in ('dispensable_gene_list is not None', 'relevant_gene_list is not None')]
    if not weight_dispatch:
        original_weights.append(ranker.body[ranker.body.index(original_weights[-1]) + 1])
    assert len(original_weights) == 3
    entry_line = next(node.lineno for node in ranker.body if isinstance(node, ast.Assign)
                      and any(isinstance(target, ast.Name) and target.id == 'native_weights' for target in node.targets)) if weight_dispatch else original_weights[0].lineno
    checkpoint_lines = {original_weights[1].lineno: 'after_dispensable',
                        original_weights[2].lineno: 'after_relevant'}
    if weight_dispatch:
        checkpoint_lines.update({weight_dispatch[0].orelse[1].lineno: 'after_dispensable',
                                 weight_dispatch[0].orelse[2].lineno: 'after_relevant'})
    following_line = next(node.lineno for node in ranker.body if isinstance(node, ast.Assign)
                          and any(isinstance(target, ast.Name) and target.id == 'pext_col_to_use' for target in node.targets))
    parser_lines = {branch.body[0].lineno: role for branch, role in
                    zip(original_weights[:2], ('dispensable', 'relevant'))}
    read_table_code = scoring.pd.read_table.__code__
    def snapshot(local, stage):
        keys = ('df', 'dispensable_gene_list', 'relevant_gene_list',
                'compatibility_penalty_rows', 'dispensable_genes', 'relevant_genes')
        # Observe without dtype encoding inside the source execution. In
        # particular, encoding nullable metadata here could warm pandas state
        # before the exact multiplication being compared.
        values = {name: local[name] for name in keys if name in local}
        values['mask_index_is_frame_index'] = local['compatibility_penalty_rows'].index is local['df'].index
        r5_snapshots.append((stage, pickle.dumps(values, protocol=5)))
    def trace(pyframe, event, arg):
        nonlocal r5_entry_seen, r6_finished
        if pyframe.f_code is ranker_code:
            if event == 'line' and pyframe.f_lineno == entry_line and not r5_entry_seen:
                r5_entry_seen = True
                gene_enabled.extend(role for role, name in (('dispensable', 'dispensable_gene_list'),
                    ('relevant', 'relevant_gene_list')) if pyframe.f_locals[name] is not None)
                snapshot(pyframe.f_locals, 'entry')
            elif event == 'line' and r5_entry_seen and pyframe.f_lineno in checkpoint_lines:
                snapshot(pyframe.f_locals, checkpoint_lines[pyframe.f_lineno])
            elif event == 'line' and pyframe.f_lineno == following_line and r5_entry_seen:
                snapshot(pyframe.f_locals, 'after_compatibility'); r6_finished = True
            elif event == 'return' and r5_entry_seen and not r6_finished:
                snapshot(pyframe.f_locals, 'failure')
        return trace
    def profile(pyframe, event, arg):
        nonlocal r4_entry_seen
        if event == 'call' and pyframe.f_code is read_table_code and pyframe.f_back is not None and (
                pyframe.f_back.f_code is ranker_code):
            assert pyframe.f_back.f_lineno in parser_lines
            gene_parser_calls.append(parser_lines[pyframe.f_back.f_lineno])
        if event == 'call' and pyframe.f_code.co_name in (
                'compute_control_common', '_variant_mechanism_masks', 'coerce_pext_value'):
            source_calls.append(pyframe.f_code.co_name)
            if pyframe.f_code.co_name == 'compute_control_common':
                local = pyframe.f_back.f_locals
                keys = ['proband', 'proband_is_male', 'proband_is_female', 'proband_het']
                if local['ped_df'] is not None and local['fam_name'] is not None:
                    keys.append('proband_row')
                r2_contexts.append({key: helper.encode(local[key]) for key in keys})
            elif pyframe.f_code.co_name == '_variant_mechanism_masks' and (
                    pyframe.f_back.f_code.co_name == 'sort_and_rank_variants'):
                # Observe actual R4 entry after R3's real assignment, without
                # warming pandas metadata before the original multiplication.
                if not r4_entry_seen:
                    r3_outputs.append(helper.encode(pyframe.f_back.f_locals['df']))
                    r4_entry_seen = True
        if event == 'call' and pyframe.f_code.co_name == 'try_ranking_compatibility' and (
                Path(pyframe.f_code.co_filename).name == 'priva_rust_ranking_compatibility.py'):
            assert not r4_entry_seen
            r3_outputs.append(helper.encode(pyframe.f_locals['frame']))
            r4_entry_seen = True
        if event == 'call' and pyframe.f_code.co_name == 'info' and pyframe.f_back is not None and (
                pyframe.f_back.f_code.co_name == 'sort_and_rank_variants'
                and pyframe.f_back.f_lineno == info_lines[0]):
            local = pyframe.f_back.f_locals
            roles = ('df', 'proband_is_male', 'proband_is_female', 'proband_het',
                'mechanism_masks', 'chrom', 'x_locus', 'has_single_alt_lof_compatible', 'exact_gof',
                'asserted_lof_effect', 'dn_without_hi_or_ambiguous', 'heterozygous_lof_incompatible',
                'high_acmg_lof', 'compatibility_penalty_rows', 'pext_modulation_rows')
            r4_outputs.append(dict(locals={name: helper.encode(local[name]) for name in roles}, aliases=[
                local['exact_gof'] is local['mechanism_masks']['is_exact_gof'],
                local['asserted_lof_effect'] is local['mechanism_masks']['is_predicted_lof'],
                local['pext_modulation_rows'] is local['high_acmg_lof']]))
        if event == 'return' and pyframe.f_code.co_name == 'try_terminal_order' and (
                Path(pyframe.f_code.co_filename).name == 'priva_rust_terminal_order.py'):
            native_calls.append(arg is not None)
        if event == 'return' and pyframe.f_code.co_name == 'try_rank_setup' and (
                Path(pyframe.f_code.co_filename).name == 'priva_rust_rank_setup.py'):
            setup_calls.append(arg is not None)
        if event == 'return' and pyframe.f_code.co_name == 'try_pedigree_context' and (
                Path(pyframe.f_code.co_filename).name == 'priva_rust_pedigree_context.py'):
            pedigree_calls.append(arg is not None)
        if event == 'return' and pyframe.f_code.co_name == 'try_control_multiply' and (
                Path(pyframe.f_code.co_filename).name == 'priva_rust_control_multiply.py'):
            product_calls.append(arg is not None)
        if event == 'return' and pyframe.f_code.co_name == 'try_ranking_compatibility' and (
                Path(pyframe.f_code.co_filename).name == 'priva_rust_ranking_compatibility.py'):
            compatibility_calls.append(arg is not None)
            if arg is not None:
                compatibility_replacements.append(dict(position=len(source_calls), name='_variant_mechanism_masks'))
        if event == 'return' and pyframe.f_code.co_name == 'try_gene_weights' and (
                Path(pyframe.f_code.co_filename).name == 'priva_rust_gene_weights.py'):
            gene_calls.append(arg is not None)
            if arg is not None:
                gene_parser_replacements.extend(gene_enabled)
    if mode == 'native':
        # Observe returns through the existing profiler. Do not import the
        # adapter early: the real ranker first imports it at its R8 entry after
        # original C1-C5/R1-R7, including in the first cold-context native call.
        kwargs['ranking_binary' if native_stages != 'r8' else 'terminal_order_binary'] = str(binary)
    # Explicit test control only. The integrated public argument always enables
    # every implemented stage; this labelled mode forces R2's original fallback
    # to keep the preceding R1/R8 control available without a production switch.
    force_r2 = mode == 'native' and native_stages == 'r1-r8'
    control = patch('priva_rust_pedigree_context.try_pedigree_context', return_value=None) if force_r2 else nullcontext(None)
    # Keep each earlier staged control genuine, without a production switch or
    # claiming the newly integrated R3 algorithm ran in its preceding modes.
    force_r3 = mode == 'native' and native_stages in ('r1-r8', 'r1-r2-r8')
    product_control = patch('priva_rust_control_multiply.try_control_multiply', return_value=None) if force_r3 else nullcontext(None)
    force_r4 = mode == 'native' and native_stages in ('r1-r8', 'r1-r2-r8', 'r1-r2-r3-r8')
    compatibility_control = patch('priva_rust_ranking_compatibility.try_ranking_compatibility', return_value=None) if force_r4 else nullcontext(None)
    force_genes = mode == 'native' and native_stages not in ('r8', 'r1-r2-r3-r4-r5-r6-r8')
    gene_control = patch('priva_rust_gene_weights.try_gene_weights', return_value=None) if force_genes else nullcontext(None)
    log = io.StringIO()
    handler = logging.StreamHandler(log)
    handler.setFormatter(logging.Formatter('%(levelname)s %(name)s %(message)s'))
    root_logger = logging.getLogger()
    old_level = root_logger.level
    root_logger.setLevel(logging.INFO)
    root_logger.addHandler(handler)
    previous, previous_trace = sys.getprofile(), sys.gettrace()
    output, error = None, None
    try:
        with warnings.catch_warnings(record=True) as caught, control as forced, product_control as forced_product, compatibility_control as forced_compatibility, gene_control as forced_genes:
            warnings.simplefilter('always')
            sys.setprofile(profile)
            sys.settrace(trace)
            try:
                output = scoring.sort_and_rank_variants(data, pedigree, **kwargs)
            except Exception as failure:
                error = dict(type=helper.typename(failure), message=str(failure), args=helper.encode(failure.args),
                    source_trace=[dict(function=row.name, statement=row.line) for row in
                        traceback.extract_tb(failure.__traceback__) if Path(row.filename).name in
                        ['acmg_scoring.py', 'acmg_consequence.py', 'acmg_variant_mechanism.py']])
            finally:
                sys.settrace(previous_trace)
                sys.setprofile(previous)
    finally:
        root_logger.removeHandler(handler)
        root_logger.setLevel(old_level)
    forced_calls = forced.call_count if force_r2 else 0
    forced_product_calls = forced_product.call_count if force_r3 else 0
    forced_compatibility_calls = forced_compatibility.call_count if force_r4 else 0
    return dict(name=case['name'], output=helper.encode(output), exception=error,
        input_after=helper.encode(data), pedigree_after=helper.encode(pedigree), source_calls=source_calls,
        warnings=[dict(type=helper.typename(value.message), message=str(value.message)) for value in caught],
        log=log.getvalue(), native_calls=native_calls, setup_calls=setup_calls,
        pedigree_calls=pedigree_calls, forced_r2_fallback_calls=forced_calls, r2_contexts=r2_contexts,
        product_calls=product_calls, forced_r3_fallback_calls=forced_product_calls, r3_outputs=r3_outputs,
        compatibility_calls=compatibility_calls, forced_r4_fallback_calls=forced_compatibility_calls,
        compatibility_replacements=compatibility_replacements, r4_outputs=r4_outputs,
        gene_calls=gene_calls, forced_gene_fallback_calls=forced_genes.call_count if force_genes else 0,
        gene_parser_calls=gene_parser_calls, gene_parser_replacements=gene_parser_replacements,
        gene_enabled=gene_enabled,
        r5_states=[dict(stage=stage, locals=helper.encode(pickle.loads(value))) for stage, value in r5_snapshots])


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
    if args.native_stages != 'r8':
        setup = load('unchanged_R1_fixture_builder', args.candidate / 'rust/tests/ranking_setup_oracle.py')
        cases.extend(dict(name='r1_complete_' + name, kind='ranking', frame=frame, ped=None,
                          kwargs=dict(fam_name=None, dispensable_gene_list=None),
                          tags=['complete-ranker-r1-integration']) for name, frame in setup.extras(helper))
    if args.native_stages in ('r1-r2-r8', 'r1-r2-r3-r8', 'r1-r2-r3-r4-r8', 'r1-r2-r3-r4-r5-r6-r8'):
        pedigree = load('unchanged_R2_fixture_builder', args.candidate / 'rust/tests/pedigree_context_oracle.py')
        for name, entry, ped, family in pedigree.extras():
            # Lift the pre-existing standalone R2 inputs into a full original
            # ranker frame; retain every R2 column, dtype, axis and scalar.
            data = helper.base_frame(len(entry))
            data.index = entry.index
            data = data.drop(columns=[column for column in data.columns if column in entry.columns])
            data = pd.concat([data, entry], axis=1)
            data.columns.name = entry.columns.name
            data.attrs = copy.deepcopy(entry.attrs)
            data.flags.allows_duplicate_labels = entry.flags.allows_duplicate_labels
            cases.append(dict(name='r2_complete_' + name, kind='ranking', frame=data, ped=ped,
                kwargs=dict(fam_name=family, dispensable_gene_list=None),
                tags=['complete-ranker-r2-integration']))
    if args.native_stages in ('r1-r2-r3-r8', 'r1-r2-r3-r4-r8', 'r1-r2-r3-r4-r5-r6-r8'):
        multiply = load('unchanged_R3_fixture_builder', args.candidate / 'rust/tests/control_multiply_oracle.py')
        for name, entry in multiply.examples():
            # Lift the existing standalone seeds through the complete ranker.
            # R1 derives sort_index from the seed's original operand dtype;
            # actual C1-C5 may overwrite control_common_index as originally
            # specified. Never substitute a precomputed control return.
            data = helper.base_frame(len(entry))
            data.index = entry.index
            data = data.drop(columns=[column for column in data.columns if column in entry.columns])
            data = pd.concat([data, entry], axis=1)
            if 'sort_index' in entry.columns:
                position = list(entry.columns).index('sort_index')
                data['ACMG_quant_score'] = entry.iloc[:, position]
            data.columns.name = entry.columns.name
            data.attrs = copy.deepcopy(entry.attrs)
            data.flags.allows_duplicate_labels = entry.flags.allows_duplicate_labels
            cases.append(dict(name='r3_complete_' + name, kind='ranking', frame=data, ped=None,
                kwargs=dict(fam_name=None, dispensable_gene_list=None),
                tags=['complete-ranker-r3-seed-through-original-R1-and-C1-C5']))
    if args.native_stages in ('r1-r2-r3-r4-r8', 'r1-r2-r3-r4-r5-r6-r8'):
        compatibility = load('unchanged_R4_fixture_builder', args.candidate / 'rust/tests/ranking_compatibility_oracle.py')
        for name, (entry, male, female, heterozygosity) in compatibility.extras():
            # Lift the historical frame seeds through actual R1/R2/C1-C5/R3.
            # The prior standalone mask is a seed for a real sample column;
            # R2 recreates its own current mask, including its raw null data.
            data = helper.base_frame(len(entry))
            data.index = entry.index
            data = data.drop(columns=[column for column in data.columns if column in entry.columns])
            data = pd.concat([data, entry], axis=1)
            sample = ['0/1' if value is True or isinstance(value, np.bool_) and value else
                      None if value is None or value is pd.NA else '0/0' for value in heterozygosity]
            data['R4_PROBAND'] = pd.Series(sample, index=data.index,
                dtype='string' if str(heterozygosity.dtype) == 'boolean' else 'object')
            data.columns.name = entry.columns.name
            data.attrs = copy.deepcopy(entry.attrs)
            data.flags.allows_duplicate_labels = entry.flags.allows_duplicate_labels
            sex = '1' if male and not female else '2' if female and not male else 'unknown'
            ped = pd.DataFrame({'#FamilyID': ['R4'], 'IndividualID': ['R4_PROBAND'],
                                'Phenotype': ['2'], 'Sex': [sex]})
            cases.append(dict(name='r4_complete_' + name, kind='ranking', frame=data, ped=ped,
                kwargs=dict(fam_name='R4', dispensable_gene_list=None),
                tags=['complete-ranker-r4-seed-through-original-upstream-stages',
                      'current-R2-mask-recomputed-by-original-source']))
    if args.native_stages == 'r1-r2-r3-r4-r5-r6-r8':
        genes = load('unchanged_R5_R6_fixture_builder', args.candidate / 'rust/tests/ranking_gene_weights_oracle.py')
        for name, (entry, dispensable, relevant, penalty) in genes.extras(args.output):
            # Lift historical seeds through every actual preceding stage. The
            # old stage-only mask seeds the sample; R2/R4 recompute the actual
            # current mask. Direct nullable/IEEE storage is proved separately.
            data = helper.base_frame(len(entry))
            # R1/R4 regenerate scores and compatibility, but never SYMBOL.
            # Preserve the old gene seed's omitted consumed column.
            if 'SYMBOL' not in entry.columns:
                data = data.drop(columns='SYMBOL')
            data.index = entry.index
            data = data.drop(columns=[column for column in data.columns if column in entry.columns])
            data = pd.concat([data, entry], axis=1)
            if 'sort_index' in entry.columns:
                position = list(entry.columns).index('sort_index')
                data['ACMG_quant_score'] = entry.iloc[:, position]
            sample = ['0/1' if value is True or isinstance(value, np.bool_) and value else
                      None if value is None or value is pd.NA else '0/0' for value in penalty]
            data['R5_PROBAND'] = pd.Series(sample, index=data.index,
                dtype='string' if str(penalty.dtype) == 'boolean' else 'object')
            data.columns.name = entry.columns.name
            data.attrs = copy.deepcopy(entry.attrs)
            data.flags.allows_duplicate_labels = entry.flags.allows_duplicate_labels
            ped = pd.DataFrame({'#FamilyID': ['R5'], 'IndividualID': ['R5_PROBAND'],
                                'Phenotype': ['2'], 'Sex': ['1']})
            cases.append(dict(name='r5_complete_' + name, kind='ranking', frame=data, ped=ped,
                kwargs=dict(fam_name='R5', dispensable_gene_list=dispensable, relevant_gene_list=relevant),
                tags=['complete-ranker-r5-r6-seed-through-original-upstream-stages',
                      'current-R2-and-R4-mask-recomputed-by-original-source']))
    prepared = args.output / 'cases.pkl'
    prepared.write_bytes(pickle.dumps(cases, protocol=5))
    reports, comparisons, mismatches = {}, [], []
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
                    pedigree_calls = observed.pop('pedigree_calls')
                    forced_r2_calls = observed.pop('forced_r2_fallback_calls')
                    product_calls = observed.pop('product_calls')
                    forced_r3_calls = observed.pop('forced_r3_fallback_calls')
                    compatibility_calls = observed.pop('compatibility_calls')
                    forced_r4_calls = observed.pop('forced_r4_fallback_calls')
                    replacements = observed.pop('compatibility_replacements')
                    gene_calls = observed.pop('gene_calls')
                    forced_gene_calls = observed.pop('forced_gene_fallback_calls')
                    parser_replacements = observed.pop('gene_parser_replacements')
                    raw_parser_calls = observed['gene_parser_calls']
                    observed['gene_parser_calls'] = _gene_parser_projection(
                        raw_parser_calls, parser_replacements, gene_calls, observed['gene_enabled'])
                    raw_source_calls = observed['source_calls']
                    observed['source_calls'] = _source_call_projection(raw_source_calls, replacements, compatibility_calls)
                    wanted = {key: value for key, value in expected.items() if key not in (
                        'native_calls', 'setup_calls', 'pedigree_calls', 'forced_r2_fallback_calls',
                        'product_calls', 'forced_r3_fallback_calls', 'compatibility_calls',
                        'forced_r4_fallback_calls', 'compatibility_replacements', 'gene_calls',
                        'forced_gene_fallback_calls', 'gene_parser_replacements')}
                    exact = observed == wanted
                    if not exact:
                        mismatches.append(dict(profile=profile, mode=label, name=expected['name'],
                            fields=[key for key in wanted if observed[key] != wanted[key]],
                            reference_report='reference-' + profile + '.json',
                            candidate_report=label + '-' + profile + '.json'))
                    comparisons.append(dict(name=expected['name'], profile=profile, mode=label,
                        exact=exact, native_used=bool(native_calls and native_calls[0]),
                        fallback=bool(native_calls and not native_calls[0]),
                        setup_native_used=bool(setup_calls and setup_calls[0]),
                        setup_fallback=bool(setup_calls and not setup_calls[0]),
                        pedigree_native_used=bool(pedigree_calls and pedigree_calls[0]),
                        pedigree_fallback=bool(pedigree_calls and not pedigree_calls[0]),
                        forced_r2_fallback_calls=forced_r2_calls,
                        product_native_used=bool(product_calls and product_calls[0]),
                        product_fallback=bool(product_calls and not product_calls[0]),
                        forced_r3_fallback_calls=forced_r3_calls,
                        compatibility_native_used=bool(compatibility_calls and compatibility_calls[0]),
                        compatibility_fallback=bool(compatibility_calls and not compatibility_calls[0]),
                        forced_r4_fallback_calls=forced_r4_calls,
                        actual_source_calls=raw_source_calls, explicit_native_R4_replacements=replacements,
                        gene_native_used=bool(gene_calls and gene_calls[0]),
                        gene_fallback=bool(gene_calls and not gene_calls[0]),
                        forced_gene_fallback_calls=forced_gene_calls,
                        actual_gene_parser_calls=raw_parser_calls, explicit_native_gene_parser_replacements=parser_replacements,
                        original_exception=expected['exception'] is not None,
                        earlier_failure_did_not_reach_R8=expected['exception'] is not None and not native_calls))
        native = [item for item in comparisons if item['mode'] == 'native']
        assert any(item['native_used'] for item in native), 'no genuine native complete-function result'
        assert any(item['fallback'] for item in native), 'no genuine current-buffer Python fallback'
        assert any(item['earlier_failure_did_not_reach_R8'] for item in native)
        if args.native_stages != 'r8':
            assert any(item['setup_native_used'] for item in native), 'no genuine native R1 complete-function result'
            assert any(item['setup_fallback'] for item in native), 'no genuine current-buffer R1 fallback'
        if args.native_stages == 'r1-r8':
            assert any(item['forced_r2_fallback_calls'] for item in native), 'R2 control did not execute the original fallback'
            assert not any(item['pedigree_native_used'] for item in native), 'R1/R8 control published native R2'
        if args.native_stages in ('r1-r2-r8', 'r1-r2-r3-r8', 'r1-r2-r3-r4-r8', 'r1-r2-r3-r4-r5-r6-r8'):
            assert any(item['pedigree_native_used'] for item in native), 'no genuine native R2 complete-function result'
            assert any(item['pedigree_fallback'] for item in native), 'no genuine current-buffer R2 fallback'
            assert not any(item['forced_r2_fallback_calls'] for item in native), 'full R2 mode used a forced fallback'
        if args.native_stages in ('r1-r8', 'r1-r2-r8'):
            assert any(item['forced_r3_fallback_calls'] for item in native), 'R3 control did not execute the original fallback'
            assert not any(item['product_native_used'] for item in native), 'preceding control published native R3'
        if args.native_stages in ('r1-r2-r3-r8', 'r1-r2-r3-r4-r8', 'r1-r2-r3-r4-r5-r6-r8'):
            assert any(item['product_native_used'] for item in native), 'no genuine native R3 complete-function product'
            assert any(item['product_fallback'] for item in native), 'no genuine current-buffer R3 fallback'
            assert not any(item['forced_r3_fallback_calls'] for item in native), 'full R3 mode used a forced fallback'
        if args.native_stages in ('r1-r8', 'r1-r2-r8', 'r1-r2-r3-r8'):
            assert any(item['forced_r4_fallback_calls'] for item in native), 'R4 control did not execute the original fallback'
            assert not any(item['compatibility_native_used'] for item in native), 'preceding control published native R4'
        if args.native_stages in ('r1-r2-r3-r4-r8', 'r1-r2-r3-r4-r5-r6-r8'):
            assert any(item['compatibility_native_used'] for item in native), 'no genuine native R4 complete-function plan'
            assert any(item['compatibility_fallback'] for item in native), 'no genuine current-buffer R4 fallback'
            assert not any(item['forced_r4_fallback_calls'] for item in native), 'full R4 mode used a forced fallback'
        if args.native_stages not in ('r8', 'r1-r2-r3-r4-r5-r6-r8'):
            assert any(item['forced_gene_fallback_calls'] for item in native), 'R5/R6 control did not execute original fallback'
            assert not any(item['gene_native_used'] for item in native), 'preceding control published native R5/R6'
        if args.native_stages == 'r1-r2-r3-r4-r5-r6-r8':
            assert any(item['gene_native_used'] for item in native), 'no genuine native R5/R6 complete-function plan'
            assert any(item['gene_fallback'] for item in native), 'no genuine current-buffer R5/R6 fallback'
            assert not any(item['forced_gene_fallback_calls'] for item in native), 'full R5/R6 mode forced fallback'
        assert not mismatches, ('complete-ranker mismatches', len(mismatches), mismatches)
    except Exception as error:
        failure = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        native = [item for item in comparisons if item['mode'] == 'native']
        summary = dict(successful=failure is None, failure=failure, source_proof=proof,
            mismatches=mismatches,
            cases=len(cases), profiles=['fresh', 'source-bootstrap'], exact_comparisons=len(comparisons),
            native_outputs=sum(item['native_used'] for item in native),
            Python_fallbacks=sum(item['fallback'] for item in native),
            native_stages=args.native_stages,
            setup_native_outputs=sum(item['setup_native_used'] for item in native),
            setup_Python_fallbacks=sum(item['setup_fallback'] for item in native),
            pedigree_native_outputs=sum(item['pedigree_native_used'] for item in native),
            pedigree_Python_fallbacks=sum(item['pedigree_fallback'] for item in native),
            forced_R2_control_fallback_calls=sum(item['forced_r2_fallback_calls'] for item in native),
            control_multiply_native_outputs=sum(item['product_native_used'] for item in native),
            control_multiply_Python_fallbacks=sum(item['product_fallback'] for item in native),
            forced_R3_control_fallback_calls=sum(item['forced_r3_fallback_calls'] for item in native),
            compatibility_native_outputs=sum(item['compatibility_native_used'] for item in native),
            compatibility_Python_fallbacks=sum(item['compatibility_fallback'] for item in native),
            forced_R4_control_fallback_calls=sum(item['forced_r4_fallback_calls'] for item in native),
            explicit_native_R4_helper_replacements=sum(len(item['explicit_native_R4_replacements']) for item in native),
            gene_weights_native_outputs=sum(item['gene_native_used'] for item in native),
            gene_weights_Python_fallbacks=sum(item['gene_fallback'] for item in native),
            forced_R5_R6_control_fallback_calls=sum(item['forced_gene_fallback_calls'] for item in native),
            explicit_native_gene_parser_replacements=sum(len(item['explicit_native_gene_parser_replacements']) for item in native),
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
    command.add_argument('--native-stages', choices=['r8', 'r1-r8', 'r1-r2-r8', 'r1-r2-r3-r8', 'r1-r2-r3-r4-r8', 'r1-r2-r3-r4-r5-r6-r8'], default='r8')
    command = sub.add_parser('worker')
    for name in ['source', 'reference', 'cases', 'output']:
        command.add_argument('--' + name, required=True, type=Path)
    command.add_argument('--binary', type=Path)
    command.add_argument('--mode', required=True, choices=['python', 'native'])
    command.add_argument('--profile', required=True, choices=['fresh', 'source-bootstrap'])
    command.add_argument('--native-stages', choices=['r8', 'r1-r8', 'r1-r2-r8', 'r1-r2-r3-r8', 'r1-r2-r3-r4-r8', 'r1-r2-r3-r4-r5-r6-r8'], default='r8')
    args = parser.parse_args()
    globals()[args.procedure](args)


if __name__ == '__main__':
    main()
