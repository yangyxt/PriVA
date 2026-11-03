"""Portable admission, restoration, state and failure controls.

Real native numerical/whole-unit comparisons require the explicit binary
driver; these tests never count a default or forced fallback as native proof.
"""
import ast
import copy
from pathlib import Path
import sys
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import priva_rust_reported_criteria as bridge
import priva_rust_reported_criteria_source as source


class DefaultAndAdmission(unittest.TestCase):
    def test_all_four_defaults_call_original_without_source_or_package_reads(self):
        arguments = {
            'mavedb_interpretation_per_row': (object(), object()),
            'mavedb_score_interpretation': (object(), object()),
            'PS3_BS3_criteria': (object(), '', {}),
            'PP5_BP6_criteria': (object(), object(), object()),
        }
        for name, args in arguments.items():
            marker = object()
            with self.subTest(name=name), mock.patch.object(bridge.original, name, return_value=marker) as original:
                with mock.patch.object(bridge, 'guard', side_effect=AssertionError('default read source')):
                    self.assertIs(getattr(bridge, name)(*args), marker)
                original.assert_called_once_with(*args)

    def test_rejected_admission_calls_original_once(self):
        arguments = (object(), object())
        marker = object()
        with mock.patch.object(bridge, 'guard', side_effect=ValueError('unsupported source')):
            with mock.patch.object(bridge.original, 'mavedb_interpretation_per_row', return_value=marker) as original:
                self.assertIs(bridge.mavedb_interpretation_per_row(*arguments, ranking_binary='unused'), marker)
                original.assert_called_once_with(*arguments)

    def test_default_confidence_mapping_uses_the_original_default_binding(self):
        marker = object()
        with mock.patch.object(bridge.original, 'PS3_BS3_criteria', return_value=marker) as original:
            with mock.patch.object(bridge, 'guard', side_effect=AssertionError('default read source')):
                self.assertIs(bridge.PS3_BS3_criteria(marker), marker)
            original.assert_called_once_with(marker, '')

    def test_original_exception_is_not_retried(self):
        failure = RuntimeError('actual source callback failure')
        with mock.patch.object(bridge, 'guard', side_effect=ValueError('unsupported')):
            with mock.patch.object(bridge.original, 'PP5_BP6_criteria', side_effect=failure) as original:
                with self.assertRaisesRegex(RuntimeError, 'actual source callback failure'):
                    bridge.PP5_BP6_criteria(object(), object(), object(), ranking_binary='unused')
                self.assertEqual(original.call_count, 1)

    def test_error_during_started_function_has_no_original_replay(self):
        failure = bridge.NativeTransportError('after source mutation')
        fake = mock.Mock()
        raw = {bridge.original.__name__: b'', bridge.interpreter.__name__: b''}
        with mock.patch.object(bridge, 'guard', return_value=raw), mock.patch.object(bridge, '_admit'):
            with mock.patch.object(bridge, 'Session', return_value=fake), mock.patch.object(bridge, 'compile_effects', return_value={}):
                with mock.patch.object(bridge.HostEffects, 'run', side_effect=failure):
                    with mock.patch.object(bridge.original, 'PP5_BP6_criteria') as original:
                        with self.assertRaisesRegex(bridge.NativeTransportError, 'after source mutation'):
                            bridge.PP5_BP6_criteria(object(), object(), object(), ranking_binary='unused')
                        original.assert_not_called()
        fake.close.assert_called_once()


class SourceRestoration(unittest.TestCase):
    def test_actual_original_guard_genuinely_accepts_supported_sources(self):
        raw = source.guard(bridge.original, bridge.interpreter)
        self.assertEqual(set(raw), set(source.SOURCE_SHA))

    def test_whole_public_source_and_mature_effect_inventory_remain_protected(self):
        raws = {module.__name__: Path(module.__file__).read_bytes() for module in (bridge.original, bridge.interpreter)}
        filenames = {module.__name__: module.__file__ for module in (bridge.original, bridge.interpreter)}
        before = {name: getattr(bridge.original, name).__code__ for name in source.PUBLIC}
        compiled = source.compile_effects(raws, filenames)
        self.assertTrue(set(source.EFFECTS).issubset(compiled))
        self.assertTrue(all(getattr(bridge.original, name).__code__ is code for name, code in before.items()))
        self.assertNotIn('exec', {node.func.id for node in ast.walk(ast.parse(Path(bridge.__file__).read_bytes()))
            if type(node) is ast.Call and type(node.func) is ast.Name})
        self.assertFalse(any(type(node) is ast.BinOp and type(node.op) in (ast.BitAnd, ast.BitOr)
            for _function, text in source.EFFECTS.values() for node in ast.walk(ast.parse(text, mode='eval'))))

    def test_missing_public_function_rejected(self):
        raw = Path(bridge.original.__file__).read_bytes()
        text = raw.decode().replace('def PP5_BP6_criteria(', 'def changed_PP5_BP6_criteria(', 1)
        raws = {bridge.original.__name__: text.encode(), bridge.interpreter.__name__: Path(bridge.interpreter.__file__).read_bytes()}
        with self.assertRaisesRegex(ValueError, 'source changed'):
            source.compile_effects(raws, {bridge.original.__name__: bridge.original.__file__, bridge.interpreter.__name__: bridge.interpreter.__file__})

    def test_duplicate_function_rejected(self):
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            source._definitions(b'def f(): pass\ndef f(): pass\n')

    def test_async_function_rejected(self):
        with self.assertRaisesRegex(ValueError, 'async'):
            source._definitions(b'async def f(): pass\n')

    def test_decorated_function_rejected(self):
        with self.assertRaisesRegex(ValueError, 'decorator'):
            source._definitions(b'@decorator\ndef f(): pass\n')

    def test_mature_effect_not_in_original_source_rejected(self):
        raws = {module.__name__: Path(module.__file__).read_bytes() for module in (bridge.original, bridge.interpreter)}
        filenames = {module.__name__: module.__file__ for module in (bridge.original, bridge.interpreter)}
        with mock.patch.dict(source.EFFECTS, {'row.missing': (source.PUBLIC[0], "float(row.get('MaveDB_urn'))")}):
            with self.assertRaisesRegex(ValueError, 'absent from original'):
                source.compile_effects(raws, filenames)

    def test_repeated_getters_and_fillna_use_their_actual_original_site(self):
        raws = {module.__name__: Path(module.__file__).read_bytes() for module in (bridge.original, bridge.interpreter)}
        filenames = {module.__name__: module.__file__ for module in (bridge.original, bridge.interpreter)}
        expressions = source.compile_effects(raws, filenames)
        for key, line in source.EFFECT_LINES.items():
            with self.subTest(key=key):
                self.assertIn(line, [row[2] for row in expressions[key].co_lines()])
        self.assertNotEqual(expressions['criteria.text.Pathogenic'].co_linetable, expressions['criteria.text.Benign'].co_linetable)

    def test_missing_repeated_expression_site_is_rejected(self):
        raws = {module.__name__: Path(module.__file__).read_bytes() for module in (bridge.original, bridge.interpreter)}
        filenames = {module.__name__: module.__file__ for module in (bridge.original, bridge.interpreter)}
        with mock.patch.dict(source.EFFECT_LINES, {'row.metadata-column': 999}):
            with self.assertRaisesRegex(ValueError, 'absent from original'):
                source.compile_effects(raws, filenames)

    def test_changed_helper_rejected_before_input_getters(self):
        changed = lambda *args, **kwargs: {'call': 'Neutral'}
        with mock.patch.object(bridge.interpreter, 'interpret_score', changed):
            with self.assertRaises(ValueError):
                source.guard(bridge.original, bridge.interpreter)

    def test_changed_libm_nan_rule_is_rejected(self):
        with mock.patch.object(bridge.interpreter.math, 'isnan', lambda value: False):
            with self.assertRaisesRegex(ValueError, 'libm'):
                source.guard(bridge.original, bridge.interpreter)

    def test_exception_catch_global_is_not_silently_rebound(self):
        with mock.patch.dict(bridge.original.__dict__, {'ValueError': TypeError}):
            with self.assertRaisesRegex(ValueError, 'builtin global'):
                source.guard(bridge.original, bridge.interpreter)

    def test_changed_native_materializer_is_rejected_by_named_source_proof(self):
        with mock.patch.object(bridge.pd.Series, '_construct_result', lambda *args, **kwargs: None):
            with self.assertRaisesRegex(ValueError, 'package executable'):
                source.package_effect_guard()


class PrimitiveTransport(unittest.TestCase):
    def test_ieee_payload_preserves_sign_and_nan_bits(self):
        for bits in ('0000000000000000', '8000000000000000', '0000000000000001',
                     '7ff0000000000000', 'fff0000000000000', '7ff8000000000042'):
            self.assertEqual(bridge._bits(bridge._float(bits)), bits)

    def test_noncanonical_text_and_float_rejected(self):
        for value in ('4A', 'f', 'ff', '61 62'):
            with self.subTest(value=value), self.assertRaises(bridge.NativeTransportError):
                bridge._untext(value)
        for value in ('NaN', '0', '000000000000000G', '7FF0000000000000'):
            with self.subTest(value=value), self.assertRaises(bridge.NativeTransportError):
                bridge._float(value)

    def test_ordinals_reject_negative_and_noncanonical(self):
        self.assertEqual(bridge._positions('0,3,2'), [0, 3, 2])
        for value in ('-1', '01', '1,', '+1'):
            with self.subTest(value=value), self.assertRaises(bridge.NativeTransportError):
                bridge._positions(value)

    def test_raw_units_and_same_spec_mutation_are_preserved(self):
        spec = bridge.interpreter.AssaySpec('Binding', bridge.interpreter.BindingSubtype.KD,
            bridge.interpreter.Directionality.NON_MONOTONIC, bridge.interpreter.Scale.RAW,
            '', None, (-0.5, 0.5))
        same = spec
        native = bridge.NativeEffects(mock.Mock())
        result = native._result(['return', bridge._text('Neutral'), bridge._bits(-0.0),
            bridge._text('Binding band'), '', '1', 'units'], spec)
        self.assertIs(same, spec)
        self.assertIs(result['notes'], None)
        self.assertIs(spec.directionality, bridge.interpreter.Directionality.HIGHER_IS_FUNCTION)
        self.assertEqual(bridge._bits(result['standardized_score']), '8000000000000000')

    def test_transport_fault_is_not_swallowed_by_log2_scientific_catch(self):
        native = bridge.NativeEffects(mock.Mock())
        native.session.exchange.side_effect = bridge.NativeTransportError('transport stopped')
        with self.assertRaisesRegex(bridge.NativeTransportError, 'transport stopped'):
            native._to_log2(1.0)


class RawAdmission(unittest.TestCase):
    def test_actual_admission_accepts_supported_dataframe_and_metadata(self):
        import pandas as pd
        frame = pd.DataFrame({'MaveDB_urn': ['one'], 'chrom': ['1'], 'pos': [1], 'ref': ['A'], 'alt': ['C']}, index=[7])
        metadata = pd.DataFrame({'urn': ['one'], 'assay_family': ['Activity'], 'rationale': ['ordinary']})
        bridge._admit('mavedb_score_interpretation', (frame, metadata))

    def test_raw_inspection_keeps_cold_and_stale_cache_untouched(self):
        import pandas as pd
        frame = pd.DataFrame({'assay_family': ['Activity', 'Binding'], 'score': [1.0, 2.0]}, index=[4, 1])
        stale = frame['score']
        stale_cache = frame._item_cache
        before = dict(stale_cache)
        bridge._raw_values(frame)
        bridge._metadata_families(frame)
        self.assertIs(frame._item_cache, stale_cache)
        self.assertEqual(list(frame._item_cache), list(before))
        self.assertIs(frame._item_cache['score'], stale)

    def test_nonstring_family_rejected_before_factorization(self):
        import pandas as pd
        frame = pd.DataFrame({'assay_family': [float('nan')]})
        bridge._raw_values(frame)
        with self.assertRaisesRegex(ValueError, 'nonstring assay-family'):
            bridge._metadata_families(frame)

    def test_rangeindex_inspection_does_not_populate_data_cache(self):
        import pandas as pd
        index = pd.RangeIndex(5, 15, 2)
        before = copy.copy(index.__dict__)
        self.assertEqual(list(bridge._axis_values(index)), [5, 7, 9, 11, 13])
        self.assertEqual(index.__dict__, before)

    def test_modified_bypassed_mask_apis_reject_before_frame_access(self):
        from pandas.core.strings.accessor import StringMethods
        from pandas.core.series import Series
        import pandas as pd
        frame = pd.DataFrame({'CLNSIG': ['Pathogenic'], 'CLNREVSTAT': ['practice_guideline']})
        for owner, name in ((StringMethods, 'contains'), (Series, '__ge__'),
                            (Series, '__and__'), (Series, '__or__'), (Series, '__invert__')):
            with self.subTest(name=name), mock.patch.object(owner, name, lambda *args, **kwargs: None):
                with mock.patch.object(bridge, '_raw_values', side_effect=AssertionError('frame accessed before API proof')):
                    with self.assertRaises(ValueError):
                        bridge._admit('PS3_BS3_criteria', (frame, ''))

    def test_text_admission_observes_raw_and_stale_cache_without_repairing_either(self):
        import pandas as pd
        frame = pd.DataFrame({'CLNSIG': ['Benign', 'Pathogenic']}, index=[4, 1])
        cached = pd.Series(['Pathogenic', 'Benign'], name='CLNSIG', index=frame.index)
        frame._item_cache['CLNSIG'] = cached
        manager, cache = frame._mgr, frame._item_cache
        bridge._text_column(frame, 'CLNSIG')
        self.assertIs(frame._mgr, manager)
        self.assertIs(frame._item_cache, cache)
        self.assertIs(frame._item_cache['CLNSIG'], cached)
        self.assertEqual(frame._mgr.iget_values(0).tolist(), ['Benign', 'Pathogenic'])
        cached.iloc[0] = 3
        with self.assertRaisesRegex(ValueError, 'text comparison'):
            bridge._text_column(frame, 'CLNSIG')

    def test_unsupported_exclusion_mask_is_rejected_before_any_source_getter(self):
        import pandas as pd
        import numpy as np
        frame = pd.DataFrame({'CLNSIG': ['Pathogenic']})
        with self.assertRaisesRegex(ValueError, 'exclusion mask'):
            bridge._table_admission('PP5_BP6_criteria', (frame, np.array([1]), np.array([False])))

    def test_lone_surrogate_text_rejects_before_native_encoding_or_effects(self):
        import pandas as pd
        row = pd.Series({'MaveDB_urn': '\ud800'}, dtype=object)
        with self.assertRaisesRegex(ValueError, 'non-UTF-8'):
            bridge._raw_values(row)

    def test_signaling_rank_nan_rejects_by_bits_without_numpy_error_callbacks(self):
        import numpy as np
        callbacks = []
        before, handler = np.geterr(), np.geterrcall()
        try:
            np.seterrcall(lambda *args: callbacks.append(args))
            np.seterr(all='call')
            self.assertFalse(bridge._rank_admitted(bridge._float('7ff0000000000001')))
            self.assertTrue(bridge._rank_admitted(bridge._float('7ff8000000000001')))
            self.assertTrue(bridge._rank_admitted(float('inf')))
            self.assertEqual(callbacks, [])
        finally:
            np.seterr(**before)
            np.seterrcall(handler)


class EffectPublications(unittest.TestCase):
    def setUp(self):
        self.session = mock.Mock()
        self.host = bridge.HostEffects(self.session, {})

    def row_context(self):
        import pandas as pd
        return {'row': pd.Series({'retained': 'x'}), '_row_publications': [], '_missing_result': '1'}

    def test_reordered_or_duplicate_row_write_rejected_before_mutation(self):
        c = self.row_context()
        with self.assertRaisesRegex(bridge.NativeTransportError, 'reordered'):
            self.host.effect('row.publish', ['MaveDB_BS3', 'false'], c)
        self.assertEqual(list(c['row'].index), ['retained'])
        self.host.effect('row.publish', ['MaveDB_PS3', 'false'], c)
        with self.assertRaisesRegex(bridge.NativeTransportError, 'reordered'):
            self.host.effect('row.publish', ['MaveDB_PS3', 'false'], c)
        self.assertIs(c['row']['MaveDB_PS3'], False)

    def test_failed_row_write_is_not_recorded_as_success(self):
        import pandas as pd
        c = self.row_context()
        c['row'] = pd.Series({'MaveDB_PS3': 'existing'}, dtype='string')
        with self.assertRaises(TypeError):
            self.host.effect('row.publish', ['MaveDB_PS3', 'false'], c)
        self.assertEqual(c['_row_publications'], [])
        self.assertEqual(c['row'].iloc[0], 'existing')

    def test_failed_numpy_index_write_retains_genuine_exception_and_prefix(self):
        import numpy as np
        c = {'_program': 'PS3_BS3_criteria', '_array_publications': [],
             'ps3_array': np.zeros(2, dtype=int), 'ps3_criteria': np.array([0, 3])}
        with self.assertRaises(IndexError):
            self.host.effect('criteria.assign', ['ps3_array', 'ps3_criteria', '3'], c)
        self.assertEqual(c['_array_publications'], [])
        self.assertEqual(c['ps3_array'].tolist(), [0, 0])

    def test_unsupported_native_effect_never_runs_arbitrary_python(self):
        c = self.row_context()
        with self.assertRaisesRegex(bridge.NativeTransportError, 'unknown'):
            self.host.effect('eval', ['row.clear()'], c)
        self.assertEqual(c['row'].to_dict(), {'retained': 'x'})

    def test_pvalue_float_catches_only_valueerror_and_keeps_nan(self):
        import math
        code = compile('float(pvalue)', '<effect-test>', 'eval')
        host = bridge.HostEffects(self.session, {'row.float-pvalue': code})
        c = {'pvalue': 'malformed'}
        self.assertEqual(host.effect('row.float-pvalue', [], c), ['ValueError'])
        self.assertTrue(math.isnan(c['pvalue_float']))
        with self.assertRaises(TypeError):
            host.effect('row.float-pvalue', [], {'pvalue': None})

    def test_band_conversion_plan_keeps_original_catch_boundaries(self):
        c = {}
        self.assertEqual(self.host.effect('band.float-tokens', ['float', bridge._text('bad')], c), ['parse-error'])
        with self.assertRaises(bridge.NativeTransportError):
            self.host.effect('band.float-tokens', ['unknown', bridge._text('1')], c)

    def test_same_spec_mutation_occurs_only_after_successful_log(self):
        spec = bridge.interpreter.AssaySpec('Binding', bridge.interpreter.BindingSubtype.KD,
            bridge.interpreter.Directionality.NON_MONOTONIC, bridge.interpreter.Scale.RAW,
            '', 'Kd', (-0.5, 0.5))
        native = bridge.NativeEffects(mock.Mock())
        native.session.exchange.side_effect = [ ['standardize'], ['none'], ['0'], ['raw', 'binding-kd'], ['binding-log10'] ]
        with self.assertRaises(ValueError):
            native.interpret_score(float('inf'), spec)
        self.assertIs(spec.directionality, bridge.interpreter.Directionality.NON_MONOTONIC)

    def test_rationale_conversion_cannot_move_before_both_criteria_writes(self):
        c = {'_row_publications': ['MaveDB_PS3'], 'score_interpretation': [object()]}
        with self.assertRaisesRegex(bridge.NativeTransportError, 'preceded criteria'):
            self.host.effect('row.rationale-texts', [], c)

    def test_metadata_writes_reject_reordering_and_duplicate_final_flag(self):
        c = {'_metadata_publications': [], '_metadata_final_publications': [], 'result': {'call': 'LoF'}}
        with self.assertRaisesRegex(bridge.NativeTransportError, 'reordered'):
            self.host.effect('pick.publish', ['urn', 'source'], c)
        c['_metadata_publications'] = ['metadata_flags', 'urn', 'n_syn_used', 'band_source', 'assay_family',
            'subtype', 'neutral_band_parsed', 'MaveDB_BS3', 'MaveDB_PS3']
        self.host.effect('pick.publish', ['MaveDB_PS3', 'true'], c)
        with self.assertRaisesRegex(bridge.NativeTransportError, 'duplicated'):
            self.host.effect('pick.publish', ['MaveDB_PS3', 'true'], c)

    def test_result_dictionary_cannot_publish_before_both_strength_writes(self):
        c = {'_array_publications': ['ps3_array'], '_result_publications': ['clinvar_patho', 'clinvar_benign'],
            'result_dict': {}, 'ps3_array': object()}
        with self.assertRaisesRegex(bridge.NativeTransportError, 'preceded strength'):
            self.host.effect('criteria.result-array', ['PS3', 'ps3_array'], c)
        self.assertEqual(c['result_dict'], {})

    def test_native_evidence_ordinals_cannot_duplicate_actual_objects(self):
        with self.assertRaisesRegex(bridge.NativeTransportError, 'duplicated'):
            self.host._selected([object()], '0,0')


class NativeMaskMaterialization(unittest.TestCase):
    def setUp(self):
        raw = {module.__name__: Path(module.__file__).read_bytes() for module in (bridge.original, bridge.interpreter)}
        names = {module.__name__: module.__file__ for module in (bridge.original, bridge.interpreter)}
        self.session = mock.Mock()
        self.session.host_compatibility_counts = {}
        self.host = bridge.HostEffects(self.session, source.compile_effects(raw, names))

    def test_native_contains_materializes_actual_cached_values_index_and_metadata(self):
        import pandas as pd
        frame = pd.DataFrame({'CLNSIG': ['Benign', 'Benign']}, index=[4, 1])
        cached = pd.Series(['Pathogenic', 'Benign'], index=frame.index, name='CLNSIG')
        cached.attrs['retained'] = ['metadata']
        frame._item_cache['CLNSIG'] = cached
        self.session.exchange.return_value = ['1', '0']
        result = self.host._contains('criteria.text.Pathogenic', 'Pathogenic', {'df': frame})
        self.session.exchange.assert_called_once_with('contains', 'Pathogenic', bridge._text('Pathogenic'), bridge._text('Benign'))
        self.assertEqual(result.tolist(), [True, False])
        self.assertIs(result.index, cached.index)
        self.assertEqual(result.attrs, cached.attrs)
        self.assertIs(frame._item_cache['CLNSIG'], cached)
        self.assertEqual(result.name, 'CLNSIG')

    def test_nullable_native_mask_keeps_boolean_extension_dtype_and_clear_mask(self):
        import pandas as pd
        frame = pd.DataFrame({'CLNSIG': pd.Series(['Pathogenic', 'Benign'], dtype='string')})
        self.session.exchange.return_value = ['1', '0']
        result = self.host._contains('criteria.text.Pathogenic', 'Pathogenic', {'df': frame})
        self.assertEqual(str(result.dtype), 'boolean')
        self.assertEqual(result.array._mask.tolist(), [False, False])
        self.assertEqual(result.array._data.tolist(), [True, False])

    def test_native_urn_mask_selects_first_duplicate_without_replacing_row_values(self):
        import pandas as pd
        metadata = pd.DataFrame({'urn': ['one', 'two', 'one'], 'rationale': ['first', 'other', 'last']}, index=[4, 1, 4])
        self.session.exchange.return_value = ['0', '0', '1', '0', '1']
        c = {'metadata_df': metadata, 'urn': 'one'}
        self.host.effect('row.metadata', [], c)
        self.assertEqual(c['meta_row']['rationale'], 'first')
        self.assertEqual(c['meta_row'].name, 4)
        self.assertEqual(metadata['rationale'].tolist(), ['first', 'other', 'last'])

    def test_native_first_duplicate_and_mask_disagreement_is_rejected(self):
        import numpy as np
        self.session.exchange.return_value = ['0', '1', '1', '0']
        with self.assertRaisesRegex(bridge.NativeTransportError, 'first duplicate'):
            self.host._urn_mask(np.array(['one', 'two'], dtype=object), 'one')

    def test_native_logical_values_retain_nullable_representation_and_result_name(self):
        import pandas as pd
        left = pd.Series([True, False], dtype='boolean', name='CLNSIG')
        right = pd.Series([False, True], index=left.index, name='CLNREVSTAT')
        left.attrs['retained'] = 'left'
        self.session.exchange.return_value = ['1', '1']
        result = self.host._logical(left, right, 'or')
        self.assertEqual(result.tolist(), [True, True])
        self.assertEqual(str(result.dtype), 'boolean')
        self.assertEqual(result.attrs, left.attrs)
        self.assertIs(result.index, left.index)
        self.assertIsNone(result.name)
        self.assertEqual(self.session.host_compatibility_counts, {})


if __name__ == '__main__':
    unittest.main()
