"""Source/evaluation-order controls, independent of the saved scientific oracle.

These tests compile synthetic expressions and the admitted source, but never
invoke an original frequency criterion or calculate a reference score.
"""
import ast
import hashlib
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import priva_rust_frequency_source as source


def lowered_synthetic(text, runtime, **namespace):
    transformer = source.Lower()
    tree = transformer.visit(ast.parse(text))
    scope = {source.RUNTIME: runtime, **namespace}
    exec(compile(ast.fix_missing_locations(tree), '<frequency-source-control>', 'exec'), scope)
    return scope, transformer


class FrequencySourceTests(unittest.TestCase):
    def test_complete_original_function_inventory_and_source(self):
        raw = (ROOT / 'scripts/acmg_pm2_bs1_ba1_frequency.py').read_bytes()
        self.assertEqual(hashlib.sha256(raw).hexdigest(), source.REFERENCE_SHA)
        nodes = [node for node in ast.parse(raw).body if isinstance(node, ast.FunctionDef)]
        self.assertEqual(tuple(node.name for node in nodes), source.PROCEDURES)
        bs1 = next(node for node in nodes if node.name == 'BS1_criteria')
        self.assertEqual([node.name for node in bs1.body if isinstance(node, ast.FunctionDef)],
                         ['_series_or_empty', '_numeric_series', '_frequency_above', '_hom_bs1'])

    def test_caller_exact_ordered_twenty_three_statements(self):
        raw = (ROOT / 'scripts/acmg_criteria_assign.py').read_bytes()
        self.assertEqual(hashlib.sha256(raw).hexdigest(), source.CALLER_SHA)
        caller = next(node for node in ast.parse(raw).body if isinstance(node, ast.FunctionDef) and node.name == 'ACMG_criteria_assign')
        nodes = [node for node in caller.body if 391 <= node.lineno <= 443]
        self.assertEqual(len(nodes), 23)
        self.assertEqual(hashlib.sha256(ast.dump(ast.Module(body=nodes, type_ignores=[]), include_attributes=False).encode()).hexdigest(), source.CALLER_AST_SHA)
        calls = [node.value.func.id for node in nodes if isinstance(node, ast.Assign)
                 and isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name)]
        self.assertEqual(calls[:3], ['BA1_criteria', 'PM2_criteria', 'BS1_criteria'])

    def test_method_lookup_precedes_argument_and_runtime(self):
        events = []
        class Owner:
            @property
            def fillna(self):
                events.append('lookup')
                return lambda value: events.append(('call', value)) or value
        def argument(): events.append('argument'); return 7
        def method(line, function, operation, name, *args, **kwargs):
            events.append('runtime')
            return operation(*args, **kwargs)
        scope, _ = lowered_synthetic('def PM2_criteria(owner):\n return owner.fillna(argument())\n',
                                    SimpleNamespace(method=method), argument=argument)
        self.assertEqual(scope['PM2_criteria'](Owner()), 7)
        self.assertEqual(events, ['lookup', 'argument', 'runtime', ('call', 7)])

    def test_failed_method_lookup_does_not_evaluate_argument(self):
        events = []
        class Owner:
            @property
            def fillna(self): raise LookupError('original lookup')
        def method(*_args, **_kwargs): self.fail('runtime method must not run after failed bound lookup')
        scope, _ = lowered_synthetic('def PM2_criteria(owner):\n return owner.fillna(argument())\n',
            SimpleNamespace(method=method), argument=lambda: events.append('argument'))
        with self.assertRaisesRegex(LookupError, 'original lookup'): scope['PM2_criteria'](Owner())
        self.assertEqual(events, [])

    def test_definition_frontier_keeps_defaults_annotations_and_decorators(self):
        events = []
        def mark(name, value): events.append(name); return value
        def decorate(value):
            events.append(('decorator_factory', value))
            def apply(function): events.append('decorator_apply'); return function
            return apply
        def binary(line, function, name, left, right):
            events.append(('runtime', name)); return left + right
        text = ("@decorate(mark('decorator_left', 1) + mark('decorator_right', 2))\n"
                "def BS1_criteria(value: mark('annotation_left', int) | mark('annotation_right', str) = mark('default', 3)) -> mark('return', int):\n"
                " return value + 1\n")
        scope, transformer = lowered_synthetic(text, SimpleNamespace(binary=binary), mark=mark, decorate=decorate)
        self.assertEqual(events, ['decorator_left', 'decorator_right', ('decorator_factory', 3),
                                 'default', 'annotation_left', 'annotation_right', 'return', 'decorator_apply'])
        self.assertEqual(len(transformer.sites), 1)
        self.assertEqual(scope['BS1_criteria'](), 4)
        self.assertEqual(events[-1], ('runtime', 'add'))

    def test_argument_callback_does_not_replace_already_bound_method(self):
        events = []
        class Owner:
            def fillna(self, value): events.append('captured_original'); return value
        def argument():
            Owner.fillna = lambda self, value: events.append('replacement')
            return 7
        def method(line, function, operation, name, *args, **kwargs):
            return operation(*args, **kwargs)
        scope, _ = lowered_synthetic('def PM2_criteria(owner):\n return owner.fillna(argument())\n',
                                    SimpleNamespace(method=method), argument=argument)
        self.assertEqual(scope['PM2_criteria'](Owner()), 7)
        self.assertEqual(events, ['captured_original'])

    def test_assignment_rhs_precedes_target_mask(self):
        events = []
        def value(): events.append('value'); return 9
        def mask(): events.append('mask'); return 0
        def assign(line, function, rhs, array, index):
            events.append('publication'); array[index] = rhs
        scope, _ = lowered_synthetic('def PM2_criteria(pm2_array):\n pm2_array[mask()] = value()\n return pm2_array\n',
            SimpleNamespace(assign=assign), value=value, mask=mask)
        owned = [0]
        self.assertIs(scope['PM2_criteria'](owned), owned)
        self.assertEqual(owned, [9])
        self.assertEqual(events, ['value', 'mask', 'publication'])

    def test_failed_assignment_rhs_leaves_target_unobserved(self):
        events = []
        def fail(): raise ValueError('original RHS')
        def assign(*_args): self.fail('runtime assignment must not run after failed RHS')
        scope, _ = lowered_synthetic('def PM2_criteria(pm2_array):\n pm2_array[mask()] = fail()\n',
            SimpleNamespace(assign=assign), fail=fail, mask=lambda: events.append('mask'))
        with self.assertRaisesRegex(ValueError, 'original RHS'): scope['PM2_criteria']([0])
        self.assertEqual(events, [])

    def test_original_ast_pin_precedes_child_lowering(self):
        text = 'def BA1_criteria(x):\n return (x + 1) > 0\n'
        original = ast.parse(text).body[0].body[0].value
        expected = hashlib.sha256(ast.dump(original, include_attributes=False).encode()).hexdigest()
        transformer = source.Lower(); transformer.visit(ast.parse(text))
        self.assertEqual(transformer.sites[-1]['original_ast_sha256'], expected)
        self.assertEqual([row['operation'] for row in transformer.sites], ['binary', 'binary'])

    def test_nested_helper_has_its_own_source_site(self):
        transformer = source.Lower()
        transformer.visit(ast.parse('def BS1_criteria(x):\n def _hom_bs1(y):\n  return y > 0\n return x > 1\n'))
        self.assertEqual([row['function'] for row in transformer.sites], ['_hom_bs1', 'BS1_criteria'])

    def test_mature_binomial_positional_and_keyword_order_retained(self):
        text = 'def control_false_neg_rate(a, b, c):\n return binom.sf(k=a, n=b, p=c)\n'
        transformer = source.Lower(); tree = transformer.visit(ast.parse(text))
        call = tree.body[0].body[0].value
        self.assertEqual(call.func.attr, 'kernel')
        self.assertEqual([item.arg for item in call.keywords], ['k', 'n', 'p'])
        self.assertEqual(ast.unparse(call.args[3]), 'binom.sf')

    def test_scalar_short_circuit_boolean_control_is_not_vectorized(self):
        transformer = source.Lower()
        tree = transformer.visit(ast.parse('def BS1_criteria(x, y):\n return x and y()\n'))
        self.assertIsInstance(tree.body[0].body[0].value, ast.BoolOp)
        self.assertEqual(transformer.sites, [])

    def test_original_codecs_and_partial_mutation_order_are_retained(self):
        raw = (ROOT / 'scripts/acmg_pm2_bs1_ba1_frequency.py').read_bytes()
        pm2 = next(node for node in ast.parse(raw).body if isinstance(node, ast.FunctionDef) and node.name == 'PM2_criteria')
        transformer = source.Lower(); transformed = transformer.visit(pm2)
        dumps = [ast.dump(node, include_attributes=False) for node in transformed.body]
        self.assertIn("attr='load'", dumps[1])
        self.assertIn("attr='open'", dumps[1])
        self.assertIn("value='Gene'", dumps[4])
        self.assertTrue(any("value='exon_patho_median_af'" in item for item in dumps))

    def test_scientific_source_is_never_written_by_lowering(self):
        text = (ROOT / 'scripts/priva_rust_frequency_source.py').read_text()
        tree = ast.parse(text)
        names = {node.func.attr for node in ast.walk(tree) if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)}
        self.assertFalse(names & {'write_text', 'write_bytes', 'unlink', 'rename', 'replace'})


if __name__ == '__main__': unittest.main()
