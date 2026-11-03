"""Frequency source admission and finite lowering; no production registration.

The canonical original files remain unchanged. A derived private namespace
retains all nine functions, nested helpers, source locations and the exact
23-statement caller segment. This file performs no scientific call on import.
"""
import ast
import copy
import gc
import hashlib
import inspect
from pathlib import Path
import stat
from types import CodeType, FunctionType

import numpy as np
import pandas as pd

REFERENCE_SHA = 'ef2e70115ddfe7087b4d24002ce4af95c4f4b91fbf38a06e8a52275eecca4fb3'
CONDITION_SHA = 'c7e8da80ca0ec47fb38a5fd1d965641c1791189cb0d25a651253f44508d0bf26'
CALLER_SHA = 'a4b79e2eee11e64d5a5d24e09d287a703feaf28cdc828bbfe0495413b3cbc373'
CALLER_AST_SHA = 'da78a9c7275c73929ed8cf78685264fa8806cdbfb12386f11d292c776823b00a'
PROCEDURES = ('_gene_max_patho_hom_freq', 'control_false_neg_rate',
    'control_false_positive_rate', 'gnomAD_rare_AF', '_rare_af_binomial',
    '_rare_hom_binomial', 'PM2_criteria', 'BS1_criteria', 'BA1_criteria')
RUNTIME = '_priva_frequency_native_runtime'

def regular_bytes(path):
    path = Path(path)
    before = path.lstat()
    if not stat.S_ISREG(before.st_mode):
        raise ValueError('nonregular source/executable')
    raw = path.read_bytes()
    after = path.lstat()
    fields = ('st_dev', 'st_ino', 'st_mode', 'st_size', 'st_mtime_ns', 'st_ctime_ns', 'st_nlink')
    if any(getattr(before, key) != getattr(after, key) for key in fields):
        raise ValueError('source/executable changed while reading')
    return raw

def placements(placement, count):
    """Read only existing Cython references; never request lazy .as_array."""
    existing = []
    references = gc.get_referents(placement)
    for value in references:
        if value is None or value is type(placement):
            continue
        if type(value) is slice:
            if any(item is not None and type(item) is not int for item in (value.start, value.stop, value.step)):
                raise ValueError('custom placement slice')
            existing.append(list(range(*value.indices(count))))
        elif type(value) is np.ndarray and value.ndim == 1 and value.dtype.kind == 'i':
            existing.append([int(item) for item in value])
        else:
            raise ValueError('unknown placement representation')
    if not 1 <= len(existing) <= 2 or len(references) > 4 or any(item != existing[0] for item in existing):
        raise ValueError('ambiguous placement representation')
    return existing[0]

def raw_columns(frame):
    from pandas.core.internals.managers import BlockManager
    from pandas.core.internals.blocks import NumpyBlock, ExtensionBlock
    if type(frame._mgr) is not BlockManager:
        raise ValueError('custom frame manager')
    count = len(frame._mgr.axes[0])
    result = [None] * count
    for block in frame._mgr.blocks:
        if type(block) not in (NumpyBlock, ExtensionBlock):
            raise ValueError('custom frame block')
        locations = placements(block._mgr_locs, count)
        values = block.values
        if type(values) is np.ndarray:
            if values.ndim != 2 or values.shape[0] != len(locations):
                raise ValueError('custom block dimensions')
            columns = [values[i] for i in range(len(locations))]
        else:
            if len(locations) != 1:
                raise ValueError('ambiguous extension block')
            columns = [values]
        for location, column in zip(locations, columns):
            if not 0 <= location < count or result[location] is not None:
                raise ValueError('invalid/overlapping block placement')
            result[location] = column
    if any(value is None for value in result):
        raise ValueError('incomplete block placement')
    return result

def plain_array(values):
    if type(values) is np.ndarray:
        if values.ndim != 1 or values.dtype.kind not in 'ObifU':
            raise ValueError('unsupported raw array')
        if values.dtype.kind == 'f' and values.dtype != np.dtype('float64'):
            # NumPy float32 scalar promotion changes exact threshold behavior;
            # do not silently substitute f64 comparison for that profile.
            raise ValueError('unsupported floating precision')
        if values.dtype.kind == 'i' and values.dtype != np.dtype('int64'):
            raise ValueError('unsupported integer storage width')
        if values.dtype.kind == 'O':
            for value in values:
                if (type(value) not in (str, bool, int, float, np.int64, np.int32, np.float64, np.float32, np.bool_)
                        and value is not None and value is not pd.NA):
                    raise ValueError('unsupported object scalar')
        return
    if type(values) in (pd.arrays.StringArray, pd.arrays.BooleanArray, pd.arrays.IntegerArray, pd.arrays.FloatingArray):
        if type(values) is pd.arrays.FloatingArray and str(values.dtype) != 'Float64':
            raise ValueError('unsupported nullable floating precision')
        if type(values) is pd.arrays.IntegerArray and str(values.dtype) != 'Int64':
            raise ValueError('unsupported nullable integer width')
        return
    raise ValueError('unsupported extension storage')


def axis_values(axis):
    """Read stored labels without creating pandas lookup engines or caches."""
    if type(axis) is pd.RangeIndex:
        return list(vars(axis)['_range'])
    if type(axis) is pd.MultiIndex:
        levels = [axis_values(level) for level in vars(axis)['_levels']]
        codes = vars(axis)['_codes']
        if any(any(int(code) < 0 for code in column) for column in codes):
            raise ValueError('missing MultiIndex labels require mature indexing')
        return [tuple(levels[j][int(codes[j][i])] for j in range(len(levels)))
                for i in range(len(codes[0]))]
    if type(axis) is not pd.Index:
        raise ValueError('unsupported index class')
    values = vars(axis)['_data']
    if type(values) is not np.ndarray or values.ndim != 1:
        raise ValueError('unsupported index storage')
    result = []
    for value in values:
        if type(value) in (str, int, tuple): result.append(value)
        elif type(value) in (np.int64, np.int32): result.append(int(value))
        elif type(value) is np.str_: result.append(str(value))
        else: raise ValueError('unsupported label value')
    return result


def function_guard(module, names, expected):
    raw = regular_bytes(module.__file__)
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError('canonical scientific source changed')
    compiled = compile(raw, str(module.__file__), 'exec')
    functions = {code.co_name: code for code in compiled.co_consts if type(code) is CodeType}
    for name in names:
        actual = getattr(module, name)
        if (type(actual) is not FunctionType or actual.__code__ != functions[name]
                or actual.__globals__ is not vars(module) or actual.__closure__ is not None
                or actual.__kwdefaults__ is not None):
            raise ValueError('canonical function binding changed: ' + name)
        declaration = next(node for node in ast.parse(raw).body
                           if isinstance(node, ast.FunctionDef) and node.name == name)
        defaults = tuple(ast.literal_eval(node) for node in declaration.args.defaults)
        if actual.__defaults__ != (defaults or None):
            raise ValueError('scientific function defaults changed: ' + name)
    return raw


def reference_guard(reference, caller_path):
    import acmg_variant_mechanism as condition
    raw = function_guard(reference, PROCEDURES, REFERENCE_SHA)
    condition_raw = function_guard(condition, ('_variant_condition_masks',), CONDITION_SHA)
    if reference._variant_condition_masks is not condition._variant_condition_masks:
        raise ValueError('condition helper binding changed')
    expected_status = {'practice_guideline': 4, 'reviewed_by_expert_panel': 3,
        'criteria_provided,_multiple_submitters,_no_conflicts': 2}
    if type(reference.HIGH_CONFIDENCE_STATUS) is not dict or list(reference.HIGH_CONFIDENCE_STATUS.items()) != list(expected_status.items()):
        raise ValueError('ClinVar review mapping changed')
    caller_raw = regular_bytes(caller_path)
    if hashlib.sha256(caller_raw).hexdigest() != CALLER_SHA:
        raise ValueError('frozen caller source changed')
    caller = next(node for node in ast.parse(caller_raw).body
                  if isinstance(node, ast.FunctionDef) and node.name == 'ACMG_criteria_assign')
    segment = [node for node in caller.body if 391 <= node.lineno and node.end_lineno <= 443]
    # This is the same exact node-list digest used by the closed oracle.
    digest = hashlib.sha256(ast.dump(ast.Module(body=segment, type_ignores=[]),
                                    include_attributes=False).encode()).hexdigest()
    if (len(segment) != 23 or segment[0].lineno != 391 or segment[-1].end_lineno != 443
            or digest != CALLER_AST_SHA):
        raise ValueError('caller frequency segment changed')
    return raw, condition_raw, segment, condition


class Lower(ast.NodeTransformer):
    """Finite table lowering. Every source expression keeps a correspondence row."""
    def __init__(self):
        self.function = None
        self.sites = []
        self.original_nodes = {}

    def visit(self, node):
        # Save the incoming source expression BEFORE recursively lowering its
        # children. A lowered operand is not an original-source digest.
        if isinstance(node, ast.AST) and hasattr(node, 'lineno'):
            self.original_nodes.setdefault(id(node), hashlib.sha256(
                ast.dump(node, include_attributes=False).encode()).hexdigest())
        return super().visit(node)

    def visit_FunctionDef(self, node):
        old = self.function
        self.function = node.name
        try:
            # Definitions retain CPython's original evaluation frontier for
            # decorators, defaults and annotations. Only runtime body
            # expressions are candidates for the native column boundary.
            node.body = [self.visit(statement) for statement in node.body]
            return node
        finally:
            self.function = old

    def call(self, node, method, operands, keywords=()):
        self.sites.append(dict(function=self.function, line=node.lineno,
            column=node.col_offset, end_line=node.end_lineno,
            end_column=node.end_col_offset, operation=method,
            original_ast_sha256=self.original_nodes[id(node)]))
        result = ast.Call(func=ast.Attribute(value=ast.Name(id=RUNTIME, ctx=ast.Load()),
                          attr=method, ctx=ast.Load()),
            args=[ast.Constant(node.lineno), ast.Constant(self.function), *operands],
            keywords=list(keywords))
        return ast.copy_location(result, node)

    def visit_BinOp(self, node):
        node = self.generic_visit(node)
        operations = {ast.Add:'add', ast.Sub:'sub', ast.Mult:'mul', ast.Div:'div',
                      ast.BitAnd:'and', ast.BitOr:'or'}
        if type(node.op) in operations:
            return self.call(node, 'binary', [ast.Constant(operations[type(node.op)]), node.left, node.right])
        return node

    def visit_Compare(self, node):
        node = self.generic_visit(node)
        operations = {ast.Eq:'eq', ast.NotEq:'ne', ast.Lt:'lt', ast.LtE:'le', ast.Gt:'gt', ast.GtE:'ge'}
        if len(node.ops) == 1 and type(node.ops[0]) in operations:
            return self.call(node, 'binary', [ast.Constant(operations[type(node.ops[0])]), node.left, node.comparators[0]])
        return node

    def visit_UnaryOp(self, node):
        node = self.generic_visit(node)
        if type(node.op) is ast.Invert:
            return self.call(node, 'invert', [node.operand])
        return node

    def visit_DictComp(self, node):
        # One literal original None-filter. Exact source admission above
        # supplies its semantics; arbitrary comprehensions stay untouched.
        if self.function == 'PM2_criteria' and node.lineno == 353:
            return self.call(node, 'filter_cache', [ast.Name(id='clinvar_patho_af_dict', ctx=ast.Load())])
        return self.generic_visit(node)

    def visit_Call(self, node):
        node = self.generic_visit(node)
        function = node.func
        if type(function) is ast.Attribute and type(function.value) is ast.Name:
            owner = function.value.id
            if owner == 'logger':
                return self.call(node, 'after_callback', [node])
            if owner == 'pd' and function.attr == 'to_numeric':
                return self.call(node, 'numeric', [function, *node.args], node.keywords)
            if owner == 'np' and function.attr in ('where','round','minimum','maximum','isnan','logical_not','zeros','full'):
                return self.call(node, 'numpy', [ast.Constant(function.attr), function, *node.args], node.keywords)
            if owner == 'binom' and function.attr in ('sf','cdf'):
                # Call the bound mature function itself, with original operands.
                return self.call(node, 'kernel', [ast.Constant(function.attr), function, *node.args], node.keywords)
        if type(function) is ast.Attribute and function.attr in (
                'fillna','isna','isnull','notna','astype','clip','mask','replace','sum','any',
                'lower','upper','removeprefix','contains','isin','eq','ge','map','apply','groupby','max','to_dict'):
            # Evaluate the actual bound-method lookup before its arguments,
            # exactly as Python does. Custom accessors are never pre-read.
            return self.call(node, 'method', [function, ast.Constant(function.attr), *node.args], node.keywords)
        return node

    def visit_Assign(self, node):
        node = self.generic_visit(node)
        arrays = {'pm2_array','bs1_array','ba1_array','pm2_hom_rare','pval','pm2_criteria'}
        if len(node.targets) == 1 and type(node.targets[0]) is ast.Subscript:
            target = node.targets[0]
            if type(target.value) is ast.Name and target.value.id in arrays:
                return ast.copy_location(ast.Expr(value=self.call(node, 'assign',
                    [node.value, ast.Name(id=target.value.id, ctx=ast.Load()), target.slice])), node)
        return node


def lowered(reference, caller_path, runtime):
    raw, condition_raw, segment, condition = reference_guard(reference, caller_path)
    nodes = [copy.deepcopy(node) for node in ast.parse(raw).body
             if isinstance(node, ast.FunctionDef) and node.name in PROCEDURES]
    condition_node = copy.deepcopy(next(node for node in ast.parse(condition_raw).body
        if isinstance(node, ast.FunctionDef) and node.name == '_variant_condition_masks'))
    skeleton = ast.parse('def ACMG_criteria_assign(anno_df, clinvar_patho_af_stat, clinvar_patho_exon_af_stat, gnomAD_extreme_rare_threshold=0.0001, expected_incidence=0.0001):\n return locals()\n').body[0]
    skeleton.lineno = 390
    skeleton.end_lineno = 445
    for argument in skeleton.args.args:
        argument.lineno = argument.end_lineno = 390
    for default in skeleton.args.defaults:
        default.lineno = default.end_lineno = 390
    skeleton.body = copy.deepcopy(segment) + [ast.Return(value=ast.Call(func=ast.Name(id='locals',ctx=ast.Load()),args=[],keywords=[]),lineno=445,col_offset=4,end_lineno=445,end_col_offset=19)]
    result, sites = {}, []
    namespace = dict(vars(reference))
    namespace.update({RUNTIME:runtime, 'gc':gc, 're':vars(condition)['re']})
    for group, filename in ((nodes, str(reference.__file__)), ([condition_node],str(condition.__file__)), ([skeleton],str(caller_path))):
        lower = Lower()
        modified = [lower.visit(node) for node in group]
        tree = ast.fix_missing_locations(ast.Module(body=modified,type_ignores=[]))
        code = compile(tree,filename,'exec')
        for item in code.co_consts:
            if type(item) is CodeType and item.co_name in (*PROCEDURES,'_variant_condition_masks','ACMG_criteria_assign'):
                original = getattr(reference,item.co_name,None)
                defaults = (0.0001,0.0001) if item.co_name=='ACMG_criteria_assign' else getattr(original,'__defaults__',None)
                function=FunctionType(item,namespace,item.co_name,defaults)
                result[item.co_name]=function
                namespace[item.co_name]=function
        sites.extend(lower.sites)
    if set(result)!=set(PROCEDURES)|{'_variant_condition_masks','ACMG_criteria_assign'}:
        raise ValueError('complete frequency lowering required')
    # The caller uses its own original logger, not the frequency logger.
    caller_namespace={name: namespace[name] for name in ('BA1_criteria','PM2_criteria','BS1_criteria','gc',RUNTIME)}
    caller_namespace['__name__']='acmg_criteria_assign'
    import logging
    caller_namespace['logger']=logging.getLogger('acmg_criteria_assign')
    caller=result['ACMG_criteria_assign']
    result['ACMG_criteria_assign']=FunctionType(caller.__code__,caller_namespace,caller.__name__,caller.__defaults__)
    return result,sites
