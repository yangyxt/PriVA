"""Strict historical summary oracle with one explicitly pinned optional seam."""
import ast
import hashlib


SUMMARY_SHA256 = '137af5bb5a6f90b0270bc94d53c0d9f6ff62160d1aca28a596f96dcb99876a0a'
HOOK_SOURCE_SHA256 = '2797f9eddb9a1b9d67e135d7fdc8a1e144a94a40b22ecfe80440029a6ea7be2a'


def definition(tree):
    nodes = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
             and n.name == 'summarize_acmg_criteria']
    if len(nodes) != 1 or not isinstance(nodes[0], ast.FunctionDef) or nodes[0].decorator_list:
        raise AssertionError('expected one undecorated synchronous top-level summary definition')
    return nodes[0]


def historical_summary_source(text):
    """Return exact historical bytes; reject any unrecognized candidate edit.

    Only the precisely pinned new branch and one literal keyword insertion may
    be removed. This is not general AST normalization or a source-hash waiver.
    The full restored function still has its original historical byte pin.
    """
    function = definition(ast.parse(text))
    segment = ast.get_source_segment(text, function)
    keywords = function.args.kwonlyargs
    if keywords:
        if (len(keywords) != 1 or keywords[0].arg != 'summary_binary'
                or keywords[0].annotation is not None
                or len(function.args.kw_defaults) != 1
                or ast.dump(function.args.kw_defaults[0]) != 'Constant(value=None)'):
            raise AssertionError('unexpected summary keyword seam')
        hooks = [node for node in function.body if isinstance(node, ast.If)
                 and ast.dump(node.test) == ast.dump(ast.parse('summary_binary is not None', mode='eval').body)]
        if len(hooks) != 1:
            raise AssertionError('missing or duplicate summary dispatch')
        hook = hooks[0]
        position = function.body.index(hook)
        following = function.body[position + 1] if position + 1 < len(function.body) else None
        if (not isinstance(following, ast.Assign) or len(following.targets) != 1
                or ast.unparse(following.targets[0]) != "df['ACMG_criteria']"):
            raise AssertionError('summary dispatch moved from the original publication boundary')
        digest = hashlib.sha256(ast.get_source_segment(text, hook).encode()).hexdigest()
        if digest != HOOK_SOURCE_SHA256:
            raise AssertionError('summary dispatch/publication order changed')
        lines = segment.splitlines(keepends=True)
        del lines[hook.lineno - function.lineno:hook.end_lineno - function.lineno + 1]
        segment = ''.join(lines)
        if segment.count(', *, summary_binary=None') != 1:
            raise AssertionError('unexpected literal summary keyword declaration')
        segment = segment.replace(', *, summary_binary=None', '', 1)
    if hashlib.sha256(segment.encode()).hexdigest() != SUMMARY_SHA256:
        raise AssertionError('historical summary source changed')
    return segment
