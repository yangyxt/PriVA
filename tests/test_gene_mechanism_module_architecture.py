from __future__ import annotations

import ast
from graphlib import TopologicalSorter
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import gene_mechanism_annotation as annotation  # noqa: E402
import gene_mechanism_common as common  # noqa: E402
import gene_mechanism_conditions as conditions  # noqa: E402
import gene_mechanism_hub as hub  # noqa: E402
import gene_mechanism_resources as resources  # noqa: E402
import gene_mechanism_variants as variants  # noqa: E402


DEFINITION_MODULES = (
    "gene_mechanism_common",
    "gene_mechanism_conditions",
    "gene_mechanism_variants",
    "gene_mechanism_resources",
    "gene_mechanism_annotation",
)


def _internal_imports(module_name: str) -> set[str]:
    source = (SCRIPTS / f"{module_name}.py").read_text(encoding="utf-8")
    tree = ast.parse(source)
    dependencies: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            dependencies.update(
                alias.name for alias in node.names if alias.name in DEFINITION_MODULES
            )
        elif (
            isinstance(node, ast.ImportFrom)
            and node.module in DEFINITION_MODULES
        ):
            dependencies.add(node.module)
    return dependencies


def test_definition_modules_form_a_one_way_acyclic_dependency_graph() -> None:
    dependencies = {
        module_name: _internal_imports(module_name)
        for module_name in DEFINITION_MODULES
    }
    order = tuple(TopologicalSorter(dependencies).static_order())

    assert set(order) == set(DEFINITION_MODULES)
    assert dependencies == {
        "gene_mechanism_common": set(),
        "gene_mechanism_conditions": {"gene_mechanism_common"},
        "gene_mechanism_variants": {"gene_mechanism_common"},
        "gene_mechanism_resources": {
            "gene_mechanism_common",
            "gene_mechanism_conditions",
            "gene_mechanism_variants",
        },
        "gene_mechanism_annotation": {
            "gene_mechanism_common",
            "gene_mechanism_conditions",
            "gene_mechanism_resources",
            "gene_mechanism_variants",
        },
    }


def test_only_hub_is_an_executable_entry_point() -> None:
    for module_name in DEFINITION_MODULES:
        source = (SCRIPTS / f"{module_name}.py").read_text(encoding="utf-8")
        tree = ast.parse(source)
        top_level_functions = {
            node.name for node in tree.body if isinstance(node, ast.FunctionDef)
        }

        assert "gene_mechanism_hub" not in source
        assert "main" not in top_level_functions
        assert not any(isinstance(node, ast.If) for node in tree.body)

    hub_tree = ast.parse(
        (SCRIPTS / "gene_mechanism_hub.py").read_text(encoding="utf-8")
    )
    assert any(
        isinstance(node, ast.FunctionDef) and node.name == "main"
        for node in hub_tree.body
    )
    assert any(isinstance(node, ast.If) for node in hub_tree.body)


def test_hub_reexports_the_existing_pipeline_contract() -> None:
    assert hub.GeneMechanismHub is resources.GeneMechanismHub
    assert (
        hub.annotate_gene_mechanism_categories
        is annotation.annotate_gene_mechanism_categories
    )
    assert hub.infer_query_variant_effect is variants.infer_query_variant_effect
    assert hub.normalize_inheritance is conditions.normalize_inheritance
    assert (
        hub.condition_cache_mechanism_assertions
        is conditions.condition_cache_mechanism_assertions
    )
    assert hub.VARIANT_MECHANISM_OUTPUT_COLUMNS is common.VARIANT_MECHANISM_OUTPUT_COLUMNS
    assert len(hub.VARIANT_MECHANISM_OUTPUT_COLUMNS) == 15
