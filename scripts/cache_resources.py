#!/usr/bin/env python3
"""Inspect PriVA's resource lifecycle map without downloading or refreshing data.

python scripts/cache_resources.py list --refreshable
python scripts/cache_resources.py show --resource clinvar_gene_statistics
python scripts/cache_resources.py refresh-plan --resource alphamissense_prescores
python scripts/cache_resources.py validate --check-routes

The refresh plan lists dependent resources for review; it never runs a command
or changes scientific policy. URL templates are not availability guarantees.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys


DEFAULT = Path(__file__).resolve().parents[1] / "cache_resources.json"


def load_registry(path: Path) -> dict:
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate registry key: {key}")
            result[key] = value
        return result

    registry = json.loads(path.read_text(), object_pairs_hook=unique)
    if registry.get("schema_version") != 1:
        raise ValueError("Unsupported cache resource schema_version")
    resources = registry.get("resources")
    if not isinstance(resources, dict) or not resources:
        raise ValueError("The resource map must be a nonempty object")
    for name, entry in resources.items():
        if not re.fullmatch(r"[a-z][a-z0-9_]*", name):
            raise ValueError(f"Invalid resource ID: {name}")
        if entry.get("lifecycle") not in {"maintained", "snapshot", "derived"}:
            raise ValueError(f"Invalid lifecycle for {name}")
        for key in ("config_keys", "bundled_paths", "depends_on"):
            values = entry.get(key)
            if not isinstance(values, list) or not all(isinstance(v, str) for v in values):
                raise ValueError(f"{name}.{key} must be a list of strings")
        assembly_dependencies = entry.get("assembly_dependencies", {})
        if not isinstance(assembly_dependencies, dict) or not all(
                isinstance(values, list) and all(isinstance(v, str) for v in values)
                for values in assembly_dependencies.values()):
            raise ValueError(f"Invalid assembly dependencies for {name}")
        refresh = entry.get("refresh", {})
        for key in ("automatic", "keep_function"):
            if type(refresh.get(key)) is not bool:
                raise ValueError(f"{name}.refresh.{key} must be boolean")
        functions = refresh.get("entrypoints")
        if not isinstance(functions, list) or not all(isinstance(f, str) for f in functions):
            raise ValueError(f"Invalid refresh entrypoints for {name}")
        if refresh["keep_function"] != bool(functions):
            raise ValueError(f"Refresh-function inventory is inconsistent for {name}")
        for dependency in dependencies(entry):
            if dependency not in resources:
                raise ValueError(f"Unknown dependency {dependency} of {name}")
    visiting, visited = set(), set()

    def visit(name):
        if name in visiting:
            raise ValueError(f"Resource dependency cycle at {name}")
        if name in visited:
            return
        visiting.add(name)
        for parent in dependencies(resources[name]):
            visit(parent)
        visiting.remove(name)
        visited.add(name)

    for name in resources:
        visit(name)
    return registry


def dependencies(entry: dict) -> set[str]:
    result = set(entry["depends_on"])
    for values in entry.get("assembly_dependencies", {}).values():
        result.update(values)
    return result


def check_routes(registry: dict, source_root: Path) -> None:
    shell = (source_root / "scripts/install_utils.sh").read_text()
    functions = set(re.findall(r"^function\s+([A-Za-z_][A-Za-z0-9_]*)\s*\(", shell, re.M))
    for name, entry in registry["resources"].items():
        for function in set(entry["refresh"]["entrypoints"] + entry.get("install_entrypoints", [])):
            if function not in functions:
                raise ValueError(f"Missing installer entrypoint {function} for {name}")


def refresh_plan(resources: dict, start: str) -> dict:
    affected = {start}
    ordered = [start]
    while True:
        following = [name for name, entry in resources.items()
                     if name not in affected and dependencies(entry) & affected]
        if not following:
            break
        affected.update(following)
        ordered.extend(following)
    return {
        "mode": "read_only_inventory",
        "requested_resource": start,
        "note": "Review affected resources for the selected assembly/profile. This list is not an execution order or authorization to refresh.",
        "affected": [{"resource": name, "lifecycle": resources[name]["lifecycle"],
                      "refresh": resources[name]["refresh"]} for name in ordered],
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("list", "show", "refresh-plan", "validate"):
        command = commands.add_parser(name)
        command.add_argument("--registry", type=Path, default=DEFAULT)
        if name in ("show", "refresh-plan"):
            command.add_argument("--resource", required=True)
        if name == "list":
            command.add_argument("--refreshable", action="store_true")
        if name == "validate":
            command.add_argument("--check-routes", action="store_true")
    args = parser.parse_args(argv)
    try:
        registry = load_registry(args.registry)
        resources = registry["resources"]
        if args.command == "validate":
            if args.check_routes:
                check_routes(registry, args.registry.resolve().parent)
            print(json.dumps({"valid": True, "resources": len(resources),
                              "cache_payloads_read": False}))
        elif args.command == "list":
            print("resource\tlifecycle\tautomatic\trefresh_entrypoints")
            for name, entry in resources.items():
                refresh = entry["refresh"]
                if args.refreshable and not refresh["keep_function"]:
                    continue
                print(f"{name}\t{entry['lifecycle']}\t{str(refresh['automatic']).lower()}\t{','.join(refresh['entrypoints'])}")
        else:
            if args.resource not in resources:
                raise ValueError(f"Unknown resource: {args.resource}")
            result = (resources[args.resource] if args.command == "show"
                      else refresh_plan(resources, args.resource))
            print(json.dumps(result, indent=2))
    except (OSError, ValueError, TypeError, KeyError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
