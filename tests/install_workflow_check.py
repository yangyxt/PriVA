#!/usr/bin/env python3
"""Regenerates workflows/install.nf from conf/install_steps.json.

Default: exit 1 when the committed workflow differs from the registry, or
when priva names an installer function that no step's code digest reaches.
--write: rewrite the workflow file.
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HEADER = ROOT / "tests" / "install_workflow_header.nf"


def render() -> str:
    reg = json.loads((ROOT / "conf" / "install_steps.json").read_text())
    steps = reg["steps"]
    ids = {s["id"]: s for s in steps}
    alias = lambda s: f"STEP_{s['id']}_{s['name']}"
    for s in steps:
        for n in s["needs"]:
            if n not in ids:
                raise SystemExit(f"{s['id']} needs unknown step {n}")
            if [x["id"] for x in steps].index(n) > [x["id"] for x in steps].index(s["id"]):
                raise SystemExit(f"{s['id']} is listed before its input {n}")
    include = "; ".join(f"INSTALL_STEP as {alias(s)}" for s in steps)
    lines = [HEADER.read_text().replace("@INCLUDES@", include).rstrip("\n")]
    for s in steps:
        needs = ", ".join(f"{alias(ids[n])}.out.receipt" for n in s["needs"])
        lines.append(f"    {alias(s)}(installWire('{s['id']}', plan, [{needs}], nativeSteps))")
    receipts = ", ".join(f"{alias(s)}.out.receipt" for s in steps)
    lines.append(FOOTER.replace("@RECEIPTS@", receipts))
    return "\n".join(lines)


FOOTER = """    def all = Channel.empty().mix(@RECEIPTS@).collect()
    def toPublish = plan.combine(all.map { [it] }).map { p, receipts ->
        def names = receipts.collect { it.toString() }.sort()
        def lineages = names.collect { r -> new groovy.json.JsonSlurper().parse(new File(r)).lineage }
        [installSha(groovy.json.JsonOutput.toJson(lineages)).take(16), p, names]
    }
    INSTALL_PUBLISH(toPublish)
    emit:
    published = INSTALL_PUBLISH.out.published
}
"""

SHELL_LIBRARIES = ("install_utils.sh", "common_bash_utils.sh", "annotation_vcf.sh")
# Entry points of the run workflow, which no install step runs: priva starts
# main_workflow for a run task, and only checks that anno_VEP_data exists
# (clinvar-gene-refresh candidate-route).
RUN_ENTRIES = {"main_workflow", "anno_VEP_data"}
WORD = r"[A-Za-z_][A-Za-z0-9_]*"


def reached_functions() -> set:
    """The installer functions some step's code digest reaches.

    Follows scripts/install_task.sh (step_functions, called_words) on the
    source text: from every step's calls and runs, through every function
    those name; in a log, echo or printf line only where a command can start.
    """
    bodies, current = {}, None
    for name in SHELL_LIBRARIES:
        for line in (ROOT / "scripts" / name).read_text().splitlines():
            start = re.match(r"function\s+(" + WORD + ")", line)
            if start:
                current = start.group(1)
                bodies[current] = []
            elif current and not line.lstrip().startswith("#"):
                bodies[current].append(line)
    reg = json.loads((ROOT / "conf" / "install_steps.json").read_text())
    queue = [c[0] for s in reg["steps"] for c in s["calls"]] + [r for s in reg["steps"] for r in s.get("runs", [])]
    for name in (r for s in reg["steps"] for r in s.get("runs", [])):
        if name not in bodies:
            raise SystemExit(f"a step's runs names {name}, which is not an installer function")
    seen = set()
    while queue:
        name = queue.pop()
        if name in seen or name not in bodies:
            continue
        seen.add(name)
        for line in bodies[name]:
            if re.match(r"\s*(log|echo|printf)\b", line):
                words = re.findall(r"(?:\$\(|`|\||&|;)\s*(" + WORD + ")", line)
            else:
                words = re.findall(WORD, line)
            queue.extend(w for w in words if w in bodies and w not in seen)
    return seen


def check_priva_runs() -> None:
    """priva runs some installer functions itself (bash install_utils.sh NAME).

    Each installer function named in a Rust string must be reached by a step's
    code digest; otherwise an edit to it would not run that step again.
    """
    reached = reached_functions()
    names = set()
    for line in SHELL_LIBRARIES:
        names |= set(re.findall(r"^function\s+(" + WORD + ")", (ROOT / "scripts" / line).read_text(), re.M))
    missing = {}
    for path in sorted((ROOT / "rust" / "src").rglob("*.rs")):
        for name in set(re.findall(r'"(' + WORD + ')"', path.read_text())) & names:
            if name not in reached and name not in RUN_ENTRIES:
                missing.setdefault(name, []).append(str(path.relative_to(ROOT)))
    if missing:
        for name, paths in sorted(missing.items()):
            print(f"{name} ({', '.join(paths)}) is reached by no install step; add it to the runs of the step whose program runs it, in conf/install_steps.json")
        sys.exit(1)


if __name__ == "__main__":
    target = ROOT / "workflows" / "install.nf"
    text = render()
    if "--write" in sys.argv[1:]:
        target.write_text(text)
    elif target.read_text() != text:
        print("workflows/install.nf differs from conf/install_steps.json; run with --write")
        sys.exit(1)
    check_priva_runs()
