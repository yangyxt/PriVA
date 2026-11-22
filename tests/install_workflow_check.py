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


def process(step: dict) -> str:
    """The process a module step runs: its own name, or the process several
    steps share (registry "process")."""
    return step.get("process", step["name"])


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
    for s in steps:
        if "module" in s:
            module = ROOT / s["module"]
            if s.get("calls"):
                raise SystemExit(f"{s['id']} has a module and installer calls; a step is one or the other")
            if not s["module"].startswith("modules/local/install_steps/") or not module.is_file():
                raise SystemExit(f"{s['id']}: its module {s['module']} is not a file of modules/local/install_steps/")
            if not re.search(rf"^process {process(s)} \{{", module.read_text(), re.M):
                raise SystemExit(f"{s['id']}: {s['module']} defines no process {process(s)}")
        elif s.get("stage") or s.get("process"):
            raise SystemExit(f"{s['id']}: only a step with a module stages files or names a process")
    # The steps that still run installer functions (install-task run), if any.
    shell = "; ".join(f"INSTALL_STEP as {alias(s)}" for s in steps if "module" not in s)
    include = f"include {{ {shell} }} from '../modules/local/install'\n" if shell else ""
    modules = "\n".join(
        f"include {{ {process(s)} as {alias(s)} }} from '../{s['module'][:-len('.nf')]}'"
        for s in steps if "module" in s)
    lines = [HEADER.read_text().replace("@INCLUDES@", include).replace("@PROCESS_INCLUDES@", modules).rstrip("\n")]
    for s in steps:
        needs = ", ".join(f"{alias(ids[n])}.out.receipt" for n in s["needs"])
        if "module" in s:
            stage = ", ".join(f"'{k}'" for k in s.get("stage", []))
            lines.append(f"    {alias(s)}(installNative('{s['id']}', plan, [{needs}], [{stage}]))")
        else:
            lines.append(f"    {alias(s)}(installWire('{s['id']}', plan, [{needs}]))")
    receipts = ", ".join(f"{alias(s)}.out.receipt" for s in steps)
    lines.append(FOOTER.replace("@RECEIPTS@", receipts))
    return "\n".join(lines)


FOOTER = """    def all = Channel.empty().mix(@RECEIPTS@).collect()
    def toPublish = plan.combine(all.map { [it] }).map { p, receipts ->
        def names = receipts.collect { it.toString() }.sort()
        // Publishing writes every step's paths into one configuration, so it
        // keys on their lineages (folders), and on the publishing code.
        def lineages = names.collect { r -> new groovy.json.JsonSlurper().parse(new File(r)).lineage }
        def digest = new groovy.json.JsonSlurper().parse(p.toFile()).publish_digest
        [installSha(groovy.json.JsonOutput.toJson([lineages: lineages, publish: digest])).take(16), p, names]
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
# Rust code that only install_utils.sh runs: alpha-annotation-base, for
# AlphaMissense_anno and the motif function, which no step calls since C3 runs
# as its own process, and clinvar-gene-refresh install, deploy and the
# candidate receipts, for ClinVar_install, which no step calls since D1 runs as
# its own process. It names installer functions to fingerprint their text or
# runs them from the shell, and goes with install_utils.sh in phase 7. So does
# the shell part of install-task's digest (rust/src/install_task/digest.rs),
# which names log, echo and printf to read the text of a step's installer
# functions: since CADD (E1) runs as its own process, no step has any.
SHELL_ONLY = {"rust/src/alpha_annotation_base/profile.rs", "rust/src/alpha_annotation_base/publication.rs",
              "rust/src/clinvar_gene_refresh/deploy.rs", "rust/src/clinvar_gene_refresh/profile.rs",
              "rust/src/clinvar_gene_refresh/refresh.rs", "rust/src/install_task/digest.rs"}
WORD = r"[A-Za-z_][A-Za-z0-9_]*"


def reached_functions() -> set:
    """The installer functions some step's code digest reaches.

    Follows priva install-task (rust/src/install_task/digest.rs:
    step_functions, called_words) on the
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
    queue = [c[0] for s in reg["steps"] for c in s.get("calls", [])] + [r for s in reg["steps"] for r in s.get("runs", [])]
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
    code digest; otherwise an edit to it would not run that step again. A
    file's unit tests (from "#[cfg(test)]" followed by "mod" to its end) are
    not read: priva runs nothing that only a test names.
    """
    reached = reached_functions()
    names = set()
    for line in SHELL_LIBRARIES:
        names |= set(re.findall(r"^function\s+(" + WORD + ")", (ROOT / "scripts" / line).read_text(), re.M))
    missing = {}
    for path in sorted((ROOT / "rust" / "src").rglob("*.rs")):
        if str(path.relative_to(ROOT)) in SHELL_ONLY:
            continue
        text = re.split(r"#\[cfg\(test\)\]\s*mod\b", path.read_text(), maxsplit=1)[0]
        for name in set(re.findall(r'"(' + WORD + ')"', text)) & names:
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
