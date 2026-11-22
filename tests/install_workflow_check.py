#!/usr/bin/env python3
"""Regenerates workflows/install.nf from conf/install_steps.json.

Default: exit 1 when the committed workflow differs from the registry.
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
        if "module" not in s:
            raise SystemExit(f"{s['id']} names no module; every step is a process of modules/local/install_steps/")
        module = ROOT / s["module"]
        if not s["module"].startswith("modules/local/install_steps/") or not module.is_file():
            raise SystemExit(f"{s['id']}: its module {s['module']} is not a file of modules/local/install_steps/")
        if not re.search(rf"^process {process(s)} \{{", module.read_text(), re.M):
            raise SystemExit(f"{s['id']}: {s['module']} defines no process {process(s)}")
    modules = "\n".join(
        f"include {{ {process(s)} as {alias(s)} }} from '../{s['module'][:-len('.nf')]}'"
        for s in steps)
    lines = [HEADER.read_text().replace("@PROCESS_INCLUDES@", modules).rstrip("\n")]
    for s in steps:
        needs = ", ".join(f"{alias(ids[n])}.out.receipt" for n in s["needs"])
        stage = ", ".join(f"'{k}'" for k in s.get("stage", []))
        lines.append(f"    {alias(s)}(installNative('{s['id']}', plan, [{needs}], [{stage}]))")
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

if __name__ == "__main__":
    target = ROOT / "workflows" / "install.nf"
    text = render()
    if "--write" in sys.argv[1:]:
        target.write_text(text)
    elif target.read_text() != text:
        print("workflows/install.nf differs from conf/install_steps.json; run with --write")
        sys.exit(1)
