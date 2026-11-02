#!/usr/bin/env python3
"""Source-matched ClinVar gene-statistics fixtures and representation probes.

  python clinvar_gene_stat_compare.py probe --producer ORIGINAL.py --work-root NEW_DIR

The pinned Python producer controls values and failure outcomes. Fixture payloads,
complete pickles, logs and detailed comparisons stay in the specified work root.
"""
from __future__ import annotations
import argparse
from contextlib import redirect_stderr, redirect_stdout
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import pickle
import pickletools
import sys

PRODUCER_SHA256 = "8abdd1964299e69d375eef83bca4ed39d534d1f70f14b96dfae80426ee14f755"
FIELDS = ["SYMBOL", "Gene", "BIOTYPE", "Feature", "DOMAINS", "Consequence", "LoF", "NMD",
          "SpliceVault_top_events", "SpliceVault_SpliceAI_delta", "SpliceAI_pred_DS_DL"]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def producer(path):
    if digest(path) != PRODUCER_SHA256:
        raise RuntimeError("Python producer source differs from the pinned fixture scope")
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location("clinvar_stat_variants", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def write_inputs(directory, symbols):
    import pysam
    from clinvar_annotation_contract import TABLE_INFO_FIELDS
    directory.mkdir()
    header = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=1000000>"]
    for name in sorted(TABLE_INFO_FIELDS):
        header.append(f'##INFO=<ID={name},Number=.,Type=Float,Description="Fixture frequency">')
    for name in ("CLNSIG", "CLNREVSTAT", "CLNDN"):
        header.append(f'##INFO=<ID={name},Number=.,Type=String,Description="Fixture annotation">')
    header.extend(['##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="Fixture identifier">',
                   '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: '+"|".join(FIELDS)+'">',
                   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"])
    lines = []
    summary = ["#AlleleID\tAssembly\tName"]
    mane = ["RefSeq_nuc\tEnsembl_nuc"]
    for index, symbol in enumerate(symbols, 1):
        values = [symbol, f"ENSG_{index}", "protein_coding", f"ENST_{index}", "Pfam:PF1&InterPro:IPR1",
                  "missense_variant", "", "", "NA", "0.0", "0.0"]
        lines.append(f"chr1\t{index}\t.\tA\tT\t.\tPASS\tALLELEID={index};CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;CSQ="+"|".join(values))
        summary.append(f"{index}\tGRCh38\tNM_{index:06d}.1:c.1A>T")
        mane.append(f"NM_{index:06d}.1\tENST_{index}.1")
    plain = directory/"input.vcf"
    plain.write_text("\n".join(header+lines)+"\n")
    vcf = directory/"input.vcf.gz"
    pysam.tabix_compress(str(plain), str(vcf), force=False)
    (directory/"summary.tsv").write_text("\n".join(summary)+"\n")
    (directory/"mane.tsv").write_text("\n".join(mane)+"\n")
    return vcf, directory/"summary.tsv", directory/"mane.tsv"


def probe(args):
    import numpy as np
    args.work_root.mkdir(parents=True, exist_ok=False)
    original = producer(args.producer)
    cases = {"text": ["FIRST", "SECOND"], "all_missing": ["", ""], "mixed_missing_text": ["", "TEXT"],
             "integers": ["1", "2"], "floats": ["0.1", "0.2"], "mixed_numeric_text": ["1", "TEXT"],
             "nullable_integer": ["", "1"], "nullable_float": ["", "0.1"], "boolean_text": ["True", "False"]}
    reports = []
    for name, symbols in cases.items():
        directory = args.work_root/name
        vcf, summary, mane = write_inputs(directory, symbols)
        output = directory/"reference.pkl.gz"
        with (directory/"producer.log").open("w") as log, redirect_stdout(log), redirect_stderr(log):
            try:
                value = original.main_stat_variants(str(vcf), str(output), 1,
                    clinvar_tab=str(summary), mane_tranx_tab=str(mane), donor_loss_cutoff=.25)
            except Exception as error:
                reports.append({"case": name, "success": False, "exception": type(error).__name__, "message": str(error)})
                continue
        items = []
        for gene in value.values():
            symbol = gene["gene_symbol"]
            entry = {"type": type(symbol).__module__+"."+type(symbol).__qualname__}
            if isinstance(symbol, np.generic):
                entry.update(dtype=str(symbol.dtype), scalar_bytes=symbol.tobytes().hex())
            items.append(entry)
        with gzip.open(output, "rb") as stream:
            payload = stream.read()
        with (directory/"pickle_opcodes.txt").open("w") as out:
            pickletools.dis(payload, out=out)
        reports.append({"case": name, "success": True, "genes": len(value), "gene_symbol_types": items,
                        "pickle_sha256": digest(output), "opcodes": sorted({opcode.name for opcode, _, _ in pickletools.genops(payload)})})
    report = {"producer_sha256": PRODUCER_SHA256, "cases": reports}
    (args.work_root/"probe.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    sub = commands.add_parser("probe")
    sub.add_argument("--producer", type=Path, required=True)
    sub.add_argument("--work-root", type=Path, required=True)
    sub.set_defaults(function=probe)
    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
