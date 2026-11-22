#!/usr/bin/env python3
"""Small inputs for comparing install steps with tests/install_step_compare.py.

A step that needs another step's results reads them through that step's
receipt. This writes such receipts, and the small files they name, for the
steps whose real results take hours to build (D1's VEP-annotated ClinVar VCF
and the like), so the steps that read them can be compared both ways:

  python3 tests/install_fixtures.py clinvar --out DIR
      DIR/receipts/D1.json and DIR/D1/data/ClinVar/...: a bgzipped, indexed
      ClinVar VCF with the INFO and CSQ fields of every ClinVar stage.
  python3 tests/install_fixtures.py alphamissense --out DIR
      DIR/receipts/C4.json and DIR/C4/data/AlphaMissense/...: a bgzipped,
      indexed AlphaMissense VCF as VEP and gnomAD annotated it (C3, C4), with
      an intolerant and a tolerant Pfam domain.

Then, for example:

  python3 tests/install_step_compare.py --step D8 --input DIR/receipts/D1.json ...
"""
import argparse
import json
import subprocess
from pathlib import Path

CSQ = [
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature",
    "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "Protein_position", "DOMAINS",
    "LoF", "NMD", "SpliceVault_top_events", "SpliceVault_SpliceAI_delta",
    "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL",
]
POPS = ["nfe", "eas", "afr", "amr", "asj", "fin", "sas", "mid", "remaining"]
FREQUENCIES = (
    ["AC_joint", "AN_joint", "AF_joint", "nhomalt_joint", "AC_grpmax_joint", "AF_grpmax_joint",
     "AN_grpmax_joint", "nhomalt_grpmax_joint"]
    + [f"AF_joint_{p}_{s}" for p in POPS for s in ("XX", "XY")]
    + [f"nhomalt_joint_{p}_XX" for p in POPS]
)


def receipt(step: str, name: str, folder: Path, keys: dict) -> dict:
    return {
        "schema": 1, "id": step, "name": name, "lineage": "fixture", "folder": str(folder),
        "keys": {k: str(v) for k, v in keys.items()}, "shipped": {},
        "content": "f" * 64, "files": [], "inputs": [],
        "finished": "2026-11-22T00:00:00+00:00", "host": "fixture", "stub": False,
    }


def csq(**fields) -> str:
    return "|".join(str(fields.get(name, "")) for name in CSQ)


def clinvar_records():
    """(chrom, pos, ref, alt, info, [csq]) of a few ClinVar variants."""
    gene_a = dict(SYMBOL="GENEA", Gene="ENSG00000000001", Feature_type="Transcript",
                  Feature="ENST00000000001", BIOTYPE="protein_coding")
    gene_b = dict(SYMBOL="GENEB", Gene="ENSG00000000002", Feature_type="Transcript",
                  Feature="ENST00000000002", BIOTYPE="protein_coding")
    rows = []
    for n, (pos, ref, alt, aa, residue, sig, stars, af) in enumerate([
        (1000, "C", "T", "p.Arg34Trp", 34, "Pathogenic", "reviewed_by_expert_panel", 0.00001),
        (1003, "G", "A", "p.Arg35Gln", 35, "Likely_pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts", 0.00002),
        (1006, "A", "G", "p.Lys36Glu", 36, "Pathogenic", "criteria_provided,_single_submitter", 0.0),
        (1009, "C", "A", "p.Ser37Ter", 37, "Pathogenic", "reviewed_by_expert_panel", 0.00005),
        (1012, "T", "C", "p.Leu38Pro", 38, "Benign", "criteria_provided,_multiple_submitters,_no_conflicts", 0.02),
        (1015, "G", "T", "p.Arg34Leu", 34, "Pathogenic", "practice_guideline", 0.00003),
    ], start=1):
        consequence = "stop_gained" if aa.endswith("Ter") else "missense_variant"
        impact = "HIGH" if consequence == "stop_gained" else "MODERATE"
        fields = dict(gene_a, Allele=alt, Consequence=consequence, IMPACT=impact, EXON="2/5",
                      HGVSc=f"ENST00000000001.1:c.{residue * 3}{ref}>{alt}",
                      HGVSp=f"ENSP00000000001.1:{aa}", Protein_position=f"{residue}/120",
                      DOMAINS="Pfam:PF00001", LoF="HC" if impact == "HIGH" else "",
                      NMD="NMD_escaping_variant" if impact == "HIGH" else "",
                      SpliceAI_pred_DS_AG="0.00", SpliceAI_pred_DS_AL="0.01",
                      SpliceAI_pred_DS_DG="0.00", SpliceAI_pred_DS_DL="0.02")
        rows.append(("chr1", pos, ref, alt, n, sig, stars, af, [csq(**fields)]))
    for n, (pos, ref, alt, aa, residue, sig, stars, af) in enumerate([
        (5000, "A", "T", "p.Asp10Val", 10, "Pathogenic", "reviewed_by_expert_panel", 0.0001),
        (5003, "C", "G", "p.Pro11Arg", 11, "Pathogenic/Likely_pathogenic", "criteria_provided,_multiple_submitters,_no_conflicts", 0.00004),
        (5100, "G", "A", "", "", "Pathogenic", "reviewed_by_expert_panel", 0.0),
    ], start=101):
        consequence = "missense_variant" if aa else "splice_donor_variant"
        fields = dict(gene_b, Allele=alt, Consequence=consequence,
                      IMPACT="MODERATE" if aa else "HIGH", EXON="4/9" if aa else "",
                      INTRON="" if aa else "4/8",
                      HGVSc=f"ENST00000000002.1:c.{pos - 4970}{ref}>{alt}",
                      HGVSp=f"ENSP00000000002.1:{aa}" if aa else "", Protein_position=f"{residue}/90" if aa else "",
                      SpliceAI_pred_DS_DL="0.91" if not aa else "0.00",
                      SpliceVault_top_events="ES:2:5%" if not aa else "")
        rows.append(("chr2", pos, ref, alt, n, sig, stars, af, [csq(**fields)]))
    return rows


def write_clinvar(out: Path):
    folder = out / "D1"
    data = folder / "data" / "ClinVar"
    data.mkdir(parents=True, exist_ok=True)
    plain = data / "clinvar.GRCh38.vep.vcf"
    lines = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=248956422>",
             "##contig=<ID=chr2,length=242193529>"]
    # As in ClinVar's own VCF: CLNSIG and CLNREVSTAT have Number=., which the
    # builders join into one text (a Number=1 text would be joined letter by
    # letter and never match).
    lines.append('##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="ClinVar allele ID">')
    lines += [f'##INFO=<ID={k},Number=.,Type=String,Description="ClinVar {k}">'
              for k in ("CLNSIG", "CLNREVSTAT")]
    lines += [f'##INFO=<ID={k},Number=A,Type=Float,Description="gnomAD {k}">' for k in FREQUENCIES]
    lines.append('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from '
                 f'Ensembl VEP. Format: {"|".join(CSQ)}">')
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    for chrom, pos, ref, alt, allele, sig, stars, af, annotations in clinvar_records():
        info = [f"ALLELEID={allele}", f"CLNSIG={sig}", f"CLNREVSTAT={stars}"]
        for key in FREQUENCIES:
            if key.startswith("AF"):
                info.append(f"{key}={af:g}")
            elif key.startswith("AN"):
                info.append(f"{key}=152000")
            else:
                info.append(f"{key}={round(af * 152000)}" if key.startswith("AC") else f"{key}=0")
        info.append("CSQ=" + ",".join(annotations))
        lines.append(f"{chrom}\t{pos}\t{allele}\t{ref}\t{alt}\t.\t.\t{';'.join(info)}")
    plain.write_text("\n".join(lines) + "\n")
    subprocess.run(["bgzip", "-f", str(plain)], check=True)
    vcf = Path(str(plain) + ".gz")
    subprocess.run(["tabix", "-f", "-p", "vcf", str(vcf)], check=True)
    keys = {"clinvar_vcf": vcf}
    for key in ("clinvar_aa_stat", "clinvar_splice_stat", "clinvar_gene_stat"):
        stand_in = data / f"{key}.fixture"
        stand_in.write_text(f"{key} stand-in: D6, D8 and D9 do not read it\n")
        keys[key] = stand_in
    receipts = out / "receipts"
    receipts.mkdir(parents=True, exist_ok=True)
    (receipts / "D1.json").write_text(json.dumps(receipt("D1", "CLINVAR", folder, keys), indent=2) + "\n")
    print(receipts / "D1.json")


AM_CSQ = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature",
          "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "Protein_position", "Amino_acids",
          "DOMAINS", "CANONICAL"]
AMINO = "ACDEFGHIKLMNPQRSTVWY"


def alphamissense_records():
    """(chrom, pos, ref, alt, score, af, pvar, csq) of four genes. GENEA's
    Pfam PF00069 and GENEB's PF00001 domains have GO molecular functions in
    InterPro (DAS's functional references); GENEC's and GENED's domains are
    not in InterPro (DAS's queries). GENEA's and GENEC's domains are
    intolerant, GENEB's and GENED's are not."""
    rows = []
    genes = [("chr1", 10000, "GENEA", "ENSG00000000101", "ENST00000000101", 120, (30, 90), "Pfam:PF00069", True),
             ("chr1", 50000, "GENEC", "ENSG00000000103", "ENST00000000103", 80, (10, 60), "Pfam:PF99901", True),
             ("chr2", 20000, "GENEB", "ENSG00000000102", "ENST00000000102", 90, (20, 70), "Pfam:PF00001", False),
             ("chr2", 60000, "GENED", "ENSG00000000104", "ENST00000000104", 70, (15, 55), "Pfam:PF99902", False)]
    for chrom, start, symbol, gene, transcript, length, (first, last), domain, intolerant in genes:
        exons = 3
        for residue in range(1, length + 1):
            pos = start + 3 * residue
            exon = min(exons, 1 + (residue - 1) * exons // length)
            inside = first <= residue <= last
            for n, alt in enumerate("CGT"):
                ref_aa = AMINO[residue % 20]
                alt_aa = AMINO[(residue + n + 1) % 20]
                if inside and intolerant:
                    score = 0.80 + ((residue * 7 + n * 3) % 19) / 100
                elif inside:
                    score = 0.20 + ((residue * 5 + n) % 23) / 100
                else:
                    score = 0.05 + ((residue * 3 + n * 11) % 29) / 100
                af = 0.0 if (residue + n) % 9 else 0.00002 * (1 + n)
                if residue % 41 == 0 and n == 0:
                    af = 0.07
                pvar = f"{ref_aa}{residue}{alt_aa}"
                fields = dict(Allele=alt, Consequence="missense_variant", IMPACT="MODERATE", SYMBOL=symbol,
                              Gene=gene, Feature_type="Transcript", Feature=transcript, BIOTYPE="protein_coding",
                              EXON=f"{exon}/{exons}", HGVSc=f"{transcript}.1:c.{3 * residue}A>{alt}",
                              HGVSp=f"ENSP{transcript[4:]}.1:p.{pvar}", Protein_position=f"{residue}/{length}",
                              Amino_acids=f"{ref_aa}/{alt_aa}", DOMAINS=domain if inside else "", CANONICAL="YES")
                rows.append((chrom, pos, "A", alt, round(score, 4), af, pvar,
                             "|".join(str(fields.get(k, "")) for k in AM_CSQ)))
    return rows


def write_alphamissense(out: Path):
    folder = out / "C4"
    data = folder / "data" / "AlphaMissense"
    data.mkdir(parents=True, exist_ok=True)
    plain = data / "AlphaMissense_hg38.vep.vcf"
    lines = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=248956422>",
             "##contig=<ID=chr2,length=242193529>",
             '##INFO=<ID=AM_PATHOGENICITY,Number=1,Type=Float,Description="AlphaMissense score">',
             '##INFO=<ID=PVAR,Number=1,Type=String,Description="Protein variant">',
             '##INFO=<ID=AF_grpmax_joint,Number=A,Type=Float,Description="gnomAD grpmax AF">',
             '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
             f'Format: {"|".join(AM_CSQ)}">',
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for chrom, pos, ref, alt, score, af, pvar, annotation in alphamissense_records():
        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t"
                     f"AM_PATHOGENICITY={score:g};PVAR={pvar};AF_grpmax_joint={af:g};CSQ={annotation}")
    plain.write_text("\n".join(lines) + "\n")
    subprocess.run(["bgzip", "-f", str(plain)], check=True)
    vcf = Path(str(plain) + ".gz")
    subprocess.run(["tabix", "-f", "-p", "vcf", str(vcf)], check=True)
    receipts = out / "receipts"
    receipts.mkdir(parents=True, exist_ok=True)
    (receipts / "C4.json").write_text(json.dumps(
        receipt("C4", "ALPHAMISSENSE_GNOMAD", folder, {"alphamissense_vep_vcf": vcf}), indent=2) + "\n")
    print(receipts / "C4.json")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("what", choices=["clinvar", "alphamissense"])
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = Path(args.out).resolve()
    {"clinvar": write_clinvar, "alphamissense": write_alphamissense}[args.what](out)


if __name__ == "__main__":
    main()
