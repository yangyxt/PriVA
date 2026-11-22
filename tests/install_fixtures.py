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
  python3 tests/install_fixtures.py mechanism --out DIR
      DIR/receipts/F1.json and DIR/F1/data/...: a MONDO disease-scope registry
      and scoped HPO assertions (F1) for four genes, and DIR/net/HOST/PATH:
      the G2P and Orphadata downloads the mechanism evidence (F4) reads, for
      tests/fake_net.py.
  python3 tests/install_fixtures.py nonlof --out DIR
      DIR/receipts/F5.json and DIR/F5/data/gofcards/...: a normalized GoFCards
      cache without ClinVar conditions (what F5 builds), and DIR/net/HOST/PATH:
      a weekly ClinVar VCV release (with its .md5, README and XSD) holding one
      of its variants, G2P and ClinGen dosage tables, for the ClinVar
      injection (F7) and the non-LOF mechanism sources and cache (F6, F8).
      PanelApp and the GoFCards workbook are left out (the server's 404).
  python3 tests/install_fixtures.py gnomad --out DIR
      DIR/net/HOST/PATH: gnomAD's 24 joint sites chromosomes (bgzipped, with
      tabix indexes) holding a few variants each, placed by the shipped
      hg38ToHg19 chain so that they lift to the same GRCh37 chromosome, to
      another one, onto the reverse strand, or nowhere; DIR/receipts/A1.hg19.json
      and DIR/receipts/A1.hg38.json: references of step A1. The GRCh37 one has
      the real chromosome lengths and holds bases only where the lifted
      variants land (a sparse file with its own .fai).
  python3 tests/install_fixtures.py clinvar_release --out DIR
      DIR/net/HOST/PATH: NCBI's weekly ClinVar VCFs of both assemblies (with
      a ##fileDate, a tabix index and the publisher's .md5), variant_summary
      (with its .md5) and MANE's summary, as the ClinVar download (D0) reads
      them from the sources of cache_resources.json.
  python3 tests/install_fixtures.py vep --out DIR
      DIR/envs/priva_acmg: an environment (CONDA_PREFIX) whose bin/ holds
      stand-ins for vep, vep_install (Ensembl's INSTALL.pl), perl and conda,
      with the Ensembl API and install helper modules naming release 116.
      vep_install logs its arguments to DIR/vep_install.log and writes a
      small converted cache and the requested plugins (UTRAnnotator: the
      pinned upstream file) where it is told to, as the VEP cache step (A2)
      runs it; DIR/receipts/A2.{hg38,hg19}.json name such a cache and
      plugins. vep annotates a VCF (logging to DIR/vep.log): a ##VEP line
      naming the cache and assembly, and a CSQ field from the transcripts of
      the alphamissense fixtures' four genes.
  python3 tests/install_fixtures.py alphamissense_inputs --out DIR
      DIR/receipts/{A1,B2,B3,C1}.{hg38,hg19}.json: a small reference, a
      conservation track, the publisher's AlphaMissense table (bgzipped,
      tabix-indexed, the variants of the alphamissense fixture) and a gnomAD
      family holding some of them: the inputs of the AlphaMissense
      annotation (C3) besides A2's, and of its gnomAD step (C4).
  python3 tests/install_fixtures.py spliceai --out DIR
      DIR/spliceai/v1/{snv,indel}.vcf.gz: SpliceAI v1 VCFs with tabix
      indexes (no assembly in their headers, so either assembly takes them),
      and DIR/spliceai/native/: the 47 native SpliceAI2 shards (sizes and
      SHA-256 pinned) and their downloads-complete manifest, for the SpliceAI
      step (B4). Without --plugin the shards hold headers only and the
      manifest no gene-context resources, so the preparation stops; with
      --plugin PATH (Ensembl's SpliceAI.pm of VEP 116, which the preparation
      pins by SHA-256) the chr1 shards hold scores inside two genes and the
      manifest pins a GENCODE 48 GTF, scored regions, a GENCODE 50 GTF that
      renames one gene, and the plugin, so the preparation completes.
  python3 tests/install_fixtures.py ranking --out DIR
      DIR/net/HOST/PATH: the public inputs of the HPO gene ranking (F10), the
      small made-up files of rust/tests/hpo_gene_ranking, compressed as the
      sources serve them; one optional Human Protein Atlas table is left out.

Then, for example:

  python3 tests/install_step_compare.py --step D8 --input DIR/receipts/D1.json ...
"""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import shutil
import zipfile
import subprocess
import sys
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
GNOMAD_FIELDS = (
    ["AC_joint", "AN_joint", "AF_joint", "nhomalt_joint", "AC_grpmax_joint", "AF_grpmax_joint",
     "AN_grpmax_joint", "nhomalt_grpmax_joint", "AC_joint_XX", "AN_joint_XX", "AC_joint_XY",
     "AN_joint_XY", "AF_joint_XX", "AF_joint_XY", "nhomalt_joint_XX", "nhomalt_joint_XY"]
    + [f"AF_joint_{p}_{s}" for p in POPS for s in ("XX", "XY")]
    + [f"nhomalt_joint_{p}_XX" for p in POPS])
# (chrom, first codon - 3, symbol, gene, transcript, length, domain residues,
# domain, intolerant) of the four AlphaMissense genes.
AM_GENES = [("chr1", 10000, "GENEA", "ENSG00000000101", "ENST00000000101", 120, (30, 90), "Pfam:PF00069", True),
            ("chr1", 50000, "GENEC", "ENSG00000000103", "ENST00000000103", 80, (10, 60), "Pfam:PF99901", True),
            ("chr2", 20000, "GENEB", "ENSG00000000102", "ENST00000000102", 90, (20, 70), "Pfam:PF00001", False),
            ("chr2", 60000, "GENED", "ENSG00000000104", "ENST00000000104", 70, (15, 55), "Pfam:PF99902", False)]


def alphamissense_records():
    """(chrom, pos, ref, alt, score, af, pvar, csq) of four genes. GENEA's
    Pfam PF00069 and GENEB's PF00001 domains have GO molecular functions in
    InterPro (DAS's functional references); GENEC's and GENED's domains are
    not in InterPro (DAS's queries). GENEA's and GENEC's domains are
    intolerant, GENEB's and GENED's are not."""
    rows = []
    for chrom, start, symbol, gene, transcript, length, (first, last), domain, intolerant in AM_GENES:
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


REGISTRY = [
    "disease_id", "mondo_id", "mondo_name", "mondo_version", "hpo_mendelian_inheritance",
    "hpo_non_monogenic_inheritance", "mondo_hereditary", "mondo_neoplasm", "mondo_predisposition",
    "mondo_omim_susceptibility", "disease_scope", "priva_scope", "scope_evidence", "scope_reference",
    "scope_review_status", "scope_note",
]
ASSERTIONS = [
    "gene_symbol", "disease_id", "hpo_id", "frequency", "evidence", "reference", "mondo_id",
    "mondo_name", "disease_scope", "priva_scope", "scope_evidence", "scope_reference",
    "scope_review_status",
]
G2P = [
    "gene symbol", "disease name", "molecular mechanism", "confidence", "allelic requirement",
    "cross cutting modifier", "phenotypes", "pmids", "g2p id", "disease mim", "disease MONDO",
]
# (gene, disease id, MONDO id, name, scope, PriVA scope)
DISEASES = [
    ("GENEA", "OMIM:100100", "MONDO:0000101", "condition A", "mendelian_non_neoplastic", "include"),
    ("GENEB", "OMIM:100200", "MONDO:0000102", "condition B", "mendelian_non_neoplastic", "include"),
    ("GENEC", "ORPHA:1001", "MONDO:0000103", "condition C", "mendelian_non_neoplastic", "include"),
    ("GENED", "ORPHA:1002", "MONDO:0000104", "condition D", "mendelian_susceptibility_uncertain", "review"),
]


def write_alphamissense_inputs(out: Path):
    """The receipts C3 reads besides A2's: a small reference (A1), a
    conservation track (B2) and the publisher's AlphaMissense table (B3) for
    each assembly; the table holds the variants of alphamissense_records."""
    receipts = out / "receipts"
    receipts.mkdir(parents=True, exist_ok=True)
    for assembly in ("hg38", "hg19"):
        a1 = out / "A1" / assembly / "data" / "reference"
        a1.mkdir(parents=True, exist_ok=True)
        fasta = a1 / f"{assembly}.fa"
        fasta.write_text("".join(f">{chrom}\n" + "ACGT" * 20000 + "\n" for chrom in ("chr1", "chr2")))
        subprocess.run(["samtools", "faidx", fasta], check=True)
        (receipts / f"A1.{assembly}.json").write_text(json.dumps(
            receipt("A1", "REFERENCE_FASTA", a1.parents[1], {"ref_genome": fasta}), indent=2) + "\n")
        b2 = out / "B2" / assembly / "data" / "Conservation"
        b2.mkdir(parents=True, exist_ok=True)
        track = b2 / f"{assembly}.phyloP100way.bw"
        track.write_bytes(b"stand-in bigWig of " + assembly.encode() + b"\n")
        (receipts / f"B2.{assembly}.json").write_text(json.dumps(
            receipt("B2", "CONSERVATION", b2.parents[1], {"conservation_file": track}), indent=2) + "\n")
        b3 = out / "B3" / assembly / "data" / "AlphaMissense"
        b3.mkdir(parents=True, exist_ok=True)
        lines = ["# Copyright 2023 DeepMind Technologies Limited", "#",
                 "# Licensed under CC BY-NC-SA 4.0 license",
                 "#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class"]
        for chrom, pos, ref, alt, score, _af, pvar, annotation in alphamissense_records():
            transcript = annotation.split("|")[AM_CSQ.index("Feature")]
            label = "likely_pathogenic" if score > 0.564 else "likely_benign" if score < 0.34 else "ambiguous"
            lines.append(f"{chrom}\t{pos}\t{ref}\t{alt}\t{assembly}\tP{transcript[-5:]}\t{transcript}.1\t{pvar}\t{score:g}\t{label}")
        plain = b3 / f"AlphaMissense_{assembly}.tsv"
        plain.write_text("\n".join(lines) + "\n")
        subprocess.run(["bgzip", "-f", plain], check=True)
        table_gz = Path(f"{plain}.gz")
        subprocess.run(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", "-S", "1", table_gz], check=True)
        (receipts / f"B3.{assembly}.json").write_text(json.dumps(
            receipt("B3", "ALPHAMISSENSE_TSV", b3.parents[1], {"alphamissense_prescore": table_gz}), indent=2) + "\n")
        # gnomAD (C1): the 24 chromosomes, with joint frequencies at every
        # fifth AlphaMissense residue of chr1 and chr2 (one allele each) and
        # a record elsewhere, as C4 reads them beside gnomad_vcf_chrX.
        c1 = out / "C1" / assembly / "data" / "gnomAD"
        c1.mkdir(parents=True, exist_ok=True)
        header = ["##fileformat=VCFv4.2"]
        for field in GNOMAD_FIELDS:
            header.append(f'##INFO=<ID={field},Number={"1" if field.startswith("AN") else "A"},'
                          f'Type={"Float" if field.startswith("AF") else "Integer"},Description="{field}">')
        sites = {}
        for chrom, pos, ref, alt, _score, _af, pvar, _annotation in alphamissense_records():
            if int(pvar[1:-1]) % 5 == 0 and alt == "CGT"[int(pvar[1:-1]) % 3]:
                sites.setdefault(chrom, []).append((pos, ref, alt))
        for chrom in [f"chr{n}" for n in range(1, 23)] + ["chrX", "chrY"]:
            lines = []
            for number, (pos, ref, alt) in enumerate(sites.get(chrom, []) + [(70000, "G", "A")]):
                values = [f"{field}={(number % 7 + 1) / 1000:.5f}" if field.startswith("AF")
                          else f"{field}={(number * 37) % 5000 + 1}" for field in GNOMAD_FIELDS]
                lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{';'.join(values)}")
            name = (f"gnomad.joint.v4.1.sites.{chrom}.vcf" if assembly == "hg38"
                    else f"gnomad.joint.v4.1.sites.hg19.{chrom}.vcf")
            plain = c1 / name
            plain.write_text("\n".join(header + [f"##contig=<ID={chrom}>", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
                                       + lines) + "\n")
            target = c1 / (name + (".bgz" if assembly == "hg38" else ".gz"))
            with open(target, "wb") as handle:
                subprocess.run(["bgzip", "-c", plain], stdout=handle, check=True)
            plain.unlink()
            subprocess.run(["tabix", "-f", "-p", "vcf", target], check=True)
        chrx = c1 / ("gnomad.joint.v4.1.sites.chrX.vcf.bgz" if assembly == "hg38"
                     else "gnomad.joint.v4.1.sites.hg19.chrX.vcf.gz")
        (receipts / f"C1.{assembly}.json").write_text(json.dumps(
            receipt("C1", "GNOMAD", c1.parents[1], {"gnomad_vcf_chrX": chrx}), indent=2) + "\n")
    print(receipts)


def table(names, rows, delimiter="\t") -> str:
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, delimiter=delimiter, lineterminator="\n")
    writer.writerow(names)
    writer.writerows(rows)
    return stream.getvalue()


def write_gz(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.GzipFile(path, "wb", mtime=0) as handle:
        handle.write(text.encode())


def write_mechanism(out: Path):
    folder = out / "F1"
    registry = folder / "data/mondo/disease_scope.tsv.gz"
    assertions = folder / "data/hpo/genes_to_phenotype.assertions.tsv.gz"
    obo = folder / "data/mondo/mondo-simple.fixture.obo"
    rows = []
    for gene, disease, mondo, name, scope, priva in DISEASES:
        rows.append([disease, mondo, name, "fixture", "True", "False", "True", "False", "False", "False",
                     scope, priva, "fixture", "fixture", "automatic", ""])
    write_gz(registry, table(REGISTRY, rows))
    hpo = []
    for gene, disease, mondo, name, scope, priva in DISEASES:
        for term, frequency in (("HP:0001250", "HP:0040281"), ("HP:0000006", ""), ("HP:0003829", "")):
            hpo.append([gene, disease, term, frequency, "PCS", "PMID:1", mondo, name, scope, priva,
                        "fixture", "fixture", "automatic"])
    write_gz(assertions, table(ASSERTIONS, hpo))
    obo.write_text("format-version: 1.2\ndata-version: mondo/releases/fixture/mondo-simple.obo\n\n"
                   "[Term]\nid: MONDO:0003847\nname: hereditary disease\n")
    net = out / "net"
    g2p = net / "www.ebi.ac.uk/gene2phenotype/api/panel/all/download"
    g2p.parent.mkdir(parents=True, exist_ok=True)
    g2p.write_text(table(G2P, [
        ["GENEA", "condition A", "loss of function", "definitive", "monoallelic_autosomal", "",
         "HP:0001250;HP:0003829", "123;456", "G2P00001", "100100", "MONDO:0000101"],
        ["GENEB", "condition B", "altered gene product structure", "strong", "monoallelic_autosomal",
         "typically mosaic", "HP:0001250", "789", "G2P00002", "100200", "MONDO:0000102"],
        ["GENEC", "condition C", "dominant negative", "moderate", "biallelic_autosomal", "", "",
         "", "G2P00003", "", "MONDO:0000103"],
        ["GENED", "condition D", "undetermined", "limited", "monoallelic_X_hemizygous", "", "", "",
         "G2P00004", "", ""],
    ], delimiter=","))
    association = ("<DisorderGeneAssociation><DisorderGeneAssociationStatus><Name>{status}</Name>"
                   "</DisorderGeneAssociationStatus><DisorderGeneAssociationType><Name>{kind}</Name>"
                   "</DisorderGeneAssociationType><Gene><Symbol>{gene}</Symbol></Gene>"
                   "<SourceOfValidation>{pmids}</SourceOfValidation></DisorderGeneAssociation>")
    disorders = []
    for code, name, gene, kind, status in (
            ("1001", "condition C", "GENEC", "Disease-causing germline mutation(s) (loss of function) in", "Assessed"),
            ("1002", "condition D", "GENED", "Disease-causing germline mutation(s) (gain of function) in", "Assessed"),
            ("1003", "condition E", "GENEA", "Disease-causing germline mutation(s) (gain of function) in", "Not yet assessed")):
        disorders.append(f"<Disorder><OrphaCode>{code}</OrphaCode><Name>{name}</Name><DisorderGeneAssociationList>"
                         + association.format(status=status, kind=kind, gene=gene, pmids="111[PMID]_222[PMID]")
                         + "</DisorderGeneAssociationList></Disorder>")
    orpha = net / "www.orphadata.com/data/xml/en_product6.xml"
    orpha.parent.mkdir(parents=True, exist_ok=True)
    orpha.write_text('<?xml version="1.0" encoding="UTF-8"?>\n<JDBOR><DisorderList>'
                     + "".join(disorders) + "</DisorderList></JDBOR>\n")
    receipts = out / "receipts"
    receipts.mkdir(parents=True, exist_ok=True)
    (receipts / "F1.json").write_text(json.dumps(receipt("F1", "HPO_MONDO", folder, {
        "hpo_assertions": assertions, "mondo_obo": obo, "mondo_disease_scope_registry": registry}), indent=2) + "\n")
    print(receipts / "F1.json")


def injection_fixtures():
    """The hand-made GoFCards cache and VCV records of the injection's own comparison."""
    path = Path(__file__).resolve().parents[1] / "rust/tests/gofcards_clinvar_injection_compare.py"
    spec = importlib.util.spec_from_file_location("gofcards_clinvar_injection_compare", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_nonlof(out: Path):
    injection = injection_fixtures()
    folder = out / "F5"
    normalized = folder / "data/gofcards/build_work/gofcards_exact_gof.normalized.json.gz"
    write_gz(normalized, json.dumps(injection.base_cache(), sort_keys=True))
    net = out / "net"
    ncbi = net / "ftp.ncbi.nlm.nih.gov/pub/clinvar"
    release = ncbi / "xml/weekly_release/ClinVarVCVRelease_00-latest_weekly.xml.gz"
    write_gz(release, injection.document(injection.archive()))
    digest = hashlib.md5(release.read_bytes()).hexdigest()
    (release.parent / (release.name + ".md5")).write_text(f"{digest}  {release.name}\n")
    (ncbi / "xml/_README").write_text("ClinVar XML files (fixture)\n")
    (ncbi / "xsd_public").mkdir(parents=True, exist_ok=True)
    (ncbi / "xsd_public/_index").write_text('<html><body><a href="ClinVar_VCV_2.4.xsd">ClinVar_VCV_2.4.xsd</a>\n'
                                            '<a href="ClinVar_VCV_2.3.xsd">ClinVar_VCV_2.3.xsd</a></body></html>\n')
    (ncbi / "xsd_public/ClinVar_VCV_2.4.xsd").write_text('<?xml version="1.0"?>\n'
                                                       '<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema"/>\n')
    g2p = net / "www.ebi.ac.uk/gene2phenotype/api/panel/all/download"
    g2p.parent.mkdir(parents=True, exist_ok=True)
    g2p.write_text(table(G2P, [
        ["GENE1", "condition one", "gain of function", "definitive", "monoallelic_autosomal", "",
         "HP:0001250", "123", "G2P00011", "100100", "MONDO:0000101"],
        ["GENE2", "condition two", "dominant negative", "strong", "monoallelic_autosomal", "", "",
         "456", "G2P00012", "", ""],
    ], delimiter=","))
    dosage = net / "ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"
    dosage.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(Path(__file__).resolve().parents[1] / "data/clingen/gene_dosage_sensitivity.hg38.tsv", dosage)
    receipts = out / "receipts"
    receipts.mkdir(parents=True, exist_ok=True)
    (receipts / "F5.json").write_text(json.dumps(receipt("F5", "GOFCARDS_EXACT_GOF", folder, {
        "gofcards_exact_gof_cache": normalized}), indent=2) + "\n")
    print(receipts / "F5.json")


RANKING = [
    ("purl.obolibrary.org/obo/hp.owl", "../hp_fixture.owl", None),
    ("purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt", "phenotype_to_genes.txt", None),
    ("purl.obolibrary.org/obo/go/go-basic.obo", "go-basic.obo", None),
    ("current.geneontology.org/annotations/goa_human.gaf.gz", "goa_human.gaf", "gz"),
    ("ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz", "Homo_sapiens.gene_info", "gz"),
    ("ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz", "gene2go", "gz"),
    ("reactome.org/download/current/NCBI2Reactome_All_Levels.txt", "NCBI2Reactome_All_Levels.txt", None),
    ("reactome.org/download/current/ReactomePathwaysRelation.txt", "ReactomePathwaysRelation.txt", None),
    ("rest.uniprot.org/uniprotkb/stream", "uniprot_human_reviewed.tsv", "gz"),
    ("purl.obolibrary.org/obo/pr.obo", "pr.obo", None),
    ("www.informatics.jax.org/downloads/reports/MGI_GenePheno.rpt", "MGI_GenePheno.rpt", None),
    ("www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", "HOM_MouseHumanSequence.rpt", None),
    ("www.orphadata.com/data/xml/en_product6.xml", "orphadata_en_product6.xml", None),
    ("www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip", "rna_tissue_consensus.tsv", "zip"),
    ("www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip", "rna_single_cell_type.tsv", "zip"),
]


def write_ranking(out: Path):
    raw = Path(__file__).resolve().parents[1] / "rust/tests/hpo_gene_ranking/tests/fixtures/raw"
    for served, name, packing in RANKING:
        target = out / "net" / served
        target.parent.mkdir(parents=True, exist_ok=True)
        source = (raw / name).resolve()
        if packing == "gz":
            with gzip.GzipFile(target, "wb", mtime=0) as handle:
                handle.write(source.read_bytes())
        elif packing == "zip":
            with zipfile.ZipFile(target, "w") as archive:
                archive.writestr(zipfile.ZipInfo(name, (1980, 1, 1, 0, 0, 0)), source.read_bytes())
        else:
            shutil.copyfile(source, target)
    print(out / "net")


GNOMAD_BASE = "storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint"
GNOMAD_CHROMOSOMES = [f"chr{n}" for n in range(1, 23)] + ["chrX", "chrY"]


def chains(path: Path):
    """The chains of a UCSC chain file: (header fields, [(t, q, size)])."""
    header, blocks, t, q = None, [], 0, 0
    with gzip.open(path, "rt") as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "chain":
                if header:
                    yield header, blocks
                header, blocks = fields, []
                t, q = int(fields[5]), int(fields[10])
                continue
            size = int(fields[0])
            blocks.append((t, q, size))
            if len(fields) == 3:
                t, q = t + size + int(fields[1]), q + size + int(fields[2])
    if header:
        yield header, blocks


def lifted(header, t_block, q_block, at):
    """The GRCh37 0-based position of GRCh38 0-based position AT in a block."""
    q = q_block + (at - t_block)
    return q if header[9] == "+" else int(header[8]) - q - 1


def write_gnomad(out: Path):
    import random
    root = Path(__file__).resolve().parents[1]
    rng = random.Random(41)
    sizes38, sizes19, picks = {}, {}, {c: [] for c in GNOMAD_CHROMOSOMES}
    unmapped = {}
    for header, blocks in chains(root / "data/liftover/hg38ToHg19.over.chain.gz"):
        source, target = header[2], header[7]
        sizes38[source] = int(header[3])
        sizes19[target] = int(header[8])
        if source not in picks:
            continue
        long = [b for b in blocks if b[2] >= 400]
        if not long:
            continue
        kinds = {kind for kind, *_ in picks[source]}
        if target == source and header[9] == "+" and "same" not in kinds:
            picks[source].append(("same", header, long[len(long) // 2]))
            if len(blocks) > 1 and source not in unmapped:
                # A position in the gap after the first block lifts nowhere.
                gap_t = blocks[0][0] + blocks[0][2]
                if blocks[1][0] - gap_t > 10:
                    unmapped[source] = gap_t + 5
        elif target == source and header[9] == "-" and "reverse" not in kinds:
            picks[source].append(("reverse", header, long[0]))
        elif target in sizes19 and target != source and target in picks and "other" not in kinds \
                and source in ("chr1", "chr7", "chrX"):
            picks[source].append(("other", header, long[0]))
    bases = {}
    records = {c: [] for c in GNOMAD_CHROMOSOMES}
    for chrom, chosen in picks.items():
        assert any(kind == "same" for kind, *_ in chosen), chrom
        for kind, header, (t_block, q_block, size) in chosen:
            offsets = (100, 200, 300) if kind == "same" else (150,)
            for number, offset in enumerate(offsets):
                at = t_block + offset
                q = lifted(header, t_block, q_block, at)
                window = bases.setdefault(header[7], {})
                for spot in range(q - 60, q + 60):
                    window.setdefault(spot, rng.choice("ACGT"))
                ref = rng.choice("ACGT")
                if number == 1:
                    alt = ref + "AC"          # an insertion
                elif number == 2:
                    ref, alt = ref + "GT", ref  # a deletion
                else:
                    # A substitution whose lifted ALT differs from the GRCh37 base.
                    lifted_base = window[q]
                    choices = [b for b in "ACGT" if (b if header[9] == "+" else "TGCA"["ACGT".index(b)]) != lifted_base]
                    alt = choices[0]
                records[chrom].append((at + 1, ref, alt))
        if chrom in unmapped:
            records[chrom].append((unmapped[chrom] + 1, "A", "G"))
    net = out / "net" / GNOMAD_BASE
    net.mkdir(parents=True, exist_ok=True)
    info = [f"##INFO=<ID={name},Number={'1' if name.startswith(('AN', 'nhomalt')) else 'A'},"
            f"Type={'Float' if name.startswith('AF') else 'Integer'},Description=\"{name}\">"
            for name in FREQUENCIES + ["AC_joint_XX", "AN_joint_XX", "AC_joint_XY", "AN_joint_XY",
                                       "AF_joint_XX", "AF_joint_XY", "nhomalt_joint_XX", "nhomalt_joint_XY"]]
    contigs = [f"##contig=<ID={c},length={sizes38[c]},assembly=GRCh38>" for c in GNOMAD_CHROMOSOMES]
    for chrom in GNOMAD_CHROMOSOMES:
        lines = ["##fileformat=VCFv4.2", '##FILTER=<ID=PASS,Description="All filters passed">'] + info + contigs
        lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
        for pos, ref, alt in sorted(records[chrom]):
            values = ";".join(f"{name}={'0.01' if name.startswith('AF') else 3}" for name in
                              ["AC_joint", "AN_joint", "AF_joint", "nhomalt_joint"])
            lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{values}")
        vcf = net / f"gnomad.joint.v4.1.sites.{chrom}.vcf.bgz"
        plain = out / f"{chrom}.vcf"
        plain.write_text("\n".join(lines) + "\n")
        with open(vcf, "wb") as handle:
            subprocess.run(["bgzip", "-c", plain], stdout=handle, check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", vcf], check=True)
        plain.unlink()
    # The GRCh37 reference: real lengths, bases only around the lifted sites.
    folder = out / "A1"
    fasta = folder / "data/ucsc/hg19.fa"
    fasta.parent.mkdir(parents=True, exist_ok=True)
    index, offset = [], 0
    with open(fasta, "wb") as handle:
        for chrom in GNOMAD_CHROMOSOMES:
            length = sizes19[chrom]
            name = f">{chrom}\n".encode()
            handle.seek(offset)
            handle.write(name)
            start = offset + len(name)
            index.append(f"{chrom}\t{length}\t{start}\t60\t61")
            for spot, base in sorted(bases.get(chrom, {}).items()):
                if 0 <= spot < length:
                    handle.seek(start + (spot // 60) * 61 + spot % 60)
                    handle.write(base.encode())
            offset = start + length + (length + 59) // 60
            handle.seek(offset - 1)
            handle.write(b"\n")
    (folder / "data/ucsc/hg19.fa.fai").write_text("\n".join(index) + "\n")
    small = out / "A1-hg38" / "data/ucsc/hg38.fa"
    small.parent.mkdir(parents=True, exist_ok=True)
    small.write_text(">chr1\nACGTACGTAC\n")
    (small.parent / "hg38.fa.fai").write_text("chr1\t10\t6\t10\t11\n")
    receipts = out / "receipts"
    receipts.mkdir(parents=True, exist_ok=True)
    for label, where, path in (("hg19", folder, fasta), ("hg38", small.parents[2], small)):
        (receipts / f"A1.{label}.json").write_text(json.dumps(
            receipt("A1", "REFERENCE_FASTA", where, {"ref_genome": path}), indent=2) + "\n")
    counts = {c: len(r) for c, r in records.items()}
    print(json.dumps(counts))
    print(out / "net")


def write_clinvar_release(out: Path):
    root = Path(__file__).resolve().parents[1]
    resources = json.loads((root / "cache_resources.json").read_text())["resources"]
    template = resources["clinvar_vcf"]["source"]["latest_url_template"]
    summary_url = resources["clinvar_variant_summary"]["source"]["url"]
    mane_url = resources["mane_select"]["source"]["url"]
    net = out / "net"
    served = lambda url: net / url.split("://", 1)[1]

    def publish(url: str, data: bytes, manifest: bool):
        path = served(url)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        if manifest:
            Path(f"{path}.md5").write_text(f"{hashlib.md5(data).hexdigest()}  {path.name}\n")
        return path

    for assembly, contig_length in (("GRCh37", 249250621), ("GRCh38", 248956422)):
        lines = ["##fileformat=VCFv4.1", "##fileDate=2026-11-16", "##source=ClinVar",
                 '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">',
                 '##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Gene">',
                 f"##contig=<ID=1,length={contig_length}>",
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
        for number, pos in enumerate((925952, 930165, 935222), start=1):
            lines.append(f"1\t{pos}\t{number}\tG\tA\t.\t.\tCLNSIG=Pathogenic;GENEINFO=SAMD11:148398")
        plain = out / f"clinvar.{assembly}.vcf"
        out.mkdir(parents=True, exist_ok=True)
        plain.write_text("\n".join(lines) + "\n")
        vcf = subprocess.run(["bgzip", "-c", plain], capture_output=True, check=True).stdout
        plain.unlink()
        path = publish(template.format(assembly=assembly), vcf, manifest=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", path], check=True)
    summary = table(["#AlleleID", "Type", "Name", "GeneID", "GeneSymbol", "ClinicalSignificance", "Assembly"],
                    [["1", "single nucleotide variant", "NM_152486.4(SAMD11):c.1A>G", "148398", "SAMD11",
                      "Pathogenic", "GRCh38"]])
    publish(summary_url, gzip.compress(summary.encode(), mtime=0), manifest=True)
    mane = table(["#NCBI_GeneID", "Ensembl_Gene", "HGNC_ID", "symbol", "name", "RefSeq_nuc", "RefSeq_prot",
                  "Ensembl_nuc", "Ensembl_prot", "MANE_status", "GRCh38_chr", "chr_start", "chr_end", "chr_strand"],
                 [["GeneID:148398", "ENSG00000187634.13", "HGNC:28706", "SAMD11", "sterile alpha motif domain containing 11",
                   "NM_001385641.1", "NP_001372570.1", "ENST00000616016.5", "ENSP00000478421.2", "MANE Select",
                   "NC_000001.11", "925731", "944574", "+"]])
    publish(mane_url, gzip.compress(mane.encode(), mtime=0), manifest=False)
    print(net)


VEP_INSTALL = r"""#!/usr/bin/env python3
# A stand-in for Ensembl's INSTALL.pl (vep_install): the cache and plugins of
# --AUTO cp, small and the same on every run.
import gzip, os, shutil, sys
from pathlib import Path
args = sys.argv[1:]
here = Path(__file__).resolve().parent
with open(here.parents[2] / "vep_install.log", "a") as log:
    log.write(" ".join(a if not a.startswith("/") else "<path>" for a in args) + "\n")
def value(name):
    return args[args.index(name) + 1] if name in args else None
version, assembly = value("--CACHE_VERSION"), value("--ASSEMBLY")
if os.environ.get("FAKE_VEP_INSTALL_FAIL"):
    sys.exit("stand-in vep_install failure")
if "c" in value("--AUTO"):
    cache = Path(value("--CACHEDIR")) / "homo_sapiens_merged" / f"{version}_{assembly}"
    (cache / "1").mkdir(parents=True, exist_ok=True)
    (cache / "info.txt").write_text(f"species\thomo_sapiens\nassembly\t{os.environ.get('FAKE_VEP_ASSEMBLY', assembly)}\n"
                                    f"source_assembly\t{assembly}\nvar_type\ttabix\n")
    (cache / "chr_synonyms.txt").write_text("1\tchr1\nNC_000001.11\t1\n")
    with gzip.GzipFile(cache / "1" / "all_vars.gz", "wb", mtime=0) as handle:
        handle.write(b"1\t925952\trs1\tG/A\n")
    (cache / "1" / "all_vars.gz.csi").write_bytes(b"CSI\x01stand-in")
if "p" in value("--AUTO"):
    folder = Path(value("--PLUGINSDIR"))
    folder.mkdir(parents=True, exist_ok=True)
    for name in value("--PLUGINS").split(","):
        target = folder / f"{name}.pm"
        if name == "UTRAnnotator":
            shutil.copyfile(here / "UTRAnnotator116.original.pm", target)
        else:
            target.write_text(f"package {name};\n# stand-in plugin of VEP release {version}\n1;\n")
    (folder / "plugin_config.txt").write_text(f"# stand-in plugin_config.txt of release {version}\n")
"""


VEP = r"""#!/usr/bin/env python3
# A stand-in for Ensembl's vep 116: its --help, and an annotation of a VCF
# with the cache's transcripts (stand-in-transcripts.json beside it), the same
# on every run: a ##VEP line naming the cache and the assembly, a CSQ field
# per record (missense on a known transcript, else intergenic) and the
# statistics file. FAKE_VEP_FAIL makes it fail; FAKE_VEP_DROP makes it lose
# the first record.
import gzip, json, os, sys
from pathlib import Path
args = sys.argv[1:]
here = Path(__file__).resolve().parent
if not args or "--help" in args:
    print("#----------------------------------#\n# ENSEMBL VARIANT EFFECT PREDICTOR #\n"
          "#----------------------------------#\n\nVersions:\n  ensembl              : 116.0\n"
          "  ensembl-vep          : 116.0\n")
    sys.exit(0)
with open(here.parents[2] / "vep.log", "a") as log:
    log.write(" ".join(a if not a.startswith("/") else "<path>" for a in args) + "\n")
def value(name):
    return args[args.index(name) + 1] if name in args else None
if os.environ.get("FAKE_VEP_FAIL"):
    sys.exit("stand-in vep failure")
assembly = value("--assembly")
patch = {"GRCh38": "GRCh38.p14", "GRCh37": "GRCh37.p13"}[assembly]
cache = f"{value('--dir_cache')}/homo_sapiens_merged/116_{assembly}"
plugins = [args[i + 1] for i, a in enumerate(args) if a == "-plugin"]
conservation = any(p.startswith("Conservation,") for p in plugins)
transcripts = json.loads((here / "stand-in-transcripts.json").read_text())
fields = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE",
          "EXON", "INTRON", "HGVSc", "HGVSp", "Protein_position", "Amino_acids", "DOMAINS", "CANONICAL"]
fields += ["Conservation"] if conservation else []
source, target = value("-i"), value("-o")
opener = gzip.open if source.endswith("gz") else open
dropped = 0
with opener(source, "rt") as inp, open(target, "w") as out:
    for line in inp:
        if line.startswith("#CHROM"):
            out.write(f'##VEP="v116.0" time="2026-11-16 00:00:00" cache="{cache}" ensembl=116.a1b2c3d '
                      f'ensembl-variation=116.a1b2c3d ensembl-io=116.a1b2c3d ensembl-funcgen=116.a1b2c3d '
                      f'1000genomes="phase3" COSMIC="99" ClinVar="202504" dbSNP="156" gencode="GENCODE 48" '
                      f'genebuild="2014-07" gnomADe="v4.1" gnomADg="v4.1" assembly="{patch}"\n')
            out.write('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
                      f'Format: {"|".join(fields)}">\n')
        elif not line.startswith("#"):
            if os.environ.get("FAKE_VEP_DROP") and dropped < 1:
                dropped += 1
                continue
            row = line.rstrip("\n").split("\t")
            info = dict(item.split("=", 1) for item in row[7].split(";") if "=" in item)
            known = transcripts.get(info.get("TRANSCRIPT", "").split(".")[0])
            values = dict(Allele=row[4], Consequence="intergenic_variant", IMPACT="MODIFIER")
            if known and "PVAR" in info:
                residue = int(info["PVAR"][1:-1])
                first, last = known["domain_residues"]
                exons = known["exons"]
                values.update(Consequence="missense_variant", IMPACT="MODERATE", SYMBOL=known["symbol"],
                              Gene=known["gene"], Feature_type="Transcript", Feature=known["transcript"],
                              BIOTYPE="protein_coding", EXON=f"{min(exons, 1 + (residue - 1) * exons // known['length'])}/{exons}",
                              HGVSc=f"{known['transcript']}.1:c.{3 * residue}{row[3]}>{row[4]}",
                              HGVSp=f"ENSP{known['transcript'][4:]}.1:p.{info['PVAR']}",
                              Protein_position=f"{residue}/{known['length']}",
                              Amino_acids=f"{info['PVAR'][0]}/{info['PVAR'][-1]}",
                              DOMAINS=known["domain"] if first <= residue <= last else "", CANONICAL="YES")
            if conservation:
                values["Conservation"] = f"{(int(row[1]) % 13) / 2:.3f}"
            csq = "CSQ=" + "|".join(str(values.get(name, "")) for name in fields)
            row[7] = csq if row[7] in (".", "") else row[7] + ";" + csq
            line = "\t".join(row) + "\n"
        out.write(line)
with open(value("--stats_file"), "w") as out:
    out.write("<html>stand-in VEP statistics</html>\n")
"""


def write_vep(out: Path):
    root = Path(__file__).resolve().parents[1]
    env = out / "envs" / "priva_acmg"
    tools = env / "bin"
    site = env / "lib" / "perl5" / "site_perl"
    site.mkdir(parents=True, exist_ok=True)
    (tools / "modules" / "Bio" / "EnsEMBL" / "VEP").mkdir(parents=True, exist_ok=True)
    (tools / "modules" / "Bio" / "EnsEMBL" / "ApiVersion.pm").write_text(
        "package Bio::EnsEMBL::ApiVersion;\nsub software_version { return 116; }\n1;\n")
    (tools / "modules" / "Bio" / "EnsEMBL" / "VEP" / "Utils.pm").write_text(
        "package Bio::EnsEMBL::VEP::Utils;\nuse Exporter 'import';\nour @EXPORT_OK = ('get_version_data');\n"
        "sub get_version_data { return {'ensembl-vep' => {'release' => 'release/116'}}; }\n1;\n")
    (tools / ".version").write_text("{}\n")
    shutil.copyfile(root / "tests/fixtures/UTRAnnotator116.original.pm", tools / "UTRAnnotator116.original.pm")
    (tools / "modules" / "Bio" / "EnsEMBL" / "VEP" / "Runner.pm").write_text(
        "package Bio::EnsEMBL::VEP::Runner;\n1;\n")
    (tools / "stand-in-transcripts.json").write_text(json.dumps({
        transcript: dict(symbol=symbol, gene=gene, transcript=transcript, length=length, domain_residues=list(span),
                         domain=domain, exons=3)
        for _chrom, _start, symbol, gene, transcript, length, span, domain, _intolerant in AM_GENES}, indent=2) + "\n")
    scripts = {
        "vep": VEP,
        "vep_install": VEP_INSTALL,
        "perl": f"#!/bin/sh\nexec /usr/bin/perl -I{site} \"$@\"\n",
        "conda": "#!/bin/sh\nexit 0\n",
    }
    for name, text in scripts.items():
        (tools / name).write_text(text)
        (tools / name).chmod(0o755)
    # What step A2 records for each assembly, from the stand-in installer.
    receipts = out / "receipts"
    receipts.mkdir(parents=True, exist_ok=True)
    plugins = "SpliceAI,SpliceVault,UTRAnnotator,AlphaMissense,GO,Phenotypes,Conservation,LOEUF,SameCodon,NMD,PrimateAI,MaveDB"
    for assembly, genome in (("hg38", "GRCh38"), ("hg19", "GRCh37")):
        folder = out / "A2" / assembly
        data = folder / "data"
        subprocess.run([sys.executable, tools / "vep_install", "--AUTO", "cp", "--CACHE_VERSION", "116",
                        "--ASSEMBLY", genome, "--CACHEDIR", data / "VEP_cache", "--PLUGINS", plugins,
                        "--PLUGINSDIR", data / "VEP_plugin"], check=True)
        (receipts / f"A2.{assembly}.json").write_text(json.dumps(receipt(
            "A2", "VEP_CACHE", folder, {"vep_cache_dir": data / "VEP_cache", "vep_plugins_dir": data / "VEP_plugin"}),
            indent=2) + "\n")
    print(env)


def write_spliceai(out: Path, plugin=None):
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location("spliceai_cache_input", root / "scripts/spliceai_cache_input.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    v1 = out / "spliceai" / "v1"
    v1.mkdir(parents=True, exist_ok=True)
    header = ["##fileformat=VCFv4.0",
              '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3 variant annotation. '
              'Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">',
              "##contig=<ID=1,length=249250621>", "\t".join(module.VCF_COLUMNS)]
    for kind, rows in (("snv", ["1\t69091\t.\tA\tC\t.\t.\tSpliceAI=C|OR4F5|0.01|0.00|0.00|0.00|-15|42|2|-27"]),
                       ("indel", ["1\t69092\t.\tT\tTA\t.\t.\tSpliceAI=TA|OR4F5|0.00|0.02|0.00|0.00|-3|8|1|-2"])):
        plain = v1 / f"{kind}.vcf"
        plain.write_text("\n".join(header + rows) + "\n")
        with open(v1 / f"{kind}.vcf.gz", "wb") as handle:
            subprocess.run(["bgzip", "-c", plain], stdout=handle, check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", v1 / f"{kind}.vcf.gz"], check=True)
        plain.unlink()
    native = out / "spliceai" / "native"
    native.mkdir(parents=True, exist_ok=True)
    columns = module.native_columns()
    chroms = [f"chr{n}" for n in range(1, 23)] + ["chrX", "chrY"]

    def row(chrom, pos, ref, alt, strand, delta):
        # Top three of each effect: delta (signed), ref and alt scores, then
        # the distance (two distances for junctions); the summary is the
        # largest of the four top0 site deltas.
        values = [chrom, str(pos), ref, alt, strand]
        for effect in ("donor_gain", "donor_loss", "acceptor_gain", "acceptor_loss", "jxn_gain", "jxn_loss"):
            values += [delta, "0.00", "0.00", "0.00", "0.10", "0.10", "0.10", "0.10", "0.10"]
            values += ["-4", "3", "12"] if not effect.startswith("jxn") else ["-4", "3", "12", "5", "-6", "7"]
        return "\t".join(values + [delta]) + "\n"

    # A scored variant inside the gene of each chromosome, and on chr1 one
    # in each of its two genes (ALPHA on +, BETA on -), when the gene-context
    # resources are written (a projection must hold a record).
    rows = {(kind, chrom): [row(chrom, 150, *alleles, "+", "0.25")] for chrom in chroms
            for kind, alleles in (("snv", ("A", "C")), ("indel", ("A", "AT")))}
    rows[("snv", "chr1")].append(row("chr1", 1500, "G", "T", "-", "0.05"))
    rows[("indel", "chr1")].append(row("chr1", 1600, "GC", "G", "-", "0.01"))
    files = []
    for kind in ("snv", "indel"):
        for chrom in chroms:
            if kind == "indel" and chrom == "chrY":
                continue
            shard = native / f"SpliceAI2_GRCh38_{kind}_{chrom}.tsv.gz"
            text = "\t".join(columns) + "\n" + ("".join(rows.get((kind, chrom), [])) if plugin else "")
            data = gzip.compress(text.encode(), mtime=0)
            shard.write_bytes(data)
            files.append({"kind": kind, "chrom": chrom, "path": str(shard), "bytes": len(data),
                          "sha256": hashlib.sha256(data).hexdigest()})
    manifest = {"schema_version": 1, "status": "downloads_complete", "confirmed_by": "fixture",
                "resource": "SpliceAI2", "version": "2.0", "assembly": "GRCh38",
                "source_revision": "cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9", "files": files}
    if plugin:
        def pinned(name, text, **extra):
            path = native / name
            path.write_text(text)
            data = path.read_bytes()
            return {"path": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(), **extra}

        def gene(gene_id, name, chrom, start, end, strand):
            return (f'{chrom}\tHAVANA\tgene\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gene_id}"; '
                    f'gene_type "protein_coding"; gene_name "{name}"; level 2;\n')

        # A scored gene on each chromosome (each is projected on its own),
        # and a second on chr1 on the reverse strand; VEP 116's GENCODE 50
        # names that one anew.
        genes = [("ENSG00000000001", "ALPHA", "chr1", 100, 900, "+"),
                 ("ENSG00000000002", "BETA", "chr1", 1000, 2000, "-")]
        genes += [(f"ENSG{100 + n:011d}", f"GENE{chrom[3:]}", chrom, 100, 900, "+")
                  for n, chrom in enumerate(chroms[1:])]
        renamed = {"BETA": "BETA2"}
        version = "##description: evidence-based annotation of the human genome (GRCh38), version {} (Ensembl {})\n"
        manifest["scored_regions"] = pinned("scored_regions.tsv", "gene_id\tchrom\tgene_start\tgene_end\tstrand\n" + "".join(
            f"{i}.5\t{c}\t{a}\t{b}\t{st}\n" for i, _, c, a, b, st in genes))
        manifest["gencode_gtf"] = pinned("gencode.v48.annotation.gtf", version.format(48, 114) + "".join(
            gene(f"{i}.5", n, c, a, b, st) for i, n, c, a, b, st in genes), release="GENCODE48")
        manifest["vep_gene_annotation"] = pinned("gencode.v50.annotation.gtf", version.format(50, 116) + "".join(
            gene(f"{i}.6", renamed.get(n, n), c, a, b, st) for i, n, c, a, b, st in genes), release="GENCODE50")
        manifest["vep_spliceai_plugin"] = pinned("SpliceAI.pm", Path(plugin).read_text())
    (native / "SpliceAI2_GRCh38.downloads-complete.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(out / "spliceai")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("what", choices=["clinvar", "alphamissense", "mechanism", "nonlof", "ranking", "gnomad", "clinvar_release", "vep", "spliceai", "alphamissense_inputs"])
    parser.add_argument("--out", required=True)
    parser.add_argument("--plugin", help="spliceai: Ensembl's SpliceAI.pm of VEP 116 (VEP_plugins release/116)")
    args = parser.parse_args()
    out = Path(args.out).resolve()
    if args.what == "spliceai":
        return write_spliceai(out, args.plugin)
    {"clinvar": write_clinvar, "alphamissense": write_alphamissense,
     "mechanism": write_mechanism, "nonlof": write_nonlof,
     "ranking": write_ranking, "gnomad": write_gnomad,
     "clinvar_release": write_clinvar_release, "vep": write_vep,
     "alphamissense_inputs": write_alphamissense_inputs}[args.what](out)


if __name__ == "__main__":
    main()
