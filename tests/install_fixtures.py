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


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("what", choices=["clinvar", "alphamissense", "mechanism", "nonlof", "ranking", "gnomad", "clinvar_release"])
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = Path(args.out).resolve()
    {"clinvar": write_clinvar, "alphamissense": write_alphamissense,
     "mechanism": write_mechanism, "nonlof": write_nonlof,
     "ranking": write_ranking, "gnomad": write_gnomad,
     "clinvar_release": write_clinvar_release}[args.what](out)


if __name__ == "__main__":
    main()
