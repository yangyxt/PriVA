#!/usr/bin/env python3
"""Runs the annotation steps both ways on synthetic data and compares them.

Old: the functions of scripts/annotation_vcf.sh the annotation tasks ran up to
phase 4, through the task's step wrapper as it ran them
(tests/fixtures/annotation-steps-cd91ec2/annotation_step.sh). New:
`priva task step`, the Rust steps (rust/src/pipeline_steps) `priva task run`
calls now.

The data are made up: a small reference named as hg38, a cohort VCF of three
samples (in UCSC and in NCBI names, with records off the primary contigs,
multi-allelic and duplicate records, `*` alleles, a wrong REF and records
without an alternative genotype), per-chromosome gnomAD files, a ClinVar file,
control VCFs and PED files. The gnomAD, ClinVar, control, VEP and CADD steps
start from the preprocess step's output of the old route.

bcftools and tabix are the real programs behind a logging wrapper; `vep`,
CADD.sh, its Conda and every priva command but `variant-ids` (the real one)
are stand-ins that log their arguments and write the files the program
writes. So a run compares the command lines the steps give their tools, and
what bcftools and tabix make of them, not VEP's or CADD's science.

For each case both routes start from the same files in their own folder and
must give the same exit status, the same set of tool commands that write a
file (bcftools with an output option, tabix indexing, vep, CADD.sh, priva;
paths and temporary names masked), the same records and header of the step's
VCF (the ##bcftools_ lines left out) and the same bytes of the CADD table.
Files a route leaves behind besides its outputs are listed: the new route
must leave none the old one did not. Exit 1 on any difference.

  python3 tests/annotation_step_compare.py --priva rust/target/release/priva \\
      --work W [--case NAME]...
"""
import argparse
import gzip
import json
import os
import random
import re
import shutil
import stat
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "annotation-steps-cd91ec2" / "annotation_step.sh"
SEP = "\x1f"

FAKE_PRIVA = r'''#!/usr/bin/env python3
import os, sys
args = sys.argv[1:]
if args == ["--version"]:
    print("priva 0.0.0-stand-in")
    sys.exit(0)
with open(os.environ["TOOL_LOG"], "a") as log:
    log.write("\x1f".join(["priva"] + args) + "\n")
if args and args[0] == "variant-ids":
    os.execv(os.environ["REAL_PRIVA"], [os.environ["REAL_PRIVA"]] + args)
if args and args[0] == "spliceai-cutoffs":
    for path in args[2:]:
        if "MISMATCH" in path:
            print("The SpliceAI SNV and indel files are of different versions", file=sys.stderr)
            sys.exit(1)
    print("SpliceAI cutoffs: stand-in")
'''

FAKE_VEP = r'''#!/usr/bin/env python3
import gzip, os, sys
args = sys.argv[1:]
with open(os.environ["TOOL_LOG"], "a") as log:
    log.write("\x1f".join(["vep"] + args) + "\n")
def value(name):
    return args[args.index(name) + 1] if name in args else None
if os.path.exists(os.path.join(value("--dir_cache") or "", "fail")):
    print("stand-in VEP failure", file=sys.stderr)
    sys.exit(3)
source, target = value("-i"), value("-o")
opener = gzip.open if source.endswith("gz") else open
with opener(source, "rt") as inp, open(target, "w") as out:
    for line in inp:
        if line.startswith("#CHROM"):
            out.write('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations '
                      'from Ensembl VEP. Format: Allele|Consequence|SYMBOL">\n')
        elif not line.startswith("#"):
            fields = line.rstrip("\n").split("\t")
            csq = f"CSQ={fields[4]}|missense_variant|GENE{int(fields[1]) % 7}"
            fields[7] = csq if fields[7] == "." else fields[7] + ";" + csq
            line = "\t".join(fields) + "\n"
        out.write(line)
with open(value("--stats_file"), "w") as out:
    out.write("<html>stand-in VEP statistics</html>\n")
'''

FAKE_CADD = r'''#!/usr/bin/env python3
import gzip, os, sys
args = sys.argv[1:]
env = [f"{key}={os.environ.get(key, '<unset>')}" for key in
       ("PATH", "TMPDIR", "PRIVA_CADD_MMSPLICE_EXON_CACHE", "PRIVA_CADD_PRESCORE_BIN")]
with open(os.environ["TOOL_LOG"], "a") as log:
    log.write("\x1f".join(["CADD.sh"] + args + env) + "\n")
output, genome, sites = args[args.index("-o") + 1], args[args.index("-g") + 1], args[-1]
with gzip.open(sites, "rt") as inp, gzip.open(output, "wt") as out:
    out.write(f"## CADD {genome}-v1.7 (c) University of Washington, stand-in\n")
    out.write("#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n")
    for line in inp:
        if not line.startswith("#"):
            chrom, pos, _, ref, alt = line.split("\t")[:5]
            score = (int(pos) % 97) / 10
            out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{score:.3f}\t{score * 3:.3f}\n")
'''

FAKE_CONDA = r'''#!/usr/bin/env bash
case "$1" in
    --version) echo "conda @VERSION@" ;;
    info) echo '{"platform": "linux-64", "conda_prefix": "@PREFIX@"}' ;;
    config) echo '{"get": {}, "rcs": {}, "success": true}' ;;
    *) exit 1 ;;
esac
'''

LOGGING_TOOL = '''#!/usr/bin/env bash
line=@NAME@
for word in "$@"; do line+=$'\\x1f'"$word"; done
printf '%s\\n' "$line" >> "$TOOL_LOG"
exec @REAL@ "$@"
'''

POPULATIONS = ["nfe", "eas", "afr", "amr", "asj", "fin", "sas", "mid", "remaining"]
GNOMAD_FIELDS = (
    ["AC_joint", "AN_joint", "AF_joint", "nhomalt_joint", "AC_grpmax_joint", "AF_grpmax_joint",
     "AN_grpmax_joint", "nhomalt_grpmax_joint", "AC_joint_XX", "AN_joint_XX", "AC_joint_XY",
     "AN_joint_XY", "AF_joint_XX", "AF_joint_XY", "nhomalt_joint_XX", "nhomalt_joint_XY"]
    + [f"AF_joint_{p}_{s}" for p in POPULATIONS for s in ("XX", "XY")]
    + [f"nhomalt_joint_{p}_XX" for p in POPULATIONS])
CONTIG_LENGTHS = {"chr1": 248956422, "chr2": 242193529, "chrX": 156040895,
                  "chrY": 57227415, "chrM": 16569}
NCBI = {"chr1": "1", "chr2": "2", "chrX": "X", "chrY": "Y", "chrM": "MT"}
SAMPLES = ["S1", "S2", "S3"]


def run(command, **kwargs):
    return subprocess.run(command, check=True, **kwargs)


def make_executable(path: Path, text: str):
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def write_vcf(path: Path, header: list, rows: list):
    """A bgzipped, tabix-indexed VCF from header lines and record lines."""
    path.parent.mkdir(parents=True, exist_ok=True)
    plain = path.parent / (path.name + ".plain")
    plain.write_text("\n".join(header + rows) + "\n")
    with open(path, "wb") as out:
        run(["bgzip", "-c", str(plain)], stdout=out)
    plain.unlink()
    run(["tabix", "-f", "-p", "vcf", str(path)])


def records_of(vcf: Path) -> list:
    text = run(["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(vcf)],
               capture_output=True, text=True).stdout
    return [line.split("\t") for line in text.splitlines()]


def build_reference(data: Path, rng: random.Random) -> dict:
    sizes = {"chr1": 4000, "chr2": 2000, "chrX": 2000, "chrY": 1500, "chrM": 600}
    sequences = {name: "".join(rng.choice("ACGT") for _ in range(size))
                 for name, size in sizes.items()}
    fasta = data / "ref" / "genome.hg38.fa"
    fasta.parent.mkdir(parents=True, exist_ok=True)
    index = []
    offset = 0
    with open(fasta, "w") as out:
        for name, sequence in sequences.items():
            header = f">{name}\n"
            out.write(header)
            offset += len(header)
            index.append(f"{name}\t{len(sequence)}\t{offset}\t60\t61")
            for at in range(0, len(sequence), 60):
                line = sequence[at:at + 60] + "\n"
                out.write(line)
                offset += len(line)
    (data / "ref" / "genome.hg38.fa.fai").write_text("\n".join(index) + "\n")
    return {"fasta": str(fasta), "sequences": sequences}


def cohort_rows(sequences: dict, rng: random.Random) -> list:
    """Records of the cohort, in UCSC names, sorted as the header lists them."""
    def genotypes(alts: int) -> str:
        cells = []
        for _ in SAMPLES:
            gt = rng.choice(["0/0", "0/1", "1/1", "./.", "0|1"] + (["1/2"] if alts > 1 else []))
            depth = [rng.randint(0, 30) for _ in range(alts + 1)]
            cells.append(f"{gt}:{','.join(map(str, depth))}")
        return "\t".join(cells)

    rows = []
    for chrom, count in [("chr1", 30), ("chr2", 8), ("chrX", 8), ("chrY", 4), ("chrM", 4)]:
        sequence = sequences[chrom]
        positions = sorted(rng.sample(range(30, len(sequence) - 30), count))
        special = {positions[i]: kind for i, kind in enumerate(
            ["multi", "dup", "star", "del", "ins", "badref", "homref"])} if chrom == "chr1" else {}
        for pos in positions:
            ref = sequence[pos - 1]
            alt = rng.choice([b for b in "ACGT" if b != ref])
            kind = special.get(pos)
            if kind == "multi":
                second = [b for b in "ACGT" if b not in (ref, alt)][0]
                rows.append([chrom, pos, ".", ref, f"{alt},{second}", genotypes(2)])
            elif kind == "star":
                rows.append([chrom, pos, ".", ref, f"{alt},*", genotypes(2)])
            elif kind == "del":
                rows.append([chrom, pos, ".", sequence[pos - 1:pos + 2], ref, genotypes(1)])
            elif kind == "ins":
                rows.append([chrom, pos, ".", ref, ref + "AG", genotypes(1)])
            elif kind == "badref":
                wrong = [b for b in "ACGT" if b not in (ref, alt)][0]
                rows.append([chrom, pos, ".", wrong, alt, genotypes(1)])
            elif kind == "homref":
                rows.append([chrom, pos, ".", ref, alt, "\t".join("0/0:9,0" for _ in SAMPLES)])
            else:
                rows.append([chrom, pos, ".", ref, alt, genotypes(1)])
                if kind == "dup":
                    rows.append(list(rows[-1]))
    rows.append(["chr1_KI270706v1_random", 100, ".", "A", "G", genotypes(1)])
    return rows


def vcf_lines(rows: list, names: dict = None, sites_only: bool = False) -> list:
    out = []
    for chrom, pos, vid, ref, alt, cells in rows:
        chrom = (names or {}).get(chrom, chrom)
        fields = [chrom, str(pos), vid, ref, alt, "50", "PASS", "."]
        if not sites_only:
            fields += ["GT:AD", cells]
        out.append("\t".join(fields))
    return out


def cohort_header(contigs: list, sites_only: bool = False, extra: list = ()) -> list:
    header = ["##fileformat=VCFv4.2"]
    header += contigs
    header += list(extra)
    if not sites_only:
        header += ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                   '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">']
    columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    if not sites_only:
        columns += "\tFORMAT\t" + "\t".join(SAMPLES)
    return header + [columns]


def build_inputs(work: Path, priva: str, stand_in: Path, tools: Path, checkout: Path,
                 path: str) -> dict:
    data = work / "data"
    done = data / "done.json"
    if done.exists():
        return json.loads(done.read_text())
    shutil.rmtree(data, ignore_errors=True)
    data.mkdir(parents=True)
    rng = random.Random(7)
    reference = build_reference(data, rng)
    rows = cohort_rows(reference["sequences"], rng)
    random_contig = "##contig=<ID=chr1_KI270706v1_random,length=175055,assembly=hg38>"

    def contigs(assembly: str, names: dict = None, attribute: bool = True) -> list:
        lengths = dict(CONTIG_LENGTHS)
        if assembly == "hg19":
            lengths["chr1"] = 249250621
        lines = []
        for name, length in lengths.items():
            tail = f",assembly={assembly}" if attribute else ""
            lines.append(f"##contig=<ID={(names or {}).get(name, name)},length={length}{tail}>")
        return lines

    inputs = {"reference": reference["fasta"]}
    ucsc = data / "cohort.vcf.gz"
    write_vcf(ucsc, cohort_header(contigs("hg38") + [random_contig]), vcf_lines(rows))
    inputs["ucsc"] = str(ucsc)
    ncbi = data / "cohort_ncbi.vcf.gz"
    ncbi_rows = [r for r in rows if r[0] in NCBI]
    write_vcf(ncbi, cohort_header(contigs("hg38", NCBI)), vcf_lines(ncbi_rows, NCBI))
    inputs["ncbi"] = str(ncbi)
    # bcftools-style contig lines (the length ends the line) and an hg38-only
    # contig: the assembly comes from the contig.
    sites = data / "sites.vcf.gz"
    alt_contig = "##contig=<ID=chr6_GL000256v2_alt,length=4929269>"
    write_vcf(sites, cohort_header(contigs("hg38", attribute=False) + [alt_contig], True),
              vcf_lines([r for r in rows if r[0] in NCBI], sites_only=True))
    inputs["sites"] = str(sites)
    hg19 = data / "cohort_hg19.vcf.gz"
    write_vcf(hg19, cohort_header(contigs("hg19")), vcf_lines(ncbi_rows))
    inputs["hg19"] = str(hg19)
    unknown = data / "cohort_unknown.vcf.gz"
    write_vcf(unknown, cohort_header(contigs("hg38", attribute=False)), vcf_lines(ncbi_rows))
    inputs["unknown"] = str(unknown)

    # The preprocessed cohort, as the old route makes it, for the other steps.
    route = data / "prepare"
    route.mkdir()
    (route / "tmp").mkdir()
    config = route / "task.yaml"
    config.write_text(yaml_text({"ref_genome": reference["fasta"], "threads": 2,
                                 "tmp_dir": str(route / "tmp")}))
    staged = route / "cohort.vcf.gz"
    shutil.copyfile(ucsc, staged)
    shutil.copyfile(str(ucsc) + ".tbi", str(staged) + ".tbi")
    anno = data / "cohort.anno.vcf.gz"
    env = route_env(route, stand_in, priva, tools, path)
    with open(route / "stderr.txt", "w") as err:
        run(["bash", str(FIXTURE), str(checkout), str(config), str(staged), "preprocess", str(anno)],
            env=env, stdout=err, stderr=err)
    inputs["anno"] = str(anno)
    records = records_of(anno)

    empty = data / "cohort.empty.anno.vcf.gz"
    run(["bcftools", "view", "-h", "-Oz", "-o", str(empty), str(anno)])
    run(["tabix", "-f", "-p", "vcf", str(empty)])
    inputs["empty"] = str(empty)
    chrm = data / "cohort.chrM.anno.vcf.gz"
    run(["bcftools", "view", "-r", "chrM", "-Oz", "-o", str(chrm), str(anno)])
    run(["tabix", "-f", "-p", "vcf", str(chrm)])
    inputs["chrM"] = str(chrm)
    plain = data / "cohort.plain.anno.vcf"
    run(["bcftools", "view", "-Ov", "-o", str(plain), str(anno)])
    inputs["plain"] = str(plain)

    # gnomAD: one file per chromosome, most cohort records and some others.
    gnomad = data / "gnomad"
    header_fields = []
    for field in GNOMAD_FIELDS:
        kind = "Float" if field.startswith("AF") else "Integer"
        number = "1" if field.startswith("AN") else "A"
        header_fields.append(f'##INFO=<ID={field},Number={number},Type={kind},Description="{field}">')
    header_fields.append('##INFO=<ID=vep,Number=.,Type=String,Description="not copied">')
    for chrom in ["chr1", "chr2", "chrX", "chrY"]:
        mine = sorted({(int(p), r, a) for c, p, r, a in records if c == chrom}
                      | {(rng.randint(30, 1400), "A", "C") for _ in range(4)})
        lines = []
        for pos, ref, alt in mine:
            if rng.random() < 0.25:
                continue
            values = []
            for field in GNOMAD_FIELDS:
                if field.startswith("AF"):
                    values.append(f"{field}={rng.random() / 10:.5f}")
                else:
                    values.append(f"{field}={rng.randint(0, 5000)}")
            values.append("vep=x")
            lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{';'.join(values)}")
        header = (["##fileformat=VCFv4.2", f"##contig=<ID={chrom},length={CONTIG_LENGTHS[chrom]}>"]
                  + header_fields + ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"])
        write_vcf(gnomad / f"gnomad.joint.v4.1.sites.{chrom}.vcf.bgz", header, lines)
    inputs["gnomad_x"] = str(gnomad / "gnomad.joint.v4.1.sites.chrX.vcf.bgz")
    # A copy without chrY, for a missing gnomAD file.
    partial = data / "gnomad_partial"
    partial.mkdir()
    for path in gnomad.iterdir():
        if ".chrY." not in path.name:
            shutil.copyfile(path, partial / path.name)
    inputs["gnomad_partial_x"] = str(partial / "gnomad.joint.v4.1.sites.chrX.vcf.bgz")
    already = data / "cohort.gnomad_tagged.anno.vcf.gz"
    header_file = data / "gnomad.hdr"
    header_file.write_text("\n".join(header_fields[:-1]) + "\n")
    run(["bcftools", "annotate", "-h", str(header_file), "-Oz", "-o", str(already), str(anno)])
    run(["tabix", "-f", "-p", "vcf", str(already)])
    inputs["gnomad_tagged"] = str(already)

    # ClinVar, in UCSC names with VEP's CSQ as the installer leaves it.
    clinvar_lines = []
    for chrom, pos, ref, alt in records + [["chr1", "3990", "A", "T"]]:
        if rng.random() < 0.4:
            continue
        info = (f"CLNDN=Disease_{pos};CLNHGVS=NC_000001.11:g.{pos}{ref}>{alt};"
                f"CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Pathogenic;"
                f"GENEINFO=GENE{int(pos) % 5}:{pos};CSQ={alt}|missense_variant|GENE{int(pos) % 5};"
                f"ALLELEID={pos}")
        clinvar_lines.append((chrom, int(pos), f"{chrom}\t{pos}\t{pos}\t{ref}\t{alt}\t.\t.\t{info}"))
    order = list(CONTIG_LENGTHS)
    clinvar_lines.sort(key=lambda item: (order.index(item[0]), item[1]))
    clinvar_header = ["##fileformat=VCFv4.1"] + [
        f"##contig=<ID={c},length={n}>" for c, n in CONTIG_LENGTHS.items()] + [
        f'##INFO=<ID={tag},Number=.,Type=String,Description="{tag}">'
        for tag in ["CLNDN", "CLNHGVS", "CLNREVSTAT", "CLNSIG", "GENEINFO", "CSQ"]] + [
        '##INFO=<ID=ALLELEID,Number=1,Type=Integer,Description="not copied">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    clinvar = data / "clinvar" / "clinvar.vcf.gz"
    write_vcf(clinvar, clinvar_header, [line for _, _, line in clinvar_lines])
    inputs["clinvar"] = str(clinvar)
    tagged = data / "cohort.clinvar_tagged.anno.vcf.gz"
    run(["bcftools", "annotate", "-a", str(clinvar), "-c",
         ".INFO/CLNDN,.INFO/CLNHGVS,.INFO/CLNREVSTAT,.INFO/CLNSIG,.INFO/GENEINFO,.INFO/CLNCSQ:=INFO/CSQ",
         "-Oz", "-o", str(tagged), str(anno)])
    run(["tabix", "-f", "-p", "vcf", str(tagged)])
    inputs["clinvar_tagged"] = str(tagged)

    # Control cohorts: one with its counts, one with genotypes only.
    control_header = ["##fileformat=VCFv4.2"] + [
        f"##contig=<ID={c},length={n}>" for c, n in CONTIG_LENGTHS.items()]
    counted, genotyped = [], []
    for chrom, pos, ref, alt in records:
        if rng.random() < 0.3:
            continue
        an = rng.randint(100, 400)
        ac = rng.randint(0, an // 2)
        counted.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t"
                       f"AC={ac};AN={an};AF={ac / an:.6f};AC_Hom={ac // 4}")
        cells = "\t".join(rng.choice(["0/0", "0/1", "1/1", "./."]) for _ in range(4))
        genotyped.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{cells}")
    write_vcf(data / "control" / "control_counts.vcf.gz", control_header + [
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="AC">',
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="AN">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">',
        '##INFO=<ID=AC_Hom,Number=A,Type=Integer,Description="AC_Hom">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"], counted)
    write_vcf(data / "control" / "control_genotypes.vcf.gz", control_header + [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tC1\tC2\tC3\tC4"], genotyped)
    # The genotype-only control comes without its index: the step makes one.
    os.unlink(data / "control" / "control_genotypes.vcf.gz.tbi")
    inputs["control_counts"] = str(data / "control" / "control_counts.vcf.gz")
    inputs["control_genotypes"] = str(data / "control" / "control_genotypes.vcf.gz")
    tagged = data / "cohort.control_tagged.anno.vcf.gz"
    run(["bcftools", "annotate", "-a", inputs["control_counts"], "-c",
         "INFO/control_AC:=INFO/AC,INFO/control_AN:=INFO/AN,INFO/control_AF:=INFO/AF,"
         "INFO/control_nhomalt:=INFO/AC_Hom", "-Oz", "-o", str(tagged), str(anno)])
    run(["tabix", "-f", "-p", "vcf", str(tagged)])
    inputs["control_tagged"] = str(tagged)
    peds = {
        "ped": "#Family\tIndividual\tFather\tMother\tSex\tPhenotype\n"
               "FAM\tS1\t0\t0\t1\t2\nFAM\tS2\t0\t0\t2\t1\nFAM\tS3\t0\t0\t1\t1\n",
        "ped_affected": "FAM\tS1\t0\t0\t1\t2\nFAM\tS2\t0\t0\t2\t2\nFAM\tS3\t0\t0\t1\t2\n",
        "ped_elsewhere": "FAM\tX1\t0\t0\t1\t1\nFAM\tS1\t0\t0\t1\t2\n",
    }
    for name, text in peds.items():
        path = data / f"{name}.ped"
        path.write_text(text)
        inputs[name] = str(path)
    spliced = data / "cohort.vep_tagged.anno.vcf.gz"
    header_file = data / "csq.hdr"
    header_file.write_text('##INFO=<ID=CSQ,Number=.,Type=String,Description="VEP">\n')
    run(["bcftools", "annotate", "-h", str(header_file), "-Oz", "-o", str(spliced), str(anno)])
    run(["tabix", "-f", "-p", "vcf", str(spliced)])
    inputs["vep_tagged"] = str(spliced)
    done.write_text(json.dumps(inputs))
    return inputs


def vep_resources(work: Path) -> dict:
    folder = work / "vep"
    config = {}
    for key, kind in [("vep_cache_dir", "dir"), ("vep_plugins_dir", "dir"),
                      ("vep_plugins_cachedir", "dir"), ("loftee_repo", "dir"),
                      ("utr_annotator_file", "file"), ("alphamissense_prescore", "file"),
                      ("spliceai_snv_prescore", "file"), ("spliceai_indel_prescore", "file"),
                      ("primateai_prescore", "file"), ("conservation_file", "file"),
                      ("splicevault_prescore", "file"), ("human_ancestor_fasta", "file"),
                      ("loftee_conservation_file", "file"), ("gerp_bigwig", "file"),
                      ("mavedb_file", "file"), ("pext_mean_bw", "file")]:
        path = folder / key
        if kind == "dir":
            path.mkdir(parents=True, exist_ok=True)
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("stand-in\n")
        config[key] = str(path)
    tissues = folder / "pext" / "tissues"
    tissues.mkdir(parents=True, exist_ok=True)
    for tissue in ["Brain_Cortex", "Cells-EBV-transformed_lymphocytes"]:
        (tissues / f"{tissue}.bw").write_text("stand-in\n")
    config["pext_tissue_bw_dir"] = str(folder / "pext")
    config["loeuf_prescore"] = "/downloads/loeuf_dataset_hg38.tsv.gz"
    failing = folder / "failing_cache"
    failing.mkdir(exist_ok=True)
    (failing / "fail").write_text("")
    mismatch = folder / "spliceai_MISMATCH.vcf.gz"
    mismatch.write_text("stand-in\n")
    return {"config": config, "failing_cache": str(failing), "mismatch": str(mismatch)}


def cadd_setups(work: Path) -> dict:
    """CADD.sh and its Conda frontend: a current Conda, an old one, none, and
    a current one without CADD.sh."""
    out = {}
    for name, version, script in [("ok", "24.11.3", True), ("old", "23.10.0", True),
                                  ("none", None, True), ("noscript", "24.11.3", False)]:
        top = work / "cadd" / name
        scripts = top / "CADD-scripts"
        scripts.mkdir(parents=True, exist_ok=True)
        if script:
            make_executable(scripts / "CADD.sh", FAKE_CADD)
        if version:
            frontend = top / "conda_frontend"
            (frontend / "bin").mkdir(parents=True, exist_ok=True)
            (frontend / "bin" / "activate").write_text("")
            text = FAKE_CONDA.replace("@VERSION@", version).replace("@PREFIX@", str(frontend))
            make_executable(frontend / "bin" / "conda", text)
        out[name] = str(scripts)
    # Only GRCh38 has PriVA's MMSplice exon tables here.
    (Path(out["ok"]) / "data" / "annotations" / "GRCh38_v1.7" / "mmsplice").mkdir(
        parents=True, exist_ok=True)
    return out


def checkout(work: Path) -> Path:
    """A checkout with copied scripts and the repository's data/ folder."""
    root = work / "checkout"
    if not root.exists():
        root.mkdir(parents=True)
        shutil.copytree(ROOT / "scripts", root / "scripts")
        (root / "data").symlink_to(ROOT / "data")
    return root


def logging_tools(work: Path, original_path: str) -> Path:
    """bcftools and tabix that log their arguments, and the stand-in vep."""
    tools = work / "bin"
    tools.mkdir(parents=True, exist_ok=True)
    for name in ["bcftools", "tabix"]:
        real = shutil.which(name, path=original_path)
        if not real:
            sys.exit(f"{name} is needed")
        make_executable(tools / name, LOGGING_TOOL.replace("@NAME@", name).replace("@REAL@", real))
    make_executable(tools / "vep", FAKE_VEP)
    return tools


def clean_path(original: str) -> str:
    """PATH without folders that hold a conda program, so only the stand-in
    frontends offer one."""
    keep = [part for part in original.split(os.pathsep)
            if part and not os.path.exists(os.path.join(part, "conda"))]
    return os.pathsep.join(keep)


def route_env(route: Path, stand_in: Path, real_priva: str, tools: Path, path: str) -> dict:
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("SLURM_", "CONDA", "PRIVA_", "_CE_"))}
    env.update(PATH=str(tools) + os.pathsep + path, PRIVA_BINARY=str(stand_in),
               REAL_PRIVA=real_priva, TOOL_LOG=str(route / "tool.log"),
               TMPDIR=str(route / "tmp"), PYTHONDONTWRITEBYTECODE="1")
    return env


def yaml_text(config: dict) -> str:
    lines = []
    for key, value in config.items():
        if value is None:
            lines.append(f"{key}: null")
        elif isinstance(value, bool):
            lines.append(f"{key}: {'true' if value else 'false'}")
        elif isinstance(value, (int, float)):
            lines.append(f"{key}: {value}")
        else:
            lines.append(f"{key}: {json.dumps(value)}")
    return "\n".join(lines) + "\n"


def base_config(inputs: dict, vep: dict, cadd: dict) -> dict:
    config = {
        "ref_genome": inputs["reference"], "threads": 2, "assembly": "hg38",
        "gnomad_vcf_chrX": inputs["gnomad_x"], "clinvar_vcf": inputs["clinvar"],
        "control_vcf": inputs["control_counts"], "ped_file": inputs["ped"],
        "cadd_base_dir": cadd["ok"], "pext_enabled": False,
        "pext_tissues": "Brain_Cortex, Liver , Cells-EBV-transformed_lymphocytes",
    }
    config.update(vep["config"])
    return config


# Each case: the step, its input VCF, and the configuration changes
# ({route} is the case's folder); `copy` stages files into the route first.
def cases(inputs: dict, vep: dict, cadd: dict) -> dict:
    anno = inputs["anno"]
    return {
        "preprocess-ucsc": dict(step="preprocess", vcf=inputs["ucsc"]),
        "preprocess-ncbi": dict(step="preprocess", vcf=inputs["ncbi"], config={"threads": 1}),
        "preprocess-sites-only": dict(step="preprocess", vcf=inputs["sites"]),
        "preprocess-assembly-unknown": dict(step="preprocess", vcf=inputs["unknown"],
                                            config={"threads": None}),
        "preprocess-assembly-mismatch": dict(step="preprocess", vcf=inputs["hg19"]),
        "preprocess-missing-reference": dict(step="preprocess", vcf=inputs["ucsc"],
                                             config={"ref_genome": "/missing/genome.fa"}),
        "gnomad": dict(step="gnomad", vcf=anno),
        "gnomad-one-thread": dict(step="gnomad", vcf=anno, config={"threads": 1}),
        "gnomad-already-tagged": dict(step="gnomad", vcf=inputs["gnomad_tagged"]),
        "gnomad-chrM-only": dict(step="gnomad", vcf=inputs["chrM"]),
        "gnomad-no-records": dict(step="gnomad", vcf=inputs["empty"]),
        "gnomad-missing-chromosome-file": dict(step="gnomad", vcf=anno, config={
            "gnomad_vcf_chrX": inputs["gnomad_partial_x"]}),
        "gnomad-plain-vcf": dict(step="gnomad", vcf=inputs["plain"], name="TEST.anno.vcf"),
        "clinvar": dict(step="clinvar", vcf=anno),
        "clinvar-no-threads": dict(step="clinvar", vcf=anno, config={"threads": None}),
        "clinvar-already-tagged": dict(step="clinvar", vcf=inputs["clinvar_tagged"]),
        "clinvar-missing-file": dict(step="clinvar", vcf=anno,
                                     config={"clinvar_vcf": "/missing/clinvar.vcf.gz"}),
        "control-counts": dict(step="control", vcf=anno, copy={"control.vcf.gz": inputs["control_counts"]},
                               config={"control_vcf": "{route}/control.vcf.gz"}),
        "control-genotypes": dict(step="control", vcf=anno,
                                  copy={"control.vcf.gz": inputs["control_genotypes"]},
                                  config={"control_vcf": "{route}/control.vcf.gz"}),
        "control-from-cohort": dict(step="control", vcf=anno, config={"control_vcf": None}),
        "control-none-optional": dict(step="control", vcf=anno, config={
            "control_vcf": None, "ped_file": inputs["ped_affected"], "control_vcf_optional": True}),
        "control-none-required": dict(step="control", vcf=anno, config={
            "control_vcf": None, "ped_file": inputs["ped_affected"]}),
        "control-ped-elsewhere": dict(step="control", vcf=anno, config={
            "control_vcf": "", "ped_file": inputs["ped_elsewhere"], "control_vcf_optional": True}),
        "control-already-tagged": dict(step="control", vcf=inputs["control_tagged"]),
        "vep-hg38": dict(step="vep", vcf=anno),
        "vep-hg19": dict(step="vep", vcf=anno, config={"assembly": "hg19", "threads": 3}),
        "vep-pext": dict(step="vep", vcf=anno, config={"pext_enabled": True}),
        "vep-pext-missing-mean": dict(step="vep", vcf=anno, config={
            "pext_enabled": True, "pext_mean_bw": "/missing/pext.bw", "pext_tissues": None}),
        "vep-assembly-from-vcf": dict(step="vep", vcf=anno, config={"assembly": ""}),
        "vep-no-threads": dict(step="vep", vcf=anno, config={"threads": None}),
        "vep-threads-text": dict(step="vep", vcf=anno, config={"threads": " +04 "}),
        "vep-missing-plugin-file": dict(step="vep", vcf=anno,
                                        config={"primateai_prescore": "/missing/primateai.tsv.gz"}),
        "vep-missing-mavedb": dict(step="vep", vcf=anno, config={"mavedb_file": None}),
        "vep-spliceai-mismatch": dict(step="vep", vcf=anno,
                                      config={"spliceai_indel_prescore": vep["mismatch"]}),
        "vep-fails": dict(step="vep", vcf=anno, config={"vep_cache_dir": vep["failing_cache"]}),
        "vep-already-annotated": dict(step="vep", vcf=inputs["vep_tagged"]),
        "cadd-hg38": dict(step="cadd", vcf=anno),
        "cadd-hg19": dict(step="cadd", vcf=anno, config={"assembly": "GRCh37"}),
        "cadd-no-records": dict(step="cadd", vcf=inputs["empty"]),
        "cadd-old-conda": dict(step="cadd", vcf=anno, config={"cadd_base_dir": cadd["old"]}),
        "cadd-no-conda": dict(step="cadd", vcf=anno, config={"cadd_base_dir": cadd["none"]}),
        "cadd-no-script": dict(step="cadd", vcf=anno, config={"cadd_base_dir": cadd["noscript"]}),
        "cadd-unknown-assembly": dict(step="cadd", vcf=anno, config={"assembly": "T2T"}),
    }


def stage(route: Path, case: dict, config: dict) -> dict:
    """The route's folder: the input VCF (in input/ for preprocess, in
    result/ for the steps that change it in place), the configuration."""
    shutil.rmtree(route, ignore_errors=True)
    result = route / "result"
    result.mkdir(parents=True)
    (route / "tmp").mkdir()
    if case["step"] == "preprocess":
        (route / "input").mkdir()
        vcf = route / "input" / "TEST.vcf.gz"
        output = result / "TEST.anno.vcf.gz"
    else:
        vcf = result / case.get("name", "TEST.anno.vcf.gz")
        output = result / "TEST.anno.cadd.tsv" if case["step"] == "cadd" else vcf
    shutil.copyfile(case["vcf"], vcf)
    if os.path.isfile(case["vcf"] + ".tbi"):
        shutil.copyfile(case["vcf"] + ".tbi", str(vcf) + ".tbi")
    for name, source in case.get("copy", {}).items():
        shutil.copyfile(source, route / name)
        if os.path.isfile(source + ".tbi"):
            shutil.copyfile(source + ".tbi", str(route / name) + ".tbi")
    config = dict(config, **case.get("config", {}))
    config = {k: v.replace("{route}", str(route)) if isinstance(v, str) else v
              for k, v in config.items()}
    config.update(output_dir=str(result), input_vcf=str(vcf), tmp_dir=str(route / "tmp"))
    path = route / "task.yaml"
    path.write_text(yaml_text(config))
    return {"config": path, "vcf": vcf, "output": output}


def old_command(root: Path, case: dict, paths: dict) -> list:
    command = ["bash", str(FIXTURE), str(root), str(paths["config"]), str(paths["vcf"]), case["step"]]
    if case["step"] in ("preprocess", "cadd"):
        command.append(str(paths["output"]))
    return command


def new_command(priva: str, root: Path, case: dict, paths: dict) -> list:
    command = [priva, "task", "step", "--step", case["step"], "--root", str(root),
               "--config", str(paths["config"]), "--vcf", str(paths["vcf"])]
    if case["step"] in ("preprocess", "cadd"):
        command += ["--output", str(paths["output"])]
    return command


TAG = re.compile(r"(?<=\.)[A-Za-z0-9_]{6,}(?=\.(?:vcf|log))|(?<=\.gnomad-)([a-z]+)\.[A-Za-z0-9_]+")


def mask(text: str, route: Path, work: Path) -> str:
    text = text.replace(str(route), "ROUTE").replace(str(work), "WORK")
    return TAG.sub(lambda m: f"{m.group(1)}.TAG" if m.group(1) else "TAG", text)


def writes(words: list) -> bool:
    """A tool command that writes a file, as opposed to one that reads."""
    tool, rest = words[0], words[1:]
    if tool == "bcftools":
        return rest[:1] == ["index"] or any(
            w == "-o" or w.startswith("-O") or w.startswith("--output") for w in rest)
    if tool == "tabix":
        return "-p" in rest
    return True


def vcf_text(path: Path) -> str:
    text = run(["bcftools", "view", str(path)], capture_output=True, text=True).stdout
    return "\n".join(line for line in text.splitlines() if not line.startswith("##bcftools_"))


def outcome(route: Path, work: Path, status: int, output: Path) -> dict:
    log = route / "tool.log"
    commands = []
    for line in (log.read_text().splitlines() if log.exists() else []):
        words = line.split(SEP)
        if writes(words):
            commands.append(" ".join(mask(w, route, work) for w in words))
    outputs = {}
    declared = {output, Path(str(output) + ".tbi")}
    for path in declared:
        if path.is_file():
            name = str(path.relative_to(route))
            if name.endswith(".tbi"):
                outputs[name] = "index"
            elif name.endswith(".vcf.gz"):
                outputs[name] = mask(vcf_text(path), route, work)
            else:
                outputs[name] = path.read_bytes()
    leftovers = set()
    for path in route.rglob("*"):
        if path.is_file() or path.is_symlink():
            name = str(path.relative_to(route))
            if path not in declared and name not in ("task.yaml", "stderr.txt", "tool.log"):
                leftovers.add(mask(name, route, work))
        elif path.is_dir() and not any(path.iterdir()) and path.name not in ("tmp", "result", "input"):
            leftovers.add(mask(str(path.relative_to(route)), route, work) + "/")
    return {"status": status == 0, "commands": sorted(commands), "outputs": outputs,
            "leftovers": leftovers, "order": commands}


def run_route(command: list, route: Path, env: dict) -> int:
    with open(route / "stderr.txt", "w") as err:
        return subprocess.run(command, cwd=route, env=env, stdout=err, stderr=err).returncode


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--priva", required=True)
    parser.add_argument("--work", required=True)
    parser.add_argument("--case", action="append")
    args = parser.parse_args()
    try:
        import packaging.version  # noqa: F401  (scripts/cadd_conda.py needs it)
    except ImportError:
        sys.exit("The old route's scripts/cadd_conda.py needs the Python package packaging")
    work = Path(args.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    priva = str(Path(args.priva).resolve())
    path = clean_path(os.environ.get("PATH", ""))
    tools = logging_tools(work, path)
    stand_in = work / "priva"
    make_executable(stand_in, FAKE_PRIVA)
    root = checkout(work)
    inputs = build_inputs(work, priva, stand_in, tools, root, path)
    vep = vep_resources(work)
    cadd = cadd_setups(work)
    config = base_config(inputs, vep, cadd)
    all_cases = cases(inputs, vep, cadd)
    failed = 0
    for name in args.case or all_cases:
        case = all_cases[name]
        results = {}
        for route_name in ("old", "new"):
            route = work / "cases" / name / route_name
            paths = stage(route, case, config)
            if route_name == "old":
                command = old_command(root, case, paths)
            else:
                command = new_command(priva, root, case, paths)
            status = run_route(command, route, route_env(route, stand_in, priva, tools, path))
            results[route_name] = outcome(route, work, status, paths["output"])
        old, new = results["old"], results["new"]
        problems = []
        if old["status"] != new["status"]:
            problems.append(f"status: old {old['status']} new {new['status']}")
        if old["commands"] != new["commands"]:
            only_old = [c for c in old["commands"] if c not in new["commands"]]
            only_new = [c for c in new["commands"] if c not in old["commands"]]
            problems.append("commands:\n      old only: " + "\n                ".join(only_old or ["-"])
                            + "\n      new only: " + "\n                ".join(only_new or ["-"]))
        if sorted(old["outputs"]) != sorted(new["outputs"]):
            problems.append(f"outputs: old {sorted(old['outputs'])} new {sorted(new['outputs'])}")
        for file in sorted(set(old["outputs"]) & set(new["outputs"])):
            if old["outputs"][file] != new["outputs"][file]:
                problems.append(f"{file} differs")
        extra = sorted(new["leftovers"] - old["leftovers"])
        if extra:
            problems.append(f"left behind by the new route only: {extra}")
        state = "ok" if new["status"] else "refused"
        if problems:
            failed += 1
            print(f"DIFF {name} ({state})")
            for problem in problems:
                print("  " + problem)
        else:
            records = ""
            for file, text in new["outputs"].items():
                if file.endswith(".vcf.gz"):
                    count = sum(1 for line in text.splitlines() if not line.startswith("#"))
                    records = f"; {count} records"
            print(f"SAME {name} ({state}; {len(new['commands'])} tool commands{records})")
        gone = sorted(old["leftovers"] - new["leftovers"])
        if gone:
            print(f"  left behind by the old route only: {gone}")
    print(f"{failed} of {len(args.case or all_cases)} cases differ")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
