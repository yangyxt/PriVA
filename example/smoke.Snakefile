# Prepared-input smoke DAG. Raw annotation/scoring qualification is never bypassed:
# these input fixtures already contain annotations and transcript-aware PHRED.
from pathlib import Path
import shutil

if config.get("smoke_fixture_only") is not True:
    raise ValueError("This workflow accepts only explicitly selected smoke fixtures.")

ROOT = Path(workflow.basedir).parent
BASE = Path(config["input_vcf"]).name.removesuffix(".vcf.gz")
OUT = str(Path(config["output_dir"]) / (BASE + ".anno"))

rule all:
    input:
        OUT + ".filtered.tsv",
        OUT + ".filtered.acmg.tsv"

rule stage_prepared_fixtures:
    input:
        vcf=config["input_vcf"],
        index=config["input_vcf"] + ".tbi",
        cadd=config["hub_cadd_file"]
    output:
        vcf=OUT + ".vcf.gz",
        index=OUT + ".vcf.gz.tbi",
        cadd=OUT + ".cadd.tsv"
    run:
        shutil.copyfile(input.vcf, output.vcf)
        shutil.copyfile(input.index, output.index)
        shutil.copyfile(input.cadd, output.cadd)

rule filter_prepared_variants:
    input:
        vcf=rules.stage_prepared_fixtures.output.vcf,
        index=rules.stage_prepared_fixtures.output.index,
        cadd=rules.stage_prepared_fixtures.output.cadd
    output:
        vcf=OUT + ".filtered.vcf.gz",
        cadd=OUT + ".filtered.cadd.tsv"
    params:
        script=str(ROOT / "scripts/filtration_vcf_per_fam.sh"),
        config=workflow.configfiles[0]
    threads: config["threads_per_fam"]
    shell:
        "bash {params.script:q} main_filtration -v {input.vcf:q} -c {params.config:q}"

rule prioritize_prepared_variants:
    input:
        vcf=rules.filter_prepared_variants.output.vcf,
        cadd=rules.filter_prepared_variants.output.cadd
    output:
        table=OUT + ".filtered.tsv",
        matrix=OUT + ".filtered.acmg.tsv",
        prepared=OUT + ".prepared.tsv"
    params:
        script=str(ROOT / "scripts/prioritization_vcf_per_fam.sh"),
        config=workflow.configfiles[0]
    threads: config["threads_per_fam"]
    shell:
        """
        bash {params.script:q} prepare_combined_tab {input.vcf:q} {input.cadd:q} {threads} {params.config:q}
        cp -- {output.table:q} {output.prepared:q}
        bash {params.script:q} main_prioritization -i {input.vcf:q} -t {input.cadd:q} -c {params.config:q}
        """
