nextflow.enable.dsl=2

include { CADD_UPSTREAM } from './modules/local/cadd_upstream'

workflow {
    if (params.help) {
        log.info '''PriVA standalone CADD: upstream CADD.sh v1.7.2 with all of CADD's own data
nextflow run cadd.nf -c conf/cadd.config --input INPUT.vcf.gz --assembly GRCh37|GRCh38
  --cadd_base_dir CADD_SCRIPTS_DIR --outdir NEW_OUTPUT_DIRECTORY -work-dir DURABLE_WORK
cadd_base_dir is the CADD-scripts folder prepared by main_install: CADD.sh, its
data and Snakemake environments, PriVA's faster copies of its prescore lookup and
MMSplice scripts (scripts/cadd_scripts) and the prescored position indexes.
CADD.sh's Snakemake needs a conda program (24.7.1 or later) on PATH.
The output INPUT.cadd.tsv is the table main.nf --mode annotated accepts as --cadd.
-stub-run writes a header-only table and runs no CADD.'''
    } else {
        ['input', 'assembly', 'cadd_base_dir', 'outdir'].each { key ->
            if (!params[key]) error "Required parameter --${key} is missing; include -c conf/cadd.config"
        }
        if (!(params.assembly in ['hg19', 'hg38', 'GRCh37', 'GRCh38'])) error 'Unsupported --assembly'
        if (workflow.profile.tokenize(',').contains('test') && !workflow.stubRun) error 'The test profile requires -stub-run'
        def base = file(params.cadd_base_dir)
        if (!workflow.stubRun && !file("${base}/CADD.sh").exists()) error "CADD.sh was not found in --cadd_base_dir ${base}"
        def input = file(params.input, checkIfExists: true)
        def meta = [id: input.name.replaceAll(/\.(vcf|bcf)(\.b?gz)?$/, ''),
                    genome: params.assembly in ['hg19', 'GRCh37'] ? 'GRCh37' : 'GRCh38',
                    cadd_base_dir: base.toString(), source_root: projectDir.toString()]
        CADD_UPSTREAM(Channel.of(tuple(meta, input)))
    }
}
