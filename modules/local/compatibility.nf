// One run covers the cohort VCF: validation and annotation run once, then
// filtration and prioritization run once per PED family, or once without a
// family when no PED is given (the shell scripts' no-pedigree mode).

process VALIDATE_INPUTS {
    tag "${meta.base}:${meta.mode}"
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: false, pattern: 'bundle/*'
    input:
    tuple val(meta), path(vcf, stageAs: 'input.vcf.gz'), path(vcf_index, stageAs: 'input.vcf.gz.tbi'), path(ped, stageAs: 'family.ped'), path(config, stageAs: 'legacy.yaml'), path(resources, stageAs: 'resources.json'), path(cadd, stageAs: 'provided.cadd.tsv')
    output:
    tuple val(meta), path(vcf), path(vcf_index), path(cadd), path('bundle/*'), emit: validated
    script:
    def caddArg = cadd ? '--cadd provided.cadd.tsv' : ''
    def pedArg = ped ? '--ped family.ped' : ''
    def familyArg = meta.family ? "--family '${meta.family}'" : ''
    """
    '${meta.priva_binary}' task validate --vcf '${vcf}' --vcf-index '${vcf_index}' ${pedArg} --config '${config}' --resources '${resources}' ${familyArg} --base '${meta.base}' --assembly '${meta.assembly}' --mode '${meta.mode}' --legacy-root '${meta.source_root}' --source-digest '${meta.source_digest}' --software-id '${meta.software_id.replace("'", "'\\''")}' --vcf-hub-mode '${meta.vcf_hub_mode}' --cadd-hub-mode '${meta.cadd_hub_mode}' --task-boundary '${meta.task_boundary}' ${caddArg}
    """
    stub:
    def caddArg = cadd ? '--cadd provided.cadd.tsv' : ''
    def pedArg = ped ? '--ped family.ped' : ''
    def familyArg = meta.family ? "--family '${meta.family}'" : ''
    """
    '${meta.priva_binary}' task validate --vcf '${vcf}' --vcf-index '${vcf_index}' ${pedArg} --config '${config}' --resources '${resources}' ${familyArg} --base '${meta.base}' --assembly '${meta.assembly}' --mode '${meta.mode}' --legacy-root '${meta.source_root}' --source-digest '${meta.source_digest}' --software-id stub --vcf-hub-mode '${meta.vcf_hub_mode}' --cadd-hub-mode '${meta.cadd_hub_mode}' --task-boundary '${meta.task_boundary}' ${caddArg} --stub
    """
}

process ANNOTATE_COHORT {
    tag "${meta.base}:${meta.mode}"
    label 'annotation'
    // Raw mode publishes the new annotation pair; annotated mode only imports it.
    publishDir "${params.outdir}/annotation", mode: 'copy', overwrite: false, pattern: 'result/*.anno.*', saveAs: { new File(it).name }, enabled: params.mode == 'raw'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: false, pattern: '*.provenance.json'
    input:
    tuple val(meta), path(vcf, stageAs: 'input.vcf.gz'), path(vcf_index, stageAs: 'input.vcf.gz.tbi'), path(cadd, stageAs: 'provided.cadd.tsv'), path(bundle, stageAs: 'bundle/*')
    output:
    tuple val(meta), path('result/*.anno.vcf.gz'), path('result/*.anno.vcf.gz.tbi'), path('result/*.anno.cadd.tsv'), path('bundle/*', includeInputs: true), emit: annotated
    path('annotation.provenance.json'), emit: provenance
    script:
    def caddArg = cadd ? '--cadd provided.cadd.tsv' : ''
    """
    '${meta.priva_binary}' task run --phase annotation --bundle bundle --vcf '${vcf}' --vcf-index '${vcf_index}' --cpus ${task.cpus} ${caddArg}
    """
    stub:
    def caddArg = cadd ? '--cadd provided.cadd.tsv' : ''
    """
    '${meta.priva_binary}' task run --phase annotation --bundle bundle --vcf '${vcf}' --vcf-index '${vcf_index}' --cpus ${task.cpus} ${caddArg} --stub
    """
}

// Raw mode with both hubs off: the annotation script's functions as their own
// tasks. Preprocessing comes first; gnomAD, ClinVar, control cohort, VEP and
// CADD then run at the same time on copies of the preprocessed VCF, and the
// merge copies each step's INFO fields onto it in the script's order.

process ANNOTATE_PREPROCESS {
    tag "${meta.base}"
    label 'annotation_step'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: false, pattern: '*.provenance.json'
    input:
    tuple val(meta), path(vcf, stageAs: 'input.vcf.gz'), path(vcf_index, stageAs: 'input.vcf.gz.tbi'), path(bundle, stageAs: 'bundle/*')
    output:
    tuple val(meta), path("result/${meta.base}.anno.vcf.gz"), path("result/${meta.base}.anno.vcf.gz.tbi"), path('bundle/*', includeInputs: true), emit: preprocessed
    path('annotation-preprocess.provenance.json'), emit: provenance
    script:
    """
    '${meta.priva_binary}' task run --phase annotation-preprocess --bundle bundle --vcf '${vcf}' --vcf-index '${vcf_index}' --cpus ${task.cpus}
    """
    stub:
    """
    '${meta.priva_binary}' task run --phase annotation-preprocess --bundle bundle --vcf '${vcf}' --vcf-index '${vcf_index}' --cpus ${task.cpus} --stub
    """
}

// step is gnomad, clinvar, control or vep. VEP gets the annotation CPUs.
process ANNOTATE_STEP {
    tag "${meta.base}:${step}"
    cpus { step == 'vep' ? params.annotation_cpus : params.annotation_step_cpus }
    memory { step == 'vep' ? params.annotation_memory : params.annotation_step_memory }
    time { step == 'vep' ? params.annotation_time : params.annotation_step_time }
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: false, pattern: '*.provenance.json'
    input:
    tuple val(meta), val(step), path(vcf, stageAs: 'preprocessed.vcf.gz'), path(vcf_index, stageAs: 'preprocessed.vcf.gz.tbi'), path(bundle, stageAs: 'bundle/*')
    output:
    tuple val(meta), val(step), path("result/${meta.base}.anno.vcf.gz"), emit: annotated
    path("annotation-${step}.provenance.json"), emit: provenance
    script:
    """
    '${meta.priva_binary}' task run --phase annotation-${step} --bundle bundle --vcf '${vcf}' --vcf-index '${vcf_index}' --cpus ${task.cpus}
    """
    stub:
    """
    '${meta.priva_binary}' task run --phase annotation-${step} --bundle bundle --vcf '${vcf}' --vcf-index '${vcf_index}' --cpus ${task.cpus} --stub
    """
}

process ANNOTATE_CADD {
    tag "${meta.base}"
    label 'annotation'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: false, pattern: '*.provenance.json'
    input:
    tuple val(meta), path(vcf, stageAs: 'preprocessed.vcf.gz'), path(vcf_index, stageAs: 'preprocessed.vcf.gz.tbi'), path(bundle, stageAs: 'bundle/*')
    output:
    tuple val(meta), path("result/${meta.base}.anno.cadd.tsv"), emit: cadd
    path('annotation-cadd.provenance.json'), emit: provenance
    script:
    """
    '${meta.priva_binary}' task run --phase annotation-cadd --bundle bundle --vcf '${vcf}' --vcf-index '${vcf_index}' --cpus ${task.cpus}
    """
    stub:
    """
    '${meta.priva_binary}' task run --phase annotation-cadd --bundle bundle --vcf '${vcf}' --vcf-index '${vcf_index}' --cpus ${task.cpus} --stub
    """
}

process ANNOTATE_MERGE {
    tag "${meta.base}"
    label 'annotation_step'
    publishDir "${params.outdir}/annotation", mode: 'copy', overwrite: false, pattern: 'result/*.anno.*', saveAs: { new File(it).name }
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: false, pattern: '*.provenance.json'
    input:
    tuple val(meta), path(vcf, stageAs: 'preprocessed.vcf.gz'), path(vcf_index, stageAs: 'preprocessed.vcf.gz.tbi'), path(gnomad, stageAs: 'gnomad.vcf.gz'), path(clinvar, stageAs: 'clinvar.vcf.gz'), path(control, stageAs: 'control.vcf.gz'), path(vep, stageAs: 'vep.vcf.gz'), path(cadd, stageAs: 'cohort.cadd.tsv'), path(bundle, stageAs: 'bundle/*')
    output:
    tuple val(meta), path('result/*.anno.vcf.gz'), path('result/*.anno.vcf.gz.tbi'), path('result/*.anno.cadd.tsv'), path('bundle/*', includeInputs: true), emit: annotated
    path('annotation-merge.provenance.json'), emit: provenance
    script:
    """
    '${meta.priva_binary}' task run --phase annotation-merge --bundle bundle --vcf '${vcf}' --vcf-index '${vcf_index}' --gnomad-vcf gnomad.vcf.gz --clinvar-vcf clinvar.vcf.gz --control-vcf control.vcf.gz --vep-vcf vep.vcf.gz --cadd cohort.cadd.tsv --cpus ${task.cpus}
    """
    stub:
    """
    '${meta.priva_binary}' task run --phase annotation-merge --bundle bundle --vcf '${vcf}' --vcf-index '${vcf_index}' --gnomad-vcf gnomad.vcf.gz --clinvar-vcf clinvar.vcf.gz --control-vcf control.vcf.gz --vep-vcf vep.vcf.gz --cadd cohort.cadd.tsv --cpus ${task.cpus} --stub
    """
}

// family is a PED family identifier, or '' for the no-pedigree run.
// Annotation hub (raw mode, --annotation_hub). The hub keeps each step's
// results per variant in HUB/SOURCE/LINEAGE/CHROM/part-SESSION.parquet; files
// are only added. A run never reads the files its own session wrote, so a
// resumed run sees what its first launch saw.
process HUB_LOOKUP {
    tag "${source}"
    label 'annotation_step'
    input:
    tuple val(meta), val(source), val(lineage), path(vcf, stageAs: 'preprocessed.vcf.gz'), path(vcf_index, stageAs: 'preprocessed.vcf.gz.tbi')
    output:
    tuple val(meta), val(source), val(lineage), path('known.parquet'), path("new/${meta.base}.vcf.gz"), path("new/${meta.base}.vcf.gz.tbi"), env('NEW_VARIANTS'), emit: split
    script:
    """
    mkdir -p new
    '${meta.priva_binary}' annotation-hub lookup --hub '${params.annotation_hub}' --source ${source} --lineage '${lineage}' --exclude-session '${workflow.sessionId}' --vcf preprocessed.vcf.gz --known known.parquet --new 'new/${meta.base}.vcf.gz'
    NEW_VARIANTS=\$(gzip -dc 'new/${meta.base}.vcf.gz' | grep -vc '^#' || true)
    """
}

// Hub refresh (--action hub_refresh): a VCF of the variants that earlier hub
// folders hold for every step whose current folder is still empty.
process HUB_REFRESH_PLAN {
    tag 'plan'
    label 'annotation_step'
    input:
    tuple val(meta), path(config, stageAs: 'legacy.yaml'), path(resources, stageAs: 'resources.json')
    output:
    tuple val(meta), path("refresh/${meta.base}.vcf.gz"), path("refresh/${meta.base}.vcf.gz.tbi"), path('refresh/plan.json'), emit: plan
    script:
    """
    mkdir -p refresh
    '${meta.priva_binary}' annotation-hub refresh-plan --hub '${params.annotation_hub}' --resources resources.json --legacy-root '${meta.source_root}' --vcf 'refresh/${meta.base}.vcf.gz' --plan refresh/plan.json
    """
}

process HUB_FILL {
    tag "${source}"
    label 'annotation_step'
    input:
    tuple val(meta), val(source), path(vcf, stageAs: 'preprocessed.vcf.gz'), path(vcf_index, stageAs: 'preprocessed.vcf.gz.tbi'), path(known, stageAs: 'known.parquet'), path(output, stageAs: 'step/*')
    output:
    tuple val(meta), val(source), path(source == 'cadd' ? "filled/${meta.base}.anno.cadd.tsv" : "filled/${meta.base}.anno.vcf.gz"), emit: filled
    script:
    def name = source == 'cadd' ? "${meta.base}.anno.cadd.tsv" : "${meta.base}.anno.vcf.gz"
    """
    mkdir -p filled
    '${meta.priva_binary}' annotation-hub fill --source ${source} --vcf preprocessed.vcf.gz --known known.parquet --annotated '${output}' --output 'filled/${name}'
    """
}

process HUB_ADD {
    tag "${source}"
    label 'annotation_step'
    input:
    tuple val(meta), val(source), val(lineage), path(fresh, stageAs: 'new.vcf.gz'), path(output, stageAs: 'step/*')
    output:
    tuple val(source), path('added.log'), emit: added
    script:
    """
    '${meta.priva_binary}' annotation-hub add --hub '${params.annotation_hub}' --source ${source} --lineage '${lineage}' --session '${workflow.sessionId}' --base new.vcf.gz --annotated '${output}' 2> added.log
    """
}

process FILTER_FAMILY {
    tag "${family ?: 'no pedigree'}"
    label 'family'
    publishDir "${params.outdir}/variants", mode: 'copy', overwrite: false, pattern: 'result/*.filtered.*'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: false, pattern: '*.provenance.json'
    input:
    tuple val(meta), val(family), path(vcf, stageAs: 'input.vcf.gz'), path(vcf_index, stageAs: 'input.vcf.gz.tbi'), path(cadd, stageAs: 'input.cadd.tsv'), path(bundle, stageAs: 'bundle/*')
    output:
    tuple val(meta), val(family), path('result/*.filtered.vcf.gz'), path('result/*.filtered.vcf.gz.tbi'), path('result/*.filtered.cadd.tsv'), path('bundle/*', includeInputs: true), emit: filtered
    path('filtration*.provenance.json'), emit: provenance
    script:
    def familyArg = family ? "--family '${family}'" : ''
    """
    '${meta.priva_binary}' task run --phase filtration --bundle bundle ${familyArg} --vcf '${vcf}' --vcf-index '${vcf_index}' --cadd '${cadd}' --cpus ${task.cpus}
    """
    stub:
    def familyArg = family ? "--family '${family}'" : ''
    """
    '${meta.priva_binary}' task run --phase filtration --bundle bundle ${familyArg} --vcf '${vcf}' --vcf-index '${vcf_index}' --cadd '${cadd}' --cpus ${task.cpus} --stub
    """
}

// Prioritization, one task per step and family: combine builds the family
// table; splicing and 5' UTR add their columns to copies of it side by side;
// join puts the copies together; ACMG assigns the criteria.
process PRIORITIZE_COMBINE {
    tag "${family ?: 'no pedigree'}"
    label 'family'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: false, pattern: '*.provenance.json'
    input:
    tuple val(meta), val(family), path(vcf, stageAs: 'input.vcf.gz'), path(vcf_index, stageAs: 'input.vcf.gz.tbi'), path(cadd, stageAs: 'input.cadd.tsv'), path(bundle, stageAs: 'bundle/*')
    output:
    tuple val(meta), val(family), path('result/*.filtered.tsv'), path('bundle/*', includeInputs: true), emit: table
    path('prioritization-combine*.provenance.json'), emit: provenance
    script:
    def familyArg = family ? "--family '${family}'" : ''
    """
    '${meta.priva_binary}' task run --phase prioritization-combine --bundle bundle ${familyArg} --vcf '${vcf}' --vcf-index '${vcf_index}' --cadd '${cadd}' --cpus ${task.cpus}
    """
    stub:
    def familyArg = family ? "--family '${family}'" : ''
    """
    '${meta.priva_binary}' task run --phase prioritization-combine --bundle bundle ${familyArg} --vcf '${vcf}' --vcf-index '${vcf_index}' --cadd '${cadd}' --cpus ${task.cpus} --stub
    """
}

process PRIORITIZE_SPLICING {
    tag "${family ?: 'no pedigree'}"
    label 'family'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: false, pattern: '*.provenance.json'
    input:
    tuple val(meta), val(family), path(table, stageAs: 'input.tsv'), path(bundle, stageAs: 'bundle/*')
    output:
    tuple val(meta), val(family), path('result/*.filtered.tsv'), emit: table
    path('prioritization-splicing*.provenance.json'), emit: provenance
    script:
    def familyArg = family ? "--family '${family}'" : ''
    """
    '${meta.priva_binary}' task run --phase prioritization-splicing --bundle bundle ${familyArg} --table '${table}' --cpus ${task.cpus}
    """
    stub:
    def familyArg = family ? "--family '${family}'" : ''
    """
    '${meta.priva_binary}' task run --phase prioritization-splicing --bundle bundle ${familyArg} --table '${table}' --cpus ${task.cpus} --stub
    """
}

process PRIORITIZE_UTR {
    tag "${family ?: 'no pedigree'}"
    label 'family'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: false, pattern: '*.provenance.json'
    input:
    tuple val(meta), val(family), path(table, stageAs: 'input.tsv'), path(bundle, stageAs: 'bundle/*')
    output:
    tuple val(meta), val(family), path('result/*.filtered.tsv'), emit: table
    path('prioritization-utr*.provenance.json'), emit: provenance
    script:
    def familyArg = family ? "--family '${family}'" : ''
    """
    '${meta.priva_binary}' task run --phase prioritization-utr --bundle bundle ${familyArg} --table '${table}' --cpus ${task.cpus}
    """
    stub:
    def familyArg = family ? "--family '${family}'" : ''
    """
    '${meta.priva_binary}' task run --phase prioritization-utr --bundle bundle ${familyArg} --table '${table}' --cpus ${task.cpus} --stub
    """
}

process PRIORITIZE_JOIN {
    tag "${family ?: 'no pedigree'}"
    label 'family'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: false, pattern: '*.provenance.json'
    input:
    tuple val(meta), val(family), path(table, stageAs: 'input.tsv'), path(splicing, stageAs: 'splicing.tsv'), path(utr, stageAs: 'utr.tsv'), path(bundle, stageAs: 'bundle/*')
    output:
    tuple val(meta), val(family), path('result/*.filtered.tsv'), path('bundle/*', includeInputs: true), emit: table
    path('prioritization-join*.provenance.json'), emit: provenance
    script:
    def familyArg = family ? "--family '${family}'" : ''
    """
    '${meta.priva_binary}' task run --phase prioritization-join --bundle bundle ${familyArg} --table '${table}' --splicing-table '${splicing}' --utr-table '${utr}' --cpus ${task.cpus}
    """
    stub:
    def familyArg = family ? "--family '${family}'" : ''
    """
    '${meta.priva_binary}' task run --phase prioritization-join --bundle bundle ${familyArg} --table '${table}' --splicing-table '${splicing}' --utr-table '${utr}' --cpus ${task.cpus} --stub
    """
}

process PRIORITIZE_ACMG {
    tag "${family ?: 'no pedigree'}"
    label 'family'
    publishDir "${params.outdir}/prioritization", mode: 'copy', overwrite: false, pattern: 'result/*.tsv'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: false, pattern: '*.provenance.json'
    input:
    tuple val(meta), val(family), path(table, stageAs: 'input.tsv'), path(bundle, stageAs: 'bundle/*')
    output:
    tuple val(meta), val(family), path('result/*.filtered.tsv'), path('result/*.filtered.acmg.tsv'), emit: tables
    path('prioritization-acmg*.provenance.json'), emit: provenance
    script:
    def familyArg = family ? "--family '${family}'" : ''
    """
    '${meta.priva_binary}' task run --phase prioritization-acmg --bundle bundle ${familyArg} --table '${table}' --cpus ${task.cpus}
    """
    stub:
    def familyArg = family ? "--family '${family}'" : ''
    """
    '${meta.priva_binary}' task run --phase prioritization-acmg --bundle bundle ${familyArg} --table '${table}' --cpus ${task.cpus} --stub
    """
}
