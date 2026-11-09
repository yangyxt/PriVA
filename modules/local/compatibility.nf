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

// family is a PED family identifier, or '' for the no-pedigree run.
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

process PRIORITIZE_FAMILY {
    tag "${family ?: 'no pedigree'}"
    label 'family'
    publishDir "${params.outdir}/prioritization", mode: 'copy', overwrite: false, pattern: 'result/*.tsv'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: false, pattern: '*.provenance.json'
    input:
    tuple val(meta), val(family), path(vcf, stageAs: 'input.vcf.gz'), path(vcf_index, stageAs: 'input.vcf.gz.tbi'), path(cadd, stageAs: 'input.cadd.tsv'), path(bundle, stageAs: 'bundle/*')
    output:
    tuple val(meta), val(family), path('result/*.filtered.tsv'), path('result/*.filtered.acmg.tsv'), emit: tables
    path('prioritization*.provenance.json'), emit: provenance
    script:
    def familyArg = family ? "--family '${family}'" : ''
    """
    '${meta.priva_binary}' task run --phase prioritization --bundle bundle ${familyArg} --vcf '${vcf}' --vcf-index '${vcf_index}' --cadd '${cadd}' --cpus ${task.cpus}
    """
    stub:
    def familyArg = family ? "--family '${family}'" : ''
    """
    '${meta.priva_binary}' task run --phase prioritization --bundle bundle ${familyArg} --vcf '${vcf}' --vcf-index '${vcf_index}' --cadd '${cadd}' --cpus ${task.cpus} --stub
    """
}
