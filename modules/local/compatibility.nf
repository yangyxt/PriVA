process VALIDATE_FAMILY {
    tag "${meta.id}:${meta.mode}"
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: false, pattern: 'bundle'
    input:
    tuple val(meta), path(vcf, stageAs: 'input.vcf.gz'), path(vcf_index, stageAs: 'input.vcf.gz.tbi'), path(ped, stageAs: 'family.ped'), path(config, stageAs: 'legacy.yaml'), path(resources, stageAs: 'resources.json'), path(cadd, stageAs: 'provided.cadd.tsv')
    output:
    tuple val(meta), path(vcf), path(vcf_index), path(cadd), path('bundle'), emit: validated
    script:
    def caddArg = cadd ? '--cadd provided.cadd.tsv' : ''
    """
    python3 '${meta.source_root}/bin/priva_task.py' validate --vcf '${vcf}' --vcf-index '${vcf_index}' --ped '${ped}' --config '${config}' --resources '${resources}' --family '${meta.id}' --base '${meta.base}' --assembly '${meta.assembly}' --mode '${meta.mode}' --legacy-root '${meta.source_root}' --source-digest '${meta.source_digest}' --software-id '${meta.software_id.replace("'", "'\\''")}' --vcf-hub-mode '${meta.vcf_hub_mode}' --cadd-hub-mode '${meta.cadd_hub_mode}' --task-boundary '${meta.task_boundary}' ${caddArg}
    """
    stub:
    def caddArg = cadd ? '--cadd provided.cadd.tsv' : ''
    """
    python3 '${meta.source_root}/bin/priva_task.py' validate --vcf '${vcf}' --vcf-index '${vcf_index}' --ped '${ped}' --config '${config}' --resources '${resources}' --family '${meta.id}' --base '${meta.base}' --assembly '${meta.assembly}' --mode '${meta.mode}' --legacy-root '${meta.source_root}' --source-digest '${meta.source_digest}' --software-id stub --vcf-hub-mode '${meta.vcf_hub_mode}' --cadd-hub-mode '${meta.cadd_hub_mode}' --task-boundary '${meta.task_boundary}' ${caddArg} --stub
    """
}

process ANNOTATE_FAMILY {
    tag "${meta.id}:${meta.mode}"
    label 'annotation'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: false, pattern: '*.provenance.json'
    input:
    tuple val(meta), path(vcf, stageAs: 'input.vcf.gz'), path(vcf_index, stageAs: 'input.vcf.gz.tbi'), path(cadd, stageAs: 'provided.cadd.tsv'), path(bundle)
    output:
    tuple val(meta), path('result/*.anno.vcf.gz'), path('result/*.anno.vcf.gz.tbi'), path('result/*.anno.cadd.tsv'), path(bundle), emit: annotated
    path('annotation.provenance.json'), emit: provenance
    script:
    def caddArg = cadd ? '--cadd provided.cadd.tsv' : ''
    """
    python3 '${meta.source_root}/bin/priva_task.py' run --phase annotation --bundle '${bundle}' --vcf '${vcf}' --vcf-index '${vcf_index}' --cpus ${task.cpus} ${caddArg}
    """
    stub:
    def caddArg = cadd ? '--cadd provided.cadd.tsv' : ''
    """
    python3 '${meta.source_root}/bin/priva_task.py' run --phase annotation --bundle '${bundle}' --vcf '${vcf}' --vcf-index '${vcf_index}' --cpus ${task.cpus} ${caddArg} --stub
    """
}

process FILTER_FAMILY {
    tag "${meta.id}"
    label 'family'
    publishDir "${params.outdir}/variants", mode: 'copy', overwrite: false, pattern: 'result/*.filtered.*'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: false, pattern: '*.provenance.json'
    input:
    tuple val(meta), path(vcf, stageAs: 'input.vcf.gz'), path(vcf_index, stageAs: 'input.vcf.gz.tbi'), path(cadd, stageAs: 'input.cadd.tsv'), path(bundle)
    output:
    tuple val(meta), path('result/*.filtered.vcf.gz'), path('result/*.filtered.vcf.gz.tbi'), path('result/*.filtered.cadd.tsv'), path(bundle), emit: filtered
    path('filtration.provenance.json'), emit: provenance
    script:
    """
    python3 '${meta.source_root}/bin/priva_task.py' run --phase filtration --bundle '${bundle}' --vcf '${vcf}' --vcf-index '${vcf_index}' --cadd '${cadd}' --cpus ${task.cpus}
    """
    stub:
    """
    python3 '${meta.source_root}/bin/priva_task.py' run --phase filtration --bundle '${bundle}' --vcf '${vcf}' --vcf-index '${vcf_index}' --cadd '${cadd}' --cpus ${task.cpus} --stub
    """
}

process PRIORITIZE_FAMILY {
    tag "${meta.id}"
    label 'family'
    publishDir "${params.outdir}/prioritization", mode: 'copy', overwrite: false, pattern: 'result/*.tsv'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: false, pattern: '*.provenance.json'
    input:
    tuple val(meta), path(vcf, stageAs: 'input.vcf.gz'), path(vcf_index, stageAs: 'input.vcf.gz.tbi'), path(cadd, stageAs: 'input.cadd.tsv'), path(bundle)
    output:
    tuple val(meta), path('result/*.filtered.tsv'), path('result/*.filtered.acmg.tsv'), emit: tables
    path('prioritization.provenance.json'), emit: provenance
    script:
    """
    python3 '${meta.source_root}/bin/priva_task.py' run --phase prioritization --bundle '${bundle}' --vcf '${vcf}' --vcf-index '${vcf_index}' --cadd '${cadd}' --cpus ${task.cpus}
    """
    stub:
    """
    python3 '${meta.source_root}/bin/priva_task.py' run --phase prioritization --bundle '${bundle}' --vcf '${vcf}' --vcf-index '${vcf_index}' --cadd '${cadd}' --cpus ${task.cpus} --stub
    """
}
