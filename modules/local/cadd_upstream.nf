// Scientific ordering is upstream CADD v1.7.2; these processes only orchestrate.
def quoteCadd(value) { "'" + value.toString().replace("'", "'\\''") + "'" }

process CADD_UPSTREAM_ANNOTATE {
    tag "${meta.assembly}:annotation"
    label 'cadd_annotation'
    publishDir "${params.outdir}/annotations", mode: 'copy', overwrite: false, pattern: 'annotation_bundle'
    input:
    tuple val(meta), path(vcf), path(profile, stageAs: 'cadd-profile.json')
    output:
    tuple val(meta), path('annotation_bundle'), path(profile), emit: annotations
    script:
    """
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus} NUMEXPR_NUM_THREADS=${task.cpus}
    python3 ${quoteCadd(meta.source_root + '/bin/cadd_task.py')} annotate --profile ${quoteCadd(profile)} --assembly ${quoteCadd(meta.assembly)} --input ${quoteCadd(vcf)} --output annotation_bundle
    """
    stub:
    """
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus} NUMEXPR_NUM_THREADS=${task.cpus}
    python3 ${quoteCadd(meta.source_root + '/bin/cadd_task.py')} annotate --profile ${quoteCadd(profile)} --assembly ${quoteCadd(meta.assembly)} --input ${quoteCadd(vcf)} --output annotation_bundle --stub
    """
}

process CADD_UPSTREAM_CALCULATE {
    tag "${meta.assembly}:calculation"
    label 'cadd_calculation'
    publishDir "${params.outdir}/calculation", mode: 'copy', overwrite: false, pattern: 'calculation_bundle'
    input:
    tuple val(meta), path(annotations, stageAs: 'annotation_bundle'), path(profile, stageAs: 'cadd-profile.json')
    output:
    tuple val(meta), path('calculation_bundle'), emit: scores
    script:
    """
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus} NUMEXPR_NUM_THREADS=${task.cpus}
    python3 ${quoteCadd(meta.source_root + '/bin/cadd_task.py')} calculate --profile ${quoteCadd(profile)} --assembly ${quoteCadd(meta.assembly)} --input ${quoteCadd(annotations)} --output calculation_bundle
    """
    stub:
    """
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus} NUMEXPR_NUM_THREADS=${task.cpus}
    python3 ${quoteCadd(meta.source_root + '/bin/cadd_task.py')} calculate --profile ${quoteCadd(profile)} --assembly ${quoteCadd(meta.assembly)} --input ${quoteCadd(annotations)} --output calculation_bundle --stub
    """
}
