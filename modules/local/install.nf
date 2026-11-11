// Install tasks. Each install step of conf/install_steps.json runs as its own
// task through scripts/install_task.sh. Its results live in
// resource_root/STEP/LINEAGE/; a folder that already holds receipt.json is
// reused without running (storeDir), across runs and after the work folder
// is deleted.

process INSTALL_PLAN {
    tag 'plan'
    label 'install_light'
    input:
    tuple path(config, stageAs: 'install_config.yaml'), path(registry, stageAs: 'install_steps.json'), val(adopt_config), val(adopt), val(run_steps)
    output:
    path('plan/plan.json'), emit: plan
    script:
    def adoptArgs = adopt_config ? "--adopt-config '${adopt_config}' --adopt '${adopt}' --run '${run_steps}'" : ''
    """
    export PRIVA_BINARY='${params.priva_binary}'
    bash '${projectDir}/scripts/install_task.sh' plan --registry '${registry}' --config '${config}' --out plan ${adoptArgs}
    """
}

process INSTALL_STEP {
    tag "${id} ${lineage}"
    label 'install'
    storeDir "${params.resource_root}/${id}/${lineage}"
    input:
    tuple val(id), val(lineage), path(plan, stageAs: 'plan.json'), val(inputs)
    output:
    path('receipt.json'), emit: receipt
    script:
    def inputArgs = inputs.collect { "--input '${it}'" }.join(' ')
    def folder = "${params.resource_root}/${id}/${lineage}"
    if (lineage.startsWith('adopted-'))
        """
        bash '${projectDir}/scripts/install_task.sh' adopt --plan plan.json --step ${id} --lineage ${lineage} --receipt receipt.json
        """
    else
        """
        export PRIVA_BINARY='${params.priva_binary}' INSTALL_THREADS=${task.cpus}
        bash '${projectDir}/scripts/install_task.sh' run --registry '${projectDir}/conf/install_steps.json' --plan plan.json --step ${id} --lineage ${lineage} --folder '${folder}' --receipt receipt.json ${inputArgs}
        """
    stub:
    def inputArgs = inputs.collect { "--input '${it}'" }.join(' ')
    def folder = "${params.resource_root}/${id}/${lineage}"
    if (lineage.startsWith('adopted-'))
        """
        bash '${projectDir}/scripts/install_task.sh' adopt --plan plan.json --step ${id} --lineage ${lineage} --receipt receipt.json
        """
    else
        """
        export PRIVA_BINARY='${params.priva_binary}' INSTALL_THREADS=${task.cpus}
        bash '${projectDir}/scripts/install_task.sh' run --registry '${projectDir}/conf/install_steps.json' --plan plan.json --step ${id} --lineage ${lineage} --folder '${folder}' --receipt receipt.json ${inputArgs} --stub
        """
}

process INSTALL_PUBLISH {
    tag "publish ${lineage}"
    label 'install_light'
    storeDir "${params.resource_root}/PUBLISH/${lineage}"
    publishDir "${params.install_outdir}", mode: 'copy', overwrite: true, enabled: params.install_outdir as boolean
    input:
    tuple val(lineage), path(plan, stageAs: 'plan.json'), val(inputs)
    output:
    tuple path('config.yaml'), path('resources.lock.json'), emit: published
    script:
    def inputArgs = inputs.collect { "--input '${it}'" }.join(' ')
    """
    export PRIVA_BINARY='${params.priva_binary}'
    bash '${projectDir}/scripts/install_task.sh' publish --registry '${projectDir}/conf/install_steps.json' --plan plan.json --lineage ${lineage} --folder '${params.resource_root}/PUBLISH/${lineage}' ${inputArgs}
    """
    stub:
    def inputArgs = inputs.collect { "--input '${it}'" }.join(' ')
    """
    export PRIVA_BINARY='${params.priva_binary}'
    bash '${projectDir}/scripts/install_task.sh' publish --registry '${projectDir}/conf/install_steps.json' --plan plan.json --lineage ${lineage} --folder '${params.resource_root}/PUBLISH/${lineage}' ${inputArgs} --stub
    """
}
