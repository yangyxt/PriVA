// The install runner's tasks (priva install-task, rust/src/install_task/):
// releases, the plan, the preflight checks and the publish. Each install step
// of conf/install_steps.json is its own process in modules/local/install_steps/;
// its results live in resource_root/STEP/LINEAGE/, and a folder that already
// holds receipt.json is reused without running (storeDir), across runs and
// after the work folder is deleted.

// The current release of a changing source (registry "release"); pinned
// sources are checked only when --refresh names them.
process INSTALL_RELEASE {
    tag "${id}"
    label 'install_light'
    input:
    tuple val(id), path(config, stageAs: 'install_config.yaml')
    output:
    path("release_${id}.json"), emit: release
    script:
    """
    '${params.priva_binary}' install-task release --source-root '${projectDir}' --registry '${projectDir}/conf/install_steps.json' --config install_config.yaml --step ${id} --out 'release_${id}.json'
    """
    stub:
    """
    '${params.priva_binary}' install-task release --source-root '${projectDir}' --registry '${projectDir}/conf/install_steps.json' --config install_config.yaml --step ${id} --out 'release_${id}.json' --stub
    """
}

process INSTALL_PLAN {
    tag 'plan'
    label 'install_light'
    input:
    tuple path(config, stageAs: 'install_config.yaml'), path(registry, stageAs: 'install_steps.json'), val(adopt_config), val(adopt), val(run_steps), val(refresh), path(releases, stageAs: 'releases/*')
    output:
    path('plan/plan.json'), emit: plan
    script:
    def adoptArgs = adopt_config ? "--adopt-config '${adopt_config}' --adopt '${adopt}' --run '${run_steps}'" : ''
    def releaseArgs = (releases instanceof List ? releases : [releases]).findAll { r -> r }.collect { r -> "--releases '${r}'" }.join(' ')
    """
    export PRIVA_BINARY='${params.priva_binary}'
    '${params.priva_binary}' install-task plan --source-root '${projectDir}' --registry '${registry}' --config '${config}' --out plan ${adoptArgs} --pins '${params.resource_root}/pins.json' --refresh '${refresh}' ${releaseArgs}
    """
}

// priva's checks before any step runs (install-task preflight): the files
// you supply for the steps that will run, and the data files shipped in the
// checkout. The plan passes through it, so no step starts before they pass.
process INSTALL_PREFLIGHT {
    tag 'preflight'
    label 'install_light'
    input:
    path(plan, stageAs: 'plan.json')
    output:
    path('checked/plan.json'), emit: plan
    script:
    """
    '${params.priva_binary}' install-task preflight --source-root '${projectDir}' --registry '${projectDir}/conf/install_steps.json' --plan plan.json
    mkdir checked && cp plan.json checked/plan.json
    """
    stub:
    """
    '${params.priva_binary}' install-task preflight --source-root '${projectDir}' --registry '${projectDir}/conf/install_steps.json' --plan plan.json --stub
    mkdir checked && cp plan.json checked/plan.json
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
    '${params.priva_binary}' install-task publish --source-root '${projectDir}' --registry '${projectDir}/conf/install_steps.json' --plan plan.json --lineage ${lineage} --folder '${params.resource_root}/PUBLISH/${lineage}' ${inputArgs} --pins '${params.resource_root}/pins.json'
    """
    stub:
    def inputArgs = inputs.collect { "--input '${it}'" }.join(' ')
    """
    export PRIVA_BINARY='${params.priva_binary}'
    '${params.priva_binary}' install-task publish --source-root '${projectDir}' --registry '${projectDir}/conf/install_steps.json' --plan plan.json --lineage ${lineage} --folder '${params.resource_root}/PUBLISH/${lineage}' ${inputArgs} --pins '${params.resource_root}/pins.json' --stub
    """
}
