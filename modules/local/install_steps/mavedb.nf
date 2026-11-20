include { installAdopt; installRecord; installSetting; installStart; installStub; installWord } from './common'

// B6: Ensembl's MaveDB table for VEP's MaveDB plugin and its tabix index
// (mavedb_file), GRCh38 only: for GRCh37 the step records no table, and the
// annotation runs without the plugin. priva checks the index against the
// plugin's columns (and rebuilds it with tabix if it cannot read the first
// record).
process MAVEDB {
    tag "${id} ${lineage}"
    label 'install_download'
    storeDir "${params.resource_root}/${id}/${lineage}"
    input:
    tuple val(id), val(lineage), path(plan, stageAs: 'plan.json'), path(receipts, stageAs: 'inputs/*/receipt.json'), val(planned), val(names), path(staged, stageAs: 'staged/*')
    output:
    path('receipt.json'), emit: receipt
    path('data'), emit: data
    script:
    def assembly = installSetting(planned, 'assembly')
    def table = "data/MaveDB/MaveDB_variants.${assembly}.tsv.gz"
    if (lineage.startsWith('adopted-'))
        installAdopt(id, lineage)
    else if (assembly in ['hg19', 'GRCh37'])
        installStart(task) + """
        echo 'MaveDB has no ${assembly} table; the annotation runs without the MaveDB plugin'
        ${installRecord(id, lineage, receipts, [:])}
        """.stripIndent()
    else if (assembly in ['hg38', 'GRCh38'])
        installStart(task) + """
        mkdir -p data/MaveDB
        curl -fL --retry 5 --retry-delay 10 -o ${table} ${installWord(installSetting(planned, 'mavedb_url'))}
        curl -fL --retry 5 --retry-delay 10 -o ${table}.tbi ${installWord(installSetting(planned, 'mavedb_url') + '.tbi')}
        "\$PRIVA_BIN" indexed-plugin-cache ensure-index --input ${table} --kind mavedb
        ${installRecord(id, lineage, receipts, [mavedb_file: table])}
        """.stripIndent()
    else
        "echo 'MaveDB: unsupported assembly ${assembly}' >&2; exit 1"
    stub:
    installStub(id, lineage, receipts, task)
}
