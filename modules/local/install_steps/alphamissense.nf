include { installAdopt; installRecord; installSetting; installStart; installStub } from './common'

// B3: the publisher's original AlphaMissense table for VEP's AlphaMissense
// plugin (alphamissense_prescore). priva alphamissense-cache downloads it with
// curl, checks its size and MD5 against the pinned publisher object, indexes
// it with tabix, reads a record back through the index and writes the ready
// receipt beside it; it prints the table's path. The ready receipt names that
// path, so record rewrites it to where storeDir moves the table.
process ALPHAMISSENSE_TSV {
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
    if (lineage.startsWith('adopted-'))
        installAdopt(id, lineage)
    else
        installStart(task) + """
        here=\$(pwd -P)
        selected=\$("\$PRIVA_BIN" alphamissense-cache install --cache-dir data --assembly ${assembly})
        table=\${selected#"\$here"/}
        [[ \$table == data/* && -s \$table && -s \$table.tbi ]] || { echo "AlphaMissense: unexpected table \$selected" >&2; exit 1; }
        ${installRecord(id, lineage, receipts, [:])} --key "alphamissense_prescore=\$table" --relocate "\${table%.tsv.gz}.ready.json"
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
