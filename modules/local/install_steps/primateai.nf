include { installAdopt; installRecord; installSetting; installStart; installStub; installWord } from './common'

// B7: PrimateAI-3D for VEP's PrimateAI plugin (primateai_prescore). Illumina
// gives the table only after you sign in, so you download it and set
// primateai_prescore; the task gets it staged. priva indexed-plugin-cache
// format-primate sorts the publisher's table with the plugin's own procedure,
// compresses it with bgzip and indexes it with tabix into data/PrimateAI, and
// prints its path. A table you supply already sorted and indexed is used where
// it is; when its index is missing or older than it, priva writes a new index
// beside the staged link, and that index is kept with a link to your table.
process PRIMATEAI {
    tag "${id} ${lineage}"
    label 'install'
    storeDir "${params.resource_root}/${id}/${lineage}"
    input:
    tuple val(id), val(lineage), path(plan, stageAs: 'plan.json'), path(receipts, stageAs: 'inputs/*/receipt.json'), val(planned), val(names), path(staged, stageAs: 'staged/*')
    output:
    path('receipt.json'), emit: receipt
    path('data'), emit: data
    script:
    def assembly = installSetting(planned, 'assembly')
    def source = names.primateai_prescore ? "staged/${names.primateai_prescore}".toString() : ''
    def configured = planned.user_paths?.primateai_prescore ?: ''
    if (lineage.startsWith('adopted-'))
        installAdopt(id, lineage)
    else if (!source)
        "echo ${installWord("PrimateAI: set primateai_prescore to the assembly-matched PrimateAI-3D table you downloaded from https://basespace.illumina.com/s/yYGFdGih1rXL (missing: ${configured})")} >&2; exit 1"
    else
        installStart(task) + """
        here=\$(pwd -P)
        source=${installWord(source)}
        mkdir -p data/PrimateAI
        prepared=\$("\$PRIVA_BIN" indexed-plugin-cache format-primate --input "\$source" --cache-dir data/PrimateAI --assembly ${assembly} --threads 1 --print-path)
        if [[ \$prepared == "\$here"/data/* ]]; then
            table=\${prepared#"\$here"/}
        elif [[ \$prepared == "\$here/\$source" && -L \$source.tbi ]]; then
            table=\$source
        elif [[ \$prepared == "\$here/\$source" ]]; then
            ln -s "\$(readlink "\$source")" "data/PrimateAI/\$(basename "\$source")"
            mv "\$source.tbi" "data/PrimateAI/\$(basename "\$source").tbi"
            table="data/PrimateAI/\$(basename "\$source")"
        else
            echo "PrimateAI: unexpected table \$prepared" >&2; exit 1
        fi
        ${installRecord(id, lineage, receipts, [:])} --key "primateai_prescore=\$table"
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
