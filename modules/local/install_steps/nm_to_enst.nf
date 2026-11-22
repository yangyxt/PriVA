include { installAdopt; installRecord; installSetting; installStart; installStub; installWord } from './common'

// A5: the RefSeq-to-Ensembl transcript table (nm_to_enst_tsv) of the assembly
// and VEP release. priva installer-prerequisites picks the table: one you set
// in nm_to_enst_tsv, or the pinned table the repository ships in data/hgnc/
// for the assembly and VEP release. When the shipped table is missing, priva
// biomart fetches the pinned Ensembl BioMart export into this step's folder
// and checks it against its pinned provenance. installer-prerequisites then
// checks the table it picked. A table you set is used where it is: the step
// keeps a link to it (and to its .provenance.json).
process NM_TO_ENST {
    tag "${id} ${lineage}"
    label 'install_download'
    storeDir "${params.resource_root}/${id}/${lineage}"
    input:
    tuple val(id), val(lineage), path(plan, stageAs: 'plan.json'), path(receipts, stageAs: 'inputs/*/receipt.json'), val(planned), val(names), path(staged, stageAs: 'staged/*')
    output:
    path('receipt.json'), emit: receipt
    path('data'), emit: data
    script:
    def adopted = lineage.startsWith('adopted-')
    def assembly = adopted ? '' : installSetting(planned, 'assembly')
    def release = adopted ? '' : installSetting(planned, 'vep_installed_version')
    def source = names.nm_to_enst_tsv ? "staged/${names.nm_to_enst_tsv}".toString() : ''
    def configured = planned.user_paths?.nm_to_enst_tsv ?: ''
    def root = installWord(projectDir)
    if (adopted)
        installAdopt(id, lineage)
    else if (configured && !source)
        "echo ${installWord("NM-to-ENST: the table set in nm_to_enst_tsv does not exist: ${configured}")} >&2; exit 1"
    else
        installStart(task) + """
        here=\$(pwd -P)
        source=${installWord(source)}
        {
            printf 'assembly: "%s"\\nvep_installed_version: "%s"\\ncache_root: "%s"\\n' ${installWord(assembly)} ${installWord(release)} "\$here/data"
            [[ -z \$source ]] || printf 'nm_to_enst_tsv: "%s"\\n' "\$here/\$source"
        } > mapping.yaml
        target=\$("\$PRIVA_BIN" installer-prerequisites mapping.yaml --source-root ${root} --stage mapping-target)
        if [[ ! -e \$target ]]; then
            source_id=\$("\$PRIVA_BIN" installer-prerequisites mapping.yaml --source-root ${root} --stage mapping-source)
            workspace=\$("\$PRIVA_BIN" installer-prerequisites mapping.yaml --source-root ${root} --stage mapping-workspace)
            fetched=data/hgnc/\$(basename "\$target")
            mkdir -p data/hgnc
            if [[ \$source_id == ensembl-grch38-mart116-ftp ]]; then
                "\$PRIVA_BIN" biomart fetch --helper-root ${root} --source "\$source_id" --output "\$fetched" --resource-dir "\$workspace"
            else
                mkdir -p "\$workspace"
                "\$PRIVA_BIN" biomart fetch --helper-root ${root} --source "\$source_id" --output "\$fetched" \\
                    --raw-output "\$workspace/\$(basename "\$fetched").biomart.raw.tsv"
            fi
            "\$PRIVA_BIN" biomart verify --helper-root ${root} --mapping "\$fetched" --assembly ${installWord(assembly)}
            printf 'nm_to_enst_tsv: "%s"\\n' "\$here/\$fetched" >> mapping.yaml
        fi
        mapping=\$("\$PRIVA_BIN" installer-prerequisites mapping.yaml --source-root ${root} --stage mapping --print-mapping)
        if [[ -n \$source && \$mapping == "\$here/\$source" ]]; then
            mkdir -p data/hgnc
            for file in "\$source" "\$source.provenance.json"; do
                [[ ! -L \$file ]] || ln -s "\$(readlink "\$file")" "data/hgnc/\$(basename "\$file")"
            done
            mapping="data/hgnc/\$(basename "\$source")"
        fi
        echo "Validated NM-to-ENST mapping selected: \$mapping"
        ${installRecord(id, lineage, receipts, [:])} --key "nm_to_enst_tsv=\$mapping"
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
