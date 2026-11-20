include { installAdopt; installRecord; installSetting; installStart; installStub; installWord } from './common'

// A1: the reference FASTA of the assembly and its samtools index (ref_genome).
// priva reference-resources install downloads UCSC's pinned FASTA with curl,
// checks its MD5, decompresses it and indexes it with samtools faidx, and
// checks the index names chromosome 1 with the assembly's length; it prints the
// path it chose. A FASTA you set in ref_genome is checked the same way and used
// where it is; when its index is missing or older than it, samtools writes a
// new index beside the staged link, and that index is kept with a link to your
// FASTA. The source receipt beside a downloaded FASTA names its path, so record
// rewrites it to where storeDir moves the FASTA.
process REFERENCE_FASTA {
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
    def source = names.ref_genome ? "staged/${names.ref_genome}".toString() : ''
    def configured = planned.user_paths?.ref_genome ?: ''
    if (lineage.startsWith('adopted-'))
        installAdopt(id, lineage)
    else if (configured && !source)
        "echo ${installWord("Reference: the FASTA set in ref_genome does not exist: ${configured}")} >&2; exit 1"
    else
        installStart(task) + """
        here=\$(pwd -P)
        source=${installWord(source)}
        supplied=\${source:+\$here/\$source}
        printf 'assembly: "%s"\\ncache_root: "%s"\\nref_genome: "%s"\\n' ${installWord(assembly)} "\$here/data" "\$supplied" > reference.yaml
        "\$PRIVA_BIN" reference-resources install reference.yaml --data-dir ${installWord("${projectDir}/data")} --scope selected --use-cache-root > selected.tsv
        selected=\$(awk -F '\\t' '\$1 == "ref_genome" {print \$2}' selected.tsv)
        relocate=()
        if [[ \$selected == "\$here"/data/* ]]; then
            fasta=\${selected#"\$here"/}
            relocate=(--relocate "\$fasta.source.json")
        elif [[ -n \$supplied && \$selected == "\$supplied" ]]; then
            fasta=\$source
            for index in "\$source.fai" "\$source.gzi"; do
                [[ ! -e \$index || -L \$index ]] || fasta=
            done
            if [[ -z \$fasta ]]; then
                mkdir -p data/reference
                for file in "\$source" "\$source.fai" "\$source.gzi"; do
                    if [[ -L \$file ]]; then
                        ln -s "\$(readlink "\$file")" "data/reference/\$(basename "\$file")"
                    elif [[ -e \$file ]]; then
                        mv "\$file" data/reference/
                    fi
                done
                fasta="data/reference/\$(basename "\$source")"
            fi
        else
            echo "Reference: unexpected FASTA \$selected" >&2; exit 1
        fi
        ${installRecord(id, lineage, receipts, [:])} --key "ref_genome=\$fasta" "\${relocate[@]}"
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
