include { installAdopt; installRecord; installSetting; installStart; installStub; installWord } from './common'

// B4: the SpliceAI scores of VEP's SpliceAI plugin (spliceai_snv_prescore and
// spliceai_indel_prescore), which you download yourself (Illumina's terms):
// SpliceAI v1 VCFs for GRCh37, SpliceAI2 v2.0 for GRCh38. priva
// spliceai-cache-input install reads the bounded header of each and checks
// them as a pair: a VCF (SpliceAI v1 with its tabix index, or a prepared
// SpliceAI2 VCF with its adapter receipt) is recorded where it is; a pair of
// native SpliceAI2 tables with their downloads-complete manifest is converted
// into data/VEP_plugin_cache/SpliceAI2/GRCh38 (priva spliceai2-install check,
// spliceai2-prepare prepare, spliceai2-install register). Without the two
// files the step stops with where to get them; it asks nothing (the shell
// function prompted).
process SPLICEAI {
    tag "${id} ${lineage}"
    label 'install'
    storeDir "${params.resource_root}/${id}/${lineage}"
    input:
    tuple val(id), val(lineage), path(plan, stageAs: 'plan.json'), path(receipts, stageAs: 'inputs/*/receipt.json'), val(planned), val(names), path(staged, stageAs: 'staged/*')
    output:
    path('receipt.json'), emit: receipt
    path('data'), emit: data
    script:
    def adopted = lineage.startsWith('adopted-')
    def assembly = adopted ? '' : installSetting(planned, 'assembly')
    def snvLink = names.spliceai_snv_prescore ? "staged/${names.spliceai_snv_prescore}".toString() : ''
    def indelLink = names.spliceai_indel_prescore ? "staged/${names.spliceai_indel_prescore}".toString() : ''
    def grch37 = assembly in ['hg19', 'GRCh37']
    def guidance = grch37
        ? 'GRCh37 installation uses SpliceAI v1. Download its assembly-matched SNV and indel VCF caches through your authorized Illumina BaseSpace account (https://basespace.illumina.com/s/otSPW8hnhaZR), index each completed compressed VCF with tabix -p vcf, and set spliceai_snv_prescore and spliceai_indel_prescore to their absolute paths.'
        : 'GRCh38 installation uses SpliceAI2 v2.0. Complete the publisher-authorized downloads of all native score shards and gene-context resources (https://huggingface.co/datasets/illumina-ai/SpliceAI2-data/tree/cd1b8b82bcef9c97ac3b3b9a24f0de3a5bbf77f9), and set spliceai_snv_prescore and spliceai_indel_prescore to an SNV and an indel input of the complete cache.'
    if (adopted)
        installAdopt(id, lineage)
    else if (!(assembly in ['hg19', 'GRCh37', 'hg38', 'GRCh38']))
        "echo 'SpliceAI: not supported assembly version: ${assembly}' >&2; exit 1"
    else if (!snvLink || !indelLink)
        "echo ${installWord('SpliceAI: ' + guidance)} >&2; exit 1"
    else
        installStart(task) + """
        here=\$(pwd -P)
        # The files as you set them (the target of each staged link).
        snv=\$(readlink ${installWord(snvLink)})
        indel=\$(readlink ${installWord(indelLink)})
        printf 'assembly: %s\\nspliceai_snv_prescore: "%s"\\nspliceai_indel_prescore: "%s"\\n' \\
            ${installWord(assembly)} "\$snv" "\$indel" > spliceai.yaml
        mkdir -p data/VEP_plugin_cache
        "\$PRIVA_BIN" spliceai-cache-input install --config "\$here/spliceai.yaml" \\
            --cache-dir "\$here/data/VEP_plugin_cache" --snv "\$snv" --indel "\$indel"
        # The two paths it registered: the files you set, or the prepared VCFs.
        registered() { sed -n "s/^\$1: *//p" spliceai.yaml | sed -e 's/^"//' -e 's/"\$//'; }
        snv=\$(registered spliceai_snv_prescore)
        indel=\$(registered spliceai_indel_prescore)
        relocate=()
        receipt=data/VEP_plugin_cache/SpliceAI2/GRCh38/adapter.ready.json
        [[ ! -f \$receipt ]] || relocate=(--relocate \$receipt)
        ${installRecord(id, lineage, receipts, [:])} --key "spliceai_snv_prescore=\${snv#\$here/}" \\
            --key "spliceai_indel_prescore=\${indel#\$here/}" "\${relocate[@]}"
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
