include { installAdopt; installRecord; installSetting; installStart; installStub; installWord } from './common'

// C6: the AlphaMissense score statistics of each protein domain
// (alphamissense_pd_stat), from the VEP- and gnomAD-annotated AlphaMissense
// VCF of step C4. priva alphamissense-domain builds one shard per ##contig of
// the VCF (bcftools reads the header), in header order, merges them, and
// checks the merged statistics against the VCF. The statistics record the
// VCF's path, so the builders read it where step C4 keeps it (the target of
// its staged link), not through the link in this task's folder.
process ALPHAMISSENSE_DOMAIN {
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
    def configured = adopted ? '' : installSetting(planned, 'assembly')
    def assembly = configured in ['hg19', 'GRCh37'] ? 'hg19' : (configured in ['hg38', 'GRCh38'] ? 'hg38' : '')
    def statistics = "data/AlphaMissense/AlphaMissense_${assembly}.vep.prot.domain.stats.parquet"
    def vcf = names.alphamissense_vep_vcf ? "staged/${names.alphamissense_vep_vcf}".toString() : ''
    if (adopted)
        installAdopt(id, lineage)
    else if (!assembly)
        "echo 'AlphaMissense domains: unsupported assembly ${configured}' >&2; exit 1"
    else if (!vcf)
        "echo 'AlphaMissense domains: step C4 records no alphamissense_vep_vcf' >&2; exit 1"
    else
        installStart(task) + """
        vcf=\$(readlink -f ${installWord(vcf)})
        mapfile -t contigs < <(bcftools view -h "\$vcf" | sed -n 's/^##contig=<ID=\\([^,>]*\\).*/\\1/p')
        [[ \${#contigs[@]} -gt 0 ]] || { echo "ERROR: \$vcf declares no ##contig lines for the chromosome builders" >&2; exit 1; }
        shards=${statistics}.chromosomes
        mkdir -p "\$shards"
        inputs=()
        for index in "\${!contigs[@]}"; do
            shard=\$shards/\$(printf '%05d' "\$index").parquet
            "\$PRIVA_BIN" alphamissense-domain chromosome --vcf "\$vcf" --chrom "\${contigs[\$index]}" --output "\$shard"
            inputs+=(--input "\$shard")
        done
        "\$PRIVA_BIN" alphamissense-domain merge "\${inputs[@]}" --output ${statistics}
        rm -rf "\$shards"
        "\$PRIVA_BIN" alphamissense-domain check --input ${statistics} --vcf "\$vcf" > /dev/null
        ${installRecord(id, lineage, receipts, [alphamissense_pd_stat: statistics])}
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
