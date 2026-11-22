include { installAdopt; installRecord; installSetting; installStart; installStub; installValue; installWord } from './common'

// C10: the AlphaMissense-intolerant protein domains (DAS): the domain
// tolerance table (alphamissense_intolerant_domains), the transcript-to-domain
// map (alphamissense_tranx_domain_map) and the set of intolerant domains
// (all_intolerant_domains). From the domain statistics of step C6, the
// InterPro mapping of step C8, the AlphaMissense VCF of step C4 and the ClinVar
// VCF of step D1: priva alpha-exon-history builds the exon history per
// ##contig and merges it under priva_python_hash_seed, priva domain-metadata
// reads the domain positions, priva das-domains transcripts writes the map and
// das-domains analyze tests each domain against the functional references and
// filters them by HCSeeker coldspots and benign ClinVar variants. Each input is
// read where its step keeps it (the target of its staged link), as the
// builders record the paths they read.
process DAS_DOMAINS {
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
    def seed = adopted ? '' : installValue(planned, 'priva_python_hash_seed', '0')
    def release = adopted ? '' : installSetting(planned, 'vep_installed_version')
    def statisticsLink = names.alphamissense_pd_stat ? "staged/${names.alphamissense_pd_stat}".toString() : ''
    def interproLink = names.interpro_mapping_pickle ? "staged/${names.interpro_mapping_pickle}".toString() : ''
    def vcfLink = names.alphamissense_vep_vcf ? "staged/${names.alphamissense_vep_vcf}".toString() : ''
    def clinvarLink = names.clinvar_vcf ? "staged/${names.clinvar_vcf}".toString() : ''
    def dir = 'data/AlphaMissense'
    def tolerance = "${dir}/domain_tolerance_analysis.${assembly}.tsv"
    def map = "${dir}/transcript_pp_domain_mapping.${assembly}.rkyv"
    def intolerant = "${dir}/all_intolerant_domains.${assembly}.rkyv"
    def root = installWord(projectDir)
    if (adopted)
        installAdopt(id, lineage)
    else if (!assembly)
        "echo 'DAS: unsupported assembly ${configured}' >&2; exit 1"
    else if (!(seed ==~ /[0-9]+/))
        "echo 'ERROR: priva_python_hash_seed must be a non-negative integer: ${seed}' >&2; exit 1"
    else if (!statisticsLink || !interproLink || !vcfLink)
        "echo 'DAS: the statistics of step C6, the InterPro mapping of step C8 or the AlphaMissense VCF of step C4 are missing' >&2; exit 1"
    else
        installStart(task) + """
        here=\$(pwd -P)
        statistics=\$(readlink -f ${installWord(statisticsLink)})
        vcf=\$(readlink -f ${installWord(vcfLink)})
        clinvar=${installWord(clinvarLink)}
        [[ -z \$clinvar ]] || clinvar=\$(readlink -f "\$clinvar")
        hcseeker=${installWord("${projectDir}/data/HCSeeker/HC_spots.${assembly}.tsv")}
        normalization=${installWord("${projectDir}/data/InterPro/InterPro_domain_norm_map.json")}
        mkdir -p ${dir}
        # The InterPro mapping as a native archive: the repository's own, or
        # the mapping of step C8 converted once into this step's folder.
        printf 'assembly: "%s"\\nvep_installed_version: "%s"\\ncache_root: "%s"\\ntmp_dir: "%s"\\ninterpro_mapping_pickle: "%s"\\n' \\
            ${installWord(configured)} ${installWord(release)} "\$here/data" "\$here" "\$(readlink -f ${installWord(interproLink)})" > das.yaml
        interpro=\$("\$PRIVA_BIN" cache-install resolve --key interpro_mapping_pickle --build-config das.yaml --source-root ${root})
        filters=()
        if [[ -f \$vcf && -f \$hcseeker ]]; then
            filters+=(--metadata ${dir}/domain_metadata.${assembly}.parquet --spots "\$hcseeker")
            echo "Filter 1 (coldspot overlap) enabled: am_vcf=\$vcf, hcseeker=\$hcseeker"
        else
            echo 'Filter 1 (coldspot overlap) skipped: am_vcf or hcseeker_spots_tsv not found'
        fi
        if [[ -n \$clinvar && -f \$clinvar ]]; then
            filters+=(--clinvar-vcf "\$clinvar")
            echo "Filter 2 (benign ClinVar) enabled: clinvar_vcf=\$clinvar"
        else
            echo 'Filter 2 (benign ClinVar) skipped: clinvar_vcf not found'
        fi
        # Exon history: one shard per ##contig of the VCF, merged in header
        # order under the seed, and the serialized-set export transcripts reads.
        mapfile -t contigs < <(bcftools view -h "\$vcf" | sed -n 's/^##contig=<ID=\\([^,>]*\\).*/\\1/p')
        [[ \${#contigs[@]} -gt 0 ]] || { echo "ERROR: \$vcf declares no ##contig lines for the chromosome builders" >&2; exit 1; }
        history=${dir}/\$(basename "\$statistics" .parquet).exon_history.seed${seed}.parquet
        shards=\$history.chromosomes
        mkdir -p "\$shards"
        inputs=()
        for index in "\${!contigs[@]}"; do
            shard=\$shards/\$(printf '%05d' "\$index").parquet
            "\$PRIVA_BIN" alpha-exon-history chromosome --vcf "\$vcf" --chrom "\${contigs[\$index]}" --output "\$shard"
            inputs+=(--input "\$shard")
        done
        "\$PRIVA_BIN" alpha-exon-history merge "\${inputs[@]}" --hash-seed ${seed} --output "\$history"
        rm -rf "\$shards"
        "\$PRIVA_BIN" alpha-exon-history export --input "\$history" --output "\${history%.parquet}.lines"
        if [[ -f \$vcf && -f \$hcseeker ]]; then
            "\$PRIVA_BIN" domain-metadata build --vcf "\$vcf" --output ${dir}/domain_metadata.${assembly}.parquet \\
                --assembly ${assembly} --threads ${task.cpus}
        fi
        identity() { [[ \$1 =~ \\"builder_source_sha256\\":\\"([0-9a-f]{64})\\" ]] && printf '%s\\n' "\${BASH_REMATCH[1]}"; }
        statistics_builder=\$(identity "\$("\$PRIVA_BIN" alphamissense-domain identity --input "\$statistics")")
        history_builder=\$(identity "\$("\$PRIVA_BIN" alpha-exon-history identity --input "\$history")")
        "\$PRIVA_BIN" das-domains transcripts --statistics "\$statistics" --exon-history "\${history%.parquet}.lines" \\
            --source-builder-sha256 "\$statistics_builder" --python-hash-seed ${seed} --output ${map} \\
            --history-builder-sha256 "\$history_builder"
        "\$PRIVA_BIN" das-domains analyze --statistics "\$statistics" --source-builder-sha256 "\$statistics_builder" \\
            --interpro "\$interpro" --normalization-json "\$normalization" --python-hash-seed ${seed} \\
            --threads ${task.cpus} --assembly ${assembly} --output-dir ${dir} "\${filters[@]}"
        for produced in ${tolerance} ${map} ${intolerant}; do
            [[ -s \$produced ]] || { echo "ERROR: Missing AlphaMissense domain output: \$produced" >&2; exit 1; }
        done
        ${installRecord(id, lineage, receipts, [alphamissense_intolerant_domains: tolerance, alphamissense_tranx_domain_map: map, all_intolerant_domains: intolerant])}
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
