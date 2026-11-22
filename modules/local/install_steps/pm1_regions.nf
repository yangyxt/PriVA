include { installAdopt; installRecord; installSetting; installStart; installStub; installWord } from './common'

// C11: the PM1 hotspot regions (pm1_regions_pkl): the intolerant domains of
// step C10 with the HCSeeker hotspots and the gnomAD regional missense
// constraint (RMC) regions shipped in data/, on transcripts mapped by the
// NM-to-ENST table of step A5. priva installer-prerequisites checks the table,
// and priva mapped-artifact runs priva pm1-regions build and publishes the
// regions with the identity of the table and the inputs they were built from.
// Each input is read where its step keeps it (the target of its staged link).
process PM1_REGIONS {
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
    def release = adopted ? '' : installSetting(planned, 'vep_installed_version')
    def domainsLink = names.all_intolerant_domains ? "staged/${names.all_intolerant_domains}".toString() : ''
    def mappingLink = names.nm_to_enst_tsv ? "staged/${names.nm_to_enst_tsv}".toString() : ''
    def regions = "data/AlphaMissense/pm1_regions.${assembly}.rkyv"
    def root = installWord(projectDir)
    if (adopted)
        installAdopt(id, lineage)
    else if (!assembly)
        "echo 'PM1 regions: unsupported assembly ${configured}' >&2; exit 1"
    else if (!domainsLink || !mappingLink)
        "echo 'PM1 regions: the intolerant domains of step C10 or the NM-to-ENST table of step A5 are missing' >&2; exit 1"
    else
        installStart(task) + """
        here=\$(pwd -P)
        hcseeker=${installWord("${projectDir}/data/HCSeeker/HC_spots.${assembly}.tsv")}
        rmc=${installWord("${projectDir}/data/gnomAD_RMC/gnomad_v2.1.1_rmc_flat.${assembly}.tsv")}
        mkdir -p data/AlphaMissense
        printf 'assembly: "%s"\\nvep_installed_version: "%s"\\ncache_root: "%s"\\ntmp_dir: "%s"\\nall_intolerant_domains: "%s"\\nnm_to_enst_tsv: "%s"\\n' \\
            ${installWord(configured)} ${installWord(release)} "\$here/data" "\$here" \\
            "\$(readlink -f ${installWord(domainsLink)})" "\$(readlink -f ${installWord(mappingLink)})" > pm1.yaml
        # The DAS set as a native archive (a Python DAS pickle is converted once).
        domains=\$("\$PRIVA_BIN" cache-install resolve --key all_intolerant_domains --build-config pm1.yaml --source-root ${root})
        mapping=\$("\$PRIVA_BIN" installer-prerequisites pm1.yaml --source-root ${root} --stage mapping --print-mapping)
        identity=\$("\$PRIVA_BIN" biomart identity --helper-root ${root} --mapping "\$mapping" --input "\$rmc" --input "\$hcseeker")
        "\$PRIVA_BIN" mapped-artifact run --artifact "\$here/${regions}" --expected-identity "\$identity" \\
            --config "\$here/pm1.yaml" --native-producer pm1-regions --helper-root ${root} \\
            --mapping "\$mapping" --kind pm1 --assembly ${assembly} --cache-root "\$here/data" \\
            --output-dir "\$here/data/AlphaMissense" --input "\$rmc" --input "\$hcseeker" -- \\
            --assembly ${assembly} --nm-to-enst "\$mapping" --intolerant-domains "\$domains" \\
            --hcseeker "\$hcseeker" --rmc "\$rmc" --output "\$here/${regions}"
        echo "PM1 regions saved to ${regions}"
        ${installRecord(id, lineage, receipts, [pm1_regions_pkl: regions])}
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
