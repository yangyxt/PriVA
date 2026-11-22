include { installAdopt; installRecord; installSetting; installStart; installStub; installWord } from './common'

// D0: the latest official ClinVar release (clinvar_release_inputs, the file
// that names its downloads): NCBI's weekly ClinVar VCF of the assembly and its
// index, variant_summary and MANE's summary, from the sources of
// cache_resources.json. priva clinvar-gene-refresh download reads the
// release's ##fileDate from a bounded header read, downloads each file with
// curl, checks the publisher MD5 stayed the same around the downloads, and
// writes release/inputs.ready.json; step D1 reads that file and downloads
// nothing. The paths it names are rewritten to where storeDir keeps them.
process CLINVAR_DOWNLOAD {
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
    def release = 'data/ClinVar/release'
    if (adopted)
        installAdopt(id, lineage)
    else
        installStart(task) + """
        here=\$(pwd -P)
        printf 'assembly: %s\\n' ${installWord(assembly)} > d0.yaml
        "\$PRIVA_BIN" clinvar-gene-refresh download --config d0.yaml --source-root ${installWord(projectDir)} \\
            --output "\$here/${release}"
        relocate=(--relocate ${release}/inputs.ready.json)
        for written in ${release}/inputs/*.source.json; do
            relocate+=(--relocate "\$written")
        done
        ${installRecord(id, lineage, receipts, [clinvar_release_inputs: "${release}/inputs.ready.json".toString()])} "\${relocate[@]}"
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
