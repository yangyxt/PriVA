include { installAdopt; installRecord; installSetting; installStart; installStub; installValue; installWord } from './common'

// B2: the GERP conservation bigWig for VEP's Conservation plugin
// (conservation_file): Ensembl's 91-mammal GERP scores for GRCh38
// (conservation_GRCh38_url) or UCSC's All_hg19_RS for GRCh37, downloaded with
// curl. priva conservation-resource publish reads it with the plugin's Perl
// bigWig reader (beside the python on PATH), checks its chromosome lengths
// against the reference of step A1 and moves it into place.
process CONSERVATION {
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
    def grch = assembly in ['hg19', 'GRCh37'] ? 'GRCh37' : (assembly in ['hg38', 'GRCh38'] ? 'GRCh38' : '')
    def url = grch == 'GRCh37' ? 'https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/All_hg19_RS.bw'
        : installValue(planned, 'conservation_GRCh38_url', 'https://ftp.ensembl.org/pub/release-114/compara/conservation_scores/91_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw')
    def target = "data/Conservation/gerp_conservation_scores.homo_sapiens.${grch}.bw"
    if (lineage.startsWith('adopted-'))
        installAdopt(id, lineage)
    else if (!grch)
        "echo 'Conservation: unsupported assembly ${assembly}' >&2; exit 1"
    else if (!(url ==~ /^https?:\/\/.*/))
        "echo 'Conservation: conservation_GRCh38_url is not an http(s) URL' >&2; exit 1"
    else
        installStart(task) + """
        [[ -s staged/${names.ref_genome} ]] || { echo 'Conservation: the reference of step A1 is missing' >&2; exit 1; }
        python=\$(type -P python) || { echo 'Conservation: no python on PATH (the Perl bigWig reader is looked for beside it)' >&2; exit 1; }
        mkdir -p data/Conservation
        curl -fL --retry 5 --retry-delay 10 -o ${target}.partial ${installWord(url)}
        "\$PRIVA_BIN" conservation-resource publish --input ${target}.partial --assembly ${grch} \\
            --reference staged/${names.ref_genome} --tool-prefix "\$(dirname "\$python")" \\
            --source-url ${installWord(url)} --publish-to ${target}
        ${installRecord(id, lineage, receipts, [conservation_file: target])}
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
