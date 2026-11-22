include { installAdopt; installCommit; installRecord; installSetting; installStart; installStub } from './common'

// B1: UTRannotator's uORF tables (github.com/Ensembl/UTRannotator), which
// VEP's UTRAnnotator plugin reads (utr_annotator_file). The step keeps the
// files of the commit it is pinned to (pins.json; the current commit before
// any refresh) and that commit's id, without git's own folder, so one commit
// always gives the same bytes.
process UTR_ANNOTATOR {
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
    def commit = installCommit(planned)
    def table = "data/UTRannotator/uORF_5UTR_${grch}_PUBLIC.txt"
    if (lineage.startsWith('adopted-'))
        installAdopt(id, lineage)
    else
        installStart(task) + """
        [[ -n '${grch}' ]] || { echo 'UTRannotator: unsupported assembly ${assembly}' >&2; exit 1; }
        git clone --quiet https://github.com/Ensembl/UTRannotator data/UTRannotator
        ${commit ? "git -C data/UTRannotator checkout --quiet ${commit}" : ':'}
        git -C data/UTRannotator rev-parse HEAD > data/UTRannotator.commit
        rm -rf data/UTRannotator/.git
        [[ -s ${table} ]] || { echo 'UTRannotator: ${table} is missing or empty' >&2; exit 1; }
        ${installRecord(id, lineage, receipts, [utr_annotator_file: table])}
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
