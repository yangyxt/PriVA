include { installAdopt; installRecord; installSetting; installStart; installStub; installValue; installWord } from './common'

// B9: gnomAD's pext bigWigs from UCSC for the annotation's pext columns: the
// mean over tissues (pext_mean_bw) and the tissues listed in pext_tissues
// (pext_tissue_bw_dir/tissues/NAME.bw), downloaded with curl. With
// pext_enabled false the step records neither. A tissue UCSC does not publish
// for the assembly stops the step with its name; the names differ between hg19
// and hg38 (some hg19 names end with an underscore).
process PEXT {
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
    def ucsc = assembly in ['hg19', 'GRCh37'] ? 'hg19' : (assembly in ['hg38', 'GRCh38'] ? 'hg38' : '')
    def enabled = installValue(planned, 'pext_enabled', 'false') == 'true'
    def tissues = installValue(planned, 'pext_tissues', '').tokenize(',').collect { t -> t.trim() }.findAll { t -> t }
    def base = "https://hgdownload.soe.ucsc.edu/gbdb/${ucsc}/gnomAD/pext"
    def mean = ucsc == 'hg38' ? 'exp_prop_mean.bw' : 'mean_proportion.bw'
    def folder = "data/pext/${ucsc}"
    if (lineage.startsWith('adopted-'))
        installAdopt(id, lineage)
    else if (!enabled)
        installStart(task) + """
        echo 'pext annotation is disabled (pext_enabled); nothing to download'
        ${installRecord(id, lineage, receipts, [:])}
        """.stripIndent()
    else if (!ucsc)
        "echo 'pext: unsupported assembly ${assembly}' >&2; exit 1"
    else
        installStart(task) + """
        mkdir -p ${folder}/tissues
        curl -fL --retry 5 --retry-delay 10 -o ${folder}/${mean} ${installWord("${base}/${mean}")}
        for tissue in ${tissues.collect { t -> installWord(t) }.join(' ')}; do
            curl -fL --retry 5 --retry-delay 10 -o "${folder}/tissues/\$tissue.bw" "${base}/\$tissue.bw" || {
                echo "pext: UCSC has no ${ucsc} tissue \$tissue (pext_tissues); see ${base}/" >&2; exit 1; }
        done
        ${installRecord(id, lineage, receipts, [pext_mean_bw: "${folder}/${mean}", pext_tissue_bw_dir: folder])}
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
