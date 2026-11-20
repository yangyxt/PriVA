include { installAdopt; installRecord; installSetting; installStart; installStub; installWord } from './common'

// B5: Ensembl's SpliceVault table for VEP's SpliceVault plugin and its tabix
// index (splicevault_prescore). Ensembl publishes GRCh38 only: for GRCh37 the
// table is lifted over with priva splicevault-convert (CrossMap, the chain
// hg38_hg19_chain: a file shipped in data/liftover/ or an absolute path, and
// the reference of step A1), sorted, compressed with bgzip and indexed with
// tabix, as the installer did. priva checks the index against the plugin's
// columns.
process SPLICEVAULT {
    tag "${id} ${lineage}"
    label 'install'
    storeDir "${params.resource_root}/${id}/${lineage}"
    input:
    tuple val(id), val(lineage), path(plan, stageAs: 'plan.json'), path(receipts, stageAs: 'inputs/*/receipt.json'), val(planned), val(names), path(staged, stageAs: 'staged/*')
    output:
    path('receipt.json'), emit: receipt
    path('data'), emit: data
    script:
    def assembly = installSetting(planned, 'assembly')
    def url = installSetting(planned, 'splicevault_url')
    def selected = planned.settings?.hg38_hg19_chain ? planned.settings.hg38_hg19_chain.toString() : ''
    def chain = selected in ['', 'null', 'hg38ToHg19.over.chain.gz'] ? "${projectDir}/data/liftover/hg38ToHg19.over.chain.gz".toString()
        : (selected.startsWith('/') ? selected : "${projectDir}/data/liftover/${selected}".toString())
    def grch38 = 'data/SpliceVault/SpliceVault_data_GRCh38.tsv.gz'
    def grch37 = 'data/SpliceVault/SpliceVault_data_GRCh37.tsv.gz'
    if (lineage.startsWith('adopted-'))
        installAdopt(id, lineage)
    else if (assembly in ['hg38', 'GRCh38'])
        installStart(task) + """
        mkdir -p data/SpliceVault
        curl -fL --retry 5 --retry-delay 10 -o ${grch38} ${installWord(url)}
        curl -fL --retry 5 --retry-delay 10 -o ${grch38}.tbi ${installWord(url + '.tbi')}
        bgzip -t ${grch38}
        tabix -l ${grch38} > /dev/null
        "\$PRIVA_BIN" indexed-plugin-cache ensure-index --input ${grch38} --kind splicevault
        ${installRecord(id, lineage, receipts, [splicevault_prescore: grch38])}
        """.stripIndent()
    else if (assembly in ['hg19', 'GRCh37'])
        installStart(task) + """
        [[ -s ${installWord(chain)} ]] || { echo 'SpliceVault: the chain hg38_hg19_chain is missing or empty: ${chain}' >&2; exit 1; }
        [[ -s staged/${names.ref_genome} ]] || { echo 'SpliceVault: the reference of step A1 is missing' >&2; exit 1; }
        mkdir -p data/SpliceVault tmp
        curl -fL --retry 5 --retry-delay 10 -o SpliceVault_data_GRCh38.tsv.gz ${installWord(url)}
        bgzip -t SpliceVault_data_GRCh38.tsv.gz
        gunzip -c SpliceVault_data_GRCh38.tsv.gz > SpliceVault_data_GRCh38.tsv
        "\$PRIVA_BIN" splicevault-convert --input_tsv SpliceVault_data_GRCh38.tsv --output_tsv SpliceVault_data_GRCh37.tsv \\
            --chain_file ${installWord(chain)} --reference_fasta staged/${names.ref_genome} --temp_dir tmp
        head -1 SpliceVault_data_GRCh37.tsv > header.tsv
        tail -n +2 SpliceVault_data_GRCh37.tsv > body.tsv
        LC_ALL=C sort -k1,1V -k2,2n --parallel=8 --buffer-size=8G -T tmp body.tsv > sorted.tsv
        cat header.tsv sorted.tsv | bgzip > ${grch37}
        bgzip -t ${grch37}
        tabix -f -s 1 -b 2 -e 2 -c '#' ${grch37}
        tabix -l ${grch37} > /dev/null
        "\$PRIVA_BIN" indexed-plugin-cache ensure-index --input ${grch37} --kind splicevault
        ${installRecord(id, lineage, receipts, [splicevault_prescore: grch37])}
        """.stripIndent()
    else
        "echo 'SpliceVault: unsupported assembly ${assembly}' >&2; exit 1"
    stub:
    installStub(id, lineage, receipts, task)
}
