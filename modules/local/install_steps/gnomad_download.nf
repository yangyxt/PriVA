include { installAdopt; installRecord; installStart; installStub; installWord } from './common'

// C1_chr1 to C1_chrY: one chromosome of gnomAD's v4.1 joint sites VCF and its
// tabix index (gnomad_download_chrN), in GRCh38 coordinates whatever the
// assembly (step C1 lifts them over for GRCh37). The 24 steps share this
// process (registry "process"); the chromosome is the step's id after C1_chr.
// curl downloads the VCF and its index from gnomAD's public bucket; tabix
// writes a new index when the downloaded one is not newer than the VCF, and
// bcftools checks the VCF reads and holds a record (check_vcf_validity).
process GNOMAD_DOWNLOAD {
    tag "${id} ${lineage}"
    label 'install_download'
    storeDir "${params.resource_root}/${id}/${lineage}"
    input:
    tuple val(id), val(lineage), path(plan, stageAs: 'plan.json'), path(receipts, stageAs: 'inputs/*/receipt.json'), val(planned), val(names), path(staged, stageAs: 'staged/*')
    output:
    path('receipt.json'), emit: receipt
    path('data'), emit: data
    script:
    def chr = id.replaceFirst(/^C1_chr/, '')
    def known = (1..22).collect { n -> n.toString() } + ['X', 'Y']
    def url = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz".toString()
    def vcf = "data/gnomAD/gnomad.joint.v4.1.sites.chr${chr}.vcf.bgz".toString()
    if (lineage.startsWith('adopted-'))
        installAdopt(id, lineage)
    else if (!(chr in known))
        "echo 'gnomAD: step ${id} names no chromosome' >&2; exit 1"
    else
        installStart(task) + """
        mkdir -p data/gnomAD
        curl -fL --retry 5 --retry-delay 10 -o ${vcf} ${installWord(url)}
        curl -fL --retry 5 --retry-delay 10 -o ${vcf}.tbi ${installWord(url + '.tbi')}
        [[ ${vcf}.tbi -nt ${vcf} ]] || tabix -f -p vcf ${vcf}
        bcftools view -h ${vcf} > /dev/null
        [[ -n \$(bcftools view -H ${vcf} 2> /dev/null | head -n 1) ]] || { echo 'gnomAD chromosome ${chr}: the download holds no record' >&2; exit 1; }
        ${installRecord(id, lineage, receipts, [("gnomad_download_chr${chr}".toString()): vcf])}
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
