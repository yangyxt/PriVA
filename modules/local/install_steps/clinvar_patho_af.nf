include { installAdopt; installRecord; installSetting; installStart; installStub; installWord } from './common'

// D6: the allele frequencies of pathogenic ClinVar variants per gene
// (clinvar_patho_af_stat) and per exon (clinvar_patho_exon_af_stat), from the
// VEP-annotated ClinVar VCF of step D1. priva clinvar-contract checks the
// VCF's INFO and CSQ fields, priva clinvar-gene-af and clinvar-exon-af write
// the two archives, and clinvar-contract checks each archive. They are named
// after the managed ClinVar VCF (clinvar.hg19.vep.vcf.gz or
// clinvar.hg38.vep.vcf.gz), as the installer named them.
process CLINVAR_PATHO_AF {
    tag "${id} ${lineage}"
    label 'install'
    storeDir "${params.resource_root}/${id}/${lineage}"
    input:
    tuple val(id), val(lineage), path(plan, stageAs: 'plan.json'), path(receipts, stageAs: 'inputs/*/receipt.json'), val(planned), val(names), path(staged, stageAs: 'staged/*')
    output:
    path('receipt.json'), emit: receipt
    path('data'), emit: data
    script:
    def name = (names.clinvar_vcf ?: '').toString()
    def adopted = lineage.startsWith('adopted-')
    def configured = adopted ? '' : installSetting(planned, 'assembly')
    def assembly = configured in ['hg19', 'GRCh37'] ? 'hg19' : (configured in ['hg38', 'GRCh38'] ? 'hg38' : '')
    def gene = "data/ClinVar/clinvar.${assembly}.vep.patho_af_stat.rkyv"
    def exon = "data/ClinVar/clinvar.${assembly}.vep.patho_exon_af_stat.rkyv"
    if (adopted)
        installAdopt(id, lineage)
    else if (!assembly)
        "echo 'ClinVar AF: unsupported assembly ${configured}' >&2; exit 1"
    else if (!name)
        "echo 'ClinVar AF: step D1 records no clinvar_vcf' >&2; exit 1"
    else
        installStart(task) + """
        vcf=staged/${installWord(name)}
        "\$PRIVA_BIN" clinvar-contract af-paths-check --vcf "\$vcf" --config plan.json --gene ${gene} --exon ${exon}
        "\$PRIVA_BIN" clinvar-contract header --vcf "\$vcf" --stage af
        mkdir -p data/ClinVar
        # Without a CSQ format the builders succeed without writing, so an
        # absent archive is a failure here.
        "\$PRIVA_BIN" clinvar-gene-af build --vcf "\$vcf" --output ${gene}.partial --threads ${task.cpus}
        [[ -s ${gene}.partial ]] || { echo "ERROR: priva clinvar-gene-af wrote no archive for \$vcf" >&2; exit 1; }
        "\$PRIVA_BIN" clinvar-contract cache --input ${gene}.partial --kind gene-af
        "\$PRIVA_BIN" clinvar-exon-af build --vcf "\$vcf" --output ${exon}.partial --threads ${task.cpus}
        [[ -s ${exon}.partial ]] || { echo "ERROR: priva clinvar-exon-af wrote no archive for \$vcf" >&2; exit 1; }
        "\$PRIVA_BIN" clinvar-contract cache --input ${exon}.partial --kind exon-af
        mv ${gene}.partial ${gene}
        mv ${exon}.partial ${exon}
        ${installRecord(id, lineage, receipts, [clinvar_patho_af_stat: gene, clinvar_patho_exon_af_stat: exon])}
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
