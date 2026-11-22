include { installAdopt; installRecord; installSetting; installStart; installStub; installWord } from './common'

// D8: the PM5 index of the VEP-annotated ClinVar VCF of step D1
// (clinvar_pm5_index): each confident pathogenic ClinVar protein change with
// its consequence, NMD status and residue, which priva acmg-assign reads for
// PM5. priva clinvar-contract checks the VCF's INFO and CSQ fields, and priva
// clinvar-pm5-index writes the index, named after the managed ClinVar VCF
// (clinvar.hg19.vep.vcf.gz or clinvar.hg38.vep.vcf.gz), as the installer
// named it.
process CLINVAR_PM5 {
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
    def index = "data/ClinVar/clinvar.${assembly}.vep.pm5_index.tsv.gz"
    if (adopted)
        installAdopt(id, lineage)
    else if (!assembly)
        "echo 'ClinVar PM5: unsupported assembly ${configured}' >&2; exit 1"
    else if (!name)
        "echo 'ClinVar PM5: step D1 records no clinvar_vcf' >&2; exit 1"
    else
        installStart(task) + """
        vcf=staged/${installWord(name)}
        [[ -s \$vcf ]] || { echo "ERROR: The VEP-annotated ClinVar VCF is missing: \$vcf" >&2; exit 1; }
        "\$PRIVA_BIN" clinvar-contract header --vcf "\$vcf" --stage vep
        mkdir -p data/ClinVar
        "\$PRIVA_BIN" clinvar-pm5-index --vcf "\$vcf" --output ${index}
        ${installRecord(id, lineage, receipts, [clinvar_pm5_index: index])}
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
