include { installAdopt; installRecord; installSetting; installStart; installStub; installWord } from './common'

// A4: Ensembl's CDS FASTA of the VEP release (vep_installed_version), for
// the alternative start codon check of start_lost variants (cds_fasta_file).
// The release is pinned in the URL, so one release always gives the same
// sequences. priva checks the FASTA records and that the headers name the
// assembly. The Python detector's pickle is no longer written: the analysis
// reads the FASTA itself (rust/src/acmg_step/alternative_start_codon.rs).
process CDS_FASTA {
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
    def release = installSetting(planned, 'vep_installed_version')
    def grch = assembly in ['hg19', 'GRCh37'] ? 'GRCh37' : (assembly in ['hg38', 'GRCh38'] ? 'GRCh38' : '')
    def url = grch == 'GRCh37'
        ? "https://ftp.ensembl.org/pub/grch37/release-${release}/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.cds.all.fa.gz"
        : "https://ftp.ensembl.org/pub/release-${release}/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz"
    def fasta = "data/Ensembl_CDS/release-${release}/Homo_sapiens.${grch}.cds.all.fa.gz"
    if (lineage.startsWith('adopted-'))
        installAdopt(id, lineage)
    else
        installStart(task) + """
        [[ -n '${grch}' ]] || { echo 'CDS FASTA: unsupported assembly ${assembly}' >&2; exit 1; }
        [[ ${installWord(release)} =~ ^[1-9][0-9]*\$ ]] || { echo 'CDS FASTA: vep_installed_version must be a release number' >&2; exit 1; }
        mkdir -p data/Ensembl_CDS/release-${release}
        curl -fL --retry 5 --retry-delay 10 -o ${fasta} ${installWord(url)}
        "\$PRIVA_BIN" cds-fasta validate --fasta ${fasta} --assembly ${installWord(assembly)}
        ${installRecord(id, lineage, receipts, [cds_fasta_file: fasta])}
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
