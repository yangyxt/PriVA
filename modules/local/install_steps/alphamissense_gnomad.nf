include { installAdopt; installRecord; installSetting; installStart; installStub; installWord } from './common'

// C4: the VEP-annotated AlphaMissense VCF of step C3 with gnomAD's joint
// frequencies (alphamissense_vep_vcf), from the chromosome family of step C1
// (gnomad_vcf_chrX names its chrX member; the others are beside it). The step
// copies C3's VCF and index into its own folder and runs the annotation's
// gnomAD step on the copy (priva task step --step gnomad, the commands of
// anno_agg_gnomAD_data: split per chromosome, bcftools annotate from that
// chromosome's gnomAD VCF, chrM kept without values, concat and sort).
process ALPHAMISSENSE_GNOMAD {
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
    def assembly = adopted ? '' : installSetting(planned, 'assembly')
    def annotated = names.alphamissense_vep_vcf ? "staged/${names.alphamissense_vep_vcf}".toString() : ''
    def gnomad = names.gnomad_vcf_chrX ? "staged/${names.gnomad_vcf_chrX}".toString() : ''
    def folder = 'data/VEP_plugin_cache/AlphaMissense'
    if (adopted)
        installAdopt(id, lineage)
    else if (!annotated || !gnomad)
        "echo 'AlphaMissense gnomAD: the annotation of step C3 or the gnomAD family of step C1 is missing' >&2; exit 1"
    else
        installStart(task) + """
        here=\$(pwd -P)
        source=\$(readlink -f ${installWord(annotated)})
        # C1's chrX member where C1 keeps it, beside the other chromosomes.
        chrx=\$(readlink ${installWord(gnomad)})
        [[ -s \$source && -s \$source.tbi ]] || { echo "AlphaMissense gnomAD: \$source or its index is missing" >&2; exit 1; }
        mkdir -p ${folder}
        target=${folder}/\$(basename "\$source")
        cp -p "\$source" "\$target"
        cp -p "\$source.tbi" "\$target.tbi"
        chmod u+w "\$target" "\$target.tbi"
        printf 'assembly: %s\\nthreads: %s\\ngnomad_vcf_chrX: "%s"\\n' ${installWord(assembly)} ${task.cpus} "\$chrx" > c4.yaml
        "\$PRIVA_BIN" task step --step gnomad --root ${installWord(projectDir)} --config "\$here/c4.yaml" --vcf "\$here/\$target"
        ${installRecord(id, lineage, receipts, [:])} --key "alphamissense_vep_vcf=\$target"
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
