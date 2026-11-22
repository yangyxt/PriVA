include { installAdopt; installRecord; installSetting; installStart; installStub; installWord } from './common'

// C8: the InterPro entry of each protein domain (interpro_mapping_pickle). The
// repository ships it (data/InterPro/Interpro_entry_mapping.pkl.gz), and the
// step records that file. Without it, curl downloads InterPro's XML
// (interpro_xml_gz_url), gzip checks it, and protein_domain_mapping.py builds
// the mapping from it and writes a report of the domain names in the
// AlphaMissense VCF of step C4; the XML is kept with its SHA-256.
process INTERPRO_MAPPING {
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
    def url = adopted ? '' : installSetting(planned, 'interpro_xml_gz_url')
    def shipped = "${projectDir}/data/InterPro/Interpro_entry_mapping.pkl.gz"
    def xml = 'data/InterPro/' + url.substring(url.lastIndexOf('/') + 1)
    def vcf = names.alphamissense_vep_vcf ? "staged/${names.alphamissense_vep_vcf}".toString() : ''
    if (adopted)
        installAdopt(id, lineage)
    else
        installStart(task) + """
        shipped=${installWord(shipped)}
        if [[ -f \$shipped && ! -L \$shipped && -s \$shipped ]]; then
            echo "Reusing nonempty bundled InterPro mapping: \$shipped"
            mapping=\$shipped
        else
            vcf=${installWord(vcf)}
            [[ -n \$vcf && -s \$vcf ]] || { echo 'InterPro: the AlphaMissense VCF of step C4 is missing' >&2; exit 1; }
            mkdir -p data/InterPro
            curl -fL --retry 5 --retry-delay 10 -o ${xml}.partial.gz ${installWord(url)}
            [[ -s ${xml}.partial.gz ]] || { echo 'The InterPro XML resource is missing or empty' >&2; exit 1; }
            gzip -t ${xml}.partial.gz
            mv ${xml}.partial.gz ${xml}
            sha256sum ${xml} | awk '{print \$1}' > ${xml}.sha256
            mapping=data/InterPro/Interpro_entry_mapping.pkl.gz
            python ${installWord("${projectDir}/scripts/protein_domain_mapping.py")} --xml ${xml} \\
                --mapping_output "\$mapping" --vcf "\$vcf" --report data/InterPro/Interpro_format_report.txt
            [[ -s \$mapping ]] || { echo "ERROR: Missing or empty generated InterPro mapping: \$mapping" >&2; exit 1; }
        fi
        ${installRecord(id, lineage, receipts, [:])} --key "interpro_mapping_pickle=\$mapping"
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
