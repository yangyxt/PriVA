include { installAdopt; installRecord; installSetting; installStart; installStub; installWord } from './common'

// D1_split: D0's ClinVar release (clinvar_release_inputs names its VCF) with
// UCSC contig names, written as one plain VCF per chromosome for the VEP steps
// D1_chr1 to D1_chrOther (chr1..chr22, chrX, chrY; every other contig in
// chrOther), each with the whole header. bcftools renames the contigs with the
// shipped GRC_to_ucsc map and sorts (as liftover_from_ucsc_to_GRCh did), awk
// splits the records, and split.json (clinvar_vep_split) lists the pieces,
// their records and a hash of all records, which step D1 checks its gathered
// VEP output against.
process CLINVAR_VEP_SPLIT {
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
    def configured = adopted ? '' : installSetting(planned, 'assembly')
    def ucsc = configured in ['hg19', 'GRCh37'] ? 'hg19' : (configured in ['hg38', 'GRCh38'] ? 'hg38' : '')
    def release = names.clinvar_release_inputs ? "staged/${names.clinvar_release_inputs}".toString() : ''
    def pieces = ((1..22).collect { n -> "chr${n}".toString() } + ['chrX', 'chrY', 'chrOther']).join(' ')
    def folder = 'data/ClinVar/vep_split'
    def converted = "${folder}/clinvar.${ucsc}.vcf.gz".toString()
    if (adopted)
        installAdopt(id, lineage)
    else if (!ucsc)
        "echo 'ClinVar split: unsupported assembly ${configured}' >&2; exit 1"
    else if (!release)
        "echo 'ClinVar split: clinvar_release_inputs (step D0) is missing' >&2; exit 1"
    else
        installStart(task) + """
        set -o pipefail
        here=\$(pwd -P)
        release=\$(readlink -f ${installWord(release)})
        source=\$(jq -r '.sources.vcf.path // empty' "\$release")
        [[ -s \$source && -s \$source.tbi ]] || { echo "ClinVar split: D0's ClinVar VCF or its index is missing: \$source" >&2; exit 1; }
        contigs=${installWord("${projectDir}/data/liftover/GRC_to_ucsc.contig.map.txt")}
        [[ -f \$contigs ]] || { echo "ClinVar split: the contig map \$contigs is missing" >&2; exit 1; }
        mkdir -p ${folder}
        bcftools annotate --rename-chrs "\$contigs" "\$source" -Ou | bcftools sort -T "\$here/sort.XXXXXX" -Oz -o ${converted}
        tabix -f -p vcf ${converted}
        bcftools view -h ${converted} > /dev/null
        [[ -n \$(bcftools view -H ${converted} 2> /dev/null | head -n 1) ]] || { echo 'ClinVar split: the renamed release holds no record' >&2; exit 1; }
        bgzip -dc ${converted} | awk -v dir="\$here/${folder}" -v names='${pieces}' '
            BEGIN { n = split(names, list, " "); for (i = 1; i <= n; i++) known[list[i]] = 1 }
            /^#/ { header = header \$0 "\\n"; next }
            {
                name = (\$1 in known) ? \$1 : "chrOther"
                if (!(name in count)) { count[name] = 0; printf "%s", header > (dir "/" name ".vcf") }
                print > (dir "/" name ".vcf"); count[name]++
            }
            END {
                for (i = 1; i <= n; i++) {
                    if (!(list[i] in count)) { count[list[i]] = 0; printf "%s", header > (dir "/" list[i] ".vcf") }
                    print list[i] "\\t" count[list[i]] > (dir "/counts.tsv")
                }
            }'
        body=\$(bgzip -dc ${converted} | grep -v '^#' | sha256sum | cut -c1-64)
        pieces='[]'
        for name in ${pieces}; do
            records=\$(awk -v n="\$name" '\$1 == n { print \$2 }' ${folder}/counts.tsv)
            [[ \$records =~ ^[0-9]+\$ && -s ${folder}/\$name.vcf ]] || { echo "ClinVar split: piece \$name was not written" >&2; exit 1; }
            pieces=\$(jq -c --arg n "\$name" --arg p "\$here/${folder}/\$name.vcf" --argjson r "\$records" \\
                '. + [{name: \$n, path: \$p, records: \$r}]' <<< "\$pieces")
        done
        jq -n --arg source "\$source" --arg converted "\$here/${converted}" --arg body "\$body" --argjson pieces "\$pieces" \\
            '{kind: "priva-clinvar-vep-split-v1", source_vcf: \$source, converted_vcf: \$converted,
              body_sha256: \$body, records: ([\$pieces[].records] | add), pieces: \$pieces}' > ${folder}/split.json
        ${installRecord(id, lineage, receipts, [clinvar_vep_split: "${folder}/split.json".toString()])} --relocate ${folder}/split.json
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
