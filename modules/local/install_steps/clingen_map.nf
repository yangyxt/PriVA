include { installAdopt; installRecord; installSetting; installStart; installStub; installValue; installWord } from './common'

// D9: the ClinGen gene-disease map (clingen_map) and the ClinGen dosage
// sensitivity table (gene_dosage_sensitivity). The repository ships both for
// hg19 and hg38 (data/clingen/), and the step records those files. With
// PRIVA_FORCE_ALL_CACHES=1 it builds them instead, as the installer did:
// priva clingen-map build reads the ClinGen eRepo table (the shipped one, or a
// download) and the VEP-annotated ClinVar VCF of step D1, and curl downloads
// the dosage table, which priva mechanism-source checks and records before awk
// keeps its header and records.
process CLINGEN_MAP {
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
    def bundle = assembly in ['hg19', 'GRCh37'] ? 'hg19' : (assembly in ['hg38', 'GRCh38'] ? 'hg38' : '')
    def grch = bundle == 'hg19' ? 'GRCh37' : 'GRCh38'
    def defaultDosage = "http://ftp.clinicalgenome.org/ClinGen_gene_curation_list_${grch}.tsv"
    def dosageUrl = bundle == 'hg38' ? installValue(planned, 'clingen_dosage_grch38_url', defaultDosage) : defaultDosage
    def erepoUrl = 'https://erepo.clinicalgenome.org/evrepo/api/summary/classifications/download'
    def shipped = "${projectDir}/data/clingen"
    def vcf = names.clinvar_vcf ? "staged/${names.clinvar_vcf}" : ''
    if (adopted)
        installAdopt(id, lineage)
    else if (!bundle)
        "echo 'ClinGen: unsupported assembly ${assembly}' >&2; exit 1"
    else
        installStart(task) + """
        # curl, then priva checks the table (KIND) and writes its .source.json.
        clingen_download() {
            local url=\$1 output=\$2 kind=\$3
            mkdir -p "\$(dirname "\$output")"
            curl -fL --retry 5 --retry-delay 10 -o "\$output.download" "\$url"
            [[ -s \$output.download ]] || { echo "The downloaded ClinGen resource is missing or empty: \$url" >&2; return 1; }
            "\$PRIVA_BIN" mechanism-source validate --input "\$output.download" --kind "\$kind"
            "\$PRIVA_BIN" mechanism-source record-source --input "\$output.download" --url "\$url"
            mv -T "\$output.download.source.json" "\$output.source.json"
            mv -T "\$output.download" "\$output"
        }
        map=''
        dosage=''
        if [[ \${PRIVA_FORCE_ALL_CACHES:-0} != 1 ]]; then
            [[ -s ${installWord(shipped + "/clingen_map.${bundle}.pkl.gz")} ]] && map=${installWord(shipped + "/clingen_map.${bundle}.pkl.gz")}
            [[ -s ${installWord(shipped + "/gene_dosage_sensitivity.${bundle}.tsv")} ]] && dosage=${installWord(shipped + "/gene_dosage_sensitivity.${bundle}.tsv")}
        fi
        if [[ -n \$map ]]; then
            echo "Reusing the shipped ClinGen map: \$map"
        else
            erepo=${installWord(shipped + '/clingen_erepo.tsv')}
            if [[ ! -s \$erepo ]]; then
                erepo=data/ClinGen/clingen_erepo.tsv
                clingen_download ${installWord(erepoUrl)} "\$erepo" clingen_erepo
            fi
            vcf=${installWord(vcf)}
            [[ -n \$vcf && -s \$vcf ]] || { echo 'ClinGen: the ClinVar VCF of step D1 is missing' >&2; exit 1; }
            map=data/ClinGen/clingen_map.${bundle}.rkyv
            mkdir -p data/ClinGen
            "\$PRIVA_BIN" clingen-map build --clinvar-vcf "\$vcf" --clingen-tab "\$erepo" --output "\$map"
        fi
        if [[ -n \$dosage ]]; then
            echo "Reusing the shipped ClinGen dosage sensitivity table: \$dosage"
        else
            dosage=data/ClinGen/gene_dosage_sensitivity.${bundle}.tsv
            clingen_download ${installWord(dosageUrl)} "\$dosage" clingen_dosage_runtime
            awk -F '\\t' '\$1 !~ /^#/ {print;} \$1 ~ /^#Gene/{print;}' "\$dosage" > "\$dosage.clean"
            mv "\$dosage.clean" "\$dosage"
            "\$PRIVA_BIN" mechanism-source record-source --input "\$dosage" --url ${installWord(dosageUrl)}
        fi
        ${installRecord(id, lineage, receipts, [:])} --key "clingen_map=\$map" --key "gene_dosage_sensitivity=\$dosage"
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
