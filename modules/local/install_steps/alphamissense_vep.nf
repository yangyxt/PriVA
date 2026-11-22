include { installAdopt; installRecord; installSetting; installStart; installStub; installWord } from './common'

// C3: the AlphaMissense table of step B3 as a VCF (alphamissense_vcf) and that
// VCF annotated by VEP (alphamissense_vep_vcf, also named
// alphamissense_vep_base_vcf: the complete annotation before step C4 adds
// gnomAD; C4 writes its own copy). priva alphamissense-vcf build converts the
// table; vep annotates it with the cache and plugins of step A2, the reference
// of A1 and the Conservation plugin on the track of B2, as
// basic_vep_annotation ran it; bcftools sorts and tabix indexes the result.
// The annotation must name the VEP release and assembly
// (priva alphamissense-cache check-annotation), carry CSQ and hold the
// converted VCF's records, chromosome by chromosome (bcftools index --stats).
process ALPHAMISSENSE_VEP {
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
    def assembly = configured in ['hg19', 'GRCh37'] ? 'GRCh37' : (configured in ['hg38', 'GRCh38'] ? 'GRCh38' : '')
    def short_name = assembly == 'GRCh37' ? 'hg19' : 'hg38'
    def release = adopted ? '' : installSetting(planned, 'vep_installed_version')
    def links = ['ref_genome', 'vep_cache_dir', 'vep_plugins_dir', 'conservation_file', 'alphamissense_prescore'].collectEntries { key ->
        [(key): names[key] ? "staged/${names[key]}".toString() : '']
    }
    def folder = 'data/VEP_plugin_cache/AlphaMissense'
    def converted = "${folder}/AlphaMissense_${short_name}.vcf.gz".toString()
    def annotated = "${folder}/AlphaMissense_${short_name}.vep.vcf.gz".toString()
    def stats = "${folder}/AlphaMissense_${short_name}.vep.stats.html".toString()
    def staged_inputs = [['ref_genome', 'A1'], ['vep_cache_dir', 'A2'], ['vep_plugins_dir', 'A2'], ['alphamissense_prescore', 'B3']]
    def missing = staged_inputs.findAll { pair -> !names[pair[0]] }.collect { pair -> "${pair[0]} (step ${pair[1]})".toString() }
    if (adopted)
        installAdopt(id, lineage)
    else if (!assembly)
        "echo 'AlphaMissense annotation: select hg19/GRCh37 or hg38/GRCh38, not ${configured}' >&2; exit 2"
    else if (!(release ==~ /[1-9][0-9]*/))
        "echo 'AlphaMissense annotation: vep_installed_version must be a numeric major release, not ${release}' >&2; exit 2"
    else if (!missing.isEmpty())
        "echo ${installWord('AlphaMissense annotation: missing ' + missing.join(', '))} >&2; exit 1"
    else
        installStart(task) + """
        set -o pipefail
        here=\$(pwd -P)
        # The files of steps A1, A2, B2 and B3 where they are (VEP finds the
        # reference's index beside it, and names the cache it read).
        reference=\$(readlink -f ${installWord(links.ref_genome)})
        cache=\$(readlink -f ${installWord(links.vep_cache_dir)})
        plugins=\$(readlink -f ${installWord(links.vep_plugins_dir)})
        prescore=\$(readlink -f ${installWord(links.alphamissense_prescore)})
        conservation=()
        if [[ -n ${installWord(links.conservation_file)} ]]; then
            conservation=(-plugin "Conservation,\$(readlink -f ${installWord(links.conservation_file)}),MAX")
        fi
        for tool in vep bcftools tabix; do
            type -P \$tool > /dev/null || { echo "AlphaMissense annotation: \$tool is not on PATH" >&2; exit 1; }
        done
        actual=\$(vep --help 2>&1 | awk '/ensembl-vep[[:space:]]/ {print \$NF; exit}')
        [[ \${actual%%.*} == ${release} ]] || {
            echo "AlphaMissense annotation: vep_installed_version is ${release}, the vep on PATH is \$actual" >&2; exit 1; }
        has_record() { bcftools view -h "\$1" > /dev/null && [[ -n \$(bcftools view -H "\$1" 2> /dev/null | head -n 1) ]]; }
        # Records per chromosome, from the index (zero counts left out).
        counts() { bcftools index --stats "\$1" | awk -F '\\t' '\$3 != 0 {print \$1 "\\t" \$3}'; }
        mkdir -p ${folder} work
        "\$PRIVA_BIN" alphamissense-vcf build --input "\$prescore" --output ${converted}
        has_record ${converted} || { echo 'AlphaMissense annotation: the converted VCF holds no record' >&2; exit 1; }
        # basic_vep_annotation, without SpliceAI.
        vep -i ${converted} --format vcf --vcf --species homo_sapiens --assembly ${assembly} \\
            --cache --offline --merged --use_transcript_ref --hgvs --symbol --canonical --total_length \\
            --domains --numbers --stats_file ${stats} --fork ${task.cpus} --buffer_size 10000 \\
            --fasta "\$reference" --dir_cache "\$cache" --dir_plugins "\$plugins" "\${conservation[@]}" \\
            --force_overwrite -o work/AlphaMissense.vep.vcf
        bcftools sort -T "\$here/work/sort.XXXXXX" -Oz -o ${annotated} work/AlphaMissense.vep.vcf
        tabix -f -p vcf ${annotated}
        rm -rf work
        has_record ${annotated} || { echo 'AlphaMissense annotation: the VEP output holds no record' >&2; exit 1; }
        bcftools view -h ${annotated} | grep -q '^##INFO=<ID=CSQ,' || {
            echo 'AlphaMissense annotation: the VEP output has no CSQ field' >&2; exit 1; }
        "\$PRIVA_BIN" alphamissense-cache check-annotation --input ${annotated} \\
            --assembly ${installWord(configured)} --vep-major ${release} > /dev/null
        [[ \$(counts ${converted}) == "\$(counts ${annotated})" ]] || {
            echo 'AlphaMissense annotation: VEP lost or added records' >&2; exit 1; }
        ${installRecord(id, lineage, receipts, [alphamissense_vcf: converted, alphamissense_vep_vcf: annotated, alphamissense_vep_base_vcf: annotated])}
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
