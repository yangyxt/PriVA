include { installAdopt; installRecord; installSetting; installStart; installStub; installValue; installWord } from './common'

// D1_chr1 to D1_chrOther: VEP on one piece of step D1_split's ClinVar release
// (clinvar_vep_chrN: clinvar.chrN.vep.vcf.gz, bgzipped; step D1 stages all 25
// to gather them, so each has its own name). The 25 steps share this
// process (registry "process"); the piece is the step's id after D1_. priva
// task step --step vep-worker runs the pipeline's VEP command line
// (anno_VEP_cache_worker of annotation_vcf.sh: the plugins and files of steps
// A1, A2 and B1 to B9) with 22 forks and leaves VEP's output unsorted. VEP on
// the whole release took 34 hours (job 1540983); with 22 forks per chromosome
// it gave the same rows (jobs 1524533 and 1524532, 26 October 2026), and with
// 8 forks the order of domains differed, so the forks stay 22 (nextflow.config
// gives these steps 24 CPUs: VEP starts up to forks + 1 children besides
// itself). A piece without records has nothing to annotate; its header-only
// file stands as the output.
process CLINVAR_VEP_CHROM {
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
    def piece = id.replaceFirst(/^D1_/, '')
    def known = (1..22).collect { n -> "chr${n}".toString() } + ['chrX', 'chrY', 'chrOther']
    def assembly = adopted ? '' : installSetting(planned, 'assembly')
    def loeuf = assembly in ['hg19', 'GRCh37'] ? 'loeuf_dataset.tsv.gz' : 'loeuf_dataset_hg38.tsv.gz'
    def resources = ['ref_genome', 'vep_cache_dir', 'vep_plugins_dir', 'utr_annotator_file', 'conservation_file',
                     'alphamissense_prescore', 'spliceai_snv_prescore', 'spliceai_indel_prescore', 'splicevault_prescore',
                     'mavedb_file', 'primateai_prescore', 'loftee_repo', 'human_ancestor_fasta', 'loftee_conservation_file',
                     'gerp_bigwig', 'pext_mean_bw', 'pext_tissue_bw_dir']
    // The configuration VEP reads: each staged file where its step keeps it
    // (empty when no step gave it; VEP's checks name what it lacks).
    def lines = resources.collect { key ->
        names[key] ? "line ${key} \"\$(readlink -f ${installWord('staged/' + names[key])})\"" : "line ${key} ''"
    }.join('\n        ')
    def folder = "data/ClinVar/vep_${piece}".toString()
    def annotated = "${folder}/clinvar.${piece}.vep.vcf".toString()
    def key = "clinvar_vep_${piece}".toString()
    if (adopted)
        installAdopt(id, lineage)
    else if (!(piece in known))
        "echo 'ClinVar VEP: step ${id} names no piece of the ClinVar release' >&2; exit 1"
    else if (!names.clinvar_vep_split)
        "echo 'ClinVar VEP: clinvar_vep_split (step D1_split) is missing' >&2; exit 1"
    else
        installStart(task) + """
        here=\$(pwd -P)
        manifest=\$(readlink -f ${installWord('staged/' + names.clinvar_vep_split)})
        piece=\$(jq -r --arg n ${piece} '.pieces[] | select(.name == \$n) | .path' "\$manifest")
        records=\$(jq -r --arg n ${piece} '.pieces[] | select(.name == \$n) | .records' "\$manifest")
        [[ -s \$piece && \$records =~ ^[0-9]+\$ ]] || { echo "ClinVar VEP: no piece ${piece} in \$manifest" >&2; exit 1; }
        mkdir -p ${folder} plugin_cache
        line() { local q="'"; printf "%s: '%s'\\n" "\$1" "\${2//\$q/\$q\$q}"; }
        {
        line assembly ${installWord(assembly)}
        line pext_enabled ${installWord(installValue(planned, 'pext_enabled', 'null'))}
        line pext_tissues ${installWord(installValue(planned, 'pext_tissues', 'null'))}
        line loeuf_prescore ${installWord("${projectDir}/data/loeuf/${loeuf}")}
        line vep_plugins_cachedir "\$here/plugin_cache"
        ${lines}
        } > vep.yaml
        if ((records == 0)); then
            cp -- "\$piece" ${annotated}
        else
            "\$PRIVA_BIN" task step --step vep-worker --root ${installWord(projectDir)} --config "\$here/vep.yaml" \\
                --vcf "\$piece" --output "\$here/${annotated}" --stats "\$here/${folder}/clinvar.${piece}.vep.stats.html" --forks 22
        fi
        bgzip -f -@ ${task.cpus} ${annotated}
        ${installRecord(id, lineage, receipts, [(key): "${annotated}.gz".toString()])}
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
