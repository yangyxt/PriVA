include { installAdopt; installRecord; installStart; installStub; installValue; installWord } from './common'

// F9: the HPO-framed condition and pathogenic-mechanism cache
// (hpo_condition_mechanism_json). The repository ships it
// (data/patho_mechanism/hpo_condition_mechanism_cache.json.gz), and the step
// records that file. With PRIVA_FORCE_HPO_MECHANISM_JSON=1 or
// PRIVA_FORCE_ALL_CACHES=1 it builds it instead: priva hpo-condition-mechanism
// joins the scoped HPO assertions of step F1, the mechanism evidence of step
// F4, the GoFCards cache of step F7 and the ClinGen dosage table of step D9,
// and checks the cache against them. Each input is read where its step keeps
// it (the target of its staged link), as the cache records the paths it read.
process HPO_CONDITION_MECHANISM {
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
    def hpoRelease = adopted ? '' : installValue(planned, 'hpo_release', '')
    def mondoRelease = adopted ? '' : installValue(planned, 'mondo_release', '')
    def shipped = "${projectDir}/data/patho_mechanism/hpo_condition_mechanism_cache.json.gz"
    def assertionsLink = names.hpo_assertions ? "staged/${names.hpo_assertions}".toString() : ''
    def evidenceLink = names.ddg2p_mechanism_evidence ? "staged/${names.ddg2p_mechanism_evidence}".toString() : ''
    def dosageLink = names.gene_dosage_sensitivity ? "staged/${names.gene_dosage_sensitivity}".toString() : ''
    def gofcardsLink = names.gofcards_exact_gof_cache ? "staged/${names.gofcards_exact_gof_cache}".toString() : ''
    def cache = 'data/patho_mechanism/hpo_condition_mechanism_cache.json.gz'
    def root = installWord(projectDir)
    if (adopted)
        installAdopt(id, lineage)
    else
        installStart(task) + """
        shipped=${installWord(shipped)}
        if [[ \${PRIVA_FORCE_HPO_MECHANISM_JSON:-0} != 1 && \${PRIVA_FORCE_ALL_CACHES:-0} != 1 && \\
              -f \$shipped && ! -L \$shipped && -s \$shipped ]]; then
            echo "Reusing nonempty bundled hpo_condition_mechanism_json: \$shipped"
            cache=\$shipped
        else
            inputs=()
            for pair in --hpo-assertions=${installWord(assertionsLink)} --mechanism-evidence=${installWord(evidenceLink)} \\
                        --clingen-dosage=${installWord(dosageLink)} --gofcards-variants=${installWord(gofcardsLink)}; do
                option=\${pair%%=*}
                input=\${pair#*=}
                [[ -n \$input && -s \$input ]] || { echo "HPO condition mechanism: the input of \$option is missing (steps F1, F4, D9, F7)" >&2; exit 1; }
                inputs+=("\$option" "\$(readlink -f "\$input")")
            done
            inputs+=(--hpo-release ${installWord(hpoRelease)} --mondo-release ${installWord(mondoRelease)})
            cache=${cache}
            mkdir -p data/patho_mechanism
            echo "Building atomic HPO-framed condition mechanism JSON: \$cache"
            "\$PRIVA_BIN" hpo-condition-mechanism --source-root ${root} "\${inputs[@]}" --output "\$cache"
            "\$PRIVA_BIN" hpo-condition-mechanism --source-root ${root} --validate-only "\$cache" "\${inputs[@]}"
        fi
        ${installRecord(id, lineage, receipts, [:])} --key "hpo_condition_mechanism_json=\$cache"
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
