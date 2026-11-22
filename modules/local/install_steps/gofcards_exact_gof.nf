include { installAdopt; installRecord; installStart; installStub; installWord } from './common'

// F5: the GoFCards exact gain-of-function variant cache
// (gofcards_exact_gof_cache). The repository ships it with its ClinVar
// conditions (data/gofcards/gofcards_exact_gof.json.gz), and the step records
// that file, which step F7 then records as it is. Rebuilding it
// (PRIVA_FORCE_GOFCARDS_CACHE, PRIVA_FORCE_GOFCARDS_CLINVAR_INJECTION or
// PRIVA_FORCE_ALL_CACHES) normalizes the GoFCards variants with VEP on the
// GRCh37 and GRCh38 caches of the configured release, which come with the
// process of the VEP cache step (A2); until then the step stops and names the
// installer function that still rebuilds it.
process GOFCARDS_EXACT_GOF {
    tag "${id} ${lineage}"
    label 'install'
    storeDir "${params.resource_root}/${id}/${lineage}"
    input:
    tuple val(id), val(lineage), path(plan, stageAs: 'plan.json'), path(receipts, stageAs: 'inputs/*/receipt.json'), val(planned), val(names), path(staged, stageAs: 'staged/*')
    output:
    path('receipt.json'), emit: receipt
    path('data'), emit: data
    script:
    def shipped = "${projectDir}/data/gofcards/gofcards_exact_gof.json.gz"
    if (lineage.startsWith('adopted-'))
        installAdopt(id, lineage)
    else
        installStart(task) + """
        shipped=${installWord(shipped)}
        if [[ \${PRIVA_FORCE_GOFCARDS_CACHE:-0} == 1 || \${PRIVA_FORCE_GOFCARDS_CLINVAR_INJECTION:-0} == 1 || \\
              \${PRIVA_FORCE_ALL_CACHES:-0} == 1 ]]; then
            echo 'GoFCards: a rebuild needs the GRCh37 and GRCh38 VEP caches, which this install does not prepare yet; rebuild it with scripts/install_utils.sh gofcards_exact_gof_cache_install CONFIG' >&2
            exit 1
        fi
        [[ -f \$shipped && ! -L \$shipped && -s \$shipped ]] || { echo "GoFCards: the shipped cache is missing or empty: \$shipped" >&2; exit 1; }
        echo "Reusing nonempty bundled gofcards_exact_gof_cache: \$shipped"
        ${installRecord(id, lineage, receipts, [:])} --key "gofcards_exact_gof_cache=\$shipped"
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
