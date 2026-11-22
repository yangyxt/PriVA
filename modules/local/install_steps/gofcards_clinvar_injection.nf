include { installAdopt; installRecord; installStart; installStub; installValue; installWord } from './common'

// F7: the GoFCards cache with the ClinVar VCV conditions of each variant
// (gofcards_exact_gof_cache), the one GoFCards cache every later step reads.
// When step F5 records the shipped cache, which already carries them, the step
// records it as it is. When F5 built a normalized cache, priva
// gofcards-clinvar-injection attaches the conditions of the ClinVar VCV XML
// of step F6 (bound to its SHA-256) at clinvar_min_review_stars, writes the
// match audit and statistics beside it, and priva gofcards-cache-validate
// checks the cache and the XML it cites. Each input is read where its step
// keeps it (the target of its staged link), as the cache records the XML's
// path.
process GOFCARDS_CLINVAR_INJECTION {
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
    def stars = adopted ? '' : installValue(planned, 'clinvar_min_review_stars', '2')
    def shipped = "${projectDir}/data/gofcards/gofcards_exact_gof.json.gz"
    def sourceLink = names.gofcards_exact_gof_cache ? "staged/${names.gofcards_exact_gof_cache}".toString() : ''
    def xmlLink = names.clinvar_vcv_xml ? "staged/${names.clinvar_vcv_xml}".toString() : ''
    def work = 'data/gofcards/build_work'
    def target = 'data/gofcards/gofcards_exact_gof.json.gz'
    if (adopted)
        installAdopt(id, lineage)
    else
        installStart(task) + """
        source=${installWord(sourceLink)}
        [[ -n \$source && -s \$source ]] || { echo 'GoFCards: the cache of step F5 is missing' >&2; exit 1; }
        source=\$(readlink -f "\$source")
        shipped=\$(readlink -f ${installWord(shipped)})
        if [[ \$source == "\$shipped" ]]; then
            if [[ \${PRIVA_FORCE_GOFCARDS_CLINVAR_INJECTION:-0} == 1 || \${PRIVA_FORCE_ALL_CACHES:-0} == 1 ]]; then
                echo "ERROR: Force cannot overwrite bundled gofcards_exact_gof_cache: \$source; rebuild step F5 (PRIVA_FORCE_GOFCARDS_CACHE=1) for a cache to inject" >&2
                exit 1
            fi
            echo "Reusing nonempty bundled gofcards_exact_gof_cache: \$source"
            ${installRecord(id, lineage, receipts, [:])} --key "gofcards_exact_gof_cache=\$source"
            exit 0
        fi
        xml=${installWord(xmlLink)}
        [[ -n \$xml && -s \$xml ]] || { echo 'GoFCards: the ClinVar VCV XML of step F6 is missing' >&2; exit 1; }
        xml=\$(readlink -f "\$xml")
        identity=\$(stat -Lc '%d:%i:%s:%y:%z' "\$xml")
        sha=\$(sha256sum -- "\$xml" | awk '{print \$1}')
        [[ \$sha =~ ^[0-9a-f]{64}\$ ]] || { echo "ERROR: cannot hash \$xml" >&2; exit 1; }
        [[ \$(stat -Lc '%d:%i:%s:%y:%z' "\$xml") == "\$identity" ]] || { echo 'ERROR: Explicit ClinVar VCV XML changed while hashing; wait for its writer' >&2; exit 1; }
        mkdir -p ${work} data/patho_mechanism
        echo "Injecting ClinVar VCV conditions into GoFCards variants -> ${target}"
        set -o pipefail
        "\$PRIVA_BIN" gofcards-clinvar-injection --source-root ${installWord(projectDir)} \\
            --gofcards-cache "\$source" --clinvar-vcv "\$xml" --clinvar-vcv-sha256 "\$sha" \\
            --out-json ${target} --min-review-stars ${installWord(stars)} \\
            --audit-tsv data/patho_mechanism/clinvar_vcv_gofcards_matches.tsv \\
            --stats-json ${work}/clinvar_injection_stats.json 2>&1 | tee ${work}/clinvar_injection.log
        "\$PRIVA_BIN" gofcards-cache-validate ${target}
        [[ \$(stat -Lc '%d:%i:%s:%y:%z' "\$xml") == "\$identity" ]] || { echo 'ERROR: Explicit ClinVar VCV XML changed during cache preparation' >&2; exit 1; }
        "\$PRIVA_BIN" gofcards-cache-validate ${target} "clinvar.source=\$xml" "clinvar.source_sha256=\$sha"
        ${installRecord(id, lineage, receipts, [gofcards_exact_gof_cache: target])}
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
