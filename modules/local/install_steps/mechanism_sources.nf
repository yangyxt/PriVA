include { installAdopt; installRecord; installStart; installStub; installValue; installWord } from './common'

// F6: the sources of the non-LOF mechanism cache, above all the weekly ClinVar
// VCV XML (clinvar_vcv_xml) that the ClinVar injection of step F7 reads before
// step F8 builds the cache. With the shipped non-LOF cache and the shipped
// GoFCards cache in use, nothing reads them, and the step records no file.
// Otherwise (PRIVA_FORCE_NONLOF_MECHANISM_CACHE, PRIVA_FORCE_GOFCARDS_CACHE,
// PRIVA_FORCE_GOFCARDS_CLINVAR_INJECTION or PRIVA_FORCE_ALL_CACHES) priva
// gene-nonlof-mechanism --fetch-sources-only downloads them with curl (the
// XML checked against its .md5) and records what it fetched, and the step
// records the XML, or the clinvar_vcv_xml you supply, which priva checks.
process MECHANISM_SOURCES {
    tag "${id} ${lineage}"
    label 'install_download'
    storeDir "${params.resource_root}/${id}/${lineage}"
    input:
    tuple val(id), val(lineage), path(plan, stageAs: 'plan.json'), path(receipts, stageAs: 'inputs/*/receipt.json'), val(planned), val(names), path(staged, stageAs: 'staged/*')
    output:
    path('receipt.json'), emit: receipt
    path('data'), emit: data
    script:
    def adopted = lineage.startsWith('adopted-')
    def urls = adopted ? [] : ['g2p_ddg2p_url', 'clingen_dosage_grch38_url', 'panelapp_panels_url'].collect { key ->
        def value = installValue(planned, key, '')
        value ? "--${key.replace('_', '-')} ${installWord(value)}" : ''
    }.findAll { it }.join(' ')
    def nonlof = "${projectDir}/data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.json.gz"
    def schema = "${projectDir}/data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.schema.json"
    def gofcards = "${projectDir}/data/gofcards/gofcards_exact_gof.json.gz"
    def hgnc = "${projectDir}/data/hgnc/hgnc_non_alt_loci_set_sha256_890b18743cc8.tsv"
    def suppliedLink = names.clinvar_vcv_xml ? "staged/${names.clinvar_vcv_xml}".toString() : ''
    def managed = 'data/patho_mechanism/raw/clinvar_vcv/ClinVarVCVRelease_00-latest_weekly.xml.gz'
    def root = installWord(projectDir)
    if (adopted)
        installAdopt(id, lineage)
    else
        installStart(task) + """
        here=\$(pwd -P)
        usable() { [[ -f \$1 && ! -L \$1 && -s \$1 ]]; }
        nonlof=${installWord(nonlof)}
        gofcards=${installWord(gofcards)}
        force_nonlof=0
        [[ \${PRIVA_FORCE_NONLOF_MECHANISM_CACHE:-0} == 1 || \${PRIVA_FORCE_ALL_CACHES:-0} == 1 ]] && force_nonlof=1
        force_gofcards=0
        [[ \${PRIVA_FORCE_GOFCARDS_CACHE:-0} == 1 || \${PRIVA_FORCE_GOFCARDS_CLINVAR_INJECTION:-0} == 1 || \\
           \${PRIVA_FORCE_ALL_CACHES:-0} == 1 ]] && force_gofcards=1
        if [[ \$force_nonlof == 0 && \$force_gofcards == 0 ]] && usable "\$nonlof" && usable "\$gofcards"; then
            echo 'The shipped non-LOF and GoFCards caches are in use; no mechanism source to fetch'
            ${installRecord(id, lineage, receipts, [:])}
            exit 0
        fi
        # What the installer passed: the shipped cache paths (or this step's
        # own when forced), which a source fetch does not read or write.
        [[ \$force_nonlof == 0 ]] && usable "\$nonlof" || nonlof=\$here/data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.json.gz
        [[ \$force_gofcards == 0 ]] && usable "\$gofcards" || gofcards=\$here/data/gofcards/gofcards_exact_gof.json.gz
        options=()
        supplied=${installWord(suppliedLink)}
        if [[ -n \$supplied ]]; then
            supplied=\$(readlink -f "\$supplied")
            [[ -s \$supplied ]] || { echo "ERROR: Missing explicit ClinVar VCV XML: \$supplied; restore that selected input" >&2; exit 1; }
            options+=(--clinvar-vcv "\$supplied")
        fi
        [[ \$force_nonlof == 0 ]] || options+=(--force)
        [[ -z \${PROXY_URL:-} ]] || options+=(--proxy-url "\$PROXY_URL" --download-tool auto)
        mkdir -p data/patho_mechanism/raw
        echo 'Checking PriVA non-LOF/ClinVar mechanism sources only; no cache build in this phase'
        "\$PRIVA_BIN" gene-nonlof-mechanism --source-root ${root} \\
            --cache-dir "\$here/data/patho_mechanism" --shared-raw-dir "\$here/data/patho_mechanism/raw" \\
            --gofcards-exact-variants "\$gofcards" --hgnc-table ${installWord(hgnc)} \\
            --output-json "\$nonlof" --output-schema ${installWord(schema)} \\
            --clinvar-min-review-stars "\${CLINVAR_GOF_MIN_REVIEW_STARS:-2}" \\
            --clinvar-max-download-seconds "\${CLINVAR_VCV_MAX_DOWNLOAD_SECONDS:-86400}" \\
            --timeout "\${GENE_MECHANISM_TIMEOUT:-120}" --retries "\${GENE_MECHANISM_RETRIES:-3}" \\
            "\${options[@]}" ${urls} --fetch-sources-only
        xml=${managed}
        [[ -z \$supplied ]] || xml=\$supplied
        [[ -s \$xml ]] || { echo "ERROR: ClinVar VCV XML missing after the source fetch: \$xml" >&2; exit 1; }
        echo "PriVA non-LOF mechanism sources staged; ClinVar VCV XML ready: \$xml"
        relocate=()
        for written in data/patho_mechanism/*.json; do
            [[ -f \$written ]] && relocate+=(--relocate "\$written")
        done
        ${installRecord(id, lineage, receipts, [:])} --key "clinvar_vcv_xml=\$xml" "\${relocate[@]}"
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
