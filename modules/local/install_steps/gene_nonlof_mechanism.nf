include { installAdopt; installRecord; installStart; installStub; installValue; installWord } from './common'

// F8: the curated non-LOF mechanism assertions of each gene
// (gene_nonlof_mechanism_json). The repository ships them
// (data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.json.gz),
// and the step records that file. With PRIVA_FORCE_NONLOF_MECHANISM_CACHE=1 or
// PRIVA_FORCE_ALL_CACHES=1 priva gene-nonlof-mechanism builds them instead,
// from G2P, ClinGen dosage, PanelApp and the GoFCards workbook (downloaded
// with curl), the GoFCards cache of step F7 and the ClinVar VCV XML of step F6
// (or its own download when F6 recorded none), and checks them against the
// shipped schema. Each input is read where its step keeps it (the target of
// its staged link), as the cache records the paths it read.
process GENE_NONLOF_MECHANISM {
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
    def shipped = "${projectDir}/data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.json.gz"
    def schema = "${projectDir}/data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.schema.json"
    def hgnc = "${projectDir}/data/hgnc/hgnc_non_alt_loci_set_sha256_890b18743cc8.tsv"
    def gofcardsLink = names.gofcards_exact_gof_cache ? "staged/${names.gofcards_exact_gof_cache}".toString() : ''
    def xmlLink = names.clinvar_vcv_xml ? "staged/${names.clinvar_vcv_xml}".toString() : ''
    def cache = 'data/patho_mechanism/gene_nonlof_mechanism_curated_assertions.json.gz'
    def root = installWord(projectDir)
    if (adopted)
        installAdopt(id, lineage)
    else
        installStart(task) + """
        here=\$(pwd -P)
        shipped=${installWord(shipped)}
        if [[ \${PRIVA_FORCE_NONLOF_MECHANISM_CACHE:-0} != 1 && \${PRIVA_FORCE_ALL_CACHES:-0} != 1 && \\
              -f \$shipped && ! -L \$shipped && -s \$shipped ]]; then
            echo "Reusing nonempty bundled gene_nonlof_mechanism_json: \$shipped"
            ${installRecord(id, lineage, receipts, [:])} --key "gene_nonlof_mechanism_json=\$shipped"
            exit 0
        fi
        gofcards=${installWord(gofcardsLink)}
        [[ -n \$gofcards && -s \$gofcards ]] || { echo 'ERROR: GoFCards exact GOF cache not found (step F7)' >&2; exit 1; }
        gofcards=\$(readlink -f "\$gofcards")
        options=()
        xml=${installWord(xmlLink)}
        [[ -z \$xml ]] || options+=(--clinvar-vcv "\$(readlink -f "\$xml")")
        [[ \${PRIVA_FORCE_NONLOF_MECHANISM_CACHE:-0} != 1 && \${PRIVA_FORCE_ALL_CACHES:-0} != 1 ]] || options+=(--force)
        [[ -z \${PROXY_URL:-} ]] || options+=(--proxy-url "\$PROXY_URL" --download-tool auto)
        mkdir -p data/patho_mechanism/raw
        echo 'Checking PriVA non-LOF/ClinVar mechanism sources and rebuilding only when inputs changed'
        "\$PRIVA_BIN" gene-nonlof-mechanism --source-root ${root} \\
            --cache-dir "\$here/data/patho_mechanism" --shared-raw-dir "\$here/data/patho_mechanism/raw" \\
            --gofcards-exact-variants "\$gofcards" --hgnc-table ${installWord(hgnc)} \\
            --output-json "\$here/${cache}" --output-schema ${installWord(schema)} \\
            --clinvar-min-review-stars "\${CLINVAR_GOF_MIN_REVIEW_STARS:-2}" \\
            --clinvar-max-download-seconds "\${CLINVAR_VCV_MAX_DOWNLOAD_SECONDS:-86400}" \\
            --timeout "\${GENE_MECHANISM_TIMEOUT:-120}" --retries "\${GENE_MECHANISM_RETRIES:-3}" \\
            "\${options[@]}" ${urls}
        "\$PRIVA_BIN" gene-nonlof-mechanism --source-root ${root} --output-schema ${installWord(schema)} --validate-only ${cache}
        relocate=()
        for written in data/patho_mechanism/*.json; do
            [[ -f \$written ]] && relocate+=(--relocate "\$written")
        done
        ${installRecord(id, lineage, receipts, [gene_nonlof_mechanism_json: cache])} "\${relocate[@]}"
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
