include { installAdopt; installRecord; installStart; installStub; installValue; installWord } from './common'

// F4: the condition mechanism evidence of each gene (ddg2p_mechanism_evidence),
// from G2P and Orphadata. The repository ships it
// (data/patho_mechanism/gene_pathogenic_mechanism_evidence.tsv), and the step
// records that file. With PRIVA_FORCE_GENE_MECHANISM_CACHE=1 or
// PRIVA_FORCE_ALL_CACHES=1 it builds it instead: curl downloads G2P's table
// (g2p_ddg2p_url) and Orphadata's product 6 (orphadata_gene_disease_url),
// priva mechanism-source checks a selected (not default) source and records
// each download, priva gene-mechanism prepared joins them to the MONDO disease
// scope of step F1, and priva mechanism-evidence validate checks the table.
// curl takes the timeout, retries and proxy of the installer
// (GENE_MECHANISM_TIMEOUT, GENE_MECHANISM_RETRIES, PROXY_URL); a failed
// download fails the step.
process GENE_PATHOGENIC_MECHANISM {
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
    def g2pDefault = 'https://www.ebi.ac.uk/gene2phenotype/api/panel/all/download'
    def orphaDefault = 'https://www.orphadata.com/data/xml/en_product6.xml'
    def g2pUrl = adopted ? '' : installValue(planned, 'g2p_ddg2p_url', g2pDefault)
    def orphaUrl = adopted ? '' : installValue(planned, 'orphadata_gene_disease_url', orphaDefault)
    def shipped = "${projectDir}/data/patho_mechanism/gene_pathogenic_mechanism_evidence.tsv"
    def registryLink = names.mondo_disease_scope_registry ? "staged/${names.mondo_disease_scope_registry}".toString() : ''
    def evidence = 'data/patho_mechanism/gene_pathogenic_mechanism_evidence.tsv'
    def raw = 'data/patho_mechanism/raw'
    if (adopted)
        installAdopt(id, lineage)
    else
        installStart(task) + """
        shipped=${installWord(shipped)}
        if [[ \${PRIVA_FORCE_GENE_MECHANISM_CACHE:-0} != 1 && \${PRIVA_FORCE_ALL_CACHES:-0} != 1 && \\
              -f \$shipped && ! -L \$shipped && -s \$shipped ]]; then
            echo "Reusing nonempty bundled ddg2p_mechanism_evidence: \$shipped"
            evidence=\$shipped
        else
            registry=${installWord(registryLink)}
            [[ -n \$registry && -s \$registry ]] || { echo 'Gene mechanism: the MONDO disease scope of step F1 is missing' >&2; exit 1; }
            registry=\$(readlink -f "\$registry")
            timeout=\${GENE_MECHANISM_TIMEOUT:-120}
            retries=\${GENE_MECHANISM_RETRIES:-3}
            transfer=(-L --fail --show-error --silent --connect-timeout "\$timeout"
                --max-time \$(( timeout * 4 > timeout + 60 ? timeout * 4 : timeout + 60 ))
                --retry \$(( retries > 1 ? retries - 1 : 0 ))
                --user-agent 'PriVA/0.1 (local pathogenic-mechanism cache; contact: yangyxt)')
            [[ -z \${PROXY_URL:-} ]] || transfer+=(--proxy "\$PROXY_URL")
            # curl; a selected source is checked (KIND) before it is kept;
            # priva writes its .source.json.
            mechanism_download() {
                local url=\$1 default=\$2 output=\$3 kind=\$4
                mkdir -p "\$(dirname "\$output")"
                curl "\${transfer[@]}" --output "\$output.download.tmp" "\$url"
                [[ \$url == "\$default" ]] || "\$PRIVA_BIN" mechanism-source validate --input "\$output.download.tmp" --kind "\$kind"
                mv -T "\$output.download.tmp" "\$output"
                "\$PRIVA_BIN" mechanism-source record-source --input "\$output" --url "\$url"
            }
            mechanism_download ${installWord(g2pUrl)} ${installWord(g2pDefault)} ${raw}/AllG2P.csv g2p_ddg2p
            mechanism_download ${installWord(orphaUrl)} ${installWord(orphaDefault)} ${raw}/orphadata/en_product6.xml orphadata
            evidence=${evidence}
            "\$PRIVA_BIN" gene-mechanism prepared --raw-dir ${raw} --output "\$evidence" \\
                --disease-scope-registry "\$registry" \\
                --g2p-ddg2p-url ${installWord(g2pUrl)} --orphadata-gene-disease-url ${installWord(orphaUrl)}
        fi
        "\$PRIVA_BIN" mechanism-evidence validate --evidence "\$evidence"
        ${installRecord(id, lineage, receipts, [:])} --key "ddg2p_mechanism_evidence=\$evidence"
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
