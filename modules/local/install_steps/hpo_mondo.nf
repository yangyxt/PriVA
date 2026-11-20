include { installAdopt; installRecord; installStart; installStub; installValue; installWord } from './common'

// F1: the HPO gene-phenotype assertions scoped by MONDO disease
// (hpo_assertions), the MONDO ontology (mondo_obo) and the disease-scope
// registry (mondo_disease_scope_registry). curl downloads HPO's
// genes_to_phenotype.txt and phenotype.hpoa and MONDO's mondo-simple.obo of the
// pinned releases (and checks mondo_sha256 when it is set); priva
// hpo-assertions build collapses the HPO tables, priva mondo-scope build scopes
// them by MONDO with the overrides shipped in data/mondo, and priva mondo-scope
// validate checks the pair. Unset settings take the releases the installer
// pinned.
process HPO_MONDO {
    tag "${id} ${lineage}"
    label 'install'
    storeDir "${params.resource_root}/${id}/${lineage}"
    input:
    tuple val(id), val(lineage), path(plan, stageAs: 'plan.json'), path(receipts, stageAs: 'inputs/*/receipt.json'), val(planned), val(names), path(staged, stageAs: 'staged/*')
    output:
    path('receipt.json'), emit: receipt
    path('data'), emit: data
    script:
    def hpo = installValue(planned, 'hpo_release', 'v2026-06-23')
    def hpoBase = "https://github.com/obophenotype/human-phenotype-ontology/releases/download/${hpo}"
    def genesUrl = installValue(planned, 'hpo_genes_to_phenotype_url', "${hpoBase}/genes_to_phenotype.txt")
    def hpoaUrl = installValue(planned, 'hpo_phenotype_hpoa_url', "${hpoBase}/phenotype.hpoa")
    def mondo = installValue(planned, 'mondo_release', 'v2026-07-06')
    def mondoUrl = installValue(planned, 'mondo_url', "https://github.com/monarch-initiative/mondo/releases/download/${mondo}/mondo-simple.obo")
    def mondoSha = installValue(planned, 'mondo_sha256', '')
    def overrides = "${projectDir}/data/mondo/disease_scope_overrides.tsv"
    def raw = "data/hpo/raw/${hpo}"
    def obo = "data/mondo/mondo-simple.${mondo}.obo"
    def registry = 'data/mondo/disease_scope.tsv.gz'
    def assertions = 'data/hpo/genes_to_phenotype.assertions.tsv.gz'
    if (lineage.startsWith('adopted-'))
        installAdopt(id, lineage)
    else
        installStart(task) + """
        mkdir -p ${raw} data/mondo
        curl -fL --retry 5 --retry-delay 10 -o ${raw}/genes_to_phenotype.txt ${installWord(genesUrl)}
        curl -fL --retry 5 --retry-delay 10 -o ${raw}/phenotype.hpoa ${installWord(hpoaUrl)}
        curl -fL --retry 5 --retry-delay 10 -o ${obo} ${installWord(mondoUrl)}
        expected=${installWord(mondoSha)}
        if [[ -n \$expected ]]; then
            found=\$(sha256sum ${obo} | awk '{print \$1}')
            [[ \$found == "\$expected" ]] || { echo "MONDO SHA256 mismatch: expected \$expected, found \$found" >&2; exit 1; }
        fi
        "\$PRIVA_BIN" hpo-assertions build --genes-to-phenotype ${raw}/genes_to_phenotype.txt \\
            --phenotype-hpoa ${raw}/phenotype.hpoa --output ${raw}/genes_to_phenotype.assertions.base.tsv.gz
        overrides=()
        [[ -s ${installWord(overrides)} ]] && overrides=(--manual-overrides ${installWord(overrides)})
        "\$PRIVA_BIN" mondo-scope build --mondo-obo ${obo} --hpo-assertions ${raw}/genes_to_phenotype.assertions.base.tsv.gz \\
            --registry-output ${registry} --annotated-assertions-output ${assertions} "\${overrides[@]}"
        grep -q '^data-version: mondo/' ${obo} || { echo 'MONDO OBO has no data-version header' >&2; exit 1; }
        grep -q '^id: MONDO:0003847\$' ${obo} || { echo 'MONDO hereditary disease root MONDO:0003847 is missing' >&2; exit 1; }
        "\$PRIVA_BIN" mondo-scope validate --registry ${registry} --hpo-assertions ${assertions}
        ${installRecord(id, lineage, receipts, [hpo_assertions: assertions, mondo_obo: obo, mondo_disease_scope_registry: registry])}
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
