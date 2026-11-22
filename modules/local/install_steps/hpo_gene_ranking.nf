include { installAdopt; installRecord; installStart; installStub; installWord } from './common'

// F10: the reference folder of the phenotype gene bonus
// (hpo_gene_ranking_reference_dir): raw/ (the public inputs), prepared/
// (term-to-gene tables) and hpo_graph/ (HPO nodes and edges). curl downloads
// the inputs as llm_gene_reranker's download_public_inputs.sh does (the same
// URLs and file names; gene2go kept to its human rows; the Human Protein Atlas
// tables optional) and logs each in raw/download_log.tsv; priva
// hpo-gene-ranking prepare-gene-evidence and build-hpo-graph build the tables.
// The manifest names the files where storeDir keeps them.
process HPO_GENE_RANKING {
    tag "${id} ${lineage}"
    label 'install_download'
    storeDir "${params.resource_root}/${id}/${lineage}"
    input:
    tuple val(id), val(lineage), path(plan, stageAs: 'plan.json'), path(receipts, stageAs: 'inputs/*/receipt.json'), val(planned), val(names), path(staged, stageAs: 'staged/*')
    output:
    path('receipt.json'), emit: receipt
    path('data'), emit: data
    script:
    def uniprot = 'https://rest.uniprot.org/uniprotkb/stream?query=(organism_id:9606)%20AND%20(reviewed:true)&fields=accession,gene_primary,gene_names,protein_name,cc_function,cc_subcellular_location,cc_tissue_specificity,cc_disease,xref_geneid&format=tsv&compressed=true'
    def downloads = [
        ['http://purl.obolibrary.org/obo/hp.owl', 'hp.owl'],
        ['http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt', 'phenotype_to_genes.txt'],
        ['http://purl.obolibrary.org/obo/go/go-basic.obo', 'go-basic.obo'],
        ['https://current.geneontology.org/annotations/goa_human.gaf.gz', 'goa_human.gaf.gz'],
        ['https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz', 'Homo_sapiens.gene_info.gz'],
        ['https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz', 'gene2go.gz'],
        ['https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt', 'NCBI2Reactome_All_Levels.txt'],
        ['https://reactome.org/download/current/ReactomePathwaysRelation.txt', 'ReactomePathwaysRelation.txt'],
        [uniprot, 'uniprot_human_reviewed.tsv.gz'],
        ['http://purl.obolibrary.org/obo/pr.obo', 'pr.obo'],
        ['https://www.informatics.jax.org/downloads/reports/MGI_GenePheno.rpt', 'MGI_GenePheno.rpt'],
        ['https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt', 'HOM_MouseHumanSequence.rpt'],
        ['https://www.orphadata.com/data/xml/en_product6.xml', 'orphadata_en_product6.xml'],
    ].collect { pair -> "fetch ${installWord(pair[0])} ${installWord(pair[1])} || failed=1" }.join('\n        ')
    def reference = 'data/hpo_gene_ranking'
    if (lineage.startsWith('adopted-'))
        installAdopt(id, lineage)
    else
        installStart(task) + """
        here=\$(pwd -P)
        reference=\$here/${reference}
        raw=\$reference/raw
        mkdir -p "\$raw" "\$reference/prepared" "\$reference/hpo_graph"
        printf 'file\\turl\\tdownloaded_at\\tbytes\\tsha256\\n' > "\$raw/download_log.tsv"
        failed=0
        # curl, then a line of raw/download_log.tsv.
        fetch() {
            local url=\$1 name=\$2 target=\$raw/\$2
            echo "Downloading HPO gene ranking input \$name"
            if curl -L --fail --retry 4 --retry-delay 5 --connect-timeout 30 -sS -o "\$target.part" "\$url"; then
                mv "\$target.part" "\$target"
                printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "\$name" "\$url" "\$(date -Is)" "\$(wc -c < "\$target" | tr -d ' ')" \\
                    "\$(sha256sum "\$target" | cut -c1-64)" >> "\$raw/download_log.tsv"
            else
                rm -f "\$target.part"
                echo "ERROR: HPO gene ranking download failed: \$name (\$url)" >&2
                return 1
            fi
        }
        ${downloads}
        for table in rna_tissue_consensus rna_single_cell_type rna_immune_cell_schmiedel; do
            fetch "https://www.proteinatlas.org/download/tsv/\$table.tsv.zip" "\$table.tsv.zip" \\
                && unzip -p "\$raw/\$table.tsv.zip" > "\$raw/\$table.tsv" \\
                || echo "Optional Human Protein Atlas table \$table not available"
        done
        for archive in goa_human.gaf Homo_sapiens.gene_info uniprot_human_reviewed.tsv; do
            if [[ -s \$raw/\$archive.gz ]]; then
                gzip -dc "\$raw/\$archive.gz" > "\$raw/\$archive.part" && mv "\$raw/\$archive.part" "\$raw/\$archive" || failed=1
            fi
        done
        if [[ -s \$raw/gene2go.gz ]]; then
            gzip -dc "\$raw/gene2go.gz" | awk -F '\\t' '/^#/ || \$1 == "9606"' > "\$raw/gene2go.part" \\
                && mv "\$raw/gene2go.part" "\$raw/gene2go" || failed=1
        fi
        rm -f "\$raw"/*.gz "\$raw"/*.zip
        [[ \$failed -eq 0 ]] || { echo 'ERROR: some required HPO gene ranking downloads failed' >&2; exit 1; }
        echo 'Building the HPO gene ranking term-to-gene tables'
        "\$PRIVA_BIN" hpo-gene-ranking prepare-gene-evidence --raw-dir "\$raw" --out-dir "\$reference/prepared"
        graph=\$reference/hpo_graph/hp_owl_hpo_only
        echo 'Building the HPO gene ranking HPO graph tables'
        "\$PRIVA_BIN" hpo-gene-ranking build-hpo-graph --hpo-owl "\$raw/hp.owl" --hpo-only --output-prefix "\$graph"
        [[ -s \$reference/prepared/manifest.json && -s \$graph.nodes.tsv && -s \$graph.edges.tsv ]] || {
            echo "ERROR: HPO gene ranking reference is incomplete: \$reference" >&2
            exit 1
        }
        ${installRecord(id, lineage, receipts, [hpo_gene_ranking_reference_dir: reference])} --relocate ${reference}/prepared/manifest.json
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
