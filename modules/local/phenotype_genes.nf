// One shell word: the value in single quotes, each single quote in it written
// as '"'"'.
def shellWord(value) {
    "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

// One family's genes ranked for its HPO terms by `hpo-gene-ranking rank-genes`
// (repository llm_gene_reranker). PRIORITIZE_ACMG turns the family's
// gene_ranked.tsv into each gene's phenotype_gene_index. The terms arrive as
// values, so a task's identity counts only its own family's terms, the
// program's bytes (its sha256) and the reference folder's path.
process RANK_PHENOTYPE_GENES {
    tag "${family}"
    label 'phenotype_ranking'
    // The task writes its tables into a folder named after the family, which
    // keeps that name when published (the strict parser evaluates the
    // publishDir path once, before any task, so it cannot name the family).
    publishDir "${params.outdir}/phenotype_gene_ranking", mode: 'copy', overwrite: false
    input:
    tuple val(family), val(terms), val(ranking)
    output:
    tuple val(family), path("${family}/*.gene_ranked.tsv"), emit: ranking
    // The other tables of rank-genes, published beside the ranking.
    path("${family}/*.{frontier_p_terms,p_pairs,frontier_gene_evidence,gene_summary,gene_score_contributions}.tsv"), emit: tables
    path("${family}/run_summary.tsv"), emit: summary
    script:
    def rows = terms.collect { term -> "${shellWord(family)} ${shellWord(term)}" }.join(' ')
    """
    printf 'sample_id\\thpo_term\\n' > hpo_terms.tsv
    printf '%s\\t%s\\n' ${rows} >> hpo_terms.tsv
    '${ranking.binary}' rank-genes --reference-dir '${ranking.reference}' --hpo-term-table hpo_terms.tsv --sample-id ${shellWord(family)} --output-dir ${shellWord(family)}
    """
    stub:
    """
    mkdir ${shellWord(family)}
    printf 'sample_id\\tgene_symbol\\tncbi_gene_ids\\ttotal_score\\n' > '${family}/${family}.gene_ranked.tsv'
    for table in frontier_p_terms p_pairs frontier_gene_evidence gene_summary gene_score_contributions; do : > "${family}/${family}.\$table.tsv"; done
    : > '${family}/run_summary.tsv'
    """
}
