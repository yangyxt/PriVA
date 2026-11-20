include { installAdopt; installCommit; installRecord; installSetting; installStart; installStub } from './common'

// B8: LOFTEE, VEP's LoF plugin (loftee_repo: github.com/konradjk/loftee, the
// master branch for GRCh37 and grch38 for GRCh38, at the commit the step is
// pinned to, without git's own folder), and the files it reads, downloaded
// with curl: the ancestral FASTA and its index (human_ancestor_fasta), the
// conservation database (loftee_conservation_file, decompressed with gunzip)
// and, for GRCh38, the GERP bigWig (gerp_bigwig). priva loftee-resources
// validate checks the FASTA's header and index and opens the database with
// SQLite; for GRCh38 priva conservation-resource validate reads the bigWig
// against the reference of step A1. A FASTA compressed with bgzip gets its
// block index (bgzip -r), which VEP's lookups need on a read-only folder.
process LOFTEE {
    tag "${id} ${lineage}"
    label 'install_download'
    storeDir "${params.resource_root}/${id}/${lineage}"
    input:
    tuple val(id), val(lineage), path(plan, stageAs: 'plan.json'), path(receipts, stageAs: 'inputs/*/receipt.json'), val(planned), val(names), path(staged, stageAs: 'staged/*')
    output:
    path('receipt.json'), emit: receipt
    path('data'), emit: data
    script:
    def assembly = installSetting(planned, 'assembly')
    def grch38 = assembly in ['hg38', 'GRCh38']
    def known = grch38 || assembly in ['hg19', 'GRCh37']
    def commit = installCommit(planned)
    def broad = 'https://personal.broadinstitute.org/konradk/loftee_data'
    def ancestorUrl = grch38 ? "${broad}/GRCh38/human_ancestor.fa.gz" : 'https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz'
    def sqlUrl = grch38 ? "${broad}/GRCh38/loftee.sql.gz" : "${broad}/GRCh37/phylocsf_gerp.sql.gz"
    def gerpUrl = "${broad}/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
    def repo = "data/LoFtee_repos/loftee-${assembly}"
    def cache = "data/LoFtee/${assembly}"
    def ancestor = "${cache}/human_ancestor.fa.gz"
    def sql = "${cache}/phylocsf_gerp.sql"
    def gerp = "${cache}/gerp_conservation_scores.homo_sapiens.bw"
    def keys = [loftee_repo: repo, human_ancestor_fasta: ancestor, loftee_conservation_file: sql]
    if (grch38) keys.gerp_bigwig = gerp
    if (lineage.startsWith('adopted-'))
        installAdopt(id, lineage)
    else if (!known)
        "echo 'LOFTEE: unsupported assembly ${assembly}' >&2; exit 1"
    else
        installStart(task) + """
        mkdir -p data/LoFtee_repos ${cache}
        git clone --quiet https://github.com/konradjk/loftee.git ${repo}
        git -C ${repo} checkout --quiet ${grch38 ? 'grch38' : 'master'}
        ${commit ? "git -C ${repo} checkout --quiet ${commit}" : ':'}
        git -C ${repo} rev-parse HEAD > ${repo}.commit
        rm -rf ${repo}/.git
        [[ -s ${repo}/LoF.pm ]] || { echo 'LOFTEE: LoF.pm is missing from the checkout' >&2; exit 1; }
        curl -fL --retry 5 --retry-delay 10 -o ${ancestor} ${ancestorUrl}
        curl -fL --retry 5 --retry-delay 10 -o ${ancestor}.fai ${ancestorUrl}.fai
        curl -fL --retry 5 --retry-delay 10 -o ${sql}.gz ${sqlUrl}
        gunzip -f ${sql}.gz
        ${grch38 ? "curl -fL --retry 5 --retry-delay 10 -o ${gerp} ${gerpUrl}" : ':'}
        gzip -t ${ancestor}
        "\$PRIVA_BIN" loftee-resources validate --ancestor ${ancestor} --conservation ${sql} --assembly ${assembly}
        ${grch38 ? """python=\$(type -P python) || { echo 'LOFTEE: no python on PATH (the Perl bigWig reader is looked for beside it)' >&2; exit 1; }
        "\$PRIVA_BIN" conservation-resource validate --input ${gerp} --assembly ${assembly} --reference staged/${names.ref_genome} --tool-prefix "\$(dirname "\$python")\"""" : ':'}
        if [[ \$(head -c 14 ${ancestor} | od -An -tx1 | tr -d ' \\n') == 1f8b0804????????????????4243 ]]; then
            bgzip -r ${ancestor}
        fi
        ${installRecord(id, lineage, receipts, keys)}
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
