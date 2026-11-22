include { installAdopt; installRecord; installStart; installStub; installValue; installWord } from './common'

// F5: the GoFCards exact gain-of-function variant cache
// (gofcards_exact_gof_cache). The repository ships it with its ClinVar
// conditions (data/gofcards/gofcards_exact_gof.json.gz), and the step records
// that file, which step F7 then records as it is.
//
// Rebuilding it (PRIVA_FORCE_GOFCARDS_CACHE or PRIVA_FORCE_ALL_CACHES) writes
// the normalized cache, without ClinVar conditions, for F7 to inject them:
//
//   1. priva installer-prerequisites --stage gofcards checks the VEP cache of
//      step A2 (the install's assembly) and the FASTAs and chain you set,
//      before anything is downloaded.
//   2. priva reference-resources install --scope gofcards takes both
//      assemblies' FASTAs and UCSC's forward hg19-to-hg38 chain: the ones you
//      set (gofcards_hg19_fasta, gofcards_hg38_fasta,
//      gofcards_hg19_to_hg38_chain), else the reference of step A1 for the
//      install's assembly, else UCSC's pinned files, downloaded with their
//      MD5 checked.
//   3. vep_install (Ensembl's INSTALL.pl, --AUTO c --CONVERT) downloads the
//      VEP cache of the other assembly for vep_installed_version; its
//      release, vep's and its Ensembl API's must all be that release.
//   4. installer-prerequisites --stage gofcards checks the FASTAs, the chain
//      and both caches; priva gofcards-exact-gof downloads the GoFCards
//      workbook, asks GoFCards' summary endpoint about each allele, and runs
//      bcftools norm, vep on each assembly and CrossMap; priva
//      gofcards-cache-validate checks the cache.
//
// The step keeps the build folder (data/gofcards/build_work). It removes the
// folder of the two VEP caches (step A2's linked, the other downloaded) and
// any reference it downloaded once the cache is built; the cache names that
// folder in its provenance (vep_cache_dir), as it does a downloaded chain
// (liftover_chain).
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
    def adopted = lineage.startsWith('adopted-')
    def configured = adopted ? '' : installValue(planned, 'assembly', '')
    def selected = configured in ['hg19', 'GRCh37'] ? 'hg19' : (configured in ['hg38', 'GRCh38'] ? 'hg38' : '')
    def genome = selected == 'hg19' ? 'GRCh37' : 'GRCh38'
    def other = selected == 'hg19' ? 'GRCh38' : 'GRCh37'
    def release = adopted ? '' : installValue(planned, 'vep_installed_version', '')
    def shipped = "${projectDir}/data/gofcards/gofcards_exact_gof.json.gz"
    def review = "${projectDir}/data/gofcards/gofcards_mechanism_reviews.tsv"
    def hgnc = "${projectDir}/data/hgnc/hgnc_non_alt_loci_set_sha256_890b18743cc8.tsv"
    def link = ['ref_genome', 'vep_cache_dir', 'gofcards_hg19_fasta', 'gofcards_hg38_fasta', 'gofcards_hg19_to_hg38_chain']
        .collectEntries { key -> [(key): names[key] ? "staged/${names[key]}".toString() : ''] }
    def supplied = ['gofcards_hg19_fasta', 'gofcards_hg38_fasta', 'gofcards_hg19_to_hg38_chain']
        .collect { key -> "${key}=${installWord(link[key])}" }.join('\n        ')
    def work = 'data/gofcards/build_work'
    def target = "${work}/gofcards_exact_gof.normalized.json.gz".toString()
    if (adopted)
        installAdopt(id, lineage)
    else
        installStart(task) + """
        shipped=${installWord(shipped)}
        if [[ \${PRIVA_FORCE_GOFCARDS_CACHE:-0} != 1 && \${PRIVA_FORCE_ALL_CACHES:-0} != 1 ]]; then
            [[ -f \$shipped && ! -L \$shipped && -s \$shipped ]] || { echo "GoFCards: the shipped cache is missing or empty: \$shipped" >&2; exit 1; }
            echo "Reusing nonempty bundled gofcards_exact_gof_cache: \$shipped"
            ${installRecord(id, lineage, receipts, [:])} --key "gofcards_exact_gof_cache=\$shipped"
            exit 0
        fi
        set -o pipefail
        here=\$(pwd -P)
        root=${installWord(projectDir)}
        release=${installWord(release)}
        [[ -n ${installWord(selected)} ]] || { echo ${installWord("GoFCards: select hg19/GRCh37 or hg38/GRCh38, not ${configured}")} >&2; exit 2; }
        [[ \$release =~ ^[0-9]+\$ ]] || { echo 'GoFCards requires the configured VEP cache directory and numeric release' >&2; exit 1; }
        review=${installWord(review)}
        hgnc=${installWord(hgnc)}
        [[ -f \$review && -s \$review ]] || {
            echo "ERROR: GoFCards rebuild requires the selected nonempty mechanism review: \$review; restore that input." >&2; exit 1; }
        # Each file where its step or you keep it (the target of its staged link).
        resolved() { [[ -z \$1 ]] || readlink -f -- "\$1"; }
        reference=\$(resolved ${installWord(link.ref_genome)})
        selected_cache=\$(resolved ${installWord(link.vep_cache_dir)})
        ${supplied}
        caches=\$here/vep_cache
        quoted() { [[ -z \$1 ]] && printf '""' || printf '"%s"' "\$1"; }
        {
            printf 'assembly: %s\\n' ${installWord(selected)}
            printf 'cache_root: %s\\n' "\$(quoted "\$here/references")"
            printf 'ref_genome: %s\\n' "\$(quoted "\$reference")"
            printf 'gofcards_hg19_fasta: %s\\n' "\$(quoted "\$(resolved "\$gofcards_hg19_fasta")")"
            printf 'gofcards_hg38_fasta: %s\\n' "\$(quoted "\$(resolved "\$gofcards_hg38_fasta")")"
            printf 'gofcards_hg19_to_hg38_chain: %s\\n' "\$(quoted "\$(resolved "\$gofcards_hg19_to_hg38_chain")")"
            printf 'vep_cache_dir: %s\\n' "\$(quoted "\$caches")"
            printf 'vep_installed_version: "%s"\\n' "\$release"
        } > gofcards.yaml

        # 1. Step A2's VEP cache for ${genome}, and the FASTAs and chain you set,
        # checked before any download.
        mkdir -p "\$caches/homo_sapiens_merged"
        cache() { printf '%s/homo_sapiens_merged/%s_%s' "\$1" "\$release" "\$2"; }
        [[ -n \$selected_cache && -d \$(cache "\$selected_cache" ${genome}) ]] || {
            echo "GoFCards: step A2 holds no VEP \$release ${genome} cache" >&2; exit 1; }
        ln -s "\$(cache "\$selected_cache" ${genome})" "\$(cache "\$caches" ${genome})"
        "\$PRIVA_BIN" installer-prerequisites gofcards.yaml --stage gofcards --vep-version "\$release" \\
            --allow-pinned-reference-defaults --source-root "\$root" --allow-missing-merged-caches
        [[ -s \$hgnc ]] || { echo "ERROR: HGNC complete-set table not found: \$hgnc" >&2; exit 1; }
        for tool in vep bcftools CrossMap; do
            type -P \$tool > /dev/null || { echo "ERROR: \$tool is not on PATH" >&2; exit 1; }
        done

        # 2. Both FASTAs and the forward chain.
        "\$PRIVA_BIN" reference-resources install gofcards.yaml --data-dir "\$root/data" --scope gofcards \\
            --write-config --use-cache-root > references.tsv
        chosen() { awk -F '\\t' -v key="\$1" '\$1 == key {print \$2}' references.tsv; }
        hg19_fasta=\$(chosen gofcards_hg19_fasta)
        hg38_fasta=\$(chosen gofcards_hg38_fasta)
        chain=\$(chosen gofcards_hg19_to_hg38_chain)

        # 3. The VEP cache of ${other} from vep_install, whose release, vep's and
        # its Ensembl API's must be the configured one.
        for tool in vep vep_install perl; do
            type -P \$tool > /dev/null || { echo "GoFCards: \$tool is not on PATH" >&2; exit 1; }
        done
        vep_dir=\$(dirname "\$(realpath "\$(type -P vep)")")
        installer_dir=\$(dirname "\$(realpath "\$(type -P vep_install)")")
        actual=\$(vep --help 2>&1 | awk '/ensembl-vep[[:space:]]/ {print \$NF; exit}')
        actual=\${actual%%.*}
        api=\$(perl -I"\$vep_dir/modules" -I"\$vep_dir" -MBio::EnsEMBL::ApiVersion \\
            -e 'print Bio::EnsEMBL::ApiVersion::software_version()')
        installer=\$(perl -I"\$installer_dir/modules" -MBio::EnsEMBL::VEP::Utils=get_version_data \\
            -e 'my \$v = get_version_data(\$ARGV[0]."/.version")->{"ensembl-vep"}{"release"}; \$v =~ s{^release/}{}; print \$v' \\
            "\$installer_dir")
        [[ \$actual == "\$release" && \$api == "\$release" && \$installer == "\$release" ]] || {
            echo "GoFCards: requested VEP \$release; executable=\$actual, API=\$api, installer=\$installer" >&2; exit 1; }
        destdir=\$(perl -e 'print join("\\n", @INC);' | grep site_perl | head -1)
        [[ -n \$destdir ]] || { echo 'GoFCards: perl has no site_perl folder' >&2; exit 1; }
        vep_install -d "\$destdir" --AUTO c --NO_UPDATE --CACHE_VERSION "\$release" -s homo_sapiens_merged \\
            --NO_HTSLIB --NO_BIOPERL --CONVERT --CACHEDIR "\$caches" --ASSEMBLY ${other}

        # 4. The checks with everything in place, the build and its check.
        "\$PRIVA_BIN" installer-prerequisites gofcards.yaml --stage gofcards --vep-version "\$release" --source-root "\$root"
        mkdir -p ${work}
        echo "Building GoFCards exact GOF cache -> ${target}"
        "\$PRIVA_BIN" gofcards-exact-gof --workdir "\$here/${work}" --out-json "\$here/${target}" \\
            --hg19-fasta "\$hg19_fasta" --hg38-fasta "\$hg38_fasta" --chain "\$chain" \\
            --vep-cache-dir "\$caches" --vep-cache-version "\$release" --hgnc-table "\$hgnc" \\
            --mechanism-review-tsv "\$review" --stats-json "\$here/${work}/gofcards_build_stats.json" 2>&1 | tee ${work}/gofcards_build.log
        "\$PRIVA_BIN" gofcards-cache-validate ${target}
        rm -rf "\$caches" "\$here/references"
        ${installRecord(id, lineage, receipts, [gofcards_exact_gof_cache: target])}
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
