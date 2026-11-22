include { installAdopt; installRecord; installSetting; installStart; installStub; installWord } from './common'

// C1: the gnomAD v4.1 chromosome family of the assembly (gnomad_vcf_chrX names
// its chrX member; the other members are named alike beside it). For GRCh38
// the downloads of steps C1_chr1 to C1_chrY are linked into data/gnomAD with
// their indexes. For GRCh37 each download is lifted over with CrossMap (the
// chain hg38_hg19_chain: a file shipped in data/liftover/ or an absolute
// path, and the reference of step A1), sorted with bcftools and split by the
// chromosome its records land on; the pieces of each chromosome are
// concatenated (bcftools concat -a -d exact), sorted and indexed, as
// gnomAD_liftover did. Each member must read and hold a record.
process GNOMAD {
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
    def assembly = adopted ? '' : installSetting(planned, 'assembly')
    def chromosomes = (1..22).collect { n -> "chr${n}".toString() } + ['chrX', 'chrY']
    def downloads = chromosomes.collect { c -> names["gnomad_download_${c}".toString()] ? "staged/${names["gnomad_download_${c}".toString()]}".toString() : '' }
    def selected = planned.settings?.hg38_hg19_chain ? planned.settings.hg38_hg19_chain.toString() : ''
    def chain = selected in ['', 'null', 'hg38ToHg19.over.chain.gz'] ? "${projectDir}/data/liftover/hg38ToHg19.over.chain.gz".toString()
        : (selected.startsWith('/') ? selected : "${projectDir}/data/liftover/${selected}".toString())
    def referenceLink = names.ref_genome ? "staged/${names.ref_genome}".toString() : ''
    // The download of each chromosome, linked into data/gnomAD (GRCh38) or
    // into hg38/ for CrossMap, which reads a gzip file by its .gz name.
    def linkData = chromosomes.collect { c ->
        "fetched ${installWord(downloads[chromosomes.indexOf(c)])} data/gnomAD/gnomad.joint.v4.1.sites.${c}.vcf.bgz"
    }.join('\n        ')
    def linkLift = chromosomes.collect { c ->
        "fetched ${installWord(downloads[chromosomes.indexOf(c)])} hg38/gnomad.joint.v4.1.sites.${c}.vcf.gz"
    }.join('\n        ')
    def common = """
        # A member reads and holds a record (check_vcf_validity).
        has_record() { [[ -n \$(bcftools view -H "\$1" 2> /dev/null | head -n 1) ]]; }
        fetched() {
            local source
            [[ -n \$1 ]] && source=\$(readlink -f "\$1") && [[ -s \$source && -s \$source.tbi ]] || {
                echo "gnomAD: the download \$2 of step C1_chr* is missing" >&2; exit 1; }
            ln -sfn "\$source" "\$2"
            ln -sfn "\$source.tbi" "\$2.tbi"
            # An index not newer than its VCF is written again, here.
            [[ \$source.tbi -nt \$source ]] || { rm -f "\$2.tbi"; tabix -f -p vcf "\$2"; }
            has_record "\$2" || { echo "gnomAD: \$2 holds no record" >&2; exit 1; }
        }
        """
    if (adopted)
        installAdopt(id, lineage)
    else if (downloads.any { d -> !d })
        "echo 'gnomAD: a chromosome download of steps C1_chr1 to C1_chrY is missing' >&2; exit 1"
    else if (assembly in ['hg38', 'GRCh38'])
        installStart(task) + common.stripIndent() + """
        mkdir -p data/gnomAD
        ${linkData}
        ${installRecord(id, lineage, receipts, [gnomad_vcf_chrX: 'data/gnomAD/gnomad.joint.v4.1.sites.chrX.vcf.bgz'])}
        """.stripIndent()
    else if (assembly in ['hg19', 'GRCh37'])
        installStart(task) + common.stripIndent() + """
        set -o pipefail
        chain=${installWord(chain)}
        [[ -s \$chain ]] || { echo "gnomAD: the chain hg38_hg19_chain is missing or empty: \$chain" >&2; exit 1; }
        chain=\$(readlink -f "\$chain")
        reference=${installWord(referenceLink)}
        [[ -n \$reference ]] && reference=\$(readlink -f "\$reference") && [[ -s \$reference ]] || {
            echo 'gnomAD: the GRCh37 reference of step A1 is missing' >&2; exit 1; }
        mkdir -p hg38 hg19 tmp data/gnomAD
        ${linkLift}
        # gnomAD_liftover_per_chromosome: lift one chromosome, then split it by
        # the chromosome each record lands on.
        lift() {
            local chrom=\$1 out=hg19/\$1
            local mixed=\$out/gnomad.joint.v4.1.sites.mixed.hg19.vcf
            mkdir -p "\$out"
            CrossMap vcf --chromid l "\$chain" hg38/gnomad.joint.v4.1.sites.\$chrom.vcf.gz "\$reference" "\$mixed"
            bcftools sort -Oz -o "\$mixed.gz" "\$mixed"
            bcftools index -f -t "\$mixed.gz"
            rm "\$mixed"
            local landed
            landed=\$(bcftools query -f '%CHROM\\n' "\$mixed.gz" | sort - | uniq -)
            [[ -n \$landed ]] || { echo "gnomAD: \$chrom lifted to GRCh37 holds no record" >&2; return 1; }
            for target in \$landed; do
                bcftools view -Oz -o "\$out/gnomad.joint.v4.1.sites.\$target.hg19.vcf.gz" "\$mixed.gz" "\$target"
                bcftools index -f -t "\$out/gnomad.joint.v4.1.sites.\$target.hg19.vcf.gz"
            done
        }
        export -f lift
        export chain reference
        printf '%s\\n' ${chromosomes.join(' ')} | xargs -P ${task.cpus} -I{} bash -euo pipefail -c 'lift {}'
        # The pieces of each GRCh37 chromosome, as bcftools_concatvcfs joined them.
        for chrom in ${chromosomes.join(' ')}; do
            target=data/gnomAD/gnomad.joint.v4.1.sites.hg19.\$chrom.vcf.gz
            shopt -s nullglob
            pieces=( hg19/*/*\$chrom.*vcf.gz )
            shopt -u nullglob
            [[ \${#pieces[@]} -gt 0 ]] || { echo "gnomAD: no lifted records land on \$chrom" >&2; exit 1; }
            for piece in "\${pieces[@]}"; do
                [[ \$piece.tbi -nt \$piece ]] || bcftools index -t -f "\$piece"
                has_record "\$piece" || { echo "gnomAD: \$piece holds no record" >&2; exit 1; }
            done
            if [[ \${#pieces[@]} -eq 1 ]]; then
                bcftools view -Oz -o "\$target" "\${pieces[0]}"
                bcftools index -f -t "\$target"
            else
                printf '%s\\n' "\${pieces[@]}" > tmp/\$chrom.list
                bcftools concat -o "\$target" -a --threads 1 -d exact -Oz -f tmp/\$chrom.list
                tabix -f -p vcf "\$target"
                bcftools sort --temp-dir tmp -o "\${target%.vcf.gz}.sorted.vcf.gz" -Oz "\$target"
                mv "\${target%.vcf.gz}.sorted.vcf.gz" "\$target"
                bcftools index -f -t "\$target"
            fi
            has_record "\$target" || { echo "gnomAD: \$target holds no record" >&2; exit 1; }
        done
        rm -rf hg38 hg19 tmp
        ${installRecord(id, lineage, receipts, [gnomad_vcf_chrX: 'data/gnomAD/gnomad.joint.v4.1.sites.hg19.chrX.vcf.gz'])}
        """.stripIndent()
    else
        "echo 'gnomAD: unsupported assembly ${assembly}' >&2; exit 1"
    stub:
    installStub(id, lineage, receipts, task)
}
