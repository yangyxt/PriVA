include { installAdopt; installRecord; installSetting; installStart; installStub; installValue; installWord } from './common'

// D1: the ClinVar release of step D0 annotated by VEP and gnomAD
// (clinvar_vcf), its amino-acid and splice statistics (clinvar_aa_stat,
// clinvar_splice_stat) and the gene statistics (clinvar_gene_stat), as
// priva clinvar-gene-refresh install made them through ClinVar_VCF_deploy,
// ClinVar_AA_stat and ClinVar_Gene_stat(_regenerate):
//
//   1. bcftools renames the release's contigs (as D1_split did); its records
//      must be those D1_split split (the hash in split.json).
//   2. The VEP outputs of steps D1_chr1 to D1_chrOther are gathered: their
//      headers must be the same apart from VEP's command line and start time;
//      bcftools sorts them into one VCF, which must hold at least one record
//      and no more than the release (VEP leaves out a record whose ALT is .).
//   3. The VEP release and assembly in its header (priva clinvar-gene-refresh
//      check-header) and the fields the statistics read (priva
//      clinvar-contract header) are checked; priva task step --step gnomad
//      adds gnomAD's joint frequencies (C1's family) and the checks run again.
//   4. priva clinvar-aa-splice builds the two statistics; clinvar-contract
//      checks each archive.
//   5. priva clinvar-gene-refresh due decides whether the selected gene
//      statistics stay (the bundled ones until D0's release is three calendar
//      months newer than them) or priva clinvar-gene-stats builds new ones
//      from the release's variant_summary and MANE's summary.
process CLINVAR {
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
    def configured = adopted ? '' : installSetting(planned, 'assembly')
    def ucsc = configured in ['hg19', 'GRCh37'] ? 'hg19' : (configured in ['hg38', 'GRCh38'] ? 'hg38' : '')
    def vep_major = adopted ? '' : installSetting(planned, 'vep_installed_version')
    def seed = installValue(planned, 'priva_python_hash_seed', '0')
    def pieces = (1..22).collect { n -> "chr${n}".toString() } + ['chrX', 'chrY', 'chrOther']
    def links = ['clinvar_release_inputs', 'clinvar_vep_split', 'gnomad_vcf_chrX', 'spliceai_snv_prescore',
                 'spliceai_indel_prescore'].collectEntries { key -> [(key): names[key] ? "staged/${names[key]}".toString() : ''] }
    def outputs = pieces.collect { p ->
        def key = "clinvar_vep_${p}".toString()
        "[${p}]=${installWord(names[key] ? 'staged/' + names[key] : '')}"
    }.join(' ')
    def needed = [['clinvar_release_inputs', 'D0'], ['clinvar_vep_split', 'D1_split'], ['gnomad_vcf_chrX', 'C1']]
    def missing = needed.findAll { pair -> !names[pair[0]] }.collect { pair -> "${pair[0]} (step ${pair[1]})".toString() }
    def folder = 'data/ClinVar'
    def converted = "${folder}/clinvar.${ucsc}.vcf.gz".toString()
    def target = "${folder}/clinvar.${ucsc}.vep.vcf.gz".toString()
    def aa = "${folder}/clinvar.${ucsc}.vep.aa_change.stats.rkyv".toString()
    def splice = "${folder}/clinvar.${ucsc}.vep.splice_change.stats.rkyv".toString()
    def gene = "${folder}/clinvar_2star_stats.rkyv".toString()
    if (adopted)
        installAdopt(id, lineage)
    else if (!ucsc)
        "echo 'ClinVar: unsupported assembly ${configured}' >&2; exit 1"
    else if (!missing.isEmpty())
        "echo ${installWord('ClinVar: missing ' + missing.join(', '))} >&2; exit 1"
    else
        installStart(task) + """
        set -o pipefail
        here=\$(pwd -P)
        root=${installWord(projectDir)}
        release=\$(readlink -f ${installWord(links.clinvar_release_inputs)})
        manifest=\$(readlink -f ${installWord(links.clinvar_vep_split)})
        # C1's chrX member where C1 keeps it, beside the other chromosomes.
        chrx=\$(readlink ${installWord(links.gnomad_vcf_chrX)})
        snv=\$([[ -z ${installWord(links.spliceai_snv_prescore)} ]] || readlink -f ${installWord(links.spliceai_snv_prescore)})
        indel=\$([[ -z ${installWord(links.spliceai_indel_prescore)} ]] || readlink -f ${installWord(links.spliceai_indel_prescore)})
        source=\$(jq -r '.sources.vcf.path // empty' "\$release")
        summary=\$(jq -r '.sources.summary.path // empty' "\$release")
        mane=\$(jq -r '.sources.mane.path // empty' "\$release")
        [[ -s \$source && -s \$source.tbi ]] || { echo "ClinVar: D0's ClinVar VCF or its index is missing: \$source" >&2; exit 1; }
        mkdir -p ${folder}
        line() { local q="'"; printf "%s: '%s'\\n" "\$1" "\${2//\$q/\$q\$q}"; }
        {
        line assembly ${installWord(configured)}
        line vep_installed_version ${installWord(vep_major)}
        line threads ${task.cpus}
        line cache_root "\$here/data"
        line clinvar_vcf_dir "\$here/${folder}"
        line clinvar_vcf "\$here/${target}"
        line clinvar_gene_stat ${installWord(installValue(planned, 'clinvar_gene_stat', 'clinvar_2star_stats.pkl.gz'))}
        line clinvar_release_inputs "\$release"
        line gnomad_vcf_chrX "\$chrx"
        line spliceai_snv_prescore "\$snv"
        line spliceai_indel_prescore "\$indel"
        } > d1.yaml

        # 1. The release with UCSC contig names: the records D1_split split.
        bcftools annotate --rename-chrs "\$root/data/liftover/GRC_to_ucsc.contig.map.txt" "\$source" -Ou \\
            | bcftools sort -T "\$here/sort.XXXXXX" -Oz -o ${converted}
        tabix -f -p vcf ${converted}
        body=\$(bgzip -dc ${converted} | grep -v '^#' | sha256sum | cut -c1-64)
        [[ \$body == "\$(jq -r .body_sha256 "\$manifest")" ]] || {
            echo "ClinVar: the per-chromosome VEP outputs were made from other ClinVar records than D0's release" >&2; exit 1; }

        # 2. The per-chromosome VEP outputs, gathered and sorted.
        declare -A annotated=(${outputs})
        header_of() { ( set +o pipefail; bgzip -dc -- "\$1" | awk '/^#/ { print; next } { exit }' ); }
        first='' reference='' total=0 gathered=()
        while IFS=\$'\\t' read -r name records; do
            ((records > 0)) || continue
            piece=\${annotated[\$name]:-}
            [[ -n \$piece && -s \$piece ]] || { echo "ClinVar: the VEP output of \$name (step D1_\$name) is missing" >&2; exit 1; }
            piece=\$(readlink -f "\$piece")
            header=\$(header_of "\$piece" | grep -v '^##VEP-command-line=' | sed -E 's/ time="[^"]*"//' | sha256sum)
            if [[ -z \$first ]]; then
                first=\$piece reference=\$header
            elif [[ \$header != "\$reference" ]]; then
                echo "ClinVar: the VEP header of \$name differs from that of the first piece (\$first)" >&2; exit 1
            fi
            gathered+=("\$piece"); total=\$((total + records))
        done < <(jq -r '.pieces[] | "\\(.name)\\t\\(.records)"' "\$manifest")
        ((\${#gathered[@]})) || { echo 'ClinVar: no VEP output to gather' >&2; exit 1; }
        {
            header_of "\$first"
            for piece in "\${gathered[@]}"; do bgzip -dc -- "\$piece" | awk '!/^#/'; done
        } | bcftools sort -T "\$here/sort.XXXXXX" -Oz -o ${target} -
        tabix -f -p vcf ${target}
        written=\$(bcftools index -n ${target})
        ((written > 0 && written <= total)) || { echo "ClinVar: gathered \$written VEP records for \$total ClinVar records" >&2; exit 1; }
        echo "ClinVar: gathered \$written VEP records (from \$total ClinVar records, \${#gathered[@]} pieces)"

        # 3. VEP's release and assembly, the fields read later, then gnomAD.
        "\$PRIVA_BIN" clinvar-gene-refresh check-header --config "\$here/d1.yaml" --source-root "\$root" > /dev/null
        "\$PRIVA_BIN" clinvar-contract header --vcf ${target} --stage vep > /dev/null
        "\$PRIVA_BIN" task step --step gnomad --root "\$root" --config "\$here/d1.yaml" --vcf "\$here/${target}"
        "\$PRIVA_BIN" clinvar-contract header --vcf ${target} --stage complete > /dev/null
        "\$PRIVA_BIN" clinvar-gene-refresh check-header --config "\$here/d1.yaml" --source-root "\$root" > /dev/null

        # 4. Amino-acid and splice statistics.
        "\$PRIVA_BIN" clinvar-contract header --vcf ${target} --stage aa > /dev/null
        "\$PRIVA_BIN" clinvar-aa-splice build --vcf ${target} --aa-output ${aa}.partial \\
            --splice-output ${splice}.partial --threads ${task.cpus}
        "\$PRIVA_BIN" clinvar-contract cache --input ${aa}.partial --kind aa > /dev/null
        "\$PRIVA_BIN" clinvar-contract cache --input ${splice}.partial --kind splice > /dev/null
        mv -- ${aa}.partial ${aa}
        mv -- ${splice}.partial ${splice}

        # 5. Gene statistics: the selected ones, or new ones once they are due.
        decision=\$("\$PRIVA_BIN" clinvar-gene-refresh due --config "\$here/d1.yaml" --source-root "\$root")
        action=\$(jq -r .action <<< "\$decision")
        echo "ClinVar gene statistics: \$action; \$(jq -r .reason <<< "\$decision")"
        if [[ \$action == rebuild ]]; then
            [[ -s \$summary && -s \$mane ]] || { echo "ClinVar: D0's variant_summary or MANE summary is missing" >&2; exit 1; }
            "\$PRIVA_BIN" clinvar-contract header --vcf ${target} --stage gene > /dev/null
            "\$PRIVA_BIN" clinvar-gene-stats build --vcf ${target} --output ${folder}/clinvar_2star_stats.partial.rkyv \\
                --python-hash-seed ${installWord(seed)} --clinvar-summary "\$summary" --mane-summary "\$mane" \\
                --threads ${task.cpus} --config "\$here/d1.yaml"
            "\$PRIVA_BIN" clinvar-contract cache --input ${folder}/clinvar_2star_stats.partial.rkyv --kind gene > /dev/null
            mv -- ${folder}/clinvar_2star_stats.partial.rkyv ${gene}
            gene=${gene}
        else
            "\$PRIVA_BIN" installer-prerequisites "\$here/d1.yaml" --source-root "\$root" --stage clinvar
            gene=\$(jq -r .selected <<< "\$decision")
        fi
        ${installRecord(id, lineage, receipts, [clinvar_vcf: target, clinvar_aa_stat: aa, clinvar_splice_stat: splice, clinvar_vcf_dir: folder])} --key "clinvar_gene_stat=\$gene"
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
