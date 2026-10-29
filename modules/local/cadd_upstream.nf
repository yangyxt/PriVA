// Runs upstream CADD-scripts' own CADD.sh (its Snakemake workflow and data) on
// one VCF, as PriVA's Calculate_CADD does. main_install replaced its prescore
// lookup and MMSplice scripts with faster copies that write the same rows.
def quoteCadd(value) { "'" + value.toString().replace("'", "'\\''") + "'" }

process CADD_UPSTREAM {
    tag "${meta.id}:${meta.genome}"
    label 'cadd'
    publishDir "${params.outdir}", mode: 'copy', overwrite: false
    input:
    tuple val(meta), path(vcf)
    output:
    tuple val(meta), path("${meta.id}.cadd.tsv"), emit: scores
    script:
    def base = meta.cadd_base_dir
    """
    # CADD expects GRCh contig names; convert UCSC-style names first.
    bcftools annotate --rename-chrs ${quoteCadd(meta.source_root + '/data/liftover/ucsc_to_GRC.contig.map.tsv')} \\
        -Oz -o input.nochr.vcf.gz ${vcf}
    # PriVA's MMSplice copy reloads its exon table from beside CADD's GTF.
    exon_tables=${quoteCadd(base + '/data/annotations/' + meta.genome + '_v1.7/mmsplice')}
    [[ -d \$exon_tables ]] || exon_tables=""
    TMPDIR=\$PWD PRIVA_CADD_MMSPLICE_EXON_CACHE=\$exon_tables \\
        ${quoteCadd(base + '/CADD.sh')} -c ${task.cpus} -a -m -g ${meta.genome} \\
        -o cadd.tsv.gz input.nochr.vcf.gz
    # The consumer reads the column header as the first row.
    zcat cadd.tsv.gz | awk -F '\\t' '\$1 !~ /^##/' > ${meta.id}.cadd.tsv
    """
    stub:
    """
    printf '#Chrom\\tPos\\tRef\\tAlt\\tRawScore\\tPHRED\\tFeatureID\\n' > ${meta.id}.cadd.tsv
    """
}
