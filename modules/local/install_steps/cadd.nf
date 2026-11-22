include { installAdopt; installFolder; installRecord; installStart; installStub; installValue; installWord } from './common'

// E1: upstream CADD-scripts and its data for one assembly (cadd_base_dir), as
// CADD_conda_install, CADD_install and CADD_upstream_check installed them:
//
//   1. scripts/cadd_conda.py finds a conda new enough for CADD's Snakemake,
//      or makes CADD's own frontend (data/CADD/conda_frontend) with the first
//      of micromamba, mamba or conda it finds (MAMBA_EXE, PATH, CONDA_EXE).
//   2. curl downloads the CADD-scripts release (cadd_zip_download_url), which
//      unzip unpacks into data/CADD, then the genome annotations (unpacked
//      with tar) and the prescored SNV and indel tables with their tabix
//      indexes; each MD5 you set is checked.
//   3. scripts/cadd_scripts_install.py puts PriVA's faster extract_scored.py
//      and MMSplice.py in place of the reviewed upstream files;
//      scripts/cadd_prescore_index.py and priva cadd-prescore-compact write
//      each prescored table's position index and compact copy.
//   4. CADD.sh runs once on its bundled test input, which makes its rule
//      environments, and must reproduce the variants and PHRED scores of its
//      expected output.
//
// Conda environments name their own folder in their files, so they cannot be
// moved: this step builds in its step folder (resource_root/E1/LINEAGE/data,
// linked as data/), where the CADD annotation finds them, instead of moving
// its results there afterwards; storeDir keeps only receipt.json. The folder
// stays writable (registry "writable"): CADD adds rule environments and
// MMSplice's exon tables to it when it runs.
process CADD {
    tag "${id} ${lineage}"
    label 'install'
    storeDir "${params.resource_root}/${id}/${lineage}"
    input:
    tuple val(id), val(lineage), path(plan, stageAs: 'plan.json'), path(receipts, stageAs: 'inputs/*/receipt.json'), val(planned), val(names), path(staged, stageAs: 'staged/*')
    output:
    path('receipt.json'), emit: receipt
    script:
    def adopted = lineage.startsWith('adopted-')
    def configured = adopted ? '' : installValue(planned, 'assembly', '')
    def genome = configured in ['hg19', 'GRCh37'] ? 'GRCh37' : (configured in ['hg38', 'GRCh38'] ? 'GRCh38' : '')
    def version = installValue(planned, 'cadd_version', '')
    def setting = ['cadd_zip_download_url', "cadd_${genome}_anno_url", "cadd_${genome}_anno_md5",
                   "cadd_${genome}_snv_anno_url", "cadd_${genome}_snv_anno_md5", "cadd_${genome}_snv_anno_tbi_md5",
                   "cadd_${genome}_indel_anno_url", "cadd_${genome}_indel_anno_md5", "cadd_${genome}_indel_anno_tbi_md5"]
        .collectEntries { key -> [(key.toString()): installValue(planned, key.toString(), '')] }
    def required = ['cadd_version'] + setting.keySet().findAll { key -> key.endsWith('_url') }
    def missing = required.findAll { key -> !(key == 'cadd_version' ? version : setting[key]) }
    def tables = ['snv', 'indel'].collect { kind ->
        def url = setting["cadd_${genome}_${kind}_anno_url".toString()]
        def name = url ? url.tokenize('/').last() : ''
        def md5 = setting["cadd_${genome}_${kind}_anno_md5".toString()]
        def tbi = setting["cadd_${genome}_${kind}_anno_tbi_md5".toString()]
        "fetch ${installWord(url)} \"\$prescored/\"${installWord(name)} ${installWord(md5)}\n" +
            "        fetch ${installWord(url + '.tbi')} \"\$prescored/\"${installWord(name + '.tbi')} ${installWord(tbi)}"
    }.join('\n        ')
    def anno = "${genome}_${version}".toString()
    def folder = installFolder(id, lineage)
    if (adopted)
        installAdopt(id, lineage)
    else if (!genome)
        "echo 'CADD: unsupported assembly ${configured}' >&2; exit 1"
    else if (!missing.isEmpty())
        "echo ${installWord('CADD: set ' + missing.join(', ') + ' in the install configuration')} >&2; exit 1"
    else
        """
        set -o pipefail
        rm -rf ${installWord(folder + '/data')}
        mkdir -p ${installWord(folder + '/data')}
        ln -sfn ${installWord(folder + '/data')} data
        """.stripIndent() + installStart(task) + """
        here=\$(pwd -P)
        root=${installWord(projectDir)}
        parent=${installWord(folder + '/data/CADD')}
        mkdir -p "\$parent"
        # A download and the MD5 you set for it (none when unset).
        fetch() {
            mkdir -p "\$(dirname "\$2")"
            curl -fL --retry 5 --retry-delay 10 -o "\$2" "\$1"
            [[ -z \$3 || \$3 == null || \$(md5sum "\$2" | cut -d ' ' -f 1) == "\$3" ]] || {
                echo "CADD: the MD5 of \$1 differs from \$3" >&2; exit 1; }
        }

        # 1. A conda for CADD's Snakemake, or CADD's own frontend.
        manager=''
        for candidate in "\${MAMBA_EXE:-}" "\$(type -P micromamba)" "\$(type -P mamba)" "\${CONDA_EXE:-}" "\$(type -P conda)"; do
            [[ -n \$candidate && -x \$candidate ]] || continue
            case \${candidate##*/} in micromamba*|mamba*|conda*) manager=\$candidate; break ;; esac
        done
        [[ -n \$manager ]] || {
            echo 'CADD: install micromamba, mamba or conda and expose it through PATH, MAMBA_EXE or CONDA_EXE' >&2; exit 1; }
        python "\$root/scripts/cadd_conda.py" install --frontend-prefix "\$parent/conda_frontend" --manager "\$manager" > /dev/null

        # 2. CADD-scripts, the genome annotations and the prescored tables.
        mkdir -p unpack
        fetch ${installWord(setting.cadd_zip_download_url)} unpack/CADD-scripts.zip ''
        unzip -q unpack/CADD-scripts.zip -d unpack
        unpacked=\$(find unpack/ -mindepth 1 -maxdepth 1 -type d -name 'CADD-scripts-*' -print -quit)
        [[ -f \$unpacked/CADD.sh ]] || { echo 'CADD: the CADD-scripts release holds no CADD-scripts folder with CADD.sh' >&2; exit 1; }
        base=\$parent/\$(basename "\$unpacked")
        mv -- "\$unpacked" "\$base"
        rm -rf unpack
        fetch ${installWord(setting["cadd_${genome}_anno_url".toString()])} "\$base/data/annotations/${anno}.tar.gz" ${installWord(setting["cadd_${genome}_anno_md5".toString()])}
        mkdir -p unpack
        tar -xzf "\$base/data/annotations/${anno}.tar.gz" -C unpack
        [[ -d unpack/${anno} ]] || { echo 'CADD: the genome annotations hold no folder ${anno}' >&2; exit 1; }
        mv -- unpack/${anno} "\$base/data/annotations/${anno}"
        rm -rf unpack "\$base/data/annotations/${anno}.tar.gz"
        prescored=\$base/data/prescored/${anno}/incl_anno
        ${tables}

        # 3. PriVA's CADD scripts, and each prescored table's index and compact copy.
        python "\$root/scripts/cadd_scripts_install.py" --source "\$base"
        for table in "\$prescored"/*.tsv.gz; do
            python "\$root/scripts/cadd_prescore_index.py" build "\$table" --threads ${task.cpus}
            "\$PRIVA_BIN" cadd-prescore-compact build "\$table" --threads ${task.cpus}
        done

        # 4. CADD.sh on its bundled test input: the same variants and PHRED.
        expected=\$base/test/output_inclAnno_${genome}.tsv.gz
        [[ -x \$base/CADD.sh && -f \$base/test/input.vcf.gz && -f \$expected ]] || {
            echo "CADD: CADD.sh or its bundled test files are missing under \$base" >&2; exit 1; }
        path=\$(python "\$root/scripts/cadd_conda.py" check --frontend-prefix "\$parent/conda_frontend" --print-path)
        exons=\$base/data/annotations/${genome}_v1.7/mmsplice
        [[ -d \$exons ]] || exons=''
        mkdir -p check/tmp
        PATH=\$path TMPDIR=\$here/check/tmp PRIVA_CADD_MMSPLICE_EXON_CACHE=\$exons \\
            "\$base/CADD.sh" -c ${task.cpus} -a -m -g ${genome} -o "\$here/check/output.tsv.gz" "\$base/test/input.vcf.gz"
        scores() { zcat "\$1" | awk -F '\\t' '\$1 !~ /^#/ { print \$1, \$2, \$3, \$4, \$NF }'; }
        raw() { zcat "\$1" | awk -F '\\t' '\$1 !~ /^#/ { print \$(NF-1) }'; }
        diff <(scores "\$expected") <(scores check/output.tsv.gz) > /dev/null || {
            echo "CADD: CADD.sh's scores differ from its bundled ${genome} test output (kept in \$here/check)" >&2; exit 1; }
        echo "CADD.sh reproduced its bundled ${genome} test PHRED scores (RawScore differs on \$(diff <(raw "\$expected") <(raw check/output.tsv.gz) | grep -c '^<') rows)"
        rm -rf check
        ${installRecord(id, lineage, receipts, [:])} --key "cadd_base_dir=data/CADD/\$(basename "\$base")"
        """.stripIndent()
    stub:
    if (lineage.startsWith('adopted-'))
        installAdopt(id, lineage)
    else
        """
        rm -rf ${installWord(installFolder(id, lineage) + '/data')}
        mkdir -p ${installWord(installFolder(id, lineage) + '/data')}
        ln -sfn ${installWord(installFolder(id, lineage) + '/data')} data
        """.stripIndent() + installStub(id, lineage, receipts, task)
}
