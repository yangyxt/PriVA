include { installAdopt; installRecord; installSetting; installStart; installStub; installWord } from './common'

// A2: Ensembl's merged VEP cache of the assembly for the VEP release
// vep_installed_version (vep_cache_dir), and the release's official plugins
// (vep_plugins_dir). The vep, its Ensembl API and its install helper
// (vep_install, Ensembl's INSTALL.pl) on PATH must all be that release.
// vep_install downloads and converts the cache (--AUTO cp --CONVERT) and
// fetches the plugins into a staging folder; the cache is checked (its
// info.txt names the assembly, homo_sapiens and tabix variants), and priva
// vep-plugins install records the plugins with their hashes, then
// ensure-utr-reporting applies the approved UTRAnnotator 116 reporting profile
// (data/vep/UTRAnnotator116_all_effects.pm). The plugins' receipt names the
// patch files; record rewrites it to where storeDir keeps them.
process VEP_CACHE {
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
    def configured = adopted ? '' : installSetting(planned, 'assembly')
    def assembly = configured in ['hg19', 'GRCh37'] ? 'GRCh37' : (configured in ['hg38', 'GRCh38'] ? 'GRCh38' : '')
    def release = adopted ? '' : installSetting(planned, 'vep_installed_version')
    def plugins = 'SpliceAI,SpliceVault,UTRAnnotator,AlphaMissense,GO,Phenotypes,Conservation,LOEUF,SameCodon,NMD,PrimateAI,MaveDB'
    def cache = "data/VEP_cache/homo_sapiens_merged/${release}_${assembly}".toString()
    if (adopted)
        installAdopt(id, lineage)
    else if (!assembly)
        "echo 'VEP cache: select hg19/GRCh37 or hg38/GRCh38, not ${configured}' >&2; exit 2"
    else if (!(release ==~ /[1-9][0-9]*/))
        "echo 'VEP cache: vep_installed_version must be a numeric major release, not ${release}' >&2; exit 2"
    else
        installStart(task) + """
        here=\$(pwd -P)
        release=${installWord(release)}
        export PRIVA_DATA_DIR=${installWord("${projectDir}/data")}
        [[ -s ${installWord("${projectDir}/data/vep/UTRAnnotator116_all_effects.pm")} ]] || {
            echo 'VEP cache: the bundled UTRAnnotator 116 reporting module is missing' >&2; exit 1; }
        # The release of vep, of its Ensembl API and of its install helper.
        for tool in vep vep_install perl; do
            type -P \$tool > /dev/null || { echo "VEP cache: \$tool is not on PATH" >&2; exit 1; }
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
            echo "VEP cache: requested VEP \$release; executable=\$actual, API=\$api, installer=\$installer" >&2; exit 1; }
        # INSTALL.pl's DESTDIR: the Perl module folder of the environment.
        destdir=\$(perl -e 'print join("\\n", @INC);' | grep site_perl | head -1)
        [[ -n \$destdir ]] || { echo 'VEP cache: perl has no site_perl folder' >&2; exit 1; }
        mkdir -p data/VEP_cache data/VEP_plugin
        staging=\$(mktemp -d "\$here/VEP_plugin.install-\$release.XXXXXX")
        vep_install -d "\$destdir" --AUTO cp --NO_UPDATE --CACHE_VERSION "\$release" -s homo_sapiens_merged \\
            --NO_HTSLIB --NO_BIOPERL --CONVERT --CACHEDIR "\$here/data/VEP_cache" --ASSEMBLY ${assembly} \\
            --PLUGINS ${plugins} --PLUGINSDIR "\$staging" --PLUGINURL https://raw.githubusercontent.com/Ensembl/VEP_plugins
        # priva_validate_vep_cache
        for expected in info.txt chr_synonyms.txt 1/all_vars.gz; do
            [[ -s ${cache}/\$expected ]] || { echo "VEP cache: ${cache}/\$expected is missing" >&2; exit 1; }
        done
        [[ -s ${cache}/1/all_vars.gz.csi || -s ${cache}/1/all_vars.gz.tbi ]] || { echo 'VEP cache: the variant index is missing' >&2; exit 1; }
        awk -v assembly=${assembly} '
            \$1 == "assembly" { correct_assembly = (\$2 == assembly) }
            \$1 == "species" { correct_species = (\$2 == "homo_sapiens") }
            \$1 == "var_type" { converted = (\$2 == "tabix") }
            END { exit !(correct_assembly && correct_species && converted) }
        ' ${cache}/info.txt || { echo 'VEP cache: info.txt names another assembly or species, or variants not converted to tabix' >&2; exit 1; }
        "\$PRIVA_BIN" vep-plugins install --directory "\$here/data/VEP_plugin" --release "\$release" \\
            --plugins ${plugins} --stage "\$staging"
        rm -rf "\$staging"
        "\$PRIVA_BIN" vep-plugins ensure-utr-reporting --directory "\$here/data/VEP_plugin" --release "\$release"
        # The plugin receipt names its patch files where storeDir keeps them.
        ${installRecord(id, lineage, receipts, [vep_cache_dir: 'data/VEP_cache', vep_plugins_dir: 'data/VEP_plugin'])} --relocate data/VEP_plugin/.priva-vep-plugins.json
        """.stripIndent()
    stub:
    installStub(id, lineage, receipts, task)
}
