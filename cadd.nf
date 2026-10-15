nextflow.enable.dsl=2

include { CADD_UPSTREAM_ANNOTATE; CADD_UPSTREAM_CALCULATE } from './modules/local/cadd_upstream'

def quoteCadd(value) { "'" + value.toString().replace("'", "'\\''") + "'" }

workflow {
    if (params.help) {
        log.info '''PriVA standalone CADD stages (all calculation stages retain upstream Python)
nextflow run cadd.nf -c conf/cadd.config --cadd_mode full|annotate|calculate
  --cadd_profile /absolute/configured-profile.json --assembly GRCh37|GRCh38
  --input INPUT.vcf.gz --outdir NEW_OUTPUT_DIRECTORY -work-dir DURABLE_WORK
For calculate mode, replace --input with --annotation_bundle COMPLETE_BUNDLE_DIR.
Annotation bundles contain annotations.tsv.gz and a matching manifest.json.
Imported complete tables: python3 bin/cadd_task.py import-annotations --help
Profiles require installed resources/runtimes; analysis never installs packages.
Reference templates: conf/cadd/{GRCh37,GRCh38}.reference.json.
-stub-run requires a profile explicitly marked stub=true; it is not acceptance.'''
    } else {
        ['cadd_mode', 'cadd_profile', 'assembly', 'outdir'].each { key ->
            if (!params[key]) error "Required parameter --${key} is missing; include -c conf/cadd.config"
        }
        if (!(params.cadd_mode in ['annotate', 'calculate', 'full'])) error '--cadd_mode must be annotate, calculate or full'
        if (!(params.assembly in ['hg19', 'hg38', 'GRCh37', 'GRCh38'])) error 'Unsupported --assembly'
        if (workflow.profile.tokenize(',').contains('test') && !workflow.stubRun) error 'The test profile requires -stub-run'
        if (workflow.profile.tokenize(',').contains('apptainer') && !params.container_image) error 'The apptainer profile requires a tested --container_image'
        if (params.cadd_mode == 'calculate' && (!params.annotation_bundle || params.input)) error 'Calculate requires --annotation_bundle and no --input'
        if (params.cadd_mode != 'calculate' && (!params.input || params.annotation_bundle)) error 'Annotate/full requires --input and no --annotation_bundle'
        def root = projectDir.toString()
        def profile = file(params.cadd_profile, checkIfExists: true)
        def input = file(params.cadd_mode == 'calculate' ? params.annotation_bundle : params.input, checkIfExists: true)
        def command = ['python3', "${root}/bin/cadd_task.py", 'preflight',
            '--profile', profile.toString(), '--assembly', params.assembly.toString(),
            '--input', input.toString(), '--mode', params.cadd_mode.toString(),
            '--output', params.outdir.toString()]
        if (workflow.stubRun) command.add('--stub')
        // Runs before Nextflow's task cache; fully cached -resume still checks
        // actual model/config/resource identities and the published run identity.
        def preflight = new ProcessBuilder(command).redirectError(ProcessBuilder.Redirect.INHERIT).start()
        def identity = preflight.inputStream.text.trim()
        if (preflight.waitFor() != 0) error 'CADD launch-time contract validation failed'
        def meta = [id: 'cadd', assembly: params.assembly, source_root: root,
                    profile_identity: identity]
        def inputs = Channel.of(tuple(meta, input, profile))
        if (params.cadd_mode == 'calculate') {
            CADD_UPSTREAM_CALCULATE(inputs)
        } else {
            CADD_UPSTREAM_ANNOTATE(inputs)
            if (params.cadd_mode == 'full') CADD_UPSTREAM_CALCULATE(CADD_UPSTREAM_ANNOTATE.out.annotations)
        }
    }
}
