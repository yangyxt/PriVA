nextflow.enable.dsl=2

include { VALIDATE_FAMILY; ANNOTATE_FAMILY; FILTER_FAMILY; PRIORITIZE_FAMILY } from './modules/local/compatibility'

// The scientific source is deliberately retained during orchestration parity.
// Its content identity, rather than only its pathname, participates in every task hash.
def sourceDigest(root) {
    def digest = java.security.MessageDigest.getInstance('SHA-256')
    def files = []
    ['scripts', 'bin', 'modules', 'conf'].each { directory ->
        new File(root.toString(), directory).eachFileRecurse { candidate ->
            if (candidate.isFile() && !candidate.name.endsWith('.pyc') && !candidate.path.contains('__pycache__')) files << candidate
        }
    }
    ['main.nf', 'nextflow.config'].each { name -> files << new File(root.toString(), name) }
    files.sort { a, b -> a.path <=> b.path }.each { candidate ->
        digest.update(root.toPath().relativize(candidate.toPath()).toString().getBytes('UTF-8'))
        digest.update('\u0000'.getBytes('UTF-8'))
        digest.update(candidate.bytes)
    }
    digest.digest().encodeHex().toString()
}

workflow {
    if (params.help) {
        log.info '''PriVA single-family compatibility workflow
Required: --input VCF.vcf.gz --ped FAMILY.ped --family ID --assembly hg19|hg38
          --legacy_config CONFIG.yaml --resource_manifest RESOURCES.json
          --software_id PINNED_ENVIRONMENT_OR_IMAGE_ID --outdir NEW_OUTPUT_DIRECTORY
Optional: --mode raw|annotated --cadd ANNOTATED.cadd.tsv --vcf_index VCF.vcf.gz.tbi
          --vcf_hub_mode off|qualified-copy --cadd_hub_mode off|copy
The annotated mode requires its matching CADD TSV. The full cohort PED and VCF
are retained for annotation; only the selected family is prioritized. Preserve
the original proband-first PED row order. Use -stub-run only for wiring tests.
See README.md for configuration, outputs and scientific-profile requirements.'''
    } else {
        ['input', 'ped', 'family', 'assembly', 'legacy_config', 'resource_manifest', 'software_id', 'outdir'].each { key ->
            if (!params[key]) error "Required parameter --${key} is missing"
        }
        if (!(params.mode in ['raw', 'annotated'])) error '--mode must be raw or annotated'
        if (!(params.assembly in ['hg19', 'hg38', 'GRCh37', 'GRCh38'])) error 'Unsupported --assembly'
        if (!(params.vcf_hub_mode in ['off', 'qualified-copy'])) error 'Unsupported --vcf_hub_mode'
        if (!(params.cadd_hub_mode in ['off', 'copy'])) error 'Unsupported --cadd_hub_mode'
        def readonlyProfile = workflow.profile.tokenize(',').contains('task_readonly')
        // 24 October 2026 CT: derive identity from the selected profile itself.
        // A separate params default can override an included profile assignment
        // (observed with pinned 25.10.4), leaving its shell and identity at odds.
        def taskBoundary = readonlyProfile ? 'bwrap' : 'off'
        if (readonlyProfile && workflow.profile.tokenize(',').contains('apptainer')) error 'Combined task_readonly/apptainer execution is not qualified; choose one profile'
        if (workflow.profile.tokenize(',').contains('test') && !workflow.stubRun) error 'The test profile requires -stub-run and is never a scientific acceptance run'
        if (workflow.profile.tokenize(',').contains('apptainer') && !params.container_image) error 'The apptainer profile requires --container_image pinned to a tested image'
        if (params.mode == 'annotated' && !params.cadd) error '--mode annotated requires the matching --cadd TSV'
        if (params.mode == 'raw' && params.cadd) error '--cadd is only accepted with --mode annotated'
        if (!(params.family ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) error 'Unsafe family identifier'
        if (!(params.input.toString().endsWith('.vcf.gz'))) error '--input must be a BGZF .vcf.gz file'
        def sourceRoot = new File(projectDir.toString()).canonicalFile
        if (!(sourceRoot.toString() ==~ /[A-Za-z0-9_\.\/:+\-]+/)) error 'Legacy shell scripts require a source path without whitespace or shell metacharacters'
        def base = new File(params.input.toString()).name.replaceFirst(/\.vcf\.gz$/, '')
        if (params.mode == 'annotated') base = base.replaceFirst(/\.anno$/, '')
        if (!(base ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) error 'Unsafe input basename'
        def meta = [id: params.family, base: base, assembly: params.assembly, mode: params.mode,
                    source_root: sourceRoot.toString(), source_digest: sourceDigest(sourceRoot),
                    software_id: params.software_id.toString(), vcf_hub_mode: params.vcf_hub_mode,
                    cadd_hub_mode: params.cadd_hub_mode, task_boundary: taskBoundary]
        // This check executes outside the task cache. A fully cached -resume
        // still verifies external resources and the published run identity.
        def launchCommand = ['python3', "${sourceRoot}/bin/priva_task.py", 'launch',
            '--vcf', file(params.input, checkIfExists: true).toString(),
            '--vcf-index', file(params.vcf_index ?: "${params.input}.tbi", checkIfExists: true).toString(),
            '--ped', file(params.ped, checkIfExists: true).toString(),
            '--config', file(params.legacy_config, checkIfExists: true).toString(),
            '--resources', file(params.resource_manifest, checkIfExists: true).toString(),
            '--family', meta.id.toString(), '--assembly', meta.assembly.toString(), '--mode', meta.mode.toString(),
            '--legacy-root', meta.source_root, '--source-digest', meta.source_digest,
            '--software-id', meta.software_id, '--vcf-hub-mode', meta.vcf_hub_mode.toString(),
            '--cadd-hub-mode', meta.cadd_hub_mode.toString(), '--outdir', params.outdir.toString(),
            '--task-boundary', meta.task_boundary.toString(),
            '--annotation-cpus', params.annotation_cpus.toString(), '--family-cpus', params.family_cpus.toString()]
        if (params.cadd) launchCommand.addAll(['--cadd', file(params.cadd, checkIfExists: true).toString()])
        if (workflow.stubRun) launchCommand.add('--stub')
        def preflight = new ProcessBuilder(launchCommand).inheritIO().start()
        if (preflight.waitFor() != 0) error 'Launch-time resource/output contract validation failed'
        def inputs = Channel.of(tuple(meta,
            file(params.input, checkIfExists: true),
            file(params.vcf_index ?: "${params.input}.tbi", checkIfExists: true),
            file(params.ped, checkIfExists: true),
            file(params.legacy_config, checkIfExists: true),
            file(params.resource_manifest, checkIfExists: true),
            params.cadd ? file(params.cadd, checkIfExists: true) : []))
        VALIDATE_FAMILY(inputs)
        ANNOTATE_FAMILY(VALIDATE_FAMILY.out.validated)
        FILTER_FAMILY(ANNOTATE_FAMILY.out.annotated)
        PRIORITIZE_FAMILY(FILTER_FAMILY.out.filtered)
    }
}
