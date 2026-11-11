nextflow.enable.dsl=2

include { VALIDATE_INPUTS; ANNOTATE_COHORT; ANNOTATE_PREPROCESS; ANNOTATE_STEP; ANNOTATE_CADD; ANNOTATE_MERGE; HUB_LOOKUP; HUB_FILL; HUB_ADD; FILTER_FAMILY; PRIORITIZE_COMBINE; PRIORITIZE_SPLICING; PRIORITIZE_UTR; PRIORITIZE_JOIN; PRIORITIZE_ACMG } from './modules/local/compatibility'
include { INSTALL_RESOURCES } from './workflows/install'

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
    // --action install builds every resource as its own task (workflows/install.nf):
    //   --install_config CONFIG.yaml --resource_root FOLDER --priva_binary PRIVA
    //   [--install_outdir FOLDER] [--adopt_config OLD.yaml --adopt all --run_steps ID,ID]
    if (params.action == 'install') {
        INSTALL_RESOURCES()
    } else if (params.help) {
        log.info '''PriVA workflow
Required: --input VCF.vcf.gz --assembly hg19|hg38
          --legacy_config CONFIG.yaml --resource_manifest RESOURCES.json
          --software_id PINNED_ENVIRONMENT_OR_IMAGE_ID --outdir NEW_OUTPUT_DIRECTORY
          --priva_binary PRIVA_PROGRAM (runs the task contract: priva task)
Install:  --action install --install_config CONFIG.yaml --resource_root FOLDER
          --priva_binary PRIVA_PROGRAM [--install_outdir FOLDER]
Optional: --ped COHORT.ped --family ID --mode raw|annotated
          --cadd ANNOTATED.cadd.tsv --vcf_index VCF.vcf.gz.tbi
          --vcf_hub_mode off|qualified-copy --cadd_hub_mode off|copy
          --annotation_hub FOLDER (raw mode with both old hub modes off)
The cohort VCF is validated and annotated once. With --ped, every PED family is
then filtered and prioritized in its own tasks (--family selects one); without
--ped, one run covers all samples without family filtering. The annotated mode
requires its matching CADD TSV. Preserve the original proband-first PED row
order. Use -stub-run only for wiring tests. See README.md for configuration,
outputs and scientific-profile requirements.'''
    } else {
        ['input', 'assembly', 'legacy_config', 'resource_manifest', 'software_id', 'outdir', 'priva_binary'].each { key ->
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
        if (params.annotation_hub && !(params.mode == 'raw' && params.vcf_hub_mode == 'off' && params.cadd_hub_mode == 'off')) error '--annotation_hub needs --mode raw with --vcf_hub_mode off and --cadd_hub_mode off'
        // A value starting with '-' after a space is read by Nextflow as a
        // separate flag, which leaves this parameter set to true (a boolean,
        // or the text 'true' depending on how the launch reads it).
        ['slurm_cluster_options', 'slurm_account', 'container_run_options'].each { name ->
            if (params[name] instanceof Boolean || params[name]?.toString() == 'true') error "--${name} needs a value; write it as --${name}='VALUE'"
        }
        if (params.family && !(params.family.toString() ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) error 'Unsafe family identifier'
        if (params.family && !params.ped) error '--family selects a family of the --ped file; pass --ped'
        if (!(params.input.toString().endsWith('.vcf.gz'))) error '--input must be a BGZF .vcf.gz file'
        def sourceRoot = new File(projectDir.toString()).canonicalFile
        if (!(sourceRoot.toString() ==~ /[A-Za-z0-9_\.\/:+\-]+/)) error 'Legacy shell scripts require a source path without whitespace or shell metacharacters'
        def privaBinary = new File(params.priva_binary.toString()).canonicalFile
        if (!privaBinary.isFile() || !privaBinary.canExecute()) error '--priva_binary must be an executable PriVA program'
        if (!(privaBinary.toString() ==~ /[A-Za-z0-9_\.\/:+\-]+/)) error '--priva_binary must be a path without whitespace or shell metacharacters'
        def base = new File(params.input.toString()).name.replaceFirst(/\.vcf\.gz$/, '')
        if (params.mode == 'annotated') base = base.replaceFirst(/\.anno$/, '')
        if (!(base ==~ /[A-Za-z0-9_][A-Za-z0-9_.-]*/)) error 'Unsafe input basename'
        def meta = [family: params.family ? params.family.toString() : '', base: base, assembly: params.assembly, mode: params.mode,
                    source_root: sourceRoot.toString(), source_digest: sourceDigest(sourceRoot),
                    software_id: params.software_id.toString(), vcf_hub_mode: params.vcf_hub_mode,
                    cadd_hub_mode: params.cadd_hub_mode, task_boundary: taskBoundary,
                    // The task program's bytes, not only its path, join every task hash.
                    priva_binary: privaBinary.toString(),
                    priva_sha256: java.security.MessageDigest.getInstance('SHA-256').digest(privaBinary.bytes).encodeHex().toString()]
        // This check executes outside the task cache. A fully cached -resume
        // still verifies external resources and the published run identity.
        def launchCommand = [meta.priva_binary, 'task', 'launch',
            '--vcf', file(params.input, checkIfExists: true).toString(),
            '--vcf-index', file(params.vcf_index ?: "${params.input}.tbi", checkIfExists: true).toString(),
            '--config', file(params.legacy_config, checkIfExists: true).toString(),
            '--resources', file(params.resource_manifest, checkIfExists: true).toString(),
            '--assembly', meta.assembly.toString(), '--mode', meta.mode.toString(),
            '--legacy-root', meta.source_root, '--source-digest', meta.source_digest,
            '--software-id', meta.software_id, '--vcf-hub-mode', meta.vcf_hub_mode.toString(),
            '--cadd-hub-mode', meta.cadd_hub_mode.toString(), '--outdir', params.outdir.toString(),
            '--task-boundary', meta.task_boundary.toString(),
            '--annotation-cpus', params.annotation_cpus.toString(), '--family-cpus', params.family_cpus.toString()]
        if (params.cadd) launchCommand.addAll(['--cadd', file(params.cadd, checkIfExists: true).toString()])
        if (params.ped) launchCommand.addAll(['--ped', file(params.ped, checkIfExists: true).toString()])
        if (meta.family) launchCommand.addAll(['--family', meta.family])
        if (workflow.stubRun) launchCommand.add('--stub')
        // The check's own message joins the error, so it stays in .nextflow.log.
        def preflight = new ProcessBuilder(launchCommand).redirectErrorStream(true).start()
        def report = preflight.inputStream.text.readLines().findAll { it.trim() }
        if (preflight.waitFor() != 0) error "Launch-time resource/output contract validation failed: ${report.takeRight(5).join(' | ')}"
        report.each { log.info it }
        def inputs = Channel.of(tuple(meta,
            file(params.input, checkIfExists: true),
            file(params.vcf_index ?: "${params.input}.tbi", checkIfExists: true),
            params.ped ? file(params.ped, checkIfExists: true) : [],
            file(params.legacy_config, checkIfExists: true),
            file(params.resource_manifest, checkIfExists: true),
            params.cadd ? file(params.cadd, checkIfExists: true) : []))
        // The bundle passes between tasks as its files, not as a folder: a
        // folder's own size and time on network file systems can differ
        // between reads, which would change the task hash and defeat -resume.
        // Sorting by name keeps the list order, also part of the hash, fixed.
        def bundleLast = { row -> row.take(row.size() - 1) + [(row.last() instanceof List ? row.last() : [row.last()]).sort(false) { it.name }] }
        VALIDATE_INPUTS(inputs)
        // The validated contract lists the families this run covers; none
        // (no PED) is one run without a family, written as ''.
        def families = VALIDATE_INPUTS.out.validated.flatMap { row ->
            def contract = row.last().find { it.name == 'contract.json' }
            def listed = new groovy.json.JsonSlurper().parse(contract.toFile()).families
            listed ? listed.collect { it.toString() } : ['']
        }
        // Raw mode with both hubs off runs the annotation script's functions as
        // their own tasks; otherwise (annotated mode, or a hub in use, whose
        // logic reshapes the annotation) the whole script is one task.
        def annotated = Channel.empty()
        if (params.mode == 'raw' && params.vcf_hub_mode == 'off' && params.cadd_hub_mode == 'off') {
            ANNOTATE_PREPROCESS(VALIDATE_INPUTS.out.validated.map(bundleLast).map { run, vcf, index, cadd, bundle -> [run, vcf, index, bundle] })
            def preprocessed = ANNOTATE_PREPROCESS.out.preprocessed.map(bundleLast)
            // Each annotation step's output, as [source, file]: gnomad, clinvar,
            // control, vep (VCFs with every preprocessed record) and cadd (table).
            def stepOutputs = Channel.empty()
            if (params.annotation_hub) {
                // Hub: each step's stored results are looked up first; the step
                // runs only on new variants; fill rebuilds its full output; add
                // stores the new results (hub folders named in the contract).
                def toLookup = preprocessed.flatMap { run, vcf, index, bundle ->
                    def contract = new groovy.json.JsonSlurper().parse(bundle.find { f -> f.name == 'contract.json' }.toFile())
                    ['gnomad', 'clinvar', 'control', 'vep', 'cadd'].collect { source -> [run, source, contract.hub_lineages[source].toString(), vcf, index] }
                }
                HUB_LOOKUP(toLookup)
                def lookedUp = HUB_LOOKUP.out.split.map { run, source, lineage, known, fresh, freshIndex, count -> [run, source, lineage, known, fresh, freshIndex, count.trim().toInteger()] }
                def withNew = lookedUp.filter { item -> item[6] > 0 }
                ANNOTATE_STEP(withNew.filter { item -> item[1] != 'cadd' }.combine(preprocessed.map { item -> [item[3]] })
                    .map { run, source, lineage, known, fresh, freshIndex, count, bundle -> [run, source, fresh, freshIndex, bundle] })
                ANNOTATE_CADD(withNew.filter { item -> item[1] == 'cadd' }.combine(preprocessed.map { item -> [item[3]] })
                    .map { run, source, lineage, known, fresh, freshIndex, count, bundle -> [run, fresh, freshIndex, bundle] })
                def annotatedNew = ANNOTATE_STEP.out.annotated.map { run, source, output -> [source, output] }
                    .mix(ANNOTATE_CADD.out.cadd.map { run, output -> ['cadd', output] })
                def noNew = file("${projectDir}/assets/hub_no_new_variants")
                def stepResults = annotatedNew.mix(lookedUp.filter { item -> item[6] == 0 }.map { item -> [item[1], noNew] })
                def toFill = lookedUp.map { run, source, lineage, known, fresh, freshIndex, count -> [source, run, known] }
                    .join(stepResults)
                    .combine(preprocessed.map { item -> [item[1], item[2]] })
                    .map { source, run, known, output, vcf, index -> [run, source, vcf, index, known, output] }
                HUB_FILL(toFill)
                def toAdd = withNew.filter { item -> item[2] != 'off' }
                    .map { run, source, lineage, known, fresh, freshIndex, count -> [source, run, lineage, fresh] }
                    .join(annotatedNew)
                    .map { source, run, lineage, fresh, output -> [run, source, lineage, fresh, output] }
                HUB_ADD(toAdd)
                stepOutputs = HUB_FILL.out.filled.map { run, source, output -> [source, output] }
            } else {
                // gnomAD, ClinVar, control cohort and VEP at the same time, and CADD.
                ANNOTATE_STEP(preprocessed.combine(Channel.of('gnomad', 'clinvar', 'control', 'vep'))
                    .map { run, vcf, index, bundle, step -> [run, step, vcf, index, bundle] })
                ANNOTATE_CADD(preprocessed)
                stepOutputs = ANNOTATE_STEP.out.annotated.map { run, source, output -> [source, output] }
                    .mix(ANNOTATE_CADD.out.cadd.map { run, output -> ['cadd', output] })
            }
            def byStep = stepOutputs.branch { item ->
                gnomad: item[0] == 'gnomad'
                clinvar: item[0] == 'clinvar'
                control: item[0] == 'control'
                vep: item[0] == 'vep'
                cadd: item[0] == 'cadd'
            }
            // One cohort per run, so combine joins exactly one item of each.
            def toMerge = preprocessed
                .combine(byStep.gnomad.map { it[1] })
                .combine(byStep.clinvar.map { it[1] })
                .combine(byStep.control.map { it[1] })
                .combine(byStep.vep.map { it[1] })
                .combine(byStep.cadd.map { it[1] })
                .map { run, vcf, index, bundle, gnomad, clinvar, control, vep, cadd -> [run, vcf, index, gnomad, clinvar, control, vep, cadd, bundle] }
            ANNOTATE_MERGE(toMerge)
            annotated = ANNOTATE_MERGE.out.annotated.map(bundleLast)
        } else {
            ANNOTATE_COHORT(VALIDATE_INPUTS.out.validated.map(bundleLast))
            annotated = ANNOTATE_COHORT.out.annotated.map(bundleLast)
        }
        // One filtration task per family: (meta, family, vcf, index, cadd, bundle).
        FILTER_FAMILY(annotated.combine(families).map { run, vcf, index, cadd, bundle, family -> [run, family, vcf, index, cadd, bundle] })
        // Families prioritize with their filtered CADD table; the no-pedigree
        // run uses the whole annotation CADD table, as the shell workflow did.
        def toPrioritize = FILTER_FAMILY.out.filtered.map(bundleLast)
            .combine(annotated.map { it[3] })
            .map { run, family, vcf, index, cadd, bundle, cohortCadd -> [run, family, vcf, index, family ? cadd : cohortCadd, bundle] }
        // Prioritization per family, one task per step: combine, then splicing
        // and 5' UTR side by side on copies of the table, join, then ACMG.
        PRIORITIZE_COMBINE(toPrioritize)
        def combined = PRIORITIZE_COMBINE.out.table.map(bundleLast)
        PRIORITIZE_SPLICING(combined)
        PRIORITIZE_UTR(combined)
        def toJoin = combined.map { run, family, table, bundle -> [family, run, table, bundle] }
            .join(PRIORITIZE_SPLICING.out.table.map { run, family, table -> [family, table] })
            .join(PRIORITIZE_UTR.out.table.map { run, family, table -> [family, table] })
            .map { family, run, table, bundle, splicing, utr -> [run, family, table, splicing, utr, bundle] }
        PRIORITIZE_JOIN(toJoin)
        PRIORITIZE_ACMG(PRIORITIZE_JOIN.out.table.map(bundleLast))
    }
}
