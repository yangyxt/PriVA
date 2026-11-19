// The install workflow: one task per install step of conf/install_steps.json,
// wired by each step's needs. Generated from the registry by
// tests/install_workflow_check.py (the check fails when the two differ);
// edit the registry, then regenerate this file.
//
// A step's lineage is a hash of its code digest (installer functions, scripts,
// data files and the Rust files of the priva commands it runs), its settings,
// the files you supply, the release of a changing source (or its pin), and the
// content of every step it needs. Its results live in
// resource_root/STEP/LINEAGE/, so a change reruns the step it touches, and
// the steps after it only when that step's results changed.

include { INSTALL_RELEASE; INSTALL_PLAN; INSTALL_PREFLIGHT; INSTALL_PUBLISH } from '../modules/local/install'
include { @INCLUDES@ } from '../modules/local/install'

def installSha(String text) {
    java.security.MessageDigest.getInstance('SHA-256').digest(text.getBytes('UTF-8')).encodeHex().toString()
}

def installLineage(String id, plan, List receipts) {
    def planned = new groovy.json.JsonSlurper().parse(plan.toFile())
    def step = planned.steps[id]
    if (step.adopt != null) return 'adopted-' + installSha(groovy.json.JsonOutput.toJson([id: id, identity: step.adopt.identity])).take(16)
    // An input counts by its content: what its results hold, not how they were made.
    def inputs = receipts.collect { r -> def x = new groovy.json.JsonSlurper().parse(new File(r.toString())); [x.id, x.content ?: x.lineage] }.sort { a, b -> a[0] <=> b[0] }
    def payload = [id: id, digest: step.digest, settings: new TreeMap(step.settings ?: [:]),
                   user_files: step.user_files, inputs: inputs]
    // Only steps with a changing source carry a release, so the others keep their lineage.
    if (step.release != null) payload.release = step.release
    installSha(groovy.json.JsonOutput.toJson(payload)).take(16)
}

// One item: [id, lineage, plan, [receipt paths]]. Every channel holds one item.
def installWire(String id, plan, List needs) {
    def joined = needs.inject(plan.map { [it] }) { acc, upstream -> acc.combine(upstream) }
    joined.map { items ->
        def p = items[0]
        def receipts = items.size() > 1 ? items[1..-1].collect { it.toString() } : []
        [id, installLineage(id, p, receipts), p, receipts]
    }
}

workflow INSTALL_RESOURCES {
    main:
    ['install_config', 'resource_root', 'priva_binary'].each { key ->
        if (!params[key]) error "Required parameter --${key} is missing"
    }
    def registryPath = file(params.install_registry ?: "${projectDir}/conf/install_steps.json")
    // Changing sources: one release check each; pinned ones only with --refresh,
    // or when the pin can expire (max_age_months; the plan decides whether it has).
    def registryEntries = new groovy.json.JsonSlurper().parse(registryPath.toFile()).steps
    def refreshSteps = (params.refresh ?: '').toString().tokenize(',').collect { s -> s.trim() }
    def checked = registryEntries.findAll { s -> s.release && (!s.pinned || s.max_age_months || refreshSteps.contains(s.id.toString())) }.collect { s -> s.id.toString() }
    INSTALL_RELEASE(Channel.fromList(checked).map { id -> [id, file(params.install_config)] })
    def releases = INSTALL_RELEASE.out.release.collect().ifEmpty([])
    INSTALL_PLAN(Channel.of([file(params.install_config), registryPath, params.adopt_config ?: '', params.adopt ?: '', params.run_steps ?: '', refreshSteps.join(',')])
        .combine(releases.map { files -> [files] }))
    // Every step waits for the preflight checks (files you supply, shipped data).
    INSTALL_PREFLIGHT(INSTALL_PLAN.out.plan)
    def plan = INSTALL_PREFLIGHT.out.plan

