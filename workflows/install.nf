// The install workflow: one task per install step of conf/install_steps.json,
// wired by each step's needs. Generated from the registry by
// tests/install_workflow_check.py (the check fails when the two differ);
// edit the registry, then regenerate this file.
//
// A step's lineage is a hash of its code digest, its settings, the files you
// supply, the release of a changing source (or its pin), the priva source
// identity for native builders, and the lineage of every step it needs. Its
// results live in resource_root/STEP/LINEAGE/, so a change reruns exactly the
// steps it reaches.

include { INSTALL_RELEASE; INSTALL_PLAN; INSTALL_PUBLISH } from '../modules/local/install'
include { INSTALL_STEP as STEP_A1_REFERENCE_FASTA; INSTALL_STEP as STEP_A2_VEP_CACHE; INSTALL_STEP as STEP_A4_CDS_FASTA; INSTALL_STEP as STEP_A5_NM_TO_ENST; INSTALL_STEP as STEP_B1_UTR_ANNOTATOR; INSTALL_STEP as STEP_B2_CONSERVATION; INSTALL_STEP as STEP_B3_ALPHAMISSENSE_TSV; INSTALL_STEP as STEP_B4_SPLICEAI; INSTALL_STEP as STEP_B5_SPLICEVAULT; INSTALL_STEP as STEP_B6_MAVEDB; INSTALL_STEP as STEP_B7_PRIMATEAI; INSTALL_STEP as STEP_B8_LOFTEE; INSTALL_STEP as STEP_B9_PEXT; INSTALL_STEP as STEP_C1_GNOMAD; INSTALL_STEP as STEP_C3_ALPHAMISSENSE_VEP; INSTALL_STEP as STEP_C4_ALPHAMISSENSE_GNOMAD; INSTALL_STEP as STEP_C6_ALPHAMISSENSE_DOMAIN; INSTALL_STEP as STEP_C8_INTERPRO_MAPPING; INSTALL_STEP as STEP_D1_CLINVAR; INSTALL_STEP as STEP_C10_DAS_DOMAINS; INSTALL_STEP as STEP_C11_PM1_REGIONS; INSTALL_STEP as STEP_D6_CLINVAR_PATHO_AF; INSTALL_STEP as STEP_D8_CLINVAR_PM5; INSTALL_STEP as STEP_D9_CLINGEN_MAP; INSTALL_STEP as STEP_E1_CADD; INSTALL_STEP as STEP_F1_HPO_MONDO; INSTALL_STEP as STEP_F4_GENE_PATHOGENIC_MECHANISM; INSTALL_STEP as STEP_F5_GOFCARDS_EXACT_GOF; INSTALL_STEP as STEP_F6_MECHANISM_SOURCES; INSTALL_STEP as STEP_F7_GOFCARDS_CLINVAR_INJECTION; INSTALL_STEP as STEP_F8_GENE_NONLOF_MECHANISM; INSTALL_STEP as STEP_F9_HPO_CONDITION_MECHANISM } from '../modules/local/install'

def installSha(String text) {
    java.security.MessageDigest.getInstance('SHA-256').digest(text.getBytes('UTF-8')).encodeHex().toString()
}

def installLineage(String id, plan, List receipts, Set nativeSteps) {
    def planned = new groovy.json.JsonSlurper().parse(plan.toFile())
    def step = planned.steps[id]
    if (step.adopt != null) return 'adopted-' + installSha(groovy.json.JsonOutput.toJson([id: id, identity: step.adopt.identity])).take(16)
    def inputs = receipts.collect { r -> def x = new groovy.json.JsonSlurper().parse(new File(r.toString())); [x.id, x.lineage] }.sort { a, b -> a[0] <=> b[0] }
    def payload = [id: id, digest: step.digest, settings: new TreeMap(step.settings ?: [:]),
                   user_files: step.user_files, inputs: inputs,
                   priva: nativeSteps.contains(id) ? planned.priva_identity : null]
    // Only steps with a changing source carry a release, so the others keep their lineage.
    if (step.release != null) payload.release = step.release
    installSha(groovy.json.JsonOutput.toJson(payload)).take(16)
}

// One item: [id, lineage, plan, [receipt paths]]. Every channel holds one item.
def installWire(String id, plan, List needs, Set nativeSteps) {
    def joined = needs.inject(plan.map { [it] }) { acc, upstream -> acc.combine(upstream) }
    joined.map { items ->
        def p = items[0]
        def receipts = items.size() > 1 ? items[1..-1].collect { it.toString() } : []
        [id, installLineage(id, p, receipts, nativeSteps), p, receipts]
    }
}

workflow INSTALL_RESOURCES {
    main:
    ['install_config', 'resource_root', 'priva_binary'].each { key ->
        if (!params[key]) error "Required parameter --${key} is missing"
    }
    def registryPath = file(params.install_registry ?: "${projectDir}/conf/install_steps.json")
    // Read once here: the parsed registry is not shared with the parallel lineage calls.
    def nativeSteps = new groovy.json.JsonSlurper().parse(registryPath.toFile()).steps.findAll { s -> s.native }.collect { s -> s.id.toString() } as Set
    // Changing sources: one release check each; pinned ones only with --refresh.
    def registryEntries = new groovy.json.JsonSlurper().parse(registryPath.toFile()).steps
    def refreshSteps = (params.refresh ?: '').toString().tokenize(',').collect { s -> s.trim() }
    def checked = registryEntries.findAll { s -> s.release && (!s.pinned || refreshSteps.contains(s.id.toString())) }.collect { s -> s.id.toString() }
    INSTALL_RELEASE(Channel.fromList(checked).map { id -> [id, file(params.install_config)] })
    def releases = INSTALL_RELEASE.out.release.collect().ifEmpty([])
    INSTALL_PLAN(Channel.of([file(params.install_config), registryPath, params.adopt_config ?: '', params.adopt ?: '', params.run_steps ?: '', refreshSteps.join(',')])
        .combine(releases.map { files -> [files] }))
    def plan = INSTALL_PLAN.out.plan
    STEP_A1_REFERENCE_FASTA(installWire('A1', plan, [], nativeSteps))
    STEP_A2_VEP_CACHE(installWire('A2', plan, [], nativeSteps))
    STEP_A4_CDS_FASTA(installWire('A4', plan, [], nativeSteps))
    STEP_A5_NM_TO_ENST(installWire('A5', plan, [], nativeSteps))
    STEP_B1_UTR_ANNOTATOR(installWire('B1', plan, [], nativeSteps))
    STEP_B2_CONSERVATION(installWire('B2', plan, [STEP_A1_REFERENCE_FASTA.out.receipt], nativeSteps))
    STEP_B3_ALPHAMISSENSE_TSV(installWire('B3', plan, [STEP_A2_VEP_CACHE.out.receipt], nativeSteps))
    STEP_B4_SPLICEAI(installWire('B4', plan, [STEP_A2_VEP_CACHE.out.receipt], nativeSteps))
    STEP_B5_SPLICEVAULT(installWire('B5', plan, [STEP_A1_REFERENCE_FASTA.out.receipt], nativeSteps))
    STEP_B6_MAVEDB(installWire('B6', plan, [], nativeSteps))
    STEP_B7_PRIMATEAI(installWire('B7', plan, [STEP_A2_VEP_CACHE.out.receipt], nativeSteps))
    STEP_B8_LOFTEE(installWire('B8', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt], nativeSteps))
    STEP_B9_PEXT(installWire('B9', plan, [], nativeSteps))
    STEP_C1_GNOMAD(installWire('C1', plan, [STEP_A1_REFERENCE_FASTA.out.receipt], nativeSteps))
    STEP_C3_ALPHAMISSENSE_VEP(installWire('C3', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt], nativeSteps))
    STEP_C4_ALPHAMISSENSE_GNOMAD(installWire('C4', plan, [STEP_C3_ALPHAMISSENSE_VEP.out.receipt, STEP_C1_GNOMAD.out.receipt], nativeSteps))
    STEP_C6_ALPHAMISSENSE_DOMAIN(installWire('C6', plan, [STEP_C4_ALPHAMISSENSE_GNOMAD.out.receipt], nativeSteps))
    STEP_C8_INTERPRO_MAPPING(installWire('C8', plan, [STEP_C4_ALPHAMISSENSE_GNOMAD.out.receipt], nativeSteps))
    STEP_D1_CLINVAR(installWire('D1', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_C1_GNOMAD.out.receipt], nativeSteps))
    STEP_C10_DAS_DOMAINS(installWire('C10', plan, [STEP_C6_ALPHAMISSENSE_DOMAIN.out.receipt, STEP_C8_INTERPRO_MAPPING.out.receipt, STEP_C4_ALPHAMISSENSE_GNOMAD.out.receipt, STEP_D1_CLINVAR.out.receipt], nativeSteps))
    STEP_C11_PM1_REGIONS(installWire('C11', plan, [STEP_C10_DAS_DOMAINS.out.receipt, STEP_C6_ALPHAMISSENSE_DOMAIN.out.receipt, STEP_A5_NM_TO_ENST.out.receipt], nativeSteps))
    STEP_D6_CLINVAR_PATHO_AF(installWire('D6', plan, [STEP_D1_CLINVAR.out.receipt], nativeSteps))
    STEP_D8_CLINVAR_PM5(installWire('D8', plan, [STEP_D1_CLINVAR.out.receipt], nativeSteps))
    STEP_D9_CLINGEN_MAP(installWire('D9', plan, [STEP_D1_CLINVAR.out.receipt], nativeSteps))
    STEP_E1_CADD(installWire('E1', plan, [], nativeSteps))
    STEP_F1_HPO_MONDO(installWire('F1', plan, [], nativeSteps))
    STEP_F4_GENE_PATHOGENIC_MECHANISM(installWire('F4', plan, [STEP_F1_HPO_MONDO.out.receipt], nativeSteps))
    STEP_F5_GOFCARDS_EXACT_GOF(installWire('F5', plan, [], nativeSteps))
    STEP_F6_MECHANISM_SOURCES(installWire('F6', plan, [], nativeSteps))
    STEP_F7_GOFCARDS_CLINVAR_INJECTION(installWire('F7', plan, [STEP_F5_GOFCARDS_EXACT_GOF.out.receipt, STEP_F6_MECHANISM_SOURCES.out.receipt], nativeSteps))
    STEP_F8_GENE_NONLOF_MECHANISM(installWire('F8', plan, [STEP_F6_MECHANISM_SOURCES.out.receipt, STEP_F7_GOFCARDS_CLINVAR_INJECTION.out.receipt], nativeSteps))
    STEP_F9_HPO_CONDITION_MECHANISM(installWire('F9', plan, [STEP_F1_HPO_MONDO.out.receipt, STEP_F4_GENE_PATHOGENIC_MECHANISM.out.receipt, STEP_F7_GOFCARDS_CLINVAR_INJECTION.out.receipt, STEP_D9_CLINGEN_MAP.out.receipt], nativeSteps))
    def all = Channel.empty().mix(STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_A4_CDS_FASTA.out.receipt, STEP_A5_NM_TO_ENST.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_C1_GNOMAD.out.receipt, STEP_C3_ALPHAMISSENSE_VEP.out.receipt, STEP_C4_ALPHAMISSENSE_GNOMAD.out.receipt, STEP_C6_ALPHAMISSENSE_DOMAIN.out.receipt, STEP_C8_INTERPRO_MAPPING.out.receipt, STEP_D1_CLINVAR.out.receipt, STEP_C10_DAS_DOMAINS.out.receipt, STEP_C11_PM1_REGIONS.out.receipt, STEP_D6_CLINVAR_PATHO_AF.out.receipt, STEP_D8_CLINVAR_PM5.out.receipt, STEP_D9_CLINGEN_MAP.out.receipt, STEP_E1_CADD.out.receipt, STEP_F1_HPO_MONDO.out.receipt, STEP_F4_GENE_PATHOGENIC_MECHANISM.out.receipt, STEP_F5_GOFCARDS_EXACT_GOF.out.receipt, STEP_F6_MECHANISM_SOURCES.out.receipt, STEP_F7_GOFCARDS_CLINVAR_INJECTION.out.receipt, STEP_F8_GENE_NONLOF_MECHANISM.out.receipt, STEP_F9_HPO_CONDITION_MECHANISM.out.receipt).collect()
    def toPublish = plan.combine(all.map { [it] }).map { p, receipts ->
        def names = receipts.collect { it.toString() }.sort()
        def lineages = names.collect { r -> new groovy.json.JsonSlurper().parse(new File(r)).lineage }
        [installSha(groovy.json.JsonOutput.toJson(lineages)).take(16), p, names]
    }
    INSTALL_PUBLISH(toPublish)
    emit:
    published = INSTALL_PUBLISH.out.published
}
