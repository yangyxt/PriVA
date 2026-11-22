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
include { INSTALL_STEP as STEP_A2_VEP_CACHE; INSTALL_STEP as STEP_B4_SPLICEAI; INSTALL_STEP as STEP_C3_ALPHAMISSENSE_VEP; INSTALL_STEP as STEP_C4_ALPHAMISSENSE_GNOMAD; INSTALL_STEP as STEP_D1_split_CLINVAR_VEP_SPLIT; INSTALL_STEP as STEP_D1_chr1_CLINVAR_VEP_CHR1; INSTALL_STEP as STEP_D1_chr2_CLINVAR_VEP_CHR2; INSTALL_STEP as STEP_D1_chr3_CLINVAR_VEP_CHR3; INSTALL_STEP as STEP_D1_chr4_CLINVAR_VEP_CHR4; INSTALL_STEP as STEP_D1_chr5_CLINVAR_VEP_CHR5; INSTALL_STEP as STEP_D1_chr6_CLINVAR_VEP_CHR6; INSTALL_STEP as STEP_D1_chr7_CLINVAR_VEP_CHR7; INSTALL_STEP as STEP_D1_chr8_CLINVAR_VEP_CHR8; INSTALL_STEP as STEP_D1_chr9_CLINVAR_VEP_CHR9; INSTALL_STEP as STEP_D1_chr10_CLINVAR_VEP_CHR10; INSTALL_STEP as STEP_D1_chr11_CLINVAR_VEP_CHR11; INSTALL_STEP as STEP_D1_chr12_CLINVAR_VEP_CHR12; INSTALL_STEP as STEP_D1_chr13_CLINVAR_VEP_CHR13; INSTALL_STEP as STEP_D1_chr14_CLINVAR_VEP_CHR14; INSTALL_STEP as STEP_D1_chr15_CLINVAR_VEP_CHR15; INSTALL_STEP as STEP_D1_chr16_CLINVAR_VEP_CHR16; INSTALL_STEP as STEP_D1_chr17_CLINVAR_VEP_CHR17; INSTALL_STEP as STEP_D1_chr18_CLINVAR_VEP_CHR18; INSTALL_STEP as STEP_D1_chr19_CLINVAR_VEP_CHR19; INSTALL_STEP as STEP_D1_chr20_CLINVAR_VEP_CHR20; INSTALL_STEP as STEP_D1_chr21_CLINVAR_VEP_CHR21; INSTALL_STEP as STEP_D1_chr22_CLINVAR_VEP_CHR22; INSTALL_STEP as STEP_D1_chrX_CLINVAR_VEP_CHRX; INSTALL_STEP as STEP_D1_chrY_CLINVAR_VEP_CHRY; INSTALL_STEP as STEP_D1_chrOther_CLINVAR_VEP_CHROTHER; INSTALL_STEP as STEP_D1_CLINVAR; INSTALL_STEP as STEP_E1_CADD } from '../modules/local/install'
include { REFERENCE_FASTA as STEP_A1_REFERENCE_FASTA } from '../modules/local/install_steps/reference_fasta'
include { CDS_FASTA as STEP_A4_CDS_FASTA } from '../modules/local/install_steps/cds_fasta'
include { NM_TO_ENST as STEP_A5_NM_TO_ENST } from '../modules/local/install_steps/nm_to_enst'
include { UTR_ANNOTATOR as STEP_B1_UTR_ANNOTATOR } from '../modules/local/install_steps/utr_annotator'
include { CONSERVATION as STEP_B2_CONSERVATION } from '../modules/local/install_steps/conservation'
include { ALPHAMISSENSE_TSV as STEP_B3_ALPHAMISSENSE_TSV } from '../modules/local/install_steps/alphamissense'
include { SPLICEVAULT as STEP_B5_SPLICEVAULT } from '../modules/local/install_steps/splicevault'
include { MAVEDB as STEP_B6_MAVEDB } from '../modules/local/install_steps/mavedb'
include { PRIMATEAI as STEP_B7_PRIMATEAI } from '../modules/local/install_steps/primateai'
include { LOFTEE as STEP_B8_LOFTEE } from '../modules/local/install_steps/loftee'
include { PEXT as STEP_B9_PEXT } from '../modules/local/install_steps/pext'
include { GNOMAD_DOWNLOAD as STEP_C1_chr1_GNOMAD_DOWNLOAD_CHR1 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr2_GNOMAD_DOWNLOAD_CHR2 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr3_GNOMAD_DOWNLOAD_CHR3 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr4_GNOMAD_DOWNLOAD_CHR4 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr5_GNOMAD_DOWNLOAD_CHR5 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr6_GNOMAD_DOWNLOAD_CHR6 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr7_GNOMAD_DOWNLOAD_CHR7 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr8_GNOMAD_DOWNLOAD_CHR8 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr9_GNOMAD_DOWNLOAD_CHR9 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr10_GNOMAD_DOWNLOAD_CHR10 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr11_GNOMAD_DOWNLOAD_CHR11 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr12_GNOMAD_DOWNLOAD_CHR12 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr13_GNOMAD_DOWNLOAD_CHR13 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr14_GNOMAD_DOWNLOAD_CHR14 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr15_GNOMAD_DOWNLOAD_CHR15 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr16_GNOMAD_DOWNLOAD_CHR16 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr17_GNOMAD_DOWNLOAD_CHR17 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr18_GNOMAD_DOWNLOAD_CHR18 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr19_GNOMAD_DOWNLOAD_CHR19 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr20_GNOMAD_DOWNLOAD_CHR20 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr21_GNOMAD_DOWNLOAD_CHR21 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chr22_GNOMAD_DOWNLOAD_CHR22 } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chrX_GNOMAD_DOWNLOAD_CHRX } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD_DOWNLOAD as STEP_C1_chrY_GNOMAD_DOWNLOAD_CHRY } from '../modules/local/install_steps/gnomad_download'
include { GNOMAD as STEP_C1_GNOMAD } from '../modules/local/install_steps/gnomad'
include { ALPHAMISSENSE_DOMAIN as STEP_C6_ALPHAMISSENSE_DOMAIN } from '../modules/local/install_steps/alphamissense_domain'
include { INTERPRO_MAPPING as STEP_C8_INTERPRO_MAPPING } from '../modules/local/install_steps/interpro_mapping'
include { CLINVAR_DOWNLOAD as STEP_D0_CLINVAR_DOWNLOAD } from '../modules/local/install_steps/clinvar_download'
include { DAS_DOMAINS as STEP_C10_DAS_DOMAINS } from '../modules/local/install_steps/das_domains'
include { PM1_REGIONS as STEP_C11_PM1_REGIONS } from '../modules/local/install_steps/pm1_regions'
include { CLINVAR_PATHO_AF as STEP_D6_CLINVAR_PATHO_AF } from '../modules/local/install_steps/clinvar_patho_af'
include { CLINVAR_PM5 as STEP_D8_CLINVAR_PM5 } from '../modules/local/install_steps/clinvar_pm5'
include { CLINGEN_MAP as STEP_D9_CLINGEN_MAP } from '../modules/local/install_steps/clingen_map'
include { HPO_MONDO as STEP_F1_HPO_MONDO } from '../modules/local/install_steps/hpo_mondo'
include { GENE_PATHOGENIC_MECHANISM as STEP_F4_GENE_PATHOGENIC_MECHANISM } from '../modules/local/install_steps/gene_mechanism'
include { GOFCARDS_EXACT_GOF as STEP_F5_GOFCARDS_EXACT_GOF } from '../modules/local/install_steps/gofcards_exact_gof'
include { MECHANISM_SOURCES as STEP_F6_MECHANISM_SOURCES } from '../modules/local/install_steps/mechanism_sources'
include { GOFCARDS_CLINVAR_INJECTION as STEP_F7_GOFCARDS_CLINVAR_INJECTION } from '../modules/local/install_steps/gofcards_clinvar_injection'
include { GENE_NONLOF_MECHANISM as STEP_F8_GENE_NONLOF_MECHANISM } from '../modules/local/install_steps/gene_nonlof_mechanism'
include { HPO_CONDITION_MECHANISM as STEP_F9_HPO_CONDITION_MECHANISM } from '../modules/local/install_steps/hpo_condition_mechanism'
include { HPO_GENE_RANKING as STEP_F10_HPO_GENE_RANKING } from '../modules/local/install_steps/hpo_gene_ranking'

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

// A step that runs as its own process (registry "module",
// modules/local/install_steps/): one item [id, lineage, plan, [receipts],
// its plan entry, {key: file name in staged/}, [files to stage]]. It reads the
// keys its registry entry lists in "stage" from the receipts of the steps it
// needs (a file shipped in the checkout from this checkout, as install-task
// run reads it), with their index files and provenance (.provenance.json),
// and the files you supply for it (a file that does not exist is left out,
// and the step says it is missing).
def installNative(String id, plan, List needs, List stage) {
    installWire(id, plan, needs).map { item ->
        def planned = new groovy.json.JsonSlurper().parse(item[2].toFile()).steps[id]
        def wanted = []
        item[3].each { r ->
            def receipt = new groovy.json.JsonSlurper().parse(new File(r.toString()))
            stage.each { key ->
                def value = receipt.keys ? receipt.keys[key] : null
                if (value != null) {
                    def shipped = receipt.shipped ? receipt.shipped[key] : null
                    wanted << [key, shipped ? "${projectDir}/${shipped}".toString() : value.toString()]
                }
            }
        }
        (planned.user_paths ?: [:]).each { key, value -> if (value) wanted << [key.toString(), value.toString()] }
        def names = [:]
        def files = []
        wanted.each { pair ->
            def main = file(pair[1])
            if (main.exists()) {
                names[pair[0]] = main.name
                files << main
                ['.fai', '.gzi', '.tbi', '.csi', '.provenance.json'].each { ext ->
                    def index = file(pair[1] + ext)
                    if (index.exists()) files << index
                }
            }
        }
        [id, item[1], item[2], item[3].collect { r -> file(r) }, planned, names, files]
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
    STEP_A1_REFERENCE_FASTA(installNative('A1', plan, [], []))
    STEP_A2_VEP_CACHE(installWire('A2', plan, []))
    STEP_A4_CDS_FASTA(installNative('A4', plan, [], []))
    STEP_A5_NM_TO_ENST(installNative('A5', plan, [], []))
    STEP_B1_UTR_ANNOTATOR(installNative('B1', plan, [], []))
    STEP_B2_CONSERVATION(installNative('B2', plan, [STEP_A1_REFERENCE_FASTA.out.receipt], ['ref_genome']))
    STEP_B3_ALPHAMISSENSE_TSV(installNative('B3', plan, [STEP_A2_VEP_CACHE.out.receipt], []))
    STEP_B4_SPLICEAI(installWire('B4', plan, [STEP_A2_VEP_CACHE.out.receipt]))
    STEP_B5_SPLICEVAULT(installNative('B5', plan, [STEP_A1_REFERENCE_FASTA.out.receipt], ['ref_genome']))
    STEP_B6_MAVEDB(installNative('B6', plan, [], []))
    STEP_B7_PRIMATEAI(installNative('B7', plan, [STEP_A2_VEP_CACHE.out.receipt], []))
    STEP_B8_LOFTEE(installNative('B8', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt], ['ref_genome']))
    STEP_B9_PEXT(installNative('B9', plan, [], []))
    STEP_C1_chr1_GNOMAD_DOWNLOAD_CHR1(installNative('C1_chr1', plan, [], []))
    STEP_C1_chr2_GNOMAD_DOWNLOAD_CHR2(installNative('C1_chr2', plan, [], []))
    STEP_C1_chr3_GNOMAD_DOWNLOAD_CHR3(installNative('C1_chr3', plan, [], []))
    STEP_C1_chr4_GNOMAD_DOWNLOAD_CHR4(installNative('C1_chr4', plan, [], []))
    STEP_C1_chr5_GNOMAD_DOWNLOAD_CHR5(installNative('C1_chr5', plan, [], []))
    STEP_C1_chr6_GNOMAD_DOWNLOAD_CHR6(installNative('C1_chr6', plan, [], []))
    STEP_C1_chr7_GNOMAD_DOWNLOAD_CHR7(installNative('C1_chr7', plan, [], []))
    STEP_C1_chr8_GNOMAD_DOWNLOAD_CHR8(installNative('C1_chr8', plan, [], []))
    STEP_C1_chr9_GNOMAD_DOWNLOAD_CHR9(installNative('C1_chr9', plan, [], []))
    STEP_C1_chr10_GNOMAD_DOWNLOAD_CHR10(installNative('C1_chr10', plan, [], []))
    STEP_C1_chr11_GNOMAD_DOWNLOAD_CHR11(installNative('C1_chr11', plan, [], []))
    STEP_C1_chr12_GNOMAD_DOWNLOAD_CHR12(installNative('C1_chr12', plan, [], []))
    STEP_C1_chr13_GNOMAD_DOWNLOAD_CHR13(installNative('C1_chr13', plan, [], []))
    STEP_C1_chr14_GNOMAD_DOWNLOAD_CHR14(installNative('C1_chr14', plan, [], []))
    STEP_C1_chr15_GNOMAD_DOWNLOAD_CHR15(installNative('C1_chr15', plan, [], []))
    STEP_C1_chr16_GNOMAD_DOWNLOAD_CHR16(installNative('C1_chr16', plan, [], []))
    STEP_C1_chr17_GNOMAD_DOWNLOAD_CHR17(installNative('C1_chr17', plan, [], []))
    STEP_C1_chr18_GNOMAD_DOWNLOAD_CHR18(installNative('C1_chr18', plan, [], []))
    STEP_C1_chr19_GNOMAD_DOWNLOAD_CHR19(installNative('C1_chr19', plan, [], []))
    STEP_C1_chr20_GNOMAD_DOWNLOAD_CHR20(installNative('C1_chr20', plan, [], []))
    STEP_C1_chr21_GNOMAD_DOWNLOAD_CHR21(installNative('C1_chr21', plan, [], []))
    STEP_C1_chr22_GNOMAD_DOWNLOAD_CHR22(installNative('C1_chr22', plan, [], []))
    STEP_C1_chrX_GNOMAD_DOWNLOAD_CHRX(installNative('C1_chrX', plan, [], []))
    STEP_C1_chrY_GNOMAD_DOWNLOAD_CHRY(installNative('C1_chrY', plan, [], []))
    STEP_C1_GNOMAD(installNative('C1', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_C1_chr1_GNOMAD_DOWNLOAD_CHR1.out.receipt, STEP_C1_chr2_GNOMAD_DOWNLOAD_CHR2.out.receipt, STEP_C1_chr3_GNOMAD_DOWNLOAD_CHR3.out.receipt, STEP_C1_chr4_GNOMAD_DOWNLOAD_CHR4.out.receipt, STEP_C1_chr5_GNOMAD_DOWNLOAD_CHR5.out.receipt, STEP_C1_chr6_GNOMAD_DOWNLOAD_CHR6.out.receipt, STEP_C1_chr7_GNOMAD_DOWNLOAD_CHR7.out.receipt, STEP_C1_chr8_GNOMAD_DOWNLOAD_CHR8.out.receipt, STEP_C1_chr9_GNOMAD_DOWNLOAD_CHR9.out.receipt, STEP_C1_chr10_GNOMAD_DOWNLOAD_CHR10.out.receipt, STEP_C1_chr11_GNOMAD_DOWNLOAD_CHR11.out.receipt, STEP_C1_chr12_GNOMAD_DOWNLOAD_CHR12.out.receipt, STEP_C1_chr13_GNOMAD_DOWNLOAD_CHR13.out.receipt, STEP_C1_chr14_GNOMAD_DOWNLOAD_CHR14.out.receipt, STEP_C1_chr15_GNOMAD_DOWNLOAD_CHR15.out.receipt, STEP_C1_chr16_GNOMAD_DOWNLOAD_CHR16.out.receipt, STEP_C1_chr17_GNOMAD_DOWNLOAD_CHR17.out.receipt, STEP_C1_chr18_GNOMAD_DOWNLOAD_CHR18.out.receipt, STEP_C1_chr19_GNOMAD_DOWNLOAD_CHR19.out.receipt, STEP_C1_chr20_GNOMAD_DOWNLOAD_CHR20.out.receipt, STEP_C1_chr21_GNOMAD_DOWNLOAD_CHR21.out.receipt, STEP_C1_chr22_GNOMAD_DOWNLOAD_CHR22.out.receipt, STEP_C1_chrX_GNOMAD_DOWNLOAD_CHRX.out.receipt, STEP_C1_chrY_GNOMAD_DOWNLOAD_CHRY.out.receipt], ['ref_genome', 'gnomad_download_chr1', 'gnomad_download_chr2', 'gnomad_download_chr3', 'gnomad_download_chr4', 'gnomad_download_chr5', 'gnomad_download_chr6', 'gnomad_download_chr7', 'gnomad_download_chr8', 'gnomad_download_chr9', 'gnomad_download_chr10', 'gnomad_download_chr11', 'gnomad_download_chr12', 'gnomad_download_chr13', 'gnomad_download_chr14', 'gnomad_download_chr15', 'gnomad_download_chr16', 'gnomad_download_chr17', 'gnomad_download_chr18', 'gnomad_download_chr19', 'gnomad_download_chr20', 'gnomad_download_chr21', 'gnomad_download_chr22', 'gnomad_download_chrX', 'gnomad_download_chrY']))
    STEP_C3_ALPHAMISSENSE_VEP(installWire('C3', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt]))
    STEP_C4_ALPHAMISSENSE_GNOMAD(installWire('C4', plan, [STEP_C3_ALPHAMISSENSE_VEP.out.receipt, STEP_C1_GNOMAD.out.receipt]))
    STEP_C6_ALPHAMISSENSE_DOMAIN(installNative('C6', plan, [STEP_C4_ALPHAMISSENSE_GNOMAD.out.receipt], ['alphamissense_vep_vcf']))
    STEP_C8_INTERPRO_MAPPING(installNative('C8', plan, [STEP_C4_ALPHAMISSENSE_GNOMAD.out.receipt], ['alphamissense_vep_vcf']))
    STEP_D0_CLINVAR_DOWNLOAD(installNative('D0', plan, [], []))
    STEP_D1_split_CLINVAR_VEP_SPLIT(installWire('D1_split', plan, [STEP_D0_CLINVAR_DOWNLOAD.out.receipt]))
    STEP_D1_chr1_CLINVAR_VEP_CHR1(installWire('D1_chr1', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr2_CLINVAR_VEP_CHR2(installWire('D1_chr2', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr3_CLINVAR_VEP_CHR3(installWire('D1_chr3', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr4_CLINVAR_VEP_CHR4(installWire('D1_chr4', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr5_CLINVAR_VEP_CHR5(installWire('D1_chr5', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr6_CLINVAR_VEP_CHR6(installWire('D1_chr6', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr7_CLINVAR_VEP_CHR7(installWire('D1_chr7', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr8_CLINVAR_VEP_CHR8(installWire('D1_chr8', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr9_CLINVAR_VEP_CHR9(installWire('D1_chr9', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr10_CLINVAR_VEP_CHR10(installWire('D1_chr10', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr11_CLINVAR_VEP_CHR11(installWire('D1_chr11', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr12_CLINVAR_VEP_CHR12(installWire('D1_chr12', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr13_CLINVAR_VEP_CHR13(installWire('D1_chr13', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr14_CLINVAR_VEP_CHR14(installWire('D1_chr14', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr15_CLINVAR_VEP_CHR15(installWire('D1_chr15', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr16_CLINVAR_VEP_CHR16(installWire('D1_chr16', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr17_CLINVAR_VEP_CHR17(installWire('D1_chr17', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr18_CLINVAR_VEP_CHR18(installWire('D1_chr18', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr19_CLINVAR_VEP_CHR19(installWire('D1_chr19', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr20_CLINVAR_VEP_CHR20(installWire('D1_chr20', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr21_CLINVAR_VEP_CHR21(installWire('D1_chr21', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chr22_CLINVAR_VEP_CHR22(installWire('D1_chr22', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chrX_CLINVAR_VEP_CHRX(installWire('D1_chrX', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chrY_CLINVAR_VEP_CHRY(installWire('D1_chrY', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_chrOther_CLINVAR_VEP_CHROTHER(installWire('D1_chrOther', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt]))
    STEP_D1_CLINVAR(installWire('D1', plan, [STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_C1_GNOMAD.out.receipt, STEP_D0_CLINVAR_DOWNLOAD.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt, STEP_D1_chr1_CLINVAR_VEP_CHR1.out.receipt, STEP_D1_chr2_CLINVAR_VEP_CHR2.out.receipt, STEP_D1_chr3_CLINVAR_VEP_CHR3.out.receipt, STEP_D1_chr4_CLINVAR_VEP_CHR4.out.receipt, STEP_D1_chr5_CLINVAR_VEP_CHR5.out.receipt, STEP_D1_chr6_CLINVAR_VEP_CHR6.out.receipt, STEP_D1_chr7_CLINVAR_VEP_CHR7.out.receipt, STEP_D1_chr8_CLINVAR_VEP_CHR8.out.receipt, STEP_D1_chr9_CLINVAR_VEP_CHR9.out.receipt, STEP_D1_chr10_CLINVAR_VEP_CHR10.out.receipt, STEP_D1_chr11_CLINVAR_VEP_CHR11.out.receipt, STEP_D1_chr12_CLINVAR_VEP_CHR12.out.receipt, STEP_D1_chr13_CLINVAR_VEP_CHR13.out.receipt, STEP_D1_chr14_CLINVAR_VEP_CHR14.out.receipt, STEP_D1_chr15_CLINVAR_VEP_CHR15.out.receipt, STEP_D1_chr16_CLINVAR_VEP_CHR16.out.receipt, STEP_D1_chr17_CLINVAR_VEP_CHR17.out.receipt, STEP_D1_chr18_CLINVAR_VEP_CHR18.out.receipt, STEP_D1_chr19_CLINVAR_VEP_CHR19.out.receipt, STEP_D1_chr20_CLINVAR_VEP_CHR20.out.receipt, STEP_D1_chr21_CLINVAR_VEP_CHR21.out.receipt, STEP_D1_chr22_CLINVAR_VEP_CHR22.out.receipt, STEP_D1_chrX_CLINVAR_VEP_CHRX.out.receipt, STEP_D1_chrY_CLINVAR_VEP_CHRY.out.receipt, STEP_D1_chrOther_CLINVAR_VEP_CHROTHER.out.receipt]))
    STEP_C10_DAS_DOMAINS(installNative('C10', plan, [STEP_C6_ALPHAMISSENSE_DOMAIN.out.receipt, STEP_C8_INTERPRO_MAPPING.out.receipt, STEP_C4_ALPHAMISSENSE_GNOMAD.out.receipt, STEP_D1_CLINVAR.out.receipt], ['alphamissense_pd_stat', 'interpro_mapping_pickle', 'alphamissense_vep_vcf', 'clinvar_vcf']))
    STEP_C11_PM1_REGIONS(installNative('C11', plan, [STEP_C10_DAS_DOMAINS.out.receipt, STEP_C6_ALPHAMISSENSE_DOMAIN.out.receipt, STEP_A5_NM_TO_ENST.out.receipt], ['all_intolerant_domains', 'nm_to_enst_tsv']))
    STEP_D6_CLINVAR_PATHO_AF(installNative('D6', plan, [STEP_D1_CLINVAR.out.receipt], ['clinvar_vcf']))
    STEP_D8_CLINVAR_PM5(installNative('D8', plan, [STEP_D1_CLINVAR.out.receipt], ['clinvar_vcf']))
    STEP_D9_CLINGEN_MAP(installNative('D9', plan, [STEP_D1_CLINVAR.out.receipt], ['clinvar_vcf']))
    STEP_E1_CADD(installWire('E1', plan, []))
    STEP_F1_HPO_MONDO(installNative('F1', plan, [], []))
    STEP_F4_GENE_PATHOGENIC_MECHANISM(installNative('F4', plan, [STEP_F1_HPO_MONDO.out.receipt], ['mondo_disease_scope_registry']))
    STEP_F5_GOFCARDS_EXACT_GOF(installNative('F5', plan, [], []))
    STEP_F6_MECHANISM_SOURCES(installNative('F6', plan, [], []))
    STEP_F7_GOFCARDS_CLINVAR_INJECTION(installNative('F7', plan, [STEP_F5_GOFCARDS_EXACT_GOF.out.receipt, STEP_F6_MECHANISM_SOURCES.out.receipt], ['gofcards_exact_gof_cache', 'clinvar_vcv_xml']))
    STEP_F8_GENE_NONLOF_MECHANISM(installNative('F8', plan, [STEP_F6_MECHANISM_SOURCES.out.receipt, STEP_F7_GOFCARDS_CLINVAR_INJECTION.out.receipt], ['gofcards_exact_gof_cache', 'clinvar_vcv_xml']))
    STEP_F9_HPO_CONDITION_MECHANISM(installNative('F9', plan, [STEP_F1_HPO_MONDO.out.receipt, STEP_F4_GENE_PATHOGENIC_MECHANISM.out.receipt, STEP_F7_GOFCARDS_CLINVAR_INJECTION.out.receipt, STEP_D9_CLINGEN_MAP.out.receipt], ['hpo_assertions', 'ddg2p_mechanism_evidence', 'gofcards_exact_gof_cache', 'gene_dosage_sensitivity']))
    STEP_F10_HPO_GENE_RANKING(installNative('F10', plan, [], []))
    def all = Channel.empty().mix(STEP_A1_REFERENCE_FASTA.out.receipt, STEP_A2_VEP_CACHE.out.receipt, STEP_A4_CDS_FASTA.out.receipt, STEP_A5_NM_TO_ENST.out.receipt, STEP_B1_UTR_ANNOTATOR.out.receipt, STEP_B2_CONSERVATION.out.receipt, STEP_B3_ALPHAMISSENSE_TSV.out.receipt, STEP_B4_SPLICEAI.out.receipt, STEP_B5_SPLICEVAULT.out.receipt, STEP_B6_MAVEDB.out.receipt, STEP_B7_PRIMATEAI.out.receipt, STEP_B8_LOFTEE.out.receipt, STEP_B9_PEXT.out.receipt, STEP_C1_chr1_GNOMAD_DOWNLOAD_CHR1.out.receipt, STEP_C1_chr2_GNOMAD_DOWNLOAD_CHR2.out.receipt, STEP_C1_chr3_GNOMAD_DOWNLOAD_CHR3.out.receipt, STEP_C1_chr4_GNOMAD_DOWNLOAD_CHR4.out.receipt, STEP_C1_chr5_GNOMAD_DOWNLOAD_CHR5.out.receipt, STEP_C1_chr6_GNOMAD_DOWNLOAD_CHR6.out.receipt, STEP_C1_chr7_GNOMAD_DOWNLOAD_CHR7.out.receipt, STEP_C1_chr8_GNOMAD_DOWNLOAD_CHR8.out.receipt, STEP_C1_chr9_GNOMAD_DOWNLOAD_CHR9.out.receipt, STEP_C1_chr10_GNOMAD_DOWNLOAD_CHR10.out.receipt, STEP_C1_chr11_GNOMAD_DOWNLOAD_CHR11.out.receipt, STEP_C1_chr12_GNOMAD_DOWNLOAD_CHR12.out.receipt, STEP_C1_chr13_GNOMAD_DOWNLOAD_CHR13.out.receipt, STEP_C1_chr14_GNOMAD_DOWNLOAD_CHR14.out.receipt, STEP_C1_chr15_GNOMAD_DOWNLOAD_CHR15.out.receipt, STEP_C1_chr16_GNOMAD_DOWNLOAD_CHR16.out.receipt, STEP_C1_chr17_GNOMAD_DOWNLOAD_CHR17.out.receipt, STEP_C1_chr18_GNOMAD_DOWNLOAD_CHR18.out.receipt, STEP_C1_chr19_GNOMAD_DOWNLOAD_CHR19.out.receipt, STEP_C1_chr20_GNOMAD_DOWNLOAD_CHR20.out.receipt, STEP_C1_chr21_GNOMAD_DOWNLOAD_CHR21.out.receipt, STEP_C1_chr22_GNOMAD_DOWNLOAD_CHR22.out.receipt, STEP_C1_chrX_GNOMAD_DOWNLOAD_CHRX.out.receipt, STEP_C1_chrY_GNOMAD_DOWNLOAD_CHRY.out.receipt, STEP_C1_GNOMAD.out.receipt, STEP_C3_ALPHAMISSENSE_VEP.out.receipt, STEP_C4_ALPHAMISSENSE_GNOMAD.out.receipt, STEP_C6_ALPHAMISSENSE_DOMAIN.out.receipt, STEP_C8_INTERPRO_MAPPING.out.receipt, STEP_D0_CLINVAR_DOWNLOAD.out.receipt, STEP_D1_split_CLINVAR_VEP_SPLIT.out.receipt, STEP_D1_chr1_CLINVAR_VEP_CHR1.out.receipt, STEP_D1_chr2_CLINVAR_VEP_CHR2.out.receipt, STEP_D1_chr3_CLINVAR_VEP_CHR3.out.receipt, STEP_D1_chr4_CLINVAR_VEP_CHR4.out.receipt, STEP_D1_chr5_CLINVAR_VEP_CHR5.out.receipt, STEP_D1_chr6_CLINVAR_VEP_CHR6.out.receipt, STEP_D1_chr7_CLINVAR_VEP_CHR7.out.receipt, STEP_D1_chr8_CLINVAR_VEP_CHR8.out.receipt, STEP_D1_chr9_CLINVAR_VEP_CHR9.out.receipt, STEP_D1_chr10_CLINVAR_VEP_CHR10.out.receipt, STEP_D1_chr11_CLINVAR_VEP_CHR11.out.receipt, STEP_D1_chr12_CLINVAR_VEP_CHR12.out.receipt, STEP_D1_chr13_CLINVAR_VEP_CHR13.out.receipt, STEP_D1_chr14_CLINVAR_VEP_CHR14.out.receipt, STEP_D1_chr15_CLINVAR_VEP_CHR15.out.receipt, STEP_D1_chr16_CLINVAR_VEP_CHR16.out.receipt, STEP_D1_chr17_CLINVAR_VEP_CHR17.out.receipt, STEP_D1_chr18_CLINVAR_VEP_CHR18.out.receipt, STEP_D1_chr19_CLINVAR_VEP_CHR19.out.receipt, STEP_D1_chr20_CLINVAR_VEP_CHR20.out.receipt, STEP_D1_chr21_CLINVAR_VEP_CHR21.out.receipt, STEP_D1_chr22_CLINVAR_VEP_CHR22.out.receipt, STEP_D1_chrX_CLINVAR_VEP_CHRX.out.receipt, STEP_D1_chrY_CLINVAR_VEP_CHRY.out.receipt, STEP_D1_chrOther_CLINVAR_VEP_CHROTHER.out.receipt, STEP_D1_CLINVAR.out.receipt, STEP_C10_DAS_DOMAINS.out.receipt, STEP_C11_PM1_REGIONS.out.receipt, STEP_D6_CLINVAR_PATHO_AF.out.receipt, STEP_D8_CLINVAR_PM5.out.receipt, STEP_D9_CLINGEN_MAP.out.receipt, STEP_E1_CADD.out.receipt, STEP_F1_HPO_MONDO.out.receipt, STEP_F4_GENE_PATHOGENIC_MECHANISM.out.receipt, STEP_F5_GOFCARDS_EXACT_GOF.out.receipt, STEP_F6_MECHANISM_SOURCES.out.receipt, STEP_F7_GOFCARDS_CLINVAR_INJECTION.out.receipt, STEP_F8_GENE_NONLOF_MECHANISM.out.receipt, STEP_F9_HPO_CONDITION_MECHANISM.out.receipt, STEP_F10_HPO_GENE_RANKING.out.receipt).collect()
    def toPublish = plan.combine(all.map { [it] }).map { p, receipts ->
        def names = receipts.collect { it.toString() }.sort()
        // Publishing writes every step's paths into one configuration, so it
        // keys on their lineages (folders), and on the publishing code.
        def lineages = names.collect { r -> new groovy.json.JsonSlurper().parse(new File(r)).lineage }
        def digest = new groovy.json.JsonSlurper().parse(p.toFile()).publish_digest
        [installSha(groovy.json.JsonOutput.toJson([lineages: lineages, publish: digest])).take(16), p, names]
    }
    INSTALL_PUBLISH(toPublish)
    emit:
    published = INSTALL_PUBLISH.out.published
}
