#!/usr/bin/env python3
"""
Build the PriVA cache/resource dependency graph and render it.

Node kinds:
  * data    -- an install-time cache file (reference/annotation resource)
  * runtime -- a cohort-specific input or pipeline runtime output
  * install -- an install_utils.sh shell function, or a python script that only
               builds a cache (no runtime output)
  * prod    -- a runtime function (annotation / filter / prioritize / ACMG)
               that reads the cohort input and emits pipeline output

Edge convention:
  * install -> prod    : "invoke"   (parent function calls child function)
  * prod    -> data    : "generate" (script produces a cache)
  * prod    -> runtime : "generate" (pipeline produces a runtime output)
  * data    -> prod    : "input"    (cache is read by a runtime function)
  * runtime -> prod    : "input"    (cohort file is read by a runtime function)

Every non-empty file path in config.yaml is a node. No node is isolated, and
every cache is reachable from a runtime (prod) function.
"""
import subprocess
import sys
import os
import xml.etree.ElementTree as ET

REPO = "/paedyl01/disk1/yangyxt/PriVA"
OUT_GRAPHML = os.path.join(REPO, "priva_cache_dependency.graphml")
OUT_DOT = os.path.join(REPO, "priva_cache_dependency.dot")
OUT_PNG = os.path.join(REPO, "priva_cache_dependency.png")

INSTALL = "install_utils.sh"

# ----------------------------------------------------------- cache data nodes
# id -> (label, path)
DATA = {
    "ref_genome":              ("ucsc ref FASTA", "ucsc.hg19.fasta"),
    "hg38_hg19_chain":         ("hg38ToHg19 chain", "hg38ToHg19.over.chain.gz"),
    "contig_map":              ("GRC->ucsc contig map", "GRC_to_ucsc.contig.map.txt"),
    "vep_cache":               ("VEP cache", "VEP_caches"),
    "vep_plugins":             ("VEP plugins", "VEP_plugins"),
    "vep_plugins_cachedir":    ("VEP plugin caches", "VEP_plugins_caches"),
    "conda_env":               ("priva_acmg conda env", "env priva_acmg"),
    "conda_env_yaml":          ("acmg_conda_rough.yml", "acmg_conda_rough.yml"),
    "gnomad_vcf_chrX":         ("gnomAD v4.1 VCFs", "gnomad.joint.v4.1.sites.hg19.chrX.vcf.gz"),
    # ClinVar
    "clinvar_vcf":             ("ClinVar VEP VCF", "clinvar.hg19.vep.vcf.gz"),
    "clinvar_gene_stat":       ("ClinVar 2-star gene stat", "clinvar_2star_stats.pkl.gz"),
    "clinvar_aa_stat":         ("ClinVar aa-change stat", "clinvar.hg19.vep.aa_change.stats.pkl"),
    "clinvar_splice_stat":     ("ClinVar splice-change stat", "clinvar.hg19.vep.splice_change.stats.pkl"),
    "clinvar_patho_af_stat":   ("ClinVar patho AF stat", "clinvar.hg19.vep.patho_af_stat.pkl"),
    "clinvar_patho_exon_af_stat": ("ClinVar patho exon-AF stat", "clinvar.hg19.vep.patho_exon_af_stat.pkl"),
    # VEP plugin caches
    "utr_annotator_file":      ("UTRAnnotator cache", "uORF_5UTR_GRCh37_PUBLIC.txt"),
    "loeuf_prescore":          ("LOEUF prescore", "loeuf_dataset.tsv.gz"),
    "loftee_repo":             ("LoFtee repo", "loftee-hg19"),
    "human_ancestor_fasta":    ("LoFtee human ancestor FA", "human_ancestor.fa.gz"),
    "loftee_conservation_file":("LoFtee GERP sql", "phylocsf_gerp.sql"),
    "spliceai_snv_prescore":   ("SpliceAI SNV prescore", "spliceai_scores.raw.snv.hg19.vcf.gz"),
    "spliceai_indel_prescore": ("SpliceAI indel prescore", "spliceai_scores.raw.indel.hg19.vcf.gz"),
    "splicevault_prescore":    ("SpliceVault prescore", "SpliceVault_data_GRCh37.tsv.gz"),
    "primateai_prescore":      ("PrimateAI prescore", "PrimateAI_scores_v0.2_GRCh37_sorted.tsv.bgz"),
    "conservation_file":       ("GERP conservation bw", "gerp_conservation_scores.homo_sapiens.GRCh37.bw"),
    "cadd_base_dir":           ("CADD scripts + cache", "CADD-scripts-1.7.2"),
    "cadd_prescore":           ("CADD prescores (SNV+indel inclAnno)", "whole_genome_SNVs_inclAnno.tsv.gz + indel_inclAnno.tsv.gz"),
    "mavedb_file":             ("MaveDB variants", "MaveDB_variants.hg38.tsv.gz"),
    "pext_mean_bw":            ("pext mean BigWig", "pext mean BigWig"),
    "pext_tissue_bw_dir":      ("pext tissue BigWigs", "pext tissue BigWig dir"),
    "cds_fasta_file":          ("Ensembl CDS FASTA", "Homo_sapiens.GRCh37.cds.all.fa.gz"),
    "interpro_mapping_pickle": ("InterPro entry mapping", "Interpro_entry_mapping.pkl.gz"),
    "clingen_map":             ("ClinGen map", "clingen_map.hg19.pkl.gz"),
    "gene_dosage_sensitivity": ("ClinGen dosage sensitivity", "gene_dosage_sensitivity.hg19.tsv"),
    # AlphaMissense / DAS / PM1
    "alphamissense_prescore":  ("AlphaMissense prescore", "AlphaMissense_hg19.tsv.gz"),
    "alphamissense_vcf":       ("AlphaMissense VCF", "AlphaMissense_hg19.vcf.gz"),
    "alphamissense_vep_vcf":   ("AlphaMissense VEP VCF", "AlphaMissense_hg19.vep.vcf.gz"),
    "alphamissense_intolerant_motifs": ("AM intolerant motifs", "AlphaMissense_hg19.kde.pkl"),
    "alphamissense_pd_stat":   ("AM protein-domain stats", "AlphaMissense_hg19.vep.prot.domain.stats.pkl"),
    "alphamissense_tranx_domain_map": ("AM transcript-domain map", "transcript_pp_domain_mapping.hg19.pkl"),
    "alphamissense_intolerant_domains": ("AM domain tolerance TSV", "domain_tolerance_analysis.hg19.tsv"),
    "all_intolerant_domains":  ("DAS all intolerant domains", "all_intolerant_domains.hg19.pkl"),
    "pm1_regions_pkl":         ("PM1 regions", "pm1_regions.hg19.pkl"),
    "hcseeker_spots_tsv":      ("HCSeeker hotspots", "HC_spots.hg19.tsv"),
    "rmc_tsv":                 ("gnomAD RMC flat TSV", "gnomad_v2.1.1_rmc_flat.hg19.tsv"),
    # gene mechanism
    "gene_mechanism_cache_dir":("gene mechanism cache", "data/patho_mechanism"),
    "gene_mechanism_raw_dir":  ("gene mechanism raw", "gene_pathogenic_mechanism/raw"),
    "ddg2p_mechanism_evidence":("DDG2P/G2P evidence TSV", "gene_pathogenic_mechanism_evidence.tsv"),
    "hgnc_table":              ("HGNC non-alt-loci table", "non_alt_loci_set.tsv"),
    "gene_nonlof_mechanism_json": ("gene non-LOF mechanism JSON", "gene_nonlof_mechanism_curated_assertions.json.gz"),
    "hpo_condition_mechanism_json": ("HPO condition-mechanism JSON", "hpo_condition_mechanism_cache.json.gz"),
    "gene_nonlof_mechanism_schema": ("gene non-LOF schema", "gene_nonlof_mechanism_curated_assertions.schema.json"),
    # GoFCards
    "gofcards_exact_gof_cache":("GoFCards exact-GOF cache", "gofcards_exact_gof.json.gz"),
    "clinvar_vcv_xml":         ("ClinVar VCV weekly XML", "ClinVarVCVRelease_00-latest_weekly.xml.gz"),
    "gofcards_mechanism_review_tsv": ("GoFCards mechanism review", "gofcards_mechanism_reviews.tsv"),
    "gofcards_hg38_fasta":     ("hg38 FASTA (GoFCards)", "ucsc.hg38.fasta"),
    "gofcards_hg19_to_hg38_chain": ("hg19->hg38 chain", "hg19ToHg38.over.chain"),
    # HPO / MONDO
    "hpo_assertions":          ("HPO gene-phenotype assertions", "genes_to_phenotype.assertions.tsv.gz"),
    "mondo_obo":               ("MONDO obo", "mondo-simple.obo"),
    "mondo_disease_scope_registry": ("MONDO disease scope", "disease_scope.tsv.gz"),
    "mondo_disease_scope_overrides": ("MONDO scope overrides", "disease_scope_overrides.tsv"),
}

# subset of DATA that are directories (not single files).
# USER_DIRS: the user specifies the full path directly (top-level mount points)
#   -> roots, no upstream function. The remaining directories are installer-
#   derived: the install function creates a named subdirectory under a
#   user-specified parent (e.g. cadd_parent_dir/CADD-scripts-1.7.2), so they
#   DO have an upstream "gen" function.
DIRS = {"vep_cache", "vep_plugins", "vep_plugins_cachedir", "cadd_base_dir",
        "loftee_repo", "pext_tissue_bw_dir", "gene_mechanism_cache_dir",
        "gene_mechanism_raw_dir"}
USER_DIRS = {"vep_cache", "vep_plugins", "vep_plugins_cachedir"}

# data nodes that are not runtime caches (legitimate leaves): the conda env is
# the runtime environment itself, its yaml is an install-time build config, and
# the AM intolerant-motifs (kde) cache is a legacy orphan superseded by the DAS
# intolerant-domain caches.
NON_CACHE = {"conda_env", "conda_env_yaml", "alphamissense_intolerant_motifs"}

# ----------------------------------------------------------- runtime io nodes
# id -> (label, description)
RUNTIME = {
    "input_vcf":      ("input cohort VCF", "cohort input VCF"),
    "ped_file":       ("PED file", "cohort pedigree"),
    "alt_disease_vcf":("alt-disease VCF", "alternative disease patients"),
    "annotated_vcf":  ("annotated cohort VCF", "VEP+CADD+gnomAD annotated output"),
    "filtered_vcf":   ("filtered cohort VCF", "per-family filtered output"),
    "acmg_result":    ("ACMG classification", "final prioritized ACMG table"),
}

# ------------------------------------------------------------- function nodes
# id -> (label, script, kind)   kind in {"install", "prod"}
FUNC = {
    # --- install_utils.sh shell functions ---
    "reference_genome_install":       ("reference_genome_install", INSTALL, "install"),
    "conda_install_vep":              ("conda_install_vep", INSTALL, "install"),
    "vep_cache_api_install":          ("vep_cache_api_install", INSTALL, "install"),
    "VEP_plugins_install":            ("VEP_plugins_install", INSTALL, "install"),
    "basic_vep_annotation":           ("basic_vep_annotation", INSTALL, "install"),
    "LoFtee_install":                 ("LoFtee_install", INSTALL, "install"),
    "UTRAnnotator_install":           ("UTRAnnotator_install", INSTALL, "install"),
    "LOEUF_install":                  ("LOEUF_install", INSTALL, "install"),
    "SpliceAI_install":               ("SpliceAI_install", INSTALL, "install"),
    "SpliceVault_install":            ("SpliceVault_install", INSTALL, "install"),
    "convert_splicevault_hg38_to_hg19": ("convert_splicevault_hg38_to_hg19", INSTALL, "install"),
    "PrimateAI_install":              ("PrimateAI_install", INSTALL, "install"),
    "Conservation_install":           ("Conservation_install", INSTALL, "install"),
    "CADD_install":                   ("CADD_install", INSTALL, "install"),
    "MaveDB_install":                 ("MaveDB_install", INSTALL, "install"),
    "MaveDB_deployment":              ("MaveDB_deployment", INSTALL, "install"),
    "pext_install":                   ("pext_install", INSTALL, "install"),
    "gnomAD_install":                 ("gnomAD_install", INSTALL, "install"),
    "gnomAD_liftover":                ("gnomAD_liftover", INSTALL, "install"),
    "ClinVar_VCF_deploy":             ("ClinVar_VCF_deploy", INSTALL, "install"),
    "ClinVar_Gene_stat":              ("ClinVar_Gene_stat", INSTALL, "install"),
    "ClinVar_AA_stat":                ("ClinVar_AA_stat", INSTALL, "install"),
    "ClinVar_patho_AF_stat":          ("ClinVar_patho_AF_stat", INSTALL, "install"),
    "AlphaMissense_install":          ("AlphaMissense_install", INSTALL, "install"),
    "AlphaMissense_anno":             ("AlphaMissense_anno", INSTALL, "install"),
    "AlphaMissense_anno_gnomAD":      ("AlphaMissense_anno_gnomAD", INSTALL, "install"),
    "AlphaMissense_pick_missense_constrained_segments": ("AlphaMissense_pick_missense_constrained_segments", INSTALL, "install"),
    "AlphaMissense_stat":             ("AlphaMissense_stat", INSTALL, "install"),
    "AlphaMissense_pick_intolerant_domains": ("AlphaMissense_pick_intolerant_domains", INSTALL, "install"),
    "prepare_pm1_regions_fn":         ("prepare_pm1_regions", INSTALL, "install"),
    "InterPro_parsing":               ("InterPro_parsing", INSTALL, "install"),
    "ClinGen_deploy":                 ("ClinGen_deploy", INSTALL, "install"),
    "CDS_FASTA_install":              ("CDS_FASTA_install", INSTALL, "install"),
    "gene_pathogenic_mechanism_cache_install": ("gene_pathogenic_mechanism_cache_install", INSTALL, "install"),
    "gene_nonlof_mechanism_cache_install":     ("gene_nonlof_mechanism_cache_install", INSTALL, "install"),
    "hpo_condition_mechanism_cache_install":   ("hpo_condition_mechanism_cache_install", INSTALL, "install"),
    "gofcards_exact_gof_cache_install":        ("gofcards_exact_gof_cache_install", INSTALL, "install"),
    "gofcards_clinvar_injection_install":      ("gofcards_clinvar_injection_install", INSTALL, "install"),
    "mondo_hpo_scope_install":                 ("mondo_hpo_scope_install", INSTALL, "install"),
    "mechanism_resource_install":              ("mechanism_resource_install", INSTALL, "install"),
    # --- install-time python cache builders ---
    "alphmis_tsv2vcf.py":                     ("alphmis_tsv2vcf", "alphmis_tsv2vcf.py", "install"),
    "am_pick_missense_constrained_segments.py": ("am_pick_missense_constrained_segments", "am_pick_missense_constrained_segments.py", "install"),
    "stat_protein_domain_amscores.py":        ("stat_protein_domain_amscores", "stat_protein_domain_amscores.py", "install"),
    "am_pick_intolerant_domains.py":          ("am_pick_intolerant_domains", "am_pick_intolerant_domains.py", "install"),
    "prepare_pm1_regions.py":                 ("prepare_pm1_regions", "prepare_pm1_regions.py", "install"),
    "clinvar_stat_variants.py":               ("clinvar_stat_variants", "clinvar_stat_variants.py", "install"),
    "stat_aachange_clinvar.py":               ("stat_aachange_clinvar", "stat_aachange_clinvar.py", "install"),
    "stat_gene_patho_afs.py":                 ("stat_gene_patho_afs", "stat_gene_patho_afs.py", "install"),
    "stat_exon_patho_afs.py":                 ("stat_exon_patho_afs", "stat_exon_patho_afs.py", "install"),
    "build_gene_pathogenic_mechanism_cache.py": ("build_gene_pathogenic_mechanism_cache", "build_gene_pathogenic_mechanism_cache.py", "install"),
    "build_gene_nonlof_mechanism_cache.py":   ("build_gene_nonlof_mechanism_cache", "build_gene_nonlof_mechanism_cache.py", "install"),
    "build_hpo_condition_mechanism_cache.py": ("build_hpo_condition_mechanism_cache", "build_hpo_condition_mechanism_cache.py", "install"),
    "build_gofcards_exact_gof_cache.py":      ("build_gofcards_exact_gof_cache", "build_gofcards_exact_gof_cache.py", "install"),
    "inject_clinvar_into_gofcards.py":        ("inject_clinvar_into_gofcards", "inject_clinvar_into_gofcards.py", "install"),
    "collapse_HPO_anno.py":                   ("collapse_HPO_anno", "collapse_HPO_anno.py", "install"),
    "build_mondo_disease_scope.py":           ("build_mondo_disease_scope", "build_mondo_disease_scope.py", "install"),
    # --- runtime production functions ---
    "annotation_vcf.sh":              ("annotation_vcf (annotate)", "annotation_vcf.sh", "prod"),
    "filtration_vcf_per_fam.sh":      ("filtration_vcf_per_fam (filter)", "filtration_vcf_per_fam.sh", "prod"),
    "prioritization_vcf_per_fam.sh":  ("prioritization_vcf_per_fam (ACMG)", "prioritization_vcf_per_fam.sh", "prod"),
    "gene_mechanism_common.py":       ("gene_mechanism_common", "gene_mechanism_common.py", "prod"),
    "acmg_criteria_assign.py":        ("acmg_criteria_assign", "acmg_criteria_assign.py", "prod"),
    "splicing_var_analysis.py":       ("splicing_var_analysis", "splicing_var_analysis.py", "prod"),
    "utr_anno_interpret.py":          ("utr_anno_interpret", "utr_anno_interpret.py", "prod"),
    "acmg_pvs1_null_variant.py":      ("acmg_pvs1_null_variant", "acmg_pvs1_null_variant.py", "prod"),
    "acmg_pm2_bs1_ba1_frequency.py":  ("acmg_pm2_bs1_ba1_frequency", "acmg_pm2_bs1_ba1_frequency.py", "prod"),
    "acmg_pm1_pp2_bp1_missense.py":   ("acmg_pm1_pp2_bp1_missense", "acmg_pm1_pp2_bp1_missense.py", "prod"),
}

# ------------------------------------------------------------------- edges
# (source, target, kind)  kind in {"gen", "in", "invoke"}
EDGES = [
    # -- genome / env / VEP --
    ("reference_genome_install", "ref_genome", "gen"),
    ("conda_env_yaml", "conda_install_vep", "in"),
    ("conda_install_vep", "conda_env", "gen"),
    ("vep_cache", "vep_cache_api_install", "in"),
    ("vep_plugins", "VEP_plugins_install", "in"),
    ("vep_plugins_cachedir", "VEP_plugins_install", "in"),
    ("VEP_plugins_install", "UTRAnnotator_install", "invoke"),
    ("VEP_plugins_install", "LOEUF_install", "invoke"),
    ("VEP_plugins_install", "AlphaMissense_install", "invoke"),
    ("VEP_plugins_install", "AlphaMissense_anno", "invoke"),
    ("VEP_plugins_install", "AlphaMissense_pick_missense_constrained_segments", "invoke"),
    ("VEP_plugins_install", "SpliceAI_install", "invoke"),
    ("VEP_plugins_install", "SpliceVault_install", "invoke"),
    ("VEP_plugins_install", "MaveDB_install", "invoke"),
    ("VEP_plugins_install", "PrimateAI_install", "invoke"),
    ("VEP_plugins_install", "Conservation_install", "invoke"),
    ("ref_genome", "basic_vep_annotation", "in"),
    ("vep_cache", "basic_vep_annotation", "in"),
    ("vep_plugins", "basic_vep_annotation", "in"),
    # -- LoFtee / plugin downloads --
    ("LoFtee_install", "loftee_repo", "gen"),
    ("LoFtee_install", "human_ancestor_fasta", "gen"),
    ("LoFtee_install", "loftee_conservation_file", "gen"),
    ("UTRAnnotator_install", "utr_annotator_file", "gen"),
    ("LOEUF_install", "loeuf_prescore", "gen"),
    ("SpliceAI_install", "spliceai_snv_prescore", "gen"),
    ("SpliceAI_install", "spliceai_indel_prescore", "gen"),
    ("SpliceVault_install", "convert_splicevault_hg38_to_hg19", "invoke"),
    ("convert_splicevault_hg38_to_hg19", "splicevault_prescore", "gen"),
    ("PrimateAI_install", "primateai_prescore", "gen"),
    ("Conservation_install", "conservation_file", "gen"),
    ("CADD_install", "cadd_base_dir", "gen"),
    ("CADD_install", "cadd_prescore", "gen"),
    ("MaveDB_install", "MaveDB_deployment", "invoke"),
    ("MaveDB_deployment", "mavedb_file", "gen"),
    ("pext_install", "pext_mean_bw", "gen"),
    ("pext_install", "pext_tissue_bw_dir", "gen"),
    # -- gnomAD --
    ("gnomAD_install", "gnomAD_liftover", "invoke"),
    ("hg38_hg19_chain", "gnomAD_liftover", "in"),
    ("gnomAD_liftover", "gnomad_vcf_chrX", "gen"),
    # -- InterPro / ClinGen / CDS --
    ("InterPro_parsing", "interpro_mapping_pickle", "gen"),
    ("ClinGen_deploy", "clingen_map", "gen"),
    ("ClinGen_deploy", "gene_dosage_sensitivity", "gen"),
    ("CDS_FASTA_install", "cds_fasta_file", "gen"),
    # -- ClinVar chain --
    ("contig_map", "ClinVar_VCF_deploy", "in"),
    ("ref_genome", "ClinVar_VCF_deploy", "in"),
    ("vep_cache", "ClinVar_VCF_deploy", "in"),
    ("spliceai_snv_prescore", "ClinVar_VCF_deploy", "in"),
    ("spliceai_indel_prescore", "ClinVar_VCF_deploy", "in"),
    ("gnomad_vcf_chrX", "ClinVar_VCF_deploy", "in"),
    ("ClinVar_VCF_deploy", "basic_vep_annotation", "invoke"),
    ("ClinVar_VCF_deploy", "clinvar_vcf", "gen"),
    ("clinvar_vcf", "ClinVar_Gene_stat", "in"),
    ("ClinVar_Gene_stat", "clinvar_stat_variants.py", "invoke"),
    ("clinvar_vcf", "clinvar_stat_variants.py", "in"),
    ("clinvar_stat_variants.py", "clinvar_gene_stat", "gen"),
    ("clinvar_vcf", "ClinVar_AA_stat", "in"),
    ("ClinVar_AA_stat", "stat_aachange_clinvar.py", "invoke"),
    ("clinvar_vcf", "stat_aachange_clinvar.py", "in"),
    ("stat_aachange_clinvar.py", "clinvar_aa_stat", "gen"),
    ("stat_aachange_clinvar.py", "clinvar_splice_stat", "gen"),
    ("clinvar_vcf", "ClinVar_patho_AF_stat", "in"),
    ("ClinVar_patho_AF_stat", "stat_gene_patho_afs.py", "invoke"),
    ("ClinVar_patho_AF_stat", "stat_exon_patho_afs.py", "invoke"),
    ("clinvar_vcf", "stat_gene_patho_afs.py", "in"),
    ("clinvar_vcf", "stat_exon_patho_afs.py", "in"),
    ("stat_gene_patho_afs.py", "clinvar_patho_af_stat", "gen"),
    ("stat_exon_patho_afs.py", "clinvar_patho_exon_af_stat", "gen"),
    # -- AlphaMissense chain --
    ("AlphaMissense_install", "alphamissense_prescore", "gen"),
    ("alphamissense_prescore", "AlphaMissense_anno", "in"),
    ("ref_genome", "AlphaMissense_anno", "in"),
    ("vep_cache", "AlphaMissense_anno", "in"),
    ("conservation_file", "AlphaMissense_anno", "in"),
    ("AlphaMissense_anno", "alphmis_tsv2vcf.py", "invoke"),
    ("alphamissense_prescore", "alphmis_tsv2vcf.py", "in"),
    ("alphmis_tsv2vcf.py", "alphamissense_vcf", "gen"),
    ("AlphaMissense_anno", "basic_vep_annotation", "invoke"),
    ("alphamissense_vcf", "basic_vep_annotation", "in"),
    ("basic_vep_annotation", "alphamissense_vep_vcf", "gen"),
    ("alphamissense_vep_vcf", "AlphaMissense_pick_missense_constrained_segments", "in"),
    ("AlphaMissense_pick_missense_constrained_segments", "am_pick_missense_constrained_segments.py", "invoke"),
    ("alphamissense_vep_vcf", "am_pick_missense_constrained_segments.py", "in"),
    ("hcseeker_spots_tsv", "am_pick_missense_constrained_segments.py", "in"),
    ("rmc_tsv", "am_pick_missense_constrained_segments.py", "in"),
    ("am_pick_missense_constrained_segments.py", "alphamissense_intolerant_motifs", "gen"),
    ("alphamissense_vep_vcf", "AlphaMissense_stat", "in"),
    ("AlphaMissense_stat", "stat_protein_domain_amscores.py", "invoke"),
    ("alphamissense_vep_vcf", "stat_protein_domain_amscores.py", "in"),
    ("stat_protein_domain_amscores.py", "alphamissense_pd_stat", "gen"),
    ("alphamissense_vep_vcf", "AlphaMissense_anno_gnomAD", "in"),
    ("gnomad_vcf_chrX", "AlphaMissense_anno_gnomAD", "in"),
    # -- DAS pick intolerant domains --
    ("AlphaMissense_pick_intolerant_domains", "am_pick_intolerant_domains.py", "invoke"),
    ("alphamissense_pd_stat", "am_pick_intolerant_domains.py", "in"),
    ("alphamissense_vep_vcf", "am_pick_intolerant_domains.py", "in"),
    ("hcseeker_spots_tsv", "am_pick_intolerant_domains.py", "in"),
    ("clinvar_vcf", "am_pick_intolerant_domains.py", "in"),
    ("am_pick_intolerant_domains.py", "alphamissense_tranx_domain_map", "gen"),
    ("am_pick_intolerant_domains.py", "alphamissense_intolerant_domains", "gen"),
    ("am_pick_intolerant_domains.py", "all_intolerant_domains", "gen"),
    # -- PM1 regions --
    ("prepare_pm1_regions_fn", "prepare_pm1_regions.py", "invoke"),
    ("all_intolerant_domains", "prepare_pm1_regions.py", "in"),
    ("hcseeker_spots_tsv", "prepare_pm1_regions.py", "in"),
    ("rmc_tsv", "prepare_pm1_regions.py", "in"),
    ("alphamissense_pd_stat", "prepare_pm1_regions.py", "in"),
    ("prepare_pm1_regions.py", "pm1_regions_pkl", "gen"),
    # -- gene mechanism --
    ("mechanism_resource_install", "gene_pathogenic_mechanism_cache_install", "invoke"),
    ("mechanism_resource_install", "gene_nonlof_mechanism_cache_install", "invoke"),
    ("mechanism_resource_install", "hpo_condition_mechanism_cache_install", "invoke"),
    ("mechanism_resource_install", "gofcards_exact_gof_cache_install", "invoke"),
    ("mechanism_resource_install", "gofcards_clinvar_injection_install", "invoke"),
    ("mechanism_resource_install", "mondo_hpo_scope_install", "invoke"),
    ("gene_pathogenic_mechanism_cache_install", "gene_mechanism_raw_dir", "gen"),
    ("gene_pathogenic_mechanism_cache_install", "build_gene_pathogenic_mechanism_cache.py", "invoke"),
    ("gene_mechanism_raw_dir", "build_gene_pathogenic_mechanism_cache.py", "in"),
    ("build_gene_pathogenic_mechanism_cache.py", "gene_mechanism_cache_dir", "gen"),
    ("build_gene_pathogenic_mechanism_cache.py", "ddg2p_mechanism_evidence", "gen"),
    ("hgnc_table", "gene_nonlof_mechanism_cache_install", "in"),
    ("gene_nonlof_mechanism_cache_install", "build_gene_nonlof_mechanism_cache.py", "invoke"),
    ("hgnc_table", "build_gene_nonlof_mechanism_cache.py", "in"),
    ("ddg2p_mechanism_evidence", "build_gene_nonlof_mechanism_cache.py", "in"),
    ("gene_nonlof_mechanism_schema", "build_gene_nonlof_mechanism_cache.py", "in"),
    ("build_gene_nonlof_mechanism_cache.py", "gene_nonlof_mechanism_json", "gen"),
    ("hpo_condition_mechanism_cache_install", "build_hpo_condition_mechanism_cache.py", "invoke"),
    ("gene_nonlof_mechanism_json", "build_hpo_condition_mechanism_cache.py", "in"),
    ("hpo_assertions", "build_hpo_condition_mechanism_cache.py", "in"),
    ("build_hpo_condition_mechanism_cache.py", "hpo_condition_mechanism_json", "gen"),
    # -- GoFCards --
    ("gofcards_exact_gof_cache_install", "build_gofcards_exact_gof_cache.py", "invoke"),
    ("gofcards_mechanism_review_tsv", "build_gofcards_exact_gof_cache.py", "in"),
    ("hgnc_table", "build_gofcards_exact_gof_cache.py", "in"),
    ("gofcards_hg19_to_hg38_chain", "build_gofcards_exact_gof_cache.py", "in"),
    ("gofcards_hg38_fasta", "build_gofcards_exact_gof_cache.py", "in"),
    ("build_gofcards_exact_gof_cache.py", "gofcards_exact_gof_cache", "gen"),
    ("gofcards_clinvar_injection_install", "inject_clinvar_into_gofcards.py", "invoke"),
    ("clinvar_vcv_xml", "inject_clinvar_into_gofcards.py", "in"),
    ("gofcards_exact_gof_cache", "inject_clinvar_into_gofcards.py", "in"),
    ("inject_clinvar_into_gofcards.py", "gofcards_exact_gof_cache", "gen"),
    # -- HPO / MONDO --
    ("mondo_hpo_scope_install", "collapse_HPO_anno.py", "invoke"),
    ("mondo_hpo_scope_install", "build_mondo_disease_scope.py", "invoke"),
    ("mondo_hpo_scope_install", "mondo_obo", "gen"),
    ("collapse_HPO_anno.py", "hpo_assertions", "gen"),
    ("mondo_obo", "build_mondo_disease_scope.py", "in"),
    ("mondo_disease_scope_overrides", "build_mondo_disease_scope.py", "in"),
    ("build_mondo_disease_scope.py", "mondo_disease_scope_registry", "gen"),
    # -- runtime pipeline flow --
    ("input_vcf", "annotation_vcf.sh", "in"),
    ("ped_file", "annotation_vcf.sh", "in"),
    ("annotation_vcf.sh", "annotated_vcf", "gen"),
    ("annotated_vcf", "filtration_vcf_per_fam.sh", "in"),
    ("filtration_vcf_per_fam.sh", "filtered_vcf", "gen"),
    ("filtered_vcf", "prioritization_vcf_per_fam.sh", "in"),
    ("prioritization_vcf_per_fam.sh", "acmg_result", "gen"),
    ("prioritization_vcf_per_fam.sh", "acmg_criteria_assign.py", "invoke"),
    ("acmg_criteria_assign.py", "splicing_var_analysis.py", "invoke"),
    ("acmg_criteria_assign.py", "utr_anno_interpret.py", "invoke"),
    ("acmg_criteria_assign.py", "acmg_pvs1_null_variant.py", "invoke"),
    ("acmg_criteria_assign.py", "gene_mechanism_common.py", "invoke"),
    # -- caches consumed by the annotation pipeline --
    ("conservation_file", "annotation_vcf.sh", "in"),
    ("alphamissense_prescore", "annotation_vcf.sh", "in"),
    ("spliceai_snv_prescore", "annotation_vcf.sh", "in"),
    ("spliceai_indel_prescore", "annotation_vcf.sh", "in"),
    ("cadd_base_dir", "annotation_vcf.sh", "in"),
    ("cadd_prescore", "annotation_vcf.sh", "in"),
    ("utr_annotator_file", "annotation_vcf.sh", "in"),
    ("loeuf_prescore", "annotation_vcf.sh", "in"),
    ("loftee_repo", "annotation_vcf.sh", "in"),
    ("human_ancestor_fasta", "annotation_vcf.sh", "in"),
    ("loftee_conservation_file", "annotation_vcf.sh", "in"),
    ("splicevault_prescore", "annotation_vcf.sh", "in"),
    ("primateai_prescore", "annotation_vcf.sh", "in"),
    ("gnomad_vcf_chrX", "annotation_vcf.sh", "in"),
    ("clinvar_vcf", "annotation_vcf.sh", "in"),
    ("mavedb_file", "annotation_vcf.sh", "in"),
    ("cds_fasta_file", "annotation_vcf.sh", "in"),
    ("vep_cache", "annotation_vcf.sh", "in"),
    ("vep_plugins", "annotation_vcf.sh", "in"),
    ("vep_plugins_cachedir", "annotation_vcf.sh", "in"),
    ("ref_genome", "annotation_vcf.sh", "in"),
    ("pext_mean_bw", "annotation_vcf.sh", "in"),
    ("pext_tissue_bw_dir", "annotation_vcf.sh", "in"),
    # -- caches consumed directly by the ACMG python scripts (the prioritization
    #    shell only invokes these scripts; the scripts read the caches) --
    ("all_intolerant_domains", "acmg_criteria_assign.py", "in"),
    ("alphamissense_intolerant_domains", "acmg_criteria_assign.py", "in"),
    ("alphamissense_vcf", "acmg_criteria_assign.py", "in"),
    ("alt_disease_vcf", "acmg_criteria_assign.py", "in"),
    ("clingen_map", "acmg_criteria_assign.py", "in"),
    ("clinvar_aa_stat", "acmg_criteria_assign.py", "in"),
    ("clinvar_gene_stat", "acmg_criteria_assign.py", "in"),
    ("clinvar_patho_af_stat", "acmg_criteria_assign.py", "in"),
    ("clinvar_patho_exon_af_stat", "acmg_criteria_assign.py", "in"),
    ("clinvar_splice_stat", "acmg_criteria_assign.py", "in"),
    ("gene_dosage_sensitivity", "acmg_criteria_assign.py", "in"),
    ("gene_nonlof_mechanism_json", "acmg_criteria_assign.py", "in"),
    ("hpo_condition_mechanism_json", "acmg_criteria_assign.py", "in"),
    ("interpro_mapping_pickle", "acmg_criteria_assign.py", "in"),
    ("gene_mechanism_cache_dir", "acmg_criteria_assign.py", "in"),
    ("pm1_regions_pkl", "acmg_pm1_pp2_bp1_missense.py", "in"),
    ("clinvar_patho_af_stat", "acmg_pm2_bs1_ba1_frequency.py", "in"),
    ("gene_dosage_sensitivity", "acmg_pvs1_null_variant.py", "in"),
    ("cds_fasta_file", "acmg_pvs1_null_variant.py", "in"),
    ("alphamissense_tranx_domain_map", "splicing_var_analysis.py", "in"),
    ("alphamissense_tranx_domain_map", "utr_anno_interpret.py", "in"),
    ("alphamissense_tranx_domain_map", "acmg_pvs1_null_variant.py", "in"),
    ("gofcards_exact_gof_cache", "gene_mechanism_common.py", "in"),
    ("gofcards_exact_gof_cache", "acmg_pvs1_null_variant.py", "in"),
    ("mondo_disease_scope_registry", "gene_mechanism_common.py", "in"),
    ("hgnc_table", "gene_mechanism_common.py", "in"),
    ("hpo_assertions", "gene_mechanism_common.py", "in"),
    ("loeuf_prescore", "gene_mechanism_common.py", "in"),
]


def esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
             .replace('"', "&quot;"))


def node_kind(node_id):
    if node_id in DIRS:
        return "dir"
    if node_id in DATA:
        return "data"
    if node_id in RUNTIME:
        return "runtime"
    if node_id in FUNC:
        return FUNC[node_id][2]  # install / prod
    raise KeyError(node_id)


def node_label(node_id):
    if node_id in DATA:
        return DATA[node_id][0]
    if node_id in RUNTIME:
        return RUNTIME[node_id][0]
    return FUNC[node_id][0]


def write_graphml():
    lines = ['<?xml version="1.0" encoding="UTF-8"?>',
             '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">',
             '  <desc>PriVA cache dependency graph. '
             'NODE kinds: data=cache file, dir=user-specified directory, '
             'runtime=cohort input/output, install=install-time function '
             '(shell or python that builds a cache), prod=runtime function '
             '(annotation/filter/ACMG). '
             'EDGE kinds: generates=function->data, input=data->function, '
             'invokes=function->function (parent calls child).</desc>',
             '  <key id="d_kind"   for="node" attr.name="kind"   attr.type="string"/>',
             '  <key id="d_script" for="node" attr.name="script" attr.type="string"/>',
             '  <key id="d_label"  for="node" attr.name="label"  attr.type="string"/>',
             '  <key id="d_desc"   for="node" attr.name="desc"   attr.type="string"/>',
             '  <key id="e_label"  for="edge" attr.name="label"  attr.type="string"/>',
             '  <graph id="G" edgedefault="directed">']
    for nid in sorted(DATA):
        k = "dir" if nid in DIRS else "data"
        lines.append(f'    <node id="data:{nid}"><data key="d_kind">{k}</data>'
                     f'<data key="d_label">{esc(DATA[nid][0])}</data>'
                     f'<data key="d_desc">{esc(DATA[nid][1])}</data></node>')
    for nid in sorted(RUNTIME):
        lines.append(f'    <node id="runtime:{nid}"><data key="d_kind">runtime</data>'
                     f'<data key="d_label">{esc(RUNTIME[nid][0])}</data>'
                     f'<data key="d_desc">{esc(RUNTIME[nid][1])}</data></node>')
    for nid in sorted(FUNC):
        lines.append(f'    <node id="func:{nid}"><data key="d_kind">{FUNC[nid][2]}</data>'
                     f'<data key="d_label">{esc(FUNC[nid][0])}</data>'
                     f'<data key="d_script">{esc(FUNC[nid][1])}</data></node>')
    for src, dst, kind in EDGES:
        sid = ("func:" if src in FUNC else "runtime:" if src in RUNTIME else "data:") + src
        did = ("func:" if dst in FUNC else "runtime:" if dst in RUNTIME else "data:") + dst
        lab = {"gen": "generates", "in": "input", "invoke": "invokes"}[kind]
        lines.append(f'    <edge source="{sid}" target="{did}"><data key="e_label">{lab}</data></edge>')
    lines.append('  </graph>')
    lines.append('</graphml>')
    with open(OUT_GRAPHML, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {OUT_GRAPHML} ({len(DATA)} data, {len(RUNTIME)} runtime, "
          f"{len(FUNC)} func, {len(EDGES)} edges)")


def cluster_of(node_id):
    am = {"alphamissense_prescore","alphamissense_vcf","alphamissense_vep_vcf",
          "alphamissense_intolerant_motifs","alphamissense_pd_stat",
          "alphamissense_tranx_domain_map","alphamissense_intolerant_domains",
          "all_intolerant_domains","pm1_regions_pkl","hcseeker_spots_tsv","rmc_tsv",
          "AlphaMissense_install","AlphaMissense_anno","AlphaMissense_anno_gnomAD",
          "AlphaMissense_pick_missense_constrained_segments","AlphaMissense_stat",
          "AlphaMissense_pick_intolerant_domains","prepare_pm1_regions_fn",
          "alphmis_tsv2vcf.py","am_pick_missense_constrained_segments.py",
          "stat_protein_domain_amscores.py","am_pick_intolerant_domains.py",
          "prepare_pm1_regions.py"}
    if node_id in am:
        return "AlphaMissense / DAS / PM1"
    clinvar = {"clinvar_vcf","clinvar_gene_stat","clinvar_aa_stat","clinvar_splice_stat",
               "clinvar_patho_af_stat","clinvar_patho_exon_af_stat","ClinVar_VCF_deploy",
               "ClinVar_Gene_stat","ClinVar_AA_stat","ClinVar_patho_AF_stat",
               "clinvar_stat_variants.py","stat_aachange_clinvar.py",
               "stat_gene_patho_afs.py","stat_exon_patho_afs.py"}
    if node_id in clinvar:
        return "ClinVar"
    mech = {"gene_mechanism_cache_dir","gene_mechanism_raw_dir","ddg2p_mechanism_evidence",
            "hgnc_table","gene_nonlof_mechanism_json","hpo_condition_mechanism_json",
            "gene_nonlof_mechanism_schema","gofcards_exact_gof_cache","clinvar_vcv_xml",
            "gofcards_mechanism_review_tsv","gofcards_hg38_fasta","gofcards_hg19_to_hg38_chain",
            "hpo_assertions","mondo_obo","mondo_disease_scope_registry",
            "mondo_disease_scope_overrides","gene_pathogenic_mechanism_cache_install",
            "gene_nonlof_mechanism_cache_install","hpo_condition_mechanism_cache_install",
            "gofcards_exact_gof_cache_install","gofcards_clinvar_injection_install",
            "mondo_hpo_scope_install","mechanism_resource_install",
            "build_gene_pathogenic_mechanism_cache.py","build_gene_nonlof_mechanism_cache.py",
            "build_hpo_condition_mechanism_cache.py","build_gofcards_exact_gof_cache.py",
            "inject_clinvar_into_gofcards.py","collapse_HPO_anno.py",
            "build_mondo_disease_scope.py","gene_mechanism_common.py"}
    if node_id in mech:
        return "Mechanism / GoFCards / HPO / MONDO"
    cons = {"acmg_criteria_assign.py","splicing_var_analysis.py","utr_anno_interpret.py",
            "acmg_pvs1_null_variant.py","acmg_pm2_bs1_ba1_frequency.py",
            "acmg_pm1_pp2_bp1_missense.py","annotation_vcf.sh",
            "filtration_vcf_per_fam.sh","prioritization_vcf_per_fam.sh",
            "input_vcf","ped_file","alt_disease_vcf","annotated_vcf","filtered_vcf","acmg_result"}
    if node_id in cons:
        return "Runtime pipeline"
    plugin = {"utr_annotator_file","loeuf_prescore","loftee_repo","human_ancestor_fasta",
              "loftee_conservation_file","spliceai_snv_prescore","spliceai_indel_prescore",
              "splicevault_prescore","primateai_prescore","conservation_file","cadd_base_dir",
              "mavedb_file","pext_mean_bw","pext_tissue_bw_dir","cds_fasta_file",
              "interpro_mapping_pickle","clingen_map","gene_dosage_sensitivity",
              "LoFtee_install","UTRAnnotator_install","LOEUF_install","SpliceAI_install",
              "SpliceVault_install","convert_splicevault_hg38_to_hg19","PrimateAI_install",
              "Conservation_install","CADD_install","MaveDB_install","MaveDB_deployment",
              "pext_install","InterPro_parsing","ClinGen_deploy","CDS_FASTA_install"}
    if node_id in plugin:
        return "VEP plugins / CADD / ClinGen"
    if node_id in {"gnomad_vcf_chrX","gnomAD_install","gnomAD_liftover","hg38_hg19_chain"}:
        return "gnomAD"
    return "Genome / VEP / env"


def write_dot():
    clusters = {}
    for nid in DATA:
        clusters.setdefault(cluster_of(nid), []).append(("data", nid))
    for nid in RUNTIME:
        clusters.setdefault(cluster_of(nid), []).append(("runtime", nid))
    for nid in FUNC:
        clusters.setdefault(cluster_of(nid), []).append((FUNC[nid][2], nid))

    lines = ["digraph PriVA {", '  rankdir=LR;', '  compound=true;',
             '  node [fontname="Helvetica", fontsize=9];',
             '  edge [fontname="Helvetica", fontsize=7, color="#999999"];']
    # legend
    lines += [
        '  subgraph cluster_legend {',
        '    label="Legend"; fontsize=14; style="rounded,filled"; fillcolor="#f8fafc"; color="#94a3b8";',
        '    node [fontsize=11];',
        '    lg_data [shape=box, style=filled, fillcolor="#dbeafe", color="#3b82f6", label="data  = cache file"];',
        '    lg_dir [shape=folder, style=filled, fillcolor="#e0f2fe", color="#0ea5e9", label="dir  = user-specified directory"];',
        '    lg_rt [shape=box, style=filled, fillcolor="#cffafe", color="#0891b2", label="runtime  = cohort input / output"];',
        '    lg_inst [shape=ellipse, style=filled, fillcolor="#fed7aa", color="#ea580c", label="install  = install-time function (builds a cache)"];',
        '    lg_prod [shape=hexagon, style=filled, fillcolor="#fbbf24", color="#92400e", label="prod  = runtime function (annotation / ACMG)"];',
        '    lg_a [style=invis, width=0, height=0, label=""];',
        '    lg_b [style=invis, width=0, height=0, label=""];',
        '    lg_c [style=invis, width=0, height=0, label=""];',
        '    lg_d [style=invis, width=0, height=0, label=""];',
        '    lg_e [style=invis, width=0, height=0, label=""];',
        '    lg_f [style=invis, width=0, height=0, label=""];',
        '    lg_a -> lg_b [color="#16a34a", label="generates  (function -> file)", fontsize=11];',
        '    lg_c -> lg_d [color="#6366f1", label="input  (file -> function)", fontsize=11];',
        '    lg_e -> lg_f [color="#d97706", label="invokes  (parent -> child function)", fontsize=11];',
        '  }',
    ]
    for ci, cname in enumerate(sorted(clusters)):
        lines.append(f'  subgraph cluster_{ci} {{')
        lines.append(f'    label="{esc(cname)}"; fontsize=11; style="rounded,filled"; fillcolor="#f1f5f9"; color="#cbd5e1";')
        for kind, nid in clusters[cname]:
            if kind == "data":
                if nid in DIRS:
                    lines.append(f'    "data:{nid}" [shape=folder, style=filled, fillcolor="#e0f2fe", color="#0ea5e9", label="{esc(DATA[nid][0])}"];')
                else:
                    lines.append(f'    "data:{nid}" [shape=box, style=filled, fillcolor="#dbeafe", color="#3b82f6", label="{esc(DATA[nid][0])}"];')
            elif kind == "runtime":
                lines.append(f'    "runtime:{nid}" [shape=box, style=filled, fillcolor="#cffafe", color="#0891b2", label="{esc(RUNTIME[nid][0])}"];')
            elif kind == "install":
                lines.append(f'    "func:{nid}" [shape=ellipse, style=filled, fillcolor="#fed7aa", color="#ea580c", label="{esc(FUNC[nid][0])}"];')
            else:  # prod
                lines.append(f'    "func:{nid}" [shape=hexagon, style=filled, fillcolor="#fbbf24", color="#92400e", label="{esc(FUNC[nid][0])}\\n({esc(FUNC[nid][1])})"];')
        lines.append('  }')
    for src, dst, kind in EDGES:
        sid = ("func:" if src in FUNC else "runtime:" if src in RUNTIME else "data:") + src
        did = ("func:" if dst in FUNC else "runtime:" if dst in RUNTIME else "data:") + dst
        color = {"gen": "#16a34a", "in": "#6366f1", "invoke": "#d97706"}[kind]
        lines.append(f'  "{sid}" -> "{did}" [color="{color}"];')
    lines.append("}")
    with open(OUT_DOT, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {OUT_DOT}")


def check_graph():
    nodes = set(DATA) | set(RUNTIME) | set(FUNC)
    deg = {n: 0 for n in nodes}
    for src, dst, kind in EDGES:
        deg[src] += 1
        deg[dst] += 1
    iso = sorted(n for n, d in deg.items() if d == 0)
    if iso:
        print("WARNING isolated nodes:", iso, file=sys.stderr)
    else:
        print("OK: no isolated nodes")
    # reachability: every cache (data) must reach a prod node
    from collections import defaultdict, deque
    adj = defaultdict(list)
    for src, dst, kind in EDGES:
        adj[src].append(dst)
    prod = {n for n in FUNC if FUNC[n][2] == "prod"}
    dangling = []
    for d in sorted(DATA):
        if d in NON_CACHE:
            continue
        seen = {d}
        q = deque([d])
        reached = False
        while q:
            x = q.popleft()
            if x in prod:
                reached = True
                break
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    q.append(y)
        if not reached:
            dangling.append(d)
    if dangling:
        print("WARNING caches NOT reaching any prod function:", dangling, file=sys.stderr)
    else:
        print("OK: every cache reaches a prod function")


def render_png():
    rc = subprocess.run(["dot", "-Tpng", OUT_DOT, "-o", OUT_PNG],
                        capture_output=True, text=True)
    if rc.returncode != 0:
        print("dot render failed:", rc.stderr, file=sys.stderr)
        return
    print(f"wrote {OUT_PNG}")


if __name__ == "__main__":
    check_graph()
    write_graphml()
    write_dot()
    render_png()
