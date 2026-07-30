#!/usr/bin/env python
"""Assign ACMG/AMP criteria to an annotated variant table.

This is the top-level orchestrator. Every criterion lives in a sibling module
named after the criteria it carries; this file threads the data between them in
the one order their dependencies allow, then hands the result to the scoring
step. Run it as a script -- see the argument parser at the bottom.

WHERE EACH CRITERION LIVES
==========================

    acmg_pvs1_null_variant          PVS1
    acmg_ps1_pm5_known_variant      PS1  PM5
    acmg_ps2_pm6_denovo             PS2  PM6
    acmg_ps3_bs3_pp5_bp6_reported   PS3  BS3  PP5  BP6
    acmg_pm1_pp2_bp1_missense       PM1  PP2  BP1
    acmg_pm2_bs1_ba1_frequency      PM2  BS1  BA1
    acmg_pm4_bp3_protein_length     PM4  BP3
    acmg_pp1_bs4_bp2_pm3_family     PP1  BS4  BP2  PM3
    acmg_pp3_bp4_bp7_insilico       PP3  BP4  BP7
    acmg_bs2_bp5_observation        BS2  BP5

    acmg_variant_mechanism          shared -- exact curated allele match, and
                                    the fifteen mechanism masks
    acmg_consequence                shared -- VEP consequence, NMD, pext
    acmg_inheritance                shared -- inheritance mode and age of onset
    acmg_scoring                    combine, score, rank

WHY THE ORDER BELOW IS NOT ARBITRARY
====================================

Several criteria consume another criterion's output, so the sequence is forced:

    annotate_exact_nonlof_variants   must precede the mechanism hub, because
                                     the hub's precedence ladder needs to know
                                     whether a curated allele match exists
    annotate_gene_mechanism_categories  must precede every mechanism-reading
                                     criterion; it writes the fifteen columns
                                     they all require, and a missing column is
                                     an error rather than a silent default
    PVS1                             produces loc_intol_domain, which PM1 and
                                     PM4/BP3 both read, and pvs1_criteria,
                                     which PS1/PM5, PM1, PM4/BP3, PP3/BP4 and
                                     BP2/PM3 all read
    PS3/BS3                          produces clinvar_patho and clinvar_benign,
                                     which PS1/PM5, PM1 and PP5/BP6 read
    PM2                              produces pm2_criteria, which BS1, BS2 and
                                     BP2/PM3 read
    identify_inheritance_mode        produces the six inheritance arrays that
                                     PP1, BS1, BS2, BS4, BP2/PM3 and BP5 read

THE CLINVAR CACHES NO LONGER DEPEND ON WHO IS __main__
======================================================

Worth knowing, because it constrained this file until 2026-07-30 and the
constraint is now gone.

``pickle`` does not store a function; it stores a module name and an attribute
name, and looks the pair up on load. The ClinVar amino-acid and splice caches
held nested ``defaultdict`` objects, and because their builder is normally run
rather than imported, the byte stream recorded the factory as
``__main__.nested_defaultdict``. Loading them evaluated
``getattr(sys.modules["__main__"], "nested_defaultdict")`` -- asking whichever
file was the ENTRY POINT, not the file calling ``pickle.load``. So the caches
read correctly only from a program that defined that name itself:

    python acmg_criteria_assign.py ...                      worked
    python my_benchmark.py                                  AttributeError
      from acmg_criteria_assign import ACMG_criteria_assign

The four deployed caches were re-serialized as plain dictionaries, which
removes the reference. Their contents were verified unchanged: an
order-independent digest of every key and value matched before and after, for
all four files. See ``to_plain_dict`` in scripts/stat_aachange_clinvar.py,
which keeps future rebuilds clean.
"""

import logging
import pickle
import gzip
import mmap
import gc
import os
import argparse as ap
from typing import Tuple
import pandas as pd
import numpy as np
from gene_mechanism_hub import (
    DEFAULT_DDG2P_MECHANISM_EVIDENCE,
    DEFAULT_HPO_CONDITION_MECHANISM_CACHE,
    DEFAULT_MECHANISM_JSON,
    annotate_gene_mechanism_categories,
)

from acmg_variant_mechanism import annotate_exact_nonlof_variants
from acmg_consequence import vep_consq_interpret
from acmg_inheritance import identify_inheritance_mode
from acmg_pvs1_null_variant import (
    PVS1_criteria,
    summarize_clinvar_gene_pathogenicity,
)
from acmg_ps1_pm5_known_variant import PS1_PM5_criteria
from acmg_ps2_pm6_denovo import PS2_PM6_criteria
from acmg_ps3_bs3_pp5_bp6_reported import PP5_BP6_criteria, PS3_BS3_criteria
from acmg_pm1_pp2_bp1_missense import PM1_criteria, PP2_BP1_criteria
from acmg_pm2_bs1_ba1_frequency import BA1_criteria, BS1_criteria, PM2_criteria
from acmg_pm4_bp3_protein_length import PM4_BP3_criteria
from acmg_pp1_bs4_bp2_pm3_family import (
    BP2_PM3_criteria,
    BS4_criteria,
    PP1_criteria,
)
from acmg_pp3_bp4_bp7_insilico import BP7_criteria, PP3_BP4_criteria
from acmg_bs2_bp5_observation import BP5_criteria, BS2_criteria
from acmg_scoring import (
    DEFAULT_DISPENSABLE_GENE_LIST,
    calculate_posterior_probability,
    sort_and_rank_variants,
    summarize_acmg_criteria,
)

# The handler goes on the ROOT logger, not on this module's own logger.
# Every criteria module logs to logging.getLogger(__name__) without a handler
# of its own, so their records reach the console only by propagating to root.
# Attaching the handler here instead would silence all of them -- 111 of the
# pipeline's 152 log calls live in those modules.
logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s",
)
logger = logging.getLogger(__name__)


def ACMG_criteria_assign(anno_table: str,
                         am_score_table: str,
                         clinvar_patho_af_stat: str,
                         clinvar_patho_exon_af_stat: str,
                         clinvar_aa_dict_pkl: str,
                         clinvar_splice_dict_pkl: str,
                         interpro_entry_map_pkl: str,
                         intolerant_domains_pkl: str,
                         am_intol_domains_tsv: str,
                         pm1_regions_pkl: str,
                         clinvar_gene_stat_pkl: str,
                         tranx_exon_domain_map_pkl: str,
                         repeat_region_file: str,
                         gene_dosage_sensitivity: str,
                         pp1_vcf: str = "",
                         pp1_ped: str = "",
                         mavedb_metadata_tsv: str = "",
                         fam_name: str = "",
                         ped_table: str = "",
                         alt_disease_vcf: str = "",
                         clingen_map_pkl: str = "",
                         apply_clingen_override: bool = True,
                         cds_fasta_path: str = "",
                         pext_tissues: str = "",
                         pext_low_expression_cutoff: float = 0.1,
                         pext_penalty_floor: float = 0.8,
                         pext_penalty_shape: float = 0.5,
                         relevant_gene_list: str = None,
                         # The same constant the CLI defaults to, so calling this
                         # function directly ranks variants identically to running
                         # the script. When this defaulted to None it was still
                         # forwarded to sort_and_rank_variants unconditionally,
                         # which overrode that function's own default and silently
                         # skipped the dispensable-gene penalty. Pass None here to
                         # mean "deliberately no list".
                         dispensable_gene_list: str = DEFAULT_DISPENSABLE_GENE_LIST,
                         hpo_condition_mechanism_json: str = str(DEFAULT_HPO_CONDITION_MECHANISM_CACHE),
                         gene_mechanism_json: str = str(DEFAULT_MECHANISM_JSON),
                         ddg2p_mechanism_evidence: str = str(DEFAULT_DDG2P_MECHANISM_EVIDENCE),
                         gnomAD_extreme_rare_threshold: float = 0.0001,
                         expected_incidence: float = 0.001,
                         threads: int = 10) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Main function to assign ACMG criteria.
    Returns both annotated DataFrame and criteria matrix.

    The ClinGen SVI group has introduced a new dimension of evaluation for each criteria called Strength of Criteria.
    But many strengths cannot be applied at bioinformatic level.

    ===========Regarding the combining rule===========
    Regarding the combining rule, SVI from ClinGen suggests that PVS1 + 1PP = Likely Pathogenic.
    We just adopted the naive Bayesian Framework to calculate the posterior probability of pathogenicity.

    ===========Regarding the criteria, Cannot be applied===========
    1. PS2/PM6, regarding the DeNovo related evidence. The strength is related to the phenotype consistency with the disease. Such clinical information is rarely normalized for automatic interpretation. Thus Skip.
    2. PP4, Patient's phenotype or family history is highly specific for a disease with a single genetic etiology. (This cannot be applied simutaneously with PP1, because if a disease is tightly linked with only one gene, the segregation is doomed. Therefore further segregation observation does not add any more confidence because the confidence from this perspective is already reaching a ceiling.)

    ===========Regarding the criteria, Can be applied===========
    1. PVS1, LOEUF <= 0.35 should be considered as intolerant to LoF.
    2. PM2, is reduced from Moderate to Supporting.
    3. PP5/BP6, reputable source reported as benign or pathogenic. Suggested to be abandoned or at least not assigned along with PS3/BS3.
    3. PM3, can be only applied if PM2 is True (sufficiently rare in gnomAD).
    """
    anno_df = pd.read_table(anno_table, low_memory=False).drop_duplicates()
    logger.info(f"Got {threads} threads to process the input table {anno_table}, now the table looks like: \n{anno_df[:5].to_string(index=False)}")

    # Log pext tissue information if provided
    if pext_tissues:
        pext_tissue_list = [t.strip() for t in pext_tissues.split(",") if t.strip()]
        logger.info(f"pext tissues configured: {pext_tissue_list}")
        # Expected pext column names: PEXT_MEAN plus PEXT_{tissue} for each tissue
        expected_pext_cols = ["PEXT_MEAN"] + [f"PEXT_{t.replace('-', '_')}" for t in pext_tissue_list]
        present_pext_cols = [c for c in expected_pext_cols if c in anno_df.columns]
        if present_pext_cols:
            logger.info(f"Found pext columns in annotation table: {present_pext_cols}")
        else:
            logger.info(f"No pext columns found in annotation table (expected: {expected_pext_cols})")

    anno_df = vep_consq_interpret(anno_df, threads)

    am_score_df = pd.read_table(am_score_table, low_memory=False)
    ped_df = None
    fam_df = None
    proband = None
    if ped_table:
        ped_df = pd.read_table(
            ped_table, dtype=str, keep_default_na=False, low_memory=False
        )
        if '#FamilyID' not in ped_df.columns:
            ped_df = pd.read_table(
                ped_table,
                header=None,
                names=[
                    '#FamilyID',
                    'IndividualID',
                    'PaternalID',
                    'MaternalID',
                    'Sex',
                    'Phenotype',
                ],
                usecols=range(6),
                dtype=str,
                keep_default_na=False,
                low_memory=False,
            )
        fam_df = ped_df[ped_df['#FamilyID'] == fam_name]
        affected = fam_df.loc[fam_df['Phenotype'] == "2", 'IndividualID']
        if affected.empty:
            raise ValueError(f"No affected individual found for family {fam_name}")
        proband = affected.values[0]


    # Convert the am_score_df to a dictionary:
    # 1. Ensembl transcript ID (column 'transcript') to mean AM score (column 'mean_am_pathogenicity')
    am_score_dict = dict(zip(am_score_df['transcript'], am_score_df['mean_am_pathogenicity']))
    # Use anno_df to create a dict map from Ensembl transcript ID to gene ID
    transcript_to_gene_map = dict(zip(anno_df['Feature'], anno_df['Gene']))
    # Use the two dict above to create dict that maps gene ID to mean AM score
    gene_to_am_score_map = {g: am_score_dict[t] for t, g in transcript_to_gene_map.items() if t in am_score_dict}
    # Both ClinVar structures are loaded once here and handed to every consumer.
    # PS1_PM5_criteria used to load them again from disk itself, which read
    # 44 MiB and 211 MiB twice per run.
    logger.info(f"Loading ClinVar AA change dict from {clinvar_aa_dict_pkl}")
    clinvar_aa_dict = pickle.load(gzip.open(clinvar_aa_dict_pkl)) if clinvar_aa_dict_pkl.endswith(".gz") else pickle.load(open(clinvar_aa_dict_pkl, "rb"))
    logger.info(f"Loading ClinVar splice dict from {clinvar_splice_dict_pkl}")
    clinvar_splice_dict = pickle.load(gzip.open(clinvar_splice_dict_pkl)) if clinvar_splice_dict_pkl.endswith(".gz") else pickle.load(open(clinvar_splice_dict_pkl, "rb"))

    logger.info(f"gene_to_am_score_map created, {len(gene_to_am_score_map)} genes are having the AM score")
    # Both structures are keyed by transcript, so the transcript-to-gene map is
    # what turns them into a gene-level answer. Reading only the amino-acid one
    # would miss any gene whose two-star pathogenic variants are all splice-site
    # or intronic, since those carry no HGVSp and cannot appear there.
    clinvar_pathogenic_genes = summarize_clinvar_gene_pathogenicity(
        transcript_to_gene_map,
        clinvar_aa_dict=clinvar_aa_dict,
        clinvar_splice_dict=clinvar_splice_dict,
    )
    anno_df["Gene_avg_AM_score"] = anno_df["Gene"].map(gene_to_am_score_map)

    # Establish the variant ID column
    anno_df["variant_id"] = anno_df["chrom"] + ":" + anno_df["pos"].astype(str) + ":" + anno_df["ref"] + "-" + anno_df["alt"]
    annotate_exact_nonlof_variants(
        anno_df,
        context="step3_initial",
        mechanism_json=gene_mechanism_json,
    )

    # Load the intolerant domains
    if intolerant_domains_pkl.endswith(".gz"):
        with gzip.open(intolerant_domains_pkl, 'rb') as f:
            mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            intolerant_domains = pickle.load(mm)
    else:
        with open(intolerant_domains_pkl, 'rb') as f:
            mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            intolerant_domains = pickle.load(mm)
    logger.info(f"Loading the recorded intolerant domains which look alike: {intolerant_domains}")

    # These inheritance arrays remain inputs to criteria such as PM3 and PP1.
    # They are not mechanism fallbacks for the hub, PVS1, BS1, BS2, BS4, or
    # ranking.
    recessive, dominant, non_monogenic, non_mendelian, haplo_insufficient, incomplete_penetrance = identify_inheritance_mode(
        anno_df,
        gene_to_am_score_map,
        gene_dosage_sensitivity,
        threads,
    )

    anno_df = annotate_gene_mechanism_categories(
        anno_df,
        condition_cache=hpo_condition_mechanism_json,
        symbol_col="SYMBOL",
        use_hgnc_package=False,
    )
    logger.info(
        "Variant-level mechanism applicability applied before PVS1; effects: \n%s",
        anno_df["variant_effect"].value_counts(dropna=False).to_string(),
    )
    gc.collect()

    # Apply the PVS1 criteria, LoF on a gene known to to be pathogenic due to LoF
    pvs1_criteria, locate_intol_domains = PVS1_criteria(anno_df,
                                                        clinvar_patho_exon_af_stat,
                                                        interpro_entry_map_pkl,
                                                        clinvar_pathogenic_genes=clinvar_pathogenic_genes,
                                                        gene_dosage_sensitivity=gene_dosage_sensitivity,
                                                        intolerant_domains=intolerant_domains,
                                                        tranx_exon_domain_map_pkl=tranx_exon_domain_map_pkl,
                                                        proband_gt_col=proband,
                                                        cds_fasta_path=cds_fasta_path ) # When test on ClinVar variants, fam_name is set to None because no genotype information are provided
    anno_df["span_functional_domains"] = locate_intol_domains
    logger.info(f"PVS1 criteria applied, {(pvs1_criteria > 0).sum()} variants are having the PVS1 criteria")
    gc.collect()

    # Apply the PS3 and BS3 criteria
    ps3bs3_results, anno_df = PS3_BS3_criteria(anno_df, mavedb_metadata_tsv)
    ps3_criteria = ps3bs3_results['PS3']
    bs3_criteria = ps3bs3_results['BS3']
    logger.info(f"PS3 criteria applied, {(ps3_criteria > 0).sum()} variants are having the PS3 criteria")
    logger.info(f"BS3 criteria applied, {(bs3_criteria > 0).sum()} variants are having the BS3 criteria")
    gc.collect()

    # Apply the PS1 and PM5 criteria
    ps1_criteria, pm5_criteria = PS1_PM5_criteria(anno_df, clinvar_aa_dict, clinvar_splice_dict, ps3bs3_results['clinvar_patho'], pvs1_criteria, threads)
    logger.info(f"PS1 criteria applied, {(ps1_criteria > 0).sum()} variants are having the PS1 criteria")
    logger.info(f"PM5 criteria applied, {(pm5_criteria > 0).sum()} variants are having the PM5 criteria")
    gc.collect()

    # Apply the PS2 and PM6 criteria
    if not ped_df is None and not fam_name is None:
        ps2_criteria, pm6_criteria = PS2_PM6_criteria(anno_df, ped_df, fam_name, threads=threads)
    else:
        logger.warning(f"No ped_table provided, skip the PS2 and PM6 criteria")
        ps2_criteria, pm6_criteria = np.array([0] * len(anno_df)), np.array([0] * len(anno_df))
    logger.info(f"PS2 criteria applied, {(ps2_criteria > 0).sum()} variants are having the PS2 criteria")
    logger.info(f"PM6 criteria applied, {(pm6_criteria > 0).sum()} variants are having the PM6 criteria")
    gc.collect()

    '''
    PS4 cannot be applied because usually we dont have enough cases to determine the frequency of the variant
    '''
    # Apply PM1 criteria, mutational hotspot or well-established functional protein domain
    pm1_criteria, loc_intol_domain = PM1_criteria(anno_df, pvs1_criteria, locate_intol_domains, ps3bs3_results['clinvar_patho'], pm1_regions_pkl, threads)
    # pm1_criteria = pm1_criteria & ~ps1_criteria  # PS1 is already a strength to indicate the intolerance to the AA changes incurred by the variant
    logger.info(f"PM1 criteria applied, {(pm1_criteria > 0).sum()} variants are having the PM1 criteria")
    gc.collect()

    # Apply PM4 and BP3 criteria together (mutually exclusive)
    # PM4: in-frame indels NOT in repetitive region OR in functional domain (or stop-loss)
    # BP3: in-frame indels IN repetitive region WITHOUT functional domain
    pm4_criteria, bp3_criteria, in_repeat_vars = PM4_BP3_criteria(anno_df, pvs1_criteria, repeat_region_file, loc_intol_domain)
    anno_df["variant_in_repeat_region"] = in_repeat_vars
    gc.collect()
    logger.info(f"PM4 criteria applied, {(pm4_criteria > 0).sum()} variants are having the PM4 criteria")
    logger.info(f"BP3 criteria applied, {(bp3_criteria > 0).sum()} variants are having the BP3 criteria")

    # Apply PP2 criteria, missense variant in a gene/domain that not only intolerant to truncating variants but also intolerant to missense variants
    clinvar_gene_stat = pickle.load(gzip.open(clinvar_gene_stat_pkl, "rb")) if clinvar_gene_stat_pkl.endswith(".gz") else pickle.load(open(clinvar_gene_stat_pkl, "rb"))
    pp2_criteria, bp1_criteria = PP2_BP1_criteria(anno_df, clinvar_gene_stat, am_intol_domains_tsv)
    logger.info(f"PP2 criteria applied, {(pp2_criteria > 0).sum()} variants are having the PP2 criteria")
    logger.info(f"BP1 criteria applied, {(bp1_criteria > 0).sum()} variants are having the BP1 criteria")
    gc.collect()
    # Apply PP3 criteria, predicted to be deleterious by in-silico tools
    # Pass pvs1_criteria to prevent double-counting per ClinGen SVI guidelines
    pp3_criteria, bp4_criteria = PP3_BP4_criteria(anno_df, pvs1_criteria)
    # bp4_criteria = bp4_criteria & ~bs3_criteria
    logger.info(f"BP4 criteria applied, {(bp4_criteria > 0).sum()} variants are having the BP4 criteria")
    logger.info(f"PP3 criteria applied, {(pp3_criteria > 0).sum()} variants are having the PP3 criteria")
    gc.collect()
    '''
    PP4 cannot be applied, Patient's phenotype or family history is highly specific for a disease with a single genetic etiology
    '''
    # Apply PP5 criteria, reported as pathogenic by a reputable source but without to many supporting evidences
    # Apply BP6 criteria, reported as benign by a reputable source but without to many supporting evidences
    pp5_criteria, bp6_criteria = PP5_BP6_criteria(anno_df, ps3bs3_results['clinvar_patho'], ps3bs3_results['clinvar_benign'])
    logger.info(f"PP5 criteria applied, {(pp5_criteria > 0).sum()} variants are having the PP5 criteria")
    logger.info(f"BP6 criteria applied, {(bp6_criteria > 0).sum()} variants are having the BP6 criteria")
    gc.collect()

    ba1_criteria = BA1_criteria(anno_df)

    # Apply PM2 criteria, absent from gnomAD or extremely rare in gnomAD
    pm2_criteria = PM2_criteria(anno_df,
                                clinvar_patho_af_stat,
                                clinvar_patho_exon_af_stat,
                                gnomAD_extreme_rare_threshold)
    logger.info(f"PM2 criteria applied, {(pm2_criteria > 0).sum()} variants are having the PM2 criteria")
    gc.collect()

    # Apply BS1, PAF of variant is greater than expected incidence of the disease.
    bs1_criteria, bs1_frequency_components = BS1_criteria(
        anno_df,
        expected_incidence=expected_incidence,
        pm2_criteria=pm2_criteria,
        non_monogenic=non_monogenic,
        non_mendelian=non_mendelian,
        incomplete_penetrance=incomplete_penetrance,
        return_frequency_components=True,
    )
    logger.info(f"BS1 criteria applied, {(bs1_criteria > 0).sum()} variants are having the BS1 criteria")

    disease_incidence_bs1 = bs1_frequency_components["disease_incidence_bs1"]
    gene_pathogenic_af_bs1 = bs1_frequency_components["gene_pathogenic_af_bs1"]
    gene_pathogenic_af_blocked_by_pm2 = bs1_frequency_components["gene_pathogenic_af_blocked_by_pm2"]

    # PM2 blocks the gene pathogenic-AF comparator, while the disease-incidence
    # BS1 route is stronger and removes either PM2 strength.
    pm2_backdown = (pm2_criteria > 0) & disease_incidence_bs1
    pm2_backdown_supporting = int(((pm2_criteria == 1) & pm2_backdown).sum())
    pm2_backdown_moderate = int(((pm2_criteria == 2) & pm2_backdown).sum())
    pm2_criteria = pm2_criteria.copy()
    pm2_criteria[pm2_backdown] = 0

    anno_df["audit_BS1_disease_incidence_gate"] = disease_incidence_bs1.astype(int)
    anno_df["audit_BS1_gene_pathogenic_AF_gate"] = gene_pathogenic_af_bs1.astype(int)
    anno_df["audit_BS1_gene_pathogenic_AF_blocked_by_PM2"] = gene_pathogenic_af_blocked_by_pm2.astype(int)
    anno_df["audit_PM2_backdown_for_BS1_disease_incidence"] = pm2_backdown.astype(int)

    logger.info(
        "BS1/PM2 route precedence removed "
        f"{pm2_backdown_supporting} PM2_Supporting and "
        f"{pm2_backdown_moderate} PM2_Moderate assignments via the disease-incidence BS1 route; "
        f"PM2 blocked the gene pathogenic-AF BS1 route for {int(gene_pathogenic_af_blocked_by_pm2.sum())} variants"
    )

    pm2_bs1_conflicts = (pm2_criteria > 0) & (bs1_criteria > 0)
    anno_df["frequency_conflict_PM2_BS1"] = pm2_bs1_conflicts.astype(int)
    if pm2_bs1_conflicts.any():
        logger.warning(
            "Unexpected unresolved PM2/BS1 frequency conflicts remain for "
            f"{int(pm2_bs1_conflicts.sum())} variants"
        )
    else:
        logger.info("No PM2/BS1 frequency conflicts detected")
    gc.collect()

    # Summarize the inheritance mode, first prepare a df composed of recessive, dominant, non_monogenic, non_mendelian, haplo_insufficient, incomplete_penetrance
    inheritance_df = pd.DataFrame({
        'recessive': ["recessive" if i else "" for i in recessive],
        'dominant': ["dominant" if i else "" for i in dominant],
        'non_monogenic': ["non_monogenic" if i else "" for i in non_monogenic],
        'non_mendelian': ["non_mendelian" if i else "" for i in non_mendelian],
        'haplo_insufficient': ["haplo_insufficient" if i else "" for i in haplo_insufficient],
        'incomplete_penetrance': ["incomplete_penetrance" if i else "" for i in incomplete_penetrance]
    })
    logger.info(f"The inheritance_df looks like: \n{inheritance_df[:10].to_string(index=False)}")

    # Apply BS2, variant observed in a healthy adult
    bs2_criteria = BS2_criteria(
        anno_df,
        non_monogenic,
        non_mendelian,
        incomplete_penetrance,
        pm2_criteria,
    )
    logger.info(f"BS2 criteria applied, {(bs2_criteria > 0).sum()} variants are having the BS2 criteria")
    gc.collect()

    # Apply PP1 criteria, variant is cosegregating with a pathogenic variant in one or more families
    pp1_criteria = PP1_criteria(anno_df, recessive, dominant, non_monogenic, non_mendelian, incomplete_penetrance, pp1_vcf, pp1_ped)
    logger.info(f"PP1 criteria applied, {(pp1_criteria > 0).sum()} variants are having the PP1 criteria")
    gc.collect()

    # Apply BS4, lack of family segregation
    if not ped_df is None and not fam_name is None:
        bs4_criteria = BS4_criteria(
            anno_df,
            ped_df,
            fam_name,
            non_monogenic,
            non_mendelian,
            incomplete_penetrance,
        )
    else:
        logger.warning(f"No ped_table provided, skip the BS4 criteria")
        bs4_criteria = np.array([0] * len(anno_df))
    logger.info(f"BS4 criteria applied, {(bs4_criteria > 0).sum()} variants are having the BS4 criteria")

    # NOTE: BP3 is now calculated together with PM4 in PM4_BP3_criteria (they are mutually exclusive)
    # BP3 was already assigned above when PM4_BP3_criteria was called

    # Apply BP5, variant found in a sample with known alternative molecular basis for disease
    if alt_disease_vcf:
        bp5_criteria = BP5_criteria(anno_df, alt_disease_vcf, recessive, dominant, non_monogenic, non_mendelian, incomplete_penetrance, threads)
    else:
        bp5_criteria = np.array([0] * len(anno_df))
    logger.info(f"BP5 criteria applied, {(bp5_criteria > 0).sum()} variants are having the BP5 criteria")
    gc.collect()
    # Apply BP7, synonymous variant, no splicing-altering consequence, not conserved.
    bp7_criteria = BP7_criteria(anno_df)
    logger.info(f"BP7 criteria applied, {(bp7_criteria > 0).sum()} variants are having the BP7 criteria")
    gc.collect()

    # Apply BP2, observed in trans with a pathogenic variant in dominant disease, Or in-cis with a pathogenic variant with any inheritance mode
    # Apply PM3, observed in trans with a pathogenic variant in recessive disease.
    if not ped_df is None and not fam_name is None:
        bp2_criteria, pm3_criteria = BP2_PM3_criteria(anno_df,
                                                      ped_df,
                                                      pm2_criteria,
                                                      pvs1_criteria,
                                                      ps1_criteria,
                                                      ps3_criteria,
                                                      recessive,
                                                      dominant,
                                                      incomplete_penetrance)
    else:
        logger.warning(f"No ped_table provided, skip the BP2 and PM3 criteria")
        bp2_criteria, pm3_criteria = np.array([0] * len(anno_df)), np.array([0] * len(anno_df))
    logger.info(f"BP2 criteria applied, {(bp2_criteria > 0).sum()} variants are having the BP2 criteria")
    logger.info(f"PM3 criteria applied, {(pm3_criteria > 0).sum()} variants are having the PM3 criteria")
    gc.collect()

    # Collect all criteria in a dictionary
    criteria_dict = {
        'PVS1': pvs1_criteria,
        'PS1': ps1_criteria,
        'PS2': ps2_criteria,
        'PS3': ps3_criteria,
        'PM1': pm1_criteria,
        'PM2': pm2_criteria,
        'PM3': pm3_criteria,
        'PM4': pm4_criteria,
        'PM5': pm5_criteria,
        'PM6': pm6_criteria,
        'PP1': pp1_criteria,
        'PP2': pp2_criteria,
        'PP3': pp3_criteria,
        'PP5': pp5_criteria,
        'BA1': ba1_criteria,
        'BS1': bs1_criteria,
        'BS2': bs2_criteria,
        'BS3': bs3_criteria,
        'BS4': bs4_criteria,
        'BP1': bp1_criteria,
        'BP2': bp2_criteria,
        'BP3': bp3_criteria,
        'BP4': bp4_criteria,
        'BP5': bp5_criteria,
        'BP6': bp6_criteria,
        'BP7': bp7_criteria
    }

    # Create summary and matrix
    anno_df, criteria_matrix = summarize_acmg_criteria(
        anno_df, criteria_dict, clingen_map_pkl,
        apply_clingen_override=apply_clingen_override,
    )
    for audit_column in (
        "audit_BS1_disease_incidence_gate",
        "audit_BS1_gene_pathogenic_AF_gate",
        "audit_BS1_gene_pathogenic_AF_blocked_by_PM2",
        "audit_PM2_backdown_for_BS1_disease_incidence",
        "frequency_conflict_PM2_BS1",
    ):
        criteria_matrix[audit_column] = anno_df[audit_column].to_numpy()
    # Save the criteria matrix to a file which the path is based on the input anno_table
    output_matrix = ".".join(anno_table.split(".")[:-1]) + ".acmg.tsv"
    criteria_matrix.to_csv(output_matrix, sep="\t", index=False)

    # Apply quantification using ACMG_criteria column and use that to sort the variants
    posterior_probability, acmg_class = zip(*criteria_matrix.apply(calculate_posterior_probability, axis=1))
    anno_df["ACMG_quant_score"] = posterior_probability
    anno_df["ACMG_class"] = acmg_class

    # Sort and rank variants
    anno_df = sort_and_rank_variants(anno_df,
                                     ped_df,
                                     fam_name,
                                     pext_tissues=pext_tissues,
                                     pext_low_expression_cutoff=pext_low_expression_cutoff,
                                     pext_penalty_floor=pext_penalty_floor,
                                     pext_penalty_shape=pext_penalty_shape,
                                     relevant_gene_list=relevant_gene_list,
                                     dispensable_gene_list=dispensable_gene_list)
    # Save the annotated table to replace the input anno_table
    anno_df.to_csv(anno_table, sep="\t", index=False)

    return anno_df, criteria_matrix


if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("--anno_table", type=str, required=True)
    parser.add_argument("--am_score_table", type=str, required=True)
    parser.add_argument("--ped_table", type=str, required=False, default=None)
    parser.add_argument("--fam_name", type=str, required=False, default=None)
    parser.add_argument("--clinvar_patho_af_stat", type=str, required=True)
    parser.add_argument("--clinvar_patho_exon_af_stat", type=str, required=True)
    parser.add_argument("--clinvar_aa_dict_pkl", type=str, required=True)
    parser.add_argument("--clinvar_splice_dict_pkl", type=str, required=True)
    parser.add_argument("--interpro_entry_map_pkl", type=str, required=True)
    parser.add_argument("--intolerant_domains_pkl", type=str, required=True)
    parser.add_argument("--gene_dosage_sensitivity", type=str, required=True)
    parser.add_argument("--am_intol_domains_tsv", type=str, required=True)
    parser.add_argument("--pm1_regions_pkl", type=str, required=True)
    parser.add_argument("--clinvar_gene_stat_pkl", type=str, required=True)
    parser.add_argument("--tranx_exon_domain_map_pkl", type=str, required=True)
    parser.add_argument("--am_score_vcf", type=str, required=False, default=None)
    parser.add_argument("--alt_disease_vcf", type=str, required=False, default=None)
    parser.add_argument("--clingen_map_pkl", type=str, required=False, default=None)
    parser.add_argument("--ignore_clingen_override", action="store_true",
                        help="Do NOT let ClinGen's curated evidence codes replace the "
                             "criteria PriVA derived. Production runs should leave this "
                             "off: an expert panel's verdict outranks PriVA's own "
                             "reasoning about the same variant. Use it only when "
                             "benchmarking PriVA against ClinGen-curated variants, "
                             "where the override would make the comparison circular.")
    parser.add_argument("--mavedb_metadata_tsv", type=str, required=False, default=None)
    parser.add_argument("--cds_fasta_path", type=str, required=False, default=None,
                        help="Path to Ensembl CDS FASTA file for alternative start codon detection. "
                             "Download from: ftp://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/cds/")
    parser.add_argument("--pext_tissues", type=str, required=False, default="",
                        help="Comma-separated list of pext tissue names used in VEP annotation. "
                             "Used to identify pext-related columns in the annotation table.")
    parser.add_argument("--pext_low_expression_cutoff", type=float, required=False, default=0.1,
                        help="PEXT value below which LoF ranking receives a bounded expression penalty.")
    parser.add_argument("--pext_penalty_floor", type=float, required=False, default=0.8,
                        help="Minimum PEXT ranking multiplier for zero expression.")
    parser.add_argument("--pext_penalty_shape", type=float, required=False, default=0.5,
                        help="Power transform shape for sub-cutoff PEXT ranking penalty.")
    parser.add_argument("--repeat_region_file", type=str, required=True)
    parser.add_argument("--gnomAD_extreme_rare_threshold", type=float, required=False, default=0.0001)
    parser.add_argument("--expected_incidence", type=float, required=False, default=0.001)
    parser.add_argument("--relevant_gene_list", type=str, required=False, default=None)
    parser.add_argument("--dispensable_gene_list", type=str, required=False, default=DEFAULT_DISPENSABLE_GENE_LIST)
    parser.add_argument(
        "--hpo_condition_mechanism_json",
        type=str,
        required=False,
        default=str(DEFAULT_HPO_CONDITION_MECHANISM_CACHE),
        help="Integrated gene-condition-inheritance-penetrance-mechanism cache.",
    )
    parser.add_argument(
        "--gene_nonlof_mechanism_json",
        "--gene_mechanism_json",
        dest="gene_mechanism_json",
        type=str,
        required=False,
        default=str(DEFAULT_MECHANISM_JSON),
        help=(
            "Schema-v2 non-LOF mechanism assertions JSON. "
            "--gene_mechanism_json is retained as a compatibility alias."
        ),
    )
    parser.add_argument("--ddg2p_mechanism_evidence", type=str, required=False, default=str(DEFAULT_DDG2P_MECHANISM_EVIDENCE))
    parser.add_argument("--threads", type=int, required=False, default=10)
    args = parser.parse_args()

    anno_df, criteria_matrix = ACMG_criteria_assign(args.anno_table,
                                                    args.am_score_table,
                                                    args.clinvar_patho_af_stat,
                                                    args.clinvar_patho_exon_af_stat,
                                                    args.clinvar_aa_dict_pkl,
                                                    args.clinvar_splice_dict_pkl,
                                                    args.interpro_entry_map_pkl,
                                                    args.intolerant_domains_pkl,
                                                    args.am_intol_domains_tsv,
                                                    args.pm1_regions_pkl,
                                                    args.clinvar_gene_stat_pkl,
                                                    args.tranx_exon_domain_map_pkl,
                                                    args.repeat_region_file,
                                                    args.gene_dosage_sensitivity,
                                                    mavedb_metadata_tsv=args.mavedb_metadata_tsv,
                                                    fam_name=args.fam_name,
                                                    ped_table=args.ped_table,
                                                    alt_disease_vcf=args.alt_disease_vcf,
                                                    clingen_map_pkl=args.clingen_map_pkl,
                                                    apply_clingen_override=not args.ignore_clingen_override,
                                                    cds_fasta_path=args.cds_fasta_path,
                                                    pext_tissues=args.pext_tissues,
                                                    pext_low_expression_cutoff=args.pext_low_expression_cutoff,
                                                    pext_penalty_floor=args.pext_penalty_floor,
                                                    pext_penalty_shape=args.pext_penalty_shape,
                                                    relevant_gene_list=args.relevant_gene_list,
                                                    dispensable_gene_list=args.dispensable_gene_list,
                                                    hpo_condition_mechanism_json=args.hpo_condition_mechanism_json,
                                                    gene_mechanism_json=args.gene_mechanism_json,
                                                    ddg2p_mechanism_evidence=args.ddg2p_mechanism_evidence,
                                                    gnomAD_extreme_rare_threshold=args.gnomAD_extreme_rare_threshold,
                                                    expected_incidence=args.expected_incidence,
                                                    threads=args.threads)
