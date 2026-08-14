#!/usr/bin/env python3
"""PVS1 -- a null variant in a gene where loss of function causes disease.

The single most consequential criterion in the framework: at Very Strong it
can carry a variant to Likely Pathogenic on its own. It therefore asks two
questions and never mixes them.

IS LOSS OF FUNCTION A DISEASE MECHANISM OF THIS GENE?
    A question about the gene, reading none of the variant's mechanism
    scores. Any one of five independent gene-level sources opens the gate:
    a curated loss-of-function entry in the condition cache, a ClinVar
    pathogenic variant at two stars or better, LOEUF below 0.35, mean
    AlphaMissense above 0.6, or a ClinGen haploinsufficiency call.
    summarize_clinvar_gene_pathogenicity supplies the second of those.

HOW CONVINCING IS THIS PARTICULAR NULL VARIANT?
    Decided from the consequence -- whether nonsense-mediated decay is
    triggered, whether an intolerant domain is spanned, what fraction of the
    protein is lost, whether an alternative start codon rescues a lost one.
    span_functional_domains, downstream_domain_impact,
    downstream_exon_patho_af and evaluate_start_lost_with_cds serve this half.

The one mechanism score PVS1 reads is is_exact_gof, far below, which
withdraws PVS1 from an allele curated as gain of function. It takes evidence
away, never grants it.
"""

import logging
import re
import pickle
import gzip
import json
import os
import pandas as pd
import numpy as np
from typing import Tuple

from alternative_start_codon import adjust_pvs1_for_start_lost, get_detector
from protein_domain_mapping import DomainNormalizer

from acmg_consequence import truncate_fraction
from acmg_variant_mechanism import _variant_mechanism_masks


logger = logging.getLogger(__name__)


def summarize_clinvar_gene_pathogenicity(transcript_to_gene_map: dict,
                                         clinvar_aa_dict: dict = None,
                                         clinvar_splice_dict: dict = None,
                                         high_confidence_status = {
        # Higher confidence (2+ stars)
        'practice_guideline': 4,                                   # 4 stars
        'reviewed_by_expert_panel': 3,                             # 3 stars
        'criteria_provided,_multiple_submitters,_no_conflicts': 2,  # 2 stars
    }) -> set:
    '''
    Which genes carry a high-confidence pathogenic variant in ClinVar?

    Returns a set of gene identifiers, using whatever key ``transcript_to_gene_map``
    maps to -- in the PriVA pipeline that is the Ensembl gene ID, matching the
    ``Gene`` column.

    WHY THIS READS TWO CLINVAR STRUCTURES AND NOT ONE
    =================================================

    PriVA's ClinVar build (scripts/stat_aachange_clinvar.py) makes one pass over
    the annotated ClinVar VCF and writes two files, because they are keyed on
    different things:

      clinvar_aa_dict     keyed on the PROTEIN consequence
        {transcript: {protein_pos: {hgvsp: {CLNSIG: [...], CLNREVSTAT: [...]}}}}
        Every entry has a non-blank HGVSp. Frameshift and nonsense variants are
        present (p.*fs, p.*Ter), but a canonical splice-site variant such as
        c.1234+1G>A produces no HGVSp in VEP and therefore cannot appear here at
        all -- nor can deep-intronic variants or whole-exon deletions.

      clinvar_splice_dict keyed on having an EXON or INTRON field
        {transcript: [{clinvar_sig, clinvar_review, exon, intron, hgvsc, ...}]}
        Those variants do appear here.

    Reading only the first silently misses any gene whose two-star pathogenic
    variants all lack a protein consequence. Measured on the hg19 build, the
    second file adds 171 genes (2,869 -> 3,040).

    THE TEST, applied identically to both structures
    ================================================

      CLNSIG split on , / | ; and lowercased, then intersected with
      {pathogenic, likely_pathogenic}
                    AND
      high_confidence_status[CLNREVSTAT] >= 2, i.e. two stars or better:
          practice_guideline                                        4
          reviewed_by_expert_panel                                  3
          criteria_provided,_multiple_submitters,_no_conflicts       2
          anything else                                             0  -> fails

    Splitting CLNSIG into tokens is what keeps Conflicting_classifications_of_
    pathogenicity out: it tokenizes to {conflicting_classifications_of_pathogenicity},
    which does not match. The splice builder's own filter is looser than this
    (it accepts any CLNSIG containing "athogenic"), so the strict test is
    re-applied here rather than trusted from the build. Measured on the hg19
    build, no gene enters through the looser match alone, so the two sources
    agree by construction.

    A transcript absent from transcript_to_gene_map contributes nothing, since
    there is no gene to attribute it to.
    '''

    def qualifies(clnsig, revstat) -> bool:
        sig_tokens = {
            token.strip().lower()
            for token in re.split(r"[,/|;]", str(clnsig))
            if token.strip()
        }
        if not ({"pathogenic", "likely_pathogenic"} & sig_tokens):
            return False
        return high_confidence_status.get(revstat, 0) >= 2

    def aa_transcript_qualifies(positions) -> bool:
        for hgvsp_dict in positions.values():
            for info in hgvsp_dict.values():
                for clnsig, revstat in zip(info.get('CLNSIG', []),
                                           info.get('CLNREVSTAT', [])):
                    if qualifies(clnsig, revstat):
                        return True
        return False

    def splice_transcript_qualifies(records) -> bool:
        for record in records:
            if qualifies(record.get('clinvar_sig'), record.get('clinvar_review')):
                return True
        return False

    pathogenic_genes = set()
    for label, source, transcript_qualifies in (
            ("protein consequence", clinvar_aa_dict or {}, aa_transcript_qualifies),
            ("exon/intron", clinvar_splice_dict or {}, splice_transcript_qualifies)):
        genes_from_source = set()
        for transcript, payload in source.items():
            gene = transcript_to_gene_map.get(transcript)
            # Already established for this gene: existence is all that is asked.
            if not gene or gene in genes_from_source:
                continue
            if transcript_qualifies(payload):
                genes_from_source.add(gene)
        logger.info(
            "ClinVar %s records: %s transcripts scanned, %s distinct genes carry "
            "a pathogenic or likely-pathogenic variant at >=2 stars, of which %s "
            "were not already found by an earlier structure",
            label, len(source), len(genes_from_source),
            len(genes_from_source - pathogenic_genes),
        )
        pathogenic_genes |= genes_from_source

    logger.info(
        "Found %s genes with a high-confidence pathogenic ClinVar variant",
        len(pathogenic_genes),
    )
    return pathogenic_genes


def identify_alternative_start_codon_genes(df: pd.DataFrame) -> set:
    """
    DEPRECATED: This function checks for alternative TRANSCRIPTS, not alternative start codons.

    Per ClinGen SVI guidelines (PMC6185798), for start_lost variants we should check for
    downstream in-frame ATG codons WITHIN THE SAME TRANSCRIPT, not alternative transcripts.

    This function is kept for backwards compatibility but should be replaced with
    CDS-based alternative start codon detection using AlternativeStartCodonDetector.

    What this function checks: Different transcripts with different TSS (transcription start sites)
    What ClinGen SVI requires: Downstream ATG in the SAME transcript (translation initiation)
    """
    if df["Consequence"].str.contains("start_lost").any():
        if df["Consequence"].str.contains("start_lost").all():
            return np.nan
        elif (np.logical_not(df["Consequence"].str.contains("start_lost")) & df["BIOTYPE"].fillna(".").str.contains("protein_coding")).any():
            return df["Gene"].values[0] if len(df["Gene"].values) > 0 else np.nan
        else:
            return np.nan
    else:
        return np.nan


def evaluate_start_lost_with_cds(
    df: pd.DataFrame,
    cds_fasta_path: str,
    intolerant_domains: np.ndarray,
    lof_intolerant: np.ndarray = None
) -> Tuple[np.ndarray, pd.DataFrame]:
    """
    Evaluate start_lost variants using CDS-based alternative start codon detection.

    Per ClinGen SVI guidelines (PMC6185798), start_lost N-terminal loss is evaluated
    similarly to in-frame deletions:

    1. NO alternative downstream ATG:
       - Gene is LoF intolerant → PVS1 (Very Strong)
       - Gene is LoF tolerant → PVS1_Strong

    2. Alternative ATG EXISTS (like in-frame deletion):
       - >10% N-terminal loss + spans critical domain → PVS1_Strong
       - >10% N-terminal loss + no critical domain → PVS1_Moderate
       - ≤10% N-terminal loss + spans critical domain → PVS1_Moderate
       - ≤10% N-terminal loss + no critical domain → PVS1_Supporting

    Args:
        df: DataFrame with variant annotations
        cds_fasta_path: Path to Ensembl CDS FASTA file
        intolerant_domains: Boolean array indicating variants in intolerant domains
        lof_intolerant: Boolean array indicating gene is LoF intolerant (LOEUF < 0.35 or high AM)

    Returns:
        Tuple of:
        - pvs1_start_lost: Array of PVS1 strength values for start_lost variants
        - df: Updated DataFrame with alternative start codon annotations
    """
    pvs1_start_lost = np.zeros(len(df), dtype=int)

    # Identify start_lost variants
    start_lost_mask = df["Consequence"].str.contains("start_lost").fillna(False)
    n_start_lost = start_lost_mask.sum()

    if n_start_lost == 0:
        logger.info("No start_lost variants found")
        return pvs1_start_lost, df

    logger.info(f"Evaluating {n_start_lost} start_lost variants for alternative start codons")

    # Default lof_intolerant to True if not provided (conservative)
    if lof_intolerant is None:
        lof_intolerant = np.ones(len(df), dtype=bool)

    # Check if CDS FASTA file exists
    if not cds_fasta_path or not os.path.exists(cds_fasta_path):
        logger.warning(f"CDS FASTA file not found: {cds_fasta_path}")
        logger.warning("Falling back to default PVS1_Moderate for all start_lost variants")
        pvs1_start_lost[start_lost_mask] = 2  # PVS1_Moderate (default for start_lost)
        return pvs1_start_lost, df

    # Initialize the detector
    try:
        detector = get_detector(cds_fasta_path)
    except Exception as e:
        logger.error(f"Failed to initialize AlternativeStartCodonDetector: {e}")
        pvs1_start_lost[start_lost_mask] = 2
        return pvs1_start_lost, df

    # Get unique transcripts with start_lost
    start_lost_df = df.loc[start_lost_mask, ['Feature', 'Gene']].copy()
    unique_transcripts = start_lost_df['Feature'].dropna().unique()

    # Batch lookup alternative start codons
    alt_start_info = {}
    for tid in unique_transcripts:
        alt_start_info[tid] = detector.find_alternative_start(tid)

    # Add alternative start info to DataFrame
    df['alt_start_exists'] = df['Feature'].map(
        lambda t: alt_start_info.get(t, {}).get('has_alternative', False) if pd.notna(t) else False
    )
    df['alt_start_truncation'] = df['Feature'].map(
        lambda t: alt_start_info.get(t, {}).get('n_terminal_truncation_fraction', None) if pd.notna(t) else None
    )
    df['alt_start_position'] = df['Feature'].map(
        lambda t: alt_start_info.get(t, {}).get('first_alt_position', None) if pd.notna(t) else None
    )

    # Calculate PVS1 strength for each start_lost variant
    for idx in df.index[start_lost_mask]:
        iloc_idx = df.index.get_loc(idx)
        transcript = df.loc[idx, 'Feature']
        info = alt_start_info.get(transcript, {})
        spans_critical = intolerant_domains[iloc_idx] if len(intolerant_domains) > 0 else False
        gene_lof_intol = lof_intolerant[iloc_idx] if len(lof_intolerant) > 0 else True

        pvs1_strength = adjust_pvs1_for_start_lost(info, spans_critical, gene_lof_intol)
        pvs1_start_lost[iloc_idx] = pvs1_strength

    # Log summary
    has_alt = df.loc[start_lost_mask, 'alt_start_exists'].sum()
    logger.info(f"start_lost variants with alternative downstream ATG: {has_alt}/{n_start_lost}")
    logger.info(f"PVS1 strength distribution for start_lost: "
                f"PVS1={np.sum(pvs1_start_lost == 4)}, "
                f"PVS1_Strong={np.sum(pvs1_start_lost == 3)}, "
                f"PVS1_Moderate={np.sum(pvs1_start_lost == 2)}, "
                f"PVS1_Supporting={np.sum(pvs1_start_lost == 1)}, "
                f"Not applicable={np.sum((pvs1_start_lost == 0) & start_lost_mask.values)}")

    return pvs1_start_lost, df


def downstream_domain_impact(protein_pos_str, tranx_id, tranx_pp_domain_map, interpro_entry_map_dict, dm_instance, domains="", intol_domains=[], ensg_id=""):
    '''
    Used for frameshift and stopgain variants to explore whether the downstream protein region involving functional domains.
    Uses protein position to check if any domain starts beyond the truncation point.
    '''
    # First check the variant's own position via VEP DOMAINS field
    if isinstance(domains, str):
        for domain in domains.split("&"):
            if ensg_id + ":" + domain in intol_domains:
                return True
            elif dm_instance.interpret_functionality(domain, interpro_entry_map_dict) == "Functional":
                return True

    # Parse protein position from VEP format "150/500" or "150-151/500"
    truncation_pos = None
    if isinstance(protein_pos_str, str) and protein_pos_str:
        try:
            truncation_pos = int(protein_pos_str.split('/')[0].split('-')[0])
        except ValueError:
            pass

    if truncation_pos is None:
        return False

    tranx_id = tranx_id.split(".", 1)[0]
    if tranx_id not in tranx_pp_domain_map:
        return False

    pp_domains = tranx_pp_domain_map[tranx_id].get('pp_domains', [])
    for aa_start, aa_end, domain_path in pp_domains:
        if aa_start > truncation_pos:
            domain = domain_path.split(":", 1)[1]
            if ensg_id + ":" + domain in intol_domains:
                return True
            elif dm_instance.interpret_functionality(domain, interpro_entry_map_dict) == "Functional":
                return True
    return False


def downstream_exon_patho_af(row, clinvar_patho_exon_af_dict, logic="any", threshold=0.01):
    exon_str = row['EXON']
    affected_exons = set([])
    if not isinstance(exon_str, str):
        return False
    elif "-" in exon_str and "/" in exon_str:
        affected_exons.update(range(int(exon_str.split("-")[0]), int(exon_str.split("/")[1]) + 1))
    elif "/" in exon_str:
        affected_exons.update(range(int(exon_str.split("/")[0]), int(exon_str.split("/")[1]) + 1))
    else:
        raise ValueError(f"Invalid exon string: {exon_str}")

    tranx_id = row['Feature']
    affected_exons_patho_af = set([])
    for exon in affected_exons:
        affected_exons_patho_af.add(clinvar_patho_exon_af_dict.get(tranx_id, {}).get(exon, (0, ))[0])

    if logic == "any":
        return any(float(epa) < threshold for epa in affected_exons_patho_af)
    elif logic == "all":
        return all(float(epa) < threshold for epa in affected_exons_patho_af)
    else:
        raise ValueError(f"Invalid logic: {logic}, it should be either 'any' or 'all', depending on your needs")


def span_functional_domains(row,
                            dm_instance=None,
                            interpro_entry_map_dict=None,
                            tranx_pp_domain_map=None,
                            clinvar_patho_exon_af_dict=None,
                            intol_domains=[],
                            exon_patho_af_threshold=0.01):
    '''
    Identify whether a truncating variant is involving functional regions on a protein
    '''
    ensg_id = str(row['Gene'])
    func_domain = False
    exon_frequent_patho = False
    if ("frameshift" in row['Consequence']) or ("stop_gained" in row['Consequence']):
        # These variants not only affect the local region of proteins, but also affect the downstream protein regions
        func_domain = downstream_domain_impact(row['Protein_position'], row['Feature'], tranx_pp_domain_map, interpro_entry_map_dict, dm_instance, domains=str(row["DOMAINS"]), intol_domains=intol_domains, ensg_id=ensg_id)
        exon_frequent_patho = downstream_exon_patho_af(row, clinvar_patho_exon_af_dict, logic="any", threshold=exon_patho_af_threshold)
        logger.info(f"For variant {row['chrom']}:{row['pos']} with transcript {row['Feature']}, the exon is: {row['EXON']}, and the pathogenic AF is: {exon_frequent_patho}")
    elif row["5UTR_lof"] or row["5UTR_frameshift"] or row["5UTR_len_changing"]:
        func_domain = row["5UTR_span_intol_domain"]
        affected_exons = {"1"}
        affected_exon_patho_af = { exon: clinvar_patho_exon_af_dict.get(row['Feature'], {}).get(exon, (0, ))[0] for exon in affected_exons }
        logger.info(f"For variant {row['chrom']}:{row['pos']} with transcript {row['Feature']}, the affected exons are: {affected_exons}, and the pathogenic AFs are: {affected_exon_patho_af}")
        exon_frequent_patho = any(float(epa) < exon_patho_af_threshold for epa in affected_exon_patho_af.values())
    elif row["splicing_lof"] or row["splicing_frameshift"] or row["splicing_len_changing"]:
        func_domain = row["splicing_span_intol_domain"]
        # Parse affected exons from splicing analysis, remove empty strings
        affected_exons_raw = set(str(row["splicing_affected_exons"]).split(","))
        affected_exons = {e.strip() for e in affected_exons_raw if e and e.strip() and e.strip().isdigit()}

        if affected_exons:
            affected_exon_patho_af = { exon: clinvar_patho_exon_af_dict.get(row['Feature'], {}).get(exon, (0, ))[0] for exon in affected_exons }
            logger.info(f"For variant {row['chrom']}:{row['pos']} with transcript {row['Feature']}, the affected exons are: {affected_exons}, and the pathogenic AFs are: {affected_exon_patho_af}")
            # For splicing-induced exon skipping: if ANY affected exon is important (rare patho variants),
            # treat as significant deletion per ClinGen PVS1 exon deletion logic
            exon_frequent_patho = any(float(epa) < exon_patho_af_threshold for epa in affected_exon_patho_af.values())
        else:
            # No affected exons identified, fall back to variant's exon/intron position
            logger.warning(f"No affected exons parsed for splicing variant at {row['chrom']}:{row['pos']}, using EXON field as fallback")
            exon_frequent_patho = float(clinvar_patho_exon_af_dict.get(row['Feature'], {}).get(row['EXON'], (0, ))[0]) < exon_patho_af_threshold
    else:
        if isinstance(row["DOMAINS"], str):
            if intol_domains:
                func_domain = any([ensg_id + ":" + domain in intol_domains for domain in str(row["DOMAINS"]).split("&")])
            else:
                func_domain = any([dm_instance.interpret_functionality(domain, interpro_entry_map_dict) == "Functional" for domain in str(row["DOMAINS"]).split("&")])
        exon_frequent_patho = float(clinvar_patho_exon_af_dict.get(row['Feature'], {}).get(row['EXON'], (0, ))[0]) < exon_patho_af_threshold
        logger.debug(f"For variant {row['chrom']}:{row['pos']} with transcript {row['Feature']}, the exon is: {row['EXON']}, and the pathogenic AF is: {exon_frequent_patho}")

    return func_domain, exon_frequent_patho


def load_non_lof_gene_set(nonlof_mechanism_json: str) -> set:
    """Load the gene symbols whose curated pathogenic mechanism is non-LoF only.

    Reads the non-LoF mechanism cache (gene_nonlof_mechanism_curated_assertions.json.gz,
    built by build_gene_nonlof_mechanism_cache.py) and returns the set of gene symbols
    it names. These are genes whose disease mechanism is gain-of-function, dominant
    negative, or triplosensitivity -- loss of function is NOT the curated mechanism, so
    PVS1 should be withheld. The list overrides the statistical LoF proxies (LOEUF /
    mean AlphaMissense) and the mechanism-agnostic ClinVar 2-star pathogenic signal;
    only an established ClinGen haploinsufficiency call (HI score 3) overrides it.
    """
    with gzip.open(nonlof_mechanism_json, "rb") as f:
        data = json.load(f)
    return {entry.get("symbol") for entry in data.values() if entry.get("symbol")}


def PVS1_criteria(df: pd.DataFrame,
                  clinvar_patho_exon_af_stat: str,
                  interpro_entry_map_pkl: str,
                  clinvar_pathogenic_genes: set = None,
                  gene_dosage_sensitivity: str = None,
                  nonlof_mechanism_json: str = None,
                  intolerant_domains: set = [],
                  tranx_exon_domain_map_pkl: str = None,
                  proband_gt_col: str = None,
                  cds_fasta_path: str = None) -> pd.DataFrame:
    '''
    Introducing varied strength levels of PVS1,
    PVS1: 4
    PVS1_strong: 3
    PVS1_moderate: 2
    PVS1_supporting: 1
    not_applicable: 0
    '''
    # NOTE:
    # PVS1 is a *variant-level* evidence category that depends on whether LoF is an established
    # disease mechanism for the gene/disease context, not on the proband's zygosity per se.
    # PriVA may still apply zygosity-/inheritance-aware prioritization elsewhere (e.g., PM3, BP2),
    # but here we keep the PVS1 "gene LoF mechanism" gate consistent across genotypes.
    # ``proband_gt_col`` is kept in the function signature for CLI/backward compatibility with
    # older PriVA calls; it is intentionally not used for PVS1 strength assignment.

    # THE GATE: IS LOSS OF FUNCTION A DISEASE MECHANISM FOR THIS GENE?
    # ================================================================
    #
    # Purely a question about the GENE. It deliberately reads none of the
    # variant's own mechanism scores: whether THIS null variant is convincing
    # is decided by the strength tree below, from the consequence, and mixing
    # the two would let the mechanism layer pre-empt PVS1's own judgement.
    #
    # Five independent sources, and ANY ONE of them opens the gate. They are
    # deliberately a union rather than a consensus: each is incomplete on its
    # own, and a gene missing from one is routinely present in another.
    #
    #   1. the HPO condition cache records a CURATED loss-of-function mechanism
    #      for one of the gene's germline conditions. Mechanisms the builder
    #      deduced from a recessive inheritance are excluded on purpose: that
    #      deduction restates the inheritance and adds no evidence here.
    #   2. the gene carries a ClinVar pathogenic or likely-pathogenic variant
    #      at a review status of two stars or better
    #   3. LOEUF below 0.35             -- intolerant to losing one copy
    #   4. mean AlphaMissense above 0.6 -- intolerant to missense generally
    #   5. ClinGen dosage curation      -- haploinsufficiency score 3, or
    #      30 / 40, which mark a gene whose phenotype is recessive and for
    #      which loss of function is therefore relevant when biallelic
    #
    #   lof_mechanism = 1 OR 2 OR 3 OR 4 OR 5
    #
    # THE STRENGTH TREE, entirely separate from the gate above
    # ========================================================
    #
    #   variant passes the LoF mechanism gate
    #             |
    #             +-- frameshift / stop_gained
    #             |         |
    #             |         +-- does NOT escape NMD ............. 4  Very Strong
    #             |         |
    #             |         +-- escapes NMD (last exon)
    #             |                   |
    #             |                   +-- spans an intolerant domain .. 3 Strong
    #             |                   +-- truncates >= 10% of protein . 3 Strong
    #             |                   +-- truncates <  10% of protein . 2 Moderate
    #             |
    #             +-- splice variant
    #             |         |
    #             |         +-- frameshift, not last exon ....... 4  Very Strong
    #             |         +-- frameshift, escapes NMD ......... 3/2 by the
    #             |         |                                     same two tests
    #             |         +-- induced inframe deletion ........ reduced strength
    #             |
    #             +-- inframe deletion .......................... reduced strength
    #             +-- 5'UTR frameshift (via UTRAnnotator) ....... assigned
    #             +-- start_lost (CDS alternative start codon) .. assigned
    #
    # is_exact_gof is read once more, far below, to withdraw PVS1 from a
    # variant the upstream hub has curated as gain of function. That is the
    # only mechanism score PVS1 reads, and it withdraws evidence rather than
    # granting it.
    mechanism_masks = _variant_mechanism_masks(df)

    # 1. the gene's own condition history, from the HPO condition cache
    hpo_lof_history = (
        df.get("gene_lof_mechanism_history", pd.Series("", index=df.index))
        .fillna("")
        .astype(str)
        .str.lower()
        .eq("true")
    )

    # 2. the gene carries ClinVar pathogenic variants at high review status
    clinvar_pathogenic = df["Gene"].isin(set(clinvar_pathogenic_genes or set()))

    # 3 and 4. constraint against losing a copy, and against missense generally
    loeuf_intolerant = df["LOEUF"].fillna(2) < 0.35
    am_intolerant = (
        df.get("Gene_avg_AM_score", pd.Series(np.nan, index=df.index)).fillna(0) > 0.6
    )

    # 5. ClinGen dosage curation. Score 3 is established haploinsufficiency;
    # 30 and 40 mark a gene whose phenotype is recessive, where loss of
    # function is still the relevant mechanism once both copies are affected.
    clingen_hi_score = pd.Series([np.nan] * len(df), index=df.index)
    if gene_dosage_sensitivity and os.path.exists(gene_dosage_sensitivity):
        try:
            clingen_dosage_df = pd.read_table(
                gene_dosage_sensitivity, low_memory=False
            ).dropna(subset=["#Gene Symbol", "Haploinsufficiency Score"])
            clingen_dosage_map = dict(
                zip(
                    clingen_dosage_df["#Gene Symbol"],
                    clingen_dosage_df["Haploinsufficiency Score"].astype(int),
                )
            )
            clingen_hi_score = df["SYMBOL"].map(clingen_dosage_map)
        except Exception as exc:
            logger.warning(
                "Failed to load ClinGen dosage sensitivity file %s: %s",
                gene_dosage_sensitivity,
                exc,
            )
    clingen_hi = clingen_hi_score == 3
    clingen_ar = clingen_hi_score.isin([30, 40])

    # Withhold PVS1 from genes whose only curated pathogenic mechanism is non-LoF
    # (GOF / dominant negative / triplosensitivity). Only the "only non-LoF" subset is
    # blocked: a gene in the non-LoF list that ALSO has a curated LoF history is mixed
    # (both mechanisms) and keeps PVS1. The non-LoF list overrides the statistical LoF
    # proxies (LOEUF / mean AlphaMissense) and the mechanism-agnostic ClinVar 2-star
    # pathogenic signal (ClinVar does not tag LoF vs non-LoF). Only an established
    # ClinGen haploinsufficiency call (HI score 3), which is an explicit LoF call,
    # overrides the list.
    non_lof_gene = pd.Series(False, index=df.index)
    if nonlof_mechanism_json and os.path.exists(nonlof_mechanism_json):
        non_lof_gene = df["SYMBOL"].isin(load_non_lof_gene_set(nonlof_mechanism_json))
    only_non_lof = non_lof_gene & ~hpo_lof_history
    non_lof_blocked = only_non_lof & ~clingen_hi
    logger.info(
        "PVS1 non-LoF gene list: %s of %s variants in non-LoF genes (%s only-non-LoF), %s blocked",
        int(non_lof_gene.sum()), len(df), int(only_non_lof.sum()), int(non_lof_blocked.sum()),
    )

    lof_mechanism = (
        hpo_lof_history
        | clinvar_pathogenic
        | loeuf_intolerant
        | am_intolerant
        | clingen_hi
        | clingen_ar
    ) & ~non_lof_blocked
    logger.info(
        "PVS1 LoF mechanism gate: %s of %s variants pass. "
        "HPO condition history %s, ClinVar pathogenic gene %s, "
        "LOEUF<0.35 %s, mean AM>0.7 %s, ClinGen HI=3 %s, ClinGen HI=30/40 %s",
        int(lof_mechanism.sum()), len(df),
        int(hpo_lof_history.sum()), int(clinvar_pathogenic.sum()),
        int(loeuf_intolerant.sum()), int(am_intolerant.sum()),
        int(clingen_hi.sum()), int(clingen_ar.sum()),
    )

    # Load the necessary dict file
    clinvar_patho_exon_af_dict = pickle.load(gzip.open(clinvar_patho_exon_af_stat)) if clinvar_patho_exon_af_stat.endswith(".gz") else pickle.load(open(clinvar_patho_exon_af_stat, 'rb'))
    interpro_entry_map_dict = pickle.load(gzip.open(interpro_entry_map_pkl)) if interpro_entry_map_pkl.endswith(".gz") else pickle.load(open(interpro_entry_map_pkl, 'rb'))
    tranx_pp_domain_map = pickle.load(gzip.open(tranx_exon_domain_map_pkl)) if tranx_exon_domain_map_pkl.endswith(".gz") else pickle.load(open(tranx_exon_domain_map_pkl, 'rb'))
    dm_instance = DomainNormalizer()
    intolerant_domains, exon_rare_patho_afs = zip(*df.apply(span_functional_domains, axis=1, dm_instance=dm_instance,
                                                                                             interpro_entry_map_dict=interpro_entry_map_dict,
                                                                                             tranx_pp_domain_map=tranx_pp_domain_map,
                                                                                             clinvar_patho_exon_af_dict=clinvar_patho_exon_af_dict,
                                                                                             intol_domains=intolerant_domains,
                                                                                             exon_patho_af_threshold=0.01))

    # Convert tuples to numpy arrays for boolean operations
    intolerant_domains = np.array(intolerant_domains)
    exon_rare_patho_afs = np.array(exon_rare_patho_afs)

    pvs1_criteria = np.zeros(len(df), dtype=int)
    # Deal with frameshift/nonsense variants
    non_fs_nmd_variants = (df['Consequence'].str.contains("stop_gained").fillna(False) | \
                           df['Consequence'].str.contains('frameshift').fillna(False)) & \
                           np.logical_not(df['NMD'].fillna(".").str.contains("escaping"))
    logger.info(f"{non_fs_nmd_variants.sum()} variants are frameshift/nonsense that predicted to cause NMD")
    pvs1_criteria[non_fs_nmd_variants & lof_mechanism] = 4

    # Deal with frameshift/nonsense variants that predicted to escape NMD
    non_fs_nmd_esp_variants = (df['Consequence'].str.contains("stop_gained").fillna(False) | df['Consequence'].str.contains('frameshift').fillna(False)) & \
                              (df['NMD'].fillna(".").str.contains("escaping") | df['LoF_filter'].fillna(".").str.contains("END_TRUNC"))
    logger.info(f"{non_fs_nmd_esp_variants.sum()} variants are frameshift/nonsense that predicted to escape NMD")

    # Check whether non-NMD escaping frameshift/nonsense variants spans functional domains
    non_fs_nmd_esp_intol_domain_vars = non_fs_nmd_esp_variants & intolerant_domains
    logger.info(f"{non_fs_nmd_esp_intol_domain_vars.sum()} variants are frameshift/nonsense that predicted to escape NMD and spans functional domains")
    pvs1_criteria[non_fs_nmd_esp_intol_domain_vars & lof_mechanism & (pvs1_criteria < 3)] = 3

    truncate_frac = truncate_fraction(df)
    pvs1_criteria[non_fs_nmd_esp_variants & lof_mechanism & ~intolerant_domains & exon_rare_patho_afs & (truncate_frac >= 0.1) & (pvs1_criteria < 3)] = 3
    pvs1_criteria[non_fs_nmd_esp_variants & lof_mechanism & ~intolerant_domains & exon_rare_patho_afs & (truncate_frac < 0.1) & (pvs1_criteria < 2)] = 2

    # Deal with inframe_deletion variants
    inframe_dels = df['Consequence'].str.contains("inframe_deletion")
    large_indels = (df["ref"].str.len() > 50) | (df["alt"].str.len() > 50) | df['EXON'].fillna("").str.match(r'^\d+-\d+/\d+$')

    pvs1_criteria[inframe_dels & intolerant_domains & large_indels] = 3
    pvs1_criteria[inframe_dels & lof_mechanism & ~intolerant_domains & exon_rare_patho_afs & large_indels & (truncate_frac >= 0.1) & (pvs1_criteria < 3)] = 3
    pvs1_criteria[inframe_dels & lof_mechanism & ~intolerant_domains & exon_rare_patho_afs & large_indels & (truncate_frac < 0.1) & (pvs1_criteria < 2)] = 2

    # Now we start to deal with splicing variants
    # Splicing frameshift variants that do NOT escape NMD -> PVS1 (Very Strong)
    # NMD escape occurs when frameshift affects only the last exon
    splicing_frameshift = df['splicing_frameshift'].fillna(False)
    splicing_last_exon_only = df['splicing_last_exon_only'].fillna(False) if 'splicing_last_exon_only' in df.columns else pd.Series(False, index=df.index)

    # Splicing frameshift with NMD (not last exon) -> PVS1=4
    splicing_fs_nmd = splicing_frameshift & ~splicing_last_exon_only
    pvs1_criteria[splicing_fs_nmd & lof_mechanism & (pvs1_criteria < 4)] = 4

    # Splicing frameshift escaping NMD (last exon only) -> reduced strength, same as regular NMD-escaping variants
    splicing_fs_nmd_escape = splicing_frameshift & splicing_last_exon_only
    # PVS1_Strong if spans intolerant domain
    pvs1_criteria[splicing_fs_nmd_escape & df['splicing_span_intol_domain'].fillna(False) & lof_mechanism & (pvs1_criteria < 3)] = 3
    # PVS1_Strong if truncates >= 10% protein
    pvs1_criteria[splicing_fs_nmd_escape & ~df['splicing_span_intol_domain'].fillna(False) & df['splicing_ten_percent_protein'].fillna(False) & lof_mechanism & exon_rare_patho_afs & (pvs1_criteria < 3)] = 3
    # PVS1_Moderate if truncates < 10% protein
    pvs1_criteria[splicing_fs_nmd_escape & ~df['splicing_span_intol_domain'].fillna(False) & ~df['splicing_ten_percent_protein'].fillna(False) & lof_mechanism & exon_rare_patho_afs & (pvs1_criteria < 2)] = 2

    # Splicing-induced inframe del variants (not frameshift)
    splicing_inframe_intol_domains = (df['splicing_lof'] | df['splicing_len_changing']) & \
                                      ~splicing_frameshift & ~splicing_last_exon_only & \
                                      df['splicing_span_intol_domain']
    pvs1_criteria[splicing_inframe_intol_domains & lof_mechanism & (pvs1_criteria < 3)] = 3
    pvs1_criteria[(df['splicing_lof'] | df['splicing_len_changing']) & \
                  ~splicing_frameshift & \
                  ~df['splicing_span_intol_domain'].fillna(False) & \
                  exon_rare_patho_afs & \
                  df['splicing_ten_percent_protein'] & \
                  lof_mechanism & \
                  (pvs1_criteria < 3)] = 3
    pvs1_criteria[(df['splicing_lof'] | df['splicing_len_changing']) & \
                  ~splicing_frameshift & \
                  ~df['splicing_span_intol_domain'].fillna(False) & \
                  exon_rare_patho_afs & \
                  ~df['splicing_ten_percent_protein'] & \
                  lof_mechanism & \
                  (pvs1_criteria < 2)] = 2

    # Thanks to UTRAnnotator, we also implement the PVS1 for frameshift 5UTR variants
    five_utr_frameshift = df['5UTR_frameshift']  #5UTR frameshift bound to NMDs
    five_utr_inframe_intol_domains = df['5UTR_len_changing'] & df['5UTR_span_intol_domain'] & ~five_utr_frameshift
    pvs1_criteria[five_utr_inframe_intol_domains & lof_mechanism & (pvs1_criteria < 3)] = 3
    pvs1_criteria[five_utr_frameshift & lof_mechanism & (pvs1_criteria < 4)] = 4

    # Deal with start_lost variants using CDS-based alternative start codon detection
    # Per ClinGen SVI guidelines (PMC6185798), we check for downstream in-frame ATG codons
    # within the same transcript, not alternative transcripts
    if cds_fasta_path and os.path.exists(cds_fasta_path):
        logger.info(f"Using CDS-based alternative start codon detection from: {cds_fasta_path}")
        pvs1_start_lost, df = evaluate_start_lost_with_cds(
            df,
            cds_fasta_path,
            intolerant_domains,
            lof_mechanism  # Pass gene LoF mechanism status for proper PVS1 strength
        )
        # Only update where start_lost PVS1 is higher than current PVS1
        start_lost_mask = df["Consequence"].str.contains("start_lost").fillna(False)
        pvs1_criteria[start_lost_mask & (pvs1_start_lost > pvs1_criteria)] = pvs1_start_lost[start_lost_mask & (pvs1_start_lost > pvs1_criteria)]
        # For start_lost without PVS1 assigned yet, use the CDS-based result
        pvs1_criteria[start_lost_mask & (pvs1_criteria == 0)] = pvs1_start_lost[start_lost_mask & (pvs1_criteria == 0)]
    else:
        # Fallback: use the old transcript-based approach (deprecated)
        logger.warning("CDS FASTA file not provided. Using deprecated transcript-based approach for start_lost variants.")
        logger.warning("For accurate start_lost evaluation per ClinGen SVI guidelines, please provide cds_fasta_path")
        alt_start_genes = set(df.groupby(['chrom', 'pos', 'ref', 'alt', 'Gene']).apply(identify_alternative_start_codon_genes).unique().tolist())
        logger.info(f"These are the {len(alt_start_genes)} genes that have functional transcripts using alternative start codons: {alt_start_genes}")
        alt_start_losts = df["Consequence"].str.contains("start_lost") & df["Gene"].isin(alt_start_genes)
        logger.info(f"{alt_start_losts.sum()} variants are having start_lost consequences to transcripts with alternative start codons")
        # Per ClinGen SVI, start_lost should be PVS1_Moderate at most (not Very Strong)
        pvs1_criteria[df['Consequence'].str.contains("start_lost") & ~alt_start_losts & (pvs1_criteria < 2)] = 2

    exact_gof = mechanism_masks["is_exact_gof"]
    suppress_pvs1 = (pvs1_criteria > 0) & exact_gof.to_numpy()
    if suppress_pvs1.any():
        logger.info(
            "Suppressing PVS1 for %s exact GoFCards variant-level GOF rows",
            int(suppress_pvs1.sum()),
        )
        pvs1_criteria[suppress_pvs1] = 0

    return pvs1_criteria, intolerant_domains
