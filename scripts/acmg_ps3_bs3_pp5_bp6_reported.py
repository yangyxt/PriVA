#!/usr/bin/env python3
"""PS3, BS3, PP5 and BP6 -- what other people have already concluded.

    PS3 / BS3  a well-established functional assay shows the variant damages
               or does not damage gene function. Read from ClinVar at two
               stars or better, and from MAVE multiplexed assay scores where
               a MaveDB study covers the residue.
    PP5 / BP6  a reputable source reports the variant as pathogenic or benign
               without supplying the evidence. These read ClinVar's
               classification directly, and deliberately exclude anything
               already counted as PS3 or BS3 -- otherwise the same assertion
               would be counted twice, once as evidence and once as reputation.

PS3_BS3_criteria also returns the clinvar_patho and clinvar_benign masks that
PS1/PM5, PM1 and PP5/BP6 all consume, which is why those four sit together.
"""

import logging
import pandas as pd
import numpy as np

from mavedb_interpreter import pick_interpreter_from_metadata


logger = logging.getLogger(__name__)


def mavedb_interpretation_per_row(row: pd.Series, metadata_df: pd.DataFrame) -> bool:
    '''
    Interpret the MaveDB scores and return the PS3 and BS3 criteria.

    When a variant has multiple MaveDB assays, we count evidence for PS3 (functional impact)
    and BS3 (no functional impact) separately. If there's conflicting evidence from different
    assays, we resolve by majority vote. If tied, we conservatively set both to False.

    Tracks URN IDs and assay types for robust logging.
    '''
    if pd.isna(row.get('MaveDB_urn', np.nan)):
        row["MaveDB_PS3"] = False
        row["MaveDB_BS3"] = False
        row["MaveDB_score_interpretation"] = np.nan
        return row

    high_confs = str(row.get('MaveDB_high_conf', '')).split('&')
    pvalues = str(row.get('MaveDB_pvalue', '')).split('&')
    scores = str(row.get('MaveDB_score', '')).split('&')
    urn_sets = str(row.get('MaveDB_urn')).split('&')

    score_interpretation = []

    # Track evidence with detailed information: list of (urn, assay_family, score, call)
    ps3_evidence = []  # URNs and assay types supporting PS3 (functional impact)
    bs3_evidence = []  # URNs and assay types supporting BS3 (no functional impact)

    for i, urn in enumerate(urn_sets):
        if urn not in metadata_df['urn'].values:
            logger.warning(f"The URN {urn} is not in the MaveDB metadata, it is ignored")
            continue

        # Get metadata for this URN
        meta_row = metadata_df.loc[metadata_df['urn'] == urn].iloc[0]
        interpretation = meta_row['rationale']
        assay_family = meta_row.get('assay_family', 'Unknown')
        assay_subtype = meta_row.get('subtype_value', 'Unknown')

        score_interpretation.append(interpretation)
        score = scores[i]
        high_conf = high_confs[i]
        pvalue = pvalues[i]

        try:
            score_float = float(score)
        except ValueError:
            continue

        try:
            pvalue_float = float(pvalue)
        except ValueError:
            pvalue_float = np.nan

        if pd.isna(pvalue_float):
            result = pick_interpreter_from_metadata(urn, score_float, metadata_df)
            call = result.get('call', 'Unknown')
            if result['MaveDB_PS3']:
                ps3_evidence.append({
                    'urn': urn,
                    'assay_family': assay_family,
                    'assay_subtype': assay_subtype,
                    'score': score_float,
                    'call': call,
                    'reason': result.get('reason', 'Unknown')
                })
            if result['MaveDB_BS3']:
                bs3_evidence.append({
                    'urn': urn,
                    'assay_family': assay_family,
                    'assay_subtype': assay_subtype,
                    'score': score_float,
                    'call': call,
                    'reason': result.get('reason', 'Unknown')
                })
        elif pvalue_float < 0.05:
            ps3_evidence.append({
                'urn': urn,
                'assay_family': assay_family,
                'assay_subtype': assay_subtype,
                'score': score_float,
                'call': 'pvalue<0.05',
                'reason': f'pvalue={pvalue_float:.4f}'
            })
        elif pvalue_float >= 0.05:
            bs3_evidence.append({
                'urn': urn,
                'assay_family': assay_family,
                'assay_subtype': assay_subtype,
                'score': score_float,
                'call': 'pvalue>=0.05',
                'reason': f'pvalue={pvalue_float:.4f}'
            })

    variant_id = f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}"

    # Format evidence for logging
    def format_evidence(evidence_list):
        return "; ".join([f"{e['urn']}({e['assay_family']}/{e['assay_subtype']}, score={e['score']:.4f}, call={e['call']})"
                         for e in evidence_list])

    # === Per-assay-family resolution ===
    # Group evidence by assay_family to handle conflicts within the same assay type
    ps3_families = set(e['assay_family'] for e in ps3_evidence)
    bs3_families = set(e['assay_family'] for e in bs3_evidence)

    # Families that have ONLY PS3 evidence (no conflicting BS3 from same family)
    ps3_only_families = ps3_families - bs3_families
    # Families that have ONLY BS3 evidence (no conflicting PS3 from same family)
    bs3_only_families = bs3_families - ps3_families
    # Families that have BOTH PS3 and BS3 (internal conflict)
    conflicted_families = ps3_families & bs3_families

    # Get evidence from non-conflicted families
    ps3_clean_evidence = [e for e in ps3_evidence if e['assay_family'] in ps3_only_families]
    bs3_clean_evidence = [e for e in bs3_evidence if e['assay_family'] in bs3_only_families]

    # For conflicted families, resolve by majority within each family
    ps3_resolved_evidence = []
    bs3_resolved_evidence = []

    for family in conflicted_families:
        family_ps3 = [e for e in ps3_evidence if e['assay_family'] == family]
        family_bs3 = [e for e in bs3_evidence if e['assay_family'] == family]

        if len(family_ps3) > len(family_bs3):
            # PS3 wins within this family
            ps3_resolved_evidence.extend(family_ps3)
            logger.debug(f"Variant {variant_id}: Family '{family}' conflict resolved to PS3 by majority "
                        f"(PS3:{len(family_ps3)} vs BS3:{len(family_bs3)})")
        elif len(family_bs3) > len(family_ps3):
            # BS3 wins within this family
            bs3_resolved_evidence.extend(family_bs3)
            logger.debug(f"Variant {variant_id}: Family '{family}' conflict resolved to BS3 by majority "
                        f"(PS3:{len(family_ps3)} vs BS3:{len(family_bs3)})")
        else:
            # Tie within family - this family contributes nothing (conservative)
            logger.debug(f"Variant {variant_id}: Family '{family}' has tied conflict "
                        f"(PS3:{len(family_ps3)} vs BS3:{len(family_bs3)}), excluded from final decision")

    # Combine clean evidence with resolved evidence
    final_ps3_evidence = ps3_clean_evidence + ps3_resolved_evidence
    final_bs3_evidence = bs3_clean_evidence + bs3_resolved_evidence

    # Final decision
    mavedb_ps3 = len(final_ps3_evidence) > 0
    mavedb_bs3 = len(final_bs3_evidence) > 0

    # Logging
    if mavedb_ps3 and mavedb_bs3:
        # Both PS3 and BS3 from different assay families - this is valid!
        logger.info(f"Variant {variant_id} has PS3 and BS3 from DIFFERENT assay families (both retained). "
                   f"PS3 families: {ps3_only_families | set(e['assay_family'] for e in ps3_resolved_evidence)}, "
                   f"BS3 families: {bs3_only_families | set(e['assay_family'] for e in bs3_resolved_evidence)}. "
                   f"PS3: [{format_evidence(final_ps3_evidence)}] | "
                   f"BS3: [{format_evidence(final_bs3_evidence)}]")
    elif conflicted_families:
        # Had conflicts but resolved
        logger.debug(f"Variant {variant_id} had conflicts in families {conflicted_families}, resolved. "
                    f"Final: PS3={mavedb_ps3}, BS3={mavedb_bs3}. "
                    f"PS3: [{format_evidence(final_ps3_evidence)}] | "
                    f"BS3: [{format_evidence(final_bs3_evidence)}]")

    row["MaveDB_PS3"] = mavedb_ps3
    row["MaveDB_BS3"] = mavedb_bs3
    row["MaveDB_score_interpretation"] = "&".join([str(x) for x in score_interpretation])
    return row


def mavedb_score_interpretation(df: pd.DataFrame, mavedb_metadata: pd.DataFrame) -> pd.DataFrame:
    '''
    Interpret the MaveDB scores and return the PS3 and BS3 criteria
    '''
    df = df.apply(mavedb_interpretation_per_row, axis=1, metadata_df=mavedb_metadata)
    mavedb_ps3_recs = df.loc[df['MaveDB_PS3'], [ "MaveDB_pvalue", "MaveDB_score", "MaveDB_high_conf", "MaveDB_score_interpretation"]]
    logger.info(f"There are {df['MaveDB_PS3'].sum()} variants with MaveDB determined PS3 criteria, they are determined by these functional assays: \n{mavedb_ps3_recs[:10].to_string()}")
    mavedb_bs3_recs = df.loc[df['MaveDB_BS3'], [ "MaveDB_pvalue", "MaveDB_score", "MaveDB_high_conf", "MaveDB_score_interpretation"]]
    logger.info(f"There are {df['MaveDB_BS3'].sum()} variants with MaveDB determined BS3 criteria, they are determined by these functional assays: \n{mavedb_bs3_recs[:10].to_string()}")
    return df


def PS3_BS3_criteria(df: pd.DataFrame, mavedb_metadata_tsv: str = "", high_confidence_status = {
        # Higher confidence (2+ stars)
        'practice_guideline': 4,                                   # 4 stars
        'reviewed_by_expert_panel': 3,                             # 3 stars
        'criteria_provided,_multiple_submitters,_no_conflicts': 2,  # 2 stars
    }) -> pd.DataFrame:
    # Basically rely on ClinVar annotations

    result_dict = {"PS3": [], "BS3": [], "clinvar_patho": [], "clinvar_benign": []}

    clinvar_lof = df['CLNSIG'].fillna("").str.contains('Pathogenic') & (df['CLNREVSTAT'].map(high_confidence_status) >= 2)
    # clinvar_lof = clinvar_lof | (df['CLNSIG'].fillna("").str.contains('athogenic') & (df['CLNREVSTAT'].map(high_confidence_status) >= 3)) # Including Likely_pathogenic
    result_dict["clinvar_patho"] = clinvar_lof

    high_conf_benign = df['CLNSIG'].fillna("").str.contains('Benign') & (df['CLNREVSTAT'].map(high_confidence_status) >= 2)
    # high_conf_benign = high_conf_benign | (df['CLNSIG'].fillna("").str.contains('enign') & (df['CLNREVSTAT'].map(high_confidence_status) >= 3)) # Including Likely_benign
    result_dict["clinvar_benign"] = high_conf_benign

    ps3_array = np.zeros(len(df), dtype=int)
    bs3_array = np.zeros(len(df), dtype=int)

    if mavedb_metadata_tsv:
        mavedb_metadata = pd.read_table(mavedb_metadata_tsv, low_memory=False)
        mavedb_metadata.drop_duplicates(subset=["urn"], inplace=True)
        logger.info(f"There are {len(mavedb_metadata)} unique URNs in the MaveDB metadata, which looks like: \n{mavedb_metadata.head().to_string(index=False)}")
        df = mavedb_score_interpretation(df, mavedb_metadata)
        ps3_criteria = clinvar_lof | df['MaveDB_PS3']
        bs3_criteria = high_conf_benign | df['MaveDB_BS3']
    else:
        ps3_criteria = clinvar_lof
        bs3_criteria = high_conf_benign

    ps3_array[ps3_criteria] = 3
    bs3_array[bs3_criteria] = 3

    result_dict["PS3"] = ps3_array
    result_dict["BS3"] = bs3_array

    return result_dict, df


def PP5_BP6_criteria(df: pd.DataFrame, clinvar_patho, clinvar_benign) -> pd.Series:
    # PP5: The variant is reported as pathogenic by a reputable source but without to many supporting evidences
    pp5_criteria = df['CLNSIG'].fillna("").str.contains('athogenic') & ~clinvar_patho
    bp6_criteria = df['CLNSIG'].fillna("").str.contains('enign') & ~clinvar_benign

    pp5_array = np.zeros(len(df), dtype=int)
    bp6_array = np.zeros(len(df), dtype=int)
    pp5_array[pp5_criteria] = 1
    bp6_array[bp6_criteria] = 1
    return pp5_array, bp6_array
