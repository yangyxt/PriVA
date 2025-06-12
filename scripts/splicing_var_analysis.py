import os
import pickle
import gzip
import pandas as pd
import numpy as np
import argparse as ap
import logging
import uuid
import multiprocessing as mp
from functools import partial
import re
import sys

from combine_annotations import na_value
from protein_domain_mapping import DomainNormalizer

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)



def splicing_interpretation(anno_table, 
                            transcript_domain_map_pkl, 
                            intolerant_domains_pkl, 
                            interpro_entry_map_pkl,
                            spliceai_cutoff = 0.8, 
                            threads = 12):
    """Process annotation table in parallel to interpret splicing effects"""
    anno_df = pd.read_table(anno_table, low_memory=False)

    intolerant_domains = pickle.load(gzip.open(intolerant_domains_pkl, 'rb')) if intolerant_domains_pkl.endswith('.gz') else pickle.load(open(intolerant_domains_pkl, 'rb'))
    transcript_domain_map = pickle.load(gzip.open(transcript_domain_map_pkl, 'rb')) if transcript_domain_map_pkl.endswith('.gz') else pickle.load(open(transcript_domain_map_pkl, 'rb'))
    interpro_entry_map_dict = pickle.load(gzip.open(interpro_entry_map_pkl, 'rb')) if interpro_entry_map_pkl.endswith('.gz') else pickle.load(open(interpro_entry_map_pkl, 'rb'))

    dm_instance = DomainNormalizer()
    # Create partial function with fixed arguments
    process_row = partial(splicing_altering_per_row, 
                         transcript_domain_map=transcript_domain_map,
                         intolerant_domains=intolerant_domains,
                         interpro_entry_map_dict=interpro_entry_map_dict,
                         spliceai_cutoff=spliceai_cutoff,
                         dm_instance=dm_instance)
    
    # Get number of CPU cores (leave some free for other processes)
    logger.info(f"Processing {len(anno_df)} variants using {threads} cores, spliceai_cutoff is {spliceai_cutoff}")
    
    # Process rows in parallel
    with mp.Pool(threads) as pool:
        results = pool.map(process_row, [row for _, row in anno_df.iterrows()])
    
    # Unpack results
    splicing_lof, splicing_len_changing, splicing_frame_shift, splicing_span_intol_domain, splicing_ten_percent_protein, splicing_canonical_gain_offset, splicing_canonical_loss_offset, splicing_affected_exons = zip(*results)
    
    # Add results to dataframe
    anno_df['splicing_lof'] = splicing_lof
    anno_df['splicing_len_changing'] = splicing_len_changing
    anno_df['splicing_frameshift'] = splicing_frame_shift
    anno_df['splicing_span_intol_domain'] = splicing_span_intol_domain
    anno_df['splicing_affected_exons'] = splicing_affected_exons
    anno_df['splicing_ten_percent_protein'] = splicing_ten_percent_protein
    anno_df['spliceai_canonical_gain_offset'] = splicing_canonical_gain_offset
    anno_df['spliceai_canonical_loss_offset'] = splicing_canonical_loss_offset

    # Save results
    tmp_tag = str(uuid.uuid4())
    output_file = f"{anno_table}.{tmp_tag}.splicing.anno.txt"
    anno_df.to_csv(output_file, sep='\t', index=False)

    # mv the temp output file to anno_table to replace the original one
    os.system(f"mv {output_file} {anno_table}")
    logger.info(f"Results saved to {anno_table}")



def splicing_altering_per_row(row, 
                              transcript_domain_map, 
                              intolerant_domains,
                              interpro_entry_map_dict,
                              dm_instance = None,
                              spliceai_cutoff = 0.8):
    '''
    Basically, we just refer to SpliceVault top-4 events and SpliceAI delta scores.

    The possible splicing altering effects can be categorized into:
    1. Frameshift effect, can be induced by exon skipping, cryptic donor/acceptor activation.
    2. In-frame retention, can be induced by exon skipping, cryptic donor/acceptor activation.
    3. In-frame truncation, can be induced by exon skipping, cryptic donor/acceptor activation.
    '''
    transcript_id = row['Feature']
    splicevault_events = row['SpliceVault_top_events']
    spliceai_delta = float(row['SpliceVault_SpliceAI_delta'])
    intron_pos = None if na_value(row['INTRON']) else row['INTRON']
    exon_pos = None if na_value(row['EXON']) else row['EXON']
    cds_pos = None if na_value(row["CDS_position"]) else row["CDS_position"]
    # Extract HGVSc if available
    hgvsc = None if na_value(row.get('HGVSc', None)) else row.get('HGVSc', None)
    
    logger.debug(f"The input transcript_id is {transcript_id}, splicevault_events is {splicevault_events}, intron_pos is {intron_pos}, exon_pos is {exon_pos}, cds_pos is {cds_pos}, hgvsc is {hgvsc}")

    # First get SpliceVault interpretation
    splicevault_lof, splicevault_len_changing, splicevault_frame_shift, splicevault_span_intol_domain, splicevault_affected_exons, splicevault_ten_percent_protein = SpliceVault_interpretation(splicevault_events, 
                                                                                                                                                                                                spliceai_delta,
                                                                                                                                                                                                transcript_id, 
                                                                                                                                                                                                transcript_domain_map, 
                                                                                                                                                                                                intolerant_domains,
                                                                                                                                                                                                interpro_entry_map_dict,
                                                                                                                                                                                                cds_pos,
                                                                                                                                                                                                intron_pos,
                                                                                                                                                                                                exon_pos,
                                                                                                                                                                                                dm_instance )
        
    
    # Then get SpliceAI interpretation
    DP_AG = float(row['SpliceAI_pred_DP_AG'])
    DP_AL = float(row['SpliceAI_pred_DP_AL'])
    DP_DG = float(row['SpliceAI_pred_DP_DG'])
    DP_DL = float(row['SpliceAI_pred_DP_DL'])
    DS_AG = float(row['SpliceAI_pred_DS_AG'])
    DS_AL = float(row['SpliceAI_pred_DS_AL'])
    DS_DG = float(row['SpliceAI_pred_DS_DG'])
    DS_DL = float(row['SpliceAI_pred_DS_DL'])
    var_pos = int(row['pos'])
    strand = str(row['STRAND'])
    hgvsc = str(row['HGVSc']) if not na_value(row['HGVSc']) else None
    spliceai_lof, spliceai_len_changing, spliceai_frame_shift, spliceai_span_intol_domain, spliceai_affected_exons, spliceai_ten_percent_protein, spliceai_canonical_gain_offset, spliceai_canonical_loss_offset = SpliceAI_interpretation( DS_AG, 
                                                                                                                                                                                                                                            DS_AL,
                                                                                                                                                                                                                                            DS_DG, 
                                                                                                                                                                                                                                            DS_DL, 
                                                                                                                                                                                                                                            DP_AG, 
                                                                                                                                                                                                                                            DP_AL, 
                                                                                                                                                                                                                                            DP_DG, 
                                                                                                                                                                                                                                            DP_DL,
                                                                                                                                                                                                                                            transcript_id, 
                                                                                                                                                                                                                                            transcript_domain_map, 
                                                                                                                                                                                                                                            intolerant_domains,
                                                                                                                                                                                                                                            consequence = row['Consequence'],
                                                                                                                                                                                                                                            cds_pos = cds_pos,
                                                                                                                                                                                                                                            intron_pos = intron_pos,
                                                                                                                                                                                                                                            exon_pos = exon_pos,
                                                                                                                                                                                                                                            hgvsc = hgvsc,  # Pass HGVSc
                                                                                                                                                                                                                                            var_pos = var_pos,
                                                                                                                                                                                                                                            strand = strand,
                                                                                                                                                                                                                                            cutoff = spliceai_cutoff,
                                                                                                                                                                                                                                            interpro_entry_map_dict = interpro_entry_map_dict,
                                                                                                                                                                                                                                            dm_instance = dm_instance)

    return splicevault_lof or spliceai_lof, \
           splicevault_len_changing or spliceai_len_changing, \
           splicevault_frame_shift or spliceai_frame_shift, \
           splicevault_span_intol_domain or spliceai_span_intol_domain, \
           splicevault_ten_percent_protein or spliceai_ten_percent_protein, \
           spliceai_canonical_gain_offset, \
           spliceai_canonical_loss_offset, \
           ",".join(set(list(map(str, splicevault_affected_exons)) + list(map(str, spliceai_affected_exons))))
    


def parse_hgvsc_splice_position(hgvsc: str, 
                                strand: str,
                                var_pos: int,
                                var_length: int = 1) -> dict:
    """
    Parse HGVSc annotation to determine variant position relative to canonical splice sites.
    Handles SNVs, single-base indels, multi-base indels, duplications, and range insertions.
    
    Args:
        hgvsc: HGVSc annotation (e.g., "ENST00000231790.2:c.116+1G>A")
        
    Returns:
        dict with keys:
        - variant_type: 'donor_region', 'acceptor_region', 'exonic', 'spanning'
        - distance_to_canonical: distance from variant to nearest canonical site
        - canonical_site_type: 'donor' or 'acceptor' 
        - exon_boundary_pos: the exon boundary position
        - intronic_offset: offset into intron (+N or -N)
        - variant_length: length of the variant change
    """
    result = {
        'variant_type': 'unknown',
        'distance_to_canonical': None,  # Should be the distance from the POS to the first base of the canonical site
        'canonical_site_type': None,
        'exon_boundary_pos': None, # Last/First base of the exon in the 5->3 direction
        "intronic_offset": None, # POS of the variant to relative to the first/last base of the exon
        "exon_skipping": False,
        "overlapping_canonical_site": False
    }

    if not isinstance(hgvsc, str):
        return None
    
    if ':c.' not in hgvsc:
        return None
    
    if (not "-" in hgvsc) and (not "+" in hgvsc):
        result["variant_type"] = "exonic"
        return None

    if "*" in hgvsc:
        logger.warning(f"The HGVSc: {hgvsc} contains a * character, which represents a 3UTR variants HGVSc annotation")
        return None
        
    # Extract the coding sequence part
    coding_part = hgvsc.split(':c.')[1]
    
    # Pattern for range deletions spanning exon-intron or within intron
    # e.g., "3646_3646+3del", "207+1_207+2del", inclusive on both ends
    range_del_pattern = r'(\d+)?(?:(\+|-)(\d+))?_(\d+)?(?:(\+|-)(\d+))?del'
    
    # Pattern for range insertions/duplications in intronic regions
    # e.g., "2746-17_2746-14dup", "4644+12_4644+13dup"
    range_dup_pattern = r'(\d+)?(?:(\+|-)(\d+))?_(\d+)?(?:(\+|-)(\d+))?dup'
    # e.g., "4644+12_4644+13insTA", "4644+12_4644+13insAAA"
    range_ins_pattern = r'(\d+)?(?:(\+|-)(\d+))?_(\d+)?(?:(\+|-)(\d+))?ins[ATCG]+'
    
    # Check for range deletions first
    range_del_match = re.search(range_del_pattern, coding_part)
    if range_del_match:
        start_offset = int(range_del_match.group(2) + range_del_match.group(3)) if range_del_match.group(2) else 0 # start pos of the variant deviation with the exon boundary
        end_offset = int(range_del_match.group(5) + range_del_match.group(6)) if range_del_match.group(5) else 0 # end pos of the variant deviation with the exon boundary
        
        if start_offset * end_offset < 0:
            # Either the HGVSc is problematic or the deletion is across the entire exon, causing exon skipping
            result["variant_type"] = "spanning"
            result["exon_skipping"] = True
            return result
        elif start_offset >= 0 and end_offset >= 0:
            # The deletion is downstream of the preceding exon, near the donor site
            result["variant_type"] = "donor_region"
            result["canonical_site_type"] = "donor"
            if strand == "1":
                result["distance_to_canonical"] = start_offset - 1 - 1 # e.g. "3646_3646+3del", start_offset is 0, the POS of this deletion should be 3645, it's 2 bp upstream of 3646 + 1
                result["exon_boundary_pos"] = var_pos - start_offset + 1
            elif strand == "-1":
                result["distance_to_canonical"] = -(end_offset + 1) + 2 # e.g. "3646_3646+3del", end_offset is +3, the POS of this deletion should be 3646+4, it's 2 bp upstream of 3646 + 2 (pos strand direction, start of the canonical donor site)
                result["exon_boundary_pos"] = var_pos + end_offset + 1
            else:
                raise ValueError(f"Invalid strand: {strand}")
        elif start_offset <= 0 and end_offset <= 0:
            # The deletion is upstream of the following exon, near the acceptor site
            result["variant_type"] = "acceptor_region"
            result["canonical_site_type"] = "acceptor"
            if strand == "1":
                result["distance_to_canonical"] = start_offset - 1 + 2 # e.g. 44-7_44-5del, start_offset is -7, the POS of this deletion should be 44-8, it's 6 bp upstream of 44-2 (start of the acceptor site, pos strand direction)
                result["exon_boundary_pos"] = var_pos + 1 - start_offset
            elif strand == "-1":
                result["distance_to_canonical"] = end_offset + 1 + 1 # e.g. 44-7_44-5del, end_offset is -5, the POS of this deletion should be 44-4, it's 3 bp downstream of 44-1 (start of the acceptor site, pos strand direction)
                result["exon_boundary_pos"] = var_pos + 1 + end_offset
            else:
                raise ValueError(f"Invalid strand: {strand}")
    
    # Check for range insertions/duplications
    range_ins_match = re.search(range_ins_pattern, coding_part)
    if range_ins_match:
        start_offset = int(range_ins_match.group(2) + range_ins_match.group(3)) if range_ins_match.group(2) else 0 # start pos of the variant deviation with the exon boundary
        end_offset = int(range_ins_match.group(5) + range_ins_match.group(6)) if range_ins_match.group(5) else 0 # end pos of the variant deviation with the exon boundary
        
        if start_offset * end_offset < 0:
            # Either the HGVSc is problematic or the deletion is across the entire exon, causing exon skipping
            result["variant_type"] = "spanning"
            result["exon_skipping"] = True
            return result
        elif start_offset >= 0 and end_offset >= 0:
            # The insertion is downstream of the preceding exon, near the donor site
            result["variant_type"] = "donor_region"
            result["canonical_site_type"] = "donor"
            if strand == "1":
                result["distance_to_canonical"] = start_offset - 1 # e.g. "4644+12_4644+13insTA", start_offset is 12, the POS of this insertion should be 4644+12, it's 11 bp downstream of 4644+1 
                result["exon_boundary_pos"] = var_pos - start_offset
            elif strand == "-1":
                result["distance_to_canonical"] = 2 - end_offset # e.g. "4644+12_4644+13insTA", end_offset is 13, the POS of this insertion should be 4644+13, it's 11 bp upstream of 4644+2 (pos strand direction, start of the canonical donor site)
                result["exon_boundary_pos"] = var_pos + end_offset
            else:
                raise ValueError(f"Invalid strand: {strand}")
        elif start_offset <= 0 and end_offset <= 0:
            # The deletion is upstream of the following exon, near the acceptor site
            result["variant_type"] = "acceptor_region"
            result["canonical_site_type"] = "acceptor"
            if strand == "1":
                result["distance_to_canonical"] = start_offset + 2 # e.g. 248-4_248-3insAAGTTTT, start_offset is -4, the POS of this insertion should be 248-4, it's 2 bp upstream of 248-2 (start of the acceptor site)
                result["exon_boundary_pos"] = var_pos - start_offset
            elif strand == "-1":
                result["distance_to_canonical"] = - 1 - end_offset # e.g. 248-4_248-3insAAGTTTT, end_offset is -3, the POS of this insertion should be 248-3, it's 2 bp downstream of 248-1 (start of the acceptor site, pos strand direction)
                result["exon_boundary_pos"] = var_pos + end_offset
            else:
                raise ValueError(f"Invalid strand: {strand}")

    range_dup_match = re.search(range_dup_pattern, coding_part)
    if range_dup_match:
        start_offset = int(range_dup_match.group(2) + range_dup_match.group(3)) if range_dup_match.group(2) else 0 # start pos of the variant deviation with the exon boundary
        end_offset = int(range_dup_match.group(5) + range_dup_match.group(6)) if range_dup_match.group(5) else 0 # end pos of the variant deviation with the exon boundary
        
        if start_offset * end_offset < 0:
            # Either the HGVSc is problematic or the deletion is across the entire exon, causing exon skipping
            result["variant_type"] = "spanning"
            result["exon_skipping"] = True
            return result
        elif start_offset >= 0 and end_offset >= 0:
            # The insertion is downstream of the preceding exon, near the donor site
            result["variant_type"] = "donor_region"
            result["canonical_site_type"] = "donor"
            if strand == "1":
                result["distance_to_canonical"] = start_offset - 1 - 1 # e.g. "4644+12_4644+13dup", start_offset is 12, the POS of this duplication should be 4644+11, it's 10 bp downstream of 4644+1 
                result["exon_boundary_pos"] = var_pos - start_offset + 1
            elif strand == "-1":
                result["distance_to_canonical"] = -(end_offset + 1) + 2 # e.g. "4644+12_4644+13dup", end_offset is 13, the POS of this duplication should be 4644+14, it's 12 bp upstream of 4644+2 (pos strand direction, start of the canonical donor site)
                result["exon_boundary_pos"] = var_pos + end_offset + 1
            else:
                raise ValueError(f"Invalid strand: {strand}")
        elif start_offset <= 0 and end_offset <= 0:
            # The deletion is upstream of the following exon, near the acceptor site
            result["variant_type"] = "acceptor_region"
            result["canonical_site_type"] = "acceptor"
            if strand == "1":
                result["distance_to_canonical"] = start_offset - 1 + 2 # e.g. "2746-17_2746-14dup", start_offset is -17, the POS of this duplication should be 2746-18, it's 16 bp upstream of 2746-2 (start of the acceptor site)
                result["exon_boundary_pos"] = var_pos - start_offset + 1
            elif strand == "-1":
                result["distance_to_canonical"] = -(end_offset + 1) - 1 # e.g. "2746-17_2746-14dup", end_offset is -14, the POS of this duplication should be 2746-13, it's 12 bp downstream of 2746-1 (start of the acceptor site, pos strand direction)
                result["exon_boundary_pos"] = var_pos + end_offset + 1
            else:
                raise ValueError(f"Invalid strand: {strand}")
    
    # Pattern for single position indels with intronic offset
    # e.g., "3354-1del", "459+2dup"
    single_indel_pattern = r'(\d+)?(\+|-)(\d+)(del|dup)'

    single_indel_match = re.search(single_indel_pattern, coding_part)
    if single_indel_match:
        offset = int(single_indel_match.group(2) + single_indel_match.group(3))
        if offset > 0:
            result["variant_type"] = "donor_region"
            result["canonical_site_type"] = "donor"
            if strand == "1":
                result["distance_to_canonical"] = offset - 1 - 1   # "459+2dup", offset is 2, the POS is 459+1, it is at the position of start of donor at 459+1
                result["exon_boundary_pos"] = var_pos - offset + 1
            elif strand == "-1":
                result["distance_to_canonical"] = - (offset + 1) + 2   # "459+2dup", offset is 2, the POS is 459+3, it is 1bp upstream of start of donor at 459+2 (pos strand direction, start of the canonical donor site)
                result["exon_boundary_pos"] = var_pos + 1 + offset
        elif offset < 0:
            result["variant_type"] = "acceptor_region"
            result["canonical_site_type"] = "acceptor"
            if strand == "1":
                result["distance_to_canonical"] = -(offset + 1)  # 3354-1del, offset is -1, the POS is 3354-2, it is at the position of start of acceptor at 3354-2
                result["exon_boundary_pos"] = var_pos + 1 - offset
            elif strand == "-1":
                result["distance_to_canonical"] = -(offset + 1) - 1   # "3354-1del", offset is -1, the POS is 3354, it is 1bp upstream of start of acceptor at 3354-1 (pos strand direction, start of the canonical acceptor site)
                result["exon_boundary_pos"] = var_pos - (offset + 1)
        else:
            raise ValueError(f"Invalid offset: {offset}")

    # Pattern for single position SNVs with intronic offset
    # e.g., "116+1G>A", "1730+3A>C", "1897-2A>G"
    snv_pattern = r'(\d+)?(\+|-)(\d+)([ATCG]>[ATCG])'

    # Check for single position SNVs
    snv_match = re.search(snv_pattern, coding_part)
    if snv_match:
        offset = int(snv_match.group(2) + snv_match.group(3))
        if offset > 0:
            result["variant_type"] = "donor_region"
            result["canonical_site_type"] = "donor"
            if strand == "1":
                result["distance_to_canonical"] = offset - 1  # 1730+3A>C, offset is 3, POS is at 1730+3, it is 2bp downstream of donor 1730+1
                result["exon_boundary_pos"] = var_pos - offset
            elif strand == "-1":
                result["distance_to_canonical"] = -offset + 2  # 1730+3A>C, offset is 3, POS is at 1730+3, it is 1bp upstream of donor 1730+2 (pos strand direction, start of the canonical donor site)
                result["exon_boundary_pos"] = var_pos + offset
        elif offset < 0:
            result["variant_type"] = "acceptor_region"
            result["canonical_site_type"] = "acceptor"
            if strand == "1":
                result["distance_to_canonical"] = offset + 2  # 1897-4A>G, offset is -4, POS is at 1897-4, it is 2bp upstream of acceptor 1897-2
                result["exon_boundary_pos"] = var_pos - offset
            elif strand == "-1":
                result["distance_to_canonical"] = - offset - 1  # 1897-4A>G, offset is -4, POS is at 1897-4, it is 3bp downstream of acceptor 1897-1 (pos strand direction, start of the canonical acceptor site)
                result["exon_boundary_pos"] = var_pos + offset
        else:
            raise ValueError(f"Invalid offset: {offset}")
    
    if result['distance_to_canonical'] is None:
        logger.warning(f"No match found for {hgvsc}, this maybe a 5UTR variant")
        return None
    else:
        result["overlapping_canonical_site"] = (result["distance_to_canonical"] + var_length) >= 0
        logger.info(f"The parse result from the HGVSc {hgvsc}, at strand {strand}, var_pos {var_pos}, is {result}")
        return result



def SpliceAI_interpretation(DS_AG, 
                            DS_AL, 
                            DS_DG, 
                            DS_DL, 
                            DP_AG, 
                            DP_AL, 
                            DP_DG, 
                            DP_DL, 
                            transcript_id: str, 
                            transcript_domain_map: dict, 
                            intolerant_domains: set, 
                            consequence = "",
                            cds_pos = None,
                            intron_pos = None,  # 1/8 or 1-3/8 for indels
                            exon_pos = None,  # 1/8 or 1-3/8 for indels
                            hgvsc = None,  # Add HGVSc for distance calculation
                            var_pos = None, 
                            strand = "1",
                            cutoff = 0.8,
                            interpro_entry_map_dict = None,
                            dm_instance = None) -> tuple:
    '''
    Expect to return a tuple of boolean values(is_harmful, length_changing)
    Regarding DP values, negative values means the site is upstream of the variant site,
    and positive values means the site is downstream of the variant site.
    '''
    length_changing = False
    is_harmful = False
    frame_shift = False
    loc_intol_domain = False
    affected_exons = set()
    ten_percent_protein = False
    canonical_loss_offset = 0
    canonical_gain_offset = 0
    logger.debug(f"For {transcript_id}, the SpliceAI predictions are DS_AG: {DS_AG}, DS_AL: {DS_AL}, DS_DG: {DS_DG}, DS_DL: {DS_DL}, DP_AG: {DP_AG}, DP_AL: {DP_AL}, DP_DG: {DP_DG}, DP_DL: {DP_DL}")
    logger.debug(f"The other input arguments are transcript_id: {transcript_id}, intron_pos: {intron_pos}, exon_pos: {exon_pos}, cutoff: {cutoff}")

    # First look at acceptor-gain
    if DS_AG > cutoff and (not "UTR" in consequence):
        length_changing = True
        # The variant is located in the intron and the new cryptic acceptor site is located upstream, which becomes an early acceptor site.
        # The downstream exon is getting longer at the 5' end.
        parse_result = parse_hgvsc_splice_position(hgvsc, strand, var_pos)
        canonical_as_offset = 0
        if parse_result is not None:
            if parse_result["variant_type"] == "acceptor_region":
                canonical_acceptor_site_distance = parse_result["distance_to_canonical"] # negative means var_pos is upstream of the canonical acceptor site, pos value means downstream
                # DP_AG negative means the new acceptor site is upstream of the var_pos
                canonical_as_offset = DP_AG + canonical_acceptor_site_distance
                canonical_gain_offset = canonical_as_offset

        if abs(canonical_as_offset) < 3:
            length_changing = False
        elif abs(canonical_as_offset) % 3 != 0:
            logger.info(f"Given the HGVSc: {hgvsc}, strand {strand}, var_pos {var_pos}, the var_pos to canonical AS distance is {canonical_acceptor_site_distance}. Combined with the DP_AG: {DP_AG}, the canonical AS offset is {canonical_as_offset} (the offset from the canonical AS to the new AS), causing a frameshift effect")
            logger.info(f"The parse result from the HGVSc is {parse_result}")
            # Causing a frameshift effect
            frame_shift = True
            is_harmful = True

        if cds_pos:
            if abs(canonical_as_offset) / int(cds_pos.split('/')[1]) > 0.1:
                # Large truncation relative to cds size
                ten_percent_protein = True
                is_harmful = True

    # Then look at donor-gain
    if DS_DG > cutoff:
        length_changing = True
        # Since the variant is at the upstream of the variant which located on an exon. The early cryptic donor site might get this exon truncated at the 3' end.
        if exon_pos:
            affected_exons.add(int(exon_pos.split('/')[0]))
        
        parse_result = parse_hgvsc_splice_position(hgvsc, strand, var_pos)
        canonical_ds_offset = 0
        if parse_result is not None:
            if parse_result["variant_type"] == "donor_region":
                canonical_donor_site_distance = parse_result["distance_to_canonical"] # negative means var_pos is upstream of the canonical donor site, pos value means downstream
                # DP_DG negative means the new donor site is upstream of the var_pos
                canonical_ds_offset = DP_DG + canonical_donor_site_distance
                canonical_gain_offset = canonical_ds_offset

        if abs(canonical_ds_offset) < 3:
            length_changing = False
        elif abs(canonical_ds_offset) % 3 != 0:
            logger.info(f"Given the HGVSc: {hgvsc}, strand {strand}, var_pos {var_pos}, the var_pos to canonical DS distance is {canonical_donor_site_distance}. Combined with the DP_DG: {DP_DG}, the canonical DS offset is {canonical_ds_offset} (the offset from the canonical DS to the new DS), causing a frameshift effect")
            logger.info(f"The parse result from the HGVSc is {parse_result}")
            # Causing a frameshift effect
            frame_shift = True
            is_harmful = True

        if cds_pos:
            if abs(canonical_ds_offset) / int(cds_pos.split('/')[1]) > 0.1:
                # Large truncation relative to cds size
                ten_percent_protein = True
                is_harmful = True

    # Then look at acceptor-loss
    if DS_AL > cutoff and (not "UTR" in consequence):    
        # For acceptor loss, check if it's disrupting canonical site
        parse_result = parse_hgvsc_splice_position(hgvsc, strand, var_pos)
        canonical_as_offset = 0
        if parse_result is not None:
            if parse_result["variant_type"] == "acceptor_region":
                canonical_acceptor_site_distance = parse_result["distance_to_canonical"] # negative means var_pos is upstream of the canonical acceptor site, pos value means downstream
                # DP_AL negative means the new acceptor site is upstream of the var_pos
                canonical_as_offset = DP_AL + canonical_acceptor_site_distance
                canonical_loss_offset = canonical_as_offset
                logger.info(f"Given the HGVSc: {hgvsc}, strand {strand}, var_pos {var_pos}, the var_pos to canonical AS distance is {canonical_acceptor_site_distance}. Combined with the DP_AL: {DP_AL}, the canonical AS offset is {canonical_as_offset} (the offset from the canonical AS to the predicted loss AS)")
                logger.info(f"The parse result from the HGVSc is {parse_result}")

        if abs(canonical_as_offset) <= 1:
            length_changing = True
            if DS_AG > 0.5:
                # See whether SpliceAI can find a nearby acceptor site that can rescue the lost acceptor site
                canonical_as_offset = 0
                if parse_result is not None:
                    if parse_result["variant_type"] == "acceptor_region":
                        canonical_acceptor_site_distance = parse_result["distance_to_canonical"] # negative means var_pos is upstream of the canonical acceptor site, pos value means downstream
                        # DP_AG negative means the new acceptor site is upstream of the var_pos
                        canonical_as_offset = DP_AG + canonical_acceptor_site_distance
                        canonical_gain_offset = canonical_as_offset

                if abs(canonical_as_offset) < 3:
                    length_changing = False
                elif abs(canonical_as_offset) % 3 != 0:
                    logger.info(f"Given the HGVSc: {hgvsc}, strand {strand}, var_pos {var_pos}, the var_pos to canonical AS distance is {canonical_acceptor_site_distance}. Combined with the DP_AG: {DP_AG}, the canonical AS offset is {canonical_as_offset} (the offset from the canonical AS to the new AS), causing a frameshift effect")
                    logger.info(f"The parse result from the HGVSc is {parse_result}")
                    # Causing a frameshift effect
                    frame_shift = True
                    is_harmful = True

                if abs(canonical_as_offset) >= 30:
                    if intron_pos:
                        affected_exons.add(int(intron_pos.split('/')[0]) + 1)
                    if exon_pos:
                        affected_exons.add(int(exon_pos.split('/')[0]))

                if cds_pos:
                    if abs(canonical_as_offset) / int(cds_pos.split('/')[1]) > 0.1:
                        # Large truncation relative to cds size
                        ten_percent_protein = True
                        is_harmful = True
            else:
                # meaning no nearby acceptor site can rescue the lost acceptor site
                if exon_pos:
                    affected_exons.add(int(exon_pos.split('/')[0]))
                if intron_pos:
                    affected_exons.add(int(intron_pos.split('/')[0]) + 1)
            
        
        # Based on the input annotations, we cannot determine whether:
        # 1. There is a downstream cryptic acceptor site that can rescue some part of the downstream exon
        # 2. Whether the truncated part will cause a frameshift effect
    
    # Then look at donor-loss
    if DS_DL > cutoff:
        parse_result = parse_hgvsc_splice_position(hgvsc, strand, var_pos)
        canonical_ds_offset = 2
        if parse_result is not None:
            if parse_result["variant_type"] == "donor_region":
                canonical_donor_site_distance = parse_result["distance_to_canonical"] # negative means var_pos is upstream of the canonical donor site, pos value means downstream
                # DP_DL negative means the new donor site is upstream of the var_pos
                canonical_ds_offset = DP_DL + canonical_donor_site_distance
                canonical_loss_offset = canonical_ds_offset
                logger.info(f"Given the HGVSc: {hgvsc}, strand {strand}, var_pos {var_pos}, the var_pos to canonical DS distance is {canonical_donor_site_distance}. Combined with the DP_DL: {DP_DL}, the canonical DS offset is {canonical_ds_offset} (the offset from the canonical DS to the predicted loss DS), causing a frameshift effect")
                logger.info(f"The parse result from the HGVSc is {parse_result}")
        if abs(canonical_ds_offset) <= 1:
            length_changing = True
            is_harmful = True
            frame_shift = True
            ten_percent_protein = True
        # Usually the canonical donor site loss is considered as null variant for the transcript.

    # If there are affected exons, we take a look at the transcript domain annotations to see if any of them are intolerant domains.
    if affected_exons:
        for exon in affected_exons:
            if transcript_id in transcript_domain_map and str(exon) in transcript_domain_map[transcript_id]:
                domains = transcript_domain_map[transcript_id][str(exon)]
                func_domains = [domain_path for domain_path in domains if dm_instance.interpret_functionality(domain_path.split(':', 1)[-1], interpro_entry_map_dict) == "Functional"]
                intolerant_affected = set(domains).intersection(intolerant_domains)
                if (len(intolerant_affected) > 0) or (len(func_domains) > 0):
                    loc_intol_domain = True
                    is_harmful = True

    return is_harmful, length_changing, frame_shift, loc_intol_domain, list(affected_exons), ten_percent_protein, canonical_gain_offset, canonical_loss_offset
    



def SpliceVault_interpretation(value, 
                               spliceai_delta,
                               transcript_id, 
                               transcript_domain_map, 
                               intolerant_domains,
                               interpro_entry_map_dict,
                               cds_pos = None,
                               intron_pos = None,
                               exon_pos = None,
                               dm_instance = None) -> tuple:
    '''
    The input value is the top-4 SpliceVault events:
    Top1:CD:-10:2%:Frameshift&Top2:CD:-4:0.6%:Frameshift&Top3:CD:+491:0.05%:Frameshift&Top4:CD:+21:0.006%:Frameshift
    Top1:CD;-51;0.5%;inFrame&Top2:CD;+394;0.1%;inFrame&Top3:CD;+389;0.09%;Frameshift&Top4:CD;+440;0.04%;Frameshift

    Expect to return a tuple of boolean values(is_harmful, length_changing)
    Return False, False if the input value is not a string
    '''
    if not isinstance(value, str):
        logger.debug(f"The input value for SpliceVault interpretation is not a string: {value}, return False, False, False, False, set()")
        return False, False, False, False, set(), False

    logger.info(f"The input transcript_id is {transcript_id}, cds_pos is {cds_pos}, intron_pos is {intron_pos}, exon_pos is {exon_pos}")
    
    events = value.split('&')
    analysis_results = []
    for event in events:
        event_analysis = analyze_splice_event(event, 
                                              transcript_id, 
                                              transcript_domain_map, 
                                              intolerant_domains, 
                                              interpro_entry_map_dict,
                                              cds_pos,
                                              intron_pos, 
                                              exon_pos,
                                              dm_instance)
        analysis_results.append(event_analysis)

    # There are in total 4 events, if all of them are harmful, then return "LoF"
    # Aggregate the fraction of harmful events
    # Aggregate the fraction of not harmful events
    # If the odds ratio is greater than 20, then return "LoF"
    # Otherwise, return "VOUS"
    harmful_frac = sum([result['fraction_of_samples'] for result in analysis_results if result['is_harmful']])
    non_harmful_frac = sum([result['fraction_of_samples'] for result in analysis_results if not result['is_harmful']])
    len_change_frac = sum([result['fraction_of_samples'] for result in analysis_results if result['length_changing']])
    no_len_change_frac = sum([result['fraction_of_samples'] for result in analysis_results if not result['length_changing']])
    frame_shift_frac = sum([result['fraction_of_samples'] for result in analysis_results if result['frame_shift']])
    no_frame_shift_frac = sum([result['fraction_of_samples'] for result in analysis_results if not result['frame_shift']])
    span_intol_domain_frac = sum([result['fraction_of_samples'] for result in analysis_results if result['span_intol_domain']])
    no_span_intol_domain_frac = sum([result['fraction_of_samples'] for result in analysis_results if not result['span_intol_domain']])
    ten_percent_protein_frac = sum([result['fraction_of_samples'] for result in analysis_results if result['ten_percent_protein']])
    no_ten_percent_protein_frac = sum([result['fraction_of_samples'] for result in analysis_results if not result['ten_percent_protein']])
    affected_exon_frac = {e: 0 for e in set([e for result in analysis_results for e in result['affected_exons']])}
    for result in analysis_results:
        for exon in result['affected_exons']:
            affected_exon_frac[exon] += result['fraction_of_samples']
    harmful_odds_ratio = harmful_frac / non_harmful_frac if non_harmful_frac > 0 else 20
    len_change_odds_ratio = len_change_frac / no_len_change_frac if no_len_change_frac > 0 else 20
    frameshift_odds_ratio = frame_shift_frac / no_frame_shift_frac if no_frame_shift_frac > 0 else 20
    span_intol_domain_odds_ratio = span_intol_domain_frac / no_span_intol_domain_frac if no_span_intol_domain_frac > 0 else 20
    ten_percent_protein_odds_ratio = ten_percent_protein_frac / no_ten_percent_protein_frac if no_ten_percent_protein_frac > 0 else 20
    harmful = False
    length_changing = False
    frame_shift = False
    span_intol_domain = False
    ten_percent_protein = False
    affected_exons = [e for e, f in affected_exon_frac.items() if f > 0.3 and spliceai_delta >= 0.5]
    if harmful_odds_ratio >= 20 and harmful_frac > 0.3 and spliceai_delta >= 0.5:
        logger.info(f"For {transcript_id}, there are {harmful_frac} fraction of samples with harmful events and {non_harmful_frac} fraction of samples with non-harmful events, the fraction ratio is {harmful_odds_ratio}, so return LoF, the interpretations are {analysis_results}\n")
        harmful = True
  
    if len_change_odds_ratio >= 20 and len_change_frac > 0.3 and spliceai_delta >= 0.5:
        logger.info(f"For {transcript_id}, there are {len_change_frac} fraction of samples with length changing events and {no_len_change_frac} fraction of samples with no length changing events, the fraction ratio is {len_change_odds_ratio}, so return LoF, the interpretations are {analysis_results}\n")
        length_changing = True

    if span_intol_domain_odds_ratio >= 20 and span_intol_domain_frac > 0.3 and spliceai_delta >= 0.5:
        logger.info(f"For {transcript_id}, there are {span_intol_domain_frac} fraction of samples with span intolerant domain events and {no_span_intol_domain_frac} fraction of samples with no span intolerant domain events, the fraction ratio is {span_intol_domain_odds_ratio}, so return LoF, the interpretations are {analysis_results}\n")
        span_intol_domain = True

    if frameshift_odds_ratio >= 20 and frame_shift_frac > 0.3 and spliceai_delta >= 0.5:
        frame_shift = True
        harmful = True
        logger.info(f"For {transcript_id}, there are {frame_shift_frac} fraction of samples with frame shift events and {no_frame_shift_frac} fraction of samples with no frame shift events, the fraction ratio is {frameshift_odds_ratio}, so return LoF, the interpretations are {analysis_results}\n")

    if ten_percent_protein_odds_ratio >= 20 and ten_percent_protein_frac > 0.3 and spliceai_delta >= 0.5:
        ten_percent_protein = True
        harmful = True
        logger.info(f"For {transcript_id}, there are {ten_percent_protein_frac} fraction of samples with ten percent protein events and {no_ten_percent_protein_frac} fraction of samples with no ten percent protein events, the fraction ratio is {ten_percent_protein_odds_ratio}, so return LoF, the interpretations are {analysis_results}\n")
    return harmful, length_changing, frame_shift, span_intol_domain, affected_exons, ten_percent_protein
    
            
    

def analyze_splice_event(event_str: str, 
                        transcript_id: str,
                        transcript_domain_map: dict,
                        intolerant_domains: set,
                        interpro_entry_map_dict: dict,
                        cds_pos = None,
                        intron_pos = None,  # Format: "2/3" meaning 2nd intron out of 3
                        exon_pos = None,    # Format: "2-3" meaning exons 2 and 3
                        dm_instance: DomainNormalizer = None
                        ) -> dict:
    """
    Analyze a single SpliceVault event to determine if it affects intolerant domains.
    
    Args:
        event_str (str): String describing splice event (e.g. "ES:2-3;12%;Frameshift")
        transcript_id (str): Ensembl transcript ID
        intron_pos (str): String indicating which intron is affected ("current/total")
        exon_pos (str): String indicating which exons are involved ("start-end")
        transcript_domain_map (dict): Mapping of transcripts to exon-domain relationships
                                      {transcript_id: {exon_number: [domain_paths]}}
        intolerant_domains (set): Set of domain paths that are intolerant to structural changes
        
    Returns:
        dict: Dictionary containing:
            - affected_exons (list): List of exon numbers affected
            - affected_domains (list): List of domain paths affected
            - is_lof (bool): Boolean indicating if event likely causes LoF
            - reason (str): String explaining LoF classification
    """
    # Parse event string
    try:
        rank, event_type, pos_str, freq, frame_impact = event_str.split(':')  # e.g., "Top1:CD:-10:2%:Frameshift&Top2:CD:-4:0.6%:Frameshift&Top3:CD:+491:0.05%:Frameshift&Top4:CD:+21:0.006%:Frameshift"
    except ValueError as ve:
        raise ValueError(f"Error splitting event string {event_str}: {ve}")

    logger.debug(f"The input transcript_id is {transcript_id}, event_type is {event_type}, pos_str is {pos_str}, freq is {freq}, frame_impact is {frame_impact}, intron_pos is {intron_pos}, exon_pos is {exon_pos}, whether exon_pos is NA value? {na_value(exon_pos)}, what is the dtype of exon_pos? {type(exon_pos), }, cds_pos is {cds_pos}")
    # We need to translate 2% to 0.02, 0.6% to 0.006, 0.05% to 0.00005, 0.006% to 0.000006
    if '%' in freq:
        freq = float(freq.rstrip('%')) / 100
    elif "#" in freq:
        freq = float(freq.lstrip('#')) / 100
    else:
        raise ValueError(f"The frequency format is not recognized: {freq}, the event string is {event_str}")
    # Initialize results
    affected_exons = set()
    is_harmful = False
    length_changing = False
    frame_shift = False
    span_intol_domain = False
    ten_percent_protein = False
    reason = []

    # Process event based on type
    if event_type == 'ES':
        ten_percent_protein = True
        affected_exons, reason, length_changing, is_harmful, frame_shift, span_intol_domain = process_exon_skipping(
            pos_str, transcript_id, transcript_domain_map, intolerant_domains, interpro_entry_map_dict, dm_instance)
    elif event_type == 'CD':
        affected_exons, reason, length_changing, is_harmful, frame_shift, span_intol_domain = process_cryptic_donor(
            pos_str, transcript_id, transcript_domain_map, intolerant_domains, intron_pos, exon_pos, cds_pos, interpro_entry_map_dict, dm_instance)
    elif event_type == 'CA':
        affected_exons, reason, length_changing, is_harmful, frame_shift, span_intol_domain = process_cryptic_acceptor(
            pos_str, transcript_id, transcript_domain_map, intolerant_domains, intron_pos, exon_pos, cds_pos, interpro_entry_map_dict, dm_instance)
    
    # Look up affected domains and determine if harmful
    affected_domains = set()
    intolerant_affected = set()
    if transcript_id in transcript_domain_map:
        for exon in affected_exons:
            exon_str = str(exon)
            if exon_str in transcript_domain_map[transcript_id]:
                domains = transcript_domain_map[transcript_id][exon_str]
                affected_domains.update(domains)
                intolerant_affected.update(set(domains).intersection(intolerant_domains))
    
    # Event is harmful if it affects intolerant domains or causes frameshift
    if intolerant_affected:
        is_harmful = True
        reason.append('affects_intolerant_domains')
    if frame_impact == 'Frameshift':
        is_harmful = True
        reason.append('frameshift')
    if "large_truncation_relative_to_cds_size" in reason:
        ten_percent_protein = True
    
    return {
        'affected_exons': list(affected_exons),
        'affected_domains': list(affected_domains),
        'intolerant_domains_affected': list(intolerant_affected),
        'is_harmful': is_harmful,
        'length_changing': length_changing,
        'frame_shift': frame_shift,
        'span_intol_domain': span_intol_domain or intolerant_affected,
        "ten_percent_protein": ten_percent_protein,
        'reason': ','.join(reason),
        'event_type': event_type,
        'fraction_of_samples': freq
    }


def process_exon_skipping(pos_str: str, 
                         transcript_id: str, 
                         transcript_domain_map: dict, 
                         intolerant_domains: set,
                         interpro_entry_map_dict: dict,
                         dm_instance: DomainNormalizer = None) -> tuple:
    """Process exon skipping events."""
    affected_exons = set()
    reason = []
    length_changing = True
    frame_shift = False
    is_lof = False
    loc_intol_domain = False
    
    if '-' in pos_str:
        start, end = map(int, pos_str.split('-'))
        affected_exons.update(range(start, end + 1))
    else:
        affected_exons.add(int(pos_str))
    
    # Check for intolerant domains
    for exon in affected_exons:
        exon_str = str(exon)
        if transcript_id in transcript_domain_map and exon_str in transcript_domain_map[transcript_id]:
            domains = transcript_domain_map[transcript_id][exon_str]
            func_domains = [domain_path for domain_path in domains if dm_instance.interpret_functionality(domain_path.split(':', 1)[-1], interpro_entry_map_dict) == "Functional"]
            intolerant_affected = set(domains).intersection(intolerant_domains)
            if (len(intolerant_affected) > 0) or (len(func_domains) > 0):
                reason.append(f'exon_{exon}_overlaps_intolerant_domain')
                is_lof = True
                loc_intol_domain = True
    return affected_exons, reason, length_changing, is_lof, frame_shift, loc_intol_domain


def process_cryptic_donor(pos_str: str, 
                         transcript_id: str,
                         transcript_domain_map: dict,
                         intolerant_domains: set,
                         intron_pos = None,
                         exon_pos = None,
                         cds_pos = None,
                         interpro_entry_map_dict: dict = None,
                         dm_instance: DomainNormalizer = None
                         ) -> tuple:
    """Process cryptic donor events."""
    affected_exons = set()
    reason = []
    length_changing = False
    is_lof = False
    frame_shift = False
    loc_intol_domain = False
    logger.warning(f"The input transcript_id is {transcript_id}, intron_pos is {intron_pos}, exon_pos is {exon_pos}, cds_pos is {cds_pos}, whether intron_pos is NA value? {na_value(intron_pos)}, whether exon_pos is NA value? {na_value(exon_pos)}")
    # Parse position information if available
    if not na_value(intron_pos):
        intron_num, total_introns = map(int, intron_pos.split('/'))
    else:
        intron_num = None
        total_introns = None
    if not na_value(exon_pos):
        exon_num, total_exons = map(int, exon_pos.split('/'))
    else:
        exon_num = None
        total_exons = None

    if not na_value(cds_pos):
        cds_pos, total_cds = cds_pos.split('/')
        if "-" in cds_pos:
            cds_pos = int(cds_pos.split('-')[0]) if cds_pos.split('-')[0] else None
            total_cds = int(total_cds)
        else:
            cds_pos = int(cds_pos)
            total_cds = int(total_cds)
    else:
        cds_pos = None
        total_cds = None
    
    offset = int(pos_str.lstrip('+-'))
    if pos_str.startswith('+'):
        # Downstream cryptic donor (in intron)
        if intron_pos:
            reason.append(f'downstream_cryptic_donor_in_intron_{intron_num}')
        if abs(offset) > 3:
            length_changing = True
        if abs(offset) % 3 != 0:
            frame_shift = True
    else:
        # Upstream cryptic donor (in exon)
        if abs(offset) > 3:  # Ignore very close cryptic sites
            length_changing = True
            if exon_pos:
                affected_exons.add(exon_num)
                reason.append(f'upstream_cryptic_donor_in_exon_{exon_num}')
            elif intron_pos:
                affected_exons.add(intron_num)
                reason.append(f'upstream_cryptic_donor_in_exon_{intron_num}')

        # Check offset fraction if > 10%, then it is likely to cause LoF
        if total_cds:
            if abs(offset) / total_cds > 0.1:
                reason.append('large_truncation_relative_to_cds_size')
                is_lof = True
        
        if abs(offset) % 3 != 0:
            frame_shift = True
    
    # Check for intolerant domains
    for exon in affected_exons:
        exon_str = str(exon)
        if transcript_id in transcript_domain_map and exon_str in transcript_domain_map[transcript_id]:
            domains = transcript_domain_map[transcript_id][exon_str]
            func_domains = [domain_path for domain_path in domains if dm_instance.interpret_functionality(domain_path.split(':', 1)[-1], interpro_entry_map_dict) == "Functional"]
            intolerant_affected = set(domains).intersection(intolerant_domains)
            if (len(intolerant_affected) > 0) or (len(func_domains) > 0):
                reason.append(f'exon_{exon}_overlaps_intolerant_domain')
                is_lof = True
                loc_intol_domain = True
    return affected_exons, reason, length_changing, is_lof, frame_shift, loc_intol_domain


def process_cryptic_acceptor(pos_str: str, 
                           transcript_id: str,
                           transcript_domain_map: dict,
                           intolerant_domains: set,
                           intron_pos = None,
                           exon_pos = None,
                           cds_pos = None,
                           interpro_entry_map_dict: dict = None,
                           dm_instance: DomainNormalizer = None
                           ) -> tuple:
    """Process cryptic acceptor events."""
    affected_exons = set()
    reason = []
    length_changing = False
    is_lof = False
    frame_shift = False
    loc_intol_domain = False
    # Parse position information if available
    if not na_value(intron_pos):
        intron_num, total_introns = map(int, intron_pos.split('/'))
    else:
        intron_num = None
        total_introns = None
    if not na_value(exon_pos):
        exon_num, total_exons = map(int, exon_pos.split('/'))
    else:
        exon_num = None
        total_exons = None
    if not na_value(cds_pos):
        cds_pos, total_cds = cds_pos.split('/')
        if "-" in cds_pos:
            logger.info(f"The cds_pos is {cds_pos}, total_cds is {total_cds}, the cds_pos is a range, we need to split it into two parts")
            try:
                cds_pos = int(cds_pos.split('-')[0])
            except ValueError as ve:
                logger.warning(f"Error splitting cds_pos {cds_pos} and total_cds {total_cds}: {ve}")
                cds_pos = None
            total_cds = int(total_cds)
        else:
            cds_pos = int(cds_pos)
            total_cds = int(total_cds)
    else:
        cds_pos = None
        total_cds = None
    
    offset = int(pos_str.lstrip('+-'))
    if pos_str.startswith('-'):
        # Upstream cryptic acceptor (in intron)
        length_changing = True
        if not na_value(intron_pos):
            affected_exons.add(intron_num + 1)
            reason.append(f'upstream_cryptic_acceptor_in_intron_{intron_num}')
        if abs(offset) % 3 != 0:
            frame_shift = True
    else:
        length_changing = True
        # Downstream cryptic acceptor (in exon)
        if not na_value(exon_pos):
            reason.append(f'downstream_cryptic_acceptor_in_exon_{exon_num}')
        elif not na_value(intron_pos):
            reason.append(f'downstream_cryptic_acceptor_in_exon_{intron_num + 1}')

        if abs(offset) % 3 != 0:
            frame_shift = True
        # Check offset fraction if > 10%, then it is likely to cause LoF
        if total_cds:
            if abs(offset) / total_cds > 0.1:
                reason.append('large_truncation_relative_to_cds_size')
                is_lof = True
    
    # Check for intolerant domains
    for exon in affected_exons:
        exon_str = str(exon)
        if transcript_id in transcript_domain_map and exon_str in transcript_domain_map[transcript_id]:
            domains = transcript_domain_map[transcript_id][exon_str]
            func_domains = [domain_path for domain_path in domains if dm_instance.interpret_functionality(domain_path.split(':', 1)[-1], interpro_entry_map_dict) == "Functional"]
            intolerant_affected = set(domains).intersection(intolerant_domains)
            if (len(intolerant_affected) > 0) or (len(func_domains) > 0):
                reason.append(f'exon_{exon}_overlaps_intolerant_domain')
                is_lof = True
                loc_intol_domain = True
    return affected_exons, reason, length_changing, is_lof, frame_shift, loc_intol_domain



if __name__ == "__main__":
    parser = ap.ArgumentParser(description="Interpret splicing variants")
    parser.add_argument("--anno_table", required=True, help="Annotation table")
    parser.add_argument("--transcript_domain_map", required=True, help="Transcript domain map, pickle file")
    parser.add_argument("--intolerant_domains", required=True, help="Intolerant domains, pickle file")
    parser.add_argument("--interpro_entry_map_pkl", required=True, help="InterPro entry map, pickle file")
    parser.add_argument("--threads", default=12, type=int, help="Number of threads")
    args = parser.parse_args()
    splicing_interpretation(args.anno_table, 
                            args.transcript_domain_map, 
                            args.intolerant_domains, 
                            args.interpro_entry_map_pkl,
                            threads=args.threads)
