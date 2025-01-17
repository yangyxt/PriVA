import os
import pickle
import pandas as pd
import numpy as np
import argparse as ap
import logging
import uuid
import multiprocessing as mp
from functools import partial


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def na_value(value):
    if isinstance(value, float) and np.isnan(value):
        return True
    return True if value in [np.nan, "NaN", "nan", "na", "NA", "NAN", None, "None", "none", "", ".", "-"] else False


def splicing_interpretation(anno_table, 
                            transcript_domain_map_pkl, 
                            intolerant_domains_pkl, 
                            spliceai_cutoff = 0.8, 
                            threads = 12):
    """Process annotation table in parallel to interpret splicing effects"""
    anno_df = pd.read_table(anno_table, low_memory=False)

    intolerant_domains = pickle.load(open(intolerant_domains_pkl, 'rb'))
    transcript_domain_map = pickle.load(open(transcript_domain_map_pkl, 'rb'))
    
    # Create partial function with fixed arguments
    process_row = partial(splicing_altering_per_row, 
                         transcript_domain_map=transcript_domain_map,
                         intolerant_domains=intolerant_domains,
                         spliceai_cutoff=spliceai_cutoff)
    
    # Get number of CPU cores (leave some free for other processes)
    logger.info(f"Processing {len(anno_df)} variants using {threads} cores")
    
    # Process rows in parallel
    with mp.Pool(threads) as pool:
        results = pool.map(process_row, [row for _, row in anno_df.iterrows()])
    
    # Unpack results
    splicing_lof, splicing_len_changing = zip(*results)
    
    # Add results to dataframe
    anno_df['splicing_lof'] = splicing_lof
    anno_df['splicing_len_changing'] = splicing_len_changing
    
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
    intron_pos = None if na_value(row['INTRON']) else row['INTRON']
    exon_pos = None if na_value(row['EXON']) else row['EXON']
    cds_pos = None if na_value(row["CDS_position"]) else row["CDS_position"]
    logger.info(f"The input transcript_id is {transcript_id}, splicevault_events is {splicevault_events}, intron_pos is {intron_pos}, exon_pos is {exon_pos}, cds_pos is {cds_pos}")

    # First get SpliceVault interpretation
    splicevault_lof, splicevault_len_changing = SpliceVault_interpretation( splicevault_events, 
                                                                            transcript_id, 
                                                                            transcript_domain_map, 
                                                                            intolerant_domains,
                                                                            cds_pos,
                                                                            intron_pos,
                                                                            exon_pos )
    
    # Then get SpliceAI interpretation
    DP_AG = float(row['SpliceAI_pred_DP_AG'])
    DP_AL = float(row['SpliceAI_pred_DP_AL'])
    DP_DG = float(row['SpliceAI_pred_DP_DG'])
    DP_DL = float(row['SpliceAI_pred_DP_DL'])
    DS_AG = float(row['SpliceAI_pred_DS_AG'])
    DS_AL = float(row['SpliceAI_pred_DS_AL'])
    DS_DG = float(row['SpliceAI_pred_DS_DG'])
    DS_DL = float(row['SpliceAI_pred_DS_DL'])
    spliceai_lof, spliceai_len_changing = SpliceAI_interpretation(DS_AG, 
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
                                                                  cds_pos,
                                                                  intron_pos,
                                                                  exon_pos,
                                                                  spliceai_cutoff)
    return splicevault_lof or spliceai_lof, splicevault_len_changing or spliceai_len_changing
    


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
                            cds_pos = None,
                            intron_pos = None, 
                            exon_pos = None,
                            cutoff = 0.8) -> tuple:
    '''
    Expect to return a tuple of boolean values(is_harmful, length_changing)
    Regarding DP values, negative values means the site is upstream of the variant site,
    and positive values means the site is downstream of the variant site.
    '''
    length_changing = False
    is_harmful = False
    affected_exons = set()

    logger.debug(f"For {transcript_id}, the SpliceAI predictions are DS_AG: {DS_AG}, DS_AL: {DS_AL}, DS_DG: {DS_DG}, DS_DL: {DS_DL}, DP_AG: {DP_AG}, DP_AL: {DP_AL}, DP_DG: {DP_DG}, DP_DL: {DP_DL}")
    logger.info(f"The other input arguments are transcript_id: {transcript_id}, intron_pos: {intron_pos}, exon_pos: {exon_pos}, cutoff: {cutoff}")

    # First look at acceptor-gain
    if DS_AG > cutoff:
        if DP_AG < 0 and intron_pos:
            length_changing = True
            # The variant is located in the intron and the new cryptic acceptor site is located upstream, which becomes an early acceptor site.
            # The downstream exon is getting longer at the 5' end.
            affected_exons.add(int(intron_pos.split('/')[0]) + 1)
            if abs(DP_AG) % 3 != 0:
                # Causing a frameshift effect
                is_harmful = True
            if cds_pos:
                if abs(DP_AG) / int(cds_pos.split('/')[1]) > 0.1:
                    # Large truncation relative to cds size
                    is_harmful = True
        else:
            # We cannot determine whether the new acceptor site is upstream or downstream of the existing acceptor site. Hence not able to predict its effect.
            pass

    # Then look at donor-gain
    if DS_DG > cutoff:
        if DP_DG < 0 and exon_pos:
            length_changing = True
            # Since the variant is at the upstream of the variant which located on an exon. The early cryptic donor site might get this exon truncated at the 3' end.
            affected_exons.add(int(exon_pos.split('/')[0]))
            if abs(DP_DG) % 3 != 0:
                # Causing a frameshift effect
                is_harmful = True
            if cds_pos:
                if abs(DP_DG) / int(cds_pos.split('/')[1]) > 0.1:
                    # Large truncation relative to cds size
                    is_harmful = True
        else:
            # We cannot determine whether the new donor site is upstream or downstream of the existing donor site. Hence not able to predict its effect.
            pass

    # Then look at acceptor-loss
    if DS_AL > cutoff:
        length_changing = True
        # Based on the input annotations, we cannot determine whether:
        # 1. There is a downstream cryptic acceptor site that can rescue some part of the downstream exon
        # 2. Whether the truncated part will cause a frameshift effect
    
    # Then look at donor-loss
    if DS_DL > cutoff:
        length_changing = True
        is_harmful = True
        # Usually the canonical donor site loss is considered as null variant for the transcript.

    # If there are affected exons, we take a look at the transcript domain annotations to see if any of them are intolerant domains.
    if affected_exons:
        for exon in affected_exons:
            if transcript_id in transcript_domain_map and str(exon) in transcript_domain_map[transcript_id]:
                domains = transcript_domain_map[transcript_id][str(exon)]
                intolerant_affected = set(domains).intersection(intolerant_domains)
                if intolerant_affected:
                    is_harmful = True

    return is_harmful, length_changing
    



def SpliceVault_interpretation(value, 
                               transcript_id, 
                               transcript_domain_map, 
                               intolerant_domains,
                               cds_pos = None,
                               intron_pos = None,
                               exon_pos = None) -> tuple:
    '''
    The input value is the top-4 SpliceVault events:
    Top1:CD:-10:2%:Frameshift&Top2:CD:-4:0.6%:Frameshift&Top3:CD:+491:0.05%:Frameshift&Top4:CD:+21:0.006%:Frameshift
    Top1:CD;-51;0.5%;inFrame&Top2:CD;+394;0.1%;inFrame&Top3:CD;+389;0.09%;Frameshift&Top4:CD;+440;0.04%;Frameshift

    Expect to return a tuple of boolean values(is_harmful, length_changing)
    '''
    if not isinstance(value, str):
        logger.debug(f"The input value for SpliceVault interpretation is not a string: {value}, return False, False")
        return False, False

    logger.info(f"The input transcript_id is {transcript_id}, cds_pos is {cds_pos}, intron_pos is {intron_pos}, exon_pos is {exon_pos}")
    
    events = value.split('&')
    analysis_results = []
    for event in events:
        event_analysis = analyze_splice_event(event, 
                                              transcript_id, 
                                              transcript_domain_map, 
                                              intolerant_domains, 
                                              cds_pos,
                                              intron_pos, 
                                              exon_pos)
        analysis_results.append(event_analysis)

    # There are in total 4 events, if all of them are harmful, then return "LoF"
    if all([result['is_harmful'] for result in analysis_results]):
        logger.info(f"All 4 events are harmful for {transcript_id}, the interpretations are {analysis_results}\n")
        return True, True
    else:
        # Aggregate the fraction of harmful events
        # Aggregate the fraction of not harmful events
        # If the odds ratio is greater than 10, then return "LoF"
        # Otherwise, return "VOUS"
        harmful_frac = sum([result['fraction_of_samples'] for result in analysis_results if result['is_harmful']])
        non_harmful_frac = sum([result['fraction_of_samples'] for result in analysis_results if not result['is_harmful']])
        odds_ratio = harmful_frac / non_harmful_frac
        if odds_ratio >= 20:
            logger.info(f"For {transcript_id}, there are {harmful_frac} fraction of samples with harmful events and {non_harmful_frac} fraction of samples with non-harmful events, the fraction ratio is {odds_ratio}, so return LoF, the interpretations are {analysis_results}\n")
            return True, True
        else:
            logger.info(f"For {transcript_id}, there are {harmful_frac} fraction of samples with harmful events and {non_harmful_frac} fraction of samples with non-harmful events, the fraction ratio is {odds_ratio}, so return VOUS, the interpretations are {analysis_results}\n")
            if all([result['length_changing'] for result in analysis_results]):
                return False, True
            else:
                return False, False
            
    

def analyze_splice_event(event_str: str, 
                        transcript_id: str,
                        transcript_domain_map: dict,
                        intolerant_domains: set,
                        cds_pos = None,
                        intron_pos = None,  # Format: "2/3" meaning 2nd intron out of 3
                        exon_pos = None,    # Format: "2-3" meaning exons 2 and 3
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
    reason = []

    # Process event based on type
    if event_type == 'ES':
        affected_exons, reason, length_changing, is_harmful = process_exon_skipping(
            pos_str, transcript_id, transcript_domain_map, intolerant_domains)
    elif event_type == 'CD':
        affected_exons, reason, length_changing, is_harmful = process_cryptic_donor(
            pos_str, transcript_id, transcript_domain_map, intolerant_domains, intron_pos, exon_pos, cds_pos)
    elif event_type == 'CA':
        affected_exons, reason, length_changing, is_harmful = process_cryptic_acceptor(
            pos_str, transcript_id, transcript_domain_map, intolerant_domains, intron_pos, exon_pos, cds_pos)
    
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
    
    return {
        'affected_exons': list(affected_exons),
        'affected_domains': list(affected_domains),
        'intolerant_domains_affected': list(intolerant_affected),
        'is_harmful': is_harmful,
        'length_changing': length_changing,
        'reason': ','.join(reason),
        'event_type': event_type,
        'fraction_of_samples': freq
    }


def process_exon_skipping(pos_str: str, 
                         transcript_id: str, 
                         transcript_domain_map: dict, 
                         intolerant_domains: set) -> tuple:
    """Process exon skipping events."""
    affected_exons = set()
    reason = []
    length_changing = True
    is_lof = True
    
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
            intolerant_affected = set(domains).intersection(intolerant_domains)
            if intolerant_affected:
                reason.append(f'exon_{exon}_overlaps_intolerant_domain')
                is_lof = True
     
    return affected_exons, reason, length_changing, is_lof


def process_cryptic_donor(pos_str: str, 
                         transcript_id: str,
                         transcript_domain_map: dict,
                         intolerant_domains: set,
                         intron_pos = None,
                         exon_pos = None,
                         cds_pos = None
                         ) -> tuple:
    """Process cryptic donor events."""
    affected_exons = set()
    reason = []
    length_changing = False
    is_lof = False

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
    
    # Check for intolerant domains
    for exon in affected_exons:
        exon_str = str(exon)
        if transcript_id in transcript_domain_map and exon_str in transcript_domain_map[transcript_id]:
            domains = transcript_domain_map[transcript_id][exon_str]
            intolerant_affected = set(domains).intersection(intolerant_domains)
            if intolerant_affected:
                reason.append(f'exon_{exon}_overlaps_intolerant_domain')
                is_lof = True
    
    return affected_exons, reason, length_changing, is_lof


def process_cryptic_acceptor(pos_str: str, 
                           transcript_id: str,
                           transcript_domain_map: dict,
                           intolerant_domains: set,
                           intron_pos = None,
                           exon_pos = None,
                           cds_pos = None
                           ) -> tuple:
    """Process cryptic acceptor events."""
    affected_exons = set()
    reason = []
    length_changing = False
    is_lof = False

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
    else:
        # Downstream cryptic acceptor (in exon)
        if not na_value(exon_pos):
            reason.append(f'downstream_cryptic_acceptor_in_exon_{exon_num}')
        elif not na_value(intron_pos):
            reason.append(f'downstream_cryptic_acceptor_in_exon_{intron_num + 1}')

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
            intolerant_affected = set(domains).intersection(intolerant_domains)
            if intolerant_affected:
                reason.append(f'exon_{exon}_overlaps_intolerant_domain')
                is_lof = True
    
    return affected_exons, reason, length_changing, is_lof



if __name__ == "__main__":
    parser = ap.ArgumentParser(description="Interpret splicing variants")
    parser.add_argument("--anno_table", required=True, help="Annotation table")
    parser.add_argument("--transcript_domain_map", required=True, help="Transcript domain map, pickle file")
    parser.add_argument("--intolerant_domains", required=True, help="Intolerant domains, pickle file")
    parser.add_argument("--threads", default=12, type=int, help="Number of threads")
    args = parser.parse_args()
    splicing_interpretation(args.anno_table, 
                            args.transcript_domain_map, 
                            args.intolerant_domains, 
                            threads=args.threads)
