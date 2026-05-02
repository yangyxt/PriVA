import pickle
import gzip
import argparse as ap
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
from scipy import stats
from typing import Dict, Set, List
import os
import pandas as pd
import logging
import json
from collections import defaultdict
from statsmodels.stats.multitest import multipletests
import random

import pysam

from protein_domain_mapping import DomainNormalizer
from stat_protein_domain_amscores import nested_defaultdict, _NON_DOMAIN_SOURCES
self_directory = os.path.dirname(os.path.abspath(__file__))


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def test_beta_fit(scores: np.ndarray) -> dict:
    """
    Fit beta distribution to the scores and perform goodness-of-fit test. 
    Deprecated at 2024-12-02.
    
    Args:
        scores: Array of AM scores (should be between 0 and 1)
        
    Returns:
        Dictionary containing fit results and test statistics
    """
    # Fit beta distribution
    alpha, beta, loc, scale = stats.beta.fit(scores)
    
    # Perform Kolmogorov-Smirnov test for goodness of fit
    ks_stat, ks_p = stats.kstest(scores, 'beta', args=(alpha, beta, loc, scale))
    
    return {
        'alpha': alpha,
        'beta': beta,
        'loc': loc,
        'scale': scale,
        'ks_test': {'statistic': ks_stat, 'p_value': ks_p}
    }


def calculate_percentile_differences(query_scores: np.ndarray, reference_scores: np.ndarray) -> dict:
    """
    Calculates differences at biologically meaningful percentiles.
    
    Args:
        query_scores: Array of AM scores for query domain
        reference_scores: Array of AM scores for reference domain
        
    Returns:
        Dictionary of percentile differences
    """
    percentiles = [25, 50, 75]  # Quartiles
    query_percentiles = np.percentile(query_scores, percentiles)
    ref_percentiles = np.percentile(reference_scores, percentiles)
    
    return {
        f'p{p}_diff': float(ref - query)
        for p, ref, query in zip(percentiles, ref_percentiles, query_percentiles)
    }


def build_fixed_equal_weight_composite(ref_arrays: List[np.ndarray], T: int, rng: np.random.Generator) -> np.ndarray:
    """
    Build a single fixed, equal-weight composite of total size T from K reference arrays.
    If an array is shorter than its allocated take, sample with replacement to reach the target.
    """
    K = len(ref_arrays)
    if K == 0:
        raise ValueError("No reference arrays provided for composite.")
    m0 = T // K
    r  = T - m0 * K

    parts = []
    for k, y in enumerate(ref_arrays):
        take = m0 + (1 if k < r else 0)
        y = np.asarray(y)
        if len(y) == 0:
            raise ValueError("One of the reference arrays is empty.")
        if len(y) >= take:
            idx = rng.choice(len(y), size=take, replace=False)
            draw = y[idx]
        else:
            # pad with replacement to reach 'take'
            pad = rng.choice(y, size=take - len(y), replace=True)
            draw = np.concatenate([y, pad])
        parts.append(draw)
    return np.concatenate(parts)


def ks_vs_fixed_composite(
    query_scores: np.ndarray,
    fixed_composite: np.ndarray,
    rng: np.random.Generator,
    cap: int | None = None
) -> dict:
    """
    One-sided KS: tests if the query distribution is more tolerant (left-shifted) than the composite.
    We downsample the composite to max(query_size, 50) samples (or 'cap'), for stability with small queries.
    Returns statistic, pvalue, sample sizes, and median values for effect size estimation.
    """
    query_scores = np.asarray(query_scores)
    n_q = len(query_scores)
    if n_q < 2:
        return {'ks_stat': float('nan'), 'p_value': float('nan'), 'n_query': n_q, 'n_comp': 0, 'median_query': float('nan'), 'median_comp': float('nan')}

    M = max(n_q, 50)
    if cap is not None:
        M = min(M, cap)
    M = min(M, len(fixed_composite))
    replace_flag = M > len(fixed_composite)

    idx = rng.choice(len(fixed_composite), size=M, replace=replace_flag)
    y_mix = fixed_composite[idx]

    # One-sided KS: 'greater' → F_query(x) - F_mix(x) (query CDF larger → more left-shift → more tolerant)
    stat, p = stats.ks_2samp(query_scores, y_mix, alternative='greater')

    return {
        'ks_stat': float(stat),
        'p_value': float(p),
        'n_query': int(n_q),
        'n_comp': int(M),
        'median_query': float(np.median(query_scores)),
        'median_comp': float(np.median(y_mix))
    }



# def analyze_domain_tolerance(query_scores: np.ndarray, 
#                              reference_domains: Dict[str, np.ndarray]) -> dict:
#     """
#     Statistical test for domain tolerance comparison using Fisher's combined probability method.
    
#     Args:
#         query_scores: Array of scores for the query domain
#         reference_domains: Dictionary mapping reference domain names to their score arrays
        
#     Returns:
#         Dictionary containing test results or None if testing fails
#     """
#     if len(reference_domains) == 0:
#         return None
        
#     results = []
#     p_values = []  # Collect p-values for Fisher's method
    
#     # Randomly pick 100 domains from the reference domains
#     ref_names = list(reference_domains.keys())
#     ref_names = random.sample(ref_names, min(100, len(ref_names)))
#     for ref_name in ref_names:
#         ref_scores = reference_domains[ref_name]
#         # Perform KS test
#         stat, p_value = stats.ks_2samp(
#             query_scores,
#             ref_scores,
#             alternative='greater'
#         )
		
# 		# For alternative='greater',
# 		# CDF of ref scores < CDF of query scores if p-value < 0.05, so PDF of ref scores > PDF of query scores, meaning query domains are more tolerant
# 		# CDF of ref scores >= CDF of query scores if p-value > 0.05, so PDF of ref scores <= PDF of query scores, meaning query domains are more or equally intolerant
        
#         # Calculate effect sizes (percentile differences)
#         percentile_diffs = calculate_percentile_differences(query_scores, ref_scores)
        
#         results.append({
#             'reference': ref_name,
#             'statistic': stat,
#             'p_value': p_value,
#             'ref_sample_size': len(ref_scores),
#             'query_sample_size': len(query_scores),
#             **percentile_diffs
#         })
        
#         p_values.append(p_value)
    
#     # Apply Fisher's method
#     # chi_square = -2 * sum(ln(p)) follows chi-square distribution with 2k degrees of freedom
#     # where k is the number of tests being combined
#     chi_square_stat = -2 * np.sum(np.log(np.array(p_values)))
#     significant_proportion = sum(p_value < 0.05 for p_value in p_values) / len(p_values)
#     combined_p_value = stats.chi2.sf(chi_square_stat, df=2*len(p_values))
#     assert isinstance(combined_p_value, float), f"combined_p_value is not a float: {combined_p_value} for query scores {query_scores} and reference domains {reference_domains}"
    
#     return {
#         'is_more_tolerant': combined_p_value < 0.05,  # You can adjust this threshold
#         'fisher_combined_p_value': combined_p_value,
#         'fisher_chi_square_stat': chi_square_stat,
#         'significant_proportion': significant_proportion,
#         'n_references': len(ref_names),
#         'individual_results': results
#     }


def visualize_single_domain(args):
    """
    Process a single domain's distribution data.
    
    Args:
        args: Tuple containing (domain_path, scores, output_dir)
        
    Returns:
        Dictionary with basic statistics or None if processing fails
    """
    domain_path, scores, output_dir = args

    output_path = os.path.join(output_dir, f"distribution_{domain_path}.png")
    if os.path.exists(output_path):
        logger.debug(f"Skipping {domain_path} because the plot {output_path} already exists")
        return None

    try:    
        if len(scores) < 8:
            return None
            
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
            
        # Create plot
        plt.figure(figsize=(8, 6))
        
        # Histogram
        plt.hist(scores, bins='auto', density=True, alpha=0.7)
        plt.title(f'Score Distribution: {domain_path}')
        plt.xlabel('AM Score')
        plt.ylabel('Density')
        
        # Add basic statistics as text
        text_results = (
            f"Distribution Statistics:\n"
            f"n={len(scores)}\n"
            f"mean={np.mean(scores):.3f}\n"
            f"std={np.std(scores):.3f}"
        )
        plt.text(0.98, 0.02, text_results, 
                fontsize=8, ha='right', va='bottom', 
                transform=plt.gca().transAxes)
        
        # Save plot
        plt.tight_layout()
        plt.savefig(output_path, bbox_inches='tight', dpi=300)
        plt.close()
        
        return {
            'domain': domain_path,
            'n_samples': len(scores),
            'mean': np.mean(scores),
            'std': np.std(scores)
        }
        
    except Exception as e:
        logger.error(f"Error processing domain {domain_path}: {str(e)}")
        return None


def collect_domain_data(d, output_dir: str, min_variants: int, path=[]) -> List:
    """
    Collect domain data for parallel processing.
    
    Args:
        d: Nested dictionary containing domain data
        output_dir: Directory to save outputs
        min_variants: Minimum number of variants required
        path: Current path in the nested dictionary (used for recursion)
        
    Returns:
        List of tuples containing (domain_path, scores, output_dir)
    """
    domain_data = []
    for k, v in d.items():
        if k == 'max_distribution':
            domain_path = ':'.join(path)
            if len(v) >= min_variants:
                logger.debug(f"Processing domain {domain_path} with {len(v)} variants")
                domain_data.append((domain_path, v, output_dir))
            else:
                logger.info(f"Skipping domain {domain_path} due to insufficient variants ({len(v)} < {min_variants})")
        elif isinstance(v, dict) and not 'distribution' in k:
            current_path = path + [k]
            domain_data.extend(collect_domain_data(v, output_dir, min_variants, current_path))
    return domain_data


def identify_functional_domain(domain_path, dm_instance=None, func_map_dict = None, func_pred_dict = None, interpro_map_dict = None):
    domain_name = domain_path.split(':', 1)[-1] # Remove ENSG ids
    interpro_entry = dm_instance.query_interpro_entry_vep_anno(domain_name, interpro_map_dict).get("interpro_entries", None)
    if interpro_entry:
        interpro_id = interpro_entry[0][0]
        interpro_type = interpro_entry[0][1]
        go_terms = interpro_entry[0][4]
        functional = False if all([go_term[1] != "molecular_function" for go_term in go_terms]) else True
    else:
        return None
    if functional and interpro_type in ["Conserved_Site", "Binding_site", "Active_site", "Domain", "PTM"]:
        return "Functional"
    elif not functional and interpro_type in ["Domain", "Conserved_Site", "Binding_site", "Active_site", "PTM"]:
        return "Non-functional"
    else:
        return None


def visualize_domain_distribution(scores_dict: dict, output_dir: str, threads: int, assembly: str='hg19') -> Dict[str, np.ndarray]:
    """
    Process and visualize domain distributions, returning a dictionary of domain scores.
    
    Args:
        scores_dict: Dictionary containing domain scores
        output_dir: Directory to save outputs
        threads: Number of threads for parallel processing
        
    Returns:
        Dictionary mapping domain paths to their score arrays
    """
    # Process all domains if no intolerant domains file provided
    domain_data = []
    domain_scores = {}  # New dictionary to store scores

    output_dir = os.path.join(output_dir, 'am_domain_distributions')
    os.makedirs(output_dir, exist_ok=True)
    
    for gene_id, gene_data in scores_dict.items():
        for db_name, db_data in gene_data.items():
            # Recursively collect domain data
            domain_data.extend(collect_domain_data(
                                                    db_data, 
                                                    output_dir, 
                                                    6,  # Minimum variants required
                                                    path=[gene_id, db_name]
                                                ))
        
    # Also store scores in flat dictionary for later use
    for path, scores, _ in domain_data:
        domain_scores[path] = scores

    # Dump the domain scores to a pickle file
    with open(os.path.join(output_dir, f'domain_amscores.{assembly}.pkl'), 'wb') as f:
        pickle.dump(domain_scores, f)

    logger.info(f"Found {len(domain_data)} domains to process, saved their AM scores to {os.path.join(output_dir, f'domain_amscores.{assembly}.pkl')}")
    
    # Process domains in parallel
    n_cores = min(threads, len(domain_data))  # Don't use more cores than tasks
    logger.info(f"Processing {len(domain_data)} domains using {n_cores} cores")
    with mp.Pool(n_cores) as pool:
        results = list(pool.imap_unordered(visualize_single_domain, domain_data))
    
    # Filter out None results
    results = [r for r in results if r is not None]
    logger.info(f"Successfully visualized distribution of {len(results)} domains")
    
    # Save results to TSV if any successful results
    if results:
        results_df = pd.DataFrame(results)
        results_df.to_csv(os.path.join(output_dir, f'distribution_test_results.{assembly}.tsv'), 
                        sep='\t', index=False)
    
    return domain_scores


def process_domain_ks(args):
    """
    Args: (domain_path, scores, fixed_composite, cap, seed)
    Returns: (domain_path, result_dict)
    """
    domain_path, scores, fixed_composite, cap, seed = args
    rng = np.random.default_rng(seed)
    res = ks_vs_fixed_composite(scores, fixed_composite, rng, cap=cap)
    return (domain_path, res)



# def process_domain_tolerance(args):
#     """
#     Process a single domain's tolerance analysis.
    
#     Args:
#         args: Tuple containing (domain_path, scores, ref_domain_scores)
        
#     Returns:
#         Tuple of (domain_path, result) or None if processing fails
#     """
#     try:
#         domain_path, scores, ref_domain_scores = args
#         result = analyze_domain_tolerance(scores, ref_domain_scores)
#         if result is not None:
#             # vis_result = {k:v for k, v in result.items() if k != 'individual_results'}
#             if result['fisher_combined_p_value'] > 0.05:
#                 logger.debug(f"Completed tolerance analysis for domain {domain_path}, the result is {result}\n")
#             return (domain_path, result)
#         else:
#             logger.warning(f"Failed to complete tolerance analysis for domain {domain_path}")
#     except Exception as e:
#         logger.error(f"Error processing domain {domain_path}: {str(e)}")
#     return None


def process_domain_dict_for_mapping(d: dict, current_path: List[str], pp_map: dict, exon_map: dict) -> None:
    """
    Process a single level of the domain dictionary to extract transcript-domain relationships
    using protein position ranges and exon-domain associations.
    
    Args:
        d: Current level of the nested dictionary
        current_path: List of strings representing the current domain path
        pp_map: Dict to store {transcript: {domain_path: (aa_min, aa_max)}}
        exon_map: Dict to store {transcript: {exon_str: set(domain_paths)}}
    """
    for k, v in d.items():
        if isinstance(v, dict):
            if 'distribution' in v:
                domain_path = ':'.join(current_path + [k])
                for transcript_id, tx_data in v.items():
                    if not (isinstance(transcript_id, str) and transcript_id.startswith('ENST')):
                        continue
                    # Handle both new format (dict with exons/aa_min/aa_max) and legacy (set)
                    if isinstance(tx_data, dict) and 'exons' in tx_data:
                        aa_min = tx_data.get('aa_min', float('inf'))
                        aa_max = tx_data.get('aa_max', 0)
                        if aa_min <= aa_max:
                            pp_map.setdefault(transcript_id, {})[domain_path] = (aa_min, aa_max)
                        for exon_id in tx_data['exons']:
                            exon_idx = exon_id.split("/")[0]
                            exon_map.setdefault(transcript_id, {}).setdefault(exon_idx, set()).add(domain_path)
                    elif isinstance(tx_data, set):
                        # Legacy format: bare exon set, no protein positions
                        for exon_id in tx_data:
                            exon_idx = exon_id.split("/")[0]
                            exon_map.setdefault(transcript_id, {}).setdefault(exon_idx, set()).add(domain_path)
            else:
                process_domain_dict_for_mapping(v, current_path + [k], pp_map, exon_map)


def generate_transcript_pp_domain_map(scores_dict: dict) -> dict:
    """
    Generate a mapping from transcript_id to protein-position-based domain ranges
    and exon-domain associations.
    
    Args:
        scores_dict: Nested dictionary containing domain data with transcript info
        
    Returns:
        Dictionary with structure:
        {
            'transcript_id': {
                'pp_domains': [(aa_start, aa_end, 'domain_path'), ...],  # sorted by aa_start
                'exon_domains': {
                    '1': ['domain_path1', 'domain_path2'],
                    '2': ['domain_path3'],
                }
            }
        }
    """
    pp_map = {}   # {transcript: {domain_path: (aa_min, aa_max)}}
    exon_map = {} # {transcript: {exon_str: set(domain_paths)}}
    
    for gene_id, gene_data in scores_dict.items():
        for db_name, domain_data in gene_data.items():
            process_domain_dict_for_mapping(domain_data, [gene_id, db_name], pp_map, exon_map)
    
    result = {}
    all_transcripts = set(pp_map.keys()) | set(exon_map.keys())
    for tx in all_transcripts:
        entry = {}
        # Build sorted pp_domains list
        if tx in pp_map:
            entry['pp_domains'] = sorted(
                [(aa_min, aa_max, dp) for dp, (aa_min, aa_max) in pp_map[tx].items()],
                key=lambda x: x[0]
            )
        else:
            entry['pp_domains'] = []
        # Build exon_domains with sorted lists
        if tx in exon_map:
            entry['exon_domains'] = {
                exon_idx: sorted(list(domain_paths))
                for exon_idx, domain_paths in exon_map[tx].items()
            }
        else:
            entry['exon_domains'] = {}
        result[tx] = entry
    
    logger.info(f"Generated mapping for {len(result)} transcripts "
                f"({sum(1 for v in result.values() if v['pp_domains'])} with protein position data)")
    return result


def _collect_gnomad_leaves(d, result, path):
    """Walk scores_dict nested structure and collect domain_path -> gnomAD_distribution array."""
    for k, v in d.items():
        if k == 'gnomAD_distribution':
            domain_path = ':'.join(path)
            result[domain_path] = v
        elif isinstance(v, dict) and 'distribution' not in k:
            _collect_gnomad_leaves(v, result, path + [k])


def extract_domain_metadata(vcf_path, target_domains=None, region=None):
    """Parse VEP-annotated AlphaMissense VCF to extract per-domain metadata.

    For each domain_path, collects gene_symbol, aa_start, aa_end from CSQ fields.
    Uses the same domain hierarchy parsing as stat_protein_domain_amscores.

    Args:
        vcf_path: Path to AlphaMissense VEP-annotated VCF (.vcf.gz)
        target_domains: Optional set of domain_paths to restrict extraction.
        region: Optional chromosome/region string for region-based fetching.

    Returns:
        dict: {domain_path: {'gene_symbol': str, 'aa_start': int, 'aa_end': int}}
    """
    vcf = pysam.VariantFile(vcf_path)
    csq_format = str(vcf.header).split('Format: ')[-1].strip('>"').split('|')

    domains_idx = csq_format.index('DOMAINS')
    gene_idx = csq_format.index('Gene')
    symbol_idx = csq_format.index('SYMBOL')
    protein_pos_idx = csq_format.index('Protein_position')
    feature_type_idx = csq_format.index('Feature_type')
    feature_idx = csq_format.index('Feature')

    domain_meta = {}
    n_records = 0

    records = vcf.fetch(region) if region else vcf
    for record in records:
        n_records += 1
        if n_records % 5_000_000 == 0:
            logger.info(f"[DOMAIN_META] Processed {n_records:,} VCF records, "
                        f"metadata for {len(domain_meta)} domains so far")

        csq_value = record.info.get('CSQ')
        if not csq_value:
            continue

        if isinstance(csq_value, tuple):
            transcript_annotations = csq_value
        else:
            transcript_annotations = csq_value.split(',')

        for csq_entry in transcript_annotations:
            fields = csq_entry.split('|')

            if fields[feature_type_idx] != 'Transcript' or not fields[feature_idx].startswith('ENST'):
                continue

            domains_str = fields[domains_idx]
            ensg_id = fields[gene_idx]
            symbol = fields[symbol_idx]
            protein_pos_str = fields[protein_pos_idx]

            if not domains_str or not ensg_id or not protein_pos_str:
                continue

            try:
                pp_parts = protein_pos_str.split('/')
                # Handle range format: "2-3" -> "2"
                aa_pos = int(pp_parts[0].split('-')[0])
                prot_length = int(pp_parts[1]) if len(pp_parts) > 1 else None
            except (ValueError, IndexError):
                continue

            # Parse domains — same logic as stat_protein_domain_amscores._parse_domain_hierarchy
            for domain_entry in domains_str.split('&'):
                parts = domain_entry.split(':')
                db_name = parts[0]
                if db_name in _NON_DOMAIN_SOURCES:
                    continue
                for i in range(1, len(parts)):
                    hierarchy = parts[1:i+1]
                    domain_path = ':'.join([ensg_id, db_name] + hierarchy)

                    if target_domains is not None and domain_path not in target_domains:
                        continue

                    if domain_path not in domain_meta:
                        domain_meta[domain_path] = {
                            'gene_symbol': symbol,
                            'aa_start': aa_pos,
                            'aa_end': aa_pos,
                            'prot_length': prot_length or 0,
                        }
                    else:
                        meta = domain_meta[domain_path]
                        if aa_pos < meta['aa_start']:
                            meta['aa_start'] = aa_pos
                        if aa_pos > meta['aa_end']:
                            meta['aa_end'] = aa_pos
                        if prot_length is not None and prot_length > meta['prot_length']:
                            meta['prot_length'] = prot_length

    vcf.close()
    logger.info(f"[DOMAIN_META] Finished: {n_records:,} records processed, "
                f"metadata for {len(domain_meta)} domains")
    return domain_meta


def _extract_metadata_worker(args):
    """Multiprocessing worker: extract domain metadata for a single chromosome."""
    vcf_path, chrom, target_domains = args
    return extract_domain_metadata(vcf_path, target_domains=target_domains, region=chrom)


def _merge_domain_metadata(target, source):
    """Merge per-chromosome domain metadata dicts in place."""
    for dp, meta in source.items():
        if dp not in target:
            target[dp] = meta
        else:
            t = target[dp]
            if meta['aa_start'] < t['aa_start']:
                t['aa_start'] = meta['aa_start']
            if meta['aa_end'] > t['aa_end']:
                t['aa_end'] = meta['aa_end']
            if meta.get('prot_length', 0) > t.get('prot_length', 0):
                t['prot_length'] = meta['prot_length']


def extract_domain_metadata_parallel(vcf_path, target_domains=None, threads=1):
    """Parallel wrapper around extract_domain_metadata, split by chromosome.

    Falls back to single-threaded when threads <= 1.
    """
    if threads <= 1:
        return extract_domain_metadata(vcf_path, target_domains=target_domains)

    vcf = pysam.VariantFile(vcf_path)
    chroms = list(vcf.header.contigs)
    vcf.close()
    n_workers = min(threads, len(chroms))
    logger.info(f"[DOMAIN_META] Extracting metadata from {len(chroms)} chromosomes "
                f"with {n_workers} workers")

    with mp.Pool(n_workers) as pool:
        results = pool.map(_extract_metadata_worker,
                           [(vcf_path, c, target_domains) for c in chroms])

    merged = {}
    for result in results:
        _merge_domain_metadata(merged, result)

    logger.info(f"[DOMAIN_META] Parallel extraction complete: "
                f"metadata for {len(merged)} domains")
    return merged


def filter_coldspot_overlap(intolerant_domains, domain_metadata, hcseeker_spots_tsv,
                            overlap_threshold=0.5):
    """Remove domains with bidirectional coldspot overlap > threshold.

    A domain is removed if, for any coldspot in the same gene, EITHER:
      - overlap / coldspot_length > threshold  (domain covers most of a coldspot), OR
      - overlap / domain_length   > threshold  (coldspot covers most of the domain).

    Args:
        intolerant_domains: set of domain_path strings
        domain_metadata: dict from extract_domain_metadata()
        hcseeker_spots_tsv: path to HC_spots.{assembly}.tsv
        overlap_threshold: fraction threshold for either direction

    Returns:
        (filtered_set, excluded_set, summary_dict)
    """
    spots_df = pd.read_csv(hcseeker_spots_tsv, sep='\t')
    coldspots = spots_df[spots_df['type'] == 'coldspot']

    gene_coldspots = defaultdict(list)
    for _, row in coldspots.iterrows():
        gene_coldspots[row['gene']].append((int(row['aa_start_pos']), int(row['aa_end_pos'])))
    logger.info(f"[FILTER1:COLDSPOT] Loaded {len(coldspots)} coldspots across {len(gene_coldspots)} genes")

    excluded = set()
    summary = {}
    n_evaluated = 0

    for domain_path in intolerant_domains:
        meta = domain_metadata.get(domain_path)
        if meta is None:
            summary[domain_path] = {
                'coldspot_max_frac_of_coldspot': 0.0,
                'coldspot_max_frac_of_domain': 0.0,
                'removed': False,
            }
            continue

        gene_sym = meta['gene_symbol']
        dom_start = meta['aa_start']
        dom_end = meta['aa_end']

        if gene_sym not in gene_coldspots:
            summary[domain_path] = {
                'coldspot_max_frac_of_coldspot': 0.0,
                'coldspot_max_frac_of_domain': 0.0,
                'removed': False,
            }
            continue

        n_evaluated += 1
        max_frac_of_cs = 0.0
        max_frac_of_dom = 0.0
        dom_len = dom_end - dom_start + 1

        for cs_start, cs_end in gene_coldspots[gene_sym]:
            overlap_start = max(dom_start, cs_start)
            overlap_end = min(dom_end, cs_end)
            if overlap_start > overlap_end:
                continue
            overlap_len = overlap_end - overlap_start + 1
            coldspot_len = cs_end - cs_start + 1
            frac_of_cs = overlap_len / coldspot_len
            frac_of_dom = overlap_len / dom_len
            if frac_of_cs > max_frac_of_cs:
                max_frac_of_cs = frac_of_cs
            if frac_of_dom > max_frac_of_dom:
                max_frac_of_dom = frac_of_dom

        removed = (max_frac_of_cs > overlap_threshold) or (max_frac_of_dom > overlap_threshold)
        summary[domain_path] = {
            'coldspot_max_frac_of_coldspot': max_frac_of_cs,
            'coldspot_max_frac_of_domain': max_frac_of_dom,
            'removed': removed,
        }
        if removed:
            excluded.add(domain_path)

    filtered = intolerant_domains - excluded

    logger.info(f"[FILTER1:COLDSPOT] Input intolerant domains: {len(intolerant_domains)}")
    logger.info(f"[FILTER1:COLDSPOT] Domains evaluated (gene matched to HCSeeker): {n_evaluated}")
    logger.info(f"[FILTER1:COLDSPOT] Domains removed (bidirectional overlap > {overlap_threshold}): {len(excluded)}")
    logger.info(f"[FILTER1:COLDSPOT] Output intolerant domains: {len(filtered)}")
    if excluded:
        top_removed = sorted(excluded, key=lambda d: max(
            summary[d]['coldspot_max_frac_of_coldspot'],
            summary[d]['coldspot_max_frac_of_domain']), reverse=True)[:5]
        for d in top_removed:
            logger.info(f"[FILTER1:COLDSPOT] Removed: {d} "
                        f"(frac_of_coldspot={summary[d]['coldspot_max_frac_of_coldspot']:.3f}, "
                        f"frac_of_domain={summary[d]['coldspot_max_frac_of_domain']:.3f})")

    return filtered, excluded, summary


def filter_benign_clinvar(intolerant_domains, clinvar_vcf):
    """Remove domains containing >=1 high-confidence benign ClinVar missense variant.

    High-confidence: CLNSIG contains 'Benign' (not just Likely_benign),
    CLNREVSTAT >= 2 review stars, consequence includes missense_variant.
    Uses ClinVar VCF's CSQ DOMAINS field to map variants to domain_paths directly.

    Args:
        intolerant_domains: set of domain_path strings
        clinvar_vcf: path to clinvar.{assembly}.vep.vcf.gz

    Returns:
        (filtered_set, excluded_set, summary_dict)
    """
    REVIEW_2STAR = {
        'practice_guideline',
        'reviewed_by_expert_panel',
        'criteria_provided,_multiple_submitters,_no_conflicts',
    }

    vcf = pysam.VariantFile(clinvar_vcf)
    csq_format = str(vcf.header).split('Format: ')[-1].strip('>"').split('|')

    domains_idx = csq_format.index('DOMAINS')
    gene_idx = csq_format.index('Gene')
    consequence_idx = csq_format.index('Consequence')
    feature_type_idx = csq_format.index('Feature_type')
    feature_idx = csq_format.index('Feature')

    domain_benign_count = defaultdict(int)

    for record in vcf:
        clnsig = record.info.get('CLNSIG')
        if clnsig is None:
            continue
        if isinstance(clnsig, tuple):
            clnsig = ','.join(clnsig)
        if 'Benign' not in clnsig:
            continue

        clnrevstat = record.info.get('CLNREVSTAT')
        if clnrevstat is None:
            continue
        if isinstance(clnrevstat, tuple):
            clnrevstat = ','.join(clnrevstat)
        if clnrevstat not in REVIEW_2STAR:
            continue

        csq_value = record.info.get('CSQ')
        if not csq_value:
            continue

        if isinstance(csq_value, tuple):
            transcript_annotations = csq_value
        else:
            transcript_annotations = csq_value.split(',')

        for csq_entry in transcript_annotations:
            fields = csq_entry.split('|')

            if fields[feature_type_idx] != 'Transcript' or not fields[feature_idx].startswith('ENST'):
                continue
            if 'missense_variant' not in fields[consequence_idx]:
                continue

            domains_str = fields[domains_idx]
            ensg_id = fields[gene_idx]

            if not domains_str or not ensg_id:
                continue

            for domain_entry in domains_str.split('&'):
                parts = domain_entry.split(':')
                db_name = parts[0]
                for i in range(1, len(parts)):
                    hierarchy = parts[1:i+1]
                    domain_path = ':'.join([ensg_id, db_name] + hierarchy)
                    if domain_path in intolerant_domains:
                        domain_benign_count[domain_path] += 1

    vcf.close()

    excluded = set()
    summary = {}

    for domain_path in intolerant_domains:
        n_benign = domain_benign_count.get(domain_path, 0)
        removed = n_benign >= 1
        summary[domain_path] = {'n_benign_missense': n_benign, 'removed': removed}
        if removed:
            excluded.add(domain_path)

    filtered = intolerant_domains - excluded

    logger.info(f"[FILTER2:BENIGN] Input intolerant domains: {len(intolerant_domains)}")
    logger.info(f"[FILTER2:BENIGN] Domains with >=1 HC benign missense: {len(excluded)}")
    logger.info(f"[FILTER2:BENIGN] Domains removed: {len(excluded)}")
    logger.info(f"[FILTER2:BENIGN] Output intolerant domains: {len(filtered)}")
    if excluded:
        top_removed = sorted(excluded, key=lambda d: domain_benign_count[d], reverse=True)[:5]
        for d in top_removed:
            logger.info(f"[FILTER2:BENIGN] Removed: {d} (n_benign={domain_benign_count[d]})")

    return filtered, excluded, summary


def filter_common_missense(intolerant_domains, scores_dict, af_threshold=0.05):
    """Remove domains containing >=1 missense variant with AF_grpmax_joint > af_threshold.

    Reads gnomAD_distribution from scores_dict (numpy arrays of AF values per domain).

    Args:
        intolerant_domains: set of domain_path strings
        scores_dict: loaded domain scores pkl containing gnomAD_distribution per domain leaf
        af_threshold: AF threshold (default 0.05, BA1 criterion)

    Returns:
        (filtered_set, excluded_set, summary_dict)
    """
    gnomad_data = {}
    _collect_gnomad_leaves(scores_dict, gnomad_data, path=[])

    excluded = set()
    summary = {}
    n_with_gnomad = 0

    for domain_path in intolerant_domains:
        gnomad_arr = gnomad_data.get(domain_path)
        if gnomad_arr is None or len(gnomad_arr) == 0:
            summary[domain_path] = {'n_common_missense': 0, 'max_af': 0.0, 'removed': False}
            continue

        n_with_gnomad += 1
        n_common = int(np.sum(gnomad_arr > af_threshold))
        max_af = float(np.max(gnomad_arr))
        removed = n_common >= 1
        summary[domain_path] = {'n_common_missense': n_common, 'max_af': max_af, 'removed': removed}
        if removed:
            excluded.add(domain_path)

    filtered = intolerant_domains - excluded

    logger.info(f"[FILTER3:COMMON] Input intolerant domains: {len(intolerant_domains)}")
    logger.info(f"[FILTER3:COMMON] Domains with gnomAD data: {n_with_gnomad} / {len(intolerant_domains)}")
    logger.info(f"[FILTER3:COMMON] Domains removed (any AF > {af_threshold}): {len(excluded)}")
    logger.info(f"[FILTER3:COMMON] Output intolerant domains: {len(filtered)}")
    if excluded:
        top_removed = sorted(excluded, key=lambda d: summary[d]['max_af'], reverse=True)[:5]
        for d in top_removed:
            logger.info(f"[FILTER3:COMMON] Removed: {d} (max_AF={summary[d]['max_af']:.4f})")

    return filtered, excluded, summary


def filter_whole_protein_family(intolerant_domains, domain_metadata,
                                dm_instance, interpro_map_dict,
                                coverage_threshold=0.9, min_protein_length=200):
    """Remove Family-type InterPro entries spanning ≥coverage_threshold of proteins ≥min_protein_length AA.

    Proteins shorter than min_protein_length typically harbor a single structural
    domain (Ramírez-Sánchez et al. 2016; Xu & Nussinov 1998), so Family entries
    on such proteins are retained as legitimate single-domain annotations.

    Args:
        intolerant_domains: set of domain_path strings
        domain_metadata: dict {domain_path: {aa_start, aa_end, prot_length, ...}}
        dm_instance: DomainNormalizer instance for InterPro lookups
        interpro_map_dict: InterPro entry mapping dict
        coverage_threshold: minimum domain/protein coverage to trigger exclusion (default 0.9)
        min_protein_length: proteins shorter than this are exempt (default 200)

    Returns:
        (filtered_set, excluded_set, summary_dict)
    """
    excluded = set()
    summary = {}

    for domain_path in intolerant_domains:
        meta = domain_metadata.get(domain_path)
        if meta is None or meta.get('prot_length', 0) == 0:
            summary[domain_path] = {'interpro_type': None, 'coverage': None,
                                    'prot_length': None, 'removed': False}
            continue

        # Resolve InterPro type
        domain_name = domain_path.split(':', 1)[-1]  # strip ENSG prefix
        interpro_result = dm_instance.query_interpro_entry_vep_anno(
            domain_name, interpro_map_dict)
        interpro_type = None
        if interpro_result is not None:
            entries = interpro_result.get('interpro_entries')
            if entries:
                interpro_type = entries[0][1]

        prot_len = meta['prot_length']
        coverage = (meta['aa_end'] - meta['aa_start'] + 1) / prot_len

        removed = (interpro_type == 'Family'
                   and coverage >= coverage_threshold
                   and prot_len >= min_protein_length)

        summary[domain_path] = {'interpro_type': interpro_type,
                                'coverage': round(coverage, 4),
                                'prot_length': prot_len, 'removed': removed}
        if removed:
            excluded.add(domain_path)

    filtered = intolerant_domains - excluded

    n_family = sum(1 for s in summary.values() if s.get('interpro_type') == 'Family')
    logger.info(f"[FILTER4:FAMILY] Input intolerant domains: {len(intolerant_domains)}")
    logger.info(f"[FILTER4:FAMILY] Family-type entries found: {n_family}")
    logger.info(f"[FILTER4:FAMILY] Domains removed "
                f"(Family + coverage≥{coverage_threshold} + prot≥{min_protein_length}): "
                f"{len(excluded)}")
    logger.info(f"[FILTER4:FAMILY] Output intolerant domains: {len(filtered)}")
    if excluded:
        top_removed = sorted(excluded,
                             key=lambda d: summary[d].get('coverage', 0),
                             reverse=True)[:5]
        for d in top_removed:
            s = summary[d]
            logger.info(f"[FILTER4:FAMILY] Removed: {d} "
                        f"(coverage={s['coverage']}, prot_length={s['prot_length']})")

    return filtered, excluded, summary


def analyze_domain_data(pickle_file: str,
                       output_dir: str,
                       threads: int = 12,
                       fdr_threshold: float = 0.05,
                       assembly: str = 'hg19',
                       am_vcf: str = None,
                       hcseeker_spots_tsv: str = None,
                       clinvar_vcf: str = None,
                       coldspot_overlap_threshold: float = 0.5,
                       common_af_threshold: float = 0.05):
    """
    Analyze domain data from the pickle file and apply post-KS filters.

    Args:
        pickle_file: Path to the pickle file containing scores_dict
        output_dir: Directory to save analysis outputs
        threads: Number of threads for parallel processing
        fdr_threshold: FDR threshold for significance after correction
        assembly: Assembly version (hg19 or hg38)
        am_vcf: Path to VEP-annotated AlphaMissense VCF (needed for Filter 1)
        hcseeker_spots_tsv: Path to HC_spots.{assembly}.tsv (None = skip Filter 1)
        clinvar_vcf: Path to clinvar.{assembly}.vep.vcf.gz (None = skip Filter 2)
        coldspot_overlap_threshold: Overlap fraction threshold for Filter 1 (default 0.5)
        common_af_threshold: AF threshold for Filter 3 (default 0.05)
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load the data
    with open(pickle_file, 'rb') as f:
        scores_dict = pickle.load(f)

    # Generate transcript-protein position-domain mapping
    transcript_map = generate_transcript_pp_domain_map(scores_dict)
    
    # Save the mapping to both JSON and pickle files
    mapping_json = os.path.join(output_dir, f'transcript_pp_domain_mapping.{assembly}.json')
    mapping_pickle = os.path.join(output_dir, f'transcript_pp_domain_mapping.{assembly}.pkl')
    
    with open(mapping_json, 'w') as f:
        json.dump(transcript_map, f, indent=2)
    with open(mapping_pickle, 'wb') as f:
        pickle.dump(transcript_map, f)
    
    logger.info(f"Saved transcript-pp-domain mapping to {mapping_json} and {mapping_pickle}")
    
    # Continue with existing analysis
    domain_scores = visualize_domain_distribution(scores_dict, output_dir, threads, assembly)
    dm_instance = DomainNormalizer()

    interpro_map_pickle = os.path.join(os.path.dirname(self_directory), 'data', 'InterPro', 'Interpro_entry_mapping.pkl.gz')
    interpro_map_dict = pickle.load(gzip.open(interpro_map_pickle, 'rb'))

    functional_map = os.path.join(os.path.dirname(self_directory), 'data', 'InterPro', 'curated_InterPro_func_domains.tsv.gz')
    func_map_df = pd.read_table(functional_map, low_memory=False)
    func_map_dict = dict(zip(func_map_df['IPR_ID'], func_map_df['Molecular_Function_GO_Terms']))
    func_pred_dict = dict(zip(func_map_df['IPR_ID'], func_map_df['Functionality_Assessment']))
    
    # Pick out the reference (curated functional) domains
    ref_domain_scores = {
        domain: domain_s
        for domain, domain_s in domain_scores.items()
        if identify_functional_domain(domain, dm_instance, func_map_dict, func_pred_dict, interpro_map_dict) == "Functional"
    }
    logger.info(f"Found {len(ref_domain_scores)} curated functional domains/sites as references")

    # Prepare the list of query domains (exclude the curated references)
    query_domains = [(d, s) for d, s in domain_scores.items() if d not in ref_domain_scores]
    logger.info(f"Found {len(query_domains)} query domains to test against the composite reference")

    if len(ref_domain_scores) == 0 or len(query_domains) == 0:
        logger.warning("No references or no query domains found; skipping KS analysis.")
        return

    # --- Build a single fixed equal-weight composite reference ONCE ---
    # Choose composite size T: median of query sizes, capped to save memory
    query_sizes = np.array([len(s) for _, s in query_domains], dtype=int)
    if len(query_sizes) == 0:
        logger.warning("No query sizes found; skipping KS analysis.")
        return
    T_default = int(np.max(query_sizes))
    T_cap = 200000
    T = max(50, min(T_default, T_cap))  # keep it reasonable
    seed_global = 42  # make reproducible
    rng_global = np.random.default_rng(seed_global)

    ref_arrays = list(ref_domain_scores.values())
    fixed_composite = build_fixed_equal_weight_composite(ref_arrays, T=T, rng=rng_global)
    logger.info(f"Built fixed equal-weight composite reference of size T={len(fixed_composite)} from {len(ref_arrays)} references")

    # --- Parallel KS over query domains ---
    # We downsample the composite per-query to match n_query (capped at T)
    cap = T
    n_tasks = len(query_domains)
    n_cores = min(threads, n_tasks) if n_tasks > 0 else 1
    logger.info(f"Processing {n_tasks} query domains using {n_cores} cores for composite KS tests")

    ks_results = {}
    tasks = []
    for i, (domain_path, scores) in enumerate(query_domains):
        # domain-specific seed for deterministic subsampling
        seed = (seed_global + i) & 0xFFFFFFFF
        tasks.append((domain_path, scores, fixed_composite, cap, seed))

    if n_tasks > 0:
        with mp.Pool(n_cores) as pool:
            for res in pool.imap_unordered(process_domain_ks, tasks):
                domain_path, result = res
                ks_results[domain_path] = result

    if not ks_results:
        logger.warning("No KS results were produced.")
        return

    # --- Collect p-values and apply BH FDR across query domains ---
    domains = list(ks_results.keys())
    raw_pvals = [ks_results[d]['p_value'] for d in domains]

    rejected, pvals_bh, _, _ = multipletests(raw_pvals, alpha=fdr_threshold, method='fdr_bh') # Benjamini-Hochberg procedure

    # --- Write results table ---
    df_rows = []
    for d, adj_p, rej in zip(domains, pvals_bh, rejected):
        row = {
            'domain': d,
            'ks_stat': ks_results[d]['ks_stat'],
            'p_value': ks_results[d]['p_value'],
            'fdr_corrected_p_value': float(adj_p),
            # significant = "more tolerant than composite" at FDR
            'is_more_tolerant': bool(rej),
            'n_query': ks_results[d]['n_query'],
            'n_comp_used': ks_results[d]['n_comp'],
            'composite_T': int(len(fixed_composite)),
            'median_query': ks_results[d]['median_query'],
            'median_composite': ks_results[d]['median_comp']
        }
        df_rows.append(row)

    results_df = pd.DataFrame(df_rows)
    out_tsv = os.path.join(output_dir, f'domain_tolerance_analysis.{assembly}.tsv')
    results_df.to_csv(out_tsv, sep='\t', index=False)
    logger.info(f"Composite KS analysis results saved to {out_tsv}")
    logger.info(f"Completed KS analysis for {len(results_df)} domains; "
                f"found {int(np.sum(rejected))} domains significantly more tolerant than composite (FDR < {fdr_threshold})")

    # --- Define "intolerant" set for downstream: those NOT significantly more tolerant + all curated refs ---
    not_more_tolerant = set(results_df.loc[~results_df['is_more_tolerant'], 'domain'].unique().tolist())
    analysis_intolerant_domains = not_more_tolerant | set(ref_domain_scores.keys())
    n_ks_intolerant = len(analysis_intolerant_domains)
    logger.info(f"[KS] Intolerant set (not-tolerant + curated refs): {n_ks_intolerant}")

    # Keep a copy of the pre-filter set for the summary TSV
    ks_intolerant_raw = set(analysis_intolerant_domains)

    # --- Post-KS filters (each filter is independent and optional) ---
    filter_summaries = {}
    excluded_sets = {}

    # --- Load/cache domain metadata (needed for Filters 1 and 4) ---
    domain_metadata = None
    if am_vcf is not None:
        domain_meta_cache = os.path.join(output_dir, f'domain_metadata.{assembly}.pkl')
        if os.path.exists(domain_meta_cache):
            logger.info(f"[DOMAIN_META] Loading cached metadata from {domain_meta_cache}")
            with open(domain_meta_cache, 'rb') as f:
                domain_metadata = pickle.load(f)
        else:
            logger.info(f"[DOMAIN_META] Extracting from AM VCF: {am_vcf}")
            domain_metadata = extract_domain_metadata_parallel(am_vcf, threads=threads)
            with open(domain_meta_cache, 'wb') as f:
                pickle.dump(domain_metadata, f)
            logger.info(f"[DOMAIN_META] Cached {len(domain_metadata)} domains to {domain_meta_cache}")

    # Filter 1: HCSeeker coldspot overlap
    if hcseeker_spots_tsv is not None:
        if domain_metadata is None:
            logger.warning("[FILTER1:COLDSPOT] Skipped: --am_vcf required for domain AA range extraction")
        else:
            analysis_intolerant_domains, excl, summ = filter_coldspot_overlap(
                analysis_intolerant_domains, domain_metadata, hcseeker_spots_tsv,
                coldspot_overlap_threshold)
            filter_summaries['coldspot'] = summ
            excluded_sets['coldspot'] = excl

    # Filter 2: High-confidence benign ClinVar missense
    if clinvar_vcf is not None:
        analysis_intolerant_domains, excl, summ = filter_benign_clinvar(
            analysis_intolerant_domains, clinvar_vcf)
        filter_summaries['benign'] = summ
        excluded_sets['benign'] = excl

    # Filter 3: Common missense (gnomAD AF > threshold) — always runs
    analysis_intolerant_domains, excl, summ = filter_common_missense(
        analysis_intolerant_domains, scores_dict, common_af_threshold)
    filter_summaries['common'] = summ
    excluded_sets['common'] = excl

    # Filter 4: Whole-protein Family classifiers — requires domain_metadata
    if domain_metadata is not None:
        analysis_intolerant_domains, excl, summ = filter_whole_protein_family(
            analysis_intolerant_domains, domain_metadata,
            dm_instance, interpro_map_dict)
        filter_summaries['family'] = summ
        excluded_sets['family'] = excl
    else:
        logger.warning("[FILTER4:FAMILY] Skipped: --am_vcf required")

    # --- Summary logging ---
    n_f1 = len(excluded_sets.get('coldspot', set()))
    n_f2 = len(excluded_sets.get('benign', set()))
    n_f3 = len(excluded_sets.get('common', set()))
    n_f4 = len(excluded_sets.get('family', set()))
    logger.info(f"[SUMMARY] Total domains analyzed: {len(domain_scores)}")
    logger.info(f"[SUMMARY] Curated functional references: {len(ref_domain_scores)}")
    logger.info(f"[SUMMARY] DAS intolerant (post-KS): {n_ks_intolerant}")
    logger.info(f"[SUMMARY] Removed by Filter 1 (coldspot): {n_f1}")
    logger.info(f"[SUMMARY] Removed by Filter 2 (benign): {n_f2}")
    logger.info(f"[SUMMARY] Removed by Filter 3 (common): {n_f3}")
    logger.info(f"[SUMMARY] Removed by Filter 4 (whole-protein family): {n_f4}")
    logger.info(f"[SUMMARY] Final intolerant domains: {len(analysis_intolerant_domains)}")

    # --- Save filter summary TSV ---
    # Resolve domain metadata for aa coordinates
    dm_for_summary = domain_metadata or {}
    if filter_summaries:
        summary_rows = []
        for dp in sorted(ks_intolerant_raw):
            meta = dm_for_summary.get(dp, {})
            row = {
                'domain_path': dp,
                'gene_symbol': meta.get('gene_symbol', ''),
                'aa_start': meta.get('aa_start', np.nan),
                'aa_end': meta.get('aa_end', np.nan),
                'in_original_intolerant': True,
            }
            if 'coldspot' in filter_summaries:
                s = filter_summaries['coldspot'].get(dp, {})
                row['coldspot_max_frac_of_coldspot'] = s.get('coldspot_max_frac_of_coldspot', np.nan)
                row['coldspot_max_frac_of_domain'] = s.get('coldspot_max_frac_of_domain', np.nan)
                row['coldspot_removed'] = s.get('removed', False)
            if 'benign' in filter_summaries:
                s = filter_summaries['benign'].get(dp, {})
                row['n_benign_missense'] = s.get('n_benign_missense', np.nan)
                row['benign_removed'] = s.get('removed', False)
            if 'common' in filter_summaries:
                s = filter_summaries['common'].get(dp, {})
                row['n_common_missense'] = s.get('n_common_missense', np.nan)
                row['max_af'] = s.get('max_af', np.nan)
                row['common_removed'] = s.get('removed', False)
            if 'family' in filter_summaries:
                s = filter_summaries['family'].get(dp, {})
                row['interpro_type'] = s.get('interpro_type', '')
                row['domain_coverage'] = s.get('coverage', np.nan)
                row['prot_length'] = s.get('prot_length', np.nan)
                row['family_removed'] = s.get('removed', False)
            row['final_intolerant'] = dp in analysis_intolerant_domains
            summary_rows.append(row)

        summary_df = pd.DataFrame(summary_rows)
        summary_tsv = os.path.join(output_dir, f'domain_filter_results.{assembly}.tsv')
        summary_df.to_csv(summary_tsv, sep='\t', index=False)
        logger.info(f"Saved filter summary ({len(summary_df)} domains) to {summary_tsv}")

    # --- Save final intolerant domain set ---
    intolerant_domains_pickle = os.path.join(output_dir, f'all_intolerant_domains.{assembly}.pkl')
    with open(intolerant_domains_pickle, 'wb') as f:
        pickle.dump(analysis_intolerant_domains, f)
    logger.info(f"Saved final intolerant domains (N={len(analysis_intolerant_domains)}) to {intolerant_domains_pickle}")

    # --- Save excluded domain sets per filter ---
    if excluded_sets:
        excluded_pickle = os.path.join(output_dir, f'excluded_domains.{assembly}.pkl')
        with open(excluded_pickle, 'wb') as f:
            pickle.dump(excluded_sets, f)
        logger.info(f"Saved excluded domain sets to {excluded_pickle}")



if __name__ == '__main__':
    parser = ap.ArgumentParser(description='Analyze domain data from the pickle file.')
    parser.add_argument('--pickle_file', required=True, help='Path to the pickle file containing scores_dict')
    parser.add_argument('--output_dir', required=True, help='Directory to save analysis outputs')
    parser.add_argument('--threads', type=int, default=62, help='Number of threads for parallel processing (default: 12)')
    parser.add_argument('--assembly', type=str, default='hg19', help='Assembly version, either hg19 or hg38')
    parser.add_argument('--am_vcf', type=str, default=None,
                        help='Path to VEP-annotated AlphaMissense VCF (required for Filter 1: coldspot overlap)')
    parser.add_argument('--hcseeker_spots_tsv', type=str, default=None,
                        help='Path to HC_spots.{assembly}.tsv for Filter 1 (None = skip)')
    parser.add_argument('--clinvar_vcf', type=str, default=None,
                        help='Path to clinvar.{assembly}.vep.vcf.gz for Filter 2 (None = skip)')
    parser.add_argument('--coldspot_overlap_threshold', type=float, default=0.5,
                        help='Overlap fraction threshold for Filter 1 (default: 0.5)')
    parser.add_argument('--common_af_threshold', type=float, default=0.05,
                        help='AF threshold for Filter 3 common missense (default: 0.05)')
    parser.add_argument('--fdr_threshold', type=float, default=0.05,
                        help='FDR threshold for KS test significance (default: 0.05)')
    args = parser.parse_args()

    analyze_domain_data(args.pickle_file,
                        args.output_dir,
                        args.threads,
                        fdr_threshold=args.fdr_threshold,
                        assembly=args.assembly,
                        am_vcf=args.am_vcf,
                        hcseeker_spots_tsv=args.hcseeker_spots_tsv,
                        clinvar_vcf=args.clinvar_vcf,
                        coldspot_overlap_threshold=args.coldspot_overlap_threshold,
                        common_af_threshold=args.common_af_threshold)
