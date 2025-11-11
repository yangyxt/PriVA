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

from protein_domain_mapping import DomainNormalizer
from stat_protein_domain_amscores import nested_defaultdict
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
    We downsample the composite to match the query size (or 'cap'), for fairness and stability.
    Returns statistic, pvalue, and sample sizes, plus a simple effect size (median diff).
    """
    query_scores = np.asarray(query_scores)
    n_q = len(query_scores)
    if n_q < 2:
        return {'ks_stat': np.nan, 'p_value': np.nan, 'n_query': n_q, 'n_comp': 0, 'median_diff': np.nan}

    M = n_q
    if cap is not None:
        M = min(M, cap)
    M = min(M, len(fixed_composite))
    replace_flag = M > len(fixed_composite)

    idx = rng.choice(len(fixed_composite), size=M, replace=replace_flag)
    y_mix = fixed_composite[idx]

    # One-sided KS: 'greater' → F_query(x) - F_mix(x) (query CDF larger → more left-shift → more tolerant)
    stat, p = stats.ks_2samp(query_scores, y_mix, alternative='greater')

    median_diff = float(np.median(y_mix) - np.median(query_scores))  # positive → query more tolerant
    return {
        'ks_stat': float(stat),
        'p_value': float(p),
        'n_query': int(n_q),
        'n_comp': int(M),
        'median_diff': median_diff
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


def process_domain_dict_for_mapping(d: dict, current_path: List[str], transcript_exon_map: dict) -> None:
    """
    Process a single level of the domain dictionary to extract transcript-exon-domain relationships.
    
    Args:
        d: Current level of the nested dictionary
        current_path: List of strings representing the current domain path
        transcript_exon_map: Dictionary to store the mapping results
    """
    for k, v in d.items():
        if isinstance(v, dict):
            # Check if this level contains transcript information
            if 'distribution' in v:
                # This level contains domain information
                domain_path = ':'.join(current_path + [k])
                for transcript_id, exon_set in v.items():
                    if isinstance(transcript_id, str) and transcript_id.startswith('ENST'):
                        if isinstance(exon_set, set):  # This is the exon set
                            # Map each exon individually
                            for exon_id in exon_set:
                                exon_idx = exon_id.split("/")[0]
                                transcript_exon_map[transcript_id][exon_idx].add(domain_path)
            else:
                # Continue building the domain path
                process_domain_dict_for_mapping(v, current_path + [k], transcript_exon_map)


def generate_transcript_exon_domain_map(scores_dict: dict) -> Dict[str, Dict[str, List[str]]]:
    """
    Generate a mapping from transcript_id to individual exon indices to domain paths.
    
    Args:
        scores_dict: Nested dictionary containing domain data with transcript and exon information
        
    Returns:
        Dictionary with structure:
        {
            'transcript_id': {
                '1': ['domain_path1', 'domain_path2'],
                '2': ['domain_path1', 'domain_path3'],
                '3': ['domain_path2']
            }
        }
    """
    transcript_exon_map = defaultdict(lambda: defaultdict(set))
    
    # Process the scores dictionary
    for gene_id, gene_data in scores_dict.items():
        for db_name, domain_data in gene_data.items():
            process_domain_dict_for_mapping(domain_data, [gene_id, db_name], transcript_exon_map)
    
    # Convert sets to lists for easier handling
    result = {}
    for transcript_id, exon_map in transcript_exon_map.items():
        result[transcript_id] = {
            exon_idx: sorted(list(domain_paths))
            for exon_idx, domain_paths in exon_map.items()
        }
    
    logger.info(f"Generated mapping for {len(result)} transcripts")
    return result


def analyze_domain_data(pickle_file: str, 
                       output_dir: str,
                       threads: int = 12,
                       fdr_threshold: float = 0.05,
                       assembly: str = 'hg19'):
    """
    Analyze domain data from the pickle file.
    
    Args:
        pickle_file: Path to the pickle file containing scores_dict
        output_dir: Directory to save analysis outputs
        intolerant_domains_tsv: Optional path to TSV file containing intolerant domains
        mechanism_tsv: Path to TSV file containing mechanism analysis results
        threads: Number of threads for parallel processing
        fdr_threshold: FDR threshold for significance after correction
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load the data
    with open(pickle_file, 'rb') as f:
        scores_dict = pickle.load(f)
    
    # Generate transcript-exon-domain mapping
    transcript_map = generate_transcript_exon_domain_map(scores_dict)
    
    # Save the mapping to both JSON and pickle files
    mapping_json = os.path.join(output_dir, f'transcript_exon_domain_mapping.{assembly}.json')
    mapping_pickle = os.path.join(output_dir, f'transcript_exon_domain_mapping.{assembly}.pkl')
    
    with open(mapping_json, 'w') as f:
        json.dump(transcript_map, f, indent=2)
    with open(mapping_pickle, 'wb') as f:
        pickle.dump(transcript_map, f)
    
    logger.info(f"Saved transcript-exon-domain mapping to {mapping_json} and {mapping_pickle}")
    
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
    T_default = int(np.median(query_sizes))
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
            'median_diff': ks_results[d]['median_diff']  # + if query more tolerant
        }
        df_rows.append(row)

    results_df = pd.DataFrame(df_rows)
    out_tsv = os.path.join(output_dir, f'domain_composite_ks_analysis.{assembly}.tsv')
    results_df.to_csv(out_tsv, sep='\t', index=False)
    logger.info(f"Composite KS analysis results saved to {out_tsv}")
    logger.info(f"Completed KS analysis for {len(results_df)} domains; "
                f"found {int(np.sum(rejected))} domains significantly more tolerant than composite (FDR < {fdr_threshold})")

    # --- Define "intolerant" set for downstream: those NOT significantly more tolerant + all curated refs ---
    # (Your earlier logic merged these with reference domains.)
    not_more_tolerant = set(results_df.loc[~results_df['is_more_tolerant'], 'domain'].unique().tolist())
    analysis_intolerant_domains = not_more_tolerant | set(ref_domain_scores.keys())

    intolerant_domains_pickle = os.path.join(output_dir, f'all_intolerant_domains.{assembly}.pkl')
    with open(intolerant_domains_pickle, 'wb') as f:
        pickle.dump(analysis_intolerant_domains, f)

    logger.info(f"Saved merged intolerant domains (N={len(analysis_intolerant_domains)}) to {intolerant_domains_pickle}")



if __name__ == '__main__':
    parser = ap.ArgumentParser(description='Analyze domain data from the pickle file.')
    parser.add_argument('--pickle_file', required=True, help='Path to the pickle file containing scores_dict')
    parser.add_argument('--output_dir', required=True, help='Directory to save analysis outputs')
    parser.add_argument('--threads', type=int, default=62, help='Number of threads for parallel processing (default: 12)')
    parser.add_argument('--assembly', type=str, default='hg19', help='Assembly version, either hg19 or hg38')
    args = parser.parse_args()
    
    analyze_domain_data(args.pickle_file, 
                        args.output_dir, 
                        args.threads,
                        assembly=args.assembly)
