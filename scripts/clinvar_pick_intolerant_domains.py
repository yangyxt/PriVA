import os
import logging
import argparse as ap
import multiprocessing as mp
from collections import defaultdict
from typing import Dict, List
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import pickle

from stat_protein_domain_clinvars import nested_defaultdict

# Set up logging
logger = logging.getLogger(__name__)
logging.basicConfig(
    format='%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)d:%(message)s',
    level=logging.INFO
)


def filter_out_lowrevstatus(d, high_confidence_status: List[str]):
    if 'distribution' in d:
        # Get indices of high confidence variants
        high_conf_idx = [
            i for i, status in enumerate(d['distribution']['clinrevstatus'])
            if any(hs in status.lower() for hs in high_confidence_status)
        ]
        
        if high_conf_idx:  # Only keep domains with high confidence variants
            new_dist = {
                'clinrevstatus': [d['distribution']['clinrevstatus'][i] for i in high_conf_idx],
                'clnsig': [d['distribution']['clnsig'][i] for i in high_conf_idx],
                'mole_consq': [d['distribution']['mole_consq'][i] for i in high_conf_idx]
            }
            return {'distribution': new_dist}
    return None


def recursive_filter_revstatus(d, high_confidence_status: List[str]):
    if not isinstance(d, dict):
        return d
    
    filtered = nested_defaultdict()
    for k, v in d.items():
        if k == 'distribution':
            processed = filter_out_lowrevstatus(d, high_confidence_status)
            if processed:
                return processed
        else:
            result = recursive_filter_revstatus(v, high_confidence_status)
            if result:
                filtered[k] = result
    
    return filtered if filtered else None



def filter_high_confidence_variants(clinvar_dict: Dict) -> Dict:
    """
    Filter variants to keep only those with review status of 1 star or above.
    
    Star Rating:
    ★★★★ practice_guideline
    ★★★  reviewed_by_expert_panel
    ★★   criteria_provided,_multiple_submitters,_no_conflicts
    ★    criteria_provided,_conflicting_classifications
    ★    criteria_provided,_single_submitter
    """
    high_confidence_status = {
        # Higher confidence (2+ stars)
        'practice_guideline',                                    # 4 stars
        'reviewed_by_expert_panel',                             # 3 stars
        'criteria_provided,_multiple_submitters,_no_conflicts',  # 2 stars
    }
    
    
    return recursive_filter_revstatus(clinvar_dict, high_confidence_status)


def is_not_patho(clnsig: str) -> bool:
    return any(term in clnsig.lower() for term in ['benign', 'likely_not_patho'])


def is_pathogenic(clnsig: str) -> bool:
    return any(term in clnsig.lower() for term in ['pathogenic', 'likely_pathogenic'])


def is_missense(consq: str) -> bool:
    """Check if a variant consequence is missense."""
    return 'missense_variant' in consq.lower()


def is_protein_disrupting(consq: str) -> bool:
    """
    Check if a variant consequence likely disrupts protein function through
    truncation, length change, NMD, or essential splice site alterations.
    
    Args:
        consq: Consequence string from ClinVar
        
    Returns:
        bool: True if variant is likely protein-disrupting
    """
    protein_disrupting_terms = {
        # Truncating variants
        'stop_gained',               # Premature stop codon
        'frameshift_variant',        # Reading frame disruption
        'start_lost',               # Loss of start codon
        'stop_lost',                # Loss of stop codon
        
        # Essential splice site variants
        'splice_acceptor_variant',   # Changes 3' splice site
        'splice_donor_variant',      # Changes 5' splice site
        
        # NMD-triggering variants
        'NMD_transcript_variant',    # Nonsense-mediated decay target
        
        # Protein length changes
        'feature_truncation',        # Reduction of genomic feature
        'inframe_deletion',          # Deletion without frameshift
        'inframe_insertion',         # Insertion without frameshift
        'feature_elongation',		# Increase in genomic feature size
        
        # Complete loss
        'transcript_ablation'        # Complete removal of transcript
    }
    
    return any(term in consq.lower() for term in protein_disrupting_terms)


def vartype_corr_exam(domain_dist: Dict) -> Dict:
    """
    Create a 2x2 contingency table for a domain's variants.
    
    Args:
        domain_dist: Dictionary containing 'clnsig' and 'mole_consq' lists
        
    Returns:
        Dictionary containing counts and statistics.
        {
            'table': {
                'missense_patho': int,
                'missense_not_patho': int,
                'truncating_patho': int,
                'truncating_not_patho': int
            },
            'pvalue': float,
            'odds_ratio': float
        }
        Fisher p-value will be calculated if:
        - Total sample size >= 5
        - At least one variant in each variant type (missense/truncating)
        - At least one variant in each clinical class (patho/benign)
        
        Chi-square test will be used if:
        - Total sample size >= 20
        - All cell counts >= 5
        
        Fisher's exact test will be used if:
        - Above conditions for p-value are met but chi-square conditions aren't
    """
    # Initialize counts
    counts = {
        'missense_patho': 0,
        'missense_not_patho': 0,
        'truncating_patho': 0,
        'truncating_not_patho': 0
    }
    
    # Count variants
    for clnsig, consq in zip(domain_dist['clnsig'], domain_dist['mole_consq']):
        is_miss = is_missense(consq)
        is_truncating = is_protein_disrupting(consq)
        is_path = is_pathogenic(clnsig)

        is_miss = False if is_truncating else is_miss
        
        if is_miss and is_path:
            counts['missense_patho'] += 1
        elif is_miss and not is_path:
            counts['missense_not_patho'] += 1
        elif is_truncating and is_path:
            counts['truncating_patho'] += 1
        elif is_truncating and not is_path:
            counts['truncating_not_patho'] += 1
    
    # Create contingency table
    table = [[counts['missense_patho'], counts['missense_not_patho']],
             [counts['truncating_patho'], counts['truncating_not_patho']]]
    
    # Check conditions for statistical testing
    total_count = sum(counts.values())
    missense_total = counts['missense_patho'] + counts['missense_not_patho']
    truncating_total = counts['truncating_patho'] + counts['truncating_not_patho']
    patho_total = counts['missense_patho'] + counts['truncating_patho']
    not_patho_total = counts['missense_not_patho'] + counts['truncating_not_patho']
    min_cell_count = min(counts.values())
    
    # Set default values
    p_value = np.nan
    odds_ratio = np.nan
    test_used = 'none'
    
    # Check if any statistical test can be performed
    if (total_count >= 5 and 
        missense_total > 0 and truncating_total > 0 and 
        patho_total > 0 and not_patho_total > 0):
        
        # Determine which test to use
        if total_count >= 20 and min_cell_count >= 5:
            # Use Chi-square test
            chi2, p_value = stats.chi2_contingency(table)[0:2]
            odds_ratio = ((counts['missense_patho'] * counts['truncating_not_patho']) / 
                         (counts['missense_not_patho'] * counts['truncating_patho']))
            test_used = 'chi_square'
        else:
            # Use Fisher's exact test
            odds_ratio, p_value = stats.fisher_exact(table)
            test_used = 'fisher'
            
        # Adjust p-value of 1 to 0.99999 to avoid log10(1) = 0
        p_value = min(p_value, 0.99999)

    
    return {
        'table': counts,
        'pvalue': p_value,
        'odds_ratio': odds_ratio,
        'test_used': test_used,
        'missense_total': missense_total,
        'truncating_total': truncating_total,
        'patho_total': patho_total,
        'not_patho_total': not_patho_total,
        'total_variants': total_count
    }


def plot_clinsig_distribution(domain_path: str, clinsig_list: List[str], output_dir: str, is_intolerant: bool = False):
    """
    Create histogram of clinical significance categories for a domain.
    
    Args:
        domain_path: Name/path of the domain
        clinsig_list: List of clinical significance values
        output_dir: Base directory for output
        is_intolerant: Whether this domain is classified as intolerant
    """
    # Define the categories we care about
    categories = [
        'Pathogenic',
        'Likely_pathogenic',
        'Uncertain_significance',
        'Likely_not_patho',
        'Benign'
    ]
    
    # Count occurrences
    category_counts = defaultdict(int)
    for clinsig in clinsig_list:
        for cat in categories:
            if cat.lower() in clinsig.lower():
                category_counts[cat] += 1
    
    # Create histogram
    plt.figure(figsize=(10, 6))
    counts = [category_counts[cat] for cat in categories]
    bars = plt.bar(categories, counts)
    
    # Customize plot
    plt.title(f"Clinical Significance Distribution\n{domain_path}")
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('Count')
    
    # Add count labels on top of bars
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}',
                ha='center', va='bottom')
    
    # Determine output subdirectory
    subdir = 'intolerant_domains' if is_intolerant else 'all_domains'
    plot_dir = os.path.join(output_dir, 'clinsig_plots', subdir)
    os.makedirs(plot_dir, exist_ok=True)
    
    # Save plot
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, f"clinsig_dist_{domain_path}.png"))
    plt.close()


def process_single_domain(args):
    """
    Process a single domain's distribution data.
    
    Args:
        args: Tuple containing (domain_path, distribution, output_dir, benign_threshold, min_variants)
        
    Returns:
        DataFrame row with results or None if processing fails
    """
    domain_path, dist, output_dir, benign_threshold, min_variants, pathogenic_threshold = args
    
    try:
        clnsig_list = dist['clnsig']
        if len(clnsig_list) < min_variants:
            return None
            
        benign_count = sum(1 for sig in clnsig_list if is_not_patho(sig))
        pathogenic_count = sum(1 for sig in clnsig_list if is_pathogenic(sig))
        benign_fraction = benign_count / len(clnsig_list)
        pathogenic_fraction = pathogenic_count / len(clnsig_list)
        
        # Create distribution plot
        plot_clinsig_distribution(domain_path, clnsig_list, output_dir, 
                                is_intolerant=(benign_fraction <= benign_threshold and pathogenic_fraction >= pathogenic_threshold))
        
        benign_stat = {
                'domain': domain_path,
                'benign_fraction': benign_fraction,
                'pathogenic_fraction': pathogenic_fraction,
                'total_variants': len(clnsig_list),
                'benign_variants': benign_count
            }
        
        if benign_fraction <= benign_threshold and pathogenic_fraction >= pathogenic_threshold:
            logger.info(f"Domain {domain_path} has been identified as intolerant with {benign_count} benign variants and {pathogenic_count} pathogenic variants out of {len(clnsig_list)} total variants")
            return benign_stat
        
    except Exception as e:
        logger.error(f"Error processing domain {domain_path}: {str(e)}")
        return None


def collect_domain_data(d, output_dir, benign_threshold, min_variants, pathogenic_threshold, path=[]):
    """Collect domain data for parallel processing"""
    domain_data = []
    for k, v in d.items():
        if k == 'distribution':
            domain_path = ':'.join(path)
            domain_data.append((domain_path, v, output_dir, benign_threshold, min_variants, pathogenic_threshold))
        elif isinstance(v, dict):
            domain_data.extend(collect_domain_data(v, output_dir, benign_threshold, min_variants, pathogenic_threshold, path + [k]))
    return domain_data


def identify_intolerant_domains(clinvar_dict: Dict, 
                              output_dir: str,
                              benign_threshold: float = 0.1,
                              pathogenic_threshold: float = 0.5,
                              min_variants: int = 5,
                              threads: int = 1) -> Dict[str, float]:
    """
    Identify domains that have very few benign variants and create distribution plots.
    
    Args:
        clinvar_dict: Filtered dictionary containing only high confidence variants
        output_dir: Directory to save outputs
        benign_threshold: Maximum fraction of benign variants allowed
        min_variants: Minimum number of variants required in a domain
        threads: Number of CPU threads to use
        
    Returns:
        Dictionary mapping domain paths to their statistics
    """
    # Collect all domain data
    domain_data = collect_domain_data(clinvar_dict, output_dir, benign_threshold, min_variants, pathogenic_threshold)
    
    # Process domains in parallel
    threads = min(threads, len(domain_data), mp.cpu_count()-1)
    logger.info(f"Processing {len(domain_data)} domains using {threads} cores")
    
    with mp.Pool(threads) as pool:
        results = list(pool.imap_unordered(process_single_domain, domain_data))
    
    # Filter out None results and convert to dictionary
    results = [r for r in results if r is not None]
    intolerant_domains = {r['domain']: r for r in results}
    
    return intolerant_domains



def process_domain_mechanism(domain_path: str, dist: Dict) -> Dict:
    """
    Process a single domain's mechanism statistics.
    
    Args:
        domain_path: Path identifying the domain
        dist: Dictionary containing domain distribution data
        
    Returns:
        Dictionary containing mechanism statistics or None if invalid
    """
    if 'distribution' not in dist:
        return None
        
    # Create contingency table
    stats = vartype_corr_exam(dist['distribution'])
    table = stats['table']
    
    # Calculate additional metrics
    total_patho = table['missense_patho'] + table['truncating_patho']
    if total_patho > 0:
        truncating_fraction = table['truncating_patho'] / total_patho
    else:
        truncating_fraction = 0
    
    return {
        'domain': domain_path,
        'missense_pathogenic': table['missense_patho'],
        'missense_not_patho': table['missense_not_patho'],
        'truncating_pathogenic': table['truncating_patho'],
        'truncating_not_patho': table['truncating_not_patho'],
        'total_variants': sum(table.values()),
        'total_pathogenic': total_patho,
        'truncating_fraction': truncating_fraction,
        'pvalue': stats['pvalue'],
        'odds_ratio': stats['odds_ratio']
    }


def traverse_mechanism_dict(d: Dict, results: List, path: List = []):
    """
    Traverse the nested dictionary to collect mechanism data.
    
    Args:
        d: Nested dictionary to traverse
        results: List to collect results
        path: Current path in the dictionary
    """
    for k, v in d.items():
        if k == 'distribution':
            domain_path = ':'.join(path)
            result = process_domain_mechanism(domain_path, d)
            if result:
                results.append(result)
        elif isinstance(v, dict):
            traverse_mechanism_dict(v, results, path + [k])




def analyze_domain_mechanisms(clinvar_dict: Dict, output_dir: str, threads: int, assembly: str = 'hg19') -> pd.DataFrame:
    """
    Analyze variant type distributions for all domains and genes.
    
    Args:
        clinvar_dict: Nested dictionary containing ClinVar data
        output_dir: Directory to save outputs
        threads: Number of CPU threads to use
        
    Returns:
        DataFrame containing mechanism analysis results
    """
    results = []
    
    # Process all domains
    traverse_mechanism_dict(clinvar_dict, results)
    
    # Convert results to DataFrame and add gene-level analysis
    if results:
        df = pd.DataFrame(results)
        
        # Add gene-level statistics
        df = analyze_gene_mechanisms(df, threads)
        
        # Save results
        output_file = os.path.join(output_dir, f'domain_mechanism_analysis.{assembly}.tsv')
        df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Saved mechanism analysis results to {output_file}")
        
        # Create summary plots
        create_mechanism_plots(df, output_dir)
        
        return df
    
    return None


def analyze_gene_mechanisms(df: pd.DataFrame, threads: int) -> pd.DataFrame:
    '''
    Aggregate domain-level data to gene level and perform statistical tests.
    
    Args:
        df: DataFrame containing domain mechanism analysis results
        threads: Number of CPU threads to use
        
    Returns:
        Original DataFrame with added gene-level statistics
    '''
    # Extract gene IDs
    df['gene_id'] = df['domain'].str.split(':').str[0]
    
    # Group by gene and sum the contingency table values
    gene_stats = df.groupby('gene_id').agg({
        'missense_pathogenic': 'sum',
        'missense_not_patho': 'sum',
        'truncating_pathogenic': 'sum',
        'truncating_not_patho': 'sum',
        'total_variants': 'sum',
        'total_pathogenic': 'sum'
    }).reset_index()
    
    # Convert DataFrame rows to dictionaries for parallel processing
    gene_dicts = gene_stats.to_dict('records')
    
    # Calculate gene-level statistics using parallel processing
    threads = min(threads, len(gene_dicts), mp.cpu_count()-1)
    logger.info(f"Calculating {len(gene_dicts)} gene statistics using {threads} cores")
    with mp.Pool(threads) as pool:
        results = pool.map(calculate_gene_stats, gene_dicts)
    
    # Convert results to DataFrame
    gene_level_df = pd.DataFrame(results)
    
    # Merge gene-level statistics back to original domain DataFrame
    df = df.merge(gene_level_df[['gene_id', 'gene_pvalue', 'gene_odds_ratio', 'gene_test_used']], 
                 on='gene_id', 
                 how='left')
    
    return df


def calculate_gene_stats(row_dict: dict) -> dict:
    '''
    Calculate statistical tests for a gene's contingency table.
    
    Args:
        row_dict: Dictionary containing gene statistics
        
    Returns:
        Dictionary containing gene statistics and test results
    '''
    # Create contingency table
    table = [
        [row_dict['missense_pathogenic'], row_dict['missense_not_patho']],
        [row_dict['truncating_pathogenic'], row_dict['truncating_not_patho']]
    ]
    
    # Initialize results
    result_stats = {
        'gene_id': row_dict['gene_id'],
        'gene_pvalue': np.nan,
        'gene_odds_ratio': np.nan,
        'gene_test_used': 'none'
    }
    
    # Check if we have enough data for statistical testing
    total_count = sum(sum(row) for row in table)
    min_cell_count = min(min(row) for row in table)
    
    if (total_count >= 5 and 
        all(sum(row) > 0 for row in table) and
        all(sum(col) > 0 for col in zip(*table))):
        try:
            if total_count >= 20 and min_cell_count >= 5:
                # Use Chi-square test
                chi2, pvalue = stats.chi2_contingency(table)[0:2]
                result_stats['gene_test_used'] = 'chi_square'
            else:
                # Use Fisher's exact test
                odds_ratio, pvalue = stats.fisher_exact(table)
                result_stats['gene_test_used'] = 'fisher'
            
            # Calculate odds ratio
            odds_ratio = ((row_dict['missense_pathogenic'] * row_dict['truncating_not_patho']) / 
                         (row_dict['missense_not_patho'] * row_dict['truncating_pathogenic'])) if row_dict['truncating_pathogenic'] > 0 and row_dict['missense_not_patho'] > 0 else np.inf
            
            result_stats['gene_pvalue'] = min(pvalue, 0.99999)
            result_stats['gene_odds_ratio'] = odds_ratio
            
        except Exception as e:
            logger.warning(f"Error calculating statistics for gene {row_dict['gene_id']}: {str(e)}")
    
    return result_stats


def create_mechanism_plots(df: pd.DataFrame, output_dir: str):
    """Create domain and gene level mechanism plots."""
    # Domain-level plot
    plt.figure(figsize=(10, 6))
    plt.scatter(df['truncating_fraction'], -np.log10(df['pvalue']))
    plt.xlabel('Fraction of Pathogenic Variants that are Truncating')
    plt.ylabel('-log10(Domain P-value)')
    plt.title('Domain Mechanism Analysis')
    plt.savefig(os.path.join(output_dir, 'domain_mechanism_plot.png'))
    plt.close()
    
    # Gene-level plot
    plt.figure(figsize=(10, 6))
    gene_df = df.drop_duplicates('gene_id')
    plt.scatter(gene_df['gene_odds_ratio'], -np.log10(gene_df['gene_pvalue']))
    plt.xlabel('Gene-level Odds Ratio (log scale)')
    plt.xscale('log')
    plt.ylabel('-log10(Gene P-value)')
    plt.title('Gene Mechanism Analysis')
    plt.savefig(os.path.join(output_dir, 'gene_mechanism_plot.png'))
    plt.close()


def main(pickle_file: str, output_dir: str, threads: int, assembly: str = 'hg19'):
    """
    Analyze domain data from the pickle file.
    
    Args:
        pickle_file: Path to the pickle file containing clinvar_dict
        output_dir: Directory to save analysis outputs
        threads: Number of CPU threads to use
        assembly: Assembly version, either 'hg19' or 'hg38'
    """
    # Load the data
    with open(pickle_file, 'rb') as f:
        clinvar_dict = pickle.load(f)
    
    # Filter and analyze ClinVar data
    filtered_clinvar = filter_high_confidence_variants(clinvar_dict)
    intolerant_domains = identify_intolerant_domains(filtered_clinvar, output_dir, threads=threads)
    
    # Save intolerant domains results
    results_df = pd.DataFrame.from_dict(intolerant_domains, orient='index')
    results_df.index.name = 'domain'
    results_df.to_csv(os.path.join(output_dir, f'intolerant_domains.{assembly}.tsv'), 
                     sep='\t', index=False)
    
    logger.info(f"Found {len(intolerant_domains)} intolerant domains")

    # Analyze domain mechanisms
    mechanism_df = analyze_domain_mechanisms(filtered_clinvar, output_dir, threads=threads, assembly=assembly)


if __name__ == "__main__":
    parser = ap.ArgumentParser(description="Analyze ClinVar domain data")
    parser.add_argument("--pickle_file", required=True, help="Path to pickle file containing clinvar_dict")
    parser.add_argument("--output_dir", required=True, help="Directory to save analysis outputs")
    parser.add_argument("--threads", type=int, default=40, help="Number of CPU threads to use (default: CPU count - 2)")
    parser.add_argument("--assembly", type=str, default='hg19', help="Assembly version, either 'hg19' or 'hg38'")

    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Run main function
    main(args.pickle_file, args.output_dir, args.threads, args.assembly)
