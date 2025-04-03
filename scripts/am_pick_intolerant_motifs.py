#!/usr/bin/env python3

import os
# Limit OpenBLAS threads
os.environ["OPENBLAS_NUM_THREADS"] = "4"
os.environ["MKL_NUM_THREADS"] = "4"
os.environ["NUMEXPR_NUM_THREADS"] = "4"
os.environ["OMP_NUM_THREADS"] = "4"


import pysam
import numpy as np
from collections import defaultdict
from scipy.stats import gaussian_kde
import logging
import pickle
import argparse
import json
import multiprocessing as mp
from functools import partial
from sklearn.neighbors import KernelDensity



logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)



class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)


class AMMotifAnalyzer:
    """
    A revised class that:

      1) Parses a VCF (via pysam), collecting (position, reference_aa, am_score) for each gene.
      2) Lets you build a 1D KDE in protein-position space (bandwidth is in amino acid units).
      3) Lets you compute a 'density threshold' by choosing only variants with AM score >= 0.564
         (if you wish), combining them across genes, and using that as a global threshold.
      4) Identifies 'intolerant regions' as consecutive positions where density > your chosen threshold.

    Note: Here, "bandwidth=10" is the standard deviation of the Gaussian kernel in
    the coordinate space of amino acid positions. Roughly, each variant influences
    positions within ± (2–3)*bandwidth around it, depending on the Gaussian's shape.
    It does NOT mean we only look at ±10 aa in a strict "box"; it's a Gaussian decay
    with σ = 10.
    """

    def __init__(self, vcf_path: str, bandwidth: float = 2.0):
        """
        Args:
            vcf_path: Path to AlphaMissense VCF file.
        """
        self.vcf_path = vcf_path
        self.bandwidth = bandwidth
        # Store list of (pos, ref_aa, score) for each gene
        self.transcript_variants = defaultdict(list)

    def parse_vcf(self, n_processes=None):
        """
        Parse VCF file using pysam and collect (position, ref_aa, score) by gene.
        The VCF should be annotated by VEP with CSQ field.
        
        Parallel processing by chromosome.
        """
        logger.info(f"Parsing VCF file: {self.vcf_path}")
        if n_processes is None:
            n_processes = mp.cpu_count()
            
        try:
            # First pass to identify chromosomes and extract field indices
            vcf = pysam.VariantFile(self.vcf_path)
            # Extract the CSQ header line
            csq_fields = vcf.header.info['CSQ'].description.split('Format: ')[1].split('|')
            gene_field = csq_fields.index('Gene')
            transcript_field = csq_fields.index('Feature')
            consq_field = csq_fields.index('Consequence')
            
            # Group by chromosome
            chromosomes = list(vcf.header.contigs)
            logger.info(f"Found {len(chromosomes)} chromosomes, processing in parallel with {n_processes} processes")
            
            # Close the first pass file
            vcf.close()
            
            # Process each chromosome in parallel
            with mp.Pool(n_processes) as pool:
                chrom_results = pool.map(
                    partial(self._parse_chromosome, 
                            gene_field=gene_field, 
                            transcript_field=transcript_field,
                            consq_field=consq_field),
                    chromosomes
                )
            
            # Combine results from all chromosomes
            total_variants = 0
            for chrom_transcript_variants in chrom_results:
                for transcript, variants in chrom_transcript_variants:
                    self.transcript_variants[transcript].extend(variants)
                    total_variants += len(variants)
                    
            logger.info(f"Stored {total_variants} variant entries across {len(self.transcript_variants)} transcripts.")
            
        except Exception as e:
            logger.error(f"Error reading VCF file: {e}")
            raise
    
    def _parse_chromosome(self, chromosome, gene_field, transcript_field, consq_field):
        """
        Parse VCF records for a specific chromosome.
        
        Args:
            chromosome: Chromosome name
            gene_field: Index of gene field in CSQ
            transcript_field: Index of transcript field in CSQ
            consq_field: Index of consequence field in CSQ
            
        Returns:
            List of tuples [(transcript_id, [(pos, ref_aa, score), ...]), ...]
        """
        chrom_variants = defaultdict(list)
        try:
            vcf = pysam.VariantFile(self.vcf_path)
            count = 0
            
            # Create a fetch iterator for just this chromosome
            for record in vcf.fetch(chromosome):
                try:
                    am_score = float(record.info['AM_PATHOGENICITY'])
                    target_tranx = record.info.get('TRANSCRIPT').split('.')[0]
                    pvar = record.info.get('PVAR')
                    transcript_specific_csqs = record.info.get("CSQ", [])
                    
                    if not transcript_specific_csqs:
                        continue

                    for tranx_specific_anno in transcript_specific_csqs:
                        tranx_specific_anno = tranx_specific_anno.split('|')
                        gene = tranx_specific_anno[gene_field]
                        transcript = tranx_specific_anno[transcript_field]
                        consequence = tranx_specific_anno[consq_field]

                        if consequence != 'missense_variant':
                            continue
                        
                        if transcript != target_tranx:
                            continue
                        
                        if not gene or not pvar:
                            continue

                        # Example PVAR: "M10I" => ref_aa = 'M', aa_pos = 10, alt_aa = 'I'
                        regex_groups = re.match(r'^([A-Z]+)([0-9]+)([A-Z]+)$', pvar)
                        ref_aa = regex_groups.group(1)  # 'M'
                        aa_pos_str = regex_groups.group(2)  # '10'
                        aa_pos = int(aa_pos_str)
                        
                        chrom_variants[target_tranx].append((aa_pos, ref_aa, am_score))
                        count += 1
                except (KeyError, ValueError) as e:
                    continue
                    
            logger.info(f"Processed {count} variant entries on chromosome {chromosome}")
            vcf.close()
            
            # Convert defaultdict to list of tuples for pickling
            return list(chrom_variants.items())
            
        except Exception as e:
            logger.error(f"Error processing chromosome {chromosome}: {e}")
            return []



    def kde_local_density(self, positions: np.ndarray, scores: np.ndarray, query_positions: np.ndarray) -> np.ndarray:
        """
        Compute weighted KDE along protein positions using sklearn's KernelDensity
        
        Args:
            positions: variant positions
            scores: AM scores (weights)
            query_positions: where to evaluate density
        """
        if len(positions) < 2:
            return np.zeros_like(query_positions)
        
        # Reshape for sklearn
        X = positions.reshape(-1, 1)
        query_positions = query_positions.reshape(-1, 1)
        
        # Fit weighted KDE
        kde = KernelDensity(kernel='gaussian', bandwidth=self.bandwidth)
        kde.fit(X, sample_weight=scores)
        
        # Evaluate density at query positions
        log_density = kde.score_samples(query_positions)
        return np.exp(log_density)


    def single_data_influence_calculation(self, distance: float = 1.0) -> float:
        scaled_distance = distance / self.bandwidth
        kernel_height = 1 / np.sqrt(2 * np.pi) * np.exp(- scaled_distance**2 / 2)
        return kernel_height


    def compute_density_threshold(self, 
                                  target_score: float = 0.564,
                                  total_weight: float = np.nan,
                                  gene_name: str = None) -> float:
        """
        Compute the density for a base located at the boundary of a 5 consecutive residues (minimal motif size) with avg score of 0.564, outside the window all residues have score 0.2 (benign AM)
        
        Since the bandwidth is approximately 2 or 3, we can only consider the nearby 4 bases because the farther bases will be too little influence to be considered.
        Consider the left 4 bases (with weight of 0.564), and the right 4 bases (with weight of 0.2), and the query base itself (with weight of 0.564). 
        
        Args:
            bandwidth: Gaussian kernel width (standard deviation)
            target_score: Score threshold (e.g., 0.564)
            protein_length: Length of the protein
            gene_name: Name of the gene
        """
        # Calculate the nearby base contributions
        total_contribution = 0
        for distance in range(0, 5):
            contribution = self.single_data_influence_calculation(distance)
            if distance > 0:
                weighted_contribution = contribution * target_score + 0.2 * contribution
            else:
                weighted_contribution = contribution * target_score
            total_contribution += weighted_contribution

        denominator = total_weight * self.bandwidth
        density_threshold = total_contribution / denominator
        
        logger.info(f"Density threshold computed: {density_threshold:.4g} for {gene_name} with total weight {total_weight}")
        return density_threshold



    def process_single_transcript(self, transcript_data) -> tuple[str, dict]:
        """
        Process a single gene's data. Must be a standalone function for Pool.
        
        Args:
            transcript_data: tuple of (gene_id, variant_list)
            
        Returns:
            tuple of (gene_id, regions_dict)
        """
        transcript, variants = transcript_data
        if len(variants) < 5:
            return transcript, {}
            
        # Group variants by position
        pos_variants = defaultdict(list)
        pos_ref_aa = {}
        for pos, ref_aa, score in variants:
            pos_variants[pos].append(score)
            pos_ref_aa[int(pos)] = ref_aa
        
        # Get arrays - convert to numpy arrays immediately
        positions = np.array(sorted(pos_variants.keys()))
        max_scores = np.array([max(pos_variants[p]) for p in positions])
        min_scores = np.array([min(pos_variants[p]) for p in positions])

        assert np.all(max_scores >= min_scores), f"Error in {transcript}: max_scores < min_scores at some positions: Here is the max_scores: {max_scores}, and min_scores: {min_scores}"
        
        if not len(positions):
            return transcript, {}
            
        # Evaluate KDE at all integer positions
        query_pos = np.arange(min(positions), max(positions) + 1)
        
        # Fit KDE to positions, weighted by scores
        kde_max = self.kde_local_density(positions, max_scores, query_pos)
        kde_min = self.kde_local_density(positions, min_scores, query_pos)

        # Compute expected density for 5aa window with avg score 0.564
        max_density_threshold = self.compute_density_threshold(
            target_score=0.564,
            total_weight=np.sum(max_scores),
            gene_name=transcript
        )
        logger.info(f"density_threshold for max_score: {max_density_threshold}")

        min_density_threshold = self.compute_density_threshold(
            target_score=0.2,
            total_weight=np.sum(min_scores),
            gene_name=transcript
        )
        logger.info(f"density_threshold for min_score: {min_density_threshold}")
        
        # Find regions where local weighted density indicates avg score > 0.564
        max_regions = self._find_high_density_regions(query_pos, kde_max, pos_ref_aa, max_density_threshold)
        logger.info(f"Found {len(max_regions)} max_regions for {transcript}, they looks like {max_regions}")
        min_regions = self._find_high_density_regions(query_pos, kde_min, pos_ref_aa, min_density_threshold)
        logger.info(f"Found {len(min_regions)} min_regions for {transcript}, they looks like {min_regions}")
        
        if max_regions or min_regions:
            return transcript, {
                'max_score_regions': max_regions,
                'min_score_regions': min_regions
            }
        return transcript, {}



    def identify_intolerant_regions(self, n_processes: int = None) -> dict:
        """
        Parallel version of region identification.
        
        Args:
            bandwidth: KDE bandwidth parameter
            n_processes: number of processes to use (default: CPU count)
        """
        if not self.transcript_variants:
            self.parse_vcf()
            
        if n_processes is None:
            n_processes = mp.cpu_count()
            
        # Convert data to list of tuples for Pool
        transcript_data = list(self.transcript_variants.items())
        
        # Create Pool and process genes in parallel
        with mp.Pool(n_processes) as pool:
            results = pool.map(self.process_single_transcript, transcript_data)
            
        # Combine results into final dictionary
        return {tid: regions for tid, regions in results if regions}



    def _find_high_density_regions(self, 
                                   positions: np.ndarray, 
                                   densities: np.ndarray, 
                                   pos_ref_aa: dict,
                                   density_threshold: float) -> list:
        """
        Find consecutive positions where density exceeds threshold.
        Threshold is computed based on what density we'd expect from
        a 5aa window with average score of 0.564.
        
        Args:
            positions: array of positions (e.g., np.arange(min_pos, max_pos + 1))
            densities: array of density values at each position
        
        Returns:
            set of positions, each position is a string of REF_AA + POSITION: "A10"
        """
        
        regions = set()
        current_positions = set()
        
        for pos, density in zip(positions, densities):
            if density > density_threshold:
                if int(pos) in pos_ref_aa:
                    current_positions.add(f"{pos_ref_aa[int(pos)]}{pos}")
            else:
                if current_positions:  # end of a region
                    regions.update(current_positions)
                    current_positions = set()
                
        # Handle last region if exists
        if current_positions:
            regions.update(current_positions)
        
        return regions


def main():
    parser = argparse.ArgumentParser(
        description="Analyze local density of variant positions along proteins."
    )
    parser.add_argument('--vcf_path', required=True, help='Path to AlphaMissense VCF file')
    parser.add_argument('--bandwidth', type=float, default=2.0,
                        help='KDE bandwidth in amino acid units (default=2.0)')
    parser.add_argument('--processes', type=int, default=None,
                        help='Number of processes to use (default: CPU count)')
    parser.add_argument('--output', help='Output file (.pkl or .json) for results')
    args = parser.parse_args()

    analyzer = AMMotifAnalyzer(args.vcf_path, args.bandwidth)
    
    # Parse VCF using same number of processes
    analyzer.parse_vcf(n_processes=args.processes)

    results = analyzer.identify_intolerant_regions(n_processes=args.processes)

    # Save or print results
    if args.output:
        if args.output.endswith('.pkl'):
            with open(args.output, 'wb') as f:
                pickle.dump(results, f)
            logger.info(f"Results saved to {args.output} in pickle format.")
        elif args.output.endswith('.json'):
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2, cls=NpEncoder)
            logger.info(f"Results saved to {args.output} in JSON format.")
        else:
            logger.warning("Unknown output file format. No file saved.")
    else:
        logger.info("No output file specified; printing results.")
        logger.info(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
