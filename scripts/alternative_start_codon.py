#!/usr/bin/env python
"""
Alternative Start Codon Detection for PVS1 Evaluation

Per ClinGen SVI guidelines (PMC6185798), for start_lost variants, we should check
if there's a downstream in-frame ATG codon that could serve as an alternative
translation initiation site.

This module provides functions to:
1. Load and index Ensembl CDS FASTA files
2. Find downstream in-frame ATG codons for any transcript
3. Adjust PVS1 strength based on alternative start codon presence

Download URLs (always latest Ensembl release via "current" symlinks):
- GRCh38: https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
- GRCh37: https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.cds.all.fa.gz
"""

import os
import gzip
import pickle
import logging
from typing import Dict, Tuple, Optional, List
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


class AlternativeStartCodonDetector:
    """
    Detect downstream in-frame ATG codons for start_lost variants.
    Per ClinGen SVI guidelines (PMC6185798).
    """
    
    def __init__(self, cds_fasta_path: str, cache_pkl_path: str = None):
        """
        Initialize with Ensembl CDS FASTA file.
        
        Args:
            cds_fasta_path: Path to Homo_sapiens.GRCh38.cds.all.fa.gz or similar
            cache_pkl_path: Optional path to cache the parsed CDS sequences
        """
        self.cds_fasta_path = cds_fasta_path
        self.cache_pkl_path = cache_pkl_path or cds_fasta_path.replace('.fa.gz', '.cds_cache.pkl.gz')
        self.cds_sequences = {}
        self.transcript_to_key = {}
        
        # Load from cache if available and newer than source
        if self._load_from_cache():
            logger.info(f"Loaded CDS sequences from cache: {self.cache_pkl_path}")
        else:
            self._load_fasta()
            self._save_to_cache()
    
    def _load_from_cache(self) -> bool:
        """Try to load from cache pickle file."""
        if not os.path.exists(self.cache_pkl_path):
            return False
        
        # Check if cache is newer than source FASTA
        if os.path.exists(self.cds_fasta_path):
            if os.path.getmtime(self.cache_pkl_path) < os.path.getmtime(self.cds_fasta_path):
                logger.info("Cache is older than source FASTA, reloading...")
                return False
        
        try:
            with gzip.open(self.cache_pkl_path, 'rb') as f:
                cached_data = pickle.load(f)
                self.cds_sequences = cached_data['cds_sequences']
                self.transcript_to_key = cached_data['transcript_to_key']
            return True
        except Exception as e:
            logger.warning(f"Failed to load cache: {e}")
            return False
    
    def _save_to_cache(self):
        """Save parsed CDS sequences to cache."""
        try:
            with gzip.open(self.cache_pkl_path, 'wb') as f:
                pickle.dump({
                    'cds_sequences': self.cds_sequences,
                    'transcript_to_key': self.transcript_to_key
                }, f)
            logger.info(f"Saved CDS cache to: {self.cache_pkl_path}")
        except Exception as e:
            logger.warning(f"Failed to save cache: {e}")
    
    def _load_fasta(self):
        """Load and parse Ensembl CDS FASTA file."""
        logger.info(f"Loading CDS FASTA from: {self.cds_fasta_path}")
        
        open_func = gzip.open if self.cds_fasta_path.endswith('.gz') else open
        
        current_header = None
        current_seq = []
        
        with open_func(self.cds_fasta_path, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # Save previous sequence
                    if current_header and current_seq:
                        self._process_fasta_entry(current_header, ''.join(current_seq))
                    
                    # Parse new header
                    # Ensembl CDS FASTA header format:
                    # >ENST00000456328.2 cds chromosome:GRCh38:1:11869:14409:1 gene:ENSG00000223972.5 ...
                    current_header = line[1:]  # Remove '>'
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Don't forget the last sequence
            if current_header and current_seq:
                self._process_fasta_entry(current_header, ''.join(current_seq))
        
        logger.info(f"Loaded {len(self.cds_sequences)} CDS sequences")
    
    def _process_fasta_entry(self, header: str, sequence: str):
        """Process a single FASTA entry."""
        # Extract transcript ID from header
        # Format: ENST00000456328.2 cds chromosome:GRCh38:...
        parts = header.split()
        transcript_with_version = parts[0]  # e.g., ENST00000456328.2
        transcript_id = transcript_with_version.split('.')[0]  # e.g., ENST00000456328
        
        # Store the sequence
        self.cds_sequences[transcript_with_version] = sequence.upper()
        
        # Build mapping from transcript ID (without version) to full key
        if transcript_id not in self.transcript_to_key:
            self.transcript_to_key[transcript_id] = transcript_with_version
    
    def find_alternative_start(self, transcript_id: str) -> Dict:
        """
        Find downstream in-frame ATG codons in the CDS.
        
        The CDS FASTA from Ensembl starts with the canonical start codon (ATG).
        We search for additional in-frame ATG codons downstream.
        
        Args:
            transcript_id: Ensembl transcript ID (e.g., ENST00000123456)
            
        Returns:
            Dict with:
                - has_alternative: bool - whether there's a downstream ATG
                - first_alt_position: int - amino acid position of first alt ATG (1-based)
                - all_alt_positions: List[int] - all downstream ATG positions
                - n_terminal_truncation_fraction: float - fraction of protein lost (0-1)
                - total_protein_length: int - total protein length in amino acids
                - error: str or None - error message if lookup failed
        """
        # Remove version if present
        transcript_id_clean = transcript_id.split('.')[0]
        
        # Try to find the CDS sequence
        cds_seq = None
        
        # Try exact match first (with version)
        if transcript_id in self.cds_sequences:
            cds_seq = self.cds_sequences[transcript_id]
        # Try without version
        elif transcript_id_clean in self.transcript_to_key:
            full_key = self.transcript_to_key[transcript_id_clean]
            cds_seq = self.cds_sequences.get(full_key)
        
        if cds_seq is None:
            return {
                'has_alternative': False,
                'first_alt_position': None,
                'all_alt_positions': [],
                'n_terminal_truncation_fraction': None,
                'total_protein_length': None,
                'error': f'Transcript {transcript_id} not found in CDS FASTA'
            }
        
        # Calculate protein length (excluding stop codon)
        # CDS includes start codon but may or may not include stop codon
        total_protein_length = len(cds_seq) // 3
        
        # Verify first codon is ATG (canonical start)
        first_codon = cds_seq[0:3]
        if first_codon != 'ATG':
            logger.warning(f"Transcript {transcript_id} CDS does not start with ATG: {first_codon}")
        
        # Find all in-frame ATG positions after the first codon
        alt_atg_positions = []
        for i in range(3, len(cds_seq) - 2, 3):  # Skip first ATG, in-frame only
            codon = cds_seq[i:i+3]
            if codon == 'ATG':
                aa_position = (i // 3) + 1  # 1-based amino acid position
                alt_atg_positions.append(aa_position)
        
        if alt_atg_positions:
            first_alt = alt_atg_positions[0]
            # Truncation fraction = how much of the N-terminus would be lost
            truncation_fraction = (first_alt - 1) / total_protein_length
            return {
                'has_alternative': True,
                'first_alt_position': first_alt,
                'all_alt_positions': alt_atg_positions,
                'n_terminal_truncation_fraction': round(truncation_fraction, 4),
                'total_protein_length': total_protein_length,
                'error': None
            }
        
        return {
            'has_alternative': False,
            'first_alt_position': None,
            'all_alt_positions': [],
            'n_terminal_truncation_fraction': None,
            'total_protein_length': total_protein_length,
            'error': None
        }
    
    def batch_find_alternative_starts(self, transcript_ids: List[str]) -> Dict[str, Dict]:
        """
        Find alternative start codons for multiple transcripts.
        
        Args:
            transcript_ids: List of Ensembl transcript IDs
            
        Returns:
            Dict mapping transcript_id -> alternative start info
        """
        results = {}
        for tid in transcript_ids:
            results[tid] = self.find_alternative_start(tid)
        return results


def adjust_pvs1_for_start_lost(
    alternative_start_info: Dict,
    spans_critical_domain: bool = False,
    gene_is_lof_intolerant: bool = True
) -> int:
    """
    Determine PVS1 strength for start_lost variants based on alternative start codon.
    Per ClinGen SVI recommendations (PMC6185798).
    
    CRITICAL: PVS1 is ONLY applicable to genes that are LoF intolerant!
    If the gene is LoF tolerant, PVS1 = 0 (Not applicable) regardless of other factors.
    
    For LoF INTOLERANT genes, the decision tree treats start_lost N-terminal deletion
    (from original ATG to alternative ATG) similarly to an in-frame deletion:
    
    1. NO alternative downstream in-frame ATG:
       - Gene is LoF intolerant → PVS1 (Very Strong)
       
    2. Alternative ATG EXISTS - evaluate like in-frame deletion:
       - >10% N-terminal loss + spans critical/intolerant domain → PVS1_Strong
       - >10% N-terminal loss + no critical domain → PVS1_Moderate  
       - ≤10% N-terminal loss + spans critical domain → PVS1_Moderate
       - ≤10% N-terminal loss + no critical domain → PVS1_Supporting
    
    Note: "spans_critical_domain" can be determined by:
    - Known pathogenic variants exist between original and alternative start
    - The lost N-terminal region contains intolerant domains (AM-based)
    - The lost region contains critical functional domains
    
    Args:
        alternative_start_info: Output from find_alternative_start()
        spans_critical_domain: Whether N-terminal region contains critical functional domains
        gene_is_lof_intolerant: Whether the gene is intolerant to LoF (LOEUF < 0.35 or high AM)
        
    Returns:
        PVS1 strength level:
            0 = Not applicable (includes LoF tolerant genes)
            1 = PVS1_Supporting
            2 = PVS1_Moderate
            3 = PVS1_Strong
            4 = PVS1 (Very Strong)
    """
    # CRITICAL: If gene is LoF tolerant, PVS1 does NOT apply
    if not gene_is_lof_intolerant:
        return 0  # Not applicable - gene tolerates LoF
    
    # If we couldn't find the transcript, be conservative
    if alternative_start_info.get('error'):
        return 2  # PVS1_Moderate (default for start_lost in LoF intolerant gene)
    
    has_alt = alternative_start_info.get('has_alternative', False)
    truncation = alternative_start_info.get('n_terminal_truncation_fraction', 0)
    
    if not has_alt:
        # No alternative start codon found - complete loss of canonical translation
        # This is essentially a null variant for this transcript in LoF intolerant gene
        return 4  # PVS1 (Very Strong)
    
    # Alternative start exists - evaluate like in-frame deletion
    if truncation is None:
        truncation = 0
    
    # Apply in-frame deletion logic per ClinGen SVI
    if truncation > 0.10:  # >10% of protein lost (significant deletion)
        if spans_critical_domain:
            return 3  # PVS1_Strong (large deletion spanning critical domain)
        else:
            return 2  # PVS1_Moderate (large deletion, function uncertain)
    else:  # ≤10% of protein lost (small deletion)
        if spans_critical_domain:
            return 2  # PVS1_Moderate (small deletion but in critical region)
        else:
            return 1  # PVS1_Supporting (small deletion, likely tolerated)


def identify_transcripts_with_alternative_starts(
    df: pd.DataFrame,
    detector: AlternativeStartCodonDetector
) -> pd.DataFrame:
    """
    Add alternative start codon information to a DataFrame of start_lost variants.
    
    Args:
        df: DataFrame with 'Feature' column containing transcript IDs
        detector: AlternativeStartCodonDetector instance
        
    Returns:
        DataFrame with added columns:
            - alt_start_exists: bool
            - alt_start_position: int or None
            - alt_start_truncation_fraction: float or None
    """
    # Get unique transcript IDs
    unique_transcripts = df['Feature'].dropna().unique()
    
    # Batch lookup
    alt_start_info = detector.batch_find_alternative_starts(unique_transcripts)
    
    # Map back to DataFrame
    df['alt_start_exists'] = df['Feature'].map(
        lambda t: alt_start_info.get(t, {}).get('has_alternative', False)
    )
    df['alt_start_position'] = df['Feature'].map(
        lambda t: alt_start_info.get(t, {}).get('first_alt_position', None)
    )
    df['alt_start_truncation_fraction'] = df['Feature'].map(
        lambda t: alt_start_info.get(t, {}).get('n_terminal_truncation_fraction', None)
    )
    df['alt_start_protein_length'] = df['Feature'].map(
        lambda t: alt_start_info.get(t, {}).get('total_protein_length', None)
    )
    
    return df


# Singleton instance for reuse
_detector_instance = None


def get_detector(cds_fasta_path: str) -> AlternativeStartCodonDetector:
    """Get or create a singleton detector instance."""
    global _detector_instance
    if _detector_instance is None or _detector_instance.cds_fasta_path != cds_fasta_path:
        _detector_instance = AlternativeStartCodonDetector(cds_fasta_path)
    return _detector_instance


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Detect alternative start codons in transcripts")
    parser.add_argument("--cds_fasta", required=True, help="Path to Ensembl CDS FASTA file")
    parser.add_argument("--transcript", help="Single transcript ID to check")
    parser.add_argument("--transcript_list", help="File with list of transcript IDs (one per line)")
    parser.add_argument("--output", help="Output file for results")
    
    args = parser.parse_args()
    
    # Initialize detector
    detector = AlternativeStartCodonDetector(args.cds_fasta)
    
    # Process transcripts
    if args.transcript:
        result = detector.find_alternative_start(args.transcript)
        print(f"\nTranscript: {args.transcript}")
        for key, value in result.items():
            print(f"  {key}: {value}")
    
    if args.transcript_list:
        with open(args.transcript_list) as f:
            transcripts = [line.strip() for line in f if line.strip()]
        
        results = detector.batch_find_alternative_starts(transcripts)
        
        # Output as TSV
        output_lines = ["transcript_id\thas_alternative\tfirst_alt_position\tn_terminal_truncation_fraction\ttotal_protein_length"]
        for tid, info in results.items():
            output_lines.append(f"{tid}\t{info['has_alternative']}\t{info['first_alt_position']}\t{info['n_terminal_truncation_fraction']}\t{info['total_protein_length']}")
        
        output_text = '\n'.join(output_lines)
        
        if args.output:
            with open(args.output, 'w') as f:
                f.write(output_text)
            print(f"Results written to: {args.output}")
        else:
            print(output_text)

