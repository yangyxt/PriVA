#!/usr/bin/env python3

import argparse
import pysam
import pandas as pd
import subprocess
import os
import logging
from typing import List, Dict, Any
import uuid

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(module)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def generate_variant_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Generate a unique variant ID by concatenating position information"""
    return f"{chrom}_{pos}_{ref}_{alt}"


def parse_splice_site_loc(site_loc: str) -> tuple:
    """Parse splice site location into chrom and pos"""
    chrom, pos = site_loc.split(':')
    return chrom, int(pos)


def create_simple_vcf_header() -> pysam.VariantHeader:
    """Create a minimal VCF header for liftover with hg38 contigs"""
    header = pysam.VariantHeader()
    header.add_meta('fileformat', 'VCFv4.2')

    # Add required hg38 contigs with chr prefix
    for i in range(1, 23):
        header.add_meta('contig', items=[('ID', f'chr{i}')])

    # Add X, Y and M chromosomes
    for chrom in ['X', 'Y', 'M']:
        header.add_meta('contig', items=[('ID', f'chr{chrom}')])

    return header


def liftover_coordinates(input_vcf: str, output_vcf: str, input_bed: str,
                        output_bed: str, chain_file: str, reference_fasta: str):
    """Perform liftover for both VCF and BED files"""
    # Liftover VCF
    cmd_vcf = f"CrossMap vcf --chromid l {chain_file} {input_vcf} {reference_fasta} {output_vcf}"
    subprocess.run(cmd_vcf, shell=True, check=True)

    # Liftover BED
    unmap_bed = output_bed.replace('.bed', '.unmap.bed')
    cmd_bed = f"CrossMap bed --chromid l --unmap-file {unmap_bed} {chain_file} {input_bed} {output_bed}"
    subprocess.run(cmd_bed, shell=True, check=True)


def process_and_update_chunk(chunk: pd.DataFrame, chain_file: str, reference_fasta: str,
                           temp_dir: str) -> pd.DataFrame:
    """Process a single chunk through the entire pipeline"""

    # Generate temporary file paths for this chunk
    chunk_id = str(uuid.uuid4())
    temp_vcf = os.path.join(temp_dir, f"temp_grch38_{chunk_id}.vcf")
    lifted_vcf = os.path.join(temp_dir, f"temp_hg19_{chunk_id}.vcf")
    temp_bed = os.path.join(temp_dir, f"temp_grch38_{chunk_id}.bed")
    lifted_bed = os.path.join(temp_dir, f"temp_hg19_{chunk_id}.bed")

    try:
        # Generate variant IDs for this chunk
        chunk['variant_id'] = chunk.apply(
            lambda row: generate_variant_id(row['#chrom'], row['pos'], row['ref'], row['alt']),
            axis=1
        )

        # Create VCF for this chunk
        with pysam.VariantFile(temp_vcf, 'w', header=create_simple_vcf_header()) as vcf_out:
            for _, row in chunk.iterrows():
                record = vcf_out.new_record()
                # logger.info(f"The record is {row.to_string()}, the chrom is {row['#chrom']}")
                record.chrom = row['#chrom']
                record.pos = row['pos']
                record.id = row['variant_id']
                record.ref = row['ref']
                record.alts = (row['alt'],)
                vcf_out.write(record)

        # Create BED for splice sites
        with open(temp_bed, 'w') as bed_out:
            for _, row in chunk.iterrows():
                chrom, pos = parse_splice_site_loc(row['splicevault_site_loc'])
                bed_out.write(f"{chrom}\t{pos-1}\t{pos}\t{row['variant_id']}\n")

        # Perform liftover for this chunk
        liftover_coordinates(temp_vcf, lifted_vcf, temp_bed, lifted_bed,
                           chain_file, reference_fasta)

        # Read lifted coordinates
        vcf_coords = {}
        with pysam.VariantFile(lifted_vcf) as vcf:
            for record in vcf:
                vcf_coords[record.id] = {
                    'chrom': record.chrom,
                    'pos': record.pos,
                    'ref': record.ref,
                    'alt': record.alts[0]
                }

        bed_coords = {}
        with open(lifted_bed) as f:
            for line in f:
                chrom, start, end, var_id = line.strip().split('\t')
                bed_coords[var_id] = f"{chrom}:{end}"

        # Update coordinates in the chunk
        updated_chunk = chunk.copy()
        unmapped_ids = set()
        for idx, row in updated_chunk.iterrows():
            var_id = row['variant_id']
            matched = False
            if var_id in vcf_coords:
                coords = vcf_coords[var_id]
                updated_chunk.at[idx, '#chrom'] = coords['chrom']
                updated_chunk.at[idx, 'pos'] = coords['pos']
                updated_chunk.at[idx, 'ref'] = coords['ref']
                updated_chunk.at[idx, 'alt'] = coords['alt']
                matched = True
            if var_id in bed_coords:
                updated_chunk.at[idx, 'splicevault_site_loc'] = bed_coords[var_id]
            if not matched:
                unmapped_ids.add(var_id)

        # Remove variant_id column
        updated_chunk = updated_chunk.loc[~updated_chunk['variant_id'].isin(unmapped_ids), :]
        updated_chunk = updated_chunk.drop('variant_id', axis=1)
        return updated_chunk

    finally:
        # Clean up temporary files
        for temp_file in [temp_vcf, lifted_vcf, temp_bed, lifted_bed]:
            if os.path.exists(temp_file):
                os.remove(temp_file)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_tsv', required=True)
    parser.add_argument('--output_tsv', required=True)
    parser.add_argument('--chain_file', required=True)
    parser.add_argument('--reference_fasta', required=True)
    parser.add_argument('--temp_dir', default='./temp')
    parser.add_argument('--chunk_size', type=int, default=40000)
    args = parser.parse_args()

    os.makedirs(args.temp_dir, exist_ok=True)

    # Process header first
    header = None
    with open(args.input_tsv) as f:
        header = f.readline().strip()

    # Write header to output file
    with open(args.output_tsv, 'w') as f:
        f.write(header + '\n')

    # Process chunks
    chunk_iterator = pd.read_csv(args.input_tsv, sep='\t', chunksize=args.chunk_size)
    for i, chunk in enumerate(chunk_iterator):
        logger.info(f"Processing chunk {i+1}...")

        # Process chunk and get updated records
        updated_chunk = process_and_update_chunk(chunk, args.chain_file,
                                               args.reference_fasta, args.temp_dir)

        # Append updated chunk to output file
        if i == 0:
            updated_chunk.to_csv(args.output_tsv, sep='\t', index=False,
                           header=True)
        else:
            updated_chunk.to_csv(args.output_tsv, sep='\t', index=False,
                           header=False, mode='a')

        logger.info(f"Completed chunk {i+1}")

if __name__ == "__main__":
    main()