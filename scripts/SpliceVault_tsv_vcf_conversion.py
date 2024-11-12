#!/usr/bin/env python3

import argparse
import pysam
import pandas as pd
import subprocess
import os
import logging
from typing import List, Dict, Any


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(module)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)



def create_vcf_header(reference_fasta: str) -> pysam.VariantHeader:
    """Create VCF header with all necessary INFO fields"""
    header = pysam.VariantHeader()

    # Add reference contigs
    with pysam.FastaFile(reference_fasta) as fasta:
        for contig in fasta.references:
            header.add_meta('contig', items=[('ID', contig)])

    # Add INFO fields for all SpliceVault columns
    info_fields = {
        'TRANSCRIPT_ID': {'Number': '1', 'Type': 'String', 'Description': "Ensembl transcript ID"},
        'SITE_TYPE': {'Number': '1', 'Type': 'String', 'Description': "SpliceVault site type (Donor_loss/Acceptor_loss)"},
        'SITE_LOC': {'Number': '1', 'Type': 'String', 'Description': "SpliceVault site location"},
        'SPLICEAI_DELTA': {'Number': '1', 'Type': 'Float', 'Description': "SpliceAI delta score"},
        'OUT_OF_FRAME': {'Number': '1', 'Type': 'String', 'Description': "Out of frame events"},
        'TOP4_EVENTS': {'Number': '1', 'Type': 'String', 'Description': "Top 4 splicing events"},
        'SAMPLE_COUNT': {'Number': '1', 'Type': 'Integer', 'Description': "Sample count"},
        'MAX_DEPTH': {'Number': '1', 'Type': 'Integer', 'Description': "Maximum depth"}
    }

    for field, props in info_fields.items():
        header.add_meta('INFO', items=[
            ('ID', field),
            ('Number', props['Number']),
            ('Type', props['Type']),
            ('Description', props['Description'])
        ])

    return header

def tsv_to_vcf(input_tsv: str, output_vcf: str, reference_fasta: str) -> None:
    """Convert SpliceVault TSV to VCF format"""
    header = create_vcf_header(reference_fasta)

    with pysam.VariantFile(output_vcf, 'w', header=header) as vcf_out:
        # Process TSV file in chunks
        for chunk in pd.read_csv(input_tsv, sep='\t', chunksize=10000):
            for _, row in chunk.iterrows():
                # Create new record
                record = vcf_out.new_record()

                # Set basic variant info
                record.chrom = row['#chrom']
                record.pos = int(row['pos'])
                record.ref = row['ref']
                record.alts = (row['alt'],)
                # Set END to be equal to or greater than POS
                record.stop = record.pos + len(record.ref) - 1  # This sets END properly

                # Add INFO fields
                record.info['TRANSCRIPT_ID'] = row['transcript_id']
                record.info['SITE_TYPE'] = row['splicevault_site_type']
                record.info['SITE_LOC'] = row['splicevault_site_loc']
                record.info['SPLICEAI_DELTA'] = float(row['spliceai_delta'])
                record.info['OUT_OF_FRAME'] = row['splicevault_out_of_frame_events']
                record.info['TOP4_EVENTS'] = row['splicevault_top4_events']
                record.info['SAMPLE_COUNT'] = int(row['splicevault_site_sample_count'])
                record.info['MAX_DEPTH'] = int(row['splicevault_site_max_depth'])

                # logger.debug(f"Writing record: \n{row.to_string()}\nto VCF record: \n{record}\n")
                vcf_out.write(record)

def vcf_to_tsv(input_vcf: str, output_tsv: str) -> None:
    """Convert VCF back to SpliceVault TSV format"""
    # Open input VCF
    vcf = pysam.VariantFile(input_vcf)

    # Prepare header
    header = ['#chrom', 'pos', 'ref', 'alt',
              'transcript_id', 'splicevault_site_type', 'splicevault_site_loc',
              'spliceai_delta', 'splicevault_out_of_frame_events',
              'splicevault_top4_events', 'splicevault_site_sample_count',
              'splicevault_site_max_depth']

    with open(output_tsv, 'w') as f:
        f.write('\t'.join(header) + '\n')

        for record in vcf:
            row = [
                record.chrom,
                str(record.pos),
                record.ref,
                record.alts[0],
                record.info['TRANSCRIPT_ID'],
                record.info['SITE_TYPE'],
                record.info['SITE_LOC'],
                str(record.info['SPLICEAI_DELTA']),
                record.info['OUT_OF_FRAME'],
                record.info['TOP4_EVENTS'],
                str(record.info['SAMPLE_COUNT']),
                str(record.info['MAX_DEPTH'])
            ]
            f.write('\t'.join(row) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Convert SpliceVault TSV between GRCh38 and hg19 coordinates')
    parser.add_argument('--input_tsv', required=True, help='Input TSV file (GRCh38)')
    parser.add_argument('--output_tsv', required=True, help='Output TSV file (hg19)')
    parser.add_argument('--chain_file', required=True, help='Chain file for liftover (hg38ToHg19.over.chain.gz)')
    parser.add_argument('--reference_fasta', required=True, help='Reference FASTA file (GRCh38)')
    parser.add_argument('--temp_dir', default=None, help='Directory for temporary files')

    args = parser.parse_args()

    # Create temporary directory if not specified
    temp_dir = args.temp_dir if args.temp_dir else os.path.dirname(args.output_tsv)

    # Intermediate file paths
    grch38_vcf = os.path.join(temp_dir, "grch38.vcf")
    hg19_vcf = os.path.join(temp_dir, "hg19.vcf")

    try:
        # 1. Convert TSV to VCF
        logger.info("Converting TSV to VCF...")
        if not os.path.exists(grch38_vcf):
            tsv_to_vcf(args.input_tsv, grch38_vcf, args.reference_fasta)

        # 2. Run CrossMap for coordinate liftover
        logger.info("Lifting over coordinates...")
        cmd = f"CrossMap vcf {args.chain_file} {grch38_vcf} {args.reference_fasta} {hg19_vcf}"
        ret = subprocess.run(cmd, shell=True, check=True)
        if ret.returncode != 0:
            raise Exception("CrossMap liftover failed")

        # 3. Convert lifted VCF back to TSV
        logger.info("Converting lifted VCF back to TSV...")
        vcf_to_tsv(hg19_vcf, args.output_tsv)
        logger.info(f"Conversion complete. Output written to {args.output_tsv}")
    except Exception as e:
        logger.error(f"Error: {e}")
        raise e


if __name__ == "__main__":
    main()