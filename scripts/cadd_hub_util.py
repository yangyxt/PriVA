#!/usr/bin/env python3

import sys
import os
import argparse
import multiprocessing as mp
import polars as pl
from pathlib import Path

def match_variants(hub_file, pos_file, covered_file, uncovered_file, threads=None):
    """Match variants in position file against hub file using Polars for efficiency"""
    if threads is None:
        threads = mp.cpu_count()

    # Check if hub file exists and has content
    if not os.path.exists(hub_file) or os.path.getsize(hub_file) == 0:
        # Write all variants to uncovered file
        variants_df = pl.read_csv(pos_file, separator='\t', has_header=False)
        variants_df = variants_df.rename({
            "column_1": "Chrom",
            "column_2": "Pos"
        })
        
        # Write uncovered positions - ensuring 1-indexed positions for bcftools -R
        # The format is CHROM\tPOS, which bcftools treats as 1-indexed positions
        with open(uncovered_file, 'w') as f:
            for row in variants_df.select(["Chrom", "Pos"]).iter_rows():
                f.write(f"{row[0]}\t{row[1]}\n")
        
        # Create empty covered file with header
        with open(covered_file, 'w') as f:
            f.write("#Chrom\tPos\tRef\tAlt\tFeatureID\tPHRED\n")
        
        print(f"Hub file {hub_file} missing or empty, all variants are uncovered", file=sys.stderr)
        return False
    
    # Read position file with variant data
    variants_df = pl.read_csv(pos_file, separator='\t', has_header=False)
    variants_df = variants_df.rename({
        "column_1": "Chrom",
        "column_2": "Pos",
        "column_3": "Ref",
        "column_4": "Alt"
    })
    
    # Create lookup key
    variants_df = variants_df.with_column(
        pl.concat_str([
            pl.col("Chrom"),
            pl.col("Pos").cast(pl.Utf8),
            pl.col("Ref"),
            pl.col("Alt")
        ], separator="_").alias("var_key")
    )
    
    # Read hub file with lazy evaluation - only needed columns
    hub_df = pl.scan_csv(
        hub_file, 
        separator='\t',
        has_header=True,
        comment_char='#'
    ).select([
        pl.col("#Chrom").alias("Chrom"),
        pl.col("Pos"),
        pl.col("Ref"),
        pl.col("Alt"),
        pl.col("FeatureID"),
        pl.col("PHRED")
    ])
    
    # Create lookup key for hub
    hub_df = hub_df.with_column(
        pl.concat_str([
            pl.col("Chrom"),
            pl.col("Pos").cast(pl.Utf8),
            pl.col("Ref"),
            pl.col("Alt")
        ], separator="_").alias("var_key")
    )
    
    # Collect hub data with key for faster lookup
    hub_df = hub_df.collect()
    
    # Create sets of variant keys for fast lookup
    variant_keys = set(variants_df["var_key"].to_list())
    hub_keys = set(hub_df["var_key"].to_list())
    
    # Find covered variants
    covered_keys = variant_keys.intersection(hub_keys)
    uncovered_keys = variant_keys - covered_keys
    
    # Filter covered variants from hub
    covered_df = hub_df.filter(pl.col("var_key").is_in(covered_keys))
    
    # Write covered variants to file
    # Write to temp file first and then move
    temp_covered = f"{covered_file}.temp"
    covered_df.drop("var_key").write_csv(temp_covered, separator='\t')
    
    if os.path.exists(temp_covered) and os.path.getsize(temp_covered) > 0:
        os.rename(temp_covered, covered_file)
    else:
        print(f"Error: Failed to create valid covered file", file=sys.stderr)
        if os.path.exists(temp_covered):
            os.remove(temp_covered)
    
    # Get uncovered variants
    uncovered_df = variants_df.filter(pl.col("var_key").is_in(uncovered_keys))
    
    # Write uncovered positions to file for bcftools
    # The format is CHROM\tPOS, which bcftools treats as 1-indexed positions
    temp_uncovered = f"{uncovered_file}.temp"
    with open(temp_uncovered, 'w') as f:
        for row in uncovered_df.select(["Chrom", "Pos"]).iter_rows():
            f.write(f"{row[0]}\t{row[1]}\n")
    
    if os.path.exists(temp_uncovered) and os.path.getsize(temp_uncovered) > 0:
        os.rename(temp_uncovered, uncovered_file)
    else:
        print(f"Error: Failed to create valid uncovered file", file=sys.stderr)
        if os.path.exists(temp_uncovered):
            os.remove(temp_uncovered)
    
    # Print stats
    total = len(variant_keys)
    covered = len(covered_keys)
    uncovered = len(uncovered_keys)
    print(f"Found {covered} cached variants, {uncovered} uncached variants", file=sys.stderr)
    
    # Return True if all variants were covered
    return uncovered == 0

def update_hub(new_file, hub_file, threads=None):
    """Update hub CADD TSV with new scores using Polars, with atomic file operations"""
    if threads is None:
        threads = mp.cpu_count()
    
    # Skip if new file doesn't exist or is empty
    if not os.path.exists(new_file) or os.path.getsize(new_file) == 0:
        print("New CADD file is empty or doesn't exist, no update needed", file=sys.stderr)
        return
    
    # Ensure hub directory exists
    hub_dir = os.path.dirname(hub_file)
    os.makedirs(hub_dir, exist_ok=True)
    
    # Create temp file in same directory to ensure atomic move
    temp_file = f"{hub_file}.temp"
    
    # If hub file doesn't exist, just copy the new file with only needed columns
    if not os.path.exists(hub_file) or os.path.getsize(hub_file) == 0:
        new_df = pl.read_csv(
            new_file, 
            separator='\t', 
            has_header=True, 
            comment_char='#'
        ).select([
            pl.col("#Chrom").alias("Chrom"),
            pl.col("Pos"),
            pl.col("Ref"),
            pl.col("Alt"),
            pl.col("FeatureID"),
            pl.col("PHRED")
        ])
        
        # Check if we need to rename the first column
        if "#Chrom" not in new_df.columns and "Chrom" in new_df.columns:
            new_df = new_df.rename({"Chrom": "#Chrom"})
        
        # Write to temp file first
        new_df.write_csv(temp_file, separator='\t')
        
        # Check if temp file is valid
        if os.path.exists(temp_file) and os.path.getsize(temp_file) > 0:
            # Move temp file to hub file (atomic operation)
            os.rename(temp_file, hub_file)
            print(f"Created new CADD hub TSV with {len(new_df)} variants", file=sys.stderr)
        else:
            print(f"Error: Failed to create valid temp file", file=sys.stderr)
            if os.path.exists(temp_file):
                os.remove(temp_file)
        return
    
    # Read existing hub data
    hub_df = pl.read_csv(
        hub_file, 
        separator='\t', 
        has_header=True,
        comment_char='#'
    ).select([
        pl.col("#Chrom").alias("Chrom"),
        pl.col("Pos"),
        pl.col("Ref"),
        pl.col("Alt"),
        pl.col("FeatureID"),
        pl.col("PHRED")
    ])
    
    original_size = len(hub_df)
    
    # Create lookup key for hub
    hub_df = hub_df.with_column(
        pl.concat_str([
            pl.col("Chrom"),
            pl.col("Pos").cast(pl.Utf8),
            pl.col("Ref"),
            pl.col("Alt")
        ], separator="_").alias("var_key")
    )
    
    # Read new data
    new_df = pl.read_csv(
        new_file, 
        separator='\t', 
        has_header=True,
        comment_char='#'
    ).select([
        pl.col("#Chrom").alias("Chrom"),
        pl.col("Pos"),
        pl.col("Ref"),
        pl.col("Alt"),
        pl.col("FeatureID"),
        pl.col("PHRED")
    ])
    
    # Create lookup key for new data
    new_df = new_df.with_column(
        pl.concat_str([
            pl.col("Chrom"),
            pl.col("Pos").cast(pl.Utf8),
            pl.col("Ref"),
            pl.col("Alt")
        ], separator="_").alias("var_key")
    )
    
    # Get hub keys for filtering
    hub_keys = set(hub_df["var_key"].to_list())
    
    # Filter new variants not in hub
    new_variants = new_df.filter(~pl.col("var_key").is_in(hub_keys))
    
    # If no new variants, we're done
    if len(new_variants) == 0:
        print("No new variants to add to hub", file=sys.stderr)
        return
    
    # Prepare new variants for append (drop var_key and rename Chrom back to #Chrom)
    new_variants = new_variants.drop("var_key")
    if "Chrom" in new_variants.columns and "#Chrom" not in new_variants.columns:
        new_variants = new_variants.rename({"Chrom": "#Chrom"})
    
    # Create the combined DataFrame
    combined_df = pl.concat([hub_df.drop("var_key"), new_variants])
    
    # Write to temp file first
    combined_df.write_csv(temp_file, separator='\t')
    
    # Verify the temp file is valid and larger than original
    if os.path.exists(temp_file) and os.path.getsize(temp_file) > 0:
        # Additional verification that the file is valid and has more entries
        try:
            # Check if the new file has more entries
            temp_df = pl.read_csv(temp_file, separator='\t', has_header=True, comment_char='#')
            if len(temp_df) > original_size:
                # Move temp file to hub file (atomic operation)
                os.rename(temp_file, hub_file)
                print(f"Added {len(new_variants)} new variants to CADD hub TSV", file=sys.stderr)
            else:
                print(f"Error: Temp file doesn't contain more entries than original", file=sys.stderr)
                os.remove(temp_file)
        except Exception as e:
            print(f"Error verifying temp file: {str(e)}", file=sys.stderr)
            if os.path.exists(temp_file):
                os.remove(temp_file)
    else:
        print(f"Error: Failed to create valid temp file", file=sys.stderr)
        if os.path.exists(temp_file):
            os.remove(temp_file)

def merge_results(covered_file, new_file, output_file):
    """Merge covered and new CADD results using Polars with atomic file updates"""
    # Check if files exist and have content
    covered_exists = os.path.exists(covered_file) and os.path.getsize(covered_file) > 0
    new_exists = os.path.exists(new_file) and os.path.getsize(new_file) > 0
    
    if not covered_exists and not new_exists:
        print("ERROR: No CADD scores available (neither covered nor newly calculated)", file=sys.stderr)
        return False
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Create temp file in same directory for atomic operation
    temp_file = f"{output_file}.temp"
    
    # Handle edge cases
    if not covered_exists:
        print("No cached variants, using only newly calculated CADD scores", file=sys.stderr)
        # Just copy the new file
        new_df = pl.read_csv(new_file, separator='\t', has_header=True, comment_char='#')
        new_df.write_csv(temp_file, separator='\t')
        
        if os.path.exists(temp_file) and os.path.getsize(temp_file) > 0:
            os.rename(temp_file, output_file)
            return True
        else:
            print("Error: Failed to create valid output file", file=sys.stderr)
            if os.path.exists(temp_file):
                os.remove(temp_file)
            return False
    
    if not new_exists:
        print("No new CADD scores, using only cached variants", file=sys.stderr)
        # Just copy the covered file
        covered_df = pl.read_csv(covered_file, separator='\t', has_header=True, comment_char='#')
        covered_df.write_csv(temp_file, separator='\t')
        
        if os.path.exists(temp_file) and os.path.getsize(temp_file) > 0:
            os.rename(temp_file, output_file)
            return True
        else:
            print("Error: Failed to create valid output file", file=sys.stderr)
            if os.path.exists(temp_file):
                os.remove(temp_file)
            return False
    
    # Read both files
    covered_df = pl.read_csv(covered_file, separator='\t', has_header=True, comment_char='#')
    new_df = pl.read_csv(new_file, separator='\t', has_header=True, comment_char='#')
    
    # Standardize column names between the two DataFrames
    if "#Chrom" in covered_df.columns and "Chrom" in new_df.columns:
        new_df = new_df.rename({"Chrom": "#Chrom"})
    elif "Chrom" in covered_df.columns and "#Chrom" in new_df.columns:
        new_df = new_df.rename({"#Chrom": "Chrom"})
    
    # Concatenate the DataFrames
    combined_df = pl.concat([covered_df, new_df])
    
    # Write to temp file first
    combined_df.write_csv(temp_file, separator='\t')
    
    # Verify the temp file is valid
    if os.path.exists(temp_file) and os.path.getsize(temp_file) > 0:
        try:
            # Check if the temp file has all entries
            temp_df = pl.read_csv(temp_file, separator='\t', has_header=True, comment_char='#')
            if len(temp_df) == len(covered_df) + len(new_df):
                # Move temp file to output file (atomic operation)
                os.rename(temp_file, output_file)
                print(f"Merged {len(covered_df)} cached and {len(new_df)} newly calculated CADD scores", file=sys.stderr)
                return True
            else:
                print(f"Error: Temp file doesn't contain the expected number of entries", file=sys.stderr)
                os.remove(temp_file)
                return False
        except Exception as e:
            print(f"Error verifying temp file: {str(e)}", file=sys.stderr)
            if os.path.exists(temp_file):
                os.remove(temp_file)
            return False
    else:
        print(f"Error: Failed to create valid temp file", file=sys.stderr)
        if os.path.exists(temp_file):
            os.remove(temp_file)
        return False

def main():
    parser = argparse.ArgumentParser(description='CADD score caching utilities')
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Match variants subcommand
    match_parser = subparsers.add_parser('match', help='Match variants against hub CADD TSV')
    match_parser.add_argument('--hub', required=True, help='Hub CADD TSV file')
    match_parser.add_argument('--pos', required=True, help='Position file from VCF')
    match_parser.add_argument('--covered', required=True, help='Output file for covered variants')
    match_parser.add_argument('--uncovered', required=True, help='Output file for uncovered positions')
    match_parser.add_argument('--threads', type=int, default=mp.cpu_count(), help='Number of threads to use')
    
    # Update hub subcommand
    update_parser = subparsers.add_parser('update', help='Update hub CADD TSV with new scores')
    update_parser.add_argument('--new', required=True, help='New CADD TSV file')
    update_parser.add_argument('--hub', required=True, help='Hub CADD TSV file')
    update_parser.add_argument('--threads', type=int, default=mp.cpu_count(), help='Number of threads to use')
    
    # Merge results subcommand
    merge_parser = subparsers.add_parser('merge', help='Merge covered and new CADD results')
    merge_parser.add_argument('--covered', required=True, help='Covered CADD TSV file')
    merge_parser.add_argument('--new', required=True, help='New CADD TSV file')
    merge_parser.add_argument('--output', required=True, help='Output merged TSV file')
    
    args = parser.parse_args()
    
    if args.command == 'match':
        success = match_variants(args.hub, args.pos, args.covered, args.uncovered, args.threads)
        sys.exit(0 if success else 1)
    elif args.command == 'update':
        update_hub(args.new, args.hub, args.threads)
        sys.exit(0)
    elif args.command == 'merge':
        success = merge_results(args.covered, args.new, args.output)
        sys.exit(0 if success else 1)
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()