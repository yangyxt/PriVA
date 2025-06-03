#!/usr/bin/env python3

import sys
import os
import argparse
import multiprocessing as mp
import polars as pl
from pathlib import Path
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(module)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

# Define standard column sets
OUTPUT_HEADER_COLS = ["#Chrom", "Pos", "Ref", "Alt", "FeatureID", "PHRED"]
INTERNAL_PROCESSING_COLS = ["Chrom", "Pos", "Ref", "Alt", "FeatureID", "PHRED"]
SELF_SCRIPT_PATH = Path(__file__).parent.parent
CONTIG_MAP_PATH = os.path.join(SELF_SCRIPT_PATH, "data", "liftover", "GRC_to_ucsc.contig.map.txt")



def table_to_dict(file_path, separator=' '):
    """
    Efficiently convert a 2-column table to dictionary.
    
    Parameters:
    - file_path: Path to the table file
    - separator: Column delimiter (default: tab)
    - has_header: Whether the file has a header row
    - key_col: Name or index of key column (if None, uses first column)
    - val_col: Name or index of value column (if None, uses second column)
    
    Returns:
    - Dictionary mapping keys to values
    """
    result_dict = {}
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            line = line.strip()
            if not line:
                continue
            key, val = line.split(separator)
            result_dict[key] = val
    return result_dict



def read_and_standardize_chrom(file_path, 
                               separator='\t', 
                               has_header=True, 
                               comment_prefix='#', 
                               is_pos_file=False, 
                               null_values=["NA", "na", "NaN", "nan", "", ".", ",", ";", "NAN"]):
    """Reads a CSV/TSV and standardizes the chromosome column name."""
    if is_pos_file: # Position file from VCF has no header and fixed columns
        force_schema = {"column_1": pl.Utf8}
        df = pl.read_csv(file_path, separator=separator, has_header=False, comment_prefix=comment_prefix, null_values=null_values, infer_schema_length=10000, schema_overrides=force_schema)
        df = df.rename({
            "column_1": "Chrom",
            "column_2": "Pos",
            "column_3": "Ref",
            "column_4": "Alt"
        })
        # Select only these four columns for pos_file
        return df.select(["Chrom", "Pos", "Ref", "Alt"])

    force_schema = {"#Chrom": pl.Utf8, "Chrom": pl.Utf8, "column_1": pl.Utf8}
    df = pl.read_csv(file_path, 
                     separator=separator, 
                     has_header=has_header,
                     null_values=null_values, 
                     infer_schema_length=1000, 
                     schema_overrides=force_schema)
    if "#Chrom" in df.columns:
        df = df.rename({"#Chrom": "Chrom"})

    # pretty print the first 10 rows of df
    logger.info(f"First 10 rows of {file_path}: \n{df.head(10)}")

    # If chr prefix is absent from the first value of Chrom column, add it
    contig_map_dict = table_to_dict(CONTIG_MAP_PATH, separator=' ')
    df = df.with_columns([pl.col("Chrom").replace_strict(contig_map_dict, default=pl.col("Chrom")).alias("Chrom")])

    logger.info(f"First 10 rows of {file_path} after contig map: \n{df.head(10)}")
    # Ensure all internal processing columns (or all output header cols if more relevant) exist for hub/new files
    # For this generic function, we'll assume "Chrom" is the main one to standardize.
    # Specific column selection will happen in the calling functions.
    return df


def select_and_rename_for_output(df: pl.DataFrame) -> pl.DataFrame:
    """Selects standard columns and renames Chrom to #Chrom for output."""
    if "Chrom" not in df.columns:
        logger.error("Cannot rename 'Chrom' to '#Chrom': 'Chrom' column does not exist.")
        # Attempt to select other columns if they exist
        cols_to_select_anyway = [col for col in OUTPUT_HEADER_COLS if col in df.columns or col.replace("#","") in df.columns]
        if not cols_to_select_anyway:
             raise ValueError("DataFrame does not contain 'Chrom' or any recognizable output columns.")
        return df.select(cols_to_select_anyway) # Or handle error more gracefully

    # Ensure all internal columns are present before trying to select them with #Chrom rename
    # This expects df to have "Chrom" and other INTERNAL_PROCESSING_COLS
    select_exprs = []
    if "Chrom" in df.columns:
        select_exprs.append(pl.col("Chrom").alias("#Chrom"))
    for col_name in INTERNAL_PROCESSING_COLS:
        if col_name != "Chrom" and col_name in df.columns:
            select_exprs.append(pl.col(col_name))
    
    # If FeatureID or PHRED are missing from df, this select will fail.
    # It's better to select based on OUTPUT_HEADER_COLS, mapping from internal names.
    final_selection = []
    if "Chrom" in df.columns: # This should be the case if the function is called correctly
        final_selection.append(pl.col("Chrom").alias("#Chrom"))
    else: # Fallback if Chrom is already #Chrom or missing
        if "#Chrom" in df.columns:
            final_selection.append(pl.col("#Chrom"))

    for col_name in OUTPUT_HEADER_COLS:
        internal_name = col_name.replace("#", "")
        if col_name == "#Chrom": continue # Already handled
        if internal_name in df.columns:
            final_selection.append(pl.col(internal_name))
        elif col_name in df.columns: # If it was already #Chrom for example
             final_selection.append(pl.col(col_name))
    
    return df.select(final_selection)


def match_variants(hub_file, pos_file, covered_file, uncovered_file, threads=None, null_values=["NA", "na", "NaN", "nan", "", ".", ",", ";", "NAN"]):
    """Match variants in position file against hub file using Polars for efficiency"""
    # Threads argument is kept for signature compatibility, Polars handles its own threading.
    # if threads is None:
    #     threads = mp.cpu_count() # Polars uses this by default

    # Check if hub file exists and has content
    if not os.path.exists(hub_file) or os.path.getsize(hub_file) == 0:
        variants_df = read_and_standardize_chrom(pos_file, is_pos_file=True)
        
        with open(uncovered_file, 'w') as f:
            for row in variants_df.select(["Chrom", "Pos"]).iter_rows():
                f.write(f"{row[0]}\t{row[1]}\n")
        
        with open(covered_file, 'w') as f:
            f.write("\t".join(OUTPUT_HEADER_COLS) + "\n")
        
        logger.warning(f"Hub file {hub_file} missing or empty, all variants are uncovered")
        return False # Indicates not all variants were covered (in this case, none were from hub)
    
    variants_df = read_and_standardize_chrom(pos_file, is_pos_file=True)
    variants_df = variants_df.with_columns(
        pl.concat_str([
            pl.col("Chrom").cast(pl.Utf8),
            pl.col("Pos").cast(pl.Utf8),
            pl.col("Ref").cast(pl.Utf8),
            pl.col("Alt").cast(pl.Utf8)
        ], separator="_").alias("var_key")
    )
    
    hub_df_raw = read_and_standardize_chrom(hub_file, comment_prefix='#')
    # Select only the necessary columns for internal processing from the hub
    try:
        hub_df = hub_df_raw.select(INTERNAL_PROCESSING_COLS)
    except pl.ColumnNotFoundError as e:
        logger.error(f"Hub file {hub_file} is missing required columns: {e}. Cannot proceed with matching.")
        # Fallback: treat as if hub was empty
        with open(uncovered_file, 'w') as f:
            for row in variants_df.select(["Chrom", "Pos"]).iter_rows():
                f.write(f"{row[0]}\t{row[1]}\n")
        with open(covered_file, 'w') as f:
            f.write("\t".join(OUTPUT_HEADER_COLS) + "\n")
        return False


    hub_df = hub_df.with_columns(
        pl.concat_str([
            pl.col("Chrom").cast(pl.Utf8),
            pl.col("Pos").cast(pl.Utf8),
            pl.col("Ref").cast(pl.Utf8),
            pl.col("Alt").cast(pl.Utf8)
        ], separator="_").alias("var_key")
    )
    
    # hub_df is already a DataFrame. The .collect() call was an error and is removed.
    
    variant_keys = set(variants_df["var_key"].to_list())
    hub_keys = set(hub_df["var_key"].to_list())
    
    covered_keys = variant_keys.intersection(hub_keys)
    uncovered_keys = variant_keys - covered_keys
    
    # Filter covered variants from hub_df (which has internal column names)
    covered_df_internal = hub_df.filter(pl.col("var_key").is_in(list(covered_keys))) # Use list for safety
    
    temp_covered = f"{covered_file}.temp.{os.getpid()}"
    if not covered_df_internal.is_empty():
        covered_df_output = select_and_rename_for_output(covered_df_internal.drop("var_key"))
        covered_df_output.write_csv(temp_covered, separator='\t')
    else: # Create empty temp file if no covered variants
        with open(temp_covered, 'w') as f:
            f.write("\t".join(OUTPUT_HEADER_COLS) + "\n")


    if os.path.exists(temp_covered): # Check if temp file was created (even if empty with header)
        if os.path.getsize(temp_covered) > 0 or not covered_df_internal.is_empty(): # Ensure non-empty content if variants were expected
            os.rename(temp_covered, covered_file)
        elif covered_df_internal.is_empty() and os.path.getsize(temp_covered) == len("\t".join(OUTPUT_HEADER_COLS) + "\n"):
             os.rename(temp_covered, covered_file) # It's an empty file with just header
        else:
            logger.error(f"Error: Failed to create valid covered file. Temp file size: {os.path.getsize(temp_covered)}")
            if os.path.exists(temp_covered): os.remove(temp_covered)
            # Fallback to creating a minimal covered file
            with open(covered_file, 'w') as f:
                f.write("\t".join(OUTPUT_HEADER_COLS) + "\n")
    else: # Should not happen if logic above is correct
        logger.error("Error: Temp covered file was not created.")
        with open(covered_file, 'w') as f:
            f.write("\t".join(OUTPUT_HEADER_COLS) + "\n")

    uncovered_pos_df = variants_df.filter(pl.col("var_key").is_in(list(uncovered_keys)))
    
    temp_uncovered = f"{uncovered_file}.temp.{os.getpid()}"
    with open(temp_uncovered, 'w') as f:
        for row in uncovered_pos_df.select(["Chrom", "Pos"]).iter_rows():
            f.write(f"{row[0]}\t{row[1]}\n")
    
    # Ensure uncovered_file is always created, even if empty
    if os.path.exists(temp_uncovered) and (os.path.getsize(temp_uncovered) > 0 or uncovered_pos_df.is_empty()):
        os.rename(temp_uncovered, uncovered_file)
    else:
        logger.warning(f"No uncovered positions to write, or failed to write. Creating empty uncovered file: {uncovered_file}")
        with open(uncovered_file, 'w') as f: # Create empty file
            pass
        if os.path.exists(temp_uncovered) and os.path.getsize(temp_uncovered) == 0: # if temp was empty, remove it
            os.remove(temp_uncovered)
        elif os.path.exists(temp_uncovered): # if temp existed with error
             logger.error(f"Failed to create valid uncovered file from temp file.")
             os.remove(temp_uncovered)


    total = len(variant_keys)
    covered = len(covered_keys)
    uncovered_count = len(uncovered_keys) # Use a different name to avoid conflict with file name
    logger.info(f"Found {covered} cached variants, {uncovered_count} uncached variants from {total} input variants.")
    
    return uncovered_count == 0


def update_hub(new_file_path, hub_file_path, threads=None, null_values=["NA", "na", "NaN", "nan", "", ".", ",", ";", "NAN"]):
    """Update hub CADD TSV with new scores using Polars, with atomic file operations"""
    if not os.path.exists(new_file_path) or os.path.getsize(new_file_path) == 0:
        logger.info(f"New CADD file {new_file_path} is empty or doesn't exist, no update needed")
        return
    
    hub_dir = os.path.dirname(hub_file_path)
    if hub_dir: # Create dir only if hub_file_path includes a directory part
        os.makedirs(hub_dir, exist_ok=True)
    
    temp_hub_file = f"{hub_file_path}.temp.{os.getpid()}"
    
    new_df_raw = read_and_standardize_chrom(new_file_path, comment_prefix='#')
    try:
        new_df_internal = new_df_raw.select(INTERNAL_PROCESSING_COLS)
    except pl.ColumnNotFoundError as e:
        logger.error(f"New file {new_file_path} is missing required columns: {e}. Cannot update hub.")
        return

    if not os.path.exists(hub_file_path) or os.path.getsize(hub_file_path) == 0:
        logger.info(f"Hub file {hub_file_path} does not exist or is empty. Creating new hub from {new_file_path}.")
        final_hub_df = select_and_rename_for_output(new_df_internal)
        final_hub_df.write_csv(temp_hub_file, separator='\t')
        
        if os.path.exists(temp_hub_file) and os.path.getsize(temp_hub_file) > 0:
            os.rename(temp_hub_file, hub_file_path)
            logger.info(f"Created new CADD hub TSV {hub_file_path} with {len(final_hub_df)} variants.")
        else:
            logger.error(f"Error: Failed to create valid temp file for new hub {hub_file_path}.")
            if os.path.exists(temp_hub_file): os.remove(temp_hub_file)
        return

    hub_df_raw = read_and_standardize_chrom(hub_file_path, comment_prefix='#')
    try:
        hub_df_internal = hub_df_raw.select(INTERNAL_PROCESSING_COLS)
    except pl.ColumnNotFoundError as e:
        logger.error(f"Existing hub file {hub_file_path} is missing required columns: {e}. Cannot update. Consider re-creating the hub.")
        return
        
    original_size = len(hub_df_internal)
    
    hub_df_internal = hub_df_internal.with_columns(
        pl.concat_str([pl.col(c).cast(pl.Utf8) for c in ["Chrom", "Pos", "Ref", "Alt"]], separator="_").alias("var_key")
    )
    new_df_internal = new_df_internal.with_columns(
        pl.concat_str([pl.col(c).cast(pl.Utf8) for c in ["Chrom", "Pos", "Ref", "Alt"]], separator="_").alias("var_key")
    )
    
    hub_keys = set(hub_df_internal["var_key"].to_list())
    new_variants_to_add = new_df_internal.filter(~pl.col("var_key").is_in(list(hub_keys)))
    
    if new_variants_to_add.is_empty():
        logger.info(f"No new unique variants from {new_file_path} to add to hub {hub_file_path}.")
        return
    
    combined_df_internal = pl.concat([
        hub_df_internal.drop("var_key"), 
        new_variants_to_add.drop("var_key")
    ])
    
    final_combined_df_output = select_and_rename_for_output(combined_df_internal)
    final_combined_df_output.write_csv(temp_hub_file, separator='\t')
    
    if os.path.exists(temp_hub_file) and os.path.getsize(temp_hub_file) > 0:
        try:
            temp_df_check = pl.read_csv(temp_hub_file, separator='\t', has_header=True, comment_prefix='#', null_values=null_values, infer_schema_length=10000)
            if len(temp_df_check) >= original_size : # Should be strictly greater if new_variants_to_add was not empty
                os.rename(temp_hub_file, hub_file_path)
                logger.info(f"Added {len(new_variants_to_add)} new variants to CADD hub {hub_file_path}. New total: {len(temp_df_check)}.")
            else:
                logger.error(f"Error: Temp hub file ({len(temp_df_check)} entries) incorrect size. Original: {original_size}, Added: {len(new_variants_to_add)}. Hub not updated.")
                os.remove(temp_hub_file)
        except Exception as e:
            logger.error(f"Error verifying temp hub file: {str(e)}")
            if os.path.exists(temp_hub_file): os.remove(temp_hub_file)
    else:
        logger.error(f"Error: Failed to create valid temp file for hub update: {temp_hub_file}")
        if os.path.exists(temp_hub_file): os.remove(temp_hub_file)


def merge_results(covered_file_path, new_file_path, output_file_path, null_values=["NA", "na", "NaN", "nan", "", ".", ",", ";", "NAN"]):
    """Merge covered and new CADD results using Polars with atomic file updates"""
    covered_exists = os.path.exists(covered_file_path) and os.path.getsize(covered_file_path) > 0
    new_exists = os.path.exists(new_file_path) and os.path.getsize(new_file_path) > 0
    
    output_dir = os.path.dirname(output_file_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    temp_output_file = f"{output_file_path}.temp.{os.getpid()}"

    dfs_to_concat = []
    len_covered = 0
    len_new = 0

    if covered_exists:
        try:
            covered_df_raw = read_and_standardize_chrom(covered_file_path, comment_prefix='#')
            # Covered file should already be in output format with #Chrom
            # but we select to ensure order and existence of columns
            if "#Chrom" not in covered_df_raw.columns and "Chrom" in covered_df_raw.columns: # If it was stored with "Chrom"
                 covered_df_raw = covered_df_raw.rename({"Chrom": "#Chrom"})

            # Select based on OUTPUT_HEADER_COLS, ensuring they exist
            present_cols_covered = [col for col in OUTPUT_HEADER_COLS if col in covered_df_raw.columns]
            covered_df = covered_df_raw.select(present_cols_covered)

            if not covered_df.is_empty():
                dfs_to_concat.append(covered_df)
                len_covered = len(covered_df)
        except Exception as e:
            logger.warning(f"Could not read or process covered file {covered_file_path}: {e}. Skipping.")
    
    if new_exists:
        try:
            new_df_raw = read_and_standardize_chrom(new_file_path, comment_prefix='#')
            if "#Chrom" not in new_df_raw.columns and "Chrom" in new_df_raw.columns: # If it was stored with "Chrom"
                 new_df_raw = new_df_raw.rename({"Chrom": "#Chrom"})
            
            present_cols_new = [col for col in OUTPUT_HEADER_COLS if col in new_df_raw.columns]
            new_df = new_df_raw.select(present_cols_new)

            if not new_df.is_empty():
                dfs_to_concat.append(new_df)
                len_new = len(new_df)
        except Exception as e:
            logger.warning(f"Could not read or process new scores file {new_file_path}: {e}. Skipping.")

    if not dfs_to_concat:
        logger.error("ERROR: No CADD scores available to merge. Creating empty output file.")
        with open(output_file_path, 'w') as f: # Create empty file with header
            f.write("\t".join(OUTPUT_HEADER_COLS) + "\n")
        return False # No successful merge

    # Polars concat can handle schema differences by name alignment, filling missing with nulls.
    # By selecting common OUTPUT_HEADER_COLS, we aim for consistency.
    try:
        combined_df = pl.concat(dfs_to_concat, how="diagonal_relaxed") # Changed 'diagonal' to 'diagonal_relaxed' in recent polars
    except TypeError: # Fallback for slightly older Polars versions if 'diagonal_relaxed' isn't there
        try:
            combined_df = pl.concat(dfs_to_concat, how="diagonal")
        except TypeError: # Fallback to default if diagonal strategies are not available/suitable
            combined_df = pl.concat(dfs_to_concat)


    # Ensure final output has the exact OUTPUT_HEADER_COLS in order, if possible
    final_cols_present_in_combined = [col for col in OUTPUT_HEADER_COLS if col in combined_df.columns]
    combined_df_final_output = combined_df.select(final_cols_present_in_combined)


    combined_df_final_output.write_csv(temp_output_file, separator='\t')
    
    if os.path.exists(temp_output_file) and os.path.getsize(temp_output_file) > 0:
        try:
            temp_df_check = pl.read_csv(temp_output_file, separator='\t', has_header=True, comment_prefix='#', null_values=null_values, infer_schema_length=10000)
            # Check if the merged file seems reasonable (not empty if input wasn't empty)
            if not temp_df_check.is_empty() or (len_covered == 0 and len_new == 0) :
                os.rename(temp_output_file, output_file_path)
                logger.info(f"Merged {len_covered} cached and {len_new} newly calculated CADD scores into {output_file_path}. Total rows: {len(temp_df_check)}.")
                return True
            else:
                logger.error(f"Error: Merged temp file {temp_output_file} is unexpectedly empty. Original covered={len_covered}, new={len_new}.")
                os.remove(temp_output_file)
                # Create empty file with header as fallback
                with open(output_file_path, 'w') as f:
                    f.write("\t".join(OUTPUT_HEADER_COLS) + "\n")
                return False
        except Exception as e:
            logger.error(f"Error verifying merged temp file: {str(e)}")
            if os.path.exists(temp_output_file): os.remove(temp_output_file)
            with open(output_file_path, 'w') as f:
                f.write("\t".join(OUTPUT_HEADER_COLS) + "\n")
            return False
    else:
        logger.error(f"Error: Failed to create valid temp file for merge: {temp_output_file}")
        if os.path.exists(temp_output_file): os.remove(temp_output_file)
        with open(output_file_path, 'w') as f:
            f.write("\t".join(OUTPUT_HEADER_COLS) + "\n")
        return False


def main():
    parser = argparse.ArgumentParser(description='CADD score caching utilities')
    # Make subcommand mandatory
    subparsers = parser.add_subparsers(dest='command', help='Command to run', required=True)
    
    # Match variants subcommand
    match_parser = subparsers.add_parser('match', help='Match variants against hub CADD TSV')
    match_parser.add_argument('--hub', required=True, help='Hub CADD TSV file')
    match_parser.add_argument('--pos', required=True, help='Position file from VCF (Chrom Pos Ref Alt, no header)')
    match_parser.add_argument('--covered', required=True, help='Output file for covered variants (TSV with header)')
    match_parser.add_argument('--uncovered', required=True, help='Output file for uncovered positions (Chrom Pos, for bcftools)')
    match_parser.add_argument('--threads', type=int, default=None, help='Number of threads (Polars usually auto-detects optimal).')
    
    # Update hub subcommand
    update_parser = subparsers.add_parser('update', help='Update hub CADD TSV with new scores')
    update_parser.add_argument('--new', required=True, help='New CADD TSV file with scores to add')
    update_parser.add_argument('--hub', required=True, help='Hub CADD TSV file to update')
    update_parser.add_argument('--threads', type=int, default=None, help='Number of threads.')
    
    # Merge results subcommand
    merge_parser = subparsers.add_parser('merge', help='Merge covered and new CADD results')
    merge_parser.add_argument('--covered', required=True, help='Covered CADD TSV file (from cache)')
    merge_parser.add_argument('--new', required=True, help='New CADD TSV file (newly calculated scores)')
    merge_parser.add_argument('--output', required=True, help='Output merged TSV file')
    # Threads not used directly by merge logic's polars calls beyond default polars behavior
    
    args = parser.parse_args()

    # Polars global settings (optional, as Polars often defaults well)
    # if args.threads is not None and args.threads > 0:
    #    try:
    #        pl.set_num_threads(args.threads) # Available in modern Polars
    #    except AttributeError:
    #        logger.warning("pl.set_num_threads not available. Ensure POLARS_NUM_THREADS env var if specific thread count needed.")
    # pl.enable_string_cache(True) # Generally recommended

    if args.command == 'match':
        # The return of match_variants is True if all variants were covered.
        # Exit code 0 usually means success of the script's operation,
        # not necessarily that all variants were covered.
        # Let's adjust to exit 0 if the command ran without catastrophic error.
        try:
            match_variants(args.hub, args.pos, args.covered, args.uncovered, args.threads)
            sys.exit(0)
        except Exception as e:
            logger.critical(f"Critical error in 'match' command: {e}", exc_info=True)
            sys.exit(1)
    elif args.command == 'update':
        try:
            update_hub(args.new, args.hub, args.threads)
            sys.exit(0)
        except Exception as e:
            logger.critical(f"Critical error in 'update' command: {e}", exc_info=True)
            sys.exit(1)
    elif args.command == 'merge':
        try:
            success_flag = merge_results(args.covered, args.new, args.output)
            sys.exit(0 if success_flag else 1) # Exit 1 if merge reported an issue (e.g. no files)
        except Exception as e:
            logger.critical(f"Critical error in 'merge' command: {e}", exc_info=True)
            sys.exit(1)
    else:
        parser.print_help() # Should not be reached if subparsers are 'required'
        sys.exit(1)

if __name__ == "__main__":
    main()