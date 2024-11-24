import pandas as pd
import logging
import sys

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def read_phenotype_data(input_file):
    """
    Read phenotype annotation data from TSV file.
    
    Args:
        input_file (str): Path to input TSV file
    
    Returns:
        pd.DataFrame: Raw phenotype annotation data
    """
    try:
        df = pd.read_csv(input_file, sep='\t')
        logger.info(f"Successfully read input file with {len(df)} rows")
        return df
    except Exception as e:
        logger.error(f"Error reading input file: {str(e)}")
        raise

def separate_inheritance_data(df):
    """
    Separate inheritance patterns from other phenotype annotations.
    Use regex to find patterns containing "inheritance" at the end of the string
    
    Args:
        df (pd.DataFrame): Input phenotype annotation dataframe
    
    Returns:
        tuple: (main_df, inheritance_df) - separated dataframes
    """
    try:
        inheritance_record_bool = df['hpo_name'].str.contains(".+ inheritance$", regex=True)
        inheritance_df = df[inheritance_record_bool].copy()
        main_df = df[~inheritance_record_bool].copy()
        
        logger.info(f"Separated {len(inheritance_df)} inheritance rows from {len(main_df)} phenotype rows")
        return main_df, inheritance_df
    except Exception as e:
        logger.error(f"Error separating inheritance data: {str(e)}")
        raise

def deduplicate_paired_values(row, separator=';'):
    """
    Deduplicate values while maintaining correspondence between paired columns.
    
    Args:
        row (pd.Series): Row containing paired values to deduplicate
        separator (str): Separator used in string concatenation
        
    Returns:
        pd.Series: Row with deduplicated paired values
    """
    try:
        # Create pairs of values
        pairs = list(zip(
            row['hpo_id'].split(separator),
            row['hpo_name'].split(separator),
            row['frequency'].split(separator),
            row['disease_id'].split(separator)
        ))
        
        # Remove duplicates while preserving order
        unique_pairs = []
        seen = set()
        for pair in pairs:
            # Use first two elements (hpo_id and hpo_name) as key for deduplication
            key = (pair[0], pair[1])
            if key not in seen:
                seen.add(key)
                unique_pairs.append(pair)
        
        # Unzip the pairs back into separate lists
        if unique_pairs:
            hpo_ids, hpo_names, frequencies, disease_ids = zip(*unique_pairs)
        else:
            hpo_ids, hpo_names, frequencies, disease_ids = [], [], [], []
        
        # Update the row with deduplicated values
        row['hpo_id'] = separator.join(hpo_ids)
        row['hpo_name'] = separator.join(hpo_names)
        row['frequency'] = separator.join(frequencies)
        row['disease_id'] = separator.join(disease_ids)
        
        return row
    except Exception as e:
        logger.error(f"Error deduplicating paired values: {str(e)}")
        raise

def collapse_main_annotations(df):
    """
    Collapse main phenotype annotations by gene with deduplication.
    
    Args:
        df (pd.DataFrame): Main phenotype annotation dataframe
    
    Returns:
        pd.DataFrame: Collapsed and deduplicated main annotations
    """
    try:
        # First collapse all annotations
        collapsed = df.groupby(['ncbi_gene_id', 'gene_symbol']).agg({
            'hpo_id': ';'.join,
            'hpo_name': ';'.join,
            'frequency': ';'.join,
            'disease_id': ';'.join
        }).reset_index()
        
        # Then deduplicate while maintaining correspondence
        collapsed = collapsed.apply(deduplicate_paired_values, axis=1)
        
        logger.info(f"Collapsed and deduplicated annotations into {len(collapsed)} gene entries")
        return collapsed
    except Exception as e:
        logger.error(f"Error collapsing main annotations: {str(e)}")
        raise

def process_inheritance_data(inheritance_df):
    """
    Process inheritance patterns into separate columns with deduplication.
    
    Args:
        inheritance_df (pd.DataFrame): Inheritance pattern dataframe
    
    Returns:
        pd.DataFrame: Processed and deduplicated inheritance data
    """
    try:
        # First group by gene
        inheritance_modes = inheritance_df.groupby('ncbi_gene_id').agg({
            'hpo_name': ';'.join,
            'disease_id': ';'.join
        }).reset_index()
        
        # Deduplicate inheritance modes while maintaining disease ID correspondence
        def deduplicate_inheritance(row):
            pairs = list(zip(
                row['hpo_name'].split(';'),
                row['disease_id'].split(';')
            ))
            unique_pairs = []
            seen = set()
            for pair in pairs:
                if pair[0] not in seen:  # Use inheritance mode as key
                    seen.add(pair[0])
                    unique_pairs.append(pair)
            
            if unique_pairs:
                modes, disease_ids = zip(*unique_pairs)
            else:
                modes, disease_ids = [], []
                
            row['hpo_name'] = ';'.join(modes)
            row['disease_id'] = ';'.join(disease_ids)
            return row
        
        inheritance_modes = inheritance_modes.apply(deduplicate_inheritance, axis=1)
        inheritance_modes.columns = ['ncbi_gene_id', 'inheritance_modes', 'inheritance_disease_ids']
        
        logger.info(f"Processed and deduplicated inheritance data for {len(inheritance_modes)} genes")
        return inheritance_modes
    except Exception as e:
        logger.error(f"Error processing inheritance data: {str(e)}")
        raise

def merge_and_format_results(collapsed_df, inheritance_df):
    """
    Merge main annotations with inheritance data and format final results.
    
    Args:
        collapsed_df (pd.DataFrame): Collapsed main annotations
        inheritance_df (pd.DataFrame): Processed inheritance data
    
    Returns:
        pd.DataFrame: Final merged and formatted dataframe
    """
    try:
        final_df = pd.merge(collapsed_df, inheritance_df, on='ncbi_gene_id', how='left')
        final_df['inheritance_modes'] = final_df['inheritance_modes'].fillna('-')
        final_df['inheritance_disease_ids'] = final_df['inheritance_disease_ids'].fillna('-')
        
        logger.info(f"Successfully merged and formatted final results")
        return final_df
    except Exception as e:
        logger.error(f"Error merging and formatting results: {str(e)}")
        raise

def collapse_gene_phenotype_annotations(input_file, output_file):
    """
    Main function to collapse gene-phenotype annotations.
    
    Args:
        input_file (str): Path to input TSV file
        output_file (str): Path to output TSV file
    
    Returns:
        pd.DataFrame: Final collapsed annotations dataframe
    """
    try:
        # Read data
        df = read_phenotype_data(input_file)
        
        # Separate inheritance data
        main_df, inheritance_df = separate_inheritance_data(df)
        
        # Process main annotations
        collapsed_df = collapse_main_annotations(main_df)
        
        # Process inheritance data
        inheritance_modes = process_inheritance_data(inheritance_df)
        
        # Merge and format results
        final_df = merge_and_format_results(collapsed_df, inheritance_modes)
        
        # Save results
        final_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Successfully wrote collapsed annotations to {output_file}")
        
        return final_df
        
    except Exception as e:
        logger.error(f"Error in main processing pipeline: {str(e)}")
        raise


if __name__ == "__main__":
    collapse_gene_phenotype_annotations(sys.argv[1], sys.argv[2])