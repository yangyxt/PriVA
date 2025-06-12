import pandas as pd
import numpy as np
import multiprocessing as mp
import pysam
import argparse as ap
import logging
import sys
import traceback
import gc
import pyarrow as pa
import pyarrow.parquet as pq
import tempfile
import subprocess
import time
import os
from typing import List, Dict, Any, Tuple
from collections import deque


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


class SubprocessLogCollector:
    def __init__(self, worker_id: int):
        self.worker_id = worker_id
        self.log_buffer = deque()
        
    def write(self, message: str):
        self.log_buffer.append(message)
        
    def flush(self):
        pass


def setup_worker_logger(worker_id: str) -> logging.Logger:
    logger = logging.getLogger(f'worker_{worker_id}')
    
    # Only setup if not already configured (avoid duplicate handlers)
    if not logger.handlers:
        logger.setLevel(logging.INFO)
        
        # Create a collector for this worker
        collector = SubprocessLogCollector(worker_id)
        handler = logging.StreamHandler(collector)
        formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        
        # Store collector as logger attribute for reuse
        logger._collector = collector
    
    return logger, logger._collector



def na_value(value):
    if isinstance(value, float):
        if np.isnan(value):
            return True
        else:
            return False
    elif isinstance(value, int):
        return False
    elif value is None:
        return True
    elif isinstance(value, str):
        return value in ["NaN", "nan", "na", "NA", "NAN", "None", "none", "", ".", "-"]
    else:
        return False


def parse_csq_field(csq_field: Tuple[str, ...], csq_fields: List[str], logger: logging.Logger) -> Dict[str, Any]:
    '''
    The function is used to parse CSQ field into a dictionary
    {value_name: [value_list],
        ...
    }
    one csq field can contain multiple variant-transcript consequences
    '''

    if not csq_field:
        return None
    
    tranx_annos = csq_field
    logger.debug(f"There are {len(tranx_annos)} transcript level annotations in the CSQ field")
    if len(tranx_annos) == 0:
        logger.warning(f"There are no transcript level annotations in the CSQ field: {csq_field}")
    anno_dict = {}
    for tranx_anno in tranx_annos:
        field_values = tranx_anno.split("|")
        feature_dict = {}
        feature_name = field_values[csq_fields.index("Feature")]
        for i in range(len(csq_fields)):
            field_value = field_values[i]
            # See if the field value is numerical thus can be converted to int or float
            if field_value.isdigit():
                try:
                    field_value = int(field_value)
                except ValueError:
                    try:
                        field_value = float(field_value)
                    except ValueError:
                        pass
            field_value = np.nan if field_value in ["", ".", None, np.nan, "-", "NaN", "nan", "NA", "N/A"] else field_value
            feature_dict[csq_fields[i]] = field_value
        if feature_name:
            anno_dict[feature_name] = feature_dict
            logger.debug(f"Adding feature dict {feature_dict} to the anno_dict with key {feature_name}, now anno_dict contain {len(anno_dict)} transcript level annotations")
    return anno_dict


def extract_record_info(record, var_source_exists: bool, control_ac_exists: bool):
    """Extract necessary information from a VariantRecord object"""
    # Initialize the info dictionary
    info_dict = {
        # Core frequency fields
        'AC_joint': record.info.get('AC_joint', [np.nan])[0],
        'AN_joint': record.info.get('AN_joint', np.nan),
        'AF_joint': record.info.get('AF_joint', [np.nan])[0],
        'nhomalt_joint': record.info.get('nhomalt_joint', [np.nan])[0],
        
        # Maximum values across populations
        'AC_grpmax_joint': record.info.get('AC_grpmax_joint', [np.nan])[0],
        'AF_grpmax_joint': record.info.get('AF_grpmax_joint', [np.nan])[0],
        'AN_grpmax_joint': record.info.get('AN_grpmax_joint', [np.nan])[0],
        'nhomalt_grpmax_joint': record.info.get('nhomalt_grpmax_joint', [np.nan])[0],
        
        # Sex-specific fields (overall)
        'AF_joint_XX': record.info.get('AF_joint_XX', [np.nan])[0],
        'AF_joint_XY': record.info.get('AF_joint_XY', [np.nan])[0],
        'nhomalt_joint_XX': record.info.get('nhomalt_joint_XX', [np.nan])[0],
        'nhomalt_joint_XY': record.info.get('nhomalt_joint_XY', [np.nan])[0],
        
        # ClinVar information
        'CLNDN': ",".join(record.info.get('CLNDN', [""])),
        'CLNSIG': ",".join(record.info.get('CLNSIG', [""])),
        'CLNREVSTAT': ",".join(record.info.get('CLNREVSTAT', [""])),

        # CSQ field
        'CSQ': record.info.get('CSQ', tuple(["",])),
    }

    if var_source_exists:
        info_dict['VARIANT_SOURCE'] = record.info.get('VARIANT_SOURCE', "")
    
    if control_ac_exists:
        info_dict['control_AC'] = record.info.get('control_AC', (np.nan,))[0]
        info_dict['control_AN'] = record.info.get('control_AN', np.nan)
        info_dict['control_AF'] = record.info.get('control_AF', (np.nan,))[0]
        info_dict['control_nhomalt'] = record.info.get('control_nhomalt', np.nan)
    
    # Add population-specific fields
    pop_codes = ["nfe", "eas", "afr", "amr", "asj", "fin", "sas", "mid", "remaining"]
    
    for pop in pop_codes:
        # Add XX/XY-specific allele frequencies for each population
        info_dict[f'AF_joint_{pop}_XX'] = record.info.get(f'AF_joint_{pop}_XX', [np.nan])[0]
        info_dict[f'AF_joint_{pop}_XY'] = record.info.get(f'AF_joint_{pop}_XY', [np.nan])[0]
        
        # Add XX homozygote counts (as specified in the script)
        info_dict[f'nhomalt_joint_{pop}_XX'] = record.info.get(f'nhomalt_joint_{pop}_XX', [np.nan])[0]
    
    return {
        'chrom': record.chrom,
        'pos': record.pos,
        'ref': record.ref,
        'alts': record.alts,
        'VCF_filters': ",".join(record.filter),
        'info': info_dict,
        
        # Extract GT info for each sample, preserving phasing information
        'GT': {
            sample: '|'.join('.' if allele is None else str(allele) for allele in record.samples[sample]['GT']) 
            if record.samples[sample].phased 
            else '/'.join('.' if allele is None else str(allele) for allele in record.samples[sample]['GT'])
            for sample in record.samples.keys()
        },
        
        # Extract AD info for each sample
        'AD': {
            f"{sample}_AD": ','.join('.' if value is None else str(value) for value in record.samples[sample]['AD'])
            for sample in record.samples.keys()
        }
    }


def convert_record_to_tab(args: tuple) -> tuple[List[Dict[str, Any]], List[str]]:
    start_time = time.time()
    """Convert a record dictionary to a list of dictionaries and collect logs."""
    record_dict, worker_id, csq_fields, clinvar_csq_fields, var_source_exists, control_ac_exists, cadd_phred_dict, hpo_symbol_map = args
    logger, collector = setup_worker_logger(worker_id)
    
    try:
        rows = []
        logger.debug(f"Processing variant at {record_dict['chrom']}:{record_dict['pos']}\n")

        assert isinstance(record_dict['chrom'], str), f"The chromosome is not a valid value: {record_dict}"
        
        # Extract the variant-specific annotations
        var_dict_items = {
            "chrom": record_dict['chrom'],
            "pos": record_dict['pos'],
            "ref": record_dict['ref'],
            "alt": record_dict['alts'][0] if record_dict['alts'] else '',
            
            # Core frequency fields
            "gnomAD_joint_AN": record_dict['info']['AN_joint'],
            "gnomAD_joint_AF": record_dict['info']['AF_joint'],
            
            # Maximum values across populations
            "gnomAD_joint_AF_max": record_dict['info']['AF_grpmax_joint'],
            "gnomAD_joint_AN_max": record_dict['info']['AN_grpmax_joint'],
            "gnomAD_nhomalt_max": record_dict['info']['nhomalt_grpmax_joint'],
            
            # Sex-specific fields (overall)
            "gnomAD_joint_AF_XX": record_dict['info']['AF_joint_XX'],
            "gnomAD_joint_AF_XY": record_dict['info']['AF_joint_XY'],
            "gnomAD_nhomalt_XX": record_dict['info']['nhomalt_joint_XX'],
            "gnomAD_nhomalt_XY": record_dict['info']['nhomalt_joint_XY'],
            
            # ClinVar information
            "CLNSIG": record_dict['info']['CLNSIG'],
            "CLNREVSTAT": record_dict['info']['CLNREVSTAT'],
            "VCF_filters": record_dict['VCF_filters'],
        }

        if var_source_exists:
            var_dict_items["VARIANT_SOURCE"] = record_dict['info']['VARIANT_SOURCE']
        
        if control_ac_exists:
            var_dict_items["control_AC"] = record_dict['info']['control_AC']
            var_dict_items["control_AN"] = record_dict['info']['control_AN']
            var_dict_items["control_AF"] = record_dict['info']['control_AF']
            var_dict_items["control_nhomalt"] = record_dict['info']['control_nhomalt']
        
        # # Add population-specific fields
        # pop_codes = ["nfe", "eas", "afr", "amr", "asj", "fin", "sas", "mid", "remaining"]
        
        # for pop in pop_codes:
        #     # Add XX/XY-specific allele frequencies for each population
        #     var_dict_items[f"gnomAD_AF_{pop}_XX"] = record_dict['info'].get(f'AF_joint_{pop}_XX', np.nan)
        #     var_dict_items[f"gnomAD_AF_{pop}_XY"] = record_dict['info'].get(f'AF_joint_{pop}_XY', np.nan)
            
        #     # Add XX homozygote counts (as specified in the script)
        #     var_dict_items[f"gnomAD_nhomalt_{pop}_XX"] = record_dict['info'].get(f'nhomalt_joint_{pop}_XX', np.nan)
        
        gt_dict = record_dict['GT']
        ad_dict = record_dict['AD']
        
        # Extract the variant-transcript level annotations
        csqs = parse_csq_field(record_dict['info']['CSQ'], csq_fields, logger)

        for feature_name, feature_dict in csqs.items():
            if feature_name.startswith("ENS"):
                row_dict = {**var_dict_items, **feature_dict, **gt_dict, **ad_dict}
                cadd_phred = cadd_phred_dict.get(feature_name, {"CADD_phred": np.nan, "CADD_reg_phred": np.nan})
                if len(cadd_phred) == 1:
                    if feature_name.startswith("ENST"):
                        cadd_phred.update({"CADD_reg_phred": np.nan})
                        row_dict.update(cadd_phred)
                    else:
                        row_dict.update({"CADD_reg_phred": cadd_phred["CADD_phred"], "CADD_phred": np.nan})
                else:
                    row_dict.update(cadd_phred)

                if len(hpo_symbol_map) == 0:
                    hpo_symbol_map = {"HPO_IDs": np.nan, "HPO_terms": np.nan, "HPO_sources": np.nan, "HPO_gene_inheritance": np.nan}
                row_dict.update(hpo_symbol_map)
                rows.append(row_dict)
        
        if rows:
            logger.debug(f"Completed processing variant at {record_dict['chrom']}:{record_dict['pos']}, which contains {len(rows)} transcript level annotations\n")
        else:
            logger.debug(f"Variant at {record_dict['chrom']}:{record_dict['pos']}:{record_dict['ref']}->{record_dict['alts'][0] if record_dict['alts'] else ''} seems to be an intergenic variant")
        
        return rows, list(collector.log_buffer), time.time() - start_time
        
    except Exception as e:
        logger.error(f"Error processing variant at {record_dict['chrom']}:{record_dict['pos']}: {str(e)}\n\n")
        raise(e)
        return [], list(collector.log_buffer), time.time() - start_time


def arg_generator(vcf_file, threads, cadd_phred_dict, hpo_symbol_map):
    with pysam.VariantFile(vcf_file) as vcf_file:
        # Extract the VEP CSQ field description from the header
        csq_fields = vcf_file.header.info['CSQ'].description.split('Format: ')[1].split('|')
        symbol_field_index = csq_fields.index("SYMBOL")
        clinvar_csq_fields = vcf_file.header.info['CLNCSQ'].description.split('Format: ')[1].split('|')
        var_source_exists = "VARIANT_SOURCE" in vcf_file.header.info
        control_ac_exists = "control_AC" in vcf_file.header.info
        # Convert records to dictionaries before multiprocessing
        worker_ind = 0
        for record in vcf_file:
            worker_ind += 1
            yield ( extract_record_info(record, var_source_exists, control_ac_exists), 
                    str(worker_ind % threads),
                    tuple(csq_fields), 
                    tuple(clinvar_csq_fields), 
                    var_source_exists, 
                    control_ac_exists,
                    cadd_phred_dict.get(f"{record.chrom}:{record.pos}:{record.ref}-{record.alts[0] if record.alts else ''}", {}),
                    {k: v for k,v in hpo_symbol_map.items() if k in [ csq.split("|")[symbol_field_index] for csq in record.info.get("CSQ", ("|" * symbol_field_index,)) ]} )


def convert_vcf_to_tab(input_vcf: str, threads=4, cadd_phred_dict: dict = None, hpo_symbol_map: dict = None, chunk_size: int = 20000) -> pd.DataFrame:
    est_anno_count = subprocess.run(f'bcftools query -f "%INFO/CSQ\\n" {input_vcf} | tr \',\' \'\\n\' | wc -l', shell=True, capture_output=True, text=True).stdout.strip()
    est_anno_count = int(est_anno_count)
    logger.info(f"The VCF file {input_vcf} has {est_anno_count} transcript_level annotations")
    
    try:
        record_args = arg_generator(input_vcf, threads, cadd_phred_dict, hpo_symbol_map)
        varcount = 0
        anno_record_count = 0
        batch_run_time = np.zeros(5000, dtype=np.float16)
        batch_extend_time = np.zeros(5000, dtype=np.float16)
        all_rows = []
        col_dicts = None

        if threads == 1:
            results = map(convert_record_to_tab, record_args)
            pool = None
        else:
            pool = mp.Pool(threads)
            results = pool.imap_unordered(convert_record_to_tab, record_args)

        for rows, logs, run_time in results:
            if logs:
                sys.stderr.write("\n".join(logs) + "\n")
                sys.stderr.flush()
            varcount += 1
            before_ext_time = time.time()
            all_rows.extend(rows)
            after_ext_time = time.time()
            batch_run_time[varcount % 5000 - 1] = run_time
            batch_extend_time[varcount % 5000 - 1] = after_ext_time - before_ext_time
            if varcount % 5000 == 0:
                logger.info(f"At Variant {varcount}: total run time across past 5000 variants: {batch_run_time.sum():.6f}s, extend time: {batch_extend_time.sum():.6f}s")
            if col_dicts is None: col_dicts = {key: np.empty(est_anno_count, dtype=object) for key in rows[0].keys()}
            
            # Force garbage collection every 1000 variants to prevent memory buildup
            if varcount % 1000 == 0:
                gc.collect()
            if len(all_rows) >= chunk_size:
                # Append dict to the col dicts
                start_collect = time.time()
                num_rows = len(all_rows)
                for key in all_rows[0].keys():
                    has_digit = False
                    values = [None] * num_rows
                    for i in range(num_rows):
                        d = all_rows[i]
                        value = np.nan if na_value(d.get(key, np.nan)) else d.get(key, np.nan)
                        try:
                            value = float(value)
                        except ValueError:
                            pass
                        else:
                            if not na_value(value):
                                has_digit = True
                        values[i] = value
                    
                    if has_digit and col_dicts[key].dtype != np.float32:
                        logger.info(f"Found the column {key} has digit values in the iterating chunk, convert the np array type to float32")
                        col_dicts[key] = col_dicts[key].astype(np.float32)
                    
                    col_dicts[key][anno_record_count:anno_record_count + num_rows] = values
                anno_record_count += num_rows
                finish_collect = time.time()
                logger.info(f"Completed processing {varcount} variants, which contains {anno_record_count} transcript level annotations, collect time for this chunk: {finish_collect - start_collect:.6f}s")
                all_rows = []
                gc.collect()  # Force garbage collection after processing each chunk
        
        if all_rows:
            num_rows = len(all_rows)
            for key in all_rows[0].keys():
                has_digit = False
                values = [None] * num_rows
                for i in range(num_rows):
                    d = all_rows[i]
                    value = np.nan if na_value(d.get(key, np.nan)) else d.get(key, np.nan)
                    try:
                        value = float(value)
                    except ValueError:
                        pass
                    else:
                        if not na_value(value):
                            has_digit = True
                    values[i] = value
                
                if has_digit and col_dicts[key].dtype != np.float32:
                    logger.info(f"Found the column {key} has digit values in the iterating chunk, convert the np array type to float32")
                    col_dicts[key] = col_dicts[key].astype(np.float32)
                
                col_dicts[key][anno_record_count:anno_record_count + num_rows] = values
            anno_record_count += num_rows
            logger.info(f"Completed processing {varcount} variants, which contains {anno_record_count} transcript level annotations")
            gc.collect()  # Force garbage collection after final batch
            
        if pool:
            pool.close()
            pool.join()
            pool.terminate()
            pool = None
        
        col_dicts = {k: v[:anno_record_count] for k, v in col_dicts.items()}
        df = pd.DataFrame(col_dicts)
        for col in df.columns:
            # Convert float back to int if possible
            df.loc[:, col] = df.loc[:, col].replace("NaN", np.nan)
            try:
                df.loc[:, col] = df.loc[:, col].astype(np.float32)
            except ValueError:
                pass
            else:
                try:
                    df.loc[:, col] = df.loc[:, col].astype(np.int32)
                except ValueError:
                    pass
        
        logger.info(f"The annotation table has {df.shape[0]} rows and {df.shape[1]} columns. And it looks like this: \n{df.head().to_string(index=False)}")
        return df
        
    except Exception as e:
        # Capture full traceback
        error_traceback = traceback.format_exc()
        sys.stderr.write(f"Error in conversion process: {str(e)}\n")
        sys.stderr.write(f"Full traceback:\n{error_traceback}\n")
        
        # Optionally log to file as well
        logger.error(f"Error in conversion process: {str(e)}")
        logger.error(f"Full traceback:\n{error_traceback}")
        
        # Reraise the exception
        raise(e)
    finally:
        if pool:
            try:
                pool.close()
                pool.join()
            except:
                pass


def main_combine_annotations(input_vcf: str,
                            cadd_tab: str, 
                            hpo_tab: str, 
                            threads=4,
                            output_tab = None):
    
    # Read CADD data with memory optimization - only load required columns
    logger.info("Loading CADD data...")
    cadd_required_columns = ["#Chrom", "Pos", "Ref", "Alt", "PHRED", "FeatureID"]
    cadd_dtypes = {
        "#Chrom": "string",  # Chromosomes are truly categorical (chr1, chr2, etc.)
        "Pos": "int32",        # Positions don't need int64
        "Ref": "string",       # Ref bases - use string, not category (variable length)
        "Alt": "string",       # Alt bases - use string, not category (variable length)
        "PHRED": "float32",    # CADD scores don't need float64 precision
        "FeatureID": "string"  # Feature IDs - use string, not category (variable length)
    }
    cadd_df = pd.read_table(cadd_tab, low_memory=False, usecols=cadd_required_columns, dtype=cadd_dtypes)
    logger.info(f"CADD data loaded with {cadd_df.shape[0]} rows and {cadd_df.shape[1]} columns (only required columns)")

    # Prepare CADD data more efficiently - data types already optimized during reading
    cadd_df["chrom"] = "chr" + cadd_df["#Chrom"].astype(str)
    cadd_df["pos"] = cadd_df["Pos"]  # Already int32
    cadd_df["ref"] = cadd_df["Ref"]  # Already string dtype, no conversion needed
    cadd_df["alt"] = cadd_df["Alt"]  # Already string dtype, no conversion needed
    cadd_df['variant_id'] = cadd_df["chrom"] + ":" + cadd_df["pos"].astype(str) + ":" + cadd_df["ref"] + "-" + cadd_df["alt"]
    cadd_df["CADD_phred"] = cadd_df["PHRED"]  # Already float32
    cadd_df["Feature"] = cadd_df["FeatureID"]  # Already string dtype

    # Create subsets to reduce memory during merges
    cadd_subset = cadd_df[["variant_id", "Feature", "CADD_phred"]].copy()
    logger.info(f"CADD subset memory usage: {cadd_subset.memory_usage(deep=True).sum() / 1024**2:.1f} MB")

    cadd_phred_dict = cadd_subset.groupby("variant_id").apply(lambda x: x.set_index("Feature")[["CADD_phred"]].to_dict(orient="index")).to_dict()

    # Read HPO data - only load required columns
    logger.info("Loading HPO data...")
    hpo_required_columns = ["gene_symbol", "hpo_id", "hpo_name", "disease_id", "inheritance_modes"]
    hpo_dtypes = {
        "gene_symbol": "string",        # Gene symbols - use string (some are long/variable)
        "hpo_id": "string",            # HPO IDs as string
        "hpo_name": "string",          # HPO names as string  
        "disease_id": "string",        # Disease IDs as string
        "inheritance_modes": "string"   # Inheritance modes as string (safer than category)
    }
    hpo_df = pd.read_table(hpo_tab, low_memory=False, usecols=hpo_required_columns, dtype=hpo_dtypes)
    logger.info(f"HPO data loaded with {hpo_df.shape[0]} rows and {hpo_df.shape[1]} columns (only required columns)")

    # Prepare HPO data efficiently - data types already optimized during reading
    hpo_df.rename(columns={
        "gene_symbol": "SYMBOL",
        "hpo_id": "HPO_IDs", 
        "hpo_name": "HPO_terms",
        "disease_id": "HPO_sources",
        "inheritance_modes": "HPO_gene_inheritance"
    }, inplace=True)
    hpo_df = hpo_df.loc[hpo_df["SYMBOL"] != "-", :]
    hpo_symbol_map = hpo_df.set_index("SYMBOL").to_dict()
    
    # Convert VCF to dataframe
    converted_tab = convert_vcf_to_tab(input_vcf, threads, cadd_phred_dict, hpo_symbol_map)
    logger.info(f"Final converted table memory usage: {converted_tab.memory_usage(deep=True).sum() / 1024**2:.1f} MB")
    
    logger.info(f"Final table memory usage: {converted_tab.memory_usage(deep=True).sum() / 1024**2:.1f} MB")
    logger.info(f"Final table shape: {converted_tab.shape}")

    if output_tab:
        logger.info(f"Writing output to {output_tab}...")
        converted_tab.to_csv(output_tab, sep="\t", index=False)
    
    return converted_tab


if __name__ == "__main__":
    args = ap.ArgumentParser()
    args.add_argument("--input", "-i", type=str, required=True, help="The input VCF file")
    args.add_argument("--output", "-o", type=str, required=False, help="The output tab file", default=None)
    args.add_argument("--threads", "-t", type=int, default=4, help="The number of threads")
    args.add_argument("--cadd", "-c", type=str, required=True, help="The CADD tab file")
    args.add_argument("--hpo", "-p", type=str, required=True, help="The HPO tab file")
    args = args.parse_args()

    main_combine_annotations(args.input, args.cadd, args.hpo, threads = args.threads, output_tab = args.output)