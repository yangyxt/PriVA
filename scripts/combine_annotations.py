import pandas as pd
import numpy as np
import multiprocessing as mp
import pysam
import argparse as ap
import json
import logging
import sys
from typing import List, Dict, Any
from collections import deque


class SubprocessLogCollector:
    def __init__(self, worker_id: int):
        self.worker_id = worker_id
        self.log_buffer = deque()
        
    def write(self, message: str):
        self.log_buffer.append(message)
        
    def flush(self):
        pass


def setup_worker_logger(worker_id: int) -> logging.Logger:
    logger = logging.getLogger(f'worker_{worker_id}')
    logger.setLevel(logging.INFO)
    
    # Create a collector for this worker
    collector = SubprocessLogCollector(worker_id)
    handler = logging.StreamHandler(collector)
    formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger, collector



def parse_csq_field(csq_field: str, csq_fields: List[str], logger: logging.Logger) -> Dict[str, Any]:
    '''
    The function is used to parse CSQ field into a dictionary
    {value_name: [value_list],
        ...
    }
    one csq field can contain multiple variant-transcript consequences
    '''

    if not csq_field:
        return None
    
    tranx_annos = csq_field.split(",")
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
            field_value = np.nan if field_value == "" else field_value
            feature_dict[csq_fields[i]] = field_value
        anno_dict[feature_name] = feature_dict

    return anno_dict




def convert_record_to_tab(args: tuple) -> tuple[List[Dict[str, Any]], List[str]]:
    """Convert a VCF record to a list of dictionaries and collect logs.
    
    Returns:
        Tuple of (list of row dictionaries, list of log messages)
    """
    record, worker_id, csq_fields, clinvar_csq_fields = args
    logger, collector = setup_worker_logger(worker_id)
    
    try:
        rows = [] # List of dictionaries, very efficiently pickable
        # Process record and create rows...
        logger.info(f"Processing variant at {record.chrom}:{record.pos}\n")
        
        # Extract the variant-specific annotations
        var_dict_items = {}
        var_dict_items["chrom"] = record.chrom
        var_dict_items["pos"] = record.pos
        var_dict_items["ref"] = record.ref
        var_dict_items["alt"] = record.alts[0]
        var_dict_items["gnomAD_joint_af"] = record.info.get('AF_joint', [np.nan])[0]
        var_dict_items["gnomAD_nhomalt_XX"] = record.info.get('nhomalt_joint_XX', [np.nan])[0]
        var_dict_items["gnomAD_nhomalt_XY"] = record.info.get('nhomalt_joint_XY', [np.nan])[0]
        var_dict_items["CLNSIG"] = record.info.get('CLNSIG', [""])[0]
        var_dict_items["CLNREVSTAT"] = record.info.get('CLNREVSTAT', [""])[0]        
        
        # Extract the variant-transcript level annotations
        csqs = parse_csq_field(record.info.get('CSQ', ["",])[0], csq_fields, logger)

        for feature_name, feature_dict in csqs.items():
            if feature_name.startswith("ENS"):
                row_dict = {**var_dict_items, **feature_dict}
                rows.append(row_dict)
            
        logger.info(f"Completed processing variant at {record.chrom}:{record.pos}\n")
        return rows, list(collector.log_buffer)
        
    except Exception as e:
        logger.error(f"Error processing variant at {record.chrom}:{record.pos}: {str(e)}\n\n")
        return [], list(collector.log_buffer)


def convert_vcf_to_tab(input_vcf: str, threads=4) -> pd.DataFrame:
    all_rows = []
    
    try:
        with pysam.VariantFile(input_vcf) as vcf_file:
            # Extract the VEP CSQ (in INFO) field description from the header
            csq_fields = vcf_file.header.info['CSQ'].description.split('Format: ')[1].split('|')
            clinvar_csq_fields = vcf_file.header.info['CLNCSQ'].description.split('Format: ')[1].split('|')
            record_args = ((record, 
                            f"{record.chrom}:{record.pos}:{record.ref}->{record.alts[0]}", 
                            tuple(csq_fields), 
                            tuple(clinvar_csq_fields)) for record in vcf_file)

            with mp.Pool(threads) as pool:
                for rows, logs in pool.imap_unordered(convert_record_to_tab, record_args):
                    # Print collected logs as a group
                    if logs:
                        sys.stderr.write("\n".join(logs) + "\n")
                        sys.stderr.flush()
                    all_rows.extend(rows)
        
        if all_rows:
            return pd.DataFrame(all_rows)
        return pd.DataFrame()
        
    except Exception as e:
        sys.stderr.write(f"Error in conversion process: {str(e)}\n")
        raise


def main(input_vcf: str, output_tab: str, cadd_tab: str,threads=4):
    converted_tab = convert_vcf_to_tab(input_vcf, threads)
    cadd_tab = pd.read_table(cadd_tab, low_memory=False)

    cadd_tab["chrom"] = "chr" + cadd_tab["#Chrom"].astype(str)
    cadd_tab["pos"] = cadd_tab["Pos"].astype(int)
    cadd_tab["ref"] = cadd_tab["Ref"].astype(str)
    cadd_tab["alt"] = cadd_tab["Alt"].astype(str)
    cadd_tab["CADD_phred"] = cadd_tab["PHRED"].astype(float)

    merged_tab = pd.merge(converted_tab, cadd_tab[["chrom", "pos", "ref", "alt", "CADD_phred"]], on=["chrom", "pos", "ref", "alt"], how="left")
    merged_tab.to_csv(output_tab, sep="\t", index=False)
    return merged_tab


if __name__ == "__main__":
    args = ap.ArgumentParser()
    args.add_argument("--input", "-i", type=str, required=True, help="The input VCF file")
    args.add_argument("--output", "-o", type=str, required=True, help="The output tab file")
    args.add_argument("--threads", "-t", type=int, default=4, help="The number of threads")
    args.add_argument("--cadd", "-c", type=str, required=True, help="The CADD tab file")
    args = args.parse_args()

    main(args.input, args.output, args.cadd, args.threads)