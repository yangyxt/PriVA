import pandas as pd
import numpy as np
import multiprocessing as mp
import pysam
import argparse as ap
import logging
import sys
import traceback
from typing import List, Dict, Any
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
    logger.setLevel(logging.DEBUG)
    
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
    
    tranx_annos = csq_field
    logger.debug(f"There are {len(tranx_annos)} transcript level annotations in the CSQ field")
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
    logger.debug(f"There are {len(anno_dict)} transcript level annotations in the CSQ field after the extraction")
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
    """Convert a record dictionary to a list of dictionaries and collect logs."""
    record_dict, worker_id, csq_fields, clinvar_csq_fields, var_source_exists, control_ac_exists = args
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
                rows.append(row_dict)
            
        logger.debug(f"Completed processing variant at {record_dict['chrom']}:{record_dict['pos']}, which contains {len(rows)} transcript level annotations\n")
        return rows, list(collector.log_buffer)
        
    except Exception as e:
        logger.error(f"Error processing variant at {record_dict['chrom']}:{record_dict['pos']}: {str(e)}\n\n")
        raise(e)
        return [], list(collector.log_buffer)


def convert_vcf_to_tab(input_vcf: str, threads=4) -> pd.DataFrame:
    all_rows = []
    
    try:
        with pysam.VariantFile(input_vcf) as vcf_file:
            # Extract the VEP CSQ field description from the header
            csq_fields = vcf_file.header.info['CSQ'].description.split('Format: ')[1].split('|')
            clinvar_csq_fields = vcf_file.header.info['CLNCSQ'].description.split('Format: ')[1].split('|')
            var_source_exists = "VARIANT_SOURCE" in vcf_file.header.info
            control_ac_exists = "control_AC" in vcf_file.header.info
            # Convert records to dictionaries before multiprocessing
            record_args = ((extract_record_info(record, var_source_exists, control_ac_exists), 
                          f"{record.chrom}:{record.pos}:{record.ref}->{record.alts[0] if record.alts else ''}",
                          tuple(csq_fields), 
                          tuple(clinvar_csq_fields), 
                          var_source_exists, 
                          control_ac_exists) for record in vcf_file)

            with mp.Pool(threads) as pool:
                varcount = 0
                results = pool.imap_unordered(convert_record_to_tab, record_args)
                for rows, logs in results:
                    if logs:
                        sys.stderr.write("\n".join(logs) + "\n")
                        sys.stderr.flush()
                    varcount += 1
                    all_rows.extend(rows)
            # varcount = 0
            # for tup in record_args:
            #     rows, logs = convert_record_to_tab(tup)
            #     if logs:
            #         sys.stderr.write("\n".join(logs) + "\n")
            #         sys.stderr.flush()
            #     varcount += 1
            #     all_rows.extend(rows)
        
        if all_rows:
            logger.info(f"Completed processing {varcount} variants, which contains {len(all_rows)} transcript level annotations")
            anno_df = pd.DataFrame(all_rows).dropna(subset=["chrom"])
            logger.info(f"The annotation table has {anno_df.shape[0]} rows and {anno_df.shape[1]} columns. And it looks like this: \n{anno_df.head().to_string(index=False)}")
            return anno_df
        return pd.DataFrame()
        
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


def main_combine_annotations(input_vcf: str,
                            cadd_tab: str, 
                            hpo_tab: str, 
                            threads=4,
                            output_tab = None):
    converted_tab = convert_vcf_to_tab(input_vcf, threads)
    cadd_tab = pd.read_table(cadd_tab, low_memory=False)

    cadd_tab["chrom"] = "chr" + cadd_tab["#Chrom"].astype(str)
    cadd_tab["pos"] = cadd_tab["Pos"].astype(int)
    cadd_tab["ref"] = cadd_tab["Ref"].astype(str)
    cadd_tab["alt"] = cadd_tab["Alt"].astype(str)
    cadd_tab["CADD_phred"] = cadd_tab["PHRED"].astype(float)
    cadd_tab["Feature"] = cadd_tab["FeatureID"]

    merged_tab = pd.merge(converted_tab, cadd_tab[["chrom", "pos", "ref", "alt", "Feature", "CADD_phred"]], on=["chrom", "pos", "ref", "alt", "Feature"], how="left")
    cadd_tab["CADD_reg_phred"] = cadd_tab["CADD_phred"]
    merged_tab = pd.merge(merged_tab, cadd_tab.loc[cadd_tab["Feature"].str.contains("ENSR").fillna(False), ["chrom", "pos", "ref", "alt", "CADD_reg_phred"]], on=["chrom", "pos", "ref", "alt"], how="left")

    # Read hpo tab file, which is a tsv.gz file
    hpo_tab = pd.read_table(hpo_tab, low_memory=False)
    hpo_tab["SYMBOL"] = hpo_tab["gene_symbol"].astype(str)
    hpo_tab["HPO_IDs"] = hpo_tab["hpo_id"].astype(str)
    hpo_tab["HPO_terms"] = hpo_tab["hpo_name"].astype(str)
    hpo_tab["HPO_sources"] = hpo_tab["disease_id"].astype(str)
    hpo_tab["HPO_gene_inheritance"] = hpo_tab["inheritance_modes"].astype(str)

    merged_tab = pd.merge(merged_tab, hpo_tab[["SYMBOL", "HPO_IDs", "HPO_terms", "HPO_sources", "HPO_gene_inheritance"]], on="SYMBOL", how="left")

    if output_tab:
        merged_tab.to_csv(output_tab, sep="\t", index=False)
    return merged_tab


if __name__ == "__main__":
    args = ap.ArgumentParser()
    args.add_argument("--input", "-i", type=str, required=True, help="The input VCF file")
    args.add_argument("--output", "-o", type=str, required=False, help="The output tab file", default=None)
    args.add_argument("--threads", "-t", type=int, default=4, help="The number of threads")
    args.add_argument("--cadd", "-c", type=str, required=True, help="The CADD tab file")
    args.add_argument("--hpo", "-p", type=str, required=True, help="The HPO tab file")
    args = args.parse_args()

    main_combine_annotations(args.input, args.cadd, args.hpo, threads = args.threads, output_tab = args.output)