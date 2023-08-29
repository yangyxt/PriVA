# from ctypes.wintypes import tagSIZE
from tokenize import group
import pandas as pd
import os
import sys
import re
import numpy as np
import argparse as ap
import logging
import gc
import datetime
import subprocess


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def convert_per_row(row, table, label, delimiter):
    # Select the cell where it matches with the Gene symbol in to_be_modified_table
    selected_record = table.loc[table['symbol'].isin(str(row[label]).split(delimiter)), 'symbol']
    # Convert the series into list
    # If we found one match
    if len(selected_record.tolist()) == 1:
        # logger.info("This record (index {}) contains an updated HGNC symbol {}".format(row.name, str(row[label])))
        return selected_record.tolist()[0]
    # If we found more than one match
    elif len(selected_record.tolist()) > 1:
        logger.debug("This record (index {}) contains multiple updated HGNC symbols {}".format(row.name, str(row[label])))
        return delimiter.join(selected_record.tolist())
    # If we haven't found one match, it means the symbol in to_be_modified table is an alias.
    else:
        logger.debug("This record (index {}) does not contain an updated HGNC symbol {}".format(row.name, str(row[label])))
        # Prepare a series of booleans based on checking if there are intersection
        # bool_ser_alias = table['alias_symbol'].apply(intersection, args=str(row[label]).split(delimiter))
        bool_ser_alias = table['alias_symbol'].str.strip('"').str.split("|", expand=True).isin(str(row[label]).split(delimiter)).any(axis=1)
        # bool_ser_prev = table['prev_symbol'].apply(intersection, args=str(row[label]).split(delimiter))
        bool_ser_prev = table['prev_symbol'].str.strip('"').str.split("|", expand=True).isin(str(row[label]).split(delimiter)).any(axis=1)
        alias_or_prev_record_list = table.loc[bool_ser_alias | bool_ser_prev, 'symbol'].tolist()
        if len(alias_or_prev_record_list) == 1:
            return alias_or_prev_record_list[0]
        elif len(alias_or_prev_record_list) > 1:
            # We may found the name match with one prev symbol for one gene and one alias symbol for another gene
            # We need to determine the true gene for the current symbol
            search_symbols = str(row[label]).split(delimiter)
            if bool_ser_alias.any(): alias_symbols = table.loc[bool_ser_alias, 'alias_symbol'].str.split("|", expand=True)
            if bool_ser_prev.any(): prev_symbols = table.loc[bool_ser_prev, 'prev_symbol'].str.split("|", expand=True)
            if 'alias_symbols' in locals() and 'prev_symbols' in locals(): 
                symbol_dict = {alias_symbols.shape[1]: 'alias_symbol', prev_symbols.shape[1]:'prev_symbol'}
                if symbol_dict[min(alias_symbols.shape[1], prev_symbols.shape[1])] == "alias_symbol":
                    return delimiter.join(table.loc[bool_ser_alias,'symbol'].tolist())
                elif symbol_dict[min(alias_symbols.shape[1], prev_symbols.shape[1])] == "prev_symbol":
                    return delimiter.join(table.loc[bool_ser_prev, 'symbol'].tolist())
            else:
                return delimiter.join(alias_or_prev_record_list)
        else:
            return str(row[label])              
                

def retrieve_latest_hgnc_tab(url="http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt", 
                             local_path="/data/hgnc/non_alt_loci_set.tsv",
                             type_col="locus_type",
                             loc_type="",
                             group_col="locus_group",
                             loc_group="",
                             core=True):
    '''
    All possible locus types:
        complex locus constituent
        endogenous retrovirus
        fragile site
        gene with protein product
        immunoglobulin gene
        immunoglobulin pseudogene
        protocadherin
        pseudogene
        readthrough
        region
        RNA, cluster
        RNA, long non-coding
        RNA, micro
        RNA, misc
        RNA, ribosomal
        RNA, small nuclear
        RNA, small nucleolar
        RNA, transfer
        RNA, vault
        RNA, Y
        T cell receptor gene
        T cell receptor pseudogene
        unknown
        virus integration site
        
    All possible locus groups:
        "non-coding RNA",
        "other",
        "protein-coding gene",
        "pseudogene"
    '''
    today = datetime.datetime.today()
    mtime = datetime.datetime.fromtimestamp(os.path.getmtime(local_path))
    duration = today - mtime
    
    if duration.days > 3:
        try:
            hgnc_anno_table = pd.read_csv(url, sep='\t', low_memory=False)
            hgnc_anno_table.to_csv(local_path + ".tmp", sep="\t", index=False)
            result = subprocess.run("mv {} {} && ls -lh {}".format(local_path + ".tmp", local_path, local_path), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
            print(result.stdout, file=sys.stderr)
        except:
            logger.warning("The Request to the URL of latest HGNC symbol table does not have a valid return code. Abandon accessing online form. Use local form instead.")
            hgnc_anno_table = pd.read_csv(local_path, sep='\t', low_memory=False)
    else:
        logger.info("Local file {} is updated enough. We dont need to access to internet to get the latest one.".format(local_path))
        hgnc_anno_table = pd.read_csv(local_path, sep='\t', low_memory=False)
        
    logger.info("Before filtering, the hgnc_anno_table has shape of {} and it looks like:\n{}".format(hgnc_anno_table.shape, hgnc_anno_table[:5].to_string(index=False)))
    
    loc_types = hgnc_anno_table[type_col].drop_duplicates().dropna().tolist()
    loc_groups = hgnc_anno_table[group_col].drop_duplicates().dropna().tolist()
    
    total_filter = np.array([ True for i in range(0, len(hgnc_anno_table)) ])
    for t, ts, tcol in [ (loc_type, loc_types, type_col), (loc_group, loc_groups, group_col) ]:
        if len(t) > 0:
            if t in ts:
                filter_arr = hgnc_anno_table[tcol].astype(str) == t
            else:
                logger.error("Input loc_type: {} is not included in this table, Check the available options: {}".format(loc_type, loc_types))
                assert t in ts
        else:
            filter_arr = np.array([ True for i in range(0, len(hgnc_anno_table)) ])
        total_filter = np.logical_and(total_filter, filter_arr)

    if core:
        return hgnc_anno_table.loc[total_filter, ["symbol", "alias_symbol", "prev_symbol", type_col, group_col]].drop_duplicates()
    else:
        return hgnc_anno_table.loc[total_filter, :].drop_duplicates()

    
    
def update_hgnc(table, d=";", label=None, parallel="False", threads=8, 
                url="http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt",
                local_path="/data/hgnc/non_alt_loci_set.tsv"):
    # a stands for the absolute path to the annotation table, if not inputted, we'll use a default URL
    # b stands for the absolute path to the to_be_modified table, or it can be a file containing a list of genes with one gene per row. b can also be a DataFrame Object
    # c stands for the column label in the to_be_modified table where we compare gene symbols. If b is a gene list file, do not input c, let it be None
    # d stands for the delimiter used in gene symbol field

    hgnc_anno_table = retrieve_latest_hgnc_tab(url, local_path)
    
    if label:
        if type(table) == pd.core.frame.DataFrame:
            to_be_modified_table = table
        elif type(table) == str:
            to_be_modified_table = pd.read_csv(table, sep='\t', low_memory=False)
    else:
        if type(table) == pd.core.frame.DataFrame:
            to_be_modified_table = table
        else:
            to_be_modified_table = pd.read_csv(table, sep='\t', low_memory=False)
        gene_symbol_regex = re.compile("^[A-Z]{1}[A-Z0-9]{2,}$|^C[0-9]+orf[0-9]+$")
        label = [ label for label in to_be_modified_table.columns.tolist() if re.search(gene_symbol_regex, str(to_be_modified_table.at[1,label]))][0]

    new_table = pd.DataFrame()
    new_table['original_symbols'] = list(dict.fromkeys(to_be_modified_table[label]))

    if parallel == "True" or parallel == True:
        import multiprocessing as mp
        pool = mp.Pool(threads)
        new_table['HGNC_symbol'] = pool.starmap(convert_per_row, [(new_table.iloc[i, :], hgnc_anno_table, 'original_symbols', d) for i in range(0, len(new_table))])
    else:
        new_table['HGNC_symbol'] = new_table.apply(convert_per_row, args=(hgnc_anno_table, 'original_symbols', d), axis=1)
        
    new_table.set_index('original_symbols', inplace=True)
    to_be_modified_table = to_be_modified_table.merge(new_table, left_on=label, right_index=True)
    new_table = None
    del new_table
    gc.collect()
    
    logger.warning("Now we have a new column recording the updated HGNC_symbol in: \n" + str(to_be_modified_table[:4].to_string(index=False)))
    original_column_index = to_be_modified_table.columns.get_loc(label)
    to_be_modified_table.drop(columns=label, inplace=True)
    tobe_mov = to_be_modified_table.pop('HGNC_symbol')
    to_be_modified_table.insert(loc=original_column_index, column=label, value=tobe_mov)
    if label == "Later_defined_header":
        to_be_modified_table.to_csv(table, sep='\t', index=False, header=False)
    elif type(table) != pd.core.frame.DataFrame:
        to_be_modified_table.to_csv(table, sep='\t', index=False)
    else:
        return to_be_modified_table


if __name__ == '__main__':
    # a stands for the absolute path to the annotation table
    # b stands for the absolute path to the to_be_modified table
    # c stands for the column label in the to_be_modified table where we compare gene symbols
    # d stands for the delimiter used in

    parser = ap.ArgumentParser()
    parser.add_argument("-ap", "--anno_path", type=str, help="Input an URL link or a path directing to the latest HGNC table", required=False)
    parser.add_argument("-u", "--url", type=str, help="The URL of the HGNC annotation table", required=False, default="http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/non_alt_loci_set.txt")
    parser.add_argument("-v", "--variant_table", type=str, help="The path of variant table where you want to update the gene symbol", required=True)
    parser.add_argument("-ch", "--column_header", type=str, help="The gene symbol column header", required=False)
    parser.add_argument("-d", "--delimiter", type=str, help="The delimiter used in gene symbol column", required=False, default=";")
    parser.add_argument("-p", "--parallel", type=str, help="String indication of the boolean value abt whether to use parallel execution or not", required=False, default="True")
    parser.add_argument("-t", "--threads", type=int, help="Parallel threads", required=False, default=8)
    
    args = parser.parse_args()
    
    update_hgnc(url=args.url, 
                table=args.variant_table, 
                label=args.column_header, 
                d=args.delimiter,
                parallel=args.parallel,
                threads=args.threads,
                local_path=args.anno_path)
    # convert_to_hgnc('/Users/skyxt/Google_Drive/HKU_YANG_LAB/HKU_Academic_Issue/Project_PID_WES/Probe_Design/non_alt_loci_set.txt', '/Users/skyxt/Google_Drive/HKU_YANG_LAB/HKU_Academic_Issue/Project_PID_WES/8_samples/CNV_calling/Arshia-Ahmed.hg19_multianno.reformat.txt', 'Gene.refGene', ';')
