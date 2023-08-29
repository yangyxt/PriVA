from unittest import result
import pandas as pd
import numpy as np
import logging
import os
import re
import sqlite3
import argparse as ap
from functools import partial
from multiprocessing import Pool
from itertools import zip_longest
from concurrent.futures import ProcessPoolExecutor
import vcfpy


logger = logging.getLogger()
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(processName)s:%(funcName)s:%(lineno)s:%(name)s:%(exc_info)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)
logger.setLevel(os.environ.get("LOGLEVEL", "INFO"))


def annotate_expected_nhomo(row, 
                            chr_col="Chr", 
                            AN_col="AN_Controls", 
                            AC_col="AC_gnomAD_Controls", 
                            AF_col="AF_gnomAD_Controls",
                            new_col="Expected_nhomalt_Controls"):
    if re.search(r"^[Cc]*[Hh]*[Rr]*X$", str(row[chr_col])):
        # For male
        row[AF_col + "_XY"] = float(row[AC_col + "_XY"])/float(row[AN_col + "_XY"]) if float(row[AN_col + "_XY"]) else np.nan
        row[new_col + "_XY"] = row[AC_col + "_XY"]
        # For female
        row[AF_col + "_XX"] = float(row[AC_col + "_XX"])/float(row[AN_col + "_XX"]) if float(row[AN_col + "_XX"]) else np.nan
        row[new_col + "_XX"] = pow(float(row[AF_col + "_XX"]), 2) * float(row[AN_col + "_XX"])/2
        # For both gender
        row[AF_col] = float(row[AC_col])/float(row[AN_col]) if float(row[AN_col]) else np.nan
        row[new_col] = row[new_col + "_XY"] + row[new_col + "_XX"]
    elif re.search(r"^[Cc]*[Hh]*[Rr]*Y$", str(row[chr_col])):
        # For male
        row[AF_col + "_XY"] = float(row[AC_col + "_XY"])/float(row[AN_col + "_XY"]) if float(row[AN_col + "_XY"]) else np.nan
        row[new_col + "_XY"] = row[AC_col + "_XY"]
        # For female (NA)
        row[AF_col + "_XX"] = np.nan
        row[new_col + "_XX"] = np.nan
        # For both gender in general
        row[AF_col] = float(row[AC_col])/float(row[AN_col]) if float(row[AN_col]) else np.nan
        row[new_col] = row[new_col + "_XY"]
    else:
        # For male
        row[AF_col + "_XY"] = float(row[AC_col + "_XY"])/float(row[AN_col + "_XY"]) if float(row[AN_col + "_XY"]) else np.nan
        row[new_col + "_XY"] = pow(float(row[AF_col + "_XY"]), 2) * float(row[AN_col + "_XY"])/2
        # For female
        row[AF_col + "_XX"] = float(row[AC_col + "_XX"])/float(row[AN_col + "_XX"]) if float(row[AN_col + "_XX"]) else np.nan
        row[new_col + "_XX"] = pow(float(row[AF_col + "_XX"]), 2) * float(row[AN_col + "_XX"])/2
        # For both gender
        row[AF_col] = float(row[AC_col])/float(row[AN_col]) if float(row[AN_col]) else np.nan
        row[new_col] = row[new_col + "_XY"] + row[new_col + "_XX"]
    return row


def check_var_presence( input_variant,
                        table_name = "gnomAD_total_variants",
                        gnomAD_sql_template="/paedyl01/disk1/yangyxt/public_data/gnomAD/total_vcf_db/gnomAD_total_variants.v3.1.2.hg19.{}.sql",
                        contigs=["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                                 "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                 "chr20", "chr21", "chr22", "chrX", "chrY"],
                        target_cols = [ "AC_Controls", "nhomalt_Controls", "AN_Controls", 
                                        "AC_Controls_XX", "nhomalt_Controls_XX", "AN_Controls_XX",
                                        "AC_Controls_XY", "nhomalt_Controls_XY", "AN_Controls_XY" ], 
                        match_cols = {"CHROM": "Chr", "POS": "Start", "REF": "Ref", "ALT": "Alt"}):
    if type(input_variant) == vcfpy.Record:
        contig = input_variant.CHROM
        pos = input_variant.POS
        ref = input_variant.REF
        alt = input_variant.ALT[0].serialize()
    elif type(input_variant) == pd.core.series.Series:
        col_names = input_variant.index.tolist()
        logger.info("Take a look at the input pd series (row) column names: \n{}\n".format(col_names))
        if (all([v in col_names for k,v in match_cols.items()])):
            contig = input_variant[match_cols["CHROM"]]
            pos = input_variant[match_cols["POS"]]
            ref = input_variant[match_cols["REF"]]
            alt = input_variant[match_cols["ALT"]]
        else:
            logger.warning("The input pandas row does not match with Chr,Ref,Alt,Start. Let's take a look at the input type: {}, and content: {}".format(type(input_variant), input_variant))
            match_cols = { ["CHROM", "POS", "REF", "ALT"][i]:col_names[[0,1,3,4]][i] for i in range(0,4) }
            contig = input_variant[match_cols["CHROM"]]
            pos = input_variant[match_cols["POS"]]
            ref = input_variant[match_cols["REF"]]
            alt = input_variant[match_cols["ALT"]]
    else:
        # logger.debug("The input variant record is an instance of {}, show the content:\n{}".format(type(input_variant), input_variant))
        contig = input_variant[match_cols["CHROM"]]
        pos = input_variant[match_cols["POS"]]
        ref = input_variant[match_cols["REF"]]
        alt = input_variant[match_cols["ALT"]]
            
    if contig in contigs:
        gnomAD_sql_path = gnomAD_sql_template.format(contig)
        assert len(gnomAD_sql_path) > 0
        conn = create_index_sql(table_name=table_name,
                                index_cols=["CHROM", "POS", "REF", "ALT"],
                                sql_path = gnomAD_sql_path)
        if conn is not None:
            cursor = conn.cursor()
            query_clause = '''
                           SELECT {tc} 
                           FROM "{tn}"
                           WHERE "CHROM" = "{cg}" AND "POS" = {ps} AND "REF" = "{rf}" AND "ALT" = "{at}";  
                           '''.format(tc = ", ".join([ '"' + c + '"' for c in target_cols ]),
                                      tn = table_name,
                                      cg = contig,
                                      ps = pos,
                                      rf = ref,
                                      at = alt)
            logger.debug("Before executing this clause, check the SQL phrase: \n{}\n".format(query_clause))
            try:
                results = cursor.execute(query_clause).fetchall()
            except sqlite3.OperationalError as se:
                logger.error("This SQL clause caused error: {}\n".format(query_clause))
                raise se
            if len(results) == 0:
                if type(input_variant) == pd.core.series.Series:
                    for tc in target_cols: input_variant[tc] = np.nan
                    return input_variant
                else: 
                    return np.nan
            elif len(results) > 1:
                logger.warning("This variant input matches with multiple records in gnomAD database: {}\n The matched records are: \n{}\n".format(input_variant, results))
                if type(input_variant) == pd.core.series.Series:
                    for i in range(0, len(target_cols)): input_variant[target_cols[i]] = results[0][i]
                    return input_variant
                else:
                    return results[0]
            else:
                if type(input_variant) == pd.core.series.Series:
                    for i in range(0, len(target_cols)): input_variant[target_cols[i]] = results[0][i]
                    return input_variant
                else:
                    return results[0]
        else:
            if type(input_variant) == pd.core.series.Series:
                for tc in target_cols: input_variant[tc] = np.nan
                return input_variant
            else: 
                return np.nan
    else:
        if type(input_variant) == pd.core.series.Series:
            for tc in target_cols: input_variant[tc] = np.nan
            return input_variant
        else: 
            return np.nan



def parallel_check_var_presence(input_variants,
                                pool_object = None,
                                result_dict = dict(),
                                table_name = "gnomAD_total_variants",
                                gnomAD_sql_template="/paedyl01/disk1/yangyxt/public_data/gnomAD/total_vcf_db/gnomAD_total_variants.v3.1.2.hg19.{}.sql",
                                contigs=["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", 
                                            "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                            "chr20", "chr21", "chr22", "chrX", "chrY"],
                                target_cols = [ "AC_Controls", "nhomalt_Controls", "AN_Controls", 
                                                "AC_Controls_XX", "nhomalt_Controls_XX", "AN_Controls_XX",
                                                "AC_Controls_XY", "nhomalt_Controls_XY", "AN_Controls_XY" ], 
                                match_cols = {"CHROM": "Chr", "POS": "Start", "REF": "Ref", "ALT": "Alt"}):
    check_var_presence_bundle = partial(check_var_presence, 
                                          gnomAD_sql_template = gnomAD_sql_template, 
                                          contigs = contigs, 
                                          target_cols = target_cols, 
                                          match_cols = match_cols)
    result_iterator = pool_object.imap_unordered(check_var_presence_bundle, input_variants)
    while True:
        try:
            result_row = next(result_iterator)
        except StopIteration:
            break
        else:
            contig = result_row["Chr"]
            result_dict[contig].append(result_row)
    return result_dict
        



def on_completion(future, progress_dict, iterator_dict, result_dict, gnomAD_db_template):
    exception = future.exception()
    if exception:
        raise exception

    result = future.result()
    group_index = progress_dict[future]
    # Add the resulting row to the corresponding chromosome set
    result_dict[group_index].add(result)
    
    row = next(iterator_dict[group_index], None)
    
    if row is not None:
        # Schedule a new task from the same group
        new_future = executor.submit(check_var_presence, row, gnomAD_sql_template = gnomAD_db_template)
        new_future.add_done_callback(on_completion)
        active_futures[new_future] = group_index
    del active_futures[future]
    

def create_index_sql(table_name, index_cols, index_name = "idx_of_cols", sql_path = "", db_conn = None):
    if db_conn is None:
        if len(sql_path) > 0:
            db_conn = create_sqlite_db(sql_path)
        else:
            logger.error("The input connection is None and the input DB path is none.")
            raise ValueError
        
    create_index_clause = '''
                          CREATE UNIQUE INDEX IF NOT EXISTS "{idn}" ON {tabn}({cl});
                          '''.format(idn = index_name,
                                     tabn = table_name, 
                                     cl = ", ".join(index_cols))
                          
    check_index_clause = 'PRAGMA index_list("{}");'.format(table_name)
    
    if db_conn is not None:
        cur = db_conn.cursor()
        idx_names = [ idx[1] for idx in cur.execute(check_index_clause).fetchall() ]
        if index_name in idx_names:
            return db_conn
        else:
            cur = cur.execute(create_index_clause)
            db_conn.commit()
            cur.close()
            return db_conn       
    else:
        logger.error("The input connection is None and the input DB path is none.")
        raise ValueError
    
    

def create_sqlite_db(dbpath=":memory:", force=False):
    from pathlib import Path
    if os.path.exists(dbpath):
        if force:
            os.remove(dbpath)
            Path(dbpath).touch()
    else:
        Path(dbpath).touch()
            
    conn = None
    try:
        conn = sqlite3.connect(dbpath, timeout=360)
        logger.debug("Establish connection with database at {} with sqlite3 version {}".format(dbpath,sqlite3.sqlite_version))
    except Exception as e:
        logger.warning("Failed to establish connection with database at {} with sqlite3 version {}".format(dbpath,sqlite3.sqlite_version))
        logger.error(e)
        raise e
    finally:
        return conn
    
    

def main_get_homo(var_table, 
                  threads=8, 
                  db_dir="",
                  assembly="hg19",
                  fetched_cols={"AF_Controls":"AF_gnomAD_Controls", 
                                "AF_Controls_XX":"AF_gnomAD_Controls_XX", 
                                "AF_Controls_XY":"AF_gnomAD_Controls_XY", 
                                "nhomalt_Controls":"nhomalt_gnomAD_Controls", 
                                "nhomalt_Controls_XX":"nhomalt_gnomAD_Controls_XX", 
                                "nhomalt_Controls_XY":"nhomalt_gnomAD_Controls_XY",
                                "AC_Controls":"AC_gnomAD_Controls", 
                                "AC_Controls_XX":"AC_gnomAD_Controls_XX", 
                                "AC_Controls_XY":"AC_gnomAD_Controls_XY"}, 
                  col_index="gnomAD_genome_ALL",
                  output_table=None):
    '''
    Prepare the table storing variants, by default input_variant is a pandas series
    '''
    if type(var_table) == str:
        var_table_path = var_table
        dtypes = pd.read_table(var_table_path, low_memory=False).dtypes.to_dict()
        logger.info("Here are the dataframe column data types:\n{}\n".format(dtypes))
        var_table = pd.read_table(var_table_path, dtype=dtypes)
        
    logger.info("Take a look at the input variant table:\n{}\n".format(var_table[:5].to_string(index=False)))    
     
    gnomAD_db = db_dir + f"/{assembly}/gnomAD_total_variants.{assembly}" + ".{}.sql"
    # First split the table by chromosome
    by_chr_df = var_table.groupby("Chr", as_index=False)
    result_by_chr = {k: [] for k, it in by_chr_df}
    # Start to process the table with multiple threads
    
    '''
    The code below using concurrent futures to realize parallel execution but it's not pickle friendly. Deprecated now
    ==========================================================================================================================================
    with ProcessPoolExecutor(max_workers=threads) as executor:
        # Dictionary the progress of subprocesses for convenience of tracing
        active_futures = {}
        
        # Schedule the first task for each group
        for group_index, it in iterator_by_chr.items():
            row = next(it, None)
            if row is not None:
                future = executor.submit(check_var_presence, row, gnomAD_sql_template = gnomAD_db)
                call_back = partial(on_completion, progress_dict = active_futures, 
                                           iterator_dict = iterator_by_chr, 
                                           result_dict = result_by_chr,
                                           gnomAD_db_template = gnomAD_db)
                future.add_done_callback(call_back)
    '''
    iter_dict = { k: [ s for i, s in g.iterrows() ] for k, g in by_chr_df }
    
    with Pool(threads) as pool:
        for rows in zip_longest(*list(iter_dict.values())):
            row_series_list = [row for row in rows if row is not None ]
            if len(row_series_list) > 0:
                logger.info("Take a look at one of the row object passed to the parallel_check_var_presence function:\n{}\n".format(row_series_list[0]))
                result_by_chr = parallel_check_var_presence(row_series_list, 
                                                            pool_object= pool, 
                                                            result_dict=result_by_chr, 
                                                            gnomAD_sql_template=gnomAD_db)
        
    # Now result by chr dict should be a dict of variant records (pd.Series format) that belonging to different chromosomes
    # First establish an empty pandas dataframe, then assigning each row values using iloc
    total_rows = sum(len(series_set) for series_set in result_by_chr.values())
    one_series_set = result_by_chr[next(iter(result_by_chr))]
    columns = next(iter(one_series_set)).index.tolist()
    logger.info("The pd.series stored in result_dict acquired by parallel execution of check_var_presence looks like :\n{}\nAnd the column names are:\n{}\n".format(next(iter(one_series_set)).to_string(),
                                                                                                                                                                    columns))
    result_df = pd.DataFrame(index=range(total_rows), columns=columns)
    logger.info("The empty dataframe prepared has these columns: \n{}\n".format(result_df[:1].to_string(index=False)))
    
    row_index = 0
    for series_set in result_by_chr.values():
        for series in series_set:
            if series.index.tolist() != result_df.columns.tolist():
                raise ValueError("The series column names does not match with the empty table's column labels. Series index:\n{}\nDataframe columns:\n{}\n".format(series.index.tolist(),
                                                                                                                                                                   result_df.columns.tolist()))
            result_df.iloc[row_index] = series
            row_index += 1
            
    logger.info("The result dataframe returned by check_var_presence looks like :\n{}\n".format(result_df[:5].to_string(index=False)))
    
    for col in ["AC_Controls", "AC_Controls_XX", "AC_Controls_XY"]:
        result_df[col.replace("AC", "AF")] = np.where(result_df[col.replace("AC", "AN")] > 0, result_df[col]/result_df[col.replace("AC", "AN")].replace(0, np.nan), np.nan)
            
    result_df = result_df.rename(columns=fetched_cols)
    logger.info("After getting gnomAD nhomalt annotation, The result df looks like: \n{}".format(result_df[:4].to_string(index=False)))
    
    # Add expected nhomo column
    result_df = result_df.apply(annotate_expected_nhomo, axis=1)
    logger.info("After generating expected nhomalt on, The result df looks like: \n{}".format(result_df[:4].to_string(index=False)))
    
    # Move the columns into specified column index.
    if col_index != "-1":
        new_cols = [x for x in result_df.columns.tolist() if x not in var_table.columns.tolist()]
        new_cols.reverse()
        for col in new_cols:
            tobe_moved = result_df.pop(col)
            try:
                insert_ind = int(col_index)
            except ValueError:
                insert_ind = result_df.columns.get_loc(col_index) + 1
            result_df.insert(insert_ind, col, tobe_moved)
            result_df = result_df.copy()
            

    if 'var_table_path' in locals():
        if not output_table:
            output_table = var_table_path
        result_df.to_csv(output_table, sep='\t', index=False)
    else:
        return result_df


if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("-v", "--variant_table", type=str, help="Input path of the variant table", required=False)
    parser.add_argument("-i", "--col_index", type=str, help="Column index or column header to specify the place of new annotated columns.", required=False, default="gnomAD_genome_ALL")
    parser.add_argument("-a", "--assembly", type=str, help="hg19 or hg38", required=False, default="hg19")
    parser.add_argument("-d", "--db_dir", type=str, help="The path storing the aggregated gnomD sqlite DB for both hg19 and hg38", required=True)
    parser.add_argument("-t", "--threads", type=int, help="Threads to be used by multiprocessing", required=False, default=8)
    parser.add_argument("-o", "--output_table", type=str, help="Absolute file path of the output table", required=False, default=None)
    # parser.add_argument("-c", "--fetched_cols", type=str, help="The column headers to be extracted and to be reserved in final df", required=False, default="AF=AF_gnomAD_v3_genome_ALL,AF_male=AF_male_gnomAD_v3_genome,AF_female=AF_female_gnomAD_v3_genome,nhomalt=nhomalt_gnomAD_v3_genome,nhomalt_male=nhomalt_male_gnomAD_v3_genome,nhomalt_female=nhomalt_female_gnomAD_v3_genome")
    args = parser.parse_args()

    # fetched_cols = { x.split("=")[0]:x.split("=")[1] for x in args.fetched_cols.split(",") }
    main_get_homo(var_table = args.variant_table, 
                  threads = args.threads, 
                  db_dir = args.db_dir,
                  assembly = args.assembly,
                  col_index = args.col_index,
                  output_table = args.output_table)





    
