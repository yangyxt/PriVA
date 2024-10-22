import logging
import pandas as pd
import numpy as np
import sqlite3
import re
import time
import uuid
import os
import sys
import argparse as ap
import math
import io
import subprocess


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)



def na_value(v):
    # logger.info("Input value is {}, and its type is {}, whether its NAN? {}, {}".format(v, type(v), v in [np.nan], math.isnan(v)))
    if type(v) == str:
        if re.search(r"^[Nn][Aa][Nn]*$", v):
            return True
        elif re.search(r"^[;,\|]*[Nn][Aa][Nn]*[;,\|]*$", v):
            return True
        elif re.search(r"^([Nn][Aa][Nn]*[;,\|_\- ]+)+([Nn][Aa][Nn]*)*$", v):
            return True
        elif re.search(r"^[\.\-\*_ ]*$", v):
            return True
        else:
            return False
    elif v is None:
        return True
    elif v in [np.nan]:
        return True
    elif type(v) == pd._libs.tslibs.timestamps.Timestamp:
        return False
    elif math.isnan(v):
        return True
    else:
        return False


def deal_with_nullkey_group_return(func):
    def check_name_wrapper(*args, **kwargs):
        first_arg = [a for a in args][0]
        res = func(*args, **kwargs)
        if type(first_arg) == pd.core.frame.DataFrame:
            try:
                group_name = first_arg.name
            except AttributeError:
                logger.debug("WARNING: This function {} is supposed to be used as a groupby apply function, while we cannot access the groupdf name".format(func.__name__))
                return res
            else:
                if check_null_groupkeys(group_name):
                    print("WARNING: Input group df to function {} has all NA group keys: {}, key type {}. The group_df has shape of {}, check the head part:\n{}".format(func.__name__,
                                                                                                                                                            group_name,
                                                                                                                                                            type(group_name),
                                                                                                                                                            first_arg.shape,
                                                                                                                                                            first_arg[:3].to_string(index=False)), file=sys.stderr)
                    if type(res) == pd.core.frame.DataFrame:
                        input_cols = first_arg.columns.tolist()
                        output_cols = res.columns.tolist()
                        new_cols = [ l for l in output_cols if l not in input_cols ]
                        if len(new_cols) == 0:
                            return res
                        else:
                            res.loc[:, new_cols] = np.nan
                            return res.drop_duplicates()
                    elif res is None:
                        pass
                    else:
                        return res
                else:
                    return res
        else:
            print("WARNING: This function {} is supposed to be used as a groupby apply function, while the first arg is not a dataframe, instead its a {}".format(func.__name__, type(first_arg)), file=sys.stderr)
            return res
    return check_name_wrapper


def check_null_groupkeys(keys):
    if type(keys) == tuple or type(keys) == list:
        return all([ na_value(k) for k in keys ])
    else:
        return na_value(keys)


def convert_input_value(v):
    if type(v) == str:
        if re.search(r"^[Tt][Rr][Uu][Ee]$", v):
            return True
        elif re.search(r"^[Ff][Aa][Ll][Ss][Ee]$", v):
            return False
        elif re.search(r"^[0-9]+$", v):
            return int(v)
        elif re.search(r"^[0-9]*\.[0-9]+$", v):
            return float(v)
        elif v == "None":
            return None
        elif re.search(r"^[Nn][Aa][Nn]$", v):
            return np.nan
        else:
            return v
    else:
        return v



def deal_with_nullkey_group_return(func):
    def check_name_wrapper(*args, **kwargs):
        first_arg = [a for a in args][0]
        res = func(*args, **kwargs)
        if type(first_arg) == pd.core.frame.DataFrame:
            try:
                group_name = first_arg.name
            except AttributeError:
                logger.debug("WARNING: This function {} is supposed to be used as a groupby apply function, while we cannot access the groupdf name".format(func.__name__))
                return res
            else:
                if check_null_groupkeys(group_name):
                    print("WARNING: Input group df to function {} has all NA group keys: {}, key type {}. The group_df has shape of {}, check the head part:\n{}".format(func.__name__,
                                                                                                                                                            group_name,
                                                                                                                                                            type(group_name),
                                                                                                                                                            first_arg.shape,
                                                                                                                                                            first_arg[:3].to_string(index=False)), file=sys.stderr)
                    if type(res) == pd.core.frame.DataFrame:
                        input_cols = first_arg.columns.tolist()
                        output_cols = res.columns.tolist()
                        new_cols = [ l for l in output_cols if l not in input_cols ]
                        if len(new_cols) == 0:
                            return res
                        else:
                            res.loc[:, new_cols] = np.nan
                            return res.drop_duplicates()
                    elif res is None:
                        pass
                    else:
                        return res
                else:
                    return res
        else:
            print("WARNING: This function {} is supposed to be used as a groupby apply function, while the first arg is not a dataframe, instead its a {}".format(func.__name__, type(first_arg)), file=sys.stderr)
            return res
    return check_name_wrapper



def extract_transcript_coordinates(gtf_file, output_tab=None):
    import gffutils
    import pandas as pd
    """
    Parses the GTF file and extracts chromosome, start, end, and transcript ID for each transcript.
    Returns a pandas DataFrame with the extracted information.

    Parameters:
    - gtf_file: Path to the GTF file.

    Returns:
    - df: A pandas DataFrame with columns ['chr', 'start', 'end', 'transcript_id'].
    """
    # Create a database from the GTF file (stored in memory)
    db = create_gffutils_db(gtf_file)

    # Initialize a list to store transcript data
    data = []

    # Iterate over all transcript features in the database
    for transcript in db.features_of_type('transcript'):
        chrom = transcript.chrom
        start = transcript.start
        end = transcript.end
        # Extract the transcript ID from attributes
        transcript_id = transcript.attributes.get('transcript_id', [''])[0]
        data.append([chrom, start, end, transcript_id])

    # Create a DataFrame from the data list
    df = pd.DataFrame(data, columns=['chr', 'start', 'end', 'transcript_id'])

    if output_tab:
        df.to_csv(output_tab, sep='\t', index=False)

    return df



def create_gffutils_db(gff_file="/paedyl01/disk1/yangyxt/public_data/gene_annotation/GCF_000001405.25_GRCh37.p13_genomic.gff.gz",
                       db_file = None, **kwargs):
    import gffutils

    if not db_file:
        db_file = re.sub("(\.g[tf]f[3]*)\.gz$", "\1.db", gff_file)
        logger.info(f"The output database file is {db_file}")

    # kwargs here can have merge_strategy = "merge" or merge_strategy = "create_unique"
    db = gffutils.create_db(gff_file,
                            dbfn=db_file,
                            force=True,
                            keep_order=True,
                            sort_attribute_values=True,
                            merge_strategy = "merge",
                            **kwargs)

    return db


def update_tranx_span(target_tsv, updated_span_tsv, output_tsv=None):
    # Expect the target_tsv to be a gz tsv file with the last three columns being:
    #   chromosome, start_position, end_position
    # And the second column being the transcript_id (column name: transcript)
    # Expect the updated_span_tsv to be a plain tsv file with the first three columns being:
    #   chr, start, end
    # And the fourth column being the transcript_id (column name: transcript_id)
    import pandas as pd
    import numpy as np
    target_df = pd.read_table(target_tsv, compression="gzip", low_memory=False)
    # Rename the transcript column to transcript_id
    target_df.rename(columns={"transcript": "transcript_id"}, inplace=True)
    updated_span_df = pd.read_table(updated_span_tsv, low_memory=False)

    # Merge the two dataframes on transcript_id
    merged_df = pd.merge(target_df, updated_span_df, on="transcript_id", how="left")

    # Drop the rows where start or end is NaN
    merged_df = merged_df.dropna(subset=["start", "end"])

    # Update the start and end position
    merged_df.loc[:, "start_position"] = merged_df.loc[:, "start"].astype(np.int64)
    merged_df.loc[:, "end_position"] = merged_df.loc[:, "end"].astype(np.int64)
    # Drop the intermediate columns
    merged_df.drop(columns=["start", "end", "chr"], inplace=True)

    # Write the output to a gzipped tsv file
    if not output_tsv:
        output_tsv = re.sub("(\.tsv\.gz)$", ".grch38.tsv.gz", target_tsv)
    merged_df.to_csv(output_tsv, sep="\t", index=False, compression="gzip")
    return output_tsv



def check_remote_file_exists(url):
    import requests
    try:
        response = requests.head(url)
        if response.status_code == 200:
            return True
        else:
            return False
    except Exception as e:
        print(f"An error occurred: {e}")
        return False


def get_updated_clinvar_reformat_tab(clinvar_dir = "/home/yangyxt/public_data/clinvar",
                                     assembly = "GRCh37",
                                     **kwargs):
    from clinvar_class import UPDATED_CLINVAR

    ncbi2ucsc_as_dict = {"GRCh38": "hg38",
                        "hg38": "hg38",
                        "GRCh37": "hg19",
                        "hg19": "hg19"}

    ucsc_assembly = ncbi2ucsc_as_dict[assembly]

    clinvar_dir = os.path.join(clinvar_dir, ucsc_assembly)
    clinvar_obj = UPDATED_CLINVAR(clinvar_dir = clinvar_dir,
                                     assembly = assembly, **kwargs)

    clinvar_obj.prepare_updated_reformat_clinvar()


if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("-f", "--function", type=str, help="The function name", required=True)
    parser.add_argument("-a", "--arguments", type=str, help="The function's input arguments, delimited by semi-colon ;", required=False, default=None)
    parser.add_argument("-k", "--key_arguments", type=str, help="Keyword arguments for the function, delimited by semi-colon ;", required=False, default=None)

    args = parser.parse_args()
    try:
        fargs = [ convert_input_value(a) for a in args.arguments.split(";") ] if type(args.arguments) == str else []
        fkwargs = { t.split("=")[0]: convert_input_value(t.split("=")[1]) for t in args.key_arguments.split(";") } if type(args.key_arguments) == str else {}
        logger.info("Running function: {}, input args are {}, input kwargs are {}".format(args.function, fargs, fkwargs))
    except Exception as e:
        logger.error("Input argument does not meet the expected format, encounter Parsing error {}, Let's check the input:\n-f {}, -a {}, -k {}".format(
            e,
            args.function,
            args.arguments,
            args.key_arguments
        ))
        raise e
    globals()[args.function](*fargs, **fkwargs)