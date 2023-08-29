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