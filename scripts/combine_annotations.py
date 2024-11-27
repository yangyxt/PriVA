import pandas as pd
import numpy as np
import multiprocessing as mp
import pysam
import argparse
import json



def convert_record_to_tab(record) -> pd.DataFrame:
    return pd.DataFrame()


def convert_vcf_to_tab(input_vcf, threads):
    with pysam.VariantFile(input_vcf) as vcf_file, mp.Pool(threads) as pool:
        results = pool.map(convert_record_to_tab, vcf_file)


