import pandas as pd
from scipy.stats import binom
import numpy as np
import re
from python_utils import na_value
import logging
import argparse as ap



logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)



def main_filter(input_table,
                output_table="",
                AF_cutoff=0.05,
                AF_cols=["gnomAD_exome_ALL", "gnomAD_genome_ALL", "AF_gnomAD_Controls"],
                obs_nhomoalt_col="nhomalt_gnomAD_Controls",
                af_main_col = "AF_gnomAD_Controls",
                an_main_col = "AN_Controls",
                excl_tags=[]):
    input_df = pd.read_table(input_table, low_memory=False, na_values=[".", "", "-", " "])
    logger.info("Before filtering on allele frequency and number of homozygous alternative allele carriers, the shape of input table is {}".format(input_df.shape))
    
    AF_filter_bools = input_df.apply(filter_af_per_record, cut_off = AF_cutoff, AF_cols = AF_cols, axis=1)
    excl_bools = np.array([False for i in range(0, len(input_df))])
    if len(excl_tags) > 0:
        for tag in excl_tags:
            excl_bools = np.logical_or(excl_bools, np.logical_not(input_df.loc[:, "FILTER"].astype(str).str.contains(tag)))
    
    # Now we need to filter on the number of homozygous alternative allele carriers in the gnomAD general population (even though some patients are having complex or chronic diseases, none of them are with Mendelian diseases)
    # Genotype of individuals generally fits into three types, homoref, homoalt and hetero, the occurence of three genotypes among a population fits into a multinomial distribution
    # Given the expected count of homoalt carriers (calculated based on AF, following HWE)
    # If observed count of homoalt carriers significantly larger than the expected count of homoalt carriers, then:
        # 1. The variant does not fits into HWE because the ALT allele is a young allele, the AF will increase with generation passing, current AF is an underestimation
        # 2. The homozygous form of the ALT allele is under positive selection (unlikely given the general AF is small)
    # We perform a one-way hypothesis test to calculate the p-Value of observed number of homoalt carriers (whether it is too large)
    # Since I only care the occurence of homozygous alternative allele individuals. The multinomial distribution can be reduced to binomial distribution
    nhomalt_filter_bools = input_df.apply(filter_nhomoalt_per_record, axis=1,
                                            obs_nhomoalt_col = obs_nhomoalt_col,
                                            af_main_col = af_main_col,
                                            an_main_col = an_main_col)
    
    filtered_df = input_df.loc[np.logical_or(AF_filter_bools & nhomalt_filter_bools, excl_bools), :]
    logger.info("After the filtration on the AF and No of homoalt carriers, the shape of the filtered table is {} and will be output to {}".format(filtered_df.shape, output_table))
    filtered_df.to_csv(output_table, index=False, sep="\t")
    
    
    

def gender_col_suffix(label, 
                      gender="M",
                      male_suffix="_XY",
                      female_suffix="_XX"):
    if gender == "M":
        return label + male_suffix
    elif gender == "F":
        return label + female_suffix
    


def filter_nhomoalt_per_record(row, 
                               chromosome_col="Chr", 
                               obs_nhomoalt_col="nhomalt_gnomAD_Controls",
                               af_main_col = "AF_gnomAD_Controls",
                               an_main_col = "AN_Controls",
                               sig_level = 0.05):
    total_AN_female = float(row[gender_col_suffix(an_main_col, "F")])
    # Test whether this variant is in autosomal or sex chromosome
    if re.search(r"X$", row[chromosome_col]):
        # Need to inspect within female population to test whether O/E ratio is too large
        p = pow(float(row[gender_col_suffix(af_main_col, "F")]), 2)
        n = total_AN_female
        k = float(row[gender_col_suffix(obs_nhomoalt_col, "F")])
        if n in [np.nan]:
            p_value = np.nan
        else:
            p_value = 1 - binom.cdf(k-1, n, p) if k >= 1 else 0
    elif re.search(r"Y$", row[chromosome_col]):
        p_value = np.nan        
    else:
        # Autosomal records
        p = pow(float(row[af_main_col]), 2)
        n = total_AN_female
        k = float(row[obs_nhomoalt_col])
        if n in [np.nan]:
            p_value = np.nan
        else:
            p_value = 1 - binom.cdf(k-1, n, p) if k >= 1 else 0
        
    return not p_value < sig_level
        
    
    
def filter_af_per_record(row, cut_off = 0.05, AF_cols=[]):
    numerical_values = [ row[col] for col in AF_cols if not na_value(row[col]) ]
    if len(numerical_values) > 0:
        return any([v for v in numerical_values if v < cut_off])
    else:
        return True



    
    
if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("-i", "--input_table", type=str, help="Path to the input anno table", required=True)
    parser.add_argument("-o", "--output_table", type=str, help="Path to the output anno table", required=True)
    parser.add_argument("-c", "--AF_cutoff", type=float, help="cutoff of the allele frequency", required=False, default=0.05)
    parser.add_argument("-l", "--AF_labels", type=str, help="Delimited by comma, a list of AF column labels used for filtering", required=False, default = ",".join(["gnomAD_exome_ALL", "gnomAD_genome_ALL", "AF_gnomAD_Controls"]))
    parser.add_argument("-a", "--gnomAD_nhomoalt", type=str, help="Label of column storing the gnomAD observed number of homozygous carriers", required=False, default = "nhomalt_gnomAD_Controls")
    parser.add_argument("-f", "--gnomAD_af_col", type=str, help="Label of the main AF column storing the gnomAD AFs", required=False, default = "AF_gnomAD_Controls")
    parser.add_argument("-n", "--gnomAD_an_col", type=str, help="Label of the main column storing the gnomAD observed total number of alleles", required=False, default = "AN_Controls")
    parser.add_argument("-e", "--excl_tags", type=str, help="Delimited by comma, a list of FILTER tags that does not apply AF filtering", required=False, default = "")

    args = parser.parse_args()
    af_labels = args.AF_labels.split(",")
    exl_tags = [ t for t in args.excl_tags.split(",") if len(t) > 0 ]
    
    main_filter(input_table=args.input_table,
                output_table=args.output_table,
                AF_cutoff=args.AF_cutoff,
                AF_cols=af_labels,
                obs_nhomoalt_col=args.gnomAD_nhomoalt,
                af_main_col = args.gnomAD_af_col,
                an_main_col = args.gnomAD_an_col,
                excl_tags=exl_tags)
    