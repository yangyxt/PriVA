import pandas as pd
import numpy as np 
import pysam as ps 
import pybedtools as pb 
import subprocess
from subprocess import PIPE
import re
import argparse as ap
import logging
import sys
import os
import urllib
import pprint
import swifter
from swifter import set_defaults
set_defaults(npartitions=10)
from check_qname_count import id_generator
from python_utils import deal_with_nullkey_group_return, na_value
from prepare_updated_clinvar_recs import UPDATED_CLINVAR

'''
This script is used to annotate intervals overlapping with the same interval.
'''

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def make_grantham_dict(grantham_mat_file="/paedyl01/disk1/yangyxt/public_data/Interpro/Grantham_matrix/Gratham.matrix.txt"):

    """
    Citation:   http://www.ncbi.nlm.nih.gov/pubmed/4843792
    Provenance: http://www.genome.jp/dbget-bin/www_bget?aaindex:GRAR740104
    .	A	R	N	D	C	Q	E	G	H	I	L	K	M	F	P	S	T	W	Y	V
    A	0	112	111	126	195	91	107	60	86	94	96	106	84	113	27	99	58	148	112	64
    R	112	0	86	96	180	43	54	125	29	97	102	26	91	97	103	110	71	101	77	96
    N	111	86	0	23	139	46	42	80	68	149	153	94	142	158	91	46	65	174	143	133
    D	126	96	23	0	154	61	45	94	81	168	172	101	160	177	108	65	85	181	160	152
    C	195	180	139	154	0	154	170	159	174	198	198	202	196	205	169	112	149	215	194	192
    Q	91	43	46	61	154	0	29	87	24	109	113	53	101	116	76	68	42	130	99	96
    E	107	54	42	45	170	29	0	98	40	134	138	56	126	140	93	80	65	152	122	121
    G	60	125	80	94	159	87	98	0	98	135	138	127	127	153	42	56	59	184	147	109
    H	86	29	68	81	174	24	40	98	0	94	99	32	87	100	77	89	47	115	83	84
    I	94	97	149	168	198	109	134	135	94	0	5	102	10	21	95	142	89	61	33	29
    L	96	102	153	172	198	113	138	138	99	5	0	107	15	22	98	145	92	61	36	32
    K	106	26	94	101	202	53	56	127	32	102	107	0	95	102	103	121	78	110	85	97
    M	84	91	142	160	196	101	126	127	87	10	15	95	0	28	87	135	81	67	36	21
    F	113	97	158	177	205	116	140	153	100	21	22	102	28	0	114	155	103	40	22	50
    P	27	103	91	108	169	76	93	42	77	95	98	103	87	114	0	74	38	147	110	68
    S	99	110	46	65	112	68	80	56	89	142	145	121	135	155	74	0	58	177	144	124
    T	58	71	65	85	149	42	65	59	47	89	92	78	81	103	38	58	0	128	92	69
    W	148	101	174	181	215	130	152	184	115	61	61	110	67	40	147	177	128	0	37	88
    Y	112	77	143	160	194	99	122	147	83	33	36	85	36	22	110	144	92	37	0	55
    V	64	96	133	152	192	96	121	109	84	29	32	97	21	50	68	124	69	88	55	0
    """

    aa_map = { "A":"Ala", "R":"Arg", "N":"Asn", "D":"Asp", "C":"Cys", "Q":"Gln", "E":"Glu","G":"Gly","H":"His","I":"Ile", "L":"Leu","K":"Lys","M":"Met","F":"Phe","P":"Pro","O":"Pyl","S":"Ser","U":"Sec","T":"Thr","W":"Trp","Y":"Tyr","V":"Val","B":"Asx","Z":"Glx","X":"Xaa","J":"Xle"}
    
    with open(grantham_mat_file, mode="r") as gf:
        headers = pd.Series(gf.readline().strip().split("\t")).map(aa_map).fillna(pd.Series(gf.readline().strip().split("\t")))
        # header = f.next().strip().split('\t')
        idx_to_aa = dict(zip(range(0,len(headers)), headers))

        grantham_dict = {}
        for line in gf.readlines():
            fields = line.strip().split('\t')
            from_aa = aa_map.get(fields[0], fields[0])

            for idx, score in enumerate(fields):
                if idx == 0:
                    continue
                to_aa = idx_to_aa[idx]
                grantham_dict[(from_aa, to_aa)] = int(score)

    return grantham_dict


def return_bed_instance(table_path, outputstr=False):
    '''
    If we are dealing with a bed file with header line, then there wont be a pure digit column.
    '''
    # First detect whether table first row contain pure digit column. 
    with open(table_path, mode="r") as t:
        first_line_cols = t.readline().strip("\n ").split("\t")
    
    digit_cols = [ v for v in first_line_cols if re.search(r'^[0-9]+$', v) ]
    if len(digit_cols) > 0:
        # Then the table does not contain header, it is a direct bed format file.
        if outputstr:
            return (os.path.abspath(table_path), tuple(["" for i in range(1, len(first_line_cols)+1)]))
        else:
            return (pb.BedTool(table_path), tuple(["" for i in range(1, len(first_line_cols)+1)]))
    else:
        # The table contain a header row.
        bed_path = ".".join(table_path.split(".")[:-1]) + ".bed"
        tab_dir = os.path.dirname(table_path)
        if os.path.exists(bed_path):
            tmp_bed_path = bed_path[:-4] + id_generator() + ".bed"
            pd.read_csv(table_path, sep='\t', encoding="utf-8").to_csv(tmp_bed_path, sep='\t', index=False, header=False)
            cmd = "mv {} {}".format(tmp_bed_path, bed_path)
            subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8", check=True)
        else: 
            pd.read_csv(table_path, sep='\t', encoding="utf-8").to_csv(bed_path, sep='\t', index=False, header=False)
        if outputstr:
            return (os.path.abspath(bed_path), tuple(first_line_cols))
        else:
            return (pb.BedTool(bed_path), tuple(first_line_cols))


def find_shared_rec(left_df, right_df, mid_df, output_bed=None, left_drop_cols=[], mid_drop_cols=[], right_drop_cols=[]):
    '''
    To be dropped cols relevant parameter should be 1-based column indices
    '''
    # rec1_tup = return_bed_instance(left_table, outputstr=True)
    rec1_bed = pb.BedTool.from_dataframe(left_df)
    logger.info("The Left table looks like:\n{}".format(left_df[:5].to_string(index=False)))
    rec1_headers = left_df.columns.tolist()

    # rec2_tup = return_bed_instance(right_table, outputstr=True)
    logger.info("The Right table looks like:\n{}".format(right_df[:5].to_string(index=False)))
    rec2_headers = right_df.columns.tolist()

    # mid_tup = return_bed_instance(mid_bed, outputstr=True)
    mid_bed = pb.BedTool.from_dataframe(mid_df)
    logger.info("The Middle table looks like:\n{}".format(mid_df[:5].to_string(index=False)))
    mid_headers = mid_df.columns.tolist()

    if "" in mid_headers:
        if "" in rec1_headers and "" in rec2_headers:
            all_headers = ["feature" + str(i) for i in range(1, len(rec1_headers+mid_headers+rec2_headers)+1)]
            rec1_headers = all_headers[:len(rec1_headers)]
            rec2_headers = all_headers[-len(rec2_headers):]
            mid_headers = [ h for h in all_headers if h not in rec1_headers and h not in rec2_headers ]
        elif "" in rec1_headers:
            both_headers = ["feature" + str(i) for i in range(1, len(rec1_headers+mid_headers)+1)]
            rec1_headers = both_headers[:len(rec1_headers)]
            mid_headers = both_headers[-len(mid_headers):]
        elif "" in rec2_headers:
            both_headers = ["feature" + str(i) for i in range(1, len(rec2_headers+mid_headers)+1)]
            mid_headers = both_headers[:len(mid_headers)]
            rec2_headers = both_headers[-len(rec2_headers):]
        else:
            mid_headers = ["feature" + str(i) for i in range(1, len(mid_headers)+1)]
    elif "" in rec1_headers and "" in rec2_headers:
        both_headers = ["feature" + str(i) for i in range(1, len(rec1_headers+rec2_headers)+1)]
        rec1_headers = both_headers[:len(rec1_headers)]
        rec2_headers = both_headers[-len(rec2_headers):]
    elif "" in rec1_headers:
        rec1_headers = ["feature" + str(i) for i in range(1, len(rec1_headers)+1)]
    elif "" in rec2_headers:
        rec2_headers = ["feature" + str(i) for i in range(1, len(rec2_headers)+1)]

    rec2_headers = [ h + "_r" if h in rec1_headers else h for h in rec2_headers ]
    mid_headers = [ h + "_mid" if h in rec1_headers or h in rec2_headers else h for h in mid_headers ]

    rec1_mid_bed = rec1_bed.intersect(mid_bed, wao=True, f=0.8, F=0.8, e=True)
    # mid_rec2_bed = mid_bed.intersect(rec2_bed, wao=True, f=0.8, F=0.8, e=True)
    dtype_dict = {h:str for h in mid_headers}
    rec1_mid_table = rec1_mid_bed.to_dataframe(disable_auto_names=True, names=rec1_headers + mid_headers + ["overlap_len"], dtype=str).dropna(how="all").drop(columns="overlap_len").drop_duplicates().astype(dtype_dict, errors="ignore").replace({"nan":np.nan, ".":np.nan, "-":np.nan})
    logger.info("The Left-Middle merged table has shape of {}, looks like:\n{}".format(rec1_mid_table.shape, rec1_mid_table[:5].to_string(index=False)))
    mid_rec2_table = pd.read_table(UPDATED_CLINVAR().prepare_domain_clinvar_overlap().replace(".tsv", ".shortv.tsv"), low_memory=False).dropna(how="all").astype(dtype_dict, errors="ignore").drop(columns=right_drop_cols, errors="ignore").replace({"nan":np.nan, ".":np.nan, "-":np.nan}).drop_duplicates().dropna(subset=["Type"])
    logger.info("The Middle-Right merged table has shape of {}, looks like:\n{}".format(mid_rec2_table.shape, mid_rec2_table[:5].to_string(index=False)))
    logger.info("Merge these two dataframe on these columns (mid_headers): {}".format(mid_headers))
    
    merged_table = rec1_mid_table.merge(mid_rec2_table.dropna(subset=[ l for l in mid_headers if l not in ["chr_start", "chr_end"] ], 
                                                              how="all"), 
                                        on=mid_headers, how="left").drop_duplicates()

    merged_table.drop(columns=mid_drop_cols, errors="ignore", inplace=True)
    merged_table.drop_duplicates(inplace=True)
    logger.info("The Left-Middle-Right merged table has shape of {}, and it looks like:\n{}".format(merged_table.shape, merged_table[:5].to_string(index=False)))
    
    if output_bed:
        merged_table.to_csv(output_bed, sep='\t', index=False)
    else:
        return merged_table
    

def replace_re_group(matchobj, pattern_dict):
    return pattern_dict.get(matchobj.group(1),matchobj.group(1))


def find_same_AA_change(row):
    aa_list = [ x.split(":")[-1] for x in str(row['AAChange.refGene']).split(",") ]
    mod_aa_list = [ x[:-1]+"=" if re.search(r'p\.(?P<aa>[A-Z]+)([0-9]+)(?P=aa)$', x) else x for x in aa_list ]
    aa_map = { "A":"Ala", "R":"Arg", "N":"Asn", "D":"Asp", "C":"Cys", "Q":"Gln", "E":"Glu","G":"Gly","H":"His","I":"Ile", "L":"Leu","K":"Lys","M":"Met","F":"Phe","P":"Pro","O":"Pyl","S":"Ser","U":"Sec","T":"Thr","W":"Trp","Y":"Tyr","V":"Val","B":"Asx","Z":"Glx","X":"Xaa","J":"Xle"}
    final_aa_list = []
    for s in mod_aa_list:
        s = re.sub(r'([A-Z])', lambda m: aa_map.get(m.group(1),m.group(1)), s)
        final_aa_list.append(s)
    logger.debug(final_aa_list)
    try:
        clinvar_aa_change = re.search(r'\((p\..+)\)$', str(row['shared_domain_Clinvar_name'])).group(1)
    except AttributeError:
        clinvar_aa_change = "."
        clinvar_residue_pos = "0"
    if clinvar_aa_change in final_aa_list and clinvar_aa_change != ".":
        return "Yes"
    elif re.search(r'p\.[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}', clinvar_aa_change):
        grantham_dict = make_grantham_dict()  # Get the AA change score matrix
        matched_obj = re.search(r'p\.([A-Z][a-z]{2})([0-9]+)([A-Z][a-z]{2})', clinvar_aa_change)
        from_aa = matched_obj.group(1)
        clinvar_residue_pos = matched_obj.group(2)
        to_aa = matched_obj.group(3)
        clin_aa_change_grandis = grantham_dict.get((from_aa, to_aa), 0)
        aa_change_list = [ (re.search(r'p\.([A-Z][a-z]{2})[0-9]+([A-Z][a-z]{2})', aa_change).group(1), re.search(r'p\.([A-Z][a-z]{2})[0-9]+([A-Z][a-z]{2})', aa_change).group(2), re.search(r'p\.([A-Z][a-z]{2})([0-9]+)([A-Z][a-z]{2})', aa_change).group(2)) for aa_change in final_aa_list if re.search(r'p\.([A-Z][a-z]{2})[0-9]+([A-Z][a-z]{2})', aa_change) ]
        aa_change_gd_list = [ (grantham_dict.get(tup, 1000), tup) for tup in aa_change_list if tup[-1] == clinvar_residue_pos ]
        if len(aa_change_gd_list) == 0:
            # for tup in aa_change_list:
            #     logger.info("The clinvar record happened at residue {} for gene {}, while the current variant is at residue {}".format(clinvar_residue_pos, str(row['Gene.refGene']), tup[-1]))
            return "Not_the_same_residue"
        similar_changes = [ x for x in aa_change_gd_list if abs(int(x[0]) - int(clin_aa_change_grandis)) <= 25 ]
        if len(similar_changes) > 0:
            for x in similar_changes:
                logger.warning("The variant changes from AA {} to AA {}, the corresponding Grantham Distance is {}. Well the recorded variant in Clinvar changes AA {} to AA {}, and the corresponding Grantham Distance is {}.".format(x[1][0], x[1][1], x[0], from_aa, to_aa, clin_aa_change_grandis))
            return "Similar"
        else:
            return "No"
    else:
        return "."
    

@deal_with_nullkey_group_return
def summary_clinvar_per_domain(groupdf, one_star_label="same_domain_CLNSIG_C1", two_star_label="same_domain_CLNSIG_C2"):
    # Prepare the clinvar confidence filter
    cln_rs_qdict = {"practice guideline":4, \
                    "reviewed by expert panel":3, \
                    "criteria provided, multiple submitters, no conflicts": 2, \
                    "criteria provided, conflicting interpretations": 1, \
                    "criteria provided, single submitter": 1, \
                    "no assertion for the single variant": 0, \
                    "no assertion criteria provided": 0, \
                    "no assertion provided": 0}
    
    cln_rs_qdict_supplement = { k.replace(" ", "_"): v for k,v in cln_rs_qdict.items() }
    cln_rs_qdict = {**cln_rs_qdict, **cln_rs_qdict_supplement}
    
    rs_bools_one_star = groupdf["shared_domain_CLINREVSTAT"].map(cln_rs_qdict) >= 1
    rs_bools_two_star = groupdf["shared_domain_CLINREVSTAT"].map(cln_rs_qdict) >= 2
    
    # One critical step is that we need to filter on the variant type.
    same_vartype_bools = groupdf.apply(judge_same_vartype_with_clinvar, axis=1)

    # First summary on clinvar sig
    sig_dict = groupdf.loc[(rs_bools_one_star & same_vartype_bools).to_numpy(), "shared_domain_CLNSIG"].value_counts().to_dict()
    sum_value = ", ".join([ str(k) + ":" + str(v) for k, v in sig_dict.items() ])
    sum_value = sum_value if len(sum_value) > 0 else np.nan
    groupdf[one_star_label] = sum_value
    
    sig_dict = groupdf.loc[(rs_bools_two_star & same_vartype_bools).to_numpy(), "shared_domain_CLNSIG"].value_counts().to_dict()
    sum_value = ", ".join([ str(k) + ":" + str(v) for k, v in sig_dict.items() ])
    sum_value = sum_value if len(sum_value) > 0 else np.nan
    groupdf[two_star_label] = sum_value
    
    return groupdf


@deal_with_nullkey_group_return
def summary_AA_change_per_variant(groupdf, 
                                  new_label_one_star="same_residue_CLNSIG_C1", 
                                  new_label_two_star="same_residue_CLNSIG_C2"):
    same_residue_bools = (groupdf["shared_domain_CLIN_same_AA"] == "Yes") | \
                         (groupdf["shared_domain_CLIN_same_AA"] == "No") | \
                         (groupdf["shared_domain_CLIN_same_AA"] == "Similar")
    same_residue_df = groupdf.loc[same_residue_bools, :]
    if len(same_residue_df) == 0:
        for value in ["Same_AA_change", "Similar_AA_change", "Diff_AA_change"]:
            groupdf[value+"_"+new_label_one_star] = np.nan
            groupdf[value+"_"+new_label_two_star] = np.nan
    else:
        same_residue_df.loc[:, "shared_domain_CLIN_same_AA"] = same_residue_df["shared_domain_CLIN_same_AA"].map({"Yes": "Same_AA_change", "Similar": "Similar_AA_change", "No": "Diff_AA_change"}).fillna(same_residue_df["shared_domain_CLIN_same_AA"])
        one_star_sum_value = ""
        for value in ["Same_AA_change", "Similar_AA_change", "Diff_AA_change"]:
            value_df = summary_clinvar_per_domain(same_residue_df.loc[same_residue_df["shared_domain_CLIN_same_AA"] == value,:], new_label_one_star, new_label_two_star)
            if len(value_df[new_label_one_star].tolist()) > 0:
                one_star = value_df[new_label_one_star].tolist()[0] if type(value_df[new_label_one_star].tolist()[0]) == str else np.nan
            else:
                one_star = np.nan
            groupdf[value+"_"+new_label_one_star] = one_star
            if len(value_df[new_label_two_star].tolist()) > 0:
                two_star = value_df[new_label_two_star].tolist()[0] if type(value_df[new_label_two_star].tolist()[0]) == str else np.nan
            else:
                two_star = np.nan
            groupdf[value+"_"+new_label_two_star] = two_star
            
    return groupdf
        

def judge_same_vartype_with_clinvar(row):
    if re.search(r'^[ATCG]{1}$', str(row["Ref"])) and re.search(r'^[ATCG]{1}$', str(row["Alt"])) and row["shared_domain_Clinvar_type"] == "single_nucleotide_variant":
        return True
    elif re.search(r'^[ATCG]+$', str(row["Ref"])) and re.search(r'^[ATCG]+$', str(row["Alt"])) and (row["shared_domain_Clinvar_type"] == "Indel" or row["shared_domain_Clinvar_type"] == "Microsatellite"):
        return True
    elif re.search(r'^<DEL>$', str(row["Alt"])) and (row["shared_domain_Clinvar_type"] == "Microsatellite" or row["shared_domain_Clinvar_type"] == "Deletion" or row["shared_domain_Clinvar_type"] == "copy number loss"):
        return True
    elif re.search(r'^<INS>$', str(row["Alt"])) and (row["shared_domain_Clinvar_type"] == "Microsatellite" or row["shared_domain_Clinvar_type"] == "Insertion"):
        return True
    elif re.search(r'^<DUP>$', str(row["Alt"])) and (row["shared_domain_Clinvar_type"] == "Microsatellite" or row["shared_domain_Clinvar_type"] == "copy number gain" or row["shared_domain_Clinvar_type"] == "Duplication"):
        return True
    elif re.search(r'^<DUP:TANDEM>$', str(row["Alt"])) and (row["shared_domain_Clinvar_type"] == "Microsatellite" or row["shared_domain_Clinvar_type"] == "Duplication"):
        return True
    elif re.search(r'^<INV>$', str(row["Alt"])) and row["shared_domain_Clinvar_type"] == "Inversion":
        return True
    else:
        return False


def check_df(df, df_name="", headers=[]):
    logger.info("Dataframe {} has {} rows and {} columns, looks like this:\n{}".format(df_name, df.shape[0], df.shape[1], df[:3].to_string()))
    if len(headers) > 0:
        new_headers = [ h for h in df.columns.tolist() if h not in headers ]
        new_cols = {}
        for col in new_headers:
            new_cols[col] = np.logical_not(df[col].apply(na_value)).sum()
        logger.info("Dataframe {} has some new headers, {}".format(df_name, ", ".join(["column " + k + " has " + str(v) + " valid records" for k,v in new_cols.items() ])))


def find_shared_domain_clinvar_rec(pat_table, clinvar_table, domain_region, output_table, clinvar_confidence=1):
    domain_region_remained_cols = ["feature_name", "note", "mis_obs_num", "lof_obs_num", "mis_exp_num", "lof_exp_num", "mis_oe_pvalue", "lof_oe_pvalue"]
    # domain_region_dropped_cols = ['chrom', 'chr_start', 'chr_end', 'strand', 'gene', 'protein_ID', 'gene_ID', 'type', 'chromosome', 'start_position', 'end_position', 'cds_length', 'num_coding_exons', 'canonical', 'exp_mis', 'mu_mis', 'possible_mis', 'exp_syn', 'mu_syn']
    dr_df = pd.read_csv(domain_region, sep='\t', low_memory=False).dropna(how="all")
    domain_region_dropped_1 = [ c for c in dr_df.columns.tolist() if c not in domain_region_remained_cols and c not in ["chrom", "chr_start", "chr_end", "strand", "gene"] ]
    domain_region_dropped_2 = ["chrom", "chr_start", "chr_end", "strand", "gene"]
    dr_df.drop(columns=domain_region_dropped_1, inplace=True, errors="ignore")
    
    clinvar_table_dropped_1 = ['GeneID', 'GeneSymbol', 'HGNC_ID', 'ClinSigSimple', 'LastEvaluated', 'PhenotypeIDS', 'Origin', 'OriginSimple', 'Assembly', 'ChromosomeAccession', 'ReferenceAllele', 'AlternateAllele', 'Cytogenetic', 'NumberSubmitters', 'Guidelines', 'TestedInGTR', 'OtherIDs', 'SubmitterCategories', 'VariationID', 'PositionVCF']
    clinvar_table_dropped_2 = ['Chr', 'Start', 'End', 'Ref', 'Alt']
    ct_df = pd.read_table(clinvar_table, low_memory=False, nrows=1000).dropna(how="all").drop_duplicates() # Just need header info
    ct_df = ct_df.drop(columns=clinvar_table_dropped_1, errors="ignore").drop_duplicates()
    
    dr_df = dr_df.loc[dr_df.iloc[:, 2].astype(str).str.contains(r"^[0-9]+$"), :].drop_duplicates()
    ct_df = ct_df.loc[ct_df.iloc[:, 2].astype(str).str.contains(r"^[0-9]+$"), :].drop_duplicates()
    
    logger.info("Dataframe {} has {} rows and {} columns, looks like this:\n{}".format("domain_table", dr_df.shape[0], dr_df.shape[1], dr_df[:3].to_string()))
    logger.info("Dataframe {} has {} rows and {} columns, looks like this:\n{}".format("Clinvar_table", ct_df.shape[0], ct_df.shape[1], ct_df[:3].to_string()))
    
    rename_dict = {"feature_name":"Interpro_domain_name", \
                   "note":"Interpro_domain_note", \
                   "Type":"shared_domain_Clinvar_type", \
                   "Name":"shared_domain_Clinvar_name", \
                   "RS# (dbSNP)":"shared_domain_CLIN_dbSNP_RS#", \
                   "nsv/esv (dbVar)":"shared_domain_CLIN_dbVar_nsv_esv", \
                   "RCVaccession":"shared_domain_CLINRCVaccession", \
                   "ClinicalSignificance":"shared_domain_CLNSIG", \
                   "PhenotypeList":"shared_domain_CLINDN", \
                   "ReviewStatus":"shared_domain_CLINREVSTAT", \
                   "mis_obs_num": "Interpro_domain_mis_obs_num", \
                   "lof_obs_num": "Interpro_domain_lof_obs_num", \
                   "mis_exp_num": "Interpro_domain_mis_exp_num", \
                   "lof_exp_num": "Interpro_domain_lof_exp_num", \
                   "mis_oe_pvalue": "Interpro_domain_mis_oe_pvalue", \
                   "lof_oe_pvalue": "Interpro_domain_lof_oe_pvalue"}
    
    raw_pat_df = pd.read_csv(pat_table, sep="\t", low_memory=False, na_values=[".", "-"]).drop(columns=[v for k,v in rename_dict.items()], errors="ignore")
    pat_df = raw_pat_df.loc[:, ["Chr", "Start", "End", "Ref", "Alt", "Gene.refGene", "AAChange.refGene", "uniq_ID"]].drop_duplicates()
    pat_df = pat_df.loc[pat_df.iloc[:, 2].astype(str).str.contains(r"^[0-9]+$") & pat_df.iloc[:, 1].astype(str).str.contains(r"^[0-9]+$"), :]
    check_df(pat_df, "Patient_record_table")
    
    merged_table = find_shared_rec(pat_df, ct_df, dr_df, mid_drop_cols=domain_region_dropped_2, right_drop_cols=clinvar_table_dropped_2).dropna(subset=["Name"])
    merged_table = merged_table.rename(columns=rename_dict)
    logger.info("The merged table, before finding same AA change with known CLINVAR variants, has shape of {}. And it looks like:\n{}".format(merged_table.shape, merged_table[:10].to_string(index=False)))
    # Use CLN_shared_domain to detect whether there are similar AA change
    results = merged_table.apply(find_same_AA_change, axis=1).to_numpy()
    headers = merged_table.columns.tolist()
    merged_table = merged_table.assign(shared_domain_CLIN_same_AA = results).dropna(subset=["Interpro_domain_name"], how="all")
    check_df(merged_table, "patient_records_matched_with_clinvar_and_domain_info", headers)
    
    merged_table = merged_table.groupby("uniq_ID", as_index=False, dropna=False).apply(summary_AA_change_per_variant).drop_duplicates().dropna(subset=["shared_domain_CLINREVSTAT"])
    logger.info("Here is the merged table shape {} after summary the same AA change evidences. IT looks like :\n{}\n".format(merged_table.shape, merged_table[:10].to_string(index=False)))
    
    # After identifying similar AA change, we need to groupby domain and summary the clinvar records.
    merged_table = merged_table.groupby(["Gene.refGene", "Interpro_domain_name"], as_index=False, dropna=False).apply(summary_clinvar_per_domain).drop_duplicates()
    check_df(merged_table, "patient_records_matched_with_clinvar_and_domain_info")
    
    # Remove redundant columns
    merged_table = merged_table.drop(columns=["shared_domain_CLIN_same_AA", "shared_domain_Clinvar_type", "shared_domain_Clinvar_name", "shared_domain_CLIN_dbSNP_RS#", "shared_domain_CLIN_dbVar_nsv_esv", "shared_domain_CLINRCVaccession", "shared_domain_CLNSIG", "shared_domain_CLINDN", "shared_domain_CLINREVSTAT"]).drop_duplicates()
    check_df(merged_table, "patient_records_matched_with_clinvar_and_domain_info")
    target_labels = [ l for l in merged_table.columns.tolist() if re.search(r'Interpro_domain|_CLNSIG_C', l) ]
    logger.info("The target labels are: {}".format(target_labels))
    
    merged_pat_df = raw_pat_df.merge(merged_table.loc[:, target_labels + ["uniq_ID"]].drop_duplicates().dropna(subset=["uniq_ID"]), how="left", on="uniq_ID")
    id_post1_loc = merged_pat_df.columns.get_loc("CLNSIG") + 1
    target_labels.reverse()
    for l in target_labels:
        tobe_moved = merged_pat_df.pop(l)
        merged_pat_df.insert(loc=id_post1_loc, column=l, value=tobe_moved)
        logger.info("After moving the column {} to {}th column, the table looks like:\n{}".format(l, id_post1_loc, merged_pat_df[:4].to_string()))

    merged_pat_df.drop_duplicates().to_csv(output_table, header=True, sep='\t', index=False, encoding="utf-8")


if __name__ == "__main__":
    
    parser = ap.ArgumentParser()
    parser.add_argument("-pt", "--pat_table", type=str, help="Path to the left table", required=True)
    parser.add_argument("-ct", "--clinvar_table", type=str, help="Path to the right table", required=False, default="/paedyl01/disk1/yangyxt/public_data/clinvar/tab_delimited/Reformat_updated_Clinvar.tsv")
    parser.add_argument("-dr", "--domain_region", type=str, help="Path to the mid table", required=False, default="/paedyl01/disk1/yangyxt/public_data/Prot2HG/hg19_prot2hg_summary.gnomAD_v2_exome_overlap.count_obs.constraint_metrics.tsv")
    parser.add_argument("-ot", "--output_table", type=str, help="Path to the output table", required=True)
    parser.add_argument("-cc", "--clinvar_confidence", type=int, help="Threshold for quantified clinvar confidence", required=False, default=1)

    args = parser.parse_args()
    find_shared_domain_clinvar_rec(os.path.abspath(args.pat_table), os.path.abspath(args.clinvar_table), os.path.abspath(args.domain_region), args.output_table)








 
    
    
