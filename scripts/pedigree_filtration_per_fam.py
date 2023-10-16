import pandas as pd
import argparse as ap
import os
import sys
import numpy as np
import re
import logging


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


# Check if list contain element and return boolean
def if_contains(input, check_list):
    if input in check_list:
        return True
    else:
        return False

vec_contain = np.vectorize(if_contains)

# Check if two object equals and return booleans
def if_equal(left, right):
    if left==right:
        return True
    else:
        return False

vec_equal = np.vectorize(if_equal)

# Fetch the GT column by ID, return header str
def find_gt_header(ID, header):
    regex = re.compile("^" + str(ID) + "_.+_GT$")
    if re.match(regex, header):
        return header

vec_findgt = np.vectorize(find_gt_header)

# Filter on the homozygous ALT allele
def filter_homoalt(gt_info, *alleles):
    logger.debug("This individual's GT info is {}".format(gt_info) + ". While the ALT alleles proband have is:" + str(alleles))
    if str(gt_info).isnumeric():
        if float(gt_info) == 2:
            return False
        else:
            return True
    try:
        if 2 in [gt_info.split(gt_info[1]).count(str(x)) for x in alleles]:
            logger.debug("This healthy individual contains two copies of either ALT alleles")
            return False
        else:
            logger.debug("This healthy individual contains only one copy of any ALT alleles")
            return True
    except IndexError:
        logger.debug("Splitting GT info by / or | failed because the individual does not have GT info on this site, So let's check GT info: {}".format(gt_info))
        return True
    

# Filter on the homozygous REF allele
def filter_homoref(gt_info):
    logger.debug("This individual's GT info is {}".format(gt_info))
    if str(gt_info).isnumeric():
        if float(gt_info) == 2:
            return False
        else:
            return True
    try:
        if gt_info.split(gt_info[1]).count("0") == 2:
            logger.debug("This sick individual contains two copies of REF alleles")
            return False
        else:
            logger.debug("This sick individual contains ALT allele on this site")
            return True
    except IndexError:
        logger.debug("Splitting GT info by / or | failed because the individual does not have GT info on this site. So let's check GT info: {}".format(gt_info))
        return False
    

# Define a function that returns a boolean array based on one individual's GT information
def filter_by_person(row, mem_gt, pro_gt, ped_table, motherID, fatherID):
    mem_ID = mem_gt.split('_')[0]
    logger.debug("Family member's ID is:" + mem_ID)
    pro_ID = pro_gt.split('_')[0]
    logger.debug("Family proband's ID is:" + pro_ID)
    mem_pheno = int(ped_table.loc[ped_table['IndividualID'] == mem_ID, 'Phenotype'])
    logger.debug("Family member's phenotype is:" + str(mem_pheno))
    mem_gender = int(ped_table.loc[ped_table['IndividualID'] == mem_ID, 'Sex'])
    logger.debug("Family member's gender is:" + str(mem_gender))
    pro_gender = int(ped_table.loc[ped_table['IndividualID'] == pro_ID, 'Sex'])
    logger.debug("Family proband's gender is:" + str(pro_gender))
    
    sex_chr = ["chrX", "chrY"]
    try:
        alt_alleles = [ str(x) for x in str(row[pro_gt]).split(row[pro_gt][1]) if str(x) != "0" ]
        logger.debug("Family proband's Alt Allele on this site is:" + str(alt_alleles))
    except (IndexError, TypeError) as e:
        if str(row[pro_gt]).isnumeric():
            logger.debug("Running into error {}. Splitting GT info by / or | failed because the proband's GT info is just an integer coming from CNV record on this site. So let's check GT info: {}".format(e, row[pro_gt]))
            if str(row[pro_gt]) == "0":
                logger.warning("This proband GT info is from a CNV record and on this site, the proband has copy number of 2. Thus discard this record:\n{}\n".format(row))
                return False
            elif str(row[pro_gt]) == "1":
                alt_alleles = [ "1" ]
            elif str(row[pro_gt]) == "2":
                alt_alleles = [ "1" ]
        else:
            logger.warning("Running into error {}. Splitting GT info by / or | failed because the proband does not have GT info on this site. So let's check GT info: {}".format(e, row[pro_gt]))
            return False
    
    if mem_ID == motherID:
        if mem_pheno == 1 and row['Chr'] != "chrY":
            return filter_homoalt(row[mem_gt], *alt_alleles)
        elif mem_pheno == 2 and row['Chr'] != "chrY":
            return filter_homoref(row[mem_gt])
        else:
            # Now handling variants on chrY
            return True
        
    elif mem_ID == fatherID:
        if mem_pheno == 1:
            return filter_homoalt(row[mem_gt], *alt_alleles)
        elif mem_pheno == 2:
            return filter_homoref(row[mem_gt])
    else:
        # We are filtering on sibs now
        if mem_pheno == 1:
            return filter_homoalt(row[mem_gt], *alt_alleles)
        elif mem_gender == 2 and pro_gender == 1 and mem_pheno == 2 and row['Chr'] != 'chrY':
            return filter_homoref(row[mem_gt])
        elif mem_pheno == 2:
            return filter_homoref(row[mem_gt])   
    

# Define main_program
def pedigree_filter(var_path="", ped_path="", target_fam="", output=""):
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(name)s:%(exc_info)s:%(message)s")
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    ped_table = pd.read_csv(ped_path, sep='\t', low_memory=False)
    var_table = pd.read_csv(var_path, sep='\t', low_memory=False, nrows=1)
    
    fam_ped = ped_table.loc[ped_table['#FamilyID'] == target_fam, :]
    logger.info(str(fam_ped))
    # Check how many persons in fam_mem_no
    fam_mem_no = fam_ped.shape[0]
    logger.info("How many persons in this family:" + str(fam_mem_no))
    
    if fam_mem_no > 1:
        # Find patient's ID and parents' ID
        patient_IDs = fam_ped.loc[fam_ped['Phenotype'] == 2, 'IndividualID'].to_list()
        motherID = fam_ped.loc[fam_ped['IndividualID']==patient_IDs[0], 'MaternalID'].squeeze()
        fatherID = fam_ped.loc[fam_ped['IndividualID']==patient_IDs[0], 'PaternalID'].squeeze()
        
        # Find Sibling's IDs
        bool_array = vec_equal(fam_ped['IndividualID'], patient_IDs[0])
        logger.debug(str(bool_array) + str(type(bool_array)))
        beside_pro = fam_ped.loc[np.logical_not(bool_array), :]
        logger.debug(str(beside_pro))
        
        fa_array = vec_equal(beside_pro['PaternalID'], fatherID)
        mo_array = vec_equal(beside_pro['MaternalID'], motherID)
        logger.debug(str(fa_array) + str(type(fa_array)))
        logger.debug(str(mo_array) + str(type(mo_array)))

        sib_IDs = beside_pro.loc[fa_array | mo_array, 'IndividualID'].tolist()
        
        # Fetch the column header of each fam's GT info. 
        proband_gt = [x for x in var_table.columns.tolist() if re.search(re.compile("^" + str(patient_IDs[0]) + "_.+_GT$"), x)]
        mother_gt = [x for x in var_table.columns.tolist() if re.search(re.compile("^" + str(motherID) + "_.+_GT$"), x)]
        father_gt = [x for x in var_table.columns.tolist() if re.search(re.compile("^" + str(fatherID) + "_.+_GT$"), x)]
        sib_gts = [(lambda x: [u for u in vec_findgt(x, var_table.columns.tolist()) if u != None][0])(x) for x in sib_IDs]
        logger.info("The proband's GT column header is:" + proband_gt[0])
        try:
            logger.info("Mother's GT column header is:" + mother_gt[0])
        except IndexError:
            logger.info("The proband do hot have his/her mother's DNA")
        try:
            logger.info("Father's GT column header is:" + father_gt[0])
        except IndexError:
            logger.info("The proband do hot have his/her father's DNA")
        try:
            logger.info("sibs' GT column headers are:" + str(sib_gts))
        except IndexError:
            logger.info(("The proband do hot have his/her sib's DNA"))
        
        # Combine gt headers to one list 
        beside_pro_gts = mother_gt + father_gt + sib_gts
        logger.info(str(beside_pro_gts) + str(type(beside_pro_gts)))
        
        # Merge with count series.
        var_df = pd.read_csv(var_path, sep='\t', low_memory=False, encoding="utf-8", chunksize=100000)
        pd.read_csv(var_path, sep='\t', low_memory=False, encoding="utf-8", nrows=0).to_csv(output, sep='\t', index=False, encoding="utf-8")
        while True:
            try:
                chunk = next(var_df)
                logger.info("For this 100000 lines, we have the chunk df looks like:\n{}".format(chunk[:4].to_string(index=False)))
                regex = re.compile("^[Cc][Hh][Rr]")
                if not re.match(regex, str(chunk['Chr'].tolist()[0])):
                    chunk['backup'] = 'chr'
                    chunk['Chr'] = chunk['backup'] + chunk['Chr'].astype('str')
                    chunk.drop(columns=['backup'], inplace=True)
                    
                if len(beside_pro_gts) > 0:
                    iter_array = chunk['Chr'] != None
                    logger.info(str(iter_array) + str(type(iter_array)))
                    for beside_pro_gt in beside_pro_gts:
                        bool_array = chunk.apply(filter_by_person, args=(beside_pro_gt, proband_gt[0], fam_ped, motherID, fatherID), axis=1)
                        iter_array = np.logical_and(iter_array, bool_array)
                    logger.info(str(iter_array) + str(type(iter_array)))
                    logger.info(str(iter_array.isnull().sum()))
                    output_table = chunk.loc[iter_array, :]
                    output_table.to_csv(output, sep='\t', index=False, mode="a", encoding="utf-8", header=False)
                else:
                    chunk.to_csv(output, sep='\t', index=False, mode="a", encoding="utf-8", header=False)
            except StopIteration:
                break         
    else:
        # We only have the proband
        logger.info("We have only the proband in this fam")
        var_df = pd.read_csv(var_path, sep='\t', low_memory=False, encoding="utf-8", chunksize=100000)
        pd.read_csv(var_path, sep='\t', low_memory=False, encoding="utf-8", nrows=0).to_csv(output, sep='\t', index=False, encoding="utf-8")
        while True:
            try:
                chunk = next(var_df)
                logger.info("For this 100000 lines, we have the chunk df looks like:\n{}".format(chunk[:4].to_string(index=False)))
                regex = re.compile("^chr")
                if not re.match(regex, str(chunk['Chr'].tolist()[0])):
                    chunk['backup'] = 'chr'
                    chunk['Chr'] = chunk['backup'] + chunk['Chr'].astype('object')
                    chunk.drop(columns=['backup'], inplace=True)
                chunk.to_csv(output, sep='\t', index=False, mode="a", encoding="utf-8", header=False)
            except StopIteration:
                break


if __name__ == '__main__':
    # args.var_table should be the path to the table of variants.
    # args.ped_table should be the path to the table of pedigree info.
    # args.target_fam should be the name of family ID.
    # args.output should be the path to the output table.
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    
    parser = ap.ArgumentParser()
    parser.add_argument("-vt", "--var_table", type=str, help="The path to the table of variants.", required=True)
    parser.add_argument("-pt", "--ped_table", type=str, help="The path to the table of pedigree info", required=True)
    parser.add_argument("-tf", "--target_fam", type=str, help="The Family ID of the target fam in this batch", required=True)
    parser.add_argument("-o", "--output", type=str, help="The path to the output table", required=True)
    
    args=parser.parse_args()
    pedigree_filter(var_path=args.var_table,ped_path=args.ped_table,target_fam=args.target_fam,output=args.output)
    
    
