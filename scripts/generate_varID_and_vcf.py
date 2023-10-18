import pandas as pd 
import numpy as np 
import re 
import os
import logging
import hashlib
import sys
import argparse as ap
from check_qname_count import id_generator
from datetime import date
import subprocess
from python_utils import deal_with_nullkey_group_return


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def Insert_row_(row_number, df, row_value): 
    # Slice the upper half of the dataframe 
    df1 = df[0:row_number] 
    # Store the result of lower half of the dataframe 
    df2 = df[row_number:] 
    # Inser the row in the upper half dataframe 
    df1.loc[row_number]=row_value 
    # Concat the two dataframes 
    df_result = pd.concat([df1, df2]) 
    # Reassign the index labels 
    df_result.index = [*range(df_result.shape[0])] 
    # Return the updated dataframe 
    return df_result


@deal_with_nullkey_group_return
def prepare_ref_alt_within_groupdf(group_df, alt_alleles):
    # logger.info("By_coordinates grouped df is "+ str(group_df.iloc[:, :5]))
    var_len = int(float(group_df['End'].tolist()[0]) - float(group_df['Start'].tolist()[0])) + 1 if group_df["Alt"].tolist()[0] != '<INS>' else re.search(r'SVLEN=([0-9]+);', group_df["INFO"].tolist()[0]).group(1)
    # logger.info(var_len)
    ref_allele = group_df['Ref'].tolist()[0]
    variants = ["+" + alt_allele.replace(ref_allele, '', 1) if len(alt_allele) >= len(ref_allele) else "-" + ref_allele.replace(alt_allele, '', 1) for alt_allele in alt_alleles]
    sub_index = [alt_allele.find(ref_allele) if len(alt_allele) >= len(ref_allele) else ref_allele.find(alt_allele) for alt_allele in alt_alleles]
    # logger.info(variants)
    kicking_out_var_ind = []
    for ind, var in enumerate(variants):
        if sub_index[ind] != -1:
            if var_len == 1:
                if "-" in var and len(var) - 1 != 1:
                    logger.info("The variant {} length does not fit with the start_end interval length {}.".format(str(var), str(var_len)))
                    kicking_out_var_ind.append(ind)
                elif "+" in var:
                    pass
                else:
                    logger.info("The variant {} length fit with the start_end interval length {}.".format(str(var), str(var_len)))
            elif var_len > 1:
                if "-" in var and len(var) - 1 != var_len:
                    logger.info("The variant {} length does not fit with the start_end interval length {}.".format(str(var), str(var_len)))
                    kicking_out_var_ind.append(ind)
                elif "+" in var:
                    kicking_out_var_ind.append(ind)
                else:
                    logger.info("The variant {} length fit with the start_end interval length {}.".format(str(var), str(var_len)))
        else:
            logger.warning("Seems both Ref and Alt alleles are not part of each other. Therefore we don't know the precise variant.")
            '''
            if abs(len(alt_alleles[ind]) - len(ref_allele)) == var_len:
                logger.info("The variant {} length fit with the start_end interval length {}.".format(str(var), str(var_len)))
            else:
                logger.info("The variant {} length does not fit with the start_end interval length {}.".format(str(var), str(var_len)))
                kicking_out_var_ind.append(ind)
            '''
            kicking_out_var_ind.append(ind)
                
    remaining_alt_alleles = [alt for ind, alt in enumerate(alt_alleles) if ind not in kicking_out_var_ind]
    if len(remaining_alt_alleles) == 0:
        return
    # logger.info(remaining_alt_alleles)
    if len(group_df) > len(remaining_alt_alleles):
    # The group_df row number is larger than Alt alleles.
        logger.info("The group_df is like \n{}, \nand the remaining alt alleles are like {}".format(group_df.iloc[:5, :5].to_string(index=False), str(remaining_alt_alleles)))
        if len(group_df) >= 2:
            while len(group_df) > len(remaining_alt_alleles):
                group_df = group_df.iloc[:-1, :]
            group_df['Alt'] = np.array(remaining_alt_alleles)
            logger.info("The returned group_df looks like: "+ group_df.iloc[:5, :5].to_string(index=False))
            return group_df
        else:
            group_df['Alt'] = [alt_alleles][0]
            logger.info("After altering the Alt allele,"+ group_df.iloc[:5, :5].to_string(index=False))
            return group_df
    else:
        logger.info('The type of Alt alleles is out number the df row number. Why?')
        logger.info("Let's take a look at the group_df \n{}\n and remaining alt alleles {}.".format(group_df.iloc[:5, :5].to_string(index=False), str(remaining_alt_alleles)))
        while len(remaining_alt_alleles) > len(group_df):
            group_df = Insert_row_(1, group_df, group_df.iloc[0, :])
        group_df['Alt'] = np.array(remaining_alt_alleles)
        logger.info("The returned group_df looks like: \n"+ group_df.iloc[:5, :5].to_string(index=False))
    return group_df


@deal_with_nullkey_group_return
def input_vep_start_end_ref_alt(group_df):
    '''
    This function is only to clean the format for records with multiple ALTs
    '''
    if len(group_df['Ref'].tolist()[0]) == 1 and group_df['Alt'].str.contains(",").any():
        alt_alleles = group_df['Alt'].tolist()[0].split(",")
        
        # Remove the * if any
        try:
            alt_alleles.remove("*")
        except ValueError:
            pass

        # Check whether some records are removed before due to duplicated rows. If removed, add them back    
        while len(alt_alleles) > len(group_df):
            group_df = Insert_row_(1, group_df, group_df.iloc[0, :])
        
        # logger.info(group_df)
        by_end = group_df.groupby(['Start','End'], as_index=False, sort=False, dropna=False)
        # logger.info(by_end.groups)
        group_df = by_end.apply(prepare_ref_alt_within_groupdf, alt_alleles).reset_index(drop=True, level=0)
        logger.info("\nThis function has been applied to one group once and the group_df it received is: \n" + group_df.iloc[:5, :5].to_string(index=False))
        return group_df
    elif len(group_df['Ref'].tolist()[0]) > 1 and group_df['Alt'].str.contains(",").any():
        alt_alleles = group_df['Alt'].tolist()[0].split(",")
        
        # Check whether some records are removed before due to duplicated rows. If removed, add them back    
        while len(alt_alleles) > len(group_df):
            group_df = Insert_row_(1, group_df, group_df.iloc[0, :])
        
        # Assign the Alt alleles to the group_df Alt column
        logger.info("Ref allele is not a single base: " + str(alt_alleles))
        logger.info("Ref allele is not a single base: \n" + group_df.iloc[:5, :5].to_string(index=False))
        by_end = group_df.groupby(['Start','End'], as_index=False, sort=False, dropna=False)
        # logger.info(by_end.groups)
        group_df = by_end.apply(prepare_ref_alt_within_groupdf, alt_alleles).reset_index(drop=True, level=0)
        logger.info("The variant happens on the same segment, but the coordinates can be slightly different, we directly assign each alt allele to each record and have this returned: \n"+ group_df.iloc[:5, :5].to_string(index=False))
        return group_df
    elif group_df['Alt'].str.contains("<").any():
        # logger.info("This group of variants are all CNVs, so do not change them."+ str(group_df.iloc[:, :5]))
        return group_df
    else:
        # isplay("This group of variants are all normal variants, do not change them."+ str(group_df.iloc[:, :5]))
        return group_df
    

def prepare_vep_input(row):
    '''
    VCF 4.2 format has 8 fixed, mandatory columns:
    #CHROM POS ID REF ALT QUAL FILTER INFO
    '''
    if re.search("^<INS>$", str(row["Alt"])):
        info_str = str(row["INFO"])
        info_fields = info_str.split(";")
        sv_lens = [ int(field_value.split(":")[-1]) for field_value in info_fields if re.search(r"^SVLEN=", field_value) ]
        if len(sv_lens) > 0:
        	var_str = str(row["Chr"] + ":" + str(row["Start"]).replace(".0", "") + "-" + str(row["End"]).replace(".0","") + ":" + str(row["Ref"]) + "-->" + str(row["Alt"])) + ":" + str(sv_lens[0])
        else:
            var_str = str(row["Chr"] + ":" + str(row["Start"]).replace(".0", "") + "-" + str(row["End"]).replace(".0","") + ":" + str(row["Ref"]) + "-->" + str(row["Alt"]))
    else:
    	var_str = str(row["Chr"] + ":" + str(row["Start"]).replace(".0", "") + "-" + str(row["End"]).replace(".0","") + ":" + str(row["Ref"]) + "-->" + str(row["Alt"]))
    var_bytes = var_str.encode("utf-8")
    uniq_id = hashlib.sha256(var_bytes).hexdigest()
    if re.search("^<[A-Z]+>$", str(row["Alt"])):
        '''
        Now we are dealing with SV records
        '''
        if re.search("DUP:TANDEM", str(row['Alt'])) or re.search("SVTYPE=DUP:TANDEM", str(row['INFO'])):
            return row['Chr'], int(float(row['Start'])), uniq_id, ".", row['Alt'], ".", "PASS", "SVTYPE=TDUP;END=" + str(row['End']), uniq_id
        else:
            return row['Chr'], int(float(row['Start'])), uniq_id, ".", row['Alt'], ".", "PASS", "SVTYPE=" + row['Alt'].strip("<").strip(">") + ";END=" + str(row['End']), uniq_id
    else:
        '''
        Now we are dealing with SNV and short Indels
        '''
        row["End"] = row["End"] if float(row["End"]) >= float(row["Start"]) else row["Start"]
        if len(row['Alt']) > len(row['Ref']):
            '''
            This record is an insertion variant
            '''
            return row['Chr'], int(float(row['Start'])), uniq_id, row['Ref'], row['Alt'], ".", "PASS", ".", uniq_id
        elif len(row['Alt']) < len(row['Ref']):
            '''
            This record is a deletion variant
            '''
            return row['Chr'], int(float(row['Start'])), uniq_id, row['Ref'], row['Alt'], ".", "PASS", ".", uniq_id
        else:
            '''
            This record is just a normal SNV
            '''
            return row['Chr'], int(float(row['Start'])), uniq_id, row['Ref'], row['Alt'], ".", "PASS", ".", uniq_id




@deal_with_nullkey_group_return
def choose_first_uid(group_df, col_label="ID"):
    group_df[col_label] = group_df[col_label].tolist()[0]
    return group_df




def main_generate_vep_input(var_table, 
                            ref_genome = "/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta",
                            ref_assembly_ver = "hg19"):
    '''
    Should convert the format to VCF 4.2 format
    Because only VCF format can contain Structural variants records. 
    Also, VCF 4.2 format allows multiple records with the same POS value. This is compatible with the sites having multiple ALT alleles.
    VCF 4.2 format has 8 fixed, mandatory columns:
    #CHROM POS ID REF ALT QUAL FILTER INFO
    '''
    if type(var_table) == str:
        var_table_path = var_table
        var_table = pd.read_table(var_table_path, low_memory=False).drop(columns=["uniq_ID"], errors="ignore")
    elif type(var_table) == pd.core.frame.DataFrame:
        pass
    
    var_table_headers = var_table.columns.tolist()
    
    # Adjust the Alt column to make the whole VCF records bi-allelic
    if var_table['Alt'].str.contains(",", regex=False).any():
        by_var_group = var_table.groupby(['Chr', 'Start', 'End', 'Ref', 'Alt', 'INFO'], as_index=False, sort=False, dropna=False)
        var_table_new_refalt = by_var_group.apply(input_vep_start_end_ref_alt).reset_index(drop=True)[var_table_headers].copy()
    else:
        var_table_new_refalt = var_table.copy(deep=True)
    
    var_table_new_refalt.dropna(inplace=True, how='all')
    logger.info("That's how new vartable looks like: \n{}".format(var_table_new_refalt[:5].to_string(index=False)))
    try:
        var_table_new_refalt.loc[:, 'Start'] = var_table_new_refalt['Start'].astype(np.longdouble).astype(int)
        var_table_new_refalt.loc[:,'End'] = var_table_new_refalt['End'].astype(np.longdouble).astype(int)
    except ValueError:
        bools = var_table_new_refalt.loc[:,'Start'].astype(str).str.contains(r'^[0-9]+$')
        logger.info("These rows can't be converted to float in Start column: \n{}".format(var_table_new_refalt.loc[np.logical_not(bools), :]))
        var_table_new_refalt = var_table_new_refalt.loc[bools, :]
        var_table_new_refalt.loc[:, 'Start'] = var_table_new_refalt['Start'].astype(np.longdouble).astype(int)
        var_table_new_refalt.loc[:,'End'] = var_table_new_refalt['End'].astype(np.longdouble).astype(int)
    # var_table_new_refalt.to_csv("test.new_ref_alt.txt", sep='\t', index=False)
    
    var_table_new_refalt = var_table_new_refalt.drop_duplicates(inplace=False)
    vep_input = pd.DataFrame()
    vep_input['#CHROM'], vep_input['POS'], vep_input['ID'], vep_input['REF'], vep_input['ALT'], vep_input['QUAL'], vep_input['FILTER'], vep_input['INFO'], var_table_new_refalt.loc[:,'uniq_ID'] = zip(*var_table_new_refalt.apply(prepare_vep_input, axis=1))
    # vep_input = vep_input.groupby(['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'INFO'], as_index=False, dropna=False).apply(choose_first_uid)
    # var_table_new_refalt = var_table_new_refalt.groupby(['Chr', 'Start', 'End', 'Ref', 'Alt'], as_index=False, dropna=False).apply(choose_first_uid, col_label="uniq_ID")
    for label in ['#CHROM', 'POS']:
        try:
            vep_input[label] = vep_input[label].astype('float').astype('int')
        except ValueError:
            vep_input[label] = vep_input[label].astype('str')
    vep_input.sort_values(by=['#CHROM','POS'], inplace=True)
    vep_input.loc[:, "REF"] = vep_input["REF"].replace({".": "N"})
    vep_input.drop_duplicates(inplace=True)
    
    try:
        logger.warning("Updating {} with new ref and alt and variant ID column and now it looks like : \n{}".format(var_table_path, var_table_new_refalt.iloc[:5, :].to_string(index=False)))
        var_table_new_refalt.to_csv(var_table_path, sep='\t', index=False)
        
        with open(var_table_path[:-4] + ".vcf", mode='w') as vcf:
            vcf.write("##fileformat=VCFv4.2\n")
            vcf.write("##fileDate=" + date.today().strftime("%d/%m/%Y") + "\n")
            vcf.write(f"##reference={ref_assembly_ver}\n")
            vcf.write(f"##assembly={ref_genome}\n")
            vcf.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
            vcf.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
            vcf.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
            vcf.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
            vcf.write("##ALT=<ID=BND,Description=\"Breakend\">\n")
            vcf.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
            vcf.write("##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">\n")
            vcf.write("##ALT=<ID=INS,Description=\"Insertion of novel sequence\">\n")
            vcf.write("##ALT=<ID=INV,Description=\"Inversion\">\n")
            
        vep_input.to_csv(var_table_path[:-4] + ".vcf", sep='\t', index=False, header=True, mode="a")
        cmd = "bash /paedyl01/disk1/yangyxt/ngs_scripts/common_bash_utils.sh \
               insert_contig_header_by_genome {vcf} {ref} && \
               bcftools sort --temp-dir /paedyl01/disk1/yangyxt/test_tmp {vcf} > {vcf}.tmp && \
               mv {vcf}.tmp {vcf} && ls -lh {vcf}".format(vcf = var_table_path[:-4] + ".vcf", ref=ref_genome)
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
        logger.info("\n{}".format(result.stdout))
    except NameError:
        return vep_input, var_table_new_refalt


if __name__ == '__main__':
    parser = ap.ArgumentParser()
    parser.add_argument("-vt", "--variant_table", type=str, help="The path to the table of variants in patients", required=True)
    parser.add_argument("-as", "--assembly", type=str, help="The assembly version of the reference genome", required=False, default="hg19")
    
    args=parser.parse_args()
    main_generate_vep_input(var_table=args.variant_table, 
                            ref_assembly_ver=args.assembly)
    


    
