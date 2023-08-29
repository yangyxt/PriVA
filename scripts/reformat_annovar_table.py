import sys
import os
import re
import pandas as pd
import numpy as np
import gc
import time
import logging
import argparse as ap
import subprocess
import io
import tempfile
from python_utils import deal_with_nullkey_group_return


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


# Def a function to be applied on every row to find the index.
@deal_with_nullkey_group_return
def find_index(group_df, fields):
    start_time = time.time()
    group_df.fillna(".", inplace=True)
    format_fields = group_df['Otherinfo12'].tolist()[0]
    field_list = format_fields.split(':')
    index_list = [ inde for inde, field_n in enumerate(fields) if field_n in field_list ]
    # print(index_list)
    # print(field_list)
    # print("building list: --- %s seconds ---" % (time.time() - start_time))
    start_time = time.time()
    index_list_np = np.asarray(index_list)
    index_list_np = index_list_np + 1 + index_format + sample_number
    # cdef int len_format
    len_format = len(format_fields_list)
    dict_list = []
    # print("set len format: --- %s seconds ---" % (time.time() - start_time))
    for n in range(0, sample_number):
        # start_time = time.time()
        final_ind = n + index_format + 1
        value_list = group_df.iloc[:, final_ind].str.split(':',expand=True)
        logger.debug("The targeted df looks like:\n{}\n".format(group_df.iloc[:5, (index_list_np + n * len_format)].to_string(index=False)))
        logger.debug("The value_list df is:\n{}\n".format(value_list[:5].to_string(index=False)))
        try:
            # start_time = time.time()
            group_df.iloc[:, index_list_np + n * len_format] = value_list.to_numpy()
            # print("try normal part: --- %s seconds ---" % (time.time() - start_time))
        except ValueError as ve:
            # Sometimes, we may find that the format content might not exactly corresponds to the format field names.
            # for vcf FORMAT field PS, we do not have a matched value.
            # start_time = time.time()
            logger.error("Running into this error: {}".format(ve))
            logger.warning("some columns are missing values hence cause the value_list:\n{}\n shorter than target columns {}".format(value_list, index_list_np + n * len_format))
            logger.warning("Take a look at this group_df: \n{}".format(group_df[:5].to_string(index=False)))
            if group_df.iloc[0,:].tolist()[:5] == group_df.columns.tolist()[:5]: 
                logger.warning("Seems the groupdf contains an odd row of headers. Quit returning this group.")
                return pd.DataFrame(columns=group_df.columns.to_list())
            value_list['PS'] = "."
            group_df.iloc[:, (index_list_np + n * len_format)] = value_list.to_numpy()
            # print("try except part: --- %s seconds ---" % (time.time() - start_time))
        except IndexError as ie:
            logger.error("Running into this error: {}".format(ie))
            logger.error("Why we run into an Index Error? Let's take a look at the final column indices used to assign values: \n{}\n{}".format(index_list_np + n * len_format, field_list))
            logger.error("It seems that this group_df is only composed of header lines. Abandon the group_df.")
            return pd.DataFrame(columns=group_df.columns.tolist(), dtype=np.dtype("O"))
    print("Total time for each iteration: --- %s seconds ---" % (time.time() - start_time), file=sys.stderr)
    return group_df.astype(np.dtype("O")).drop_duplicates()


def change_ref_alt(row):
    if "," in str(row['Alt.1']):
        if "<" in row['Alt.1']:
            return row['Ref'], row['Alt']
        elif "-" == row['Ref']:
        # Means this row is an insertion event
            logger.info("The Ref is {}, Alt is {}, Ref.1 is {}, Alt.1 is {}".format(row['Ref'], row['Alt'], row['Ref.1'], row['Alt.1']))
            '''
            Here is a case saying why it is important to replace Alt from Alt.1 to judge whether the rest part equals to Ref.1
            The Ref is -, Alt is CGACCCTGCCCTGGAGGCCC, Ref.1 is ACGCCCCTGCCCTGGAGGCCC, Alt.1 is ACGACCCTGCCCTGGAGGCCCCGCCCCTGCCCTGGAGGCCC,A,ACGCCCCTGCCCTGGAGGCCCCGCCCCTGCCCTGGAGGCCC
            '''
            return row['Ref.1'], [alt for alt in row['Alt.1'].split(",") if alt.replace(row['Alt'], '', 1) == row['Ref.1']][0]
        elif "-" == row['Alt']:
        # Means this row is an deletion event
            logger.info("The Ref is {}, Alt is {}, Ref.1 is {}, Alt.1 is {}".format(row['Ref'], row['Alt'], row['Ref.1'], row['Alt.1']))
            return row['Ref.1'], [alt for alt in row['Alt.1'].split(",") if row['Ref.1'].replace(row['Ref'], '', 1) == alt][0]
        else:
        # Means multiple possible SNV
            all_alts = row['Alt.1'].split(",")
            logger.info("The alt alleles are " + ",".join(all_alts))
            logger.info("The alt allele this row records is " + row['Alt'])
            try:
                all_alts.remove("*")
                logger.debug("The * is removed from alt alleles.")
            except ValueError:
                logger.warning("The alt alleles does not contain a *")
            try:
                return row['Ref'], [alt for alt in row['Alt.1'].split(",") if alt == row['Alt']][0]
            except IndexError:
                logger.warning("This site does not have a certain alt allele.")
                return row['Ref'], row['Alt']
    else:
        return row['Ref.1'], row['Alt.1']



def reformat(a, b, c, d, ped, output, fam):
    # a corresponds to the dir of hg19_multianno txt file and vcf file, which is the annotation dir.
    # b corresponds to the file_name of hg19_multianno txt file
    # c corresponds to the file_name of hg19_multianno vcf file
    # d corresponds to the dir where contains the ped file
    # ped corresponds to the ped file name, including the extension
    # output corresponds to the output tsv file name, including the extension
    # fam corresponds to the family name, usually it is $1 in bash scripts.

    os.chdir(a)

    # start_time = time.time()
    # Extract the first row as headers
    hg19_original_head = pd.read_csv(b, sep='\t', low_memory=False, nrows=2)
    headers = hg19_original_head.columns.tolist()
    logger.info("The headers of the original table are: \n"+ str(headers)) 

    # Extract the index of column Otherinfo
    index_Otherinfo = hg19_original_head.columns.get_loc('Otherinfo1')
    index_Format = hg19_original_head.columns.get_loc('Otherinfo12')
    # # print(index_Otherinfo)

    # Extract the single column of format_fields
    format_field = pd.read_table(b, low_memory=False, usecols=[index_Format]).squeeze("columns")
    gc.collect()
    format_field = format_field.dropna(how='all').tolist()
    global format_fields_list
    format_fields_list = []
    format_field = set(format_field)
    for field in format_field:
        format_fields_list.extend(field.split(':'))
    format_fields_list = set(format_fields_list)
    format_fields_list = list(format_fields_list)
    format_fields_list.sort()
    format_fields_list.remove('GT')
    format_fields_list.insert(0, 'GT')
    logger.info("In this table, the fields in FORMAT column are: " + str(format_fields_list))
    # # print('format_fields_list = ', format_fields_list)

    # Read the ped file to extract family info.
    
    ped_df = pd.read_csv(os.path.join(d, ped), sep='\t', engine='python')
    family_df = ped_df.loc[ped_df['#FamilyID'] == fam, :]
    logger.info("Family df looks like: \n{}".format(family_df.to_string(index=False)))
    
    # Get family mems in the order of the annotation vcf file.
    with open(c, mode="r") as fv:
        while True:
            line = fv.readline().strip("\n")
            if line.split("\t")[0] == "#CHROM":           
                family_mems = line.split("\t")[9:]
                break
        fv.close()
        
    proband = family_df.loc[family_df['Phenotype'].astype('float').astype('int') == 2, 'IndividualID'].tolist()[0]
    
    ID_mem_dict = {}
    for ind, mem in enumerate(family_mems):
        common_parents = [ x for x in family_df.loc[family_df['IndividualID'] == proband, 'PaternalID':'MaternalID'].squeeze(axis=0).tolist() if x in family_df.loc[family_df['IndividualID'] == mem, 'PaternalID':'MaternalID'].squeeze(axis=0).tolist() ]
        if mem == family_df.loc[family_df['IndividualID'] == proband, 'PaternalID'].tolist()[0]:
            ID_mem_dict[mem] = '_father_'
        elif mem == family_df.loc[family_df['IndividualID'] == proband, 'MaternalID'].tolist()[0]:
            ID_mem_dict[mem] = '_mother_'
        elif len(common_parents) > 0 and (common_parents != [0] and common_parents != ['0']):
            if mem == proband:
                ID_mem_dict[mem] = '_patient_child_'
            else:
                ID_mem_dict[mem] = '_sib_'
        else:
            if mem == proband:
                ID_mem_dict[mem] = '_patient_child_'
            else:
                ID_mem_dict[mem] = '_sib_'
        

    # Calculate the number of the samples
    global sample_number
    sample_number = len(family_df)
    logger.info("The family number is: " + str(sample_number))

    '''
    # Then we take a look at the family pedigree info to see whether there is only one father as the parent(normally it's mother)
    # Very important, dealing with situation that we only have one father as the parent. 1 stands for male, 2 stands for female, 0 stands for non-existent.
    if sample_number == 2:
        # now we only have one parent+proband
        os.chdir(d)
        with open(ped, "r") as ped:
            pedigree = ped.readlines()
            # cut off the headline and retain the records
            pedigree = pedigree[1:]
            for record in pedigree:
                # fam is the string variable we import from the system
                if record.split('\t')[0] == fam:
                    # If father and mother column is both 0 when sample_n value = 2, then the gender is the parent gender
                    if record.split('\t')[2] == record.split('\t')[3]:
                        parent_gender = int(record.split('\t')[4])
    else:
        parent_gender = 0
    '''

    # Now we start to make new headers based on sample numbers.
    # Define a preset family_member list
    '''
    family_member = ['_patient_child_', '_mother_', '_father_', '_sib1_', '_sib2_', '_sib3_', '_sib4_', '_sib5_',
                     '_sib6_', '_sib7_', '_sib8_', '_sib9_']
    # Define a preset family_member list for special situation that we have only father except the proband.
    family_member_op = ['_patient_child_', '_father_', '_sib1_', '_sib2_', '_sib3_', '_sib4_', '_sib5_', '_sib6_',
                        '_sib7_', '_sib8_', '_sib9_']
    '''
    for mem in family_mems:
        for field in format_fields_list:
            new_header = mem + ID_mem_dict[mem] + field
            if new_header not in headers:
                headers.append(new_header)
    '''
    for i in range(0, sample_number):
        # If we have both parents or no parent at all with the proband, then parent_gender = 0
        for field in format_fields_list:
            if parent_gender == 0 or parent_gender == 2:
                headers.append(vcf_heads[9 + i] + family_member[i] + field)
            elif parent_gender == 1:
                headers.append(vcf_heads[9 + i] + family_member_op[i] + field)
    '''
    
    logger.info("After reformatting, the headers now are: \n" + str(headers))
    # So now we have a new header line. We can start dealing with the contents.

    # Import tsv file into pandas dataframe. Set the column names as headers
    hit_column_length_limit = True
    min_itemsize = 3500
    while hit_column_length_limit:
        logger.info("Now entering the outer while loop, change the looping condition hit_column_length_limit to False.")
        hit_column_length_limit = False  
        os.chdir(a)
        hg19_txt = pd.read_csv(b, sep='\t', low_memory=False, names=headers, chunksize=50000, dtype=np.dtype("O"))
        hg19_txt_head = pd.read_csv(b, sep='\t', low_memory=False, names=headers, nrows=1)
        gc.collect()
        # Extract the index of column FORMAT
        global index_format
        index_format = hg19_txt_head.columns.get_loc('Otherinfo12')
        # # print('index_format =', index_format)
        # print("Reading main table costs: --- %s seconds ---" % (time.time() - start_time))
        first_chunk = next(hg19_txt).dropna(how='all').fillna(".")
        by_format = first_chunk.groupby('Otherinfo12', as_index=False, dropna=False)
        new_chunk = by_format.apply(find_index, fields=format_fields_list).astype(np.dtype("O")).dropna(how='all').fillna('.')
        # new_chunk = multi_process(func=find_index, data=first_chunk, num_process=8, fields=format_fields_list).fillna('.').drop_duplicates()
        buf = io.StringIO()
        new_chunk.info(verbose=True,buf=buf)
        logger.info("The first_chunk dtypes are: \n{}".format(buf.getvalue()))
        
        tmp_txt_path = tempfile.NamedTemporaryFile(delete=False).name
        if len(new_chunk) > 0: new_chunk.to_csv(tmp_txt_path, sep='\t', index=False)
        
        done_looping=False
        while not done_looping:
            try:
                chunk = next(hg19_txt).dropna(how='all').fillna(".")
                by_format = chunk.groupby('Otherinfo12', as_index=False, dropna=False)
                new_chunk = by_format.apply(find_index, fields=format_fields_list).astype(np.dtype("O")).dropna(how='all').fillna(".")
                buf = io.StringIO()
                new_chunk.info(verbose=True,buf=buf)
                # new_chunk = multi_process(func=find_index, data=chunk, num_process=8, fields=format_fields_list).fillna('.').drop_duplicates()
                # if len(new_chunk) > 0: new_chunk.to_hdf(tmp_hdf5_path, key='table', format='table', append=True)
                if len(new_chunk) > 0: new_chunk.to_csv(tmp_txt_path, sep='\t', index=False, header=False, mode='a')
            except StopIteration:
                done_looping=True
                # backup code: new_chunk = chunk.apply(find_index, axis=1, args=format_fields_list)
            except ValueError as ve:
                if "Trying to store a string with len" in repr(ve):
                    logger.exception("Some string column contain a string length over the size limit {} of the column, which is determined by first chunk.".format(min_itemsize))
                    hit_column_length_limit = True   # Make the big while loop running again.
                    min_itemsize = min_itemsize + 500    # Increase the size limit of column length
                    # cmd = "rm " + tmp_hdf5_path
                    # subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
                    logger.exception("Since the column string length limit hitting, check the min_itemsize {} now and the looping condition {} now.".format(min_itemsize, hit_column_length_limit))
                    break   # This will get out of the inner while loop, since the hit_column_length_limit = True now, the outer while loop will run again.
                elif "cannot match existing table structure" in repr(ve):
                    logger.exception("Some column's dtype does not fit with previous chunk df, check the dtypes of current chunk \n{}\ncurrent chunk shape is {}".format(buf.getvalue(), new_chunk.shape))
                    sys.exit(1)
            except KeyError:
                logger.exception("Why there is a key error? Check current chunk shape {} and dtypes {}. As well as current chunk column labels: \n{}".format(new_chunk.shape, buf.getvalue(), new_chunk.columns.tolist()))
                continue
                
            
    # start_time = time.time()
    # Export the dataframe to a tsv file. Use the argument index=False to ignore the extra index column.
    new_chunks = pd.read_csv(tmp_txt_path, sep='\t', low_memory=False).drop_duplicates()
    cmd = "rm " + tmp_txt_path
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
    
    logger.info("The new table looks like: \n" + new_chunks[:4].to_string(index=False))

    # Change the content of column other_info to the occurence number of genes in the variants list
    count_occurence = new_chunks.loc[:, 'Gene.refGene'].value_counts()
    count_occurence_df = count_occurence.to_frame().reset_index(level=0)  # Now original index column named as Gene.refGene while the later column named as count
    new_chunks = new_chunks.merge(count_occurence_df, on='Gene.refGene', how='left')

    index_otherinfo = new_chunks.columns.get_loc('Otherinfo1')
    # This command happens inplace
    new_chunks.drop(columns=['Otherinfo1'], inplace=True)
    co = new_chunks['count']
    # This insert command also happens inplace
    new_chunks.insert(index_otherinfo, 'Variant_number', co)
    new_chunks.drop(columns=['count'], inplace=True)

    # Output the variable to the new file
    new_headers = 'unknown1\tunknown2\tChr.1\tPOS\tID\tRef.1\tAlt.1\tQUAL.1\tFILTER.1\tINFO.1\tFORMAT'.split('\t')
    rename_dict = {}
    for i in range(2, 13): rename_dict["Otherinfo" + str(i)] = new_headers[i-2] 
    new_chunks.rename(columns=rename_dict, inplace=True)
    
    # Find out the label storing the VCF pos value.
    # pos_col_ind = [ i for i,v in enumerate(new_chunks.iloc[2, :].astype(str).tolist()) if re.search(r"^chr[0-9M]+", v) ][-1] + 1
    new_chunks.loc[:, "Start"] = new_chunks["POS"]
    
    logger.info("After renaming the Otherinfos columns, the table looks like:\n {}\n".format(new_chunks[:5].to_string(index=False)))

    # Drop redundant columns
    other_regex = re.compile("^Otherinfo1[0-9]$")
    tobe_dropped = [label for label in new_chunks.columns.tolist() if re.search(other_regex, label)]
    new_chunks.drop(columns=tobe_dropped, inplace=True)
    logger.info("After dropping the redundant columns, the table looks like:\n {}\n".format(new_chunks[:5].to_string(index=False)))

    new_chunks['Ref_new'], new_chunks['Alt_new'] = zip(*new_chunks.apply(change_ref_alt, axis=1))
    logger.info("After adding the two new ref alt columns, the table looks like:\n {}\n".format(new_chunks[:5].to_string(index=False)))
    
    for label in ['Ref_new', 'Alt_new']:
        loc = new_chunks.columns.get_loc(label.split('_')[0])
        tobe_mov = new_chunks.pop(label)
        new_chunks.drop(columns=[label.split('_')[0]], inplace=True)
        new_chunks.insert(loc=loc, column=label.split('_')[0], value=tobe_mov)

    new_chunks.drop(columns=["unknown1", "unknown2", "Chr.1", "POS", "ID", "Ref.1", "Alt.1", "QUAL.1", "FILTER.1", "INFO.1"], inplace=True, errors = "ignore")
    logger.info("The final table looks like:\n {}\n".format(new_chunks[:5].to_string(index=False)))
    
    # Output the final result
    new_chunks.to_csv(output, sep='\t', index=False, encoding='utf-8')




if __name__ == '__main__':
    parser = ap.ArgumentParser()
    parser.add_argument("-at", "--annotation_table", type=str, help="The path to the annotation table", required=True)
    parser.add_argument("-av", "--annotated_vcf", type=str, help="The path to the annotated vcf", required=False)
    parser.add_argument("-ped", "--ped_file", type=str, help="The path to the table of qname counts", required=True)
    parser.add_argument("-ot", "--output_table", type=str, help="Output annotation table", required=False)
    parser.add_argument("-fam", "--familyID", type=str, help="FamilyID of the processing anno table", required=False)
    
    args = parser.parse_args()
    
    anno_table = os.path.abspath(args.annotation_table)
    
    if args.annotated_vcf is None:
        annotated_vcf = anno_table[:-4] + ".vcf"
    else:
        annotated_vcf = args.annotated_vcf
        
    if args.output_table is None:
        output_table = anno_table[:-4] + ".reformat.txt"
    else:
        output_table = args.output_table
        
    if args.familyID is None:
        famID = os.path.basename(anno_table).split(".")[0]
    else:
        famID = args.familyID
        
    anno_dir = os.path.dirname(anno_table)
    ped_dir = os.path.dirname(os.path.abspath(args.ped_file))
    
    reformat(anno_dir, anno_table, annotated_vcf, ped_dir, os.path.basename(args.ped_file), output_table, famID)


