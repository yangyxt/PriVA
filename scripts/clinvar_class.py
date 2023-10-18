import pandas as pd
import numpy as np
import requests
import os
import logging
import sys
import re
import time
import subprocess
import io
import gzip
import time
import uuid
from lxml import etree as et
import tempfile
import argparse as ap
from python_utils import get_cds_len_from_ensembl_tranxid, \
                         deal_with_nullkey_group_return


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)



ncbi2ucsc_as_dict = {"GRCh38": "hg38", 
                    "hg38": "hg38",
                    "GRCh37": "hg19",
                    "hg19": "hg19"}


class UPDATED_CLINVAR(object):
    def __init__(self, \
                 clinvar_dir="/home/yangyxt/public_data/clinvar/{assembly}", \
                 path="Reformat_updated_Clinvar.tsv", \
                 vep="Reformat_updated_Clinvar.vepoutput", \
                 raw_tab="variant_summary.txt", \
                 vcf="updated_clinvar_simple.vcf", \
                 tab_url="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz", \
                 xml_url="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/clinvar_variation/weekly_release/ClinVarVariationRelease_00-latest_weekly.xml.gz", \
                 local_xml="ClinVarVariationRelease_00-latest_weekly.xml", \
                 funcsq_tab="allele_func_consq.map.csq.tsv", \
                 domain_cln_tab="Prot2hg.Reformat_updated_Clinvar.tsv",
                 domain_mh_tab="Prot2hg.Reformat_updated_Clinvar.mutation_hotspot.tsv", \
                 kde_table="updated_Clinvar_mutation_hotspot_density.tsv", \
                 kde_region_template="mutation_hotspot.kde.{chrom}.tsv", \
                 assembly="GRCh37", \
                 type_res = "tab", \
                 update = False, \
                 expiry = 180):
        clinvar_dir = clinvar_dir.format(assembly = ncbi2ucsc_as_dict[assembly])
        self.path = os.path.join(clinvar_dir, path)
        self.vep = os.path.join(clinvar_dir, vep)
        self.tab_url = tab_url
        self.xml_url = xml_url
        self.local_xml = os.path.join(clinvar_dir, local_xml)
        self.type_res = type_res
        self.raw_tab = os.path.join(clinvar_dir, raw_tab)
        self.vcf = os.path.join(clinvar_dir, vcf)
        self.update = update
        self.expiry = expiry
        self.assembly = assembly
        self.funcsq_tab = os.path.join(clinvar_dir, funcsq_tab)
        self.dr_cln_tab = os.path.join(clinvar_dir, domain_cln_tab)
        self.dr_mh_tab = os.path.join(clinvar_dir, domain_mh_tab)
        self.kde_table = os.path.join(clinvar_dir, kde_table)
        self.kde_region_temp = os.path.join(clinvar_dir, kde_region_template)
        # Here type_res can be either tab or xml or vep
    
    
    def value(self):
        if self.type_res == "xml":
            return self.read_clinvar_xml()
        elif self.type_res == "tab":
            return self.prepare_updated_reformat_clinvar()
        elif self.type_res == "raw":
            return self.prepare_updated_rawtab_clinvar(self.raw_tab)
        elif self.type_res == "vcf":
            return self.prepare_updated_vcf_clinvar(self.vcf)
        else:
            return self.prepare_updated_reformat_clinvar()
        

    def get_updated_vep_anno(self):
        return VEP_CLINVAR(path=self.prepare_updated_reformat_clinvar(), vep=self.vep, tab_url=self.tab_url).value()
    
    
    def get_updated_prot2hg(self):    
        return self.prepare_domain_clinvar_overlap()
    
    
    def get_update_mut_kde(self):
        return self.prepare_kde()


    def get_updated_tab_summary(self,  
                                return_whole=False):
        remote_url = self.tab_url
        local_path = self.raw_tab
        res = requests.get(remote_url)
        if return_whole:
            try:
                if re.search(r'\.gz$', remote_url):
                    with gzip.open(io.BytesIO(res.content), mode="rb") as gf:
                        return pd.read_csv(io.BytesIO(gf.read()), sep='\t', low_memory=False)
                else:
                    f = io.StringIO(res.content)
                    return pd.read_csv(f, sep='\t', low_memory=False)
            except Exception as e:
                logger.warning("Remote URL {} can't be accessed, run into error: {}. Try load local file {}".format(remote_url, e, local_path))
                return pd.read_csv(local_path, sep='\t', low_memory=False)
        else:
            try:
                if re.search(r'\.gz$', remote_url):
                    with gzip.open(io.BytesIO(res.content), mode="rb") as gf:
                        df_iter = pd.read_csv(io.BytesIO(gf.read()), sep='\t', low_memory=False, chunksize=100000)  # return bytes
                        return df_iter
                else:
                    f = io.StringIO(res.content)
                    df_iter = pd.read_csv(f, sep='\t', chunksize=100000)
                    return df_iter
            except Exception as e:
                logger.warning("Remote URL {} can't be accessed, run into error: {}. Try load local file {}".format(remote_url, e, local_path))
                df_iter = pd.read_csv(local_path, sep='\t', low_memory=False, chunksize=100000)
                return df_iter
            
            
    @staticmethod
    def iterate_xml(xmlfile, tag=None, logger=logger):
        '''
        Seems it can only be used to parse plain text XML file.
        '''
        if tag:
            doc = et.iterparse(xmlfile, events=('start', 'end'), tag=tag, encoding="utf-8")
        else:
            doc = et.iterparse(xmlfile, events=('start', 'end'), encoding="utf-8")
        _, root = next(doc)
        start_tag = None
        for event, element in doc:
            if event == 'start' and start_tag is None:
                start_tag = element.tag
            if event == 'end' and element.tag == start_tag:
                yield element
                start_tag = None
                root.clear()
                
                  
    def read_clinvar_xml(self):
        current = time.time()
        local_xml_path = self.local_xml
        updated_url = self.xml_url
        if (current - os.path.getmtime(local_xml_path))/3600 > (24*self.expiry):
            logger.warning("Local table file {} is at least 1 month old. Update it.".format(local_xml_path))
            try:
                import requests
                r = requests.get(updated_url)
                if re.search(r'\.gz$', updated_url):
                    with gzip.open(io.BytesIO(r.content), "rb") as gf:
                        f = io.BytesIO(gf.read())  # return bytes
                else:
                    f = io.BytesIO(r.content)
            except Exception as e:
                logger.warning("Remote URL {} can't be accessed, run into error: {}. Try load local file {}".format(updated_url, e, local_xml_path))
                if re.search(r'\.gz$', local_xml_path):
                    with gzip.open(local_xml_path, "rb") as lx:
                        f = io.BytesIO(lx.read())
                else:
                    with open(local_xml_path, 'rb') as lx:
                        f = io.BytesIO(lx.read())
        else:
            logger.warning("Local table file {} is at least 1 month old. Update it.".format(local_xml_path))
            if re.search(r'\.gz$', local_xml_path):
                with gzip.open(local_xml_path, "rb") as lx:
                    f = io.BytesIO(lx.read())
            else:
                with open(local_xml_path, 'rb') as lx:
                    f = io.BytesIO(lx.read())

        root_iter = self.iterate_xml(f)
        logger.info("Having a generator of 1st level child element under the root element {}.".format(root_iter))
        return root_iter

    
    def get_highest_patho_afs(self):
        reformat_clinvar_vcf = self.vcf
        annovar_tsv = reformat_clinvar_vcf.replace(".vcf", ".pathogenic.hg19_multianno.txt")
        patho_clinvar_vcf = reformat_clinvar_vcf.replace(".vcf", ".pathogenic.vcf")
        reformat_clinvar_df = pd.read_table(reformat_clinvar_vcf, low_memory=False)
        pathogenic_df = reformat_clinvar_df.loc[reformat_clinvar_df["INFO"].astype(str).str.contains(r"[Pp]athogenic;", regex=True), :]
        pathogenic_df.to_csv(patho_clinvar_vcf, sep="\t", index=False)
        cmd = "if [[ {annovar_tab} -ot {raw_vcf} ]]; then bash /paedyl01/disk1/yangyxt/ngs_scripts/common_bash_utils.sh basic_annovar {patho_vcf}; fi".format(annovar_tab = annovar_tsv,
                                                                                                                                                            raw_vcf = reformat_clinvar_vcf,
                                                                                                                                                            patho_vcf = patho_clinvar_vcf)
        res = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8", check=True)
        print(res.stdout, file=sys.stderr)
        annovar_df = pd.read_table(annovar_tsv, low_memory=False, usecols=["Chr","Start","End","Ref","Alt","Func.refGene","Gene.refGene","ExonicFunc.refGene","AAChange.refGene","gnomAD_exome_ALL","gnomAD_genome_ALL"]).drop_duplicates()
        annovar_df = annovar_df.loc[annovar_df["Gene.refGene"] != ".", :].replace(".", np.nan)
        def get_max_af_pergene(groupdf):
            groupdf["max_AF"] = groupdf[["gnomAD_exome_ALL", "gnomAD_genome_ALL"]].max(axis=1)
            return groupdf
        annovar_df = annovar_df.groupby("Gene.refGene", as_index=False).apply(get_max_af_pergene)
        annovar_df.to_csv(annovar_tsv, sep="\t", index=False)
        return annovar_tsv
    
    
    def get_alleleID_funcconsq_mapping(self, output, threads=10, chunksize=1000):
        if os.path.exists(output):
            if os.path.getmtime(self.local_xml) > os.path.getmtime(output):
                logger.warning("Local table file {} is older than the XML file {}. Update it.".format(output, self.local_xml))
                update = True
            else:
                update = False
        else:
            update = True
            
        if update:
            chunk_iter = self.read_clinvar_xml()
            headers = ["#AlleleID", "#VariantID", "Name", "Type", "Func_Consq"] # The header is corresponding with the func retrieve_alleleID_and_funcconsq
            # If we set chunksize for func read_clinvar_xml, the yielding part is a big chunk of original files
            tmp_suffix = str(uuid.uuid4())
            tmp_output = output + "." + tmp_suffix
            #pool = mp.Pool(threads)
            n = 0
            while True:
                try:
                    chunk = next(chunk_iter) # The chunk is a list of root_children
                except StopIteration:
                    break
                else:
                    n += 1
                    results = list(map(self.retrieve_alleleID_and_funcconsq, chunk)) # results are list of list_of_tuple
                    rows = [t for l in results for t in l]
                    op_chunk_df = pd.DataFrame(rows, columns=headers)
                    if n <= 1:
                        logger.debug("Output the first chunk df, containing {} rows and {} columns. It looks like: \n{}\n".format(op_chunk_df.shape[0], op_chunk_df.shape[1], op_chunk_df[:4].to_string(index_names=False)))
                        op_chunk_df.drop_duplicates().to_csv(tmp_output, sep='\t', index=False)
                    else:
                        logger.debug("Output this chunk df containing {} rows and {} columns.".format(op_chunk_df.shape[0], op_chunk_df.shape[1]))
                        op_chunk_df.drop_duplicates().to_csv(tmp_output, sep='\t', index=False, header=False, mode="a")
                        
            result = subprocess.run("mv {} {} && ls -lh {}".format(tmp_output, output, output), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
            logger.info("\n{}".format(result.stdout))
            return output
        else:
            return output


    @staticmethod
    def retrieve_alleleID_and_funcconsq(root_child, logger=logger):
        if type(root_child) == str:
            root_child = et.fromstring(root_child)
        row_set = []
        allele_lst = root_child.findall(path=".//{}".format("SimpleAllele"))
        genotype_lst = root_child.findall(path=".//{}".format("Genotype"))
        haplotype_lst = root_child.findall(path=".//{}".format("Haplotype"))
        allele_lst = allele_lst + genotype_lst + haplotype_lst
        for allele in allele_lst:
            attrib_dict = dict(allele.attrib)
            alleleID = attrib_dict.get("AlleleID", "NA")
            varID = attrib_dict.get("VariationID", "NA")
            func_consq = None
            for se in allele:
                if se.tag == "FunctionalConsequence":
                    if len([ s for s in se if s.tag == "XRef"]) > 0:
                        xrefs = " [" + ", ".join([ str(s.attrib["DB"]) + " " + str(s.attrib["ID"]) for s in se if s.tag == "XRef" ]) + "]"
                    else:
                        xrefs = ""
                    str_value = str(se.attrib["Value"]) + xrefs
                    if func_consq:
                        func_consq = func_consq + "; " + str_value
                    else:
                        func_consq = str_value
                elif se.tag == "Name":
                    allele_name = se.text.strip("\n ") if type(se.text) == str else "NA"
                elif se.tag == "VariantType":
                    allele_vt = se.text.strip("\n ") if type(se.text) == str else "NA"
            row_set.append(tuple([alleleID, varID, allele_name, allele_vt, func_consq])) #Return a list of tups
        return row_set
    

    @staticmethod
    @deal_with_nullkey_group_return
    def summary_clinvar_per_gene(groupdf, 
                                 two_star_label="same_gene_CLNSIG_C2",
                                 revstat_col="ReviewStatus",
                                 clnsig_col="ClinicalSignificance",
                                 logger = logger):
        '''
        This function is used to inspect whether a gene contain a high confidence LoF Pathogenic variant
        '''
        
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
    
        rs_bools_two_star = groupdf[revstat_col].map(cln_rs_qdict) >= 2

        # First summary on clinvar sig
        sig_dict = groupdf.loc[(rs_bools_two_star).to_numpy(), clnsig_col].value_counts().to_dict()
        sum_value = ", ".join([ str(k) + ":" + str(v) for k, v in sig_dict.items() ])
        sum_value = sum_value if len(sum_value) > 0 else np.nan
        groupdf[two_star_label] = sum_value
        valid_funcsq = groupdf.loc[groupdf["Func_Consq"].str.len() > 3, "Func_Consq"].drop_duplicates().tolist()
        groupdf["Patho_funcsq_per_gene"] = " | ".join(valid_funcsq)

        return groupdf.drop_duplicates()
    
    
    @staticmethod
    def convert_to_vcf_perchunk(df, assembly="GRCh37", logger=logger):
        # First filter out the variants recorded based on GRCh38
        df = df.loc[df['Assembly'] == assembly, :]  # remove out complex variants
        # Filter on chromosomes
        df = df.loc[df['Chromosome'].str.contains(r'^[0-9XYMT]+$', regex=True), :]
        old_columns = df.columns.tolist()
        # Then create vcf columns
        df["#CHROM"] = "chr" + df["Chromosome"].astype(str).str.replace("MT", "M")
        df["POS"] = df["PositionVCF"]
        df["ID"] = df["#AlleleID"]
        df["REF"] = df["ReferenceAlleleVCF"]
        df["ALT"] = df["AlternateAlleleVCF"]
        df["QUAL"] = 99
        df["FILTER"] = "PASS"
        df["INFO"] = "END=" + df["Stop"].astype(str) + ";" + \
                     "ClinicalSignificance=" + df["ClinicalSignificance"] + ";" + \
                     "ReviewStatus" + df["ReviewStatus"] + ";" + \
                     "PhenotypeList" + df["PhenotypeList"] + ";" + \
                     "PhenotypeIDS" + df["PhenotypeIDS"]
        return df.drop(columns=old_columns).drop_duplicates()
        
        
        
    @staticmethod
    def adjust_format_per_chunk(df, assembly="GRCh37", logger=logger):
        # First filter out the variants recorded based on GRCh38
        df = df.loc[df['Assembly'] == assembly, :]  # remove out complex variants
        # Then filter on the Chromosome
        df = df.loc[df['Chromosome'].str.contains(r'^[0-9XYMT]+$', regex=True), :]
        # Filter out complex
        df = df.loc[df['Type'] != "Complex", :]
        # Filter out translocation
        df = df.loc[df['Type'] != "Translocation", :]
        # Filter out variation
        df = df.loc[df['Type'] != "Variation", :]
        # Categorize on the variant type
        # First pick out variantion type
        variations = df.loc[df['Type'] == "Variation", :]
        variations = variations.loc[variations['Name'].str.contains(r'c\..+>[ATCG]'), :]  # Pick out cds change
        variations.loc[:, "ReferenceAlleleVCF"] = variations['Name'].str.extract(r'c\.-*[0-9]+([ATCG])>[ATCG]', expand=False)
        variations.loc[:, "AlternateAlleleVCF"] = variations['Name'].str.extract(r'c\.-*[0-9]+[ATCG]>([ATCG])', expand=False)
        variations = variations.loc[variations['ReferenceAlleleVCF'] != variations['AlternateAlleleVCF'], :]
        df = df.loc[df['Type'] != "Variation", :]
        df = pd.concat([df,variations], ignore_index=True)
        # Last pick out 
        # Then filter on the pos and ref/alt
        df = df.loc[df['Start'].astype(str) != "-1", :]
        # Change column order
        chr_col = df.pop('Chromosome').str.replace("MT", "M")
        df.insert(0, "Chr", "chr" + chr_col)
        start_col = df.pop("Start")
        start_col = df.pop("PositionVCF")
        df.insert(1, "Start", start_col)
        end_col = df.pop("Stop")
        end_col = start_col.astype(int) + df["ReferenceAlleleVCF"].astype(str).str.len() - 1
        df.insert(2, "End", end_col)
        ref_col = df.pop("ReferenceAlleleVCF")
        df.insert(3, "Ref", ref_col)
        alt_col = df.pop("AlternateAlleleVCF")
        df.insert(4, "Alt", alt_col)
        # Change REF/ALT notation, first prepare notation based on the Type
        replace_dict = {'Complex': '<COMP>', \
                        "Indel":"<SUB>", \
                        "Insertion": "<INS>", \
                        "Deletion":"<DEL>", \
                        "Duplication":"<DUP>", \
                        "Microsatelite": "<DUP:TANDEM>", \
                        "copy number gain":"<DUP>", \
                        "copy number loss":"<DEL>", \
                        "fusion":"<FUS>", \
                        "Invertion":"<INV>", \
                        "single nucleotide variant":"<SNV>", \
                        "Tandem duplication":"<DUP:TANDEM>", \
                        "Translocation":"<TRA>"}
        df['SVtype'] = df['Type'].replace(replace_dict)
        df['Ref'] = np.where(ref_col == "na", "N", ref_col)
        df['Alt'] = np.where(alt_col == "na", df['SVtype'], alt_col)
        df.drop(columns = ['SVtype'], inplace=True)
        return df.drop_duplicates()


    @staticmethod
    def anno_funcsq_to_clinvar(ref_clinv_df, funcsq_df, logger=logger):
        ref_clinv_df = ref_clinv_df.drop(columns=["Func_Consq", "Patho_funcsq_per_gene"], errors="ignore").drop_duplicates()
        funcsq_df.loc[:, "Func_Consq"] = funcsq_df["Func_Consq"].fillna("")
        ref_clinv_df = ref_clinv_df.merge(funcsq_df[["#AlleleID", "Func_Consq"]].drop_duplicates().dropna(subset=["#AlleleID"]), on="#AlleleID", how="left")
        return ref_clinv_df
    
    
    def prepare_kde(self):    
        # Only select short variants with high deleterious effect. 
        # Do not limit on the ClinVar Reviewstatus, only limit on the clinvar significance
        import multiprocessing as mp
        latest_clnvar = self.path
        output_table = self.kde_table
        path_template = self.kde_region_temp
        hg19_bed = path_template.format(chrom = "hg19").replace(".tsv", ".hotspot.bed")
        
        if self.update == False:
            return hg19_bed
        elif os.path.getmtime(hg19_bed) < os.path.getmtime(self.path):
            logger.warning("The mutation hotspot region bed file {} is older than the updated clinvar record table {}".format(self.hg19_bed, self.path))
            update = True
        else:
            update = False
            
        if update == False:
            return hg19_bed
        
        clinvar_df = pd.read_table(latest_clnvar, low_memory=False)
        var_type_bools = (clinvar_df["End"] - clinvar_df["Start"]) <= 10 
        patho_bools = clinvar_df["ClinSigSimple"].astype(str) == "1"
        missense_bools = clinvar_df["Name"].str.contains(r"\(p\.[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}\)$", regex=True)
        stop_gain = clinvar_df["Name"].str.contains(r"\(p\.[A-Z][a-z]{2}[0-9]+Ter\)$", regex=True)
        frameshift = clinvar_df["Name"].str.contains(r"\(p\.[A-Z][a-z]{2}[0-9]+fs\)$", regex=True)

        patho_mis_clinvar_df = clinvar_df.loc[ var_type_bools & patho_bools & missense_bools & ~stop_gain & ~frameshift, :]
        
        # First group the dataframe by chromosome
        by_chr = patho_mis_clinvar_df.groupby("Chr", as_index=False)
        all_chrs = patho_mis_clinvar_df.loc[:, "Chr"].drop_duplicates().to_list()

        # Prepare the chromosome length
        ref_genome="/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta"

        from Bio import SeqIO
        input_file = open(ref_genome)
        chr_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))

        pool = mp.Pool(6)
        # chr_new_df = by_chr.apply(perform_kde_per_chr, chr_dict=chr_dict)
        results = pool.starmap(self.perform_kde_per_chr, [(df, chr_dict, path_template) for name, df in by_chr])
        
        pd.concat(results).reset_index(drop=True).drop_duplicates().reset_index(drop=True).to_csv(output_table, index=False, sep="\t")
        
        # Now we would like to concatenate all the bedfiles into one
        hg19_bed = path_template.format(chrom = "hg19").replace(".tsv", ".hotspot.bed")
        cmd1 = "cat " + " ".join(path_template.format(chrom = chromosome).replace(".tsv", ".hotspot.bed") for chromosome in all_chrs)
        cmd2 = "sort -k 1,1 -k2,2n > " + hg19_bed + " && ls -lh " + hg19_bed
        p1 = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr = subprocess.STDOUT)
        p2 = subprocess.run(cmd2, stdin=p1.stdout, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        logger.info("The final bed file is {}".format(p2.stdout))
        return hg19_bed
    
    
    
    @staticmethod
    def post_process_mut_hotspot_tab(path_template="/paedyl01/disk1/yangyxt/public_data/clinvar/tab_delimited/mutation_hotspot.kde.{chrom}.hotspot.tsv",
                                     chromosome="",
                                     logger=logger):
        import pybedtools as pb
        assert len(chromosome) > 0
        path = path_template.format(chrom=chromosome)
        df = pd.read_table(path, low_memory=False)
        columns = df.columns.tolist()
        # Now we import pybedtools to generate a bed file
        bed = pb.BedTool.from_dataframe(df)
        merged = bed.merge(c=4, o="mean", d=9)
        # Convert the merged bed file to dataframe
        merged_df = merged.to_dataframe(disable_auto_names=True,
                                        names = columns)
        # Output the merged dataframe
        logger.info("Before merging, the table has {} records and covers {} bp. After merging, the table has {} records and covers {} bp.".format(df.shape[0],
                                                                                                                                                (df.iloc[:, 2] - df.iloc[:, 1]).sum() + df.shape[0],
                                                                                                                                                merged_df.shape[0],
                                                                                                                                                (merged_df.iloc[:, 2] - merged_df.iloc[:, 1]).sum() + merged_df.shape[0]))
        merged_df.to_csv(path.replace(".tsv", ".bed"), index=False, header=False, sep="\t")
        result = subprocess.run("ls -lh {}".format(path.replace(".tsv", ".bed")), check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
        logger.info("The subprocess run result is: {}".format(result.stdout))
        
    
    @staticmethod
    def perform_kde_per_chr(groupdf, 
                            chr_dict=None, 
                            chr_bed_temp="/paedyl01/disk1/yangyxt/public_data/clinvar/tab_delimited/mutation_hotspot.kde.{chrom}.tsv",
                            ref_genome="/paedyl01/disk1/yangyxt/indexed_genome/ucsc.hg19.fasta"):
    
        from numpy import array, linspace
        from sklearn.neighbors import KernelDensity
        # First determine the total length of the chromosome.
        chromosome = groupdf.loc[:, "Chr"].to_list()[0]
        
        if chr_dict:
            chr_size = len(chr_dict[chromosome].seq)
        else:
            from Bio import SeqIO
            input_file = open(ref_genome)
            chr_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
            chr_size = len(chr_dict[chromosome].seq)
            
        # Then determine the coordinates of the pathogenic short variants
        mut_positions = groupdf.loc[:, ["Start", "End"]].mean(axis=1).round().astype(int)
        groupdf = groupdf.assign(middle_pos=mut_positions)

        # Prepare the input array for KDE to fit
        arr = mut_positions.sort_values().to_numpy().reshape(-1, 1)
        logger.info("The chromosome {} size is: {}, input pathogenic mutation positions are :\n{}".format(chromosome, chr_size, arr))
        
        # perform KDE on arr, first have a fitting kde model. 
        kde = KernelDensity(kernel="epanechnikov", bandwidth=100).fit(arr)
        
        # Then use the model to calculate the mutation density for every base in the chromosome
        space = linspace(1, chr_size, num=chr_size).astype(np.longdouble)
        estimate = kde.score_samples(space.reshape(-1,1)).astype(np.longdouble)
        
        # Then Find out the top 10 percentile of the mutation hotspot in the 
        chr_den_bed = pd.DataFrame({"Chr": np.array([chromosome for i in range(0, chr_size)]), 
                                    "Start": space.astype(np.int64),
                                    "End": space.astype(np.int64),
                                    "mutation_density": estimate}, 
                                columns=["Chr",
                                            "Start",
                                            "End",
                                            "mutation_density"],
                                dtype="object")
        
        logger.info("The density recording dataframe for chromosome {} has a dtypes of: \n{}\n. And it looks like: \n{}\n".format(chromosome,
                                                                                                                                chr_den_bed.dtypes, 
                                                                                                                                chr_den_bed.iloc[:5, :].to_string(index=False)))
        
        # Now need to conver the column storing mutation density to long double in case of storing negative infinite value
        chr_den_bed.loc[:, "mutation_density"] = chr_den_bed.loc[:, "mutation_density"].astype(np.longdouble)
        chr_den_bed.loc[:, "Start"] = chr_den_bed.loc[:, "Start"].astype(np.int64)
        chr_den_bed.loc[:, "End"] = chr_den_bed.loc[:, "End"].astype(np.int64)
        
        logger.info("After datatype conversion for chromosome Start, End, and mutation density, the table now has dtypes like {}. And it looks like :\n{}\n".format(chr_den_bed.dtypes, 
                                                                                                                                                                    chr_den_bed.iloc[:5, :].to_string(index=False)))
        chr_den_bed.to_csv(chr_bed_temp.format(chrom=chromosome), sep="\t", index=False)

        # Then identify the top 10 percentile
        cutoff = chr_den_bed.mutation_density.quantile(0.95, interpolation="nearest")
        chr_hotspot = chr_den_bed.loc[chr_den_bed.loc[:, "mutation_density"] > cutoff, :]
        logger.info("After selecting the top 10 percentile bases with highest mutation densities, the remaining table for chromosome {} has {}bp covered.".format(chromosome, chr_hotspot.shape[0]))
        
        # Output the top 10 percentile bases
        chr_hotspot.to_csv(chr_bed_temp.format(chrom=chromosome + ".hotspot"), sep="\t", index=False)
        
        # Convert the dataframe to bed file and save them again    
        UPDATED_CLINVAR.post_process_mut_hotspot_tab(path_template=chr_bed_temp.replace(".tsv", ".hotspot.tsv"), chromosome=chromosome)
        
        groupdf = groupdf.assign(mis_patho_density=kde.score_samples(groupdf["middle_pos"].to_numpy().reshape(-1,1))).drop_duplicates()
        return groupdf

        
    def prepare_updated_reformat_clinvar(self):
        output = self.path
        update = self.update
        if os.path.exists(output):
            if update:
                logger.warning("User specify to force update the Clinvar records.")
            elif os.path.getmtime(self.raw_tab) > os.path.getmtime(output):
                logger.warning("Local table file {} is older than the input file {}. Update it.".format(output, self.raw_tab))
                update = True
            else:
                update = False
        else:
            update = True
            
        if update:
            tmp_suffix = str(uuid.uuid4())
            tmp_output = output + "." + tmp_suffix
            df_iter = self.get_updated_tab_summary()
            funcsq_df = pd.read_table(self.get_alleleID_funcconsq_mapping(output=self.funcsq_tab), low_memory=False)
            first_chunk = next(df_iter)
            first_chunk = self.anno_funcsq_to_clinvar(self.adjust_format_per_chunk(first_chunk, assembly=self.assembly), funcsq_df)
            first_chunk.to_csv(tmp_output, sep="\t", index=False, encoding="utf-8")
            while True:
                try:
                    chunk = next(df_iter)
                except StopIteration:
                    break
                else:
                    chunk = self.anno_funcsq_to_clinvar(self.adjust_format_per_chunk(chunk, assembly = self.assembly), funcsq_df)
                    chunk.to_csv(tmp_output, sep="\t", index=False, header=False, encoding="utf-8", mode="a")
            # summary the pathogenic variants per gene
            df = pd.read_table(tmp_output, low_memory=False)
            logger.info("The updated clinvar reformat table at tmp path {} looks like:\n{}".format(tmp_output, df[:5].to_string(index=False)))
            df.groupby("GeneSymbol", as_index=False, dropna=False).apply(self.summary_clinvar_per_gene).drop_duplicates().to_csv(tmp_output, sep="\t", index=False)
            result = subprocess.run("mv {} {} && ls -lh {}".format(tmp_output, output, output), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
            logger.info("\n{}".format(result.stdout))
            return output
        else:
            return output
        
    
    def prepare_updated_vcf_clinvar(self, output):
        assembly = self.assembly
        if os.path.exists(output):
            if os.path.getmtime(self.path) > os.path.getmtime(output):
                logger.warning("Local table file {} is older than the input tab file {}. Update it.".format(output, self.path))
                update = True
            else:
                update = False
        else:
            update = True
            
        if update:
            tmp_suffix = str(uuid.uuid4())
            tmp_output = output + "." + tmp_suffix
            df_iter = self.get_updated_tab_summary()
            first_chunk = next(df_iter)
            self.convert_to_vcf_perchunk(first_chunk, assembly).to_csv(tmp_output, sep="\t", index=False, encoding="utf-8")
            while True:
                try:
                    chunk = next(df_iter)
                except StopIteration:
                    break
                else:
                    self.convert_to_vcf_perchunk(chunk, assembly).to_csv(tmp_output, sep="\t", index=False, header=False, encoding="utf-8", mode="a")
            # summary the pathogenic variants per gene
            result = subprocess.run("mv {} {} && ls -lh {}".format(tmp_output, output, output), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
            logger.info("\n{}".format(result.stdout))
            return output
        else:
            return output
        
        
    def prepare_updated_rawtab_clinvar(self, output):
        assembly = self.assembly
        if os.path.exists(output):
            if os.path.getmtime(self.path) > os.path.getmtime(output):
                logger.warning("Local table file {} is older than the input file {}. Update it.".format(output, self.path))
                update = True
            else:
                update = False
        else:
            update = True
            
        if update:
            tmp_suffix = str(uuid.uuid4())
            tmp_output = output + "." + tmp_suffix
            df_iter = self.get_updated_tab_summary()
            first_chunk = next(df_iter)
            first_chunk.loc[first_chunk["Assembly"] == assembly, :].to_csv(tmp_output, sep="\t", index=False, encoding="utf-8")
            while True:
                try:
                    chunk = next(df_iter)
                except StopIteration:
                    break
                else:
                    chunk.loc[chunk["Assembly"] == assembly, :].to_csv(tmp_output, sep="\t", index=False, header=False, encoding="utf-8", mode="a")
            # summary the pathogenic variants per gene
            result = subprocess.run("mv {} {} && ls -lh {}".format(tmp_output, output, output), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
            logger.info("\n{}".format(result.stdout))
            return output
        else:
            return output
        

    
    def prepare_domain_clinvar_overlap(self, 
                                       domain_region="/home/yangyxt/public_data/Prot2HG/hg19_prot2hg_summary.gnomAD_v2_exome_overlap.count_obs.constraint_metrics.tsv"):
        assembly = self.assembly
        if os.path.exists(self.dr_cln_tab):
            if os.path.getmtime(self.path) > os.path.getmtime(self.dr_cln_tab):
                logger.warning("Local table file {} is older than the input file {}. Update it.".format(self.dr_cln_tab, self.path))
                update = True
            else:
                update = False
        else:
            logger.warning("Local table file {} not existed. Create it.".format(self.dr_cln_tab))
            update = True
            
        if update:
            import pybedtools as pb
            
            clinvar_table = self.prepare_updated_reformat_clinvar(self.path, update=update)
            domain_region_remained_cols = ["feature_name", "note", "mis_obs_num", "lof_obs_num", "mis_exp_num", "lof_exp_num", "mis_oe_pvalue", "lof_oe_pvalue"]
            # domain_region_dropped_cols = ['chrom', 'chr_start', 'chr_end', 'strand', 'gene', 'protein_ID', 'gene_ID', 'type', 'chromosome', 'start_position', 'end_position', 'cds_length', 'num_coding_exons', 'canonical', 'exp_mis', 'mu_mis', 'possible_mis', 'exp_syn', 'mu_syn']
            dr_df = pd.read_csv(domain_region, sep='\t', low_memory=False, na_values=[".", "-"]).dropna(how="all")
            domain_region_dropped_1 = [ c for c in dr_df.columns.tolist() if c not in domain_region_remained_cols and c not in ["chrom", "chr_start", "chr_end", "strand", "gene"] ]
            # domain_region_dropped_2 = ["chrom", "chr_start", "chr_end", "strand", "gene"]
            dr_df.drop(columns=domain_region_dropped_1, inplace=True, errors="ignore")
            
            clinvar_table_dropped_1 = ['GeneID', 'GeneSymbol', 'HGNC_ID', 'ClinSigSimple', 'LastEvaluated', 'PhenotypeIDS', 'Origin', 'OriginSimple', 'Assembly', 'ChromosomeAccession', 'ReferenceAllele', 'AlternateAllele', 'Cytogenetic', 'NumberSubmitters', 'Guidelines', 'TestedInGTR', 'OtherIDs', 'SubmitterCategories', 'VariationID', 'PositionVCF']
            # clinvar_table_dropped_2 = ['Chr', 'Start', 'End', 'Ref', 'Alt']
            ct_df = pd.read_csv(clinvar_table, sep='\t', low_memory=False, na_values=[".", "-"]).dropna(how="all").drop_duplicates()
            ct_df = ct_df.loc[ct_df["Assembly"].str.contains(r"{}".format(assembly)), :]
            ct_df = ct_df.drop(columns=clinvar_table_dropped_1, errors="ignore").drop_duplicates()
            ct_df = ct_df.loc[ct_df["Type"] != "copy number gain", :]
            
            dr_df = dr_df.loc[dr_df.iloc[:, 2].astype(str).str.contains(r"^[0-9]+$") & dr_df.iloc[:, 1].astype(str).str.contains(r"^[0-9]+$"), :].drop_duplicates()
            ct_df = ct_df.loc[ct_df.iloc[:, 2].astype(str).str.contains(r"^[0-9]+$") & ct_df.iloc[:, 1].astype(str).str.contains(r"^[0-9]+$"), :].drop_duplicates()
            
            logger.info("Dataframe {} has {} rows and {} columns, looks like this:\n{}\n".format("domain_table", dr_df.shape[0], dr_df.shape[1], dr_df[:10].to_string()))
            logger.info("Dataframe {} has {} rows and {} columns, looks like this:\n{}\n".format("Clinvar_table", ct_df.shape[0], ct_df.shape[1], ct_df[:10].to_string()))
            
            clinvar_bed = pb.BedTool.from_dataframe(ct_df)
            clinvar_headers = ct_df.columns.tolist()

            domain_bed = pb.BedTool.from_dataframe(dr_df)
            domain_headers = dr_df.columns.tolist()
            
            domain_clinvar_bed = domain_bed.intersect(clinvar_bed, wo=True, f=0.8, F=0.8, e=True)
            domain_clinvar_df = domain_clinvar_bed.to_dataframe(disable_auto_names=True, names=domain_headers + clinvar_headers + ["overlap_len"], dtype=str).drop(columns="overlap_len").drop_duplicates().replace({"nan":np.nan})
            
            tmp_out = self.dr_cln_tab + "." + str(uuid.uuid4())
            domain_clinvar_df.to_csv(tmp_out, sep="\t", index=False)
            result = subprocess.run("mv {} {} && ls -lh {}".format(tmp_out, self.dr_cln_tab, self.dr_cln_tab), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
            logger.info("\n{}".format(result.stdout))
            domain_clinvar_df = domain_clinvar_df.loc[(domain_clinvar_df["End"].astype(np.longdouble) - domain_clinvar_df["Start"].astype(np.longdouble)) <= 100000, :]
            domain_clinvar_df.to_csv(tmp_out, sep="\t", index=False)
            result = subprocess.run("mv {} {} && ls -lh {}".format(tmp_out, self.dr_cln_tab.replace(".tsv", ".shortv.tsv"), self.dr_cln_tab.replace(".tsv", ".shortv.tsv")), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
            logger.info("\n{}".format(result.stdout))
            
        return self.dr_cln_tab
        
        


class VEP_CLINVAR(UPDATED_CLINVAR):
    def __init__(self, 
                 type_res="tab",
                 update = False,
                 expiry = 180,
                 threads = 4,
                 clinvar_dir = "/home/yangyxt/public_data/clinvar",
                 path="/Reformat_updated_Clinvar.tsv",
                 vep="/paedyl01/disk1/yangyxt/public_data/clinvar/tab_delimited/Reformat_updated_Clinvar.vepoutput",
                 trunc_sum="/paedyl01/disk1/yangyxt/public_data/clinvar/tab_delimited/Reformat_updated_Clinvar.vepoutput.truncsum.tsv",
                 tab_url="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"):
        super(VEP_CLINVAR, self).__init__(path=path, vep=vep, tab_url=tab_url, type_res=type_res)
        self.path = path
        self.vep = vep
        self.tab_url = tab_url
        self.update = update
        self.expiry = expiry
        self.threads = threads
        self.trunc_sum = trunc_sum
        logger.info("The updated clinvar path is {}, and the updated vep anno result is {}".format(self.path, self.vep))
    
    def value(self):
        return self.main_convert(input_tab = self.path, output_vep = self.vep)
    
    def updated_vcf(self):
        return self.convertab2vcf(self.path)
        
    @staticmethod
    def convert2vepvcf_per_row(row, all_labels=""):
        row["#CHROM"] = row["Chr"]
        row["POS"] = row["Start"]
        row["ID"] = row["#AlleleID"]
        row["REF"] = row["Ref"]
        row["ALT"] = row["Alt"]
        # Deal with other complex variants
        if type(row["Alt"]) != str:
            row["ALT"] = np.nan
        elif re.search(r"^[ATCG]+$",row["Alt"]):
            pass
        elif "<" in row["Alt"]:
            pass
        elif re.search(r"^[BD-FH-SU-Za-z]$", row["Alt"]):
            row["ALT"] = np.nan
        else:
            row["ALT"] = "<" + row["ALT"] + ">"
        row["QUAL"] = "99"
        row["FILTER"] = "PASS"
        row.replace({"<DUP:TANDEM>":"<TDUP>"}, inplace=True)
        if type(row["ALT"]) == str:
            if "<" in row["ALT"]:
                row["INFO"] = "SVTYPE=" + row["ALT"].strip("<>") + ";" + "END=" + str(row["End"])
            else:
                row["INFO"] = "END=" + str(row["End"])
        else:
            row["INFO"] = "END=" + str(row["End"])
        logger.debug("\n{}".format(row))
        added_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        row.drop(labels = [ l for l in all_labels.split(",") if l not in added_cols ], inplace=True, errors="ignore")
        if row["ALT"] == "<DUP:TANDEM>" or row["ALT"] == "<TDUP>": logger.info("\n{}".format(row.tolist()))
        if row["REF"] == row["ALT"]:
            row["ALT"] = np.nan
        return row
    
    
    @staticmethod    
    def merge_with_clinvar(vep_result="/paedyl01/disk1/yangyxt/public_data/clinvar/tab_delimited/Reformat_updated_Clinvar.vepoutput",
                           clinvar_tab="/paedyl01/disk1/yangyxt/public_data/clinvar/tab_delimited/Reformat_updated_Clinvar.tsv",
                           output_tab="/paedyl01/disk1/yangyxt/public_data/clinvar/tab_delimited/Reformat_updated_Clinvar.vep.tsv"):
        vep_df = pd.read_table(vep_result, low_memory=False)
        clinvar_df = pd.read_table(clinvar_tab, low_memory=False)
        merged = clinvar_df.merge(vep_df.rename(columns={"#Uploaded_variation": "#AlleleID", "Feature": "TranscriptID"}).drop(columns=["Location", "Allele", "Feature_type"]).dropna(subset=["#AlleleID"]).drop_duplicates(), how="left", on="#AlleleID")
        logger.info("This is the results of annotating all clinvar records with VEP:\n{}".format(merged[:5].to_string(index=False)))
        tmp_suffix = "." + str(uuid.uuid4())
        tmp_output = output_tab[:-4] + tmp_suffix + ".tsv"
        merged.to_csv(tmp_output, sep="\t", index=False)
        result = subprocess.run("mv {} {} && ls -lh {}".format(tmp_output, output_tab, output_tab), \
                                shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
        logger.info("\n{}".format(result.stdout))
        return output_tab
    
    
    @staticmethod
    @deal_with_nullkey_group_return
    def sum_trunc_var_per_tranx(groupdf, csq_col="Consequence"):
        trunc_bools = groupdf.loc[:,csq_col].astype(str).str.contains(r"stop_gained") | \
                    groupdf.loc[:,csq_col].astype(str).str.contains(r"stop_lost") | \
                    groupdf.loc[:,csq_col].astype(str).str.contains(r"start_lost") | \
                    groupdf.loc[:,csq_col].astype(str).str.contains(r"frameshift_variant") | \
                    groupdf.loc[:,csq_col].astype(str).str.contains(r"protein_altering_variant") | \
                    groupdf.loc[:,csq_col].astype(str).str.contains(r"NMD_transcript_variant") | \
                    groupdf.loc[:,csq_col].astype(str).str.contains(r"splice_[a-z_]+_variant") | \
                    groupdf.loc[:,csq_col].astype(str).str.contains(r"transcript_[a-z_]+tion") | \
                    groupdf.loc[:,csq_col].astype(str).str.contains(r"incomplete_terminal_codon_variant") | \
                    (groupdf.loc[:,csq_col].astype(str).str.contains(r"feature_truncation") & \
                    groupdf.loc[:,csq_col].astype(str).str.contains(r"coding_sequence_variant"))
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
    
        hc_patho_bools = groupdf["ReviewStatus"].map(cln_rs_qdict) >= 2 & \
                        groupdf["ClinicalSignificance"].astype(str).str.contains(r"[P]athogenic")
        hc_patho_num = groupdf.loc[hc_patho_bools, "#AlleleID"].drop_duplicates().size
        hc_patho_trunc_num = groupdf.loc[( hc_patho_bools & trunc_bools ), "#AlleleID"].drop_duplicates().size
        groupdf["HC_Patho_No_pertranx"] = hc_patho_num
        groupdf["HC_Patho_trunc_No_pertranx"] = hc_patho_trunc_num
        groupdf["HC_Patho_not_trunc_pertranx"] = hc_patho_num - hc_patho_trunc_num if hc_patho_num > 0 else np.nan
        tranx_df = groupdf[["TranscriptID"]].drop_duplicates()
        tranx_df["CDS_len"] = tranx_df["TranscriptID"].apply(get_cds_len_from_ensembl_tranxid)
        groupdf = groupdf.merge(tranx_df.dropna(subset=["TranscriptID"]), how="left", on="TranscriptID")
        groupdf["HC_Patho_not_trunc_freq_pertranx"] = groupdf["HC_Patho_not_trunc_pertranx"]/groupdf["CDS_len"]
        groupdf["HC_Patho_trunc_proportion"] = groupdf["HC_Patho_trunc_No_pertranx"]/groupdf["HC_Patho_No_pertranx"]
        return groupdf[["TranscriptID", 
                        "HC_Patho_No_pertranx", 
                        "HC_Patho_trunc_No_pertranx", 
                        "HC_Patho_not_trunc_pertranx", 
                        "HC_Patho_not_trunc_freq_pertranx",
                        "HC_Patho_trunc_proportion"]].drop_duplicates().dropna(subset=["TranscriptID"])
        
    
    def summary_trunc_patho_per_tranx(self):
        if not os.path.exists(self.trunc_sum):
            update = True
        else:
            if os.path.getmtime(self.trunc_sum) > os.path.getmtime(self.vep):
                logger.warning("Local truncating summary table {} is newer than the VEP annotated CLINVAR recs {} , so use the local version for now.".format(self.trunc_sum, self.vep))
                update = False
            else:
                logger.warning("Local truncating summary table {} is newer than the VEP annotated CLINVAR recs {} , Warn myself via email to update it separately but use the local version for now.".format(self.trunc_sum, self.vep))
                update = False
        
        if update:
            vep_cln_df = pd.read_table(self.path[:-4] + ".vep.tsv", low_memory=False)        
            condensed_df = vep_cln_df[[ "#AlleleID", 
                                        "TranscriptID",
                                        "Consequence",
                                        "ReviewStatus",
                                        "ClinicalSignificance"]].drop_duplicates().groupby("TranscriptID", as_index=False, dropna=False).apply(self.sum_trunc_var_per_tranx).drop_duplicates()
            tmp_out = self.trunc_sum + "." + str(uuid.uuid4())
            condensed_df.to_csv(tmp_out, sep="\t", index=False)
            result = subprocess.run("mv {} {} && ls -lh {}".format(tmp_out, self.trunc_sum, self.trunc_sum), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
            logger.info("\n{}".format(result.stdout))
        
        return self.trunc_sum
    
    
    def convertab2vcf(self, 
                      input_tab="/paedyl01/disk1/yangyxt/public_data/clinvar/tab_delimited/Reformat_updated_Clinvar.tsv"):
        vep_vcf = input_tab[:-4] + ".vep.vcf"
        if os.path.getmtime(input_tab) > os.path.getmtime(vep_vcf):
            logger.warning("Local table file {} older than the input table {}. Warn myself to update it but use the local version for now.".format(vep_vcf, input_tab))
            input_df = pd.read_table(input_tab, low_memory=False)
            logger.info("\n{}".format(input_df[:5].to_string(index=False)))
            tmp_suffix = "." + str(uuid.uuid4())
            tmp_vep_vcf = vep_vcf[:-4] + tmp_suffix + ".vcf"
            logger.info("Before converting to vcf tab:\n{}".format(input_df[:5].to_string(index=False)))
            vcf_tab = input_df.apply(self.convert2vepvcf_per_row, all_labels=",".join(input_df.columns.tolist()), axis=1).dropna(subset=["ALT"])
            logger.info("After converting to vcf tab:\n{}".format(vcf_tab[:5].to_string(index=False)))
            vcf_tab.to_csv(tmp_vep_vcf, sep="\t", index=False)
            result = subprocess.run("mv {} {} && ls -lh {}".format(tmp_vep_vcf, vep_vcf, vep_vcf), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
            logger.info("\n{}".format(result.stdout))
        else:
            logger.info("VEP vcf is already updated: {}".format(vep_vcf))
        return vep_vcf
        

    def main_convert(self, 
                     input_tab="/paedyl01/disk1/yangyxt/public_data/clinvar/tab_delimited/Reformat_updated_Clinvar.tsv", \
                     output_vep = "/paedyl01/disk1/yangyxt/public_data/clinvar/tab_delimited/Reformat_updated_Clinvar.vepoutput", \
                     email_addr = "xxchen99@hku.hk"):
        vep_vcf = input_tab[:-4] + ".vep.vcf"
        output_tab = input_tab[:-4] + ".vep.tsv"
        tmp_suffix = "." + str(uuid.uuid4())
        tmp_vep_vcf = vep_vcf[:-4] + tmp_suffix + ".vcf"
        if os.path.getmtime(input_tab) > os.path.getmtime(output_vep):
            logger.warning("Local table file {} is at least {} days old. Warn myself to update it but use the local version for now.".format(output_vep, self.update))
            if os.path.getmtime(vep_vcf) < os.path.getmtime(input_tab):
                input_df = pd.read_table(input_tab, low_memory=False)
                logger.info("\n{}".format(input_df[:5].to_string(index=False)))
                logger.info("Before converting to vcf tab:\n{}".format(input_df[:5].to_string(index=False)))
                vcf_tab = input_df.apply(self.convert2vepvcf_per_row, all_labels=",".join(input_df.columns.tolist()), axis=1).dropna(subset=["ALT"])
                logger.info("VCF tab:\n{}".format(vcf_tab[:5].to_string(index=False)))
                vcf_tab.to_csv(tmp_vep_vcf, sep="\t", index=False)
                result = subprocess.run("mv {} {} && ls -lh {}".format(tmp_vep_vcf, vep_vcf, vep_vcf), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
                logger.info("\n{}".format(result.stdout))
            else:
                logger.info("{} is newer than {}. We dont have to update the VEP vcf file".format(vep_vcf, input_tab))
            
            vep_cmd = "cp -f {vepvcf} {vep} && \
                       bash /paedyl01/disk1/yangyxt/ngs_scripts/annovar_filtration_per_family.sh vep_onestop_annotation {vep} {mode} && \
                       ls -lh {vep_op} && rm {vep} && \
                       mv {vep_op} {vep_fop} && ls -lh {vep_fop}".format(vepvcf = vep_vcf, \
                                                                         vep=tmp_vep_vcf, \
                                                                         mode="", \
                                                                         vep_op = tmp_vep_vcf[:-4] + ".vepoutput", \
                                                                         vep_fop = output_vep)
            email_cmd = '''
                        echo "{vepcmd}" | mail -s "WARNING: VEP annotated Clinvar table already 6 months old. Execute cmd below to update." {email} - 
                        '''.format(vepcmd = vep_cmd, email=email_addr)
            # result = subprocess.run(email_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, encoding="utf-8")
            # logger.info("\n{}".format(result.stdout))             
            return self.merge_with_clinvar(output_vep, input_tab, output_tab)
        else:
            return output_tab
    


    