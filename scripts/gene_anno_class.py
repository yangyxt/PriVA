import pandas as pd
import numpy as np
import subprocess
import sqlite3
import json
import requests
import logging
import gffutils
import os
import time
import argparse as ap
from python_utils import create_gffutils_db, check_remote_file_exists


ncbi2ucsc_as_dict = {"GRCh38": "hg38", 
					"hg38": "hg38",
					"GRCh37": "hg19",
					"hg19": "hg19"}


class GENE_ANNO(object):
    def __init__(self, 
                 wkd="/home/yangyxt/public_data/gene_annotation/{assembly}",
                 assembly="GRCh37",
                 source="NCBI",
                 ncbi_urls = {"GRCh37": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz", \
                     		  "GRCh38": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz"}, 
                 gencode_urls = {"GRCh37": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/GRCh37_mapping/gencode.v{release_ver}lift37.basic.annotation.gff3.gz", \
                     			 "GRCh38": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v{release_ver}.basic.annotation.gff3.gz"}, 
                 ensembl_urls = {"GRCh37": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/GRCh37_mapping/gencode.v{release_ver}lift37.basic.annotation.gff3.gz", \
                     			 "GRCh38": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v{release_ver}.basic.annotation.gff3.gz"}):
        self.wkd = wkd.format(assembly=ncbi2ucsc_as_dict[assembly])
        self.gff_gz=os.path.join(wkd,"{}.{}.gff.gz".format(source, ncbi2ucsc_as_dict[assembly]))
        self.gff_db=os.path.join(wkd,"{}.{}.gff.db".format(source, ncbi2ucsc_as_dict[assembly]))
        
        self.logger = logging.getLogger(__name__)
		self.setup_logger()
        for k,v in gencode_urls.items():
            release_ver = 44
            while not check_remote_file_exists(v.format(release_ver=release_ver)):
                release_ver += 1
			gencode_urls[k] = v.format(release_ver=release_ver)
			ensembl_urls[k] = v.format(release_ver=release_ver)
		if re.search(r"[Nn][Cc][Bb][Ii]", source):
			self.gff_url = ncbi_urls[assembly]
		elif re.search(r"[Gg][Ee][Nn][Cc][Oo][Dd][Ee]", source):
			self.gff_url = gencode_urls[assembly]
		else:
			self.gff_url = ensembl_urls[assembly]
   

	 def setup_logger(self):
        # Configure logger
        self.logger.setLevel(logging.INFO)
        
        # Create a console handler
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        
        # Create a formatter and set the formatter for the handler
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        
        # Add the handlers to the logger
        self.logger.addHandler(ch)
        

    def access_gff_gz(self):
        if os.path.exists(self.gff_gz):
            file_time = os.path.getmtime(self.gff_gz)
            current_time = time.time()
            file_age_days = (current_time - file_time) / (60 * 60 * 24)
            
            if file_age_days > 30:
                update = True
            else:
                update = False
        else:
            update = True
            
        if update:
			response = requests.get(self.gff_url)
			if response.status_code == 200:
				with open(self.gff_gz, 'wb') as f:
					f.write(response.content)
			else:
				raise Exception(f"Cannot download gff file from {self.gff_url}")
		return self.gff_gz


	def access_gff_db(self):
		return create_gffutils_db(gff_file = self.access_gff_gz(), db_file = self.gff_db)





	
        

        