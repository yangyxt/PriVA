import pandas as pd
import argparse as ap
import logging
import pysam
import pickle
import gzip

from combine_annotations import convert_vcf_to_tab


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(module)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


def offer_clinvar_alleleid_to_transcript_map(clinvar_tab: str = "/paedyl01/disk1/yangyxt/public_data/clinvar/tab_delimited/variant_summary.txt.gz",
                                            mane_tranx_tab: str = "/paedyl01/disk1/yangyxt/PriVA/data/MANE/MANE.GRCh38.v1.4.summary.txt.gz"):

    # URL for updated clinvar tab: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
    # URL for updated MANE tab: https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.summary.txt.gz

    clinvar_df = pd.read_table(clinvar_tab, low_memory=False)
    hg38_clinvar_df = clinvar_df[clinvar_df["Assembly"] == "GRCh38"]
    hg38_clinvar_df["Tranx_RefSeq_ID"] = hg38_clinvar_df["Name"].str.extract(r'(NM_\d+\.\d+)')[0]

    mane_tranx_df = pd.read_table(mane_tranx_tab, low_memory=False)

    refseq_to_ensembl = dict(zip(mane_tranx_df["RefSeq_nuc"], mane_tranx_df["Ensembl_nuc"]))

    hg38_clinvar_df["Tranx_Ensembl_ID"] = hg38_clinvar_df["Tranx_RefSeq_ID"].map(refseq_to_ensembl)
    allele_id_to_tranx = dict(zip(hg38_clinvar_df["#AlleleID"], hg38_clinvar_df["Tranx_Ensembl_ID"].str.split(".", expand=True).iloc[:, 0]))
    return allele_id_to_tranx


def extract_allele_id_map(clinvar_vcf: str):
    """
    Create dictionaries mapping variant IDs to ClinVar annotations
    
    Args:
        clinvar_vcf: Path to the ClinVar VCF file
        
    Returns:
        Three dictionaries mapping variant_id to:
        - CLNDN (disease names)
        - CLNSIG (clinical significance)
        - CLNREVSTAT (review status)
    """
    alleleid_map = {}   # Allele ID
    
    # Open VCF file
    vcf = pysam.VariantFile(clinvar_vcf)
    
    # Process each record
    for record in vcf:
        # Create variant ID
        variant_id = f"{record.chrom}_{record.pos}_{record.ref}_{record.alts[0]}"

        if 'ALLELEID' in record.info:
            alleleid_map[variant_id] = record.info['ALLELEID']
    
    return alleleid_map


def analyze_splice_event(event_str: str, spliceai_delta_score: float) -> dict:
    """
    Analyze a single SpliceVault event to determine if it affects intolerant domains.
    
    Args:
        event_str (str): String describing splice event (e.g. "ES:2-3;12%;Frameshift")
        transcript_id (str): Ensembl transcript ID
        intron_pos (str): String indicating which intron is affected ("current/total")
        exon_pos (str): String indicating which exons are involved ("start-end")
        transcript_domain_map (dict): Mapping of transcripts to exon-domain relationships
                                      {transcript_id: {exon_number: [domain_paths]}}
        intolerant_domains (set): Set of domain paths that are intolerant to structural changes
        
    Returns:
        dict: Dictionary containing:
            - affected_exons (list): List of exon numbers affected
            - affected_domains (list): List of domain paths affected
            - is_lof (bool): Boolean indicating if event likely causes LoF
            - reason (str): String explaining LoF classification
    """
    if event_str in ["NA", "None", "none", "None", "NONE", "Null", "null", "Null", "NULL", "nan", "NaN", "NAN", "N/A"]:
        return None

    if spliceai_delta_score < 0.3:
        # Alt splicing Not necessarily happening
        return None

    stat_dict = {"large_AA_change": 0,
                 "small_AA_change": 0,
                 "frameshift": 0}

    for event in event_str.split('&'):
        # Parse event string
        try:
            rank, event_type, pos_str, freq, frame_impact = event.split(':')  # e.g., "Top1:CD:-10:2%:Frameshift&Top2:CD:-4:0.6%:Frameshift&Top3:CD:+491:0.05%:Frameshift&Top4:CD:+21:0.006%:Frameshift"
        except ValueError as ve:
            raise ValueError(f"Error splitting event string {event_str}: {ve}")
        # We need to translate 2% to 0.02, 0.6% to 0.006, 0.05% to 0.00005, 0.006% to 0.000006
        if '%' in freq:
            freq = float(freq.rstrip('%')) / 100
        elif "#" in freq:
            freq = float(freq.lstrip('#')) / 100
        else:
            raise ValueError(f"The frequency format is not recognized: {freq}, the event string is {event_str}")


        # Process event based on type
        if event_type == 'ES':
            if frame_impact == "Frameshift":
                stat_dict["frameshift"] += freq
            else:
                stat_dict["large_AA_change"] += freq
        elif event_type == 'CD' or event_type == 'CA':
            offset = int(pos_str.lstrip('+-'))
            # Upstream cryptic donor (in exon)
            if frame_impact == "Frameshift":
                stat_dict["frameshift"] += freq
            elif abs(offset) > 15:  # Ignore very close cryptic sites
                stat_dict["large_AA_change"] += freq
            elif abs(offset) >= 3 and abs(offset) <= 9:
                stat_dict["small_AA_change"] += freq

    return stat_dict


def per_gene_stat(gene_df: pd.DataFrame):
    """
    Calculate statistics for each gene, only keep high review status variants
    
    Args:
        gene_df: DataFrame containing gene information
        
    Returns:
        A dictionary with the following keys:
        - "gene_symbol": gene symbol
        - "ensembl_id": ensembl gene id
        - "clinvar_variants": ids of variants in the gene
        - "clinvar_pathogenic": ids of pathogenic/likely pathogenic variants in the gene
        - "clinvar_not_pathogenic": ids of benign/likely benign/vous variants in the gene
        - "putative_nmd_variants": ids of variants that are putative NMD variants in the gene
        - "small_aachange": ids of variants with small amino acid changes in the gene
        - "large_aachange": ids of variants with large amino acid changes in the gene
        - "no_aachange": ids of variants with no amino acid changes in the gene
        - "transcripts": { clinvar corresponding transcripts (enst id, must be protein coding)
                            "transcript1": {
                                            "clinvar_variants": ids of variants in the gene
                                            "clinvar_pathogenic": ids of pathogenic/likely pathogenic variants in the gene
                                            "clinvar_not_pathogenic": ids of benign/likely benign/vous variants in the gene
                                            "putative_nmd_variants": ids of variants that are putative NMD variants in the gene
                                            "small_aachange": ids of variants with small amino acid changes in the gene
                                            "large_aachange": ids of variants with large amino acid changes in the gene
                                            "no_aachange": ids of variants with no amino acid changes in the gene
                                            }
                            ...
                            "transcriptN": {
                                        ...
                                        }
                          }
        - "domains": "domain1": { 
                                  "clinvar_variants": ids of variants in the domain, 
                                  "clinvar_pathogenic": ids of pathogenic/likely pathogenic variants in the domain, 
                                  "clinvar_not_pathogenic": ids of benign/likely benign/vous variants in the domain, 
                                  "putative_nmd_variants": ids of variants that are putative NMD variants in the domain, 
                                  "small_aachange": ids of variants with small amino acid changes in the domain, 
                                  "large_aachange": ids of variants with large amino acid changes in the domain, 
                                  "no_aachange": ids of variants with no amino acid changes in the domain 
                                }
                     ...
                     "domainN": {
                                  ...
                                }   
    """
    high_review_status = {
        'practice_guideline',                                    # 4 stars
        'reviewed_by_expert_panel',                              # 3 stars
        'criteria_provided,_multiple_submitters,_no_conflicts',  # 2 stars
    }
    gene_symbol = gene_df["SYMBOL"].iloc[0]
    ensembl_gene_id = gene_df["Gene"].iloc[0]
    logger.info(f"The gene_df looks like this: {gene_df.head().to_string(index=False)}")
    gene_df = gene_df.loc[gene_df["BIOTYPE"] == "protein_coding", :]
    gene_df = gene_df.loc[gene_df["CLNREVSTAT"].isin(high_review_status), :]
    if len(gene_df) == 0:
        logger.warning(f"No high quality (CLNREVSTAT >= 2 stars) protein coding transcripts found for {gene_symbol} ({ensembl_gene_id})")
        return None
    protein_coding_transcripts = set(gene_df["Feature"].unique().tolist())
    all_covered_domains = gene_df["DOMAINS"].drop_duplicates().dropna().str.split("&", expand=True).stack().dropna()
    all_covered_domains = set([ensembl_gene_id + ":" + str(domain) for domain in all_covered_domains.unique().tolist() if "nan" not in domain])
    gene_stat_dict = {"gene_symbol": gene_symbol,
                      "ensembl_id": ensembl_gene_id,
                      "clinvar_variants": set(gene_df["variant_id"].unique().tolist()),
                      "clinvar_pathogenic": set(),
                      "clinvar_not_pathogenic": set(),
                      "putative_nmd_variants": set(),
                      "small_aachange": set(),
                      "large_aachange": set(),
                      "no_aachange": set(),
                      "transcripts": {tranx: {"clinvar_variants": set(),
                                              "clinvar_pathogenic": set(),
                                              "clinvar_not_pathogenic": set(),
                                              "putative_nmd_variants": set(),
                                              "small_aachange": set(),
                                              "large_aachange": set(),
                                              "no_aachange": set()} for tranx in protein_coding_transcripts},
                      "domains": {domain: {"clinvar_variants": set(),
                                           "clinvar_pathogenic": set(),
                                           "clinvar_not_pathogenic": set(),
                                           "putative_nmd_variants": set(),
                                           "small_aachange": set(),
                                           "large_aachange": set(),
                                           "no_aachange": set()} for domain in all_covered_domains},
                      }
    
    for i in range(len(gene_df)):
        row = gene_df.iloc[i]
        variant_id = str(row["variant_id"])
        ensembl_gene_id = str(row["Gene"])
        transcript = str(row["Feature"])
        domains = [ensembl_gene_id + ":" + domain for domain in str(row["DOMAINS"]).split("&") if "nan" not in domain]
        if transcript not in protein_coding_transcripts:
            continue
        gene_stat_dict["transcripts"][transcript]["clinvar_variants"].add(variant_id)
        for domain in domains:
            gene_stat_dict["domains"][domain]["clinvar_variants"].add(variant_id)
        # Then we have to determine whether the variant is pathogenic and what kind of aa change it brings to the transcript
        if "athogenic" in str(row["CLNSIG"]):
            gene_stat_dict["clinvar_pathogenic"].add(variant_id)
            gene_stat_dict["transcripts"][transcript]["clinvar_pathogenic"].add(variant_id)
            for domain in domains:
                gene_stat_dict["domains"][domain]["clinvar_pathogenic"].add(variant_id)
        else:
            gene_stat_dict["clinvar_not_pathogenic"].add(variant_id)
            gene_stat_dict["transcripts"][transcript]["clinvar_not_pathogenic"].add(variant_id)
            for domain in domains:
                gene_stat_dict["domains"][domain]["clinvar_not_pathogenic"].add(variant_id)

        # Use SpliceAI and SpliceVault to identify splicing brought AA changes
        spliceai_delta = float(row["SpliceVault_SpliceAI_delta"])
        splicing_stat = analyze_splice_event(str(row["SpliceVault_top_events"]), spliceai_delta)
    
        # Identify putative NMD variants
        vep_consq_lof = (("HC" in str(row["LoF"])) or ("OS" in str(row["LoF"]))) and any(csq in str(row["Consequence"]) for csq in ["stop_gained", "start_lost", "frameshift"]) and (not "escape" in str(row["NMD"]))
        splic_donor_loss = float(row["SpliceAI_pred_DS_DL"]) > 0.5
        putative_nmd = vep_consq_lof or splic_donor_loss
        if putative_nmd:
            gene_stat_dict["putative_nmd_variants"].add(variant_id)
            gene_stat_dict["transcripts"][transcript]["putative_nmd_variants"].add(variant_id)
            for domain in domains:
                gene_stat_dict["domains"][domain]["putative_nmd_variants"].add(variant_id)

        # Identify large aa change variants (in-frame indels, exon skipping, upstream cryptic donor (deletion), downstream cryptic acceptor (elongation))
        large_inframe_indel = any(csq in str(row["Consequence"]) for csq in ["inframe_deletion", "inframe_insertion"]) and abs(len(row["ref"]) - len(row["alt"])) >= 15
        splicing_large_aa_change = False
        splicing_small_aa_change = False
        if splicing_stat is not None:
            if (splicing_stat["frameshift"] + splicing_stat["large_AA_change"]) >= 0.2 and spliceai_delta >= 0.5:
                splicing_large_aa_change = True
            elif (splicing_stat["frameshift"] + splicing_stat["large_AA_change"]) >= 0.4 and spliceai_delta >= 0.3:
                splicing_large_aa_change = True
            splicing_small_aa_change = (splicing_stat["small_AA_change"] >= 0.2 and spliceai_delta >= 0.3) or (splicing_stat["small_AA_change"] >= 0.4 and spliceai_delta >= 0.1)

        if (large_inframe_indel or splicing_large_aa_change) and (not putative_nmd):
            gene_stat_dict["large_aachange"].add(variant_id)
            gene_stat_dict["transcripts"][transcript]["large_aachange"].add(variant_id)
            for domain in domains:
                gene_stat_dict["domains"][domain]["large_aachange"].add(variant_id)

        # Identify small aa change variants (in-frame indels, exon skipping, upstream cryptic donor (deletion), downstream cryptic acceptor (elongation) and missense variants)
        missense = any(csq in str(row["Consequence"]) for csq in ["missense_variant"])
        small_inframe_indel = any(csq in str(row["Consequence"]) for csq in ["inframe_deletion", "inframe_insertion"]) and abs(len(row["ref"]) - len(row["alt"])) <= 9
        if (small_inframe_indel or splicing_small_aa_change or missense) and (not putative_nmd) and (not large_inframe_indel) and (not splicing_large_aa_change):
            gene_stat_dict["small_aachange"].add(variant_id)
            gene_stat_dict["transcripts"][transcript]["small_aachange"].add(variant_id)
            for domain in domains:
                gene_stat_dict["domains"][domain]["small_aachange"].add(variant_id)

        # Identify no aa change variants (synonymous variants)
        synonymous = any(csq in str(row["Consequence"]) for csq in ["synonymous_variant"])
        if synonymous and (not putative_nmd) and (not large_inframe_indel) and (not splicing_large_aa_change) and (not splicing_small_aa_change) and (not missense) and (not small_inframe_indel):
            gene_stat_dict["no_aachange"].add(variant_id)
            gene_stat_dict["transcripts"][transcript]["no_aachange"].add(variant_id)
            for domain in domains:
                gene_stat_dict["domains"][domain]["no_aachange"].add(variant_id)

    return gene_stat_dict



def main_stat_variants( clinvar_anno_vcf: str, 
                        output_pickle: str,
                        threads=10 ):
    combined_df = convert_vcf_to_tab(clinvar_anno_vcf, threads)
    combined_df["variant_id"] = combined_df["chrom"] + "_" + combined_df["pos"].astype(str) + "_" + combined_df["ref"] + "_" + combined_df["alt"]

    # First we need to identify the corresponding transcript for each variant in ClinVar
    allele_id_to_tranx = offer_clinvar_alleleid_to_transcript_map()
    var_id_to_allele_id = extract_allele_id_map(clinvar_anno_vcf)

    combined_df["allele_id"] = combined_df["variant_id"].map(var_id_to_allele_id)
    combined_df["target_transcript"] = combined_df["allele_id"].map(allele_id_to_tranx)

    refined_df = combined_df.loc[combined_df["target_transcript"] == combined_df["Feature"], :]

    groupby_gene = refined_df.groupby("Gene", as_index=False)
    total_stat_dict = {}
    for gene, gene_df in groupby_gene:
        gene_stat_dict = per_gene_stat(gene_df)
        if gene_stat_dict is not None:
            logger.info(f"Processed {gene} with {len(gene_df)} variants")
            total_stat_dict[gene] = gene_stat_dict

    logger.info(f"Total number of genes processed: {len(total_stat_dict)}")
    if not output_pickle.endswith(".gz"):
        output_pickle = output_pickle + ".gz"
    with gzip.open(output_pickle, "wb", compresslevel=9) as f:
        pickle.dump(total_stat_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
    return total_stat_dict


if __name__ == "__main__":
    parser = ap.ArgumentParser()
    parser.add_argument("-v", "--clinvar_anno_vcf", type=str, required=True)
    parser.add_argument("-o", "--output_pickle", type=str, required=True)
    parser.add_argument("-t", "--threads", type=int, default=10)
    args = parser.parse_args()

    # This script needs to first annotate the updated GRCh38 ClinVar VCF with annotation_vcf.sh to generate the annotated vcf
    main_stat_variants(args.clinvar_anno_vcf, args.output_pickle, args.threads)

