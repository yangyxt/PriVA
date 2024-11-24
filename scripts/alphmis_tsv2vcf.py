import gzip
import subprocess
import sys
import os


def tsv_to_vcf(input_tsv, output_vcf):
    """
    Convert AlphaMissense TSV file to VCF format
    
    Parameters:
    input_tsv (str): Path to input TSV file
    output_vcf (str): Path to output VCF file
    """
    # VCF header
    vcf_header = "\n".join([
        "##fileformat=VCFv4.2",
        "##INFO=<ID=UNIPROT,Number=1,Type=String,Description=\"UniProt ID\">",
        "##INFO=<ID=TRANSCRIPT,Number=1,Type=String,Description=\"Ensembl Transcript ID\">",
        "##INFO=<ID=PVAR,Number=1,Type=String,Description=\"Protein Variant\">",
        "##INFO=<ID=AM_PATHOGENICITY,Number=1,Type=Float,Description=\"Pathogenicity Score\">",
        "##INFO=<ID=AM_CLASS,Number=1,Type=String,Description=\"Classification\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    ])
    
    if not output_vcf.endswith(".gz"):
        raise ValueError("Output VCF file must end with .gz")
        
    tmp_output = output_vcf + ".tmp.vcf"
    # Open gzip file in text mode with explicit encoding
    with gzip.open(input_tsv, 'rt', encoding='utf-8') as fin, open(tmp_output, 'w') as fout:
        # Write VCF header
        fout.write(vcf_header + '\n')
        
        # Skip header in TSV
        header = fin.readline()
        
        # Process each line
        for line in fin:
            if line.startswith('#'):  # Skip comment lines
                continue
                
            fields = line.strip().split()
            chrom = fields[0]
            pos = fields[1]
            ref = fields[2]
            alt = fields[3]
            uniprot = fields[5]
            transcript = fields[6]
            pvar = fields[7]
            pathogenicity = fields[8]
            classification = fields[9]
            
            # Create INFO field
            info = f"UNIPROT={uniprot};TRANSCRIPT={transcript};PVAR={pvar};AM_PATHOGENICITY={pathogenicity};AM_CLASS={classification}"
            
            # Write VCF line
            vcf_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}"
            fout.write(vcf_line + '\n')

    cmd = f"bgzip -f -c {tmp_output} > {output_vcf} && tabix -p vcf -f {output_vcf}"
    subprocess.run(cmd, shell=True)
    
    os.remove(tmp_output)



if __name__ == "__main__":
    tsv_to_vcf(sys.argv[1], sys.argv[2])