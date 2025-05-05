import pandas as pd
import pysam
import logging
import argparse as ap
import base64


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
console_handler=logging.StreamHandler()
console_handler.setLevel(logging.INFO)
formatter = logging.Formatter("%(levelname)s:%(asctime)s:%(funcName)s:%(lineno)s:%(message)s")
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


# Now we need to generate a modularized script to add variant ID on the ID field for each variant record parsed by pysam
# The variant ID should be composed by chromosome:pos:ref_allele:alt_alleles, the string should be converted
# There should be two version of notation supported.
# One is the indel notation, the other is the symbolic notation. For the symbolic notation, we need to extract the END coordinate in the INFO field and compose the ID str like this: chromosome:pos-end:ref_allele:alt_allele to ensure uniqueness of the variant allele
# The variant representing str should not be directly used as the ID sequence. It can be hashed into a base64 str sequence and to make the ID reversely trackable.


def variant_to_id(chrom, pos, ref, alt, end=None):
    """
    Convert variant information to a unique, reversible ID.

    Args:
    chrom (str): Chromosome
    pos (int): Start position
    ref (str): Reference allele
    alt (str): Alternate allele
    end (int, optional): End position for symbolic alleles

    Returns:
    str: Unique, reversible variant ID
    """
    if end is not None:
        var_str = f"{chrom}:{pos}-{end}:{ref}:{alt}"
    else:
        var_str = f"{chrom}:{pos}:{ref}:{alt}"

    # Encode the variant string to base64
    encoded = base64.urlsafe_b64encode(var_str.encode()).decode()

    return encoded


def id_to_variant(var_id):
    """
    Convert a variant ID back to its original representation.

    Args:
    var_id (str): Unique variant ID

    Returns:
    tuple: (chrom, pos, ref, alt, end)
    """
    # Decode the base64 string
    decoded = base64.urlsafe_b64decode(var_id).decode()

    # Parse the decoded string
    parts = decoded.split(':')
    chrom = parts[0]

    if '-' in parts[1]:
        pos, end = map(int, parts[1].split('-'))
    else:
        pos = int(parts[1])
        end = None

    ref = parts[2]
    alt = parts[3]

    return chrom, pos, ref, alt, end


def add_variant_ids_to_vcf(input_vcf, output_vcf):
    """
    Add variant IDs to a VCF file.

    Args:
    input_vcf (str): Path to input VCF file
    output_vcf (str): Path to output VCF file
    """
    with pysam.VariantFile(input_vcf, 'r') as vcf_in:
        with pysam.VariantFile(output_vcf, 'w', header=vcf_in.header) as vcf_out:
            for record in vcf_in:
                chrom = record.chrom
                pos = record.pos
                ref = record.ref

                # Check if variant is structural by looking for SV-specific INFO fields
                is_sv = (
                    ('SVTYPE' in record.info and 'END' in record.info) or
                    any(alt.startswith('<') and alt.endswith('>') for alt in record.alts)
                )

                for alt in record.alts:
                    end = record.info.get('END', record.stop) if is_sv else None
                    var_id = variant_to_id(chrom, pos, ref, alt, end)
                    record.id = var_id
                vcf_out.write(record)


def decode_variant_ids_from_vcf(input_vcf, output_vcf):
    """
    Decode variant IDs from a VCF file.
    """
    with pysam.VariantFile(input_vcf, 'r') as vcf_in:
        with pysam.VariantFile(output_vcf, 'w', header=vcf_in.header) as vcf_out:
            for record in vcf_in:
                var_id = record.id
                chrom, pos, ref, alt, end = id_to_variant(var_id)
                literal_id =  chrom + ":" + str(pos) + "-" + str(end) + ":" + ref + "->" + alt
                record.id = literal_id
                vcf_out.write(record)
    


if __name__ == '__main__':
    parser = ap.ArgumentParser()
    parser.add_argument("-m", "--mode", type=str, required=True, help="The mode to run the script, either 'encode' or 'decode'")
    parser.add_argument("-i", "--input_vcf", type=str, help="The path to the VCF of a family containing patients", required=True)
    parser.add_argument("-o", "--output_vcf", type=str, help="The path to the output VCF file", required=False, default="")

    args=parser.parse_args()

    if args.mode == "encode":
        if len(args.output_vcf) > 0:
            output_vcf = args.output_vcf
        else:
            output_vcf = args.input_vcf.replace(".vcf", ".uniqID.vcf")
        add_variant_ids_to_vcf(args.input_vcf, output_vcf)
    elif args.mode == "decode":
        if len(args.output_vcf) > 0:
            output_vcf = args.output_vcf
        else:
            output_vcf = args.input_vcf.replace(".vcf", ".decodeID.vcf")
        decode_variant_ids_from_vcf(args.input_vcf, output_vcf)
    else:
        parser.print_help()
        sys.exit(1)





