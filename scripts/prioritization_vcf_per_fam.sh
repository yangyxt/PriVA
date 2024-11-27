#!/usr/bin/env bash

# This script is used to prioritize the variants in the VCF file for each family
# The script expect a Family VCF input along with a pedigree table contains the Family pedigree information.
# Notice that the family sample ID should be consistent between the pedigree table and the VCF file.
# Notice that the proband of the family should stay as the first record of the family in the pedigree table

# Maintainer: yangyxt@gmail.com, yangyxt@hku.hk

SCRIPT_DIR=$(dirname $(readlink -f $0))
source ${SCRIPT_DIR}/common_bash_utils.sh

# In this script, we have several steps:
# 1. Convert the annotated VCF to a Table with transcript specific annotations as rows
# 2. Assign the ACMG criterias to each variant-transcript annotation record
# 3. Prioritize the variants based on the ACMG criterias


function prepare_combined_tab () {
    local input_vcf=${1}  # This input VCF file should be annotated with many things
    local CADD_anno_file=${2}
    local config_file=${3}

    



}



