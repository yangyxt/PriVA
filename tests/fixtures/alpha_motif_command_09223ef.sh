    python ${pick_intolerant_motifs_py} \
    --vcf_path ${alphamissense_vep_vcf} \
    --rmc_tsv ${SCRIPT_DIR}/../data/gnomAD_RMC/gnomad_v2.1.1_rmc_flat.${assembly}.tsv \
    --hcseeker_tsv ${SCRIPT_DIR}/../data/HCSeeker/HC_spots.${assembly}.tsv \
    --output ${alphamissense_intolerant_motifs} \
    --processes ${threads}