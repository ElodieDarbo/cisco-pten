#!/bin/bash

module load hicpro/3.1.0
source ~/.bashrc

#HiC-Pro mapping
HiC-Pro -i [FASTQ_INPUT_FOLDER]
        -o [HI-C_FOLDER] \
        -c Hi-C_config_file \
        -s mapping \
        -s quality_checks

#HiC-Pro proc-hic
HiC-Pro -i [HI-C_FOLDER]/bowtie_results/bwt2/ \
        -o [HI-C_FOLDER] \
        -c Hi-C_config_file \
        -s proc_hic \
        -s quality_checks

#HiC-Pro merge_sample
HiC-Pro -i [HI-C_FOLDER]/hic_results/data/ \
        -o [HI-C_FOLDER] \
        -c Hi-C_config_file \
        -s merge_persample \
        -s quality_checks

#HiC-Pro matrix       
HiC-Pro -i [HI-C_FOLDER]/hic_results/data/ \
        -o [HI-C_FOLDER] \
        -c Hi-C_config_file \
        -s build_contact_maps 
        -s ice_norm
