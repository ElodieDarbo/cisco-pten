#!/bin/bash

module load hicpro/3.1.0
source ~/.bashrc

#HiC-Pro mapping
HiC-Pro -i [FASTQ_INPUT_FOLDER]
        -o [Tiled-C_FOLDER] \
        -c Tiled-C_config_file \
        -s mapping \
        -s quality_checks

#HiC-Pro proc-hic
HiC-Pro -i [Tiled-C_FOLDER]/bowtie_results/bwt2/ \
        -o [Tiled-C_FOLDER] \
        -c Tiled-C_config_file \
        -s proc_hic \
        -s quality_checks

#HiC-Pro merge_sample
HiC-Pro -i [Tiled-C_FOLDER]/hic_results/data/ \
        -o [Tiled-C_FOLDER] \
        -c Tiled-C_config_file \
        -s merge_persample \
        -s quality_checks

#HiC-Pro matrix       
HiC-Pro -i [Tiled-C_FOLDER]/hic_results/data/ \
        -o [Tiled-C_FOLDER] \
        -c Tiled-C_config_file \
        -s build_contact_maps 
        -s ice_norm
