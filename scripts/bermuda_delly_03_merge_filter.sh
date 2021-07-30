#!/bin/bash
#SBATCH -A snic2020-5-504
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00

module load bioinfo-tools
module load delly
module load bcftools

set -eu

## Type of analysis (DEL INS INV DUP BND)

TYPE=$1

if [ ! -d delly/merged ]; then
    mkdir delly/merged
fi

bcftools merge \
    -m id \
    -O b \
    -o delly/merged/delly_merged_raw_${TYPE}.bcf \
    $(ls delly/genotyped/*_${TYPE}.bcf)
    
bcftools index delly/merged/delly_merged_raw_${TYPE}.bcf

delly filter \
    -f germline \
    -t ${TYPE} \
    -o delly/merged/delly_merged_filtered_${TYPE}.bcf \
    delly/merged/delly_merged_raw_${TYPE}.bcf
    