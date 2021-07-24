#!/bin/bash
#SBATCH -A snic2020-5-504
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00

module load bioinfo-tools
module load delly
module load bcftools

set -eu

if [ ! -d delly/merged ]; then
    mkdir delly/merged
fi

bcftools merge \
    -m id \
    -O b \
    -o delly/merged/delly_merged_raw.bcf \
    $(ls delly/genotyped/*.bcf)

delly filter \
    -f germline \
    -o delly/mergred/delly_merged_filtered.bcf \
    delly/merged/delly_merged_raw.bcf
    