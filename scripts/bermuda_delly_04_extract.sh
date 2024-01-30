#!/bin/bash
#SBATCH -A snic2020-5-504
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00

module load bioinfo-tools
module load bcftools

set -eu

for TYPE in DEL INS INV DUP; do

    ## Total stats

    bcftools stats \
        delly/merged/delly_merged_filtered_${TYPE}.bcf \
        > delly/merged/delly_stats_filtered_${TYPE}.txt

    bcftools stats \
        delly/merged/delly_merged_raw_${TYPE}.bcf \
        > delly/merged/delly_stats_raw_${TYPE}.txt

    ## Extract all

    bcftools query \
        -f '%CHROM\t%POS\t%INFO\t%QUAL\t%FILTER[\t%GT]\n' \
        -o delly/merged/delly_merged_raw_${TYPE}.txt \
        delly/merged/delly_merged_raw_${TYPE}.bcf 
        
    bcftools query \
        -f '%CHROM\t%POS\t%INFO\t%QUAL\t%FILTER[\t%GT]\n' \
        -o delly/merged/delly_merged_filtered_${TYPE}.txt \
        delly/merged/delly_merged_filtered_${TYPE}.bcf 
        
done