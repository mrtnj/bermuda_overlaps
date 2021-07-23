#!/bin/bash
#SBATCH -A snic2020-5-504
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 1-00:00:00
#SBATCH -J bermuda_coverage


module load bioinfo-tools
module load BEDTools

set -eu

if [ ! -d coverage ]; then
    mkdir coverage
fi

bedtools makewindows \
    -g galGal4_chr_sizes.txt \
    -w 10000 > \
    10kbp_windows.bed

for FILE in /proj/sllstore2017078/nobackup/bermuda_bam/*.bam;
do
    echo $FILE
    
    SAMPLE=${FILE##*/};
    SAMPLE=${SAMPLE%_recal.bam}
    
    echo $SAMPLE
    
    bedtools coverage \
        -a 10kbp_windows.bed \
        -b $FILE \
        -sorted \
        > coverage/${SAMPLE}_coverage.txt
    
done