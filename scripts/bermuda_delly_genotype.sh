#!/bin/bash
#SBATCH -A snic2020-5-504
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00

module load bioinfo-tools
module load delly

set -eu

delly merge \
    -o delly/delly_merged.bcf \
    $(ls delly/P4806*.bcf)

if [ ! -d delly/genotyped ]; then
    mkdir delly/genotyped
fi

for FILE in /proj/sllstore2017078/nobackup/bermuda_bam/*.bam;
do
    echo $FILE
    
    SAMPLE=${FILE##*/};
    SAMPLE=${SAMPLE%_recal.bam}
    
    echo $SAMPLE
    
    if [ ! -f delly/genotyped/${SAMPLE}_delly_genotyped.bcf ]; then
    
        delly call \
            -g /proj/sllstore2017078/private/martinj/galGal4_validated.fa \
            $FILE \
            -v delly/delly_merged.bcf \
            -o delly/genotyped/${SAMPLE}_delly_genotyped.bcf
            
    fi 
    
done

