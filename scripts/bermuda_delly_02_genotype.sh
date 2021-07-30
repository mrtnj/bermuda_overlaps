#!/bin/bash
#SBATCH -A snic2020-5-504
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00

module load bioinfo-tools
module load delly

set -eu

## Type of analysis (DEL INS INV DUP BND)

TYPE=$1

delly merge \
    -t ${TYPE} \
    -o delly/delly_merged_${TYPE}.bcf \
    $(ls delly/P4806_*delly_${TYPE}.bcf)

if [ ! -d delly/genotyped ]; then
    mkdir delly/genotyped
fi

for FILE in /proj/sllstore2017078/nobackup/bermuda_bam/*.bam;
do
    echo $FILE
    
    SAMPLE=${FILE##*/};
    SAMPLE=${SAMPLE%_recal.bam}
    
    echo $SAMPLE
    
    if [ ! -f delly/genotyped/${SAMPLE}_delly_genotyped_${TYPE}.bcf ]; then
    
        delly call \
            -t ${TYPE} \
            -g /proj/sllstore2017078/private/martinj/galGal4_validated.fa \
            $FILE \
            -v delly/delly_merged_${TYPE}.bcf \
            -o delly/genotyped/${SAMPLE}_delly_genotyped_${TYPE}.bcf
            
    fi 
    
done

