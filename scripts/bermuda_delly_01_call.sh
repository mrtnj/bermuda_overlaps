#!/bin/bash
#SBATCH -A snic2020-5-504
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2-00:00:00

module load bioinfo-tools
module load delly

set -eu

## Type of analysis (DEL INS INV DUP BND)

TYPE=$1

if [ ! -d delly ]; then
    mkdir delly
fi

for FILE in /proj/sllstore2017078/nobackup/bermuda_bam/*.bam;
do
    echo $FILE
    
    SAMPLE=${FILE##*/};
    SAMPLE=${SAMPLE%_recal.bam}
    
    echo $SAMPLE
    echo $TYPE 
    
    if [ ! -f delly/${SAMPLE}_delly_${TYPE}.bcf ]; then

        delly call \
            -g /proj/sllstore2017078/private/martinj/galGal4_validated.fa \
            -t $TYPE \
            --noindels \
            -o delly/${SAMPLE}_delly_${TYPE}.bcf \
            $FILE 
            
    fi
    
done

