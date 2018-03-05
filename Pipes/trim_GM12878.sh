#!/bin/bash

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -t|--trim)
    TRIM="$2"
    shift # past argument
    ;;
    -l|--length)
    LENGTH="$2"
    shift # past argument
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
    # unknown option
    ;;
esac
shift # past argument or value
done


module load trimmomatic
module load bowtie/2-2.2.9
module load fastqc

java -Djava.io.tmpdir=. -jar $TRIMMOJAR SE \
    /data/ChioriniCompCor/Drews_ATAC_vs_DNase_H2H/data/ATAC_GM_Buenrostro_50k/SRR891271.fastq \
    $TRIM.$LENGTH.output.fastq HEADCROP:$TRIM CROP:$LENGTH
rm $TRIM.$LENGTH.output.fastq
