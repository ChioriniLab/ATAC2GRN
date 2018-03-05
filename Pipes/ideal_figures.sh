#!/bin/sh
# Thomas Pranzatelli
# 6/23/16
# NIDCR
# A master bash script for a footprinting pipeline.
#
# This part of the script accepts the four input variables:
# -d, -e, -i, -s

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -d|--data)
    DATA="$2"
    shift # past argument
    ;;
    -e|--experiment)
    EXPERIMENT="$2"
    shift # past argument
    ;;
    -g|--genome)
    GENOME="$2"
    shift # past argument
    ;;
    -fp|--footprinting)
    FOOTPRINTING="$2"
    shift # past argument
    ;;
    -bias|--bias-correction)
    BIAS="$2"
    shift # past argument
    ;;
    -fdr|--wellington-fdr-cutoff)
    WELLFDR="$2"
    shift # past argument
    ;;
    -fdrlimit|--wellington-fdr-limit)
    WELLFDRLIMIT="$2"
    shift # past argument
    ;;
    -md|--homer-mindist)
    HOMERMINDIST="$2"
    shift # past argument
    ;;
    -h|--homer-size)
    HOMERSIZE="$2"
    shift # past argument
    ;;
    -p|--picard)
    PICARD="$2"
    shift # past argument
    ;;
    -a|--alignment)
    ALIGNMENT="$2"
    shift # past argument
    ;;
    -atac|--sequencing)
    ATAC="$2"
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

# Here we load the different modules.
module load bowtie/2-2.2.9
module load samtools
module load homer/4.8.2
module load bedtools/2.25.0
module load pyDNase/0.2.4
module load R
module load deeptools/2.2.4
module load picard
module load ucsc
module load rgt
module load bedops

cp -R $DATA/ChIP-validation $EXPERIMENT/
for chip in $EXPERIMENT/ChIP-validation/*
    do bedtools intersect -a $chip/segBarozzi.bed -b $EXPERIMENT/footprints_R1.bed -wa -wb > $chip/positives.bed
    bedtools intersect -a $chip/segBarozzi.bed -b $EXPERIMENT/footprints_R1.bed -v > $chip/negatives.bed
    bedtools intersect -a $chip/positives.bed -b $chip/chipseqresults.bed -wa > $chip/true_positives.bed
    bedtools intersect -a $chip/positives.bed -b $chip/chipseqresults.bed -v > $chip/false_positives.bed
    bedtools intersect -a $chip/negatives.bed -b $chip/chipseqresults.bed -wa > $chip/false_negatives.bed
    bedtools intersect -a $chip/negatives.bed -b $chip/chipseqresults.bed -v > $chip/true_negatives.bed
done

> $EXPERIMENT/overlap.txt
for num in {1..100}
do bedops -e $num% $EXPERIMENT/footprints_R1.bed $EXPERIMENT/footprints_R2.bed > $EXPERIMENT/overlaP.bed
overlap=($(wc -l $EXPERIMENT/overlaP.bed))
original=($(wc -l $EXPERIMENT/footprints_R1.bed))
echo print $num, $overlap/$original. | python >> $EXPERIMENT/overlap.txt
done

rm overlaP.bed
