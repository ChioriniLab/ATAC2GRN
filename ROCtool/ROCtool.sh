#!/bin/sh
# Thomas Pranzatelli
# 6/23/16
# NIDCR
# Pass the footprint bed file as -fp.

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -fp|--footprints)
    FOOTPRINT="$2"
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

for chip in ChIP-validation/*
    do bedtools intersect -a $chip/segBarozzi.bed -b $FOOTPRINT -wa -wb > $chip/positives.bed
    bedtools intersect -a $chip/segBarozzi.bed -b $FOOTPRINT -v > $chip/negatives.bed
    bedtools intersect -a $chip/positives.bed -b $chip/chipseqresults.bed -wa > $chip/true_positives.bed
    bedtools intersect -a $chip/positives.bed -b $chip/chipseqresults.bed -v > $chip/false_positives.bed
    bedtools intersect -a $chip/negatives.bed -b $chip/chipseqresults.bed -wa > $chip/false_negatives.bed
    bedtools intersect -a $chip/negatives.bed -b $chip/chipseqresults.bed -v > $chip/true_negatives.bed
done

python ROCtool.py