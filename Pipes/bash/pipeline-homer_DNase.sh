#!/bin/sh
# Thomas Pranzatelli
# 6/23/16
# NIDCR
# A master bash script for a footprinting pipeline.
#
# This part of the script accepts the four input variables:
# -d, -e, -i, -s

function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -e|--encsr)
    ENCSR="$2"
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

READS=1
EXPERIMENT=/data/pranzatellitj/check_HOMER_stitching/homer_DNase/$ENCSR
DATA=/data/pranzatellitj/check_HOMER_stitching/titrated_data_DNase/211792753.fastq
GENOME=/data/ChioriniCompCor/metamachine/genomes/hg19

# Here we load the different modules.
module load bowtie/2-2.2.9
module load samtools
module load homer/4.8.2
module load bedtools/2.25.0
module load pyDNase/0.2.4
module load bedops
module load rgt

if [[ $READS == 1 ]]
    then
    # .fastq files are concatenated.
    cat $DATA > $EXPERIMENT/cat.fastq
    echo "Single-end reads were concatenated."
    wc -l $EXPERIMENT/cat.fastq
    # These files are passed to bowtie2 to be aligned to the given genome index.
    bowtie2 -x $GENOME/Bowtie2Index/genome -p 8 -t -q -U $EXPERIMENT/cat.fastq -S $EXPERIMENT/bowtie2.sam
    echo "Single-end reads were aligned end-to-end."
    rm $EXPERIMENT/cat.fastq
elif [[ $READS == 2 ]]
    then
    # .fastq files in each paired end folder are concatenated.
    cat $DATA/P1/*.fastq > $EXPERIMENT/P1.fastq
    cat $DATA/P2/*.fastq > $EXPERIMENT/P2.fastq
    echo "Paired-end reads were concatenated."
    wc -l $EXPERIMENT/P1.fastq
    # These files are passed to bowtie2 to be aligned to the given genome index.
    bowtie2 -x $GENOME/Bowtie2Index/genome -p 8 -t -q -1 $EXPERIMENT/P1.fastq -2 $EXPERIMENT/P2.fastq -S $EXPERIMENT/bowtie2.sam
    echo "Paired-end reads were aligned end-to-end."
    rm $EXPERIMENT/P1.fastq
    rm $EXPERIMENT/P2.fastq
else
    echo "Use a read count argument of -r 1 for single-end reads and a read count argument of -r 2 for paired-end reads."
    exit
fi

# This alignment is converted to binary, and then sorted and indexed.
samtools view -u -b $EXPERIMENT/bowtie2.sam > $EXPERIMENT/view.bam
samtools sort -@ 7 $EXPERIMENT/view.bam -o $EXPERIMENT/picard.bam
samtools index $EXPERIMENT/picard.bam
rm $EXPERIMENT/bowtie2.sam
rm $EXPERIMENT/view.bam

# Peaks are estimated from the .bam files by HOMER.
makeTagDirectory $EXPERIMENT/tagDirectory/ $EXPERIMENT/picard.bam
findPeaks $EXPERIMENT/tagDirectory/ -region -size $ENCSR -minDist 50 -o auto -tbp 0
echo "Open chromatin peaks found."
# These peaks are converted to .bed files and sorted and merged.
pos2bed.pl $EXPERIMENT/tagDirectory/peaks.txt > $EXPERIMENT/peaks.bed
echo "Peak .bed file produced."
bedtools sort -i $EXPERIMENT/peaks.bed > $EXPERIMENT/peaks.sorted.bed
echo "Peak .bed file sorted."
grep -v "#" $EXPERIMENT/tagDirectory/peaks.txt > $EXPERIMENT/tagDirectory/peaks2.txt
mv $EXPERIMENT/tagDirectory/peaks2.txt $EXPERIMENT/tagDirectory/peaks.txt
cut -f2-4,8 $EXPERIMENT/tagDirectory/peaks.txt > $EXPERIMENT/peaks.bedgraph
echo "Peak .bedgraph file produced."
rm $EXPERIMENT/peaks.bed
bedtools merge -i $EXPERIMENT/peaks.sorted.bed > $EXPERIMENT/peaks.merged.bed
rm $EXPERIMENT/peaks.sorted.bed
echo "Peaks merged and sorted for interval estimation."

cp /data/pranzatellitj/check_HOMER_stitching/experiment_matrix.bed $EXPERIMENT/
mkdir $EXPERIMENT/hintbc_output
cd $EXPERIMENT
rgt-hint --default-bias-correction --output-location $EXPERIMENT/hintbc_output/ $EXPERIMENT/experiment_matrix.bed
sort -k1,1 -k2,2n -k3,3n $EXPERIMENT/hintbc_output/*.bed > $EXPERIMENT/footprints.bed
rm $EXPERIMENT/experiment_matrix.bed
rm -r $EXPERIMENT/hintbc_output

cp -R /data/ChioriniCompCor/Drews_ATAC_vs_DNase_H2H/ChIP-validation $EXPERIMENT/
for chip in $EXPERIMENT/ChIP-validation/*
    do bedtools intersect -a $chip/segBarozzi.bed -b $EXPERIMENT/footprints.bed -wa -wb > $chip/positives.bed
    bedtools intersect -a $chip/segBarozzi.bed -b $EXPERIMENT/footprints.bed -v > $chip/negatives.bed
    bedtools intersect -a $chip/positives.bed -b $chip/chipseqresults.bed -wa > $chip/true_positives.bed
    bedtools intersect -a $chip/positives.bed -b $chip/chipseqresults.bed -v > $chip/false_positives.bed
    bedtools intersect -a $chip/negatives.bed -b $chip/chipseqresults.bed -wa > $chip/false_negatives.bed
    bedtools intersect -a $chip/negatives.bed -b $chip/chipseqresults.bed -v > $chip/true_negatives.bed
done

module load python
python /data/pranzatellitj/tools/ROC-validation.py -e $EXPERIMENT -f -5 > $EXPERIMENT/AUC.txt
rm -R $EXPERIMENT/ChIP-validation

