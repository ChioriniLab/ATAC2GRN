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
    -e|--experiment)
    EXP="$2"
    shift # past argument
    ;;
    -g|--genome)
    GEN="$2"
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

EXPERIMENT=$EXP
DATA=data
GENOME=genomes/$GEN

# Here we load the different modules.
module load bowtie
module load samtools
module load homer
module load bedtools
module load rgt
module load bedops
module load trimmomatic
module load picard

# .fastq files in each paired end folder are concatenated.
cat $DATA/*.1.fastq.gz > $EXPERIMENT/P1.fastq.gz
cat $DATA/*.2.fastq.gz > $EXPERIMENT/P2.fastq.gz
java -Djava.io.tmpdir=. -jar $TRIMMOJAR PE \
    $EXPERIMENT/P1.fastq.gz $EXPERIMENT/P2.fastq.gz \
    $EXPERIMENT/P1.fq.gz /lscratch/$SLURM_JOBID/output_forward_unpaired.fq.gz \
    $EXPERIMENT/P2.fq.gz /lscratch/$SLURM_JOBID/output_reverse_unpaired.fq.gz \
    HEADCROP:20
echo "Paired-end reads were concatenated."
# These files are passed to bowtie2 to be aligned to the given genome index.
bowtie2 -x $GENOME/Bowtie2Index/genome --very-sensitive -p 4 -t -q -1 $EXPERIMENT/P1.fq.gz -2 $EXPERIMENT/P2.fq.gz -S $EXPERIMENT/bowtie2.sam
echo "Paired-end reads were aligned end-to-end."
rm $EXPERIMENT/*fastq.gz
rm $EXPERIMENT/*fq.gz

# This alignment is converted to binary, and then sorted and indexed.
java -Xmx10g -XX:ParallelGCThreads=3 -jar $PICARDJARPATH/picard.jar SortSam INPUT=$EXPERIMENT/bowtie2.sam OUTPUT=$EXPERIMENT/sorted.bam SORT_ORDER=coordinate
> $EXPERIMENT/picard.bam
java -Xmx10g -XX:ParallelGCThreads=3 -jar $PICARDJARPATH/picard.jar MarkDuplicates INPUT=$EXPERIMENT/sorted.bam OUTPUT=$EXPERIMENT/picard.bam METRICS_FILE=/lscratch/$SLURM_JOBID/marked_metrics.txt REMOVE_DUPLICATES=true
java -Xmx10g -XX:ParallelGCThreads=3 -jar $PICARDJARPATH/picard.jar BuildBamIndex INPUT=$EXPERIMENT/picard.bam
rm $EXPERIMENT/bowtie2.sam
rm $EXPERIMENT/sorted.bam

# Peaks are estimated from the .bam files by HOMER.
makeTagDirectory $EXPERIMENT/tagDirectory/ $EXPERIMENT/picard.bam
findPeaks $EXPERIMENT/tagDirectory/ -region -size 500 -minDist 50 -o auto -tbp 0
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

# The .bed and .bam files are used by HINT to estimate footprints.
cp /data/ChioriniCompCor/experiments/2017-09-01_Metamachine_Footprinting_Pipelines/patroloffice/experiment_matrix.bed $EXPERIMENT/
rm -R $EXPERIMENT/footprints
mkdir /lscratch/$SLURM_JOBID/footprints
cd $EXPERIMENT
rgt-hint --output-location /lscratch/$SLURM_JOBID/footprints/ --organism $GEN experiment_matrix.bed || fail "HINT-BC failed to run to completion."
cd -
cp -R /lscratch/$SLURM_JOBID/footprints $EXPERIMENT || fail "Transmitting the footprints from lscratch failed to run to completion."
echo "Footprint intervals produced."
rm -r $EXPERIMENT/tagDirectory
rm $EXPERIMENT/picard.bam
rm $EXPERIMENT/experiment_matrix.bed
rm $EXPERIMENT/picard.bai