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
module load bedops
module load pyDNase/0.2.4
module load R
module load deeptools/2.2.4
module load picard
module load ucsc
module load rgt
module load trimmomatic

# .fastq files are concatenated.
if [[ $ATAC == 0 ]]
    then cat $DATA/data/DNase_GM_Stam/*.fastq > /lscratch/$SLURM_JOBID/cat.fastq
fi
if [[ $ATAC == 1 ]]
    then cat $DATA/data/ATAC_GM_Buenrostro_50k_1/*.fastq > /lscratch/$SLURM_JOBID/cat.fastq
    java -Djava.io.tmpdir=. -jar $TRIMMOJAR SE \
    /lscratch/$SLURM_JOBID/cat.fastq /lscratch/$SLURM_JOBID/trimmed.fastq CROP:20
    mv /lscratch/$SLURM_JOBID/trimmed.fastq /lscratch/$SLURM_JOBID/cat.fastq
fi
echo "Single-end reads were concatenated."
# These files are passed to bowtie2 to be aligned to the given genome index.
if [[ $ALIGNMENT == 1 ]]
    then bowtie2 -x $GENOME/Bowtie2Index/genome -p 16 -t -q -U /lscratch/$SLURM_JOBID/cat.fastq -S /lscratch/$SLURM_JOBID/bowtie2.sam
    echo "Single-end reads were aligned end-to-end."
elif [[ $ALIGNMENT == 2 ]]
    then bowtie2 -x $GENOME/Bowtie2Index/genome -p 16 -t -q --sensitive-local -U /lscratch/$SLURM_JOBID/cat.fastq -S /lscratch/$SLURM_JOBID/bowtie2.sam
    echo "Single-end reads were aligned locally."
elif [[ $ALIGNMENT == 3 ]]
    then bowtie2 -x $GENOME/Bowtie2Index/genome -p 16 -t -q --very-sensitive -U /lscratch/$SLURM_JOBID/cat.fastq -S /lscratch/$SLURM_JOBID/bowtie2.sam
    echo "Single-end reads were very sensitively aligned end-to-end."
elif [[ $ALIGNMENT == 4 ]]
    then bowtie2 -x $GENOME/Bowtie2Index/genome -p 16 -t -q --very-sensitive-local -U /lscratch/$SLURM_JOBID/cat.fastq -S /lscratch/$SLURM_JOBID/bowtie2.sam
    echo "Single-end reads were very sensitively aligned locally."
else
    echo "Use -a 1 for end-to-end behavior, -a 2 for local behavior, -a 3 for sensitive end-to-end behavior, and -a 4 for sensitive local behavior."
    exit
fi
rm /lscratch/$SLURM_JOBID/cat.fastq


if [[ $PICARD == 1 ]]
    then # This alignment is converted to binary, and then sorted and indexed.
    java -Xmx10g -XX:ParallelGCThreads=15 -jar $PICARDJARPATH/picard.jar SortSam INPUT=/lscratch/$SLURM_JOBID/bowtie2.sam OUTPUT=/lscratch/$SLURM_JOBID/sorted.bam SORT_ORDER=coordinate
    > /lscratch/$SLURM_JOBID/picard.bam
    java -Xmx10g -XX:ParallelGCThreads=15 -jar $PICARDJARPATH/picard.jar MarkDuplicates INPUT=/lscratch/$SLURM_JOBID/sorted.bam OUTPUT=/lscratch/$SLURM_JOBID/picard_R1.bam METRICS_FILE=/lscratch/$SLURM_JOBID/marked_metrics.txt REMOVE_DUPLICATES=true
    rm /lscratch/$SLURM_JOBID/marked_metrics.txt
    if [[ FOOTPRINTING == 2 && ATAC == 1 ]]
        then bedtools bamtobed -i /lscratch/$SLURM_JOBID/picard_R1.bam -bedpe | awk -v OFS="\t" '{if($9=="+"){print $1,$2+4,$6+4}else if($9=="-"){print $1,$2-5,$6-5}}' > /lscratch/$SLURM_JOBID/picard_R1.bed
        bedtools bedpetobam -i /lscratch/$SLURM_JOBID/picard_R1.bed -g human.hg19.genome > /lscratch/$SLURM_JOBID/picard_R1.bam
    fi
    java -Xmx10g -XX:ParallelGCThreads=15 -jar $PICARDJARPATH/picard.jar BuildBamIndex INPUT=/lscratch/$SLURM_JOBID/picard_R1.bam
    rm /lscratch/$SLURM_JOBID/bowtie2.sam
    rm /lscratch/$SLURM_JOBID/sorted.bam
else
    samtools view -b /lscratch/$SLURM_JOBID/bowtie2.sam > /lscratch/$SLURM_JOBID/view.bam
    samtools sort /lscratch/$SLURM_JOBID/view.bam -o /lscratch/$SLURM_JOBID/picard_R1.bam
    if [[ FOOTPRINTING == 2 && ATAC == 1 ]]
        then bedtools bamtobed -i /lscratch/$SLURM_JOBID/picard_R1.bam -bedpe | awk -v OFS="\t" '{if($9=="+"){print $1,$2+4,$6+4}else if($9=="-"){print $1,$2-5,$6-5}}' > /lscratch/$SLURM_JOBID/picard_R1.bed
        bedtools bedpetobam -i /lscratch/$SLURM_JOBID/picard_R1.bed -g human.hg19.genome > /lscratch/$SLURM_JOBID/picard_R1.bam
    fi
    samtools index /lscratch/$SLURM_JOBID/picard_R1.bam
    rm /lscratch/$SLURM_JOBID/bowtie2.sam
    rm /lscratch/$SLURM_JOBID/view.bam
fi


# .fastq files are concatenated.
if [[ $ATAC == 0 ]]
    then cat $DATA/data/DNase_GM_Crawford/*.fastq > /lscratch/$SLURM_JOBID/cat.fastq
fi
if [[ $ATAC == 1 ]]
    then cat $DATA/data/ATAC_GM_Buenrostro_50k_2/*.fastq > /lscratch/$SLURM_JOBID/cat.fastq
    java -Djava.io.tmpdir=. -jar $TRIMMOJAR SE \
    /lscratch/$SLURM_JOBID/cat.fastq /lscratch/$SLURM_JOBID/trimmed.fastq CROP:20
    mv /lscratch/$SLURM_JOBID/trimmed.fastq /lscratch/$SLURM_JOBID/cat.fastq
fi
echo "Single-end reads were concatenated."
# These files are passed to bowtie2 to be aligned to the given genome index.
if [[ $ALIGNMENT == 1 ]]
    then bowtie2 -x $GENOME/Bowtie2Index/genome -p 16 -t -q -U /lscratch/$SLURM_JOBID/cat.fastq -S /lscratch/$SLURM_JOBID/bowtie2.sam
    echo "Single-end reads were aligned end-to-end."
elif [[ $ALIGNMENT == 2 ]]
    then bowtie2 -x $GENOME/Bowtie2Index/genome -p 16 -t -q --sensitive-local -U /lscratch/$SLURM_JOBID/cat.fastq -S /lscratch/$SLURM_JOBID/bowtie2.sam
    echo "Single-end reads were aligned locally."
elif [[ $ALIGNMENT == 3 ]]
    then bowtie2 -x $GENOME/Bowtie2Index/genome -p 16 -t -q --very-sensitive -U /lscratch/$SLURM_JOBID/cat.fastq -S /lscratch/$SLURM_JOBID/bowtie2.sam
    echo "Single-end reads were very sensitively aligned end-to-end."
elif [[ $ALIGNMENT == 4 ]]
    then bowtie2 -x $GENOME/Bowtie2Index/genome -p 16 -t -q --very-sensitive-local -U /lscratch/$SLURM_JOBID/cat.fastq -S /lscratch/$SLURM_JOBID/bowtie2.sam
    echo "Single-end reads were very sensitively aligned locally."
else
    echo "Use -a 1 for end-to-end behavior, -a 2 for local behavior, -a 3 for sensitive end-to-end behavior, and -a 4 for sensitive local behavior."
    exit
fi
rm /lscratch/$SLURM_JOBID/cat.fastq


if [[ $PICARD == 1 ]]
    then # This alignment is converted to binary, and then sorted and indexed.
    java -Xmx10g -XX:ParallelGCThreads=15 -jar $PICARDJARPATH/picard.jar SortSam INPUT=/lscratch/$SLURM_JOBID/bowtie2.sam OUTPUT=/lscratch/$SLURM_JOBID/sorted.bam SORT_ORDER=coordinate
    > /lscratch/$SLURM_JOBID/picard.bam
    java -Xmx10g -XX:ParallelGCThreads=15 -jar $PICARDJARPATH/picard.jar MarkDuplicates INPUT=/lscratch/$SLURM_JOBID/sorted.bam OUTPUT=/lscratch/$SLURM_JOBID/picard_R2.bam METRICS_FILE=/lscratch/$SLURM_JOBID/marked_metrics.txt REMOVE_DUPLICATES=true
    rm /lscratch/$SLURM_JOBID/marked_metrics.txt
    if [[ FOOTPRINTING == 2 && ATAC == 1 ]]
        then bedtools bamtobed -i /lscratch/$SLURM_JOBID/picard_R2.bam -bedpe | awk -v OFS="\t" '{if($9=="+"){print $1,$2+4,$6+4}else if($9=="-"){print $1,$2-5,$6-5}}' > /lscratch/$SLURM_JOBID/picard_R2.bed
        bedtools bedpetobam -i /lscratch/$SLURM_JOBID/picard_R2.bed -g human.hg19.genome > /lscratch/$SLURM_JOBID/picard_R2.bam
    fi
    java -Xmx10g -XX:ParallelGCThreads=15 -jar $PICARDJARPATH/picard.jar BuildBamIndex INPUT=/lscratch/$SLURM_JOBID/picard_R2.bam
    rm /lscratch/$SLURM_JOBID/bowtie2.sam
    rm /lscratch/$SLURM_JOBID/sorted.bam
else
    samtools view -b /lscratch/$SLURM_JOBID/bowtie2.sam > /lscratch/$SLURM_JOBID/view.bam
    samtools sort /lscratch/$SLURM_JOBID/view.bam -o /lscratch/$SLURM_JOBID/picard_R2.bam
    if [[ FOOTPRINTING == 2 && ATAC == 1 ]]
        then bedtools bamtobed -i /lscratch/$SLURM_JOBID/picard_R2.bam -bedpe | awk -v OFS="\t" '{if($9=="+"){print $1,$2+4,$6+4}else if($9=="-"){print $1,$2-5,$6-5}}' > /lscratch/$SLURM_JOBID/picard_R2.bed
        bedtools bedpetobam -i /lscratch/$SLURM_JOBID/picard_R2.bed -g human.hg19.genome > /lscratch/$SLURM_JOBID/picard_R2.bam
    fi
    samtools index /lscratch/$SLURM_JOBID/picard_R2.bam
    rm /lscratch/$SLURM_JOBID/bowtie2.sam
    rm /lscratch/$SLURM_JOBID/view.bam
fi

# Here, the .bam files generated are compared to the data copy .bam files.
multiBamSummary bins --bamfiles /lscratch/$SLURM_JOBID/picard_R1.bam /lscratch/$SLURM_JOBID/picard_R2.bam -out /lscratch/$SLURM_JOBID/bamresults.npz -p 16
plotCorrelation -in /lscratch/$SLURM_JOBID/bamresults.npz -o /lscratch/$SLURM_JOBID/bamp.pdf -c pearson -p scatterplot --removeOutliers --outFileCorMatrix $EXPERIMENT/bam_pearson.tab
rm /lscratch/$SLURM_JOBID/bamp.pdf
plotCorrelation -in /lscratch/$SLURM_JOBID/bamresults.npz -o /lscratch/$SLURM_JOBID/bams.pdf -c spearman -p scatterplot --removeOutliers --outFileCorMatrix $EXPERIMENT/bam_spearman.tab
rm /lscratch/$SLURM_JOBID/bams.pdf
rm /lscratch/$SLURM_JOBID/bamresults.npz

# Peaks are estimated from the .bam files by HOMER.
makeTagDirectory /lscratch/$SLURM_JOBID/tagDirectory_R1/ /lscratch/$SLURM_JOBID/picard_R1.bam
findPeaks /lscratch/$SLURM_JOBID/tagDirectory_R1/ -region -size $HOMERSIZE -minDist $HOMERMINDIST -o auto -tbp 0
echo "Open chromatin peaks found."
# These peaks are converted to .bed files and sorted and merged.
pos2bed.pl /lscratch/$SLURM_JOBID/tagDirectory_R1/peaks.txt > /lscratch/$SLURM_JOBID/peaks.bed
echo "Peak .bed file produced."
bedtools sort -i /lscratch/$SLURM_JOBID/peaks.bed > /lscratch/$SLURM_JOBID/peaks.sorted.bed
echo "Peak .bed file sorted."
rm /lscratch/$SLURM_JOBID/peaks.bed
bedtools merge -i /lscratch/$SLURM_JOBID/peaks.sorted.bed > /lscratch/$SLURM_JOBID/R1peaks.merged.bed
rm /lscratch/$SLURM_JOBID/peaks.sorted.bed
echo "Peaks merged and sorted for interval estimation."

# Peaks are estimated from the .bam files by HOMER.
makeTagDirectory /lscratch/$SLURM_JOBID/tagDirectory_R2/ /lscratch/$SLURM_JOBID/picard_R2.bam
findPeaks /lscratch/$SLURM_JOBID/tagDirectory_R2/ -region -size $HOMERSIZE -minDist $HOMERMINDIST -o auto -tbp 0
echo "Open chromatin peaks found."
# These peaks are converted to .bed files and sorted and merged.
pos2bed.pl /lscratch/$SLURM_JOBID/tagDirectory_R2/peaks.txt > /lscratch/$SLURM_JOBID/peaks.bed
echo "Peak .bed file produced."
bedtools sort -i /lscratch/$SLURM_JOBID/peaks.bed > /lscratch/$SLURM_JOBID/peaks.sorted.bed
echo "Peak .bed file sorted."
rm /lscratch/$SLURM_JOBID/peaks.bed
bedtools merge -i /lscratch/$SLURM_JOBID/peaks.sorted.bed > /lscratch/$SLURM_JOBID/R2peaks.merged.bed
rm /lscratch/$SLURM_JOBID/peaks.sorted.bed
echo "Peaks merged and sorted for interval estimation."

# For comparison with the data copy, the output .bed files are converted to .bw.
bedtools genomecov -i /lscratch/$SLURM_JOBID/R1peaks.merged.bed -g $GENOME/*.sizes -bga > /lscratch/$SLURM_JOBID/peaks.bg
LC_COLLATE=C sort -k1,1 -k2,2n /lscratch/$SLURM_JOBID/peaks.bg > /lscratch/$SLURM_JOBID/peaks.collate.bg
rm /lscratch/$SLURM_JOBID/peaks.bg
bedGraphToBigWig /lscratch/$SLURM_JOBID/peaks.collate.bg $GENOME/*.sizes /lscratch/$SLURM_JOBID/peaks_R1.bw
rm /lscratch/$SLURM_JOBID/peaks.collate.bg
# For comparison with the data copy, the output .bed files are converted to .bw.
bedtools genomecov -i /lscratch/$SLURM_JOBID/R2peaks.merged.bed -g $GENOME/*.sizes -bga > /lscratch/$SLURM_JOBID/peaks.bg
LC_COLLATE=C sort -k1,1 -k2,2n /lscratch/$SLURM_JOBID/peaks.bg > /lscratch/$SLURM_JOBID/peaks.collate.bg
rm /lscratch/$SLURM_JOBID/peaks.bg
bedGraphToBigWig /lscratch/$SLURM_JOBID/peaks.collate.bg $GENOME/*.sizes /lscratch/$SLURM_JOBID/peaks_R2.bw
rm /lscratch/$SLURM_JOBID/peaks.collate.bg

# As before, these files are compared to the data .bw files.
multiBigwigSummary bins -b /lscratch/$SLURM_JOBID/peaks_R1.bw /lscratch/$SLURM_JOBID/peaks_R2.bw -out /lscratch/$SLURM_JOBID/bigresults.npz -p 16
plotCorrelation --corData /lscratch/$SLURM_JOBID/bigresults.npz --plotFile /lscratch/$SLURM_JOBID/bigp.pdf --corMethod pearson --whatToPlot scatterplot --removeOutliers --outFileCorMatrix $EXPERIMENT/bed_pearson.tab
rm /lscratch/$SLURM_JOBID/bigp.pdf
plotCorrelation --corData /lscratch/$SLURM_JOBID/bigresults.npz --plotFile /lscratch/$SLURM_JOBID/bigs.pdf --corMethod spearman --whatToPlot scatterplot --removeOutliers --outFileCorMatrix $EXPERIMENT/bed_spearman.tab
rm /lscratch/$SLURM_JOBID/bigs.pdf
rm /lscratch/$SLURM_JOBID/bigresults.npz
rm /lscratch/$SLURM_JOBID/peaks_R1.bw
rm /lscratch/$SLURM_JOBID/peaks_R2.bw

if [[ $FOOTPRINTING == 2 ]]
    then
    sh /data/pranzatellitj/test_piq/PIQ/piq.sh -e $EXPERIMENT
fi

if [[ $FOOTPRINTING == 1 ]]
    then
    cp /data/pranzatellitj/experiments/2016-10-17-JUBAL_Cotillion_DNase_Grid_Search/experiment_matrix_R1.bed /lscratch/$SLURM_JOBID/
    cp /data/pranzatellitj/experiments/2016-10-17-JUBAL_Cotillion_DNase_Grid_Search/experiment_matrix_R2.bed /lscratch/$SLURM_JOBID/
    mkdir /lscratch/$SLURM_JOBID/R1_hintbc_output
    mkdir /lscratch/$SLURM_JOBID/R2_hintbc_output
    cd /lscratch/$SLURM_JOBID
    if [[ $BIAS == 0 ]]
        then rgt-hint --output-location /lscratch/$SLURM_JOBID/R1_hintbc_output/ /lscratch/$SLURM_JOBID/experiment_matrix_R1.bed
        rgt-hint --output-location /lscratch/$SLURM_JOBID/R2_hintbc_output/ /lscratch/$SLURM_JOBID/experiment_matrix_R2.bed
    fi
    if [[ $BIAS == 1 ]]
        then rgt-hint --default-bias-correction --output-location /lscratch/$SLURM_JOBID/R1_hintbc_output/ /lscratch/$SLURM_JOBID/experiment_matrix_R1.bed
        rgt-hint --default-bias-correction --output-location /lscratch/$SLURM_JOBID/R2_hintbc_output/ /lscratch/$SLURM_JOBID/experiment_matrix_R2.bed
    fi
    if [[ $BIAS == 2 ]]
        then rgt-hint --estimate-bias-correction --output-location /lscratch/$SLURM_JOBID/R1_hintbc_output/ /lscratch/$SLURM_JOBID/experiment_matrix_R1.bed
        rgt-hint --estimate-bias-correction --output-location /lscratch/$SLURM_JOBID/R2_hintbc_output/ /lscratch/$SLURM_JOBID/experiment_matrix_R2.bed
    fi
    sort -k1,1 -k2,2n -k3,3n /lscratch/$SLURM_JOBID/R1_hintbc_output/*.bed > /lscratch/$SLURM_JOBID/footprints_R1.bed
    sort -k1,1 -k2,2n -k3,3n /lscratch/$SLURM_JOBID/R2_hintbc_output/*.bed > /lscratch/$SLURM_JOBID/footprints_R2.bed
    rm -r /lscratch/$SLURM_JOBID/tagDirectory_R1
    rm /lscratch/$SLURM_JOBID/R1peaks.merged.bed
    rm /lscratch/$SLURM_JOBID/picard_R1.bam
    rm /lscratch/$SLURM_JOBID/*.bai
    rm -r /lscratch/$SLURM_JOBID/tagDirectory_R2
    rm /lscratch/$SLURM_JOBID/R2peaks.merged.bed
    rm /lscratch/$SLURM_JOBID/picard_R2.bam
    rm $EXPERIMENT/experiment_matrix_R1.bed
    rm $EXPERIMENT/experiment_matrix_R2.bed

    cp -R $DATA/ChIP-validation /lscratch/$SLURM_JOBID/ChIP-validation
    for chip in /lscratch/$SLURM_JOBID/ChIP-validation/*
        do bedtools intersect -a $chip/segBarozzi.bed -b /lscratch/$SLURM_JOBID/footprints_R1.bed -wa -wb > $chip/positives.bed
        bedtools intersect -a $chip/segBarozzi.bed -b /lscratch/$SLURM_JOBID/footprints_R1.bed -v > $chip/negatives.bed
        bedtools intersect -a $chip/positives.bed -b $chip/chipseqresults.bed -wa > $chip/true_positives.bed
        bedtools intersect -a $chip/positives.bed -b $chip/chipseqresults.bed -v > $chip/false_positives.bed
        bedtools intersect -a $chip/negatives.bed -b $chip/chipseqresults.bed -wa > $chip/false_negatives.bed
        bedtools intersect -a $chip/negatives.bed -b $chip/chipseqresults.bed -v > $chip/true_negatives.bed
    done

    python /data/pranzatellitj/tools/ROC-validation.py -e /lscratch/$SLURM_JOBID -f -10 > $EXPERIMENT/AUC.txt
    rm -R /lscratch/$SLURM_JOBID/ChIP-validation

fi

if [[ $FOOTPRINTING == 0 ]]
    then
    # The .bed and .bam files are used by Wellington to estimate footprints.
    rm -R /lscratch/$SLURM_JOBID/footprints_R1
    mkdir /lscratch/$SLURM_JOBID/footprints_R1
    if [[ $ATAC == 1 ]]
        then wellington_footprints.py -A -p 15 -fdr $WELLFDR -fdrlimit $WELLFDRLIMIT /lscratch/$SLURM_JOBID/R1peaks.merged.bed /lscratch/$SLURM_JOBID/picard_R1.bam /lscratch/$SLURM_JOBID/footprints_R1/
    fi
    if [[ $ATAC == 0 ]]
        then wellington_footprints.py -p 15 -fdr $WELLFDR -fdrlimit $WELLFDRLIMIT /lscratch/$SLURM_JOBID/R1peaks.merged.bed /lscratch/$SLURM_JOBID/picard_R1.bam /lscratch/$SLURM_JOBID/footprints_R1/
    fi
    sort -k1,1 -k2,2n -k3,3n /lscratch/$SLURM_JOBID/footprints_R1/*.bed > /lscratch/$SLURM_JOBID/footprints_R1.bed
    echo "Footprint intervals produced."
    rm -r /lscratch/$SLURM_JOBID/tagDirectory_R1
    rm /lscratch/$SLURM_JOBID/R1peaks.merged.bed
    rm /lscratch/$SLURM_JOBID/picard_R1.bam
    rm /lscratch/$SLURM_JOBID/picard_R1.bai
    # The .bed and .bam files are used by Wellington to estimate footprints.
    rm -R /lscratch/$SLURM_JOBID/footprints_R2
    mkdir /lscratch/$SLURM_JOBID/footprints_R2
    if [[ $ATAC == 1 ]]
        then wellington_footprints.py -A -p 15 -fdr $WELLFDR -fdrlimit $WELLFDRLIMIT /lscratch/$SLURM_JOBID/R2peaks.merged.bed /lscratch/$SLURM_JOBID/picard_R2.bam /lscratch/$SLURM_JOBID/footprints_R2/
    fi
    if [[ $ATAC == 0 ]]
        then wellington_footprints.py -p 15 -fdr $WELLFDR -fdrlimit $WELLFDRLIMIT /lscratch/$SLURM_JOBID/R2peaks.merged.bed /lscratch/$SLURM_JOBID/picard_R2.bam /lscratch/$SLURM_JOBID/footprints_R2/
    fi
    sort -k1,1 -k2,2n -k3,3n /lscratch/$SLURM_JOBID/footprints_R2/*.bed > /lscratch/$SLURM_JOBID/footprints_R2.bed
    echo "Footprint intervals produced."
    rm -r /lscratch/$SLURM_JOBID/tagDirectory_R2
    rm /lscratch/$SLURM_JOBID/R2peaks.merged.bed
    rm /lscratch/$SLURM_JOBID/picard_R2.bam
    rm /lscratch/$SLURM_JOBID/picard_R2.bai
    rm -R /lscratch/$SLURM_JOBID/footprints_R1
    rm -R /lscratch/$SLURM_JOBID/footprints_R2

    cp -R $DATA/ChIP-validation /lscratch/$SLURM_JOBID/ChIP-validation
    for chip in /lscratch/$SLURM_JOBID/ChIP-validation/*
        do bedtools intersect -a $chip/segBarozzi.bed -b /lscratch/$SLURM_JOBID/footprints_R1.bed -wa -wb > $chip/positives.bed
        bedtools intersect -a $chip/segBarozzi.bed -b /lscratch/$SLURM_JOBID/footprints_R1.bed -v > $chip/negatives.bed
        bedtools intersect -a $chip/positives.bed -b $chip/chipseqresults.bed -wa > $chip/true_positives.bed
        bedtools intersect -a $chip/positives.bed -b $chip/chipseqresults.bed -v > $chip/false_positives.bed
        bedtools intersect -a $chip/negatives.bed -b $chip/chipseqresults.bed -wa > $chip/false_negatives.bed
        bedtools intersect -a $chip/negatives.bed -b $chip/chipseqresults.bed -v > $chip/true_negatives.bed
    done

    python /data/pranzatellitj/tools/ROC-validation.py -e /lscratch/$SLURM_JOBID -f $WELLFDRLIMIT > $EXPERIMENT/AUC.txt
    rm -R /lscratch/$SLURM_JOBID/ChIP-validation

fi

bedops -e 1 /lscratch/$SLURM_JOBID/footprints_R1.bed /lscratch/$SLURM_JOBID/footprints_R2.bed > /lscratch/$SLURM_JOBID/overlap.bed
overlap=($(wc -l /lscratch/$SLURM_JOBID/overlap.bed))
original=($(wc -l /lscratch/$SLURM_JOBID/footprints_R1.bed))
echo print $overlap/$original. | python > $EXPERIMENT/overlap.txt

