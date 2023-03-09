#!/bin/bash

ref=/reference/Mus_musculus/Ensembl/GRCm39/Sequence/BOWTIE2Index/genome
adapters=/conda/share/trimmomatic-0.38-1/adapters/TruSeqPE.fa

mkdir -p trimmed
mkdir -p aligned
mkdir -p tracks

for f in raw/*_1.fastq.gz; do # begin of loop

base=`basename $f _1.fastq.gz`
bam=aligned/"$base".bam
echo $base

# TRIMMING

trimdir=trimmed

if [ -f "$trimdir/"$base"_1.paired.fastq.gz" ]; then
  echo $trimdir/$base"_*fastq.gz found; skipping"
else
trimmomatic PE -threads 6 -summary $trimdir/"$base".stats -phred33 \
  raw/"$base"_1.fastq.gz raw/"$base"_2.fastq.gz \
  $trimdir/"$base"_1.paired.fastq.gz $trimdir/"$base"_1.unpaired.fastq.gz \
  $trimdir/"$base"_2.paired.fastq.gz $trimdir/"$base"_2.unpaired.fastq.gz \
  ILLUMINACLIP:$adapters:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25
fi

# ALIGNMENT

if [ -f "$bam" ]; then
    echo "$bam found; skipping"
else
(bowtie2 -p 6 --dovetail --no-mixed --no-discordant -I 15 -X 2000 -x $ref \ 
  -1 $trimdir/"$base"_1.paired.fastq.gz -2 $trimdir/"$base"_2.paired.fastq.gz) 2> aligned/"$base".bowtie2 |\
  samtools view -bS - | samtools sort -@4 -m 2G - > $bam

java -jar /common/picard.jar MarkDuplicates I=$bam O=$base.bam.2 \
  M=aligned/$base.picard.dupMetrics.txt && mv $base.bam.2 $bam

samtools index $bam
fi

# COVERAGE TRACKS

if [ -f "tracks/$base.bw" ]; then
    echo "tracks/$base.bw found; skipping"
else
  bamCoverage -p 6 --ignoreDuplicates --effectiveGenomeSize 2652783500 --normalizeUsing CPM -b aligned/"$base".bam -o tracks/$base.bw
fi

done # end of loop





