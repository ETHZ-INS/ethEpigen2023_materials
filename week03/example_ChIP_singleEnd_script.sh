#!/bin/bash

ref=/reference/Mus_musculus/Ensembl/GRCm39/Sequence/BOWTIE2Index/genome
adapters=/conda/share/trimmomatic-0.38-1/adapters/TruSeqPE.fa

mkdir -p trimmed
mkdir -p aligned
mkdir -p tracks

for f in raw/SRR*.fastq.gz; do # begin of loop

base=`basename $f .fastq.gz`
bam=aligned/"$base".bam
echo $base

# TRIMMING

trimdir=trimmed

if [ -f "$trimdir/$base.fastq.gz" ]; then
  echo $trimdir/$base".fastq.gz found; skipping"
else
trimmomatic SE -threads 6 -summary $trimdir/"$base".stats -phred33 $f \
  $trimdir/$base.fastq.gz ILLUMINACLIP:$adapters:2:15:4:4:true \
  LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25
fi

# ALIGNMENT

if [ -f "$bam" ]; then
    echo "$bam found; skipping"
else
(bowtie2 -p 6 -x $ref -U $trimdir/$base.fastq.gz) 2> aligned/$base.bowtie2 |\
  samtools view -bS - | samtools sort -@4 -m 2G - > $bam
samtools index $bam
fi

# COVERAGE TRACKS

if [ -f "tracks/$base.bw" ]; then
    echo "tracks/$base.bw found; skipping"
else
  bamCoverage -p 6 --ignoreDuplicates --effectiveGenomeSize 2652783500 --normalizeUsing CPM \
    -b aligned/"$base".bam -o tracks/$base.bw
fi

done  # end of loop

