#!/bin/bash
#tag read heads with UMIs, trim and fastqc with Trim Galore

indir=/scratch/Users/nila7826/4suseq/raw
outdir=/scratch/Users/nila7826/4suseq/trim
read1=$1
read2=$2
umi=$3
name=$4
#read1 and read2 are paired-end reads
#umi is the read file
#software versions:
#python v3.7.4
#UMI-tools v1.1.2
#Trim Galore v0.6.6

umi_tools extract \
--bc-pattern=NNNNNNNN \
--stdin=${indir}/${umi} \
--read2-in=${indir}/${read1} \
--stdout=${indir}/${name}_R1.fastq.gz \
--read2-stdout

umi_tools extract \
--bc-pattern=NNNNNNNN \
--stdin=${indir}/${umi} \
--read2-in=${indir}/${read2} \
--stdout=${indir}/${name}_R2.fastq.gz \
--read2-stdout

/Users/nila7826/apps/TrimGalore-0.6.6/trim_galore \
-j 8 \
--fastqc \
--2colour 20 \
--paired \
-o ${outdir} \
${indir}/${name}_R1.fastq.gz \
${indir}/${name}_R2.fastq.gz