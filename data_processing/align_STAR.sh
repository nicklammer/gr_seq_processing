#!/bin/bash
#align reads to genome with STAR

indir=/scratch/Users/nila7826/4suseq/trim
outdir=/scratch/Users/nila7826/4suseq/align
genome=hg38
#genome index was generated using GRCh38.primary_assembly.genome.fa
#from gencode
index=gencode.v39.pri.index
indexdir=/scratch/Users/nila7826/genome/${genome}/${index}
read1=$1
read2=$2
name=$3
#read1 and read2 are paired-end reads
#software versions:
#STAR v2.7.3a
#samtools v1.10

STAR \
--genomeDir ${indexdir} \
--runThreadN 8 \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${outdir}/${name}. \
--readFilesIn ${indir}/${read1} ${indir}/${read2}

samtools index \
-@8 \
${outdir}/${name}.Aligned.sortedByCoord.out.bam