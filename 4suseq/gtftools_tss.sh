#!/bin/bash
#use GTFtools to generate TSS bed file from gencode gtf annotation

gtf=gencode.v39.annotation.gtf
outname=gencode.v39.tss.bed
indir=/scratch/Users/nila7826/genome/hg38
outdir=/scratch/Users/nila7826/genome/hg38
gtfpath=${indir}/${gtf}
outfile=${outdir}/${outname}

#software versions:
#python v3.6.3
#GTFtools v0.9.0

. /Users/nila7826/pyenv/gtftools/bin/activate
gtftools=/Users/nila7826/pyenv/gtftools/bin/gtftools

${gtftools} \
-t ${outfile} \
-w 0-1 \
${gtfpath}