#!/bin/bash
#bedtools window and reldist function to analyze overlap
#between gene TSSs and chip peaks
#return TSS that has nearby chip peak

genedir=/scratch/Users/nila7826/4suseq/chip/beds
chipdir=/scratch/Users/nila7826/4suseq/chip/beds
genebed=$1
shuffbed=$2
window=$3
outname=$4
outdir=/scratch/Users/nila7826/4suseq/chip/results/${outname}.${window}
genepath=${genedir}/${genebed}
shuffpath=${genedir}/${shuffbed}
chip=${chipdir}/stallcup_lis_dex.intersect.bed
genomebed=/scratch/Users/nila7826/genome/hg38/hg38.genome
#genebed is the gene TSS bed file
#shuffbed is the shuffled version of genebed
#window is an integer window size
#chip is the bed file with chip peaks
#software versions:
#bedtools v2.30.0

bedtools window \
-a ${genepath} \
-b ${chip} \
-w ${window} \
-u \
> ${outdir}/${outname}.chip.${window}.bed

bedtools window \
-a ${shuffpath} \
-b ${chip} \
-w ${window} \
-u \
> ${outdir}/${outname}.chip.shuffle.${window}.bed

bedtools reldist \
-a ${genepath} \
-b ${chip} \
> ${outdir}/${outname}.chip.dist

bedtools reldist \
-a ${shuffpath} \
-b ${chip} \
> ${outdir}/${outname}.chip.shuffle.dist