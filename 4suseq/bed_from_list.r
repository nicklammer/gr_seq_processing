#get overlapping regions from list of genes (ensembl id with version number)
library("dplyr")

args = commandArgs(trailingOnly=TRUE)
#bed file with gene regions and matching gene names
#typically used with a TSS bed file, gene name in fourth field
bedfilepath <- args[1]
#list of genes with names matching bed file sep by new line
genelistpath <- args[2]
outname <- args[3]

bedfile <- read.table(bedfilepath,sep='\t',header = FALSE)
genelist <- read.csv(genelistpath, header=FALSE)
newbed <- bedfile %>% filter(V4 %in% genelist$V1)

write.table(newbed,file=outname,sep='\t',row.names = FALSE,col.names = FALSE,quote=FALSE)