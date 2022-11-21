#deseq for GR mutant 4sU-seq data
library("DESeq2")

#file with metadata
metadata_file <- "./4suseq/samples.all.csv"
#counts from featureCounts
counts_file <- "./4suseq/featureCounts_4su_exon_stranded.coverage.csv"
#file with ensembl ID to gene symbol matches (optional)
ids <- read.csv("./geneids.gencode.v39.csv")

metadata <- read.csv(metadata_file)
metadata$reps <- factor(rep(c(1,1,1,1,2,2,2,2),3))#assign replicate numbers
metadata$time <- factor(metadata$time)
metadata$mutant <- factor(metadata$mutant)
metadata$mutant <- relevel(metadata$mutant,"wt")


counts <- read.csv(counts_file, row.names = 1)
metadata$sample <- colnames(counts)
rownames(metadata) <- metadata$sample

full_model <- ~mutant+time+mutant:time
dds <- DESeqDataSetFromMatrix(counts,metadata,design=full_model)
DEdds <- DESeq(dds)
#results for diff expression over time
#change dex time for desired hour of treatment (1-3)
dextime <- 1 #treatment time vs 0
timename <- paste0("time_",dextime,"_vs_0")
softimename <- paste0("mutantrna.time",dextime)
ctrltimename <- paste0("mutantdna.time",dextime)

wtres <- results(DEdds, contrast=c("time",dextime,0))
sofres <- results(DEdds, list(c(timename,softimename)))
ctrlres <- results(DEdds, list(c(timename,ctrltimename)))

#assign gene symbols
wtres$ens <- rownames(wtres)
sofres$ens <- rownames(sofres)
ctrlres$ens <- rownames(ctrlres)
idx <- match(wtres$ens, ids$ensembl_gene_id_version)
wtres$name <- ids$hgnc_symbol[idx]
sofres$name <- ids$hgnc_symbol[idx]
ctrlres$name <- ids$hgnc_symbol[idx]
#save results
write.csv(wtres, file=paste0("./final files/res_wt_",dextime,"h.csv"))
write.csv(sofres, file=paste0("./final files/res_sof_",dextime,"h.csv"))
write.csv(ctrlres, file=paste0("./final files/res_ctrl_",dextime,"h.csv"))

#results for anything up/down relative to wt
sofvswt <- results(DEdds, name="mutant_rna_vs_wt")
ctrlvswt <- results(DEdds, name="mutant_dna_vs_wt")
#assign gene symbols
sofvswt$ens <- rownames(sofvswt)
ctrlvswt$ens <- rownames(ctrlvswt)
idx <- match(sofvswt$ens, ids$ensembl_gene_id_version)
sofvswt$name <- ids$hgnc_symbol[idx]
ctrlvswt$name <- ids$hgnc_symbol[idx]
#save results
write.csv(sofvswt, file="./final files/res_sofvswt.csv")
write.csv(ctrlvswt, file="./final files/res_ctrlvswt.csv")

#save rlog transformed counts for plotting
rld <- rlog(DEdds)
rld_mat <- assay(rld)
write.csv(rld_mat, file="./final files/4su_counts_rlog.csv")

#plots pca with rlog counts
plotPCA(rld,intgroup="time")
plotPCA(rld,intgroup="mutant")
plotPCA(rld,intgroup="reps")

#plot p-value histograms
hist(wtres$pvalue)
hist(sofres$pvalue)
hist(ctrlres$pvalue)
hist(sofvswt$pvalue)
hist(ctrlvswt$pvalue)


