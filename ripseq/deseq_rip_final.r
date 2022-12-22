#deseq for GR RIP data
library("DESeq2")
library("ggplot2")

metadata_file <- "./ripseq/samples.all.csv"
counts_file <- "./ripseq/featureCounts_rip_gene_stranded.coverage.csv"
#file with ensembl ID to gene symbol matches
ids <- read.csv("./geneids.gencode.v39.csv")

metadata <- read.csv(metadata_file)
metadata$reps <- factor(rep(c(1,1,1,1,2,2,2,2),2))
metadata$time <- factor(metadata$time)
metadata$fraction <- factor(metadata$fraction)

counts <- read.csv(counts_file, row.names = 1)
metadata$sample <- colnames(counts)
rownames(metadata) <- metadata$sample

full_model <- ~fraction+time+fraction:time
reduced_model <- ~fraction+time

dds <- DESeqDataSetFromMatrix(counts,metadata,design=full_model)
DEdds <- DESeq(dds)
#results for enrichment at each time point
res0 <- results(DEdds, contrast=c("fraction","ip","input"))
res1 <- results(DEdds, list(c("fraction_ip_vs_input","fractionip.time1")))
res2 <- results(DEdds, list(c("fraction_ip_vs_input","fractionip.time2")))
res3 <- results(DEdds, list(c("fraction_ip_vs_input","fractionip.time3")))

#assign gene symbols
res0$ens <- rownames(res0)
res1$ens <- rownames(res1)
res2$ens <- rownames(res2)
res3$ens <- rownames(res3)
idx <- match(res0$ens, ids$ensembl_gene_id_version)
res0$name <- ids$hgnc_symbol[idx]
res1$name <- ids$hgnc_symbol[idx]
res2$name <- ids$hgnc_symbol[idx]
res3$name <- ids$hgnc_symbol[idx]

write.csv(res0, file=paste0("./final files/rip_res_0_20221109.csv"))
write.csv(res1, file=paste0("./final files/rip_res_1_20221109.csv"))
write.csv(res2, file=paste0("./final files/rip_res_2_20221109.csv"))
write.csv(res3, file=paste0("./final files/rip_res_3_20221109.csv"))

#save rlog transformed counts for plotting
rld <- rlog(DEdds)
rld_mat <- assay(rld)
write.csv(rld_mat, file="./final files/rip_counts_rlog_20221109.csv")

#plots pca with rlog counts
mm <- 0.0393701
fontsize <- 12
figwidth <- 150*mm
figheight <- 100*mm
pca_time <- plotPCA(rld,intgroup="time")
pca_time <- pca_time +
  theme(text = element_text(size = fontsize)) +
  labs(color = "Dex time\n(hours)") +
  theme_bw()
ggsave(filename="./final files/pca_rip_time.png",plot=pca_time,
       width=figwidth,height=figheight,device='png',dpi=300)
pca_mutant <- plotPCA(rld,intgroup="fraction")
pca_mutant <- pca_mutant +
  theme(text = element_text(size = fontsize)) +
  scale_color_discrete(labels=c('input','IP')) +
  labs(color = "RIP\nfraction") +
  theme_bw()
ggsave(filename="./final files/pca_rip_fraction.png",plot=pca_mutant,
       width=figwidth,height=figheight,device='png',dpi=300)
pca_reps <- plotPCA(rld,intgroup="reps")
pca_reps <- pca_reps +
  theme(text = element_text(size = fontsize)) +
  labs(color = "Replicate") +
  theme_bw()
ggsave(filename="./final files/pca_rip_reps.png",plot=pca_reps,
       width=figwidth,height=figheight,device='png',dpi=300)

#plot p-value histograms
mm <- 0.0393701
fontsize <- 12
figwidth <- 60*mm
figheight <- 40*mm
hist_0 <- ggplot(data.frame(res0), aes(x=pvalue)) + 
  geom_histogram(color="black", fill="grey") +
  ylab('Count') +
  xlab('p-value') +
  theme(text = element_text(size = fontsize)) +
  theme_bw()
ggsave(filename=paste0("./final files/hist_rip_0.png"),plot=hist_0,
       width=figwidth,height=figheight,device='png',dpi=300)
hist_1 <- ggplot(data.frame(res1), aes(x=pvalue)) + 
  geom_histogram(color="black", fill="grey") +
  ylab('Count') +
  xlab('p-value') +
  theme(text = element_text(size = fontsize)) +
  theme_bw()
ggsave(filename=paste0("./final files/hist_rip_1.png"),plot=hist_1,
       width=figwidth,height=figheight,device='png',dpi=300)
hist_2 <- ggplot(data.frame(res2), aes(x=pvalue)) + 
  geom_histogram(color="black", fill="grey") +
  ylab('Count') +
  xlab('p-value') +
  theme(text = element_text(size = fontsize)) +
  theme_bw()
ggsave(filename=paste0("./final files/hist_rip_2.png"),plot=hist_2,
       width=figwidth,height=figheight,device='png',dpi=300)
hist_3 <- ggplot(data.frame(res3), aes(x=pvalue)) + 
  geom_histogram(color="black", fill="grey") +
  ylab('Count') +
  xlab('p-value') +
  theme(text = element_text(size = fontsize)) +
  theme_bw()
ggsave(filename="./final files/hist_rip_3.png",plot=hist_3,
       width=figwidth,height=figheight,device='png',dpi=300)
