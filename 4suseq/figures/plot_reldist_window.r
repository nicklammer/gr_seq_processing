#plot window and reldist analyses together
library("dplyr")
library("ggplot2")
library("gridExtra")

#reldist files in a csv with sample names
distfiles <- read.csv("/scratch/Users/nila7826/4suseq/chip/results/samples_reldist.csv",stringsAsFactors=FALSE)
#bed files from bedtools window analysis in csv with sample names and shuffle files
beds <- read.csv("/scratch/Users/nila7826/4suseq/chip/results/samples_window_30000.csv",stringsAsFactors=FALSE)
windowsize <- "30 kb"
filename <- "/Users/nila7826/figs/reldist_window_30000.png"

#reldist
#concatenate reldist data for each sample and assign sample names
distdata <- data.frame(reldist=character(),fraction=character(),geneSet=character())
for (i in 1:length(distfiles$reldist)) {
  distdf <- read.table(distfiles$reldist[i],sep='\t',header = TRUE)
  newdata <- distdf %>% select(reldist,fraction) %>% mutate(geneSet=distfiles$set[i])
  distdata <- rbind(distdata, newdata)
}

#window
allratios <- c()
ngenes <- c()
#calculate ratio of genes in window bed file vs all gene bed file
for (i in 1:length(beds$tss)) {
  reg <- length(readLines(beds$regular[i]))
  shuff <- length(readLines(beds$shuffle[i]))
  tss <- length(readLines(beds$tss[i]))
  regratio <- reg/tss
  shuffratio <- shuff/tss
  allratios <- append(allratios, regratio)
  allratios <- append(allratios, shuffratio)
  ngenes <- append(ngenes, paste0("total=",toString(tss)))
  ngenes <- append(ngenes, '')
}
geneSets <- c()
for (i in 1:length(beds$set)) {
  geneSets <- append(geneSets, rep(beds$set[i],2))
}
type <- rep(c("Regular","Shuffle"),length(beds$set))
windata <- data.frame(ratio=allratios,geneSet=geneSets,ngenes=ngenes,type=type)
obj <- ggplot(windata, aes(fill=type, y=ratio, x=geneSet)) +
  geom_bar(position="dodge", stat="identity")
yrange <- ggplot_build(obj)$layout$panel_params[[1]]$y.range

#plot
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
textsize=10
pdist <- ggplot(data=distdata,aes(x=reldist,y=fraction)) +
  geom_line(aes(color=geneSet)) +
  geom_point(aes(color=geneSet)) +
  scale_color_manual(values=cbPalette) +
  xlab("Relative distance between gene TSSs and nearest ChIP peaks") +
  ylab("Distribution") +
  theme_bw() +
  labs(color="") +
  theme(legend.title = element_blank(),
        text = element_text(size = textsize),
        axis.text = element_text(color="black"),
        legend.position=c(.97,.95),
        legend.justification=c("right","top"),
        legend.spacing.y = unit(0, "mm"),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "gray"),
        legend.text=element_text(size=textsize-1)) +
  guides(color = guide_legend(byrow = TRUE, keyheight=0.7)) 
pwin <- ggplot(windata, aes(fill=type, y=ratio, x=geneSet)) + 
  geom_bar(position="dodge", stat="identity") +
  #geom_text(aes(label=ngenes), vjust = -0.2, hjust='right') +
  ylim(yrange[1],yrange[2]+0.01) +
  xlab("") +
  ylab(paste0("Ratio of genes with nearby\nGR ChIP peak (TSS ± ",windowsize,")")) +
  theme_bw() +
  labs(fill="") +
  theme(legend.title = element_blank(),
        text = element_text(size = textsize),
        axis.text = element_text(color="black"),
        legend.position=c(.97,.95),
        legend.justification=c("right","top"),
        legend.spacing.y = unit(0, "mm"),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "gray"),
        legend.text=element_text(size=textsize-1))
pboth <- grid.arrange(pdist, pwin, nrow = 2)
mm <- 0.0393701
figwidth <- 120*mm
figheight <- 120*mm
ggsave(file=filename, plot=pboth, dpi=300, width=figwidth, height=figheight, scale=1)