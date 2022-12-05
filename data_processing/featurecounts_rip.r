library("Rsubread")

args = commandArgs(trailingOnly=TRUE)
samplesheet <- args[1] #csv with metadata including file paths
sample_name <- args[2]
gtf <- args[3] #gene annotation file (gencode.v39.annotation.gtf)
outdir <- args[4]

annot_file <- gtf
filetable <- read.csv(samplesheet, header = TRUE)
filelist <- as.vector(filetable$filename)

GTFfeatureType="gene" #count over gene bodies
GTFattrType = "gene_id" 
strandSpec = 1 #0 is unstranded (default), 1 is stranded, 2 is reverse stranded

coverage <- featureCounts(files=filelist,
                          annot.ext=annot_file,
                          isGTFAnnotationFile=TRUE,
                          useMetaFeatures=TRUE,
                          GTF.featureType=GTFfeatureType,
                          GTF.attrType=GTFattrType,
                          allowMultiOverlap=TRUE,
                          largestOverlap=TRUE,
                          isPairedEnd=TRUE,
                          requireBothEndsMapped=FALSE,
                          strandSpecific=strandSpec,
                          nthreads=8
                          )

colnames(coverage$counts) <- filetable$sampID

#save results
fileroot <- paste0(outdir, "/featureCounts_",sample_name)

save.image(paste0(fileroot, ".RData"))

write.csv(coverage$counts, paste(fileroot,".coverage.csv", sep=""))
write.csv(coverage$stat, paste(fileroot,".stat.csv", sep=""))
write.csv(coverage$annotation, paste(fileroot,".annotation.csv", sep=""))
write.csv(coverage$targets, paste(fileroot,".targets.csv", sep=""))