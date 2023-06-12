rm(list=ls())
library(DESeq2)
library(pheatmap)
library(RUVSeq)
library(bitops)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(ggrepel)

setwd('C:/Users/cmorc/OneDrive - Harvard University/Unshared Compositions/Wisconsin Compositions/Fibroblast Paper/RNA Sequencing Analysis/Repository/')

options(stringsAsFactors = FALSE)

# library is a second strand first (RF) format
##Load raw counts with ERCC counts and perform normalizations and clean up data##

# get output file prefixes
my.outprefix <- paste(Sys.Date(),"Moorelab_Fibroblasts_DESeq2_model",sep="_")

# read kallisto counts
my.rna.data <- read.table("C:/Users/cmorc/OneDrive - Harvard University/Unshared Compositions/Wisconsin Compositions/Fibroblast Paper/RNA Sequencing Analysis/Repository/2021-01-08_Chris_RNAseq_Fibroblasts_with_ERCC.txt", sep = "\t", header = T)
my.fibro <- my.rna.data[,-c(2:3)]

# aggregate over transcripts to get counts per gene and round values, need integers for DEseq2
my.fibro.agg <- aggregate(my.fibro[,-1], by = list(my.fibro$GeneSymbol), FUN = 'sum')
my.fibro.agg[,-1] <- round(my.fibro.agg[,-1])

# see deseq2 vignette, remove genes without consistent expression
my.keep <- apply(my.fibro.agg[,-1]> 0, 1, sum) > 6

# Now pull out the null/low expressed genes
my.filtered.matrix           <- my.fibro.agg[my.keep,-1]
rownames(my.filtered.matrix) <- my.fibro.agg[my.keep,1]

# clean up column names
colnames(my.filtered.matrix) <- unlist(strsplit(colnames(my.filtered.matrix), ".kallisto_res.abundance.tsv"))

## Run RUV normalization using spike ins
my.fibro.RUV <- RUVg(as.matrix(my.filtered.matrix), grep("ERCC-",rownames(my.filtered.matrix) ), k=2)

# check the RUV values
my.fibro.RUV$W
cbind(colnames(my.filtered.matrix),my.fibro.RUV$W)

####################################################################
my.status <- factor(c(rep('ID',3),rep('IM',3),rep('PD',3),rep('PM',3),rep('QD',3),rep('QM',3),rep('SD',3),rep('SM',3)))
my.colors <- c(rep('tomato',3),rep('gold',3), rep('pink',3),rep('purple',3),rep('dodgerblue1',3),rep('dodgerblue4',3), rep('olivedrab1',3),rep('olivedrab4',3))

# build design matrix
dataDesign = data.frame(row.names = colnames( my.fibro.RUV$normalizedCounts ), 
                        status = my.status) 

dds <- DESeqDataSetFromMatrix(countData = my.fibro.RUV$normalizedCounts,
                              colData = dataDesign,
                              design = ~ status )

# run DESeq normalizations and export results
dds.deseq <- DESeq(dds)

# get DESeq2 normalized expression value
tissue.cts <- getVarianceStabilizedData(dds.deseq)
colnames(tissue.cts) <- unlist(strsplit(colnames(tissue.cts),".kallisto_res.abundance.tsv"))
resultsNames(dds.deseq)

num <- "PD"
dem <- "QD"

condition <- paste(num, " over ", dem, sep="")
res.act <- results(dds.deseq, contrast = c("status",num,dem))

res.act <- res.act[!is.na(res.act$padj),]
genes.act <- rownames(res.act)[res.act$padj < 0.05]

my.out.stats <- paste(my.outprefix,condition,".txt",sep = "_")
write.table(res.act, file = my.out.stats , sep = "\t" , row.names = T, quote=F)

my.out.fdr5 <- paste(my.outprefix,condition,"sig.txt",sep = "_")
write.table(res.act[genes.act,], file = my.out.fdr5, sep = "\t" , row.names = T, quote=F)
