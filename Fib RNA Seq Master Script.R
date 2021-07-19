setwd('C:/Users/cmorrow5/Desktop/Fib RNA Seq/2021-01-08_Chris_FibroRNAseq/DEseq2_RUV/')

options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library(RUVSeq)
library(bitops)
#library(org.Mm.eg.db)
#library(clusterProfiler)

# 2021-01-08
# analyze Chris data

# library is a second strand first (RF) format

####################################################################
# get output file prefixes
my.outprefix <- paste(Sys.Date(),"Moorelab_Fibroblasts_DESeq2_model",sep="_")

# read kallisto counts
my.rna.data <- read.table("../kallisto/2021-01-08_Chris_RNAseq_Fibroblasts_with_ERCC.txt", sep = "\t", header = T)
my.fibro <- my.rna.data[,-c(2:3)]

# aggregate over trafibroripts to get counts per gene and round values
my.fibro.agg <- aggregate(my.fibro[,-1], by = list(my.fibro$GeneSymbol), FUN = 'sum')
my.fibro.agg[,-1] <- round(my.fibro.agg[,-1]) # need integers for DEseq2

colSums(my.fibro.agg[,-1])
#ID1_fibroblasts.kallisto_res.abundance.tsv ID2_fibroblasts.kallisto_res.abundance.tsv 
#21866383                                   21359592 
#ID3_fibroblasts.kallisto_res.abundance.tsv IM1_fibroblasts.kallisto_res.abundance.tsv 
#18710969                                   22778339 
#IM2_fibroblasts.kallisto_res.abundance.tsv IM3_fibroblasts.kallisto_res.abundance.tsv 
#21229449                                   22064456 
#PD1_fibroblasts.kallisto_res.abundance.tsv PD2_fibroblasts.kallisto_res.abundance.tsv 
#19870274                                   20335478 
#PD3_fibroblasts.kallisto_res.abundance.tsv PM1_fibroblasts.kallisto_res.abundance.tsv 
#20749155                                   17817653 
#PM2_fibroblasts.kallisto_res.abundance.tsv PM3_fibroblasts.kallisto_res.abundance.tsv 
#17728202                                   20668309 
#QD1_fibroblasts.kallisto_res.abundance.tsv QD2_fibroblasts.kallisto_res.abundance.tsv 
#16688167                                   19670770 
#QD3_fibroblasts.kallisto_res.abundance.tsv QM1_fibroblasts.kallisto_res.abundance.tsv 
#20012479                                   24036496 
#QM2_fibroblasts.kallisto_res.abundance.tsv QM3_fibroblasts.kallisto_res.abundance.tsv 
#22916787                                   21228783 
#SD1_fibroblasts.kallisto_res.abundance.tsv SD2_fibroblasts.kallisto_res.abundance.tsv 
#23390491                                   23757453 
#SD3_fibroblasts.kallisto_res.abundance.tsv SM1_fibroblasts.kallisto_res.abundance.tsv 
#21201862                                   23430989 
#SM2_fibroblasts.kallisto_res.abundance.tsv SM3_fibroblasts.kallisto_res.abundance.tsv 
#21297131                                   20336591 

# see deseq2 vignette, remove genes without consistent expression
my.keep <- apply(my.fibro.agg[,-1]> 0, 1, sum) > 6

# Now pull out the null/low expressed genes
my.filtered.matrix           <- my.fibro.agg[my.keep,-1] # 20045
rownames(my.filtered.matrix) <- my.fibro.agg[my.keep,1]

# clean up column names
colnames(my.filtered.matrix) <- unlist(strsplit(colnames(my.filtered.matrix), ".kallisto_res.abundance.tsv"))


a <- log2(my.filtered.matrix[grep("ERCC-",rownames(my.filtered.matrix)),]+0.1)
head(a)
b <- a[,c(7,8,9,10,11,12,19,20,21,22,23,24,13,14,15,16,17,18,1,2,3,4,5,6)]          
head(b)

#setwd("C:/Users/cmorrow5/Desktop/")
pdf(paste0(my.outprefix,"_ERCC_counts_preNorm.pdf"))
par(oma=c(3,0.5,0.5,0.5))

boxplot(b,
        ylab = "ERCC counts (log2)",
        col=c(rep('tomato',3),rep('gold',3), rep('pink',3),rep('purple',3),rep('dodgerblue1',3),rep('dodgerblue4',3), rep('olivedrab1',3),rep('olivedrab4',3)),
        las = 2
        )
dev.off()

## Run RUV normalization using spike ins
my.fibro.RUV <- RUVg(as.matrix(my.filtered.matrix), grep("ERCC-",rownames(my.filtered.matrix) ), k=2)

# check the RUV values
my.fibro.RUV$W
cbind(colnames(my.filtered.matrix),my.fibro.RUV$W) # QM has lower RNA for sure


#setwd("C:/Users/cmorrow5/Desktop/")
pdf(paste0(my.outprefix,"_ERCC_counts_postNorm_RUVg.pdf"))
par(oma=c(3,0.5,0.5,0.5))
boxplot(log2(my.fibro.RUV$normalizedCounts[,c(7,8,9,10,11,12,19,20,21,22,23,24,13,14,15,16,17,18,1,2,3,4,5,6)]+0.1),
        col=c(rep('tomato',3),rep('gold',3), rep('pink',3),rep('purple',3),rep('dodgerblue1',3),rep('dodgerblue4',3), rep('olivedrab1',3),rep('olivedrab4',3)),
        cex=0.5,ylab="Log2 RUV Normalized counts", las = 2)  
dev.off()

# plot correlation for all pairs
jpeg(paste0(my.outprefix,"_Correlations_samples.jpg"), height = 20, width = 20, units = "cm", res = 100)
pairs(log2(my.fibro.RUV$normalizedCounts+0.1), log = 'xy')
dev.off()

write.table(my.fibro.RUV$normalizedCounts, file = paste0(my.outprefix,"_RUVg_ERCC_normalized_counts.txt"), sep = "\t",col.names = T, quote = F)

####################################################################
my.status <- factor(c(rep('ID',3),rep('IM',3),rep('PD',3),rep('PM',3),rep('QD',3),rep('QM',3),rep('SDA',3),rep('SM',3)))
my.colors <- c(rep('tomato',3),rep('gold',3), rep('pink',3),rep('purple',3),rep('dodgerblue1',3),rep('dodgerblue4',3), rep('olivedrab1',3),rep('olivedrab4',3))

# build design matrix
dataDesign = data.frame(row.names = colnames( my.fibro.RUV$normalizedCounts ), 
                        status = my.status) 

# get matrix using age as a modeling covariate
dds <- DESeqDataSetFromMatrix(countData = my.fibro.RUV$normalizedCounts,
                              colData = dataDesign,
                              design = ~ status )

# run DESeq normalizations and export results
dds.deseq <- DESeq(dds)

# plot dispersion
my.disp.out <- paste(my.outprefix,"_dispersion_plot.pdf")

pdf(my.disp.out)
plotDispEsts(dds.deseq)
dev.off()

# get DESeq2 normalized expression value
tissue.cts <- getVarianceStabilizedData(dds.deseq)

#setwd("C:/Users/cmorrow5/Desktop/")
#write.table(tissue.cts, file = "123flamingo", sep = "\t",col.names = T, quote = F)

# expression range
pdf(paste(my.outprefix,"_ERCC_VST_Normalized_counts_boxplot.pdf"))
boxplot(tissue.cts,col=my.colors,cex=0.5,ylab="Log2 DESeq2 Normalized counts", las = 2)  
dev.off()

# check sex of animals
if( sum(rownames(tissue.cts)  %in% "Xist")  != 0) {
  # plot Xist expression
  pdf(paste(my.outprefix,"_Normalized_counts_Xist_expression_barplot.pdf"))
  barplot(tissue.cts["Xist",], ylab = "Normalized log2(counts) Xist expression", las = 2, col = my.colors)
  dev.off()
}


my.pv <- pvclust::pvclust(tissue.cts,nboot = 100)

pdf(paste(my.outprefix,"_clustering_pvclust_100boots.pdf"))
plot(my.pv)
dev.off()

# MDS analysis
mds.result <- cmdscale(1-cor(tissue.cts,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

pdf(paste(my.outprefix,"_MDS_plot.pdf"))
plot(x, y, 
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="Multi-dimensional Scaling", 
     cex=3, col= my.colors, pch= 16,
     cex.lab = 1.5,
     cex.axis = 1.5)
text(x, y, labels =colnames(tissue.cts))
dev.off()

# PCA analysis
my.pos.var <- apply(tissue.cts,1,var) > 0
my.pca <- prcomp(t(tissue.cts[my.pos.var,]),scale = TRUE)
x <- my.pca$x[,1]
y <- my.pca$x[,2]

my.summary <- summary(my.pca)

pdf(paste(my.outprefix,"_PCA_plot.pdf",sep=""))
plot(x,y, 
     pch = 16, cex=3, col= my.colors, 
     xlab = paste('PC1 (', round(100*my.summary$importance[,1][2],1),"%)", sep=""),
     ylab = paste('PC2 (', round(100*my.summary$importance[,2][2],1),"%)", sep=""),
     cex.lab = 1.5,
     cex.axis = 1.5) 
dev.off()
##########################################################################################

sessionInfo()



################################################################################################
################################################################################################
################################################################################################
################################################################################################
##Pathway Analyses##

##########################################################
######## A. GO Gene Set Enrichment Analysis
go.all.gsea <- gseGO(geneList     = nscAct.geneList,
                     OrgDb        = org.Mm.eg.db,
                     keyType      = "ENTREZID",
                     ont          = "ALL",
                     nPerm        = 1000,
                     minGSSize    = 100,   
                     maxGSSize    = 2500,  
                     pvalueCutoff = 0.05,
                     verbose      = FALSE)



# GSEA plot of interesting pathways
gseaplot(go.all.gsea, geneSetID = "GO:0000165", title = "GO:0000165 MAPK")                      
gseaplot(go.all.gsea, geneSetID = "GO:0000502", title = "GO:0000502 Proteasome Complex") 
gseaplot(go.all.gsea, geneSetID = "GO:0051087", title = "GO:0051087 Chaperone binding")
gseaplot(go.all.gsea, geneSetID = "GO:0030424", title = "GO:0030424 Axon")





# write results to file
write.table(go.all.gsea@result, file = paste(Sys.Date()," ",condition, " ", "Activation_GO_GSEA_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")
my.go.all.gsea.res <- go.all.gsea@result

# make some plots
#pdf(paste0(my.outprefix," ", condition, " ", "GO BP-Dot-Plot.pdf"))
par(oma=c(3,0.5,0.5,0.5))
dotplot(go.all.gsea,  x = "Count", paste("Significant GO Terms", condition))
#dev.off()

#cnetplot(go.all.gsea,  categorySize="pvalue", node_label = 'category')
#emapplot(go.all.gsea)

##########################################################
######## B. KEGG Gene Set Enrichment Analysis
kegg.gsea <- gseKEGG(geneList     = nscAct.geneList,
                     organism     = 'mmu',
                     keyType      = 'ncbi-geneid',
                     nPerm        = 1000,
                     pvalueCutoff = 0.05,
                     verbose      = FALSE)

# write results to file
write.table(kegg.gsea@result, file = paste(Sys.Date()," ",condition, " ", "Activation_KEGG_GSEA_Analysis_FDR5.txt", sep = "_"), quote = F, sep = "\t")
my.kegg.gsea.res <- kegg.gsea@result


# GSEA plot of interesting pathways
gseaplot(kegg.gsea, geneSetID = "mmu03050", title = "Proteasome")                 
gseaplot(kegg.gsea, geneSetID = "mmu04142", title = "Lysosome") 

#pdf(paste0(my.outprefix," ", condition, " ", "KEGG-Dot-Plot.pdf"))
par(oma=c(3,0.5,0.5,0.5))
dotplot(kegg.gsea,  x = "Count", title = paste("Significant KEGG Terms", condition))
#dev.off()

#emapplot(kegg.gsea)
