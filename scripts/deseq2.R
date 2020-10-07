# library(here)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)

# import snakemake variables ------------------------------------------
catalog <- snakemake@input[['catalog']]
md <- snakemake@params[['metadata']]

print(paste("Peaks catalog: ", catalog))

print(paste("Metadata file: ", md))
# import data ---------------------------------------------------------
catalog <- read.table(catalog, sep = "\t", header = T, stringsAsFactors = F, check.names = F)
md <- read.table(md, sep = "\t", header = T, stringsAsFactors = F)

print(head(catalog))
print(head(md))

counts <- catalog[,md$SampleID]
rownames(counts) <- catalog$V4
rownames(md) <- md$SampleID

# deseq2 --------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = md, design = ~Condition)
dds <- DESeq(dds)
res <- results(dds)

# export normalized data
NormCounts <-counts(dds, normalized=TRUE)
NormCounts <- 0.00001+(NormCounts) #this makes it so there are no 0 values#
LogNormCounts <- log2(NormCounts)
LogNormCounts <- as.data.frame(LogNormCounts)
LogNormCounts$ID <- rownames(LogNormCounts)

write.csv(NormCounts, snakemake@output[['norm_counts']], quote = F)
write.csv(LogNormCounts, snakemake@output[['log_norm_counts']], quote = F)
print("Exported normalized counts tables")
# contrasts

# clustering
# pca
vsd <- vst(dds, blind = FALSE)
pca <- plotPCA(vsd, intgroup = "Condition") + geom_label_repel(aes(label=name))
ggsave(snakemake@output[['pca']], pca, width = 16, height = 9, dpi = 300, units = "in")

# heatmap
# sampleDists <- dist(t(assay(vsd)))
# sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- paste(vsd$Condition, vsd$type, sep="-")
# colnames(sampleDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# pheatmap(sampleDistMatrix,
# 		clustering_distance_rows=sampleDists,
# 		clustering_distance_cols=sampleDists,
# 		col=colors,
# 		filename = snakemake@output[['heatmap']],
# 		width = 9,
# 		height = 9)
print("Exported clustering figures")