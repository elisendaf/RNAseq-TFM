## ----EXPLORATORY ANALYSIS-----------------------------------------------------


# 1. Load the counts file
counts_data <- read.csv("./gene_counts2.csv", row.names = 1)


# Modify the matrix of the counts
conditions <- factor(c(rep("Control", 3), rep("Salinity", 3), rep("S+H", 3), rep("Heat", 3)))
coldata <- data.frame(condition = conditions, row.names = colnames(counts_data))


# 2. Creation of a DESeqDataSet object
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = counts_data, 
                              colData = coldata, 
                              design = ~ condition)

# 3. Apply a VST transformation to normalize data
vsd <- vst(dds, blind = TRUE)

# A. PCA plot
plotPCA(vsd, intgroup = "condition")

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

library(ggplot2)
ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("#DDB1F3", "#F89C74", "#F6CF71","#66C5CC")) +
  labs(color = "Condition")+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("Principal Component Analysis") +
  theme(plot.title = element_text(hjust = 0.5, size = 16)) 
  


# B. Heatmap of most variable genes
library(pheatmap)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 30)
mat <- assay(vsd)[topVarGenes, ]
custom_colors <- list(condition = c("Control" = "#DDB1F3", "Heat" = "#F89C74", 
                "S+H" = "#F6CF71", "Salinity" = "#66C5CC"))

pheatmap(mat, 
         scale = "row", 
         annotation_col = coldata, 
         annotation_colors = custom_colors,
         main = "Heatmap of most variable genes")


# C. Dendrogram
dists <- dist(t(assay(vsd)))  #compute Euclidean distances based on VST transformed data
plot(hclust(dists), main = "Dendrogram of sample distances")


# D. Heatmap of sample to sample distances
sampleDistMatrix <- as.matrix(dists)  #Converts dists object into a matrix

pheatmap(sampleDistMatrix, 
         clustering_distance_rows = dists, 
         clustering_distance_cols = dists, 
         main = "Sample-to-Sample Distance Heatmap",
         annotation_col = coldata,
         annotation_colors = custom_colors)





## ----DIFFERENTIAL EXPRESSION ANALYSIS-----------------------------------------

# Use the DESeq function
ddsDE <- dds[rowSums(counts(dds)) > 10,]  #Filter low number of counts
ddsDE <- DESeq(dds)


# Extracion of results from the DESeq analysis

# List of comparisons between conditions
comparisons <- list(
  c("Salinity", "Control"),
  c("S+H", "Control"),
  c("Heat", "Control"),
  c("S+H", "Salinity"),
  c("Heat", "Salinity"),
  c("Heat", "S+H")
)

# Define the directories for results
DEG_dir <- "./DEGs2"
filteredDEG_dir <- "./filtered_DEGs2"
volcano_dir <- "./volcano_plots2"


# Iterate through each pair of comparisons
for (comp in comparisons) {
  contrast <- c("condition", comp[1], comp[2])
  
# Extraction of the results of DESeq2
counts_res <- results(ddsDE, contrast = contrast)
counts_res <- counts_res[order(counts_res$padj, na.last = NA),]
  
# Filtering of differentially expressed genes
filtered_degs <- counts_res[!is.na(counts_res$padj) & counts_res$padj <= 0.05 & abs(counts_res$log2FoldChange) > 1, ]
  
# Save the DEGs
deg_file <- file.path(DEG_dir, paste0("DEG2_", comp[1], "_vs_", comp[2], ".csv"))
write.csv(as.data.frame(counts_res), deg_file)

# Save the filtered DEGs
filtered_file <- file.path(filteredDEG_dir, paste0("Filtered_DEG2_", comp[1], "_vs_", comp[2], ".csv"))
write.csv(as.data.frame(filtered_degs), filtered_file)
  
# Volcano plot
volcano <- ggplot(as.data.frame(counts_res), aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), size = 1) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  theme(legend.position = "none") +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted p-value") +
  ggtitle(paste(comp[1], "vs", comp[2]))
  


# Guardar el gràfic individualment
volcano_plots <- file.path(volcano_dir, paste0("Volcano2_", comp[1], "_vs_", comp[2], ".png"))
ggsave(volcano_plots, plot = volcano, width = 6, height = 4)


}


