# Step 9: Differential expression analysis using DESeq2 and GEOquery

# Load the required libraries
library(DESeq2)
library(GEOquery)

# Load counts data
counts <- read.table("counts.txt", header = TRUE, row.names = 1)

# Remove the first five columns which are not sample columns
counts <- counts[ , -(1:5)]

# Download metadata from GEO using accession number GSE161199
geo_data <- getGEO("GSE161199", GSEMatrix = TRUE)
metadata <- pData(phenoData(geo_data[[1]]))

# Ensure that the sample names in the metadata match the column names in counts
colnames(counts) <- metadata$geo_accession

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ condition)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results with adjusted p-values
res <- results(dds)

# Step 10: Output results
# Order results by adjusted p-value
resOrdered <- res[order(res$padj), ]

# Write results to a CSV file
write.csv(as.data.frame(resOrdered), file = "deseq2_results.csv")

# Step 11: Create volcano plot and MA plot for DEGs

# Volcano plot
library(ggplot2)
res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
res$padj[is.na(res$padj)] <- 1
ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  ggtitle("Volcano Plot")

# MA plot
plotMA(res, main="DESeq2 MA Plot")

# Step 12: List top 10 Differentially Expressed Genes
top_genes <- head(resOrdered, 10)
print(top_genes)
