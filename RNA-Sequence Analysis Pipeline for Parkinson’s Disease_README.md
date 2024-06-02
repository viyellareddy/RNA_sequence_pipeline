
# RNA-Sequence Analysis Pipeline for Parkinson’s Disease

## Introduction

This project implements an RNA-sequencing (RNA-seq) analysis pipeline to study differential gene expression in Parkinson's Disease (PD). The analysis utilizes RNA-seq data from the SRA database and follows a series of steps to identify differentially expressed genes, perform pathway enrichment analysis, and visualize the results.

## Abstract

Parkinson’s disease (PD) is a neurodegenerative disorder characterized by the loss of dopaminergic neurons and intracellular Lewy body deposits. This project leverages RNA-seq data from peripheral blood mononuclear cell (PBMC) samples of PD patients and healthy controls to elucidate the molecular mechanisms underlying PD. Using a comprehensive RNA-seq pipeline, 1,757 differentially expressed genes were identified, with significant pathway enrichment highlighting protein processing, mRNA processing, and mitochondrial dysfunction.

## Dataset

- Dataset Name: SRX9470334 (SRR13019394-SRR13019404)
- Source: SRA database
- Description: RNA-seq data from 6 PD patients and 6 healthy controls.

## Requirements

To run the RNA-seq pipeline, ensure you have the following tools and libraries installed:

```bash
module load sra-toolkit/2.10.8
module load samtools/1.9
module load fastqc/0.11.9
module load hisat2/2.1.0
module load subread/2.0.6
```

## Steps

### Step 1: open the termoinal
- .

### Step 2: Install Required Libraries
- Load the required libraries using the module commands provided above.

### Step 3: Download SRR Files
```bash
for i in {13019394..13019404}
do
  fastq-dump SRR${i}
done
```

### Step 4: Quality Control using FastQC
```bash
for file in *.fastq
do
  fastqc $file
done
```

### Step 5: Download Genome Indices for HISAT2
```bash
curl -O https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xvzf grch38_genome.tar.gz
```

### Step 6: Align Reads using HISAT2
```bash
for file in *.fastq
do
  hisat2 -x grch38_genome -U $file -S ${file%.fastq}.sam
  samtools view -bS ${file%.fastq}.sam > ${file%.fastq}.bam
  samtools sort ${file%.fastq}.bam -o ${file%.fastq}_sorted.bam
  rm ${file%.fastq}.sam
  rm ${file%.fastq}.bam
done
```

### Step 7: View BAM Files
```bash
for file in *_sorted.bam
do
  samtools view $file
done
```

### Step 8: Quantification using FeatureCounts
```bash
curl -O http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz
gunzip Homo_sapiens.GRCh38.106.gtf.gz
featureCounts -a Homo_sapiens.GRCh38.106.gtf -o counts.txt *_sorted.bam
```

### Step 9: Differential Expression Analysis using DESeq2
```r
# Load libraries
library(DESeq2)
library(GEOquery)

# Load counts data
counts <- read.table("counts.txt", header = TRUE, row.names = 1)
counts <- counts[ , -(1:5)]

# Download metadata
geo_data <- getGEO("GSE161199", GSEMatrix = TRUE)
metadata <- pData(phenoData(geo_data[[1]]))

# Ensure matching sample names
colnames(counts) <- metadata$geo_accession

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ condition)

# Run DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds)

# Order and save results
resOrdered <- res[order(res$padj), ]
write.csv(as.data.frame(resOrdered), file = "deseq2_results.csv")
```

### Step 10: Visualization and Analysis in R
#### Volcano Plot
```r
library(ggplot2)
res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
res$padj[is.na(res$padj)] <- 1
ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(alpha = 0.4) +
    theme_minimal() +
    ggtitle("Volcano Plot")
```

#### MA Plot
```r
plotMA(res, main="DESeq2 MA Plot")
```

### Step 11: Pathway Enrichment Analysis using RShinyGO
- Perform Gene Ontology (GO) and KEGG pathway enrichment analysis using the RShinyGO web interface: http://bioinformatics.sdstate.edu/go/

### Step 12: Top Differentially Expressed Genes
- The top differentially expressed genes identified include MATR3, MCUB, PLAUR, CLK1, UFL1, BIRC3, YAF2, RGPD5, ZRANB1, STRAP.

### Conclusion
This RNA-seq analysis pipeline effectively identifies differentially expressed genes and pathways associated with Parkinson's Disease. The findings highlight the utility of RNA sequencing in understanding the molecular disruptions underlying PD progression, providing insights into potential diagnostic and prognostic biomarkers.

