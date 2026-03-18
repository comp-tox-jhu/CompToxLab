# Pathway Analysis

## What is pathway or functional enrichment analysis?

Pathways analysis helps us understand the biological functions and networks affected by your differentially expressed genes (DEGs). Instead of looking at individual genes in isolation, we look at groups of genes that work together in known biological processes.


## How is this done?

Typically, you examine the overlap between your DEGs and sets of genes with known biological functions. If the overlap is greater than expected by chance you can say that your DEGs are enriched for genes related to that pathway. You can compare your DEGs to different databases of gene sets (KEGG, GO, Reactome, etc.). At it's simplist implementation you can just examine the overlap between your DEGs and your own gene set if you have one!

In the differential expression tutorial we looked at DEGs between AD and control, then we ended with GO enrichment results. Let's dig a little deeper into those results.

## Overview of what the data looks like

Before we dive into analysis, let's understand what we're looking at. From your differential expression analysis, you should have a results table (`res_mapped`) with the following columns:

- **GeneID**: The NCBI gene identifier
- **Symbol**: The gene name (e.g., GABRB2, GRIN1)
- **log2FoldChange**: How much the gene is up or down regulated (positive = up in AD, negative = down in AD)
- **padj**: The adjusted p-value (< 0.05 is significant)
- **baseMean**: Average expression level

Your significant DEGs likely show one of two patterns:

**Upregulated genes** (log2FoldChange > 0): More active in your condition of interest (e.g., AD)
**Downregulated genes** (log2FoldChange < 0): Less active in your condition of interest

## Loading/Cleaning Data

```{r}
library(tidyverse)         # data manipulation/plotting
library(DESeq2)            # differential expression
library(EnhancedVolcano)   # creating volcano plots
library(ggpubr)            # publication ready plotting
library(clusterProfiler)   # enrichment analysis 
```
