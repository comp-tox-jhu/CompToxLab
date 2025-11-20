## Selecting PCR Candidates from RNA-Seq Results

# Overview

Here we will describe how to select **qPCR validation targets** from RNA‑seq differential expression analysis. 

We will cover:

1. Loading RNA‑seq data  
2. Differential expression  
3. Enrichment analysis  
4. A step‑by‑step PCR‑target selection workflow:
   - Identify significant DEGs  
   - Ensure genes are well‑expressed  
   - Prioritize genes in enriched pathways  
   - Filter to protein coding genes
   - Check for adequate number exons and transcript length
   - Remove outlier‑driven changes  

Here we look to provide a biologically informed method to validate RNA‑seq findings.

---

# 1. Loading & Cleaning Data

```{r}
library(tidyverse)         
library(tidybulk)          
library(DESeq2)             
library(ggpubr)            
library(clusterProfiler)   
library(org.Hs.eg.db)      
library(biomaRt)
```

Download and load the SummarizedExperiment object:

```{r}
download.file("https://github.com/comp-tox-jhu/CompToxLab/raw/refs/heads/main/docs/omics/transcriptomics/pcr_validation/data/se.rds",
              destfile = "../data/se.rds")

se <- readRDS("../data/se.rds")
```

### What is a SummarizedExperiment

A SummarizedExperiment is a way of storing omics data:

![](img/se.svg){height=250px}

### Filtering and Stabilizing Expression

Before assessing our data for differentially expressed genes we should filter out genes that are lowly expressed to avoid noise.

```{r}
se <- se |> 
  identify_abundant() |> 
  keep_abundant()
```

Now to prepare for when we plot our data, let's add another assay, or counts matrix, that is variance stabilized. In this way we not only put our data in log scale but address issues with high variance genes tending to have higher means which can distort the expression trends we vizualize later.

```{r}
assay(se,"vst") <- vst(assay(se))
```

---

# 2. Differential Expression Analysis

We analyze females separately to detect sex‑specific transcriptional patterns. We then join the results with our `rowData`, which contains metadata about our genes, like the HGNC gene symbol.

```{r}
female_res <- se[,se$Sex=="F"] |> 
    test_differential_abundance(~condition, 
                                method = "deseq2",
                                fitType = "local",
                                action = "get") |> 
    as.data.frame() |> 
    rownames_to_column("GeneID") |> 
    inner_join(rowData(se) |> 
                 as.data.frame() |> 
                 mutate(GeneID=as.character(GeneID)), 
               by = "GeneID") |> 
    mutate(direction = ifelse(log2FoldChange > 0, "Up", "Down"))

```

### Significant DEGs

We first apply a significance filter (Adjusted p-value < 0.01, |log2FC| > log2(1)).

```{r}
sig_df <- female_res |> 
  filter(padj<0.01 & abs(log2FoldChange) > log2(1) ) |> 
  filter(!is.na(Symbol))
```

---

# 3. Enrichment Analysis

Genes in enriched pathways make strong PCR validation targets because they support biological interpretation.

```{r}
# enrich significant degs
enrich <- compareCluster(
  Symbol~direction,                # formula for clustering
  data = sig_df,                   # deg data frame
  fun = "enrichGO",                # enrichment function
  OrgDb = 'org.Hs.eg.db',          # organism your samples belong to
  keyType = "SYMBOL",              # type of gene name
  ont = "ALL",                     # enrichment ontology
  pvalueCutoff = 0.1,              # p-value cutoff
  pAdjustMethod = "fdr",           # p-value adjustment method
  universe = female_res$Symbol,    # what other genes did you test?
  qvalueCutoff = 0.1,              # q-value cutoff
  minGSSize = 2,                   # minimum number of genes per term
  maxGSSize = 1000                 # maximum number of genes per term
)
```

---

# 4. PCR Validation Target Selection Workflow 

Selecting qPCR targets is not simply choosing the most significant DEGs. It is a multi‑step evaluation process balancing statistics, biology, and assay design.

---

## 4.1 Step 1 - Start With Significant, Robust DEGs

**Criteria:**

- padj < 0.01  
- |log2FC| > 1.5  
- baseMean > 20  

PCR requires **consistent, strong, highly-expressed** differences.

```{r}
pcr_candidates <- sig_df |>
  filter(
    padj < 0.01,
    abs(log2FoldChange) > log2(1.5),
    baseMean > 20
  )
```

---

## 4.2 Step 2 - Prioritize Genes in Enriched Pathways

Genes belonging to enriched biological processes or pathways are more likely to reflect true biology and provide meaningful validation.

```{r}
enriched_genes <- enrich@compareClusterResult |>
               separate_rows(geneID,sep = "/") |>
               group_by(geneID) |> 
               reframe(pathway_hit=n()) |> 
               arrange(desc(pathway_hit))

pcr_candidates <- pcr_candidates |>
  inner_join(enriched_genes,
             by=c("Symbol"="geneID"))
```


---

## 4.3 Step 3 - Filtering by Gene Metadata

Now we will pull additional gene meta data from bioMart to understand:

- Is the gene protein coding?

- Does the gene have more than 2 exons?

- Is the minimum transcript length over 100 bp?


```{r}
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_meta <- getBM(
  attributes = c(
    "hgnc_symbol",
    "gene_biotype",
    "ensembl_exon_id",
    "transcript_length",
    "gene_biotype"
  ),
  filters = "hgnc_symbol",
  values = unique(ranked$Symbol),
  mart = mart
) |> 
  group_by(hgnc_symbol) |> 
  reframe(count_exon=length(unique(ensembl_exon_id)),
          protein_coding="protein_coding" %in% gene_biotype,
          transcript_length=min(transcript_length)) |> 
  distinct() |> 
  filter(count_exon>2,
         protein_coding==TRUE,
         transcript_length>100)
```

Now we will filter by these metrics.

```{r}

pcr_candidates <- pcr_candidates |> 
  filter(Symbol %in% gene_meta$hgnc_symbol)
```

## 4.4 Step 4 - Rank Genes by Biological and Statistical Strength

We sort the genes so the most qPCR‑ready candidates rise to the top.

```{r}
ranked <- pcr_candidates |>
  arrange(
    desc(abs(log2FoldChange)),  # then effect size
    desc(pathway_hit),          # biological relevance first
    desc(baseMean),             # then expression level
    desc(-log10(padj))          # then the p-value
  )
```

Now let's take a look at the top 10 upregulated genes!

```{r}
ranked |> 
  filter(direction == "Up") |> 
  dplyr::select(Symbol,log2FoldChange,baseMean,pathway_hit,padj) |> 
  slice_head(n=10)
```

```
Symbol log2FoldChange baseMean pathway_hit padj
<chr>  <dbl>          <dbl>    <int>       <dbl>
HSPA1A	2.455791	35096.8118	23	2.064107e-04
HSPA1B	2.443789	33144.9815	18	1.303694e-04
HSPA1L	1.887385	758.4465	1	3.737354e-05
ATF3	1.837993	632.1215	14	6.419790e-08
DNAJB1	1.830695	9915.0189	6	3.032908e-04
CXCR4	1.717091	516.1339	14	2.460203e-04
FOXJ1	1.643526	391.6230	19	3.994640e-03
SERPINH1	1.642161	1144.9492	1	4.483712e-04
RGS1	1.616057	905.4349	1	6.102832e-03
BCL6B	1.548744	207.3333	5	9.924136e-05
```



And the top 10 downregulated genes!

```{r}
ranked |> 
  filter(direction == "Down") |> 
  dplyr::select(Symbol,log2FoldChange,baseMean,pathway_hit,padj) |> 
  slice_head(n=10)
```


```
Symbol log2FoldChange baseMean pathway_hit padj
<chr>  <dbl>          <dbl>    <int>       <dbl>
PCSK1	-2.205614	258.76353	21	5.346873e-04
PCDH8	-2.163541	190.43384	19	4.758336e-05
ADCYAP1	-2.158613	158.49512	18	8.329149e-05
VGF	-2.104217	562.34778	28	6.733200e-06
BDNF	-2.073503	156.44470	34	1.000697e-06
CARTPT	-2.065815	129.88813	26	2.415623e-04
CHGB	-2.035318	1504.91848	1	4.556266e-04
GABRA4	-1.963389	295.29849	67	5.375256e-04
PPEF1	-1.935684	64.45886	1	1.592429e-05
OLFM3	-1.929110	308.58967	3	1.136850e-03
```

---


## 4.5 Step 5 - Check for Outlier‑Driven Effects 

RNA‑seq log2FC values can be inflated by outlier samples. Before selecting a PCR target, visualize gene counts:

```{r}
#| fig-width: 5
#| fig-height: 5

se[rowData(se)$Symbol %in% "PCSK1",] |> 
  tidybulk() |> 
  ggplot(aes(
    x=condition,
    y=vst,
    fill=condition,
    color=condition
  ))+
  geom_violin(alpha=.7)+
  geom_jitter(alpha=.3)+
  geom_pwc()+
  stat_summary(fun = "median",size=1)+
  theme_pubr(legend = "right")+
  labs(
    x="",
    y="Vst Exp.",
    color="",
    fill=""
    
  )
```

![](img/bad_pcr.png)

You want:

- clean separation between conditions  
- no single sample driving the difference  
- replicate consistency  

Let's try the another candidate!

```{r}
#| fig-width: 5
#| fig-height: 5

se[rowData(se)$Symbol %in% "VGF",] |> 
  tidybulk() |> 
  ggplot(aes(
    x=condition,
    y=vst,
    fill=condition,
    color=condition
  ))+
  geom_violin(alpha=.7)+
  geom_jitter(alpha=.3)+
  geom_pwc()+
  stat_summary(fun = "median",size=1)+
  theme_pubr(legend = "right")+
  labs(
    x="",
    y="Vst Exp.",
    color="",
    fill=""
    
  )
```

![](img/good_pcr.png)

This one has much cleaner separation between control and AD!

---



# 5. Summary

Here is our final checklist for selecting PCR candidats from RNA-seq results.

**Final checklist:**

- padj < 0.05  
- |log2FC| > log2(1.5)  
- baseMean > 20  
- in enriched pathways  
- protein coding
- adequate exons and transcript length
- not outlier‑driven effects  

# References

1. [Practical: Organizing data with SummarizedExperiment](https://bioconductor.org/help/course-materials/2019/BSS2019/04_Practical_CoreApproachesInBioconductor.html)
2. [tidybulk: an R tidy framework for modular transcriptomic data analysis](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02233-7)
3. [Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
4. [clusterProfiler: an R Package for Comparing Biological Themes Among Gene Clusters](https://pmc.ncbi.nlm.nih.gov/articles/PMC3339379/)
5. [Genome wide annotation for Human](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)
