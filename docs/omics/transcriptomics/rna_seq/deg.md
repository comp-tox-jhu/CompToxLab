# Differential Expression Tutorial

## Setup

First thing is first we need a folder to work with! Open RStudio and create a folder by clicking the folder button:

![](./img/folder_button.png)

Now call this folder `rna_seq_tutorial`! Click inside this folder and create a few more folders:

1. `data`
2. `scripts`
3. `results`

Great now in console we will download our count data and meta data!

```R
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE125583&format=file&file=GSE125583_raw_counts_GRCh38.p13_NCBI.tsv.gz",destfile = "data/count_data.tsv.gz")
```
Now we will unzip this file using the `gunzip` function from `R.utils`. If you don't have `R.utils` go ahead and install using `install.packages("R.utils")`.

```R
R.utils::gunzip("data/count_data.tsv.gz")
```
Great, now let's get our meta data:

```R
download.file("https://raw.githubusercontent.com/comp-tox-jhu/CompToxLab/refs/heads/main/docs/omics/transcriptomics/rna_seq/data/meta.csv",destfile = "data/meta.csv")
```

Now our annotation data:

```R
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz",destfile = "data/annot.tsv.gz")
```

And we will uncompress it:

```R
R.utils::gunzip("data/annot.tsv.gz")
```

Wonderful, now that we have both we can get to work! Let's create a script by going to create a new script by going to `File` > `New File` > `R Notebook`:

![](img/new_script.png)

When we create a new R markdown there is some helpful text to get you going, let's delete that until all you see is:

```
---
title: "R Notebook"
output: html_notebook
---
```

Now let's change that title to "RNA-seq Tutorial". To start we need to make some headers to guide what we will do. We make headers by using hashtags `#` and the more hashtags the smaller the header. Let's make the following:

```
## Loading/Cleaning Data
## Pre-Processing
## Principal Component Analysis 
## Differential Expression
## Volcano Plots
## Functional Enrichment
```

Great! now under `## Setup` we will make a code chunk and we can do this using the short cut `Cntrl+Alt+I` or `Cmd+Alt+I` on a Mac. in that chunk we will load our libraries:

```R
library(tidyverse)         # data manipulation/plotting
library(DESeq2)            # differential expression
library(EnhancedVolcano)   # creating volcano plots
library(ggpubr)            # publication ready plotting
library(clusterProfiler)   # enrichment analysis 
```
Wonderful with this we can start bringing in data to play with!

## Loading/Cleaning Data

To load our data we will be pulling from the `readr` package from the `tidyverse`. So make another code chunk under `## Setup`:

```R
counts <- readr::read_tsv("../data/count_data.tsv")
meta <- readr::read_csv("../data/meta.csv")
annot <- readr::read_tsv("../data/annot.tsv")
```

If you'll notice, we don't have gene names in our first column of our counts data frame. So we need to map it to our annotation file:

```R
# remove the first column
counts_clean <- counts |> 
  column_to_rownames("GeneID") 
```

So what we did here is we took the counts data frame, and made the `GeneID` column the row names. Now let's clean up our meta data:

```R
# remove the first column
meta_clean <- meta |> 
  column_to_rownames("geo_accession")
```
Now we have our GEO accession number as our rownames, this ensures that the row names of our meta data match the column names of our counts data! But we still need to make sure everything is in the same order:

```R
# reoder meta data to match the order of the counts data
meta_clean <- meta_clean[match(colnames(counts_clean), rownames(meta_clean)), ]

# ensure the order is the same
all(rownames(meta_clean) == colnames(counts_clean))
```

Great! if this says `TRUE` that means the order of our samples is the same for both the count and meta data! Now let's make a column in our meta data to define our groups:

```R
# add in a variable for condition
meta_clean <- meta_clean |> 
  mutate(condition = case_when(
    grepl("AD", title) ~ "AD",
    grepl("CON ", title) ~ "Control"
  )) |> 
  mutate(condition = factor(condition,levels=c("Control", "AD")))
```

Here we say we are creating a new column `condition`, and when another column, `title` has the pattern `AD` we put `AD` as the value, and if the pattern `CON` is found, we put `Control`.

## Differential Expression

Ok onto the fun part, let's do differential expression! We will be using DESeq2 to perform differential expression and find genes that are differentially expressed between Alzheimer's disease and controls.

```R
# create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=counts_clean, 
                              colData=meta_clean, 
                              design=~condition)

# run DESeq2 on data
dds <- DESeq(dds)

# convert results to a data frame
res <- data.frame(results(dds))
```

Nice what we have done is taken count data and a meta data file, defined our formula (~ condition), performed differential expression, and collected our results as a data.frame! However, if you click on the results, you'll notice that `GeneID` isn't super informative! They are just numbers. To get gene names we will need to map it back to our annotation file:

```R
# add gene names to the results
res_mapped <- res |>
  mutate(GeneID=rownames(res)) |>
  mutate(GeneID=as.character(GeneID)) |>
  inner_join(annot |> 
               mutate(GeneID=as.character(GeneID)),
             by="GeneID")
```

Here we took our results made a column named `GeneID` converted it to a character for mapping,  joined this data frame with the annotation data frame where we also converted the `GeneID` to a character value, and mapped on the `GeneID` column. 

## Principal Component Analysis

Before we check out what these genes do, let's examine how our samples group using principal component analysis. Luckily for us, DESeq2 has wrapper functions to take care of this!



