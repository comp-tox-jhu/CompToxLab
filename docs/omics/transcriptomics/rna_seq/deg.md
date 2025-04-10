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

To load our data we will be pulling from the `readr` package from the `tidyverse`:

So in that code chunk under 

