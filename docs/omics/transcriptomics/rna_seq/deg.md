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

```R

```

