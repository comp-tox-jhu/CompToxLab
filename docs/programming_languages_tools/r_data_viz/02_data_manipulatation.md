## Loading Libraries

To manipulate our data and ultimately plot it, we often use libraries of functions to reduce the amount of coding we need to do. For the purposes of this tutorial we will load the following:

```R
library(tidyverse)          # data manipulation/plotting
library(janitor)            # data tabulation
library(reshape2)           # reshape to plug into ggplot
library(ggfortify)          # pca plot
library(ggplotify)          # convert to ggplot
library(png)                # use png images
library(grid)               # modify png images
library(pheatmap)           # create heatmaps
library(patchwork)          # combine multiple plots
library(DESeq2)             # run differential expression
source("./de_wrappers.R")   # source custom differential expression function
```

## Loading Data

To load our data we will be using functions from the `readr` packages instead of the base functions provided through base R. The `readr` collection of functions are advantageous as they are typically more user friendly and are faster at importing data. Let's load our sample and meta data:

```R
# load sample data
meta <- read_delim(file="../data/meta_data.tsv",
             delim="\t")

# load mRNA data
counts <- read_delim(file="../data/expression_data.tsv",
             delim="\t")
```

## Inspecting Data

We can examine the structure of our data with the `str` function:

```{r,message=F,warning=F}
str(meta)
```

!!! info "output"

    ```
        spc_tbl_ [113 × 6] (S3: spec_tbl_df/tbl_df/tbl/data.frame)
     $ patient_id     : chr [1:113] "GSM3561843" "GSM3561844" "GSM3561845" "GSM3561846" ...
     $ braak.score:ch1: chr [1:113] "III" "IV" "II" "V" ...
     $ cell type:ch1  : chr [1:113] "myeloid" "endothelial" "neuron" "neuron" ...
     $ diagnosis:ch1  : chr [1:113] "Control" "Control" "Control" "AD" ...
     $ expired_age:ch1: chr [1:113] ">90" "88" "79" "84" ...
     $ Sex:ch1        : chr [1:113] "Female" "Female" "Male" "Female" ...
     - attr(*, "spec")=
      .. cols(
      ..   patient_id = col_character(),
      ..   `braak.score:ch1` = col_character(),
      ..   `cell type:ch1` = col_character(),
      ..   `diagnosis:ch1` = col_character(),
      ..   `expired_age:ch1` = col_character(),
      ..   `Sex:ch1` = col_character()
      .. )
     - attr(*, "problems")=<externalptr> 
    ```
    
Here we see the overall structure of our meta data. We have a data frame with 113 rows, 6 columns, the first few observations per column, and the data type of each column. Now if this is a bit much to look at we can also just look at the top 6 rows of our data:

```R
head(counts[,1:5])
```

!!! info "output"

    ```
        # A tibble: 6 × 5
      SYMBOL   GSM3561843 GSM3561844 GSM3561845 GSM3561846
      <chr>         <dbl>      <dbl>      <dbl>      <dbl>
    1 AKT3           1391       1918        472        773
    2 NR2E3             2         50          9         20
    3 NAALADL1         10         18          2          1
    4 SIGLEC14        512          0          1          1
    5 MIR708            4          0          0          0
    6 NAV2-AS6         24         44         47         35
    ```
    
## Convert A Column To Rownames

Sometimes when we import our data, we may not get it in the form we need it in. In that case we will need to do some data munging (AKA data cleaning). A common operation is take a column of unique identifiers and turn that column into the rownames for the data frame:

!!! info "Column to Rownames"

    ![](../images/column_to_rownames.png)
    
We will need to make the gene names our row names we can do this with the `column_to_rownames` funtion:

```R
# convert the SYMBOL column to rownames
counts <- counts %>%
  column_to_rownames("SYMBOL")

# now let's inspect the new counts data frame
head(counts[,1:5])
```


!!! info "output"

    ```
                 GSM3561843 GSM3561844 GSM3561845 GSM3561846 GSM3561847
    AKT3           1391       1918        472        773        140
    NR2E3             2         50          9         20          9
    NAALADL1         10         18          2          1          0
    SIGLEC14        512          0          1          1         31
    MIR708            4          0          0          0          0
    NAV2-AS6         24         44         47         35          0
    ```
    
Let's do the same for patient id in our meta data. This way we can better map our counts data frame to our meta data data frame:

```R
# make the patient_id column the rownames
meta <- meta %>%
  column_to_rownames("patient_id")

# now let's inspect the new counts data frame
head(meta[,1:5])
```

!!! info "output"

    ```
                   braak.score:ch1 cell type:ch1 diagnosis:ch1 expired_age:ch1 Sex:ch1
    GSM3561843             III       myeloid       Control             >90  Female
    GSM3561844              IV   endothelial       Control              88  Female
    GSM3561845              II        neuron       Control              79    Male
    GSM3561846               V        neuron            AD              84  Female
    GSM3561847              II       myeloid       Control              73    <NA>
    GSM3561848               V       myeloid            AD              84  Female
    ```
    
## Changing Column Names

Data won't always come in with user friendly column names. We can change column names in a data frame with the `rename`/`rename_with` functions:


!!! info  "Renaming Column Names"

    ![](../images/rename.png)
    
In our meta data frame all our column names have extra characters on the end. We can change all these column names with the `rename_with` function:

```{r,message=F,warning=F}
# rename columns that have the pattern :ch1
meta <- meta %>%
  rename_with(~ gsub(":ch1","",.))

# we can check our data frame to see our changes
head(meta)
```

!!! info "output"

    ```
                   braak.score   cell type diagnosis expired_age    Sex
    GSM3561843         III     myeloid   Control         >90 Female
    GSM3561844          IV endothelial   Control          88 Female
    GSM3561845          II      neuron   Control          79   Male
    GSM3561846           V      neuron        AD          84 Female
    GSM3561847          II     myeloid   Control          73   <NA>
    GSM3561848           V     myeloid        AD          84 Female
    ```

If we were interested in changing the name of one column we could use the `rename` function to do so:

```{r,message=F,warning=F}
# rename the cell type column
meta <- meta %>%
  dplyr::rename("cell_type"="cell type")

# let's check the names of our data frame now!
head(meta)
```


!!! info "output"

    ```
                   braak.score   cell_type diagnosis expired_age    Sex
    GSM3561843         III     myeloid   Control         >90 Female
    GSM3561844          IV endothelial   Control          88 Female
    GSM3561845          II      neuron   Control          79   Male
    GSM3561846           V      neuron        AD          84 Female
    GSM3561847          II     myeloid   Control          73   <NA>
    GSM3561848           V     myeloid        AD          84 Female
    ```
    
## Filtering Rows

Not all rows in your data will be worth keeping we can remove rows with the `filter` function:


!!! info "Filtering Data Frames"

    ![](../images/filter.png)
    
    
We will be assessing differentially expressed genes between Alzheimer's disease patients and control patients. However, we will be filtering out any non-neuronal tissue samples:

```R
# only keep neuronal samples
meta_filt <- meta %>%
  filter(cell_type=="neuron")

# we can check how many rows were filtered out by comparing the filtered
# data frame with the previous data frame
nrow(meta)
nrow(meta_filt)
```


!!! info "output"

    ```
    [1] 113
    [1] 42
    ```
    
## Selecting Columns

we have just removed 71 rows (so 71 samples) from our meta data. To run DESeq2 we will need to remove those samples from our counts data frame as well. We can remove columns using the `select` function:

!!! info "Selecting Columns"

    ![](../images/select.png)
    
    
```{r,message=F,warning=F}
# select only columns whose names are present in the rownames of the meta data
counts_filt <- counts %>%
  select(rownames(meta_filt))

# let's check how many columns are present in the filtered and unfiltered data frames
ncol(counts)
ncol(counts_filt)
```

!!! info "output"

    ```
    [1] 113
    [1] 42
    ```
    
## Add A Column

Right now our condition variable is `Control` v. `AD`. We may want to add a column that spells out that `AD` represents Alzheimer's Disease. We can add columns with the `mutate` function:

!!! info "Adding Columns To A Data Frame"

    ![](../images/mutate.png)
    
```R
# let's create a column for Alzheimer's Disease Status using the ifelse function to insert a value if a test is TRUE
meta_filt <- meta_filt %>%
  mutate(diagnosis = gsub("AD","Alzheimer's Disease",diagnosis))
# let's inspect this new column!
meta_filt$diagnosis[1:10]
```

!!! info "output"

    ```
    [1] "Control"             "Alzheimer's Disease" "Control"             "Alzheimer's Disease"
    [5] "Alzheimer's Disease" "Control"             "Control"             "Alzheimer's Disease"
    [9] "Control"             "Alzheimer's Disease"
    ```
    
## Changing A Column

We can also use mutate to update an existing column as well let's update the disease status column to be a factor setting this as a factor will help with DESeq2 later on:

```R
# factor the diagnosis column to set the `Control` value as the reference
meta_filt <- meta_filt %>%
  mutate(diagnosis=factor(diagnosis,
                               levels = c(
                                 "Control",
                                 "Alzheimer's Disease"
                               )))

# let's inspect our updated column!
meta_filt$diagnosis[1:10]
```

!!! info "output"

    ```
    [1] Control             Alzheimer's Disease Control             Alzheimer's Disease
    [5] Alzheimer's Disease Control             Control             Alzheimer's Disease
    [9] Control             Alzheimer's Disease
    Levels: Control Alzheimer's Disease
    ```
    
## Merging Data Frames

Sometimes we need to take two data frames and merge them into one data frame. We can do this with the `join` functions. Here we will do an inner join (meaning we merge two data frames an only keep rows that match by some column):

!!! info "Merging Data Frames"
    
    ![](../images/inner_join.png)
    
```R
# flip the counts data so patients are the columns
# make the rownames a patient id column
counts_flip <- counts_filt %>%
  t() %>%
  data.frame() %>%
  rownames_to_column("patient_id")

# make the rownames a patient id column
meta_patient_to_col <- meta_filt %>%
  rownames_to_column("patient_id")

# inner join the meta data and the counts data
merged <- meta_patient_to_col %>%
  inner_join(
    .,
    counts_flip,
    by="patient_id"
  )%>%
  `rownames<-`(.[["patient_id"]])

nrow(merged)
ncol(merged)
```

!!! info "output"

    ```
    [1] 42
    [1] 3006
    ```
    
## Creating a Table For Publication

When getting reading for a publication we may want to include a summary of our data. The `janitor` package contains many helpful functions for that purpose! Let's try to create a table of sex by disease status:

```{r,message=F,warning=F}
# often times you may want to tabulate your data to prepare a 
# table 1 for your paper. An easy way of doing this is to
# use the tabyl/adorn functions from the janitor package
table_1 <- merged %>%                     # data frame to work with
  tabyl(Sex, diagnosis) %>%               # tabulate sex by diagnosis
  adorn_totals("row") %>%                 # adorn totals by row
  adorn_percentages("all") %>%            # adorn percentages to all cells
  adorn_pct_formatting(digits = 1) %>%    # round the percentage to the first digit
  adorn_ns                                # keep both count and percentage

# let's take a look at our table 1!
table_1
```

!!! info "output"

    ```
            Sex    Control Alzheimer's Disease
     Female 14.3%  (6)          16.7%  (7)
       Male 28.6% (12)          19.0%  (8)
       <NA>  7.1%  (3)          14.3%  (6)
      Total 50.0% (21)          50.0% (21)
    ```
