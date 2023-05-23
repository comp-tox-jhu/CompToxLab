## Introduction To RStudio

RStudio is what is known as an Integrated Development Environment or IDE. Here you can write scripts, run R code, use R packages, view plots, and manage projects. This pane is broken up into three panels:

- **The Interactive R console/Terminal (left)**
- **Environment/History/Connections (upper right)**
- **Files/Plots/Packages/Help/Viewer (lower right)**

![](images/rstudio1.png)

## Project Management

Before we dive into R it is worth taking a moment to talk about project management. Often times data analysis is incremental and files build up over time resulting in messy directories:

![](images/messy.png)

Sifting through a non-organized file system can make it difficult to find files, share data/scripts, and identify different versions of scripts. To remedy this, It is reccomended to work within an R Project. Before we make this project, we should make sure you are in your home directory. To do this click on the three dots in the files tab:

![](images/three_dots.png)

Then enter in a ~ symbol to go home!

![](images/getting_home.png)


## R Project

For the following intro to R tutorial we will be using glioblastoma data from [cBioPortal](https://www.cbioportal.org/study/summary?id=gbm_cptac_2021). When working within R it is useful to set up an R project. R projects will set your working directory relative to the project directory. This can help ensure you are only working with files within this project space. To create a new project:

1. Go to `File` > `New Project`
2. `New Directory`
3. `New Project`
4. Create a name for your project (e.g. `intro_to_r`)
5. `Create Project`
     
When analyzing data it is useful to create a folder to house your raw data, scripts and results. We can do this by clicking the `New Folder` icon to create these folders:

1. Click `New Folder` > Enter `data` > Click OK
2. Click `New Folder` > Enter `scripts` > Click OK
3. Click `New Folder` > Enter `results` > Click OK
    
Now that we have our project set up we will need to download our data. In the `data` folder we will download our data and decompress it:

``` R
download.file(url = "https://github.com/BioNomad/omicsTrain/tree/main/docs/programming_languages_tools/intro_to_r/data/data.zip",destfile = "./data/initial_data.zip" )
unzip(tarfile = "./data/data",exdir = "./data/")
```

## Data Principles

- Treat data as read-only
- Store raw data separately from cleaned data if you do need to manipulate it
- Ensure scripts to clean data are kept in a separate `scripts` folder
- Treat reproducible results as disposable

!!! tip
    Result files are good candidate files to cut if you are getting low on storage.


## Getting Data

![](images/data_summary.png)

- Today we will be using a fake dataset assessing the taxa count on the mouse microbiome before and after antibiotic usage.
- To copy over this data we will use an R function called file.copy. 
- A function takes some input and delivers an output. 
- In this case we specify two inputs the location of our file and where we want to copy it to. 
- The function's output is copying over this file. So let's try it copy over using the following commands:



```{r data.copy,warning=F,message=F}
file.copy(from="/cluster/tufts/bio/tools/training/intro-to-r/data/meta.tsv", to="./data/")
file.copy(from="/cluster/tufts/bio/tools/training/intro-to-r/data/meta2.tsv", to="./data/")
```

So here you'll note we copied over the file metadata.tsv to the data folder. Let's copy over our script:

```{r script.copy,warning=F,message=F}
file.copy(from="/cluster/tufts/bio/tools/training/intro-to-r/scripts/intro-to-r.Rmd", to="./scripts")
```

Here we copy over our script intro-to-r to the scripts folder.



## Opening the Script

Now let's start by opening our script. Go to scripts and then double click on intro-to-r.Rmd!
