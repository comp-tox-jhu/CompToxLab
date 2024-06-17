## Introduction To RStudio

**Author:** Jason Laird, PhD Student, M.Sc.

!!! info "Learning Objectives"

    - [Introduction To RStudio](programming_languages_tools/intro_to_r/01_introduction.md)
    - [Data Types/Variables/Vectors](programming_languages_tools/intro_to_r/02_data_types_variables_vectors.md)
    - [Data Structures](programming_languages_tools/intro_to_r/03_data-structures.md)
    - [Functions/Flow](programming_languages_tools/intro_to_r/04_functions-flow.md)
    - [Inspecting/Manipulating Data](programming_languages_tools/intro_to_r/05_inspecting-manipulating-data.md)
    - [Data Vizualization](programming_languages_tools/intro_to_r/06_visualization.md)


RStudio is what is known as an Integrated Development Environment or IDE. Here you can write scripts, run R code, use R packages, view plots, and manage projects. To download RStudio on your own computer, use the link below:

!!! abstract "[Download RStudio](https://posit.co/download/rstudio-desktop/)"

Once downloaded, open RStudio. You'll see the following three panels:

- **The Interactive R console/Terminal (left)**
- **Environment/History/Connections (upper right)**
- **Files/Plots/Packages/Help/Viewer (lower right)**

!!! info "RStudio Layout"

    ![](images/rstudio1.png)

## Project Management

Before we dive into R it is worth taking a moment to talk about project management. Often times data analysis is incremental and files build up over time resulting in messy directories:

!!! info "Example of a Messy Directory"

    ![](images/messy.png)

Sifting through a non-organized file system can make it difficult to find files, share data/scripts, and identify different versions of scripts. To remedy this, It is reccomended to work within an R Project. Before we make this project, we should make sure you are in your home directory. To do this click on the three dots in the files tab:

!!! info "Navigating Folders"

    ![](images/three_dots.png)

Then enter in a ~ symbol to go home!

!!! info "Getting Home"

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
download.file("https://raw.githubusercontent.com/BioNomad/omicsTrain/main/docs/programming_languages_tools/intro_to_r/data/metadata.csv",destfile = "data/metadata.csv")
download.file("https://raw.githubusercontent.com/BioNomad/omicsTrain/main/docs/programming_languages_tools/intro_to_r/data/metadata.tsv",destfile = "data/metadata.tsv")
download.file("https://raw.githubusercontent.com/BioNomad/omicsTrain/main/docs/programming_languages_tools/intro_to_r/data/test.xlsx",destfile = "data/test.xlsx")
```

## Data Principles

- Treat data as read-only
- Store raw data separately from cleaned data if you do need to manipulate it
- Ensure scripts to clean data are kept in a separate `scripts` folder
- Treat reproducible results as disposable

!!! tip
    Result files are good candidate files to cut if you are getting low on storage.

## New R script

Now we will create an R script. R commands can be entered into the console, but saving these commands in a script will allow us to rerun these commands at a later date. To create an R script we will need to either:

- Go to `File > New File > R script`
- Click the `New File` icon and select R script

!!! info "Creating a New R Script"

    ![](images/newFile.png)

## Running R Code

When running R code you have a few options:

  Running One Line/Chunk:
  
  - Put your cursor at the beginning of the line of code and hit `Ctrl + Enter` on Windows or  &#8984; + `Enter` on MacOSX.
    
  - Highlight the line/chunk of code and hit `Ctrl + Enter` or &#8984; + `Enter`.
    
  Running The Entire Script:
  
  - Clicking `Source` at the top of the script window.
  
## References

1. [R for Reproducible Scientific Analysis](https://swcarpentry.github.io/r-novice-gapminder/)
2. [Base R Cheat Sheet](https://iqss.github.io/dss-workshops/R/Rintro/base-r-cheat-sheet.pdf)
