# Introduction to Biostatistics

Biostatistics attempts to use statiscal methods to solve biological problems. This involves data so for the purpose of our biostatistics tutorials we will need to do some setup.

## Setup 

For the following machine learning tutorials we will be using glioblastoma data from [cBioPortal](https://www.cbioportal.org/study/summary?id=gbm_cptac_2021). When working within R it is useful to set up an R project. R projects will set your working directory relative to the project directory. This can help ensure you are only working with files within this project space. To create a new project:
    
1. Go to `File` > `New Project`
2. `New Directory`
3. `New Project`
4. Create a name for your project (e.g. `biostatistics`)
5. `Create Project`
     
When analyzing data it is useful to create a folder to house your raw data, scripts and results. We can do this by clicking the `New Folder` icon to create these folders:

1. Click `New Folder` > Enter `data` > Click OK
2. Click `New Folder` > Enter `scripts` > Click OK
3. Click `New Folder` > Enter `results` > Click OK
    
Now that we have our project set up we will need to download our data. In the `data` folder we will download our data and decompress it:

``` R
download.file(url = "https://cbioportal-datahub.s3.amazonaws.com/gbm_cptac_2021.tar.gz",destfile = "./data/gbm_cptac_2021.tar.gz" )
untar(tarfile = "./data/gbm_cptac_2021.tar.gz",exdir = "./data/")
```
