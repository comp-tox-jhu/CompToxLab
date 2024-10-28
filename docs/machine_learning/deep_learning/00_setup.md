## Software

First things first you'll need an integrated development environment or IDE to start coding. We recomend VS code for coding in python and really anything other than R. Use the following link to install VS code:

[VS Code Installation](https://code.visualstudio.com/download){:target="_blank" .md-button .md-button--primary}

Additionally, you'll need Anaconda for managing packages, use the following link to install Anaconda:

 [Anaconda Installation](https://docs.anaconda.com/anaconda/install/){:target="_blank" .md-button .md-button--primary}

## Extensions

Once both are installed you'll need to install some extensions in VS code. We will be working with Jupyter notebooks so we need to install the Jupyter extension:

1. Click `Extensions`
2. Type `Jupyter` in the search bar and download the Jupyter extension

## Project Setup

Now let's set up a space for us to work in. We will start by creating a virtual environment:

1. Click `Open Folder`
2. Create a folder and name it `deep_learning`
3. Go inside the deep learning folder `cd deep_learning`
4. Create 3 folders inside `deep_learning`:
5. `mkdir data`
6. `mkdir scripts`
7. `mkdir results`

To download the data we will use:

1. Go into the `data` folder: `cd data`
2. `curl --output gbm_cptac_2021.tar.gz https://cbioportal-datahub.s3.amazonaws.com/gbm_cptac_2021.tar.gz`

Now to set up the virtual environment:

1. Go to `View` > `Terminal`
2. In terminal enter: `conda create -n dl python=3.12 pytorch torchvision plotly`
3. Then enter `y`
4. Then at the top of VS code select "Select Kernel"
5. Select "dl"

[Now Let's Start Coding!](./01_basics.md){:target="_blank" .md-button .md-button--primary}
