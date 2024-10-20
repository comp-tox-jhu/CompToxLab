## Software

First things first you'll need an integrated development environment or IDE to start coding. We recomend VS code for coding in python and really anything other than R. Use the following link to install VS code:

[VS Code Installation](https://code.visualstudio.com/download){ .md-button .md-button--primary :target="_blank"}

Additionally, you'll need Anaconda for managing packages, use the following link to install Anaconda:

 [Anaconda Installation](https://docs.anaconda.com/anaconda/install/){ .md-button .md-button--primary :target="_blank"}

## Extensions

Once both are installed you'll need to install some extensions in VS code. We will be working with Jupyter notebooks so we need to install the Jupyter extension:

1. Click `Extensions`
2. Type `Jupyter` in the search bar and download the Jupyter extension

## Project Setup

Now let's set up a space for us to work in. We will start by creating a virtual environment:

1. Click `Open Folder`
2. Create a folder and name it `deep_learning`

Now to set up the virtual environment:

1. Go to `View` > `Terminal`
2. In terminal enter: `conda create -n dl python=3.12 pytorch torchvision plotly`
3. Then enter `y`
4. Then at the top of VS code select "Select Kernel"
5. Select "dl"

[Now Let's Start Coding!](/machine_learning/deep_learning/01_basics.md){ .md-button .md-button--primary}

Now we will use the command palette in VS code!

1. Use `cmd+shift+p`
2. Type and choose "Create: New Jupyter Notebook"
3. In the notebook select `Select Kernel`
4. Then select `Python Environments..`
5. Then "Create Python Environment"
6. Then "Conda" and select "Python 3.12"

Now in the notebook, create a code chunk and run:

```
!conda install pytorch torchvision
!conda install plotly
```


## Project Setup

1. Go to `View` > `Terminal`
2. In terminal enter: `python3 -m venv dl_venv`
3. Then enter: `source dl_venv/bin/activate`
4. Then enter: `pip3 install ipykernel`
5. Then enter: `python3 -m ipykernel install --user --name=dl_kernel`

## Notebook Setup

Now we will use the command palette in VS code!

1. Use `cmd+shift+p`
2. Type and choose "Create: New Jupyter Notebook"
3. Use `cmd+shift+p` again
4. Type and choose "Notebook: Select Notebook Kernel"
5. Select "Jupyter Kernal.."
6. Choose `dl_kernel`
7. In the notebook enter the following in a code cell: `pip install python=3.12 pytorch torchvision plotly`
