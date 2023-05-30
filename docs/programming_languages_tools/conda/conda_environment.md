## Conda Environments

For reproducible research it is advisable to keep the software versions you use consistent. An easy way of ensuring this is by creating a Conda environment, where you can explicitly say which versions are necessary to run your  pipeline:

!!! info "[Conda Environments](https://nbisweden.github.io/workshop-scRNAseq/conda_instructions.html)"

    ![](images/conda_envs.png)
    
## Downloading Anaconda

To build a conda environment, you will need the Anaconda Python distribution. You can download this at the following link:

!!! info "[Download Anaconda](https://www.anaconda.com/download)"


## Creating a Conda Environment

Once you have downloaded Anaconda, go to your terminal (Terminal for Macs and [Mobaxterm](https://mobaxterm.mobatek.net/download.html) is recommended for PC). Now in the terminal create your environment with the desired version of python:

```sh
conda create -n yourenvname python=3.8
```

Once created you can activate your environment like so:


```sh
conda acvitate yourenvname
```

## Installing Tools

To install tools in your conda environment use the following syntax:

```sh
conda install yourpackage
```

Check out a list of available tools here:

!!! info "[Anaconda Packages](https://anaconda.org/)


Or if there is no available conda package you can also try installing these tools with pip:

```sh
pip install yourpackage
```


To see what's installed in your conda environment use:

```sh
conda list
```

## Deactivating An Environment

When you are finished using your environment, deactivate the environment with:

```sh
conda deactivate
```


## Installing A Conda Environment From a File

When you run an analysis using a conda environment, you can pass that environment onto others by exporting the tool versions to a yml file with the following command:

```sh
conda env export -n yourenvname -f yourenvname.yml --no-builds
```

Here we add no `--no-builds` because the environment you build is often specific to the machine you built it on. If you build the environment on Mac and want others to use it on a PC or linux platform it is advisable you export it without the "build" information your specific machine used.


## References

1. [Manage Environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
2. [Conda Instructions](https://nbisweden.github.io/workshop-scRNAseq/conda_instructions.html)
