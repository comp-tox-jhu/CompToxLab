# Navigating Files/Directories

File systems follow a heirarchy. Think of the files in your computer, you have folders which can contain folders or files. Now that have our tutorial data we will enter or *change* into our `unix_tutorial` folder or *directory* with the `cd` command (`cd` being short for change directory):

```bash
cd unix_tutorial
```

Now we have entered a directory! But how do we understand where we are? We can do this with the `pwd` command (print working directory):

```bash
pwd
```

!!! info "output"

    ```bash
    ~/Documents/unix_tutorial
    ```
    
This output will be a little different for every user based on your computer's file structure. However, you will note here that we are within the root folder (`~`) which contains all other folders, the `Documents` folder, and finally the `unix_tutorial` folder. This string of folders showing where we are is called our **file path** or just **path** for short. 

Now this `unix_tutorial` folder has other folders inside it. We can _list_ these folders with the `ls` command:

```bash
ls
```


!!! info "output"

    ```bash
    data_folder	results_folder	scripts_folder
    ```
    
Here we see that we have three folders: `data_folder`,	`results_folder`, and	`scripts_folder`. If we want more information on our files we can add *options* to our command. Let's do this by adding the long option to our `ls` command:

```bash
ls -l
```

!!! info "output"

    ```bash
    total 0
    drwxr-xr-x  4 username  groupname  128 May 13 09:33 data_folder
    drwxr-xr-x  2 username  groupname   64 May 13 09:26 results_folder
    drwxr-xr-x  2 username  groupname   64 May 13 09:26 scripts_folder
    ```
The long format will usually have 8 columns:

- permissions: split by us

The Unix file system is organized in a hierarchical structure, starting from the root directory ("/") and branching out into various directories and subdirectories. Here are some essential commands for navigating the file system:
pwd: Print the current working directory.
ls: List files and directories in the current directory.
-l: Long format, providing detailed information about files.
-a: Include hidden files and directories.
cd: Change directory.
cd [directory]: Move to the specified directory.
cd ..: Move up one directory level.
cd ~: Move to the home directory.
