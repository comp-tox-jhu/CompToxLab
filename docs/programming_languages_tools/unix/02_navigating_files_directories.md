# Navigating Files/Directories

## Changing and Checking Directories

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

## Listing Directory Contents

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
- permissions of the file
    
    !!! info "permissions"
    
        permissions are split into a line of 10 characters. The first character will be a `d` if it is a directory and a `-` if it is a file. Then the characters are read in groups of three - the first three characters are the user's permissions, the next three are the group permissions and the last three are other user permissions. `r` means someone can read the file/folder, `x` means someone can execute (like as in executing a script), and `w` means someone can write to the file/folder. 
    
- number of files (folders count as 2, files count as 1, so for example data_folder has 2 files in it)
- the username
- the group name
- size of the file/folder
- day of modification
- time of modification
- file/folder name

## Navigation Shortcuts

Now let's practice navigating directories by entering our `data_folder`:

```bash
cd data_folder
pwd
```

!!! info "output"

    ```bash
    ~/Documents/unix_tutorial/data_folder
    ```
   

To move around we have a few shortcuts. If we want to move to our root folder (where all other folders live) we can use the `cd ~` command. To move one folder up we can use the `cd ..` command:

```bash
cd ..
pwd
```

!!! info "output"

    ```bash
    ~/Documents/unix_tutorial
    ```
