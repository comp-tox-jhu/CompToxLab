# Managing Files/Directories

## Making Files/Directories

Just like your Folders or Finder app you can use the command line to manipulate your files/directories. Let's start by making a new directory with the `mkdir` command  (short for make directory):

```bash
mkdir new_dir
ls
```

!!! info "output"

    ```bash
    data_folder	new_dir		results_folder	scripts_folder
    ```
Great! We've made a new directory! Now how about files? We can make empty files with the `touch` command:

```bash
touch new_file
ls
```

!!! info "output"

    ```bash
    data_folder	new_dir		new_file	results_folder	scripts_folder
    ```
## Copying/Renaming Files/Directories

Sometimes we may want to copy a file to another folder, which we can do with the cp command (short for copy)!

```bash
cp new_file data_folder
ls data_folder
```

!!! info "output"

    ```bash
    accList.txt	meta.txt	new_file
    ```
We can see that we have successfully copied over our file! If we wanted to copy a directory we can use the `-r` option:

```bash
cp new_dir data_folder
ls data_folder
```

!!! info "output"

    ```bash
    accList.txt	meta.txt	new_file	new_dir
    ```
    
If we just wanted to move our file or directory without copying, we can use the `mv` commmand (short for move). However, the `mv` command can also be used to rename a file/folder as well! Let's use it to rename our file!

```bash
mv new_file brand_new_file
ls
```

!!! info "output"

    ```bash
    data_folder	new_dir		brand_new_file	results_folder	scripts_folder
    ```

## Removing Files/Folders

We can also delete/remove files/folders - However, do this with extreme caution as this is a permanent deletion! We will remove our empty files that we just created with the `rm` (short for remove) command!

```bash
rm brand_new_file
ls
```

!!! info "output"

    ```bash
    data_folder	new_dir	results_folder	scripts_folder
    ```
To remove our directory we will need the `-r` option:

```bash
rm -r new_dir
ls
```

!!! info "output"

    ```bash
    data_folder	results_folder	scripts_folder
    ```

## Viewing Files

To view the contents of a file, you can use the `cat` command. So let's enter our `data_folder` to view our data!

```bash
cd data_folder
cat accList.txt
```

!!! info "output"

    ```bash
    SRR1219879	
    SRR1219880
    ```
    
Sometimes a file can be thousands of lines and you may only want to view  a certain portion of it. We can view the begnnning of a file with `more` and the end of a file with the `less` command. To exit out of this viewing mode, just hit `q`! To get a quick idea of what your file contains you can also use the `head` command to grab the first 6 lines!

## Editing File Content

To edit the content of a file we need to use a text editor like nano or vim to edit on the command line. The most user-friendly editor is typically nano. So let's edit `new_file` which we have copied over to the `data_folder`:

```bash
nano new_file
```

!!! info "output"

    ```bash
    UW PICO 5.09                         File: new_file                           



    ^G Get Help  ^O WriteOut  ^R Read File ^Y Prev Pg   ^K Cut Text  ^C Cur Pos   
    ^X Exit      ^J Justify   ^W Where is  ^V Next Pg   ^U UnCut Text^T To Spell  
    ```
    
Here we will write the words `Hello World` then hit `Control` and then `X` where you will be prompted to save:

!!! info "output"

    ```bash
    UW PICO 5.09                         File: new_file 
    Hello World



    Save modified buffer (ANSWERING "No" WILL DESTROY CHANGES) ?                    
                 Y Yes                                                            
    ^C Cancel    N No  
    ```
    
Hit `Y` to save then `Enter` to go back to the command line! Let's check the contents of our file:

```bash
cat new_file
```


!!! info "output"

    ```bash
    Hello World
    ```
