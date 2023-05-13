# Introduction To Shell/Unix

Many software don't have a graphical user interface (or GUI) that users can use to run a pipeline. Instead many of these tools are only available through the command line interface. To use the command line interface you will need to know shell commands - which is exactly what we will cover here!

!!! abstract "Learning Objectives"

    1. [Navigating Files & Directories](./02_navigating_files_directories.md)
    2. [Managing Files & Directories](./03_managing_files_directories.md)
    3. [Redirection and Pipes](./04_redirection_pipes.md)
    4. [Regular Expressions](./05_regular_expressions.md)
    5. [Shell Scripting](./06_shell_scripts.md)
    
    
## Getting To the Terminal

Getting to a shell terminal will be different based on whether or not you have a Mac or a Windows machine:

- **Mac:** Click the Terminal icon on your application bar at the bottom of your screen to open a shell terminal
- **Windows:** I recommend using a program like [MobaXterm](https://mobaxterm.mobatek.net/download.html){:target="_blank" rel="noopener"}. MobaXterm is an enhanced terminal for Windows machines that allows you to use Unix commands. This is useful because many bioinformatics programs are written for Unix shells! [Follow the download instructions here](https://mobaxterm.mobatek.net/download.html){:target="_blank" rel="noopener"} to get started!

## What is a Command?

When working with command line you will need to use "commands" to perform tasks. You can think of a command as the same thing as a click using a mouse - you can open a file/folder, move between folders, edit a file, etc.. However, commands on command line are far more powerful and can automate all this clicking. A command has the following structure:

```bash
command [options] [arguments]
```

!!! example "Explanation of Terms"

    - `command` name of the command you'd like to use
    - `[options]` options you can input to modify the original command
    - `[arguments]` the inputs you'd like the command to work with
    
## Data For Today's Tutorial

To get started let's copy some data to get started! So in command line enter the following code to copy data to your computer:

```bash
wget https://github.com/BioNomad/omicsTrain/blob/main/docs/programming_languages_tools/unix/unix_tutorial.zip
```

The `wget` command can be used to pull most files from the internet! However this file is compressed (ends in either `.zip`,`.tar`, or `.tar.gz`), meaning we can't access the files yet. To uncompress our data let's use the unzip command:

```bash
unzip unix_tutorial.zip
```
Great, now we are ready to get started!

## References

1. [Introduction to the Command Line for Genomics](https://datacarpentry.org/shell-genomics/01-introduction.html)
