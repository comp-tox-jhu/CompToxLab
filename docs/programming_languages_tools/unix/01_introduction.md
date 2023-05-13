# Introduction To Shell/Unix

Many software don't have a graphical user interface (or GUI) that users can use to run a pipeline. Instead many of these tools are only available through the command line interface. To use the command line interface you will need to know shell commands - which is exactly what we will cover here!

!!! abstract "Learning Objectives"

    1. Navigating Files & Directories
    2. Managing Files & Directories
    3. Redirection and Pipes
    4. Regular Expressions
    5. Project Organization
    
    
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
    



