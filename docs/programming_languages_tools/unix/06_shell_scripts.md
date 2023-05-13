# Writing Shell Scripts

## Creating a Script

Shell scripts are convenient ways of running multiple commands. These scripts usually end with `.sh`. Let's use nano to create a shell script that says the date. First move to our `scripts_folder`:

```bash 
cd ../scripts_folder/
```

Now create a script called "first_script.sh" with the following contents:


!!! info "first_script.sh"

    ```bash
    #!/bin/bash
    
    # run the date command
    date
    ```

You'll note we use `#!/bin/bash` to make sure this is a bash shell script. When we start using other programming languages (like python) this will be swapped out with `#!/bin/ptyhon`. We also have a line that starts with a `#`. This is called a comment. Comments aren't read as code, they are just plain text we can add to give our commands/code context. Comments are very useful in making your code readable. We should also note that we can't run our script just yet. To run our script we will need to change the permissions on it with the `chmod` command:

```bash
chmod +x first_script.sh 
```

This command allows us to *execute* our script. To execute our script we can enter the following:

```bash
./first_script.sh 
```

!!! info "first_script.sh"

    ```bash
    Sat May 13 13:00:45 EDT 2023
    ```
    
## Variables

Now let's make our script more complicated with a variable (basically a word that represents some value):

!!! info "first_script.sh"

    ```bash
    #!/bin/bash
    
    # define a variable
    initial="Today is:"
    
    # echo the variable and date
    echo $initial 
    date
    ```
    
Here we see that we have text ("Today is:") that is assigned to `initial`, which we reference with a `$` sign. To run this script we will enter the following:


```bash
./first_script.sh 
```

!!! info "first_script.sh"

    ```bash
    Today is:
    Sat May 13 13:00:45 EDT 2023
    ```

## Functions

We can create a function using shell scripting with the following syntax:

```bash
function_name (){
    command
}
```

To call that function we write out function_name(). Let's modify our script to add a function!

!!! info "first_script.sh"

    ```bash
    #!/bin/bash
    
    # define a variable
    initial="Today is:"
    
    # echo the variable and date
    echo $initial 
    date
    
    # define a function to print the first column of the first argument you give the script
    first_col() {
        awk '{print $1}' $1
    }
    
    # call this function
    first_col $1
    ```
    
In our function, we ask it to print the first column of the first argument passed to the function. When we call the function, we ask it to apply this function to the first argument we give the script itself. So let's run this script with the meta.txt file as our argument:

```bash
./first_script.sh
```


!!! info "output"

   ```bash
   Today is:
   Sat May 13 13:30:30 EDT 2023
   Run
   SRR1219879
   SRR1219880
   ```
   
## Loops

We can use loops to perform repetitive tasks. Loops follow the following syntax:

```bash
for variable in list
do
    command $variable
done
```

Let's modify our script to use a loop to examine the first line of all files in working directory:

!!! info "first_script.sh"

    ```bash
    #!/bin/bash
    
    # define a variable
    initial="Today is:"
    
    # echo the variable and date
    echo $initial 
    date
    
    # define a function to print the first column of the first argument you give the script
    first_col() {
        awk '{print $1}' $1
    }
    
    # call this function
    first_col $1
    
    # make a loop to examine the first line of all files in working directory
    for file in *
    do 
        head -n 1 $file
    done
    ```

To run the script we will enter the following:

```bash
./first_script.sh
```

!!! info "output"

   ```bash
   Today is:
   Sat May 13 13:30:30 EDT 2023
   Run
   SRR1219879
   SRR1219880
   #!/bin/bash
   ```


## References

1. [Introduction to the Command Line for Genomics](https://datacarpentry.org/shell-genomics/05-writing-scripts.html)
