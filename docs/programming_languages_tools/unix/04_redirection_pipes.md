# Redirection and Pipes

## Redirection

We can move around the output of our commands to be the input of another command. For instance we can take the output of one file and add it to another with the `>` symbol:

```bash
cat new_file > super_new_file
cat super_new_file
```

!!! info "output" 

    ```bash
    Hello World
    ```
    
Here we looked at the contents of `new_file`  and put those contents in a new file which we called `super_new_file`. However, be careful with the `>` symbol as it will overwrite the contents of files. If we just wanted to append the contents we can use `>>` symbol (much safer option):

```bash
cat new_file >> super_new_file
cat super_new_file
```

!!! info "output" 

    ```bash
    Hello World
    Hello World
    ```

## Piping

So far we have entered commands separately, but we can actually take the output of one command and use it as the input for another with the `|` symbol! Here we will use the `cat` command to look at the contents of `meta.txt` then use the command `wc`  (word count) to get the number of lines, words and characters

```bash
cat new_file | wc 
```


!!! info "output" 

    ```bash
    1       2      12
    ```
    
Just to note, that you can use `wc new_file` to get the same result, this is just an example of how to pipe the output of one command into another. 

## References

1. [Introduction to the Command Line for Genomics](https://datacarpentry.org/shell-genomics/04-redirection.html)
