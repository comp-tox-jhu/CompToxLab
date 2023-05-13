Unix allows you to redirect input and output streams and use pipes to chain commands together:
>: Redirect output to a file (overwrite).
command > [file]: Redirect the output of a command to the specified file.
>>: Redirect output to a file (append).
command >> [file]: Append the output of a command to the specified file.
|: Use the output of one command as the input for another.
command1 | command2: Redirect the output of command1 as input to command2.
