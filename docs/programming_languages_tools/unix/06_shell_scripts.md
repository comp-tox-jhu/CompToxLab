# Writing Shell Scripts

Step 1: Create a new file
Open a text editor and create a new file. Save it with a .sh extension, for example, myscript.sh. This will be your Unix script.

Step 2: Set up the script
To begin, add the following line at the beginning of your script to specify the shell you want to use:

bash
Copy code
#!/bin/bash
Step 3: Define variables
You can declare variables to store values that you'll use throughout your script. Variables in bash don't have a specific type and are defined without any prefixes. Here's an example:

bash
Copy code
#!/bin/bash

# Define a variable
name="John"
age=30
Step 4: Use variables
You can access the values stored in variables using the $ symbol. Here's an example of printing the values of the variables we defined:

bash
Copy code
#!/bin/bash

name="John"
age=30

# Print the values
echo "Name: $name"
echo "Age: $age"
Step 5: Define functions
Functions allow you to group a set of commands together and reuse them in your script. Here's an example of defining a function:

bash
Copy code
#!/bin/bash

# Define a function
greet() {
  echo "Hello, $1!"
}
In this example, $1 represents the first argument passed to the function.

Step 6: Call functions
To execute a function, simply call it by its name followed by parentheses. Here's an example of calling the greet function we defined earlier:

bash
Copy code
#!/bin/bash

greet() {
  echo "Hello, $1!"
}

# Call the function
greet "Alice"
Step 7: Process command-line arguments
You can pass arguments to your script when you run it from the command line. The arguments are accessible using special variables: $1, $2, $3, and so on. Here's an example:

bash
Copy code
#!/bin/bash

# Process command-line arguments
echo "First argument: $1"
echo "Second argument: $2"
Step 8: Implement loops
Loops allow you to repeat a set of commands multiple times. There are various types of loops in bash, but let's focus on two common ones: the for loop and the while loop.

Example of a for loop:

bash
Copy code
#!/bin/bash

# Iterate over a list using a for loop
for item in "apple" "banana" "orange"
do
  echo "Fruit: $item"
done
Example of a while loop:

bash
Copy code
#!/bin/bash

# Use a while loop to print numbers from 1 to 5
counter=1
while [ $counter -le 5 ]
do
  echo "Number: $counter"
  counter=$((counter + 1))
done
Step 9: Make the script executable
Before you can run the script, you need to make it executable. Open a terminal and navigate to the directory where your script is located. Then, run the following command:

bash
Copy code
chmod +x myscript.sh
Step 10: Run the script
To execute your script, navigate to the directory where it's located and run the following command in the terminal:

bash
Copy code
./myscript.sh
