# Regular Expressions

Introduction to Regular Expressions:
Regular expressions are patterns used to match and manipulate text. They consist of a sequence of characters that define a search pattern. Unix provides several tools that support regular expressions, such as grep, sed, and awk.

Basic Regular Expression Syntax:
Before diving into the tools, let's cover the basic syntax of regular expressions:

Literal Characters: Most characters in a regular expression match themselves. For example, the regular expression "cat" matches the string "cat" exactly.

Character Classes: Square brackets ([]) define a character class, which matches any single character within the brackets. For example, the regular expression "[aeiou]" matches any vowel.

Metacharacters: Certain characters have special meaning in regular expressions. Some commonly used metacharacters are:

. (dot): Matches any single character.
* (asterisk): Matches zero or more occurrences of the preceding character or group.
+ (plus): Matches one or more occurrences of the preceding character or group.
? (question mark): Matches zero or one occurrence of the preceding character or group.
| (pipe): Acts as an OR operator, matching either the expression before or after it.
Anchors: Anchors are used to match the position within a line. Two commonly used anchors are:

^ (caret): Matches the start of a line.
$ (dollar sign): Matches the end of a line.
Using grep:
grep is a powerful command-line tool for searching files for lines that match a specific pattern. Here's how you can use grep with regular expressions:
perl
Copy code
grep pattern file
pattern represents the regular expression pattern you want to search for.
file specifies the file(s) you want to search in.
Example:
Let's say we have a file called "data.txt" containing the following lines:

bash
Copy code
apple
banana
cat
dog
To search for lines containing the word "cat," you can run the following command:

kotlin
Copy code
grep "cat" data.txt
Using sed:
sed is a stream editor that allows you to perform text transformations based on regular expressions. Here's a simple example of using sed with regular expressions:
arduino
Copy code
sed 's/pattern/replacement/g' file
pattern represents the regular expression pattern you want to search for.
replacement represents the text you want to replace the matched pattern with.
file specifies the file(s) you want to perform the replacement on.
Example:
Let's say we have a file called "data.txt" containing the following lines:

bash
Copy code
apple
banana
cat
dog
To replace all occurrences of "cat" with "lion," you can run the following command:

kotlin
Copy code
sed 's/cat/lion/g' data.txt
Using awk:
awk is a versatile programming language for manipulating structured data. It supports regular expressions to select and transform text. Here's a basic example of using awk with regular expressions:
arduino
Copy code
awk '/pattern/ { action }' file
/pattern/ represents the regular expression pattern you want to search for.
{ action } represents the action you want to perform on lines matching the pattern.
file specifies the file(s) you want to process.
