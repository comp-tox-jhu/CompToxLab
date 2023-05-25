# R Cheatsheet 

!!! abstract "This Cheatsheet is taken Directly From [Base R Cheat Sheet](https://iqss.github.io/dss-workshops/R/Rintro/base-r-cheat-sheet.pdf)"

## Getting Help

**Accessing the help files**

```R
# Get help of a particular function.
?mean
# Search the help files for a word or phrase.
help.search(‘weighted mean’)
# Find help for a package. 
help(package = ‘dplyr’)
```

**More about an object**

```R 
# Get a summary of an object’s structure.
str(iris)
# Find the class an object belongs to.
class(iris)
```

## Using Libraries

```R
install.packages(‘dplyr’)
# Download and install a package from CRAN.
library(dplyr)
# Load the package into the session, making all its functions available to use.
dplyr::select
# Use a particular function from a package.
data(iris)
# Load a built-in dataset into the environment. 
```


## Working Directory

```R
getwd()
# Find the current working directory (where inputs are found and outputs are sent).
setwd(‘C://file/path’)
# Change the current working directory. Use projects in RStudio to set the working directory to the folder you are working in. 
```

## Vectors

### Creating Vectors


||||
|-|-|-|
| `c(2, 4, 6)`| 2 4 6 | Join elements into a vector| 
|`2:6` |2 3 4 5 6 |An integer sequence |
|`seq(2, 3, by=0.5)` |2.0 2.5 3.0 |A complex sequence|
|`rep(1:2, times=3)`| 1 2 1 2 1 2| Repeat a vector |
|`rep(1:2, each=3)`| 1 1 1 2 2 2 |Repeat elements of a vector |

### Vector Functions

```R
sort(x)
# Return x sorted.
rev(x)
# Return x reversed.
table(x)
# See counts of values.
unique(x)
# See unique values.
```

### Selecting Vector Elements 

**By Position** 

```R
x[4] # The fourth element.
x[-4] # All but the fourth.
x[2:4] # Elements two to four.
x[-(2:4)] # All elements except two to four.
x[c(1, 5)] # Elements one and five.
```

**By Value**

```R
x[x == 10] # Elements which are equal to 10.
x[x < 0] # All elements less than zero.
x[x %in% c(1, 2, 5)] # Elements in the set 1, 2, 5.
```
**Named Vectors**

```R
x[‘apple’] Element with name ‘apple’.
```

## Programming

### For Loop

```R
for (variable in sequence){
Do something
}
```

**Example**

```R
for (i in 1:4){
j <- i + 10
print(j)
}
```

### While Loop

```R 
while (condition){
Do something
}
```

**Example**

```R
while (i < 5){
print(i)
i <- i + 1
}
```

### If Statements

```R
if (condition){
Do something
} else {
Do something different
}
```

**Example**

```R
if (i > 3){
print(‘Yes’)
} else {
print(‘No’)
}
```

### Functions

```R
function_name <- function(var){
Do something
return(new_variable)
}
```

**Example**

```R
square <- function(x){
squared <- x*x
return(squared)
}
```

## Reading and Writing Data


|Input |Ouput| Description|
|-|-|-|
|`df <- read.table(‘file.txt’)`| `write.table(df, ‘file.txt’)` |Read and write a delimited text file.|
| `df <- read.csv(‘file.csv’)`| `write.csv(df, ‘file.csv’)` |Read and write a comma separated value file. This is a special case of read.table and write.table. |
| `load(‘file.RData’`) |`save(df, file = ’file.Rdata’)` |Read and write an R data file, a file type special for R. |

## Conditions

|||
|-|-|
|`a == b`| Are equal |
|`a > b`| Greater than |
|`a >= b` |Greater than or equal to |
|`is.na(a)`| Is missing|
|`a != b`| Not equal |
|`a < b`| Less than |
|`a <= b`| Less than or equal to |
|`is.null(a)`| Is null |

## Data Types

||||
|-|-|-|
|`as.logical`| TRUE, FALSE, TRUE |Boolean values (TRUE or FALSE).|
|`as.numeric`| 1, 0, 1 |Integers or floating point numbers.|
|`as.character`| '1', '0', '1' |Character strings. Generally preferred to factors.|
|`as.factor`| '1', '0', '1', levels: '1', '0' |Character strings with preset levels. Needed for some statistical models |

## Maths Functions 

```R
log(x) # Natural log. sum(x) Sum.
exp(x) # Exponential. mean(x) Mean.
max(x) # Largest element. median(x) Median.
min(x) # Smallest element. quantile(x) Percentage quantiles.
round(x, n) # Round to n decimal places.
rank(x) # Rank of elements.
signif(x, n) # Round to n significant figures.
var(x) # The variance.
cor(x, y) # Correlation. 
sd(x) # The standard deviation
```

## Variable Assignment

```R
a <- 'apple'
a
```

!!! info "output"

    ```R
    [1] 'apple'
    ```
    
## The Environment

```R
ls() # List all variables in the environment.
rm(x) # Remove x from the environment.
rm(list = ls()) # Remove all variables from the environment.
```

## Matrixes

```R
m <- matrix(x, nrow = 3, ncol = 3) # Create a matrix from x.
m[2, ] #  Select a row
m[ , 1] #  Select a column
m[2, 3] #  Select an element
t(m) # Transpose
m %*% n#  Matrix Multiplication
solve(m, n) # Find x in: m * x = n
```

## Lists

```R
l <- list(x = 1:5, y = c('a', 'b')) #  A list is collection of elements which can be of different types.
l[[2]] # Second element of l.
l[1] # New list with only the first element. 
l$x  # Element named x. 
l['y'] # New list with only element named y.
```

## Data Frames

```R
df <- data.frame(x = 1:3, y = c('a', 'b', 'c')) # A special case of a list where all elements are the same length.
df$x # grab column of name x
df[[2]] # grab second column values
df[ , 2] # grab second column
df[2, ] # grab second row
df[2, 2] # grab the value in the second column and second row
nrow(df) # number of rows
ncol(df) # number of columns
dim(df) # number of rows and number of columns
cbind # bind by columns
rbind # bind by rows
```

## Strings 


```R 
paste(x, y, sep = ' ')     # Join multiple vectors together.
paste(x, collapse = ' ')   # Join elements of a vector together.
grep(pattern, x)           #  Find regular expression matches in x.
gsub(pattern, replace, x). #  Replace matches in x with a string.
toupper(x)                 # Convert to uppercase.
tolower(x)                 # Convert to lowercase.
nchar(x)                   #  Number of characters in a string. 
```

## Factors 

```R 
factor(x) # Turn a vector into a factor. Can set the levels of the factor and the order.
cut(x, breaks = 4) # Turn a numeric vector into a factor but ‘cutting’ into sections. 
```


## Statistics 

```R
lm(x ~ y, data=df) # Linear model.
glm(x ~ y, data=df) # Generalised linear model.
summary # Get more detailed information out a model.
t.test(x, y) # Preform a t-test for difference between means.
pairwise.t.test # Preform a t-test for paired data.
prop.test # Test for a difference between proportions.
aov # Analysis of variance. 
```


## Distributions

||Random Variates |Density Function |Cumulative Distribution| Quantile|
|-|-|-|-|-|
|Normal |`rnorm `|dnorm` |`pnorm` |`qnorm`|
|Poison |`rpois` |`dpois`| `ppois` |`qpois`|
|Binomial |`rbinom` |`dbinom` |`pbinom` |`qbinom`|
|Uniform |`runif` |`dunif `|`punif` |`qunif`|

