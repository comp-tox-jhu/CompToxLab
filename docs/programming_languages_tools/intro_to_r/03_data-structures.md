## Data Structures

So we have all this lovely data to play with and in R we typically organize in a few ways:

## Vectors

Vectors are collections of data, like a:

  collection of numbers - ```c(1,2,3)```
  collection of characters -  ```c("1","2","3")```
  collection of logical values - ```c(TRUE,FALSE,TRUE)```

!!! note
   It should be noted that a vector needs to be a collection of the **same type** of data. You will also note that each list is separated by commas and surrounded by ```c()```. This is necessary to create vectors so make sure to remember the ```c()```!

On the previous page we mentioned that values in vectors can be selected by their value (e.g. `someVector[4]` to grab the forth element) or their name if they are named (e.g. `someVector["a"]` to grab the element named `a`)

## Strings

When dealing with character data there are a few functions that are helpful when manipulating them. Let's create a few character vectors to play with:

```R
x=c("a","b","c")
y=c("d","e","f")
```


```R
paste(x, y, sep = '_') # join two or more vectors together
```

!!! info "output"

    ```R
    [1] "a_d" "b_e" "c_f"
    ```
    
```R
paste(x, collapse = '_') # collapes values in a vector together
```

!!! info "output"

    ```R
    [1] "a_b_c"
    ```

```R    
grep("a", x) # Find regular expression matches in x
```

!!! info "output"

    ```R
    [1] 1
    ```

Here we note that the pattern `a` was found 1 time. We can also replace patterns as well

```R
gsub("a", "A", x) # Replace pattern matches in x 
```

!!! info "output"

    ```R
    [1] "A" "b" "c"
    ```
    
```R    
toupper(x) # Convert to uppercase
```

!!! info "output"

    ```R
    [1] "A" "B" "C"
    ```
    
```R    
tolower(x) # Convert to lowercase.
```

!!! info "output"

    ```R
    [1] "a" "b" "c"
    ```

```R
nchar(x) # Number of characters in a string
```

!!! info "output"

    ```R
    [1] 3
    ```

### Matrices

A matrix can be created by combining vectors of the **same length and same data type**. They are used frequently when performing operations on numeric data but can include other data types. In R we can create a matrix with the `matrix()` function:

```R
m <- matrix(data=1:9,nrow = 3,ncol=3)
m
```

!!! info "output"

    ```
         [,1] [,2] [,3]
    [1,]    1    4    7
    [2,]    2    5    8
    [3,]    3    6    9
    ```

Here we take a vector and specify how many columns and how many rows we'd like. To select elements in a matrix we can use the following:

```R
m[2, ] # Select a row
```

!!! info "output"

    ```
     [1] 2 5 8
    ```
    
```R
m[ , 1] # Select a column
```

!!! info "output"

    ```
    [1] 1 2 3
    ```
    
```R
m[2, 3] # Select an element
```

!!! info "output"

    ```
    [1] 8
    ```
    
```R
nrow(m) # number of rows 
```

!!! info "output"

    ```
    [1] 3
    ```
    
```R
ncol(m) # number of columns 
```

!!! info "output"

    ```
    [1] 3
    ```
    
```R
dim(m) # number of rows then columns 
```

!!! info "output"

    ```
    [1] 3 3
    ```    
    
### Data Frames

Data frames are also collections of vectors of the **same length**. However, they do not need to be the same data type. Here we create a data.frame with the `data.frame()` function:

```R
df <- data.frame(
  characters=c("past","present","future"),
  numbers=c(1,2,3),
  logical=c(TRUE,FALSE,TRUE),
  integer=c(1L,2L,3L)
)
```

!!! info "output"

    ```
      characters numbers logical integer
    1       past       1    TRUE       1
    2    present       2   FALSE       2
    3     future       3    TRUE       3
    ```

We can manipulate data frames the same way we manipulate matrices. However, we have a short hand option to grab a column:

```R 
df$characters # select the column with the title characters
```

!!! info "output"

    ```R
    [1] "past"    "present" "future" 
    ```

### Lists

Lists are collections of data that **do not** need to be the same type or length. We can create lists with the `list()` function:

```R
l <- list(
  data.frame=data.frame(numbers=1:3,characters=c("past","present","future")),
  numbers=1:5,
  characters=c("past","present","future")
)
```

!!! info "output"

    ```
    $data.frame
      numbers characters
    1       1       past
    2       2    present
    3       3     future
    
    $numbers
    [1] 1 2 3 4 5
    
    $characters
    [1] "past"    "present" "future" 
    ```

```R
l[2] # list with just the second thing in the list
```


!!! info "output"

    ```R
    $numbers
    [1] 1 2 3 4 5
    ```

```R
l[[2]] # select the second thing in the list
```

!!! info "output"

    ```R
    [1] 1 2 3 4 5
    ```
    
    
```R
l$characters # grab the element named characters
```


!!! info "output"

    ```R
    [1] "past"    "present" "future" 
    ```
    
    
```R
l['numbers'] # list with just the item named numbers


!!! info "output"

    ```R
    $numbers
    [1] 1 2 3 4 5
    ```
    

