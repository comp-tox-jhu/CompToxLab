# Regular Expressions

## Syntax

Regular expressions are a great way to find/manipulate patterns in files. Regular expressions follow a specific syntax:

- `.` used for matching any character
- `*` used for mathcing zero or more times
- `+`  used for matching one or more times
- `?` used for matching zero or one times
- `|` or operator
- `^` used for matching the start of a line
- `$` used for matching at the end of a line

So if we wanted to list all files that end with `txt` we would use:

```bash
ls *txt
```

!!! info "output"

```bash
accList.txt	meta.txt
```

However, to use the rest of the symbols we need to use special commands which we will describe below.

## Search Files

Regular expressions are a great way to find/manipulate patterns in files. Let's learn how to search for files starting with `a` with `grep`:

```bash
ls | grep ^a
```

!!! info "output"

    ```bash
    accList.txt
    ```

We can also search inside files for lines with certain patterns. Let's try to find the line that contains the word "Run" in our meta data file:


```bash
grep "Run" meta.txt 
``` 

!!! info "output"

    ```bash
    Run	analyte_type	Assay.Type	body_site
    ```
    
## Replace Patterns

Sometimes you may want to replace a pattern in a file and we can accomplish this with the `sed` command! Let's replace the word "body_site" with "tissue" in our `meta.txt` file:

```bash
sed s/body_site/tissue/g meta.txt
```

!!! info "output"

    ```bash
    Run	analyte_type	Assay.Type	tissue
    SRR1219879	DNA	WGS	Peripheral blood
    SRR1219880	DNA	WGS	Peripheral blood
    ```
    
The `sed` command follows the pattern `s/pattern/replacement/g`.

## Manipulating Structured Data

If we have structured data, like a data frame (think csv, tsv, excel file), then we can manipulate the data with `awk` which follows the following pattern:

```bash
awk '/pattern/ { action }' file
```

Let's print out the first column of our `meta.txt` file with `awk` as an example without a pattern:

```bash
awk '{print $1}' meta.txt
```

!!! info "output"

    ```bash
    Run
    SRR1219879
    SRR1219880
    ```
    
You can pick other columns with the `$` and the number of the column. If you choose `$0` the entire file will be printed. For example if we wanted to print any line that had the pattern `SRR` we would use:


```bash
awk '/SRR/ {print $0}' meta.txt 
```

!!! info "output"

    ```bash
    Run
    SRR1219879	DNA	WGS	Peripheral blood
    SRR1219880	DNA	WGS	Peripheral blood
    ```

## References

1. [Regular Expressions Tutorial](https://ryanstutorials.net/regular-expressions-tutorial/)
2. [AWK Quick Guide](https://www.tutorialspoint.com/awk/awk_quick_guide.htm)
