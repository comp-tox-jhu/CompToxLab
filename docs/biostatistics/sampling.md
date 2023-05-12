# Sampling

When we try to assess an underlying population we often take samples of that population. Let's try and take a sample using the `sample()` function in R:

```R
library(tidyverse)
# load meta data
meta <- read.table("./data/gbm_cptac_2021/data_clinical_patient.txt",
                   header = T,
                   sep="\t")

## defined some population of ages
ages <- sample(meta$AGE,20)
```

If we wanted to take the same random sample we could use the `set.seed()` function:

```R
## grab the same sample
set.seed(123)
ages1 <- sample(meta$AGE,20)
ages2 <- sample(meta$AGE,20)
```

## Sampling Error

Not every sample is going to be a true approximation of the underline population. This difference is known as the sampling error.
What's assess our sample and see how it stacks up against our population:

```R
data.frame(
  Sample_Mean=mean(ages,na.rm = T),
  Population_Mean=mean(meta$AGE,na.rm = T)
)
```

```
  Sample_Mean Population_Mean
1        57.4        57.88889
```

Here we note that while similar to our true meta data mean, it is not exact. When we don't know the actual population mean we can get a whole range (or distribution) of means. The standard error of the mean is the measure of that sampling distribution:

$$\frac{\sigma}{\sqrt{N}}$$

!!! example "Explanation of Terms"
    - $\sigma$ Standard deviation of the sample
    - $N$ Number of observations in the sample

!!! tip "Math Tip"
    We can see that increasing the size of the sample, decreases the standard error of the mean.

## References

- [BIOL202 Tutorials](https://ubco-biology.github.io/BIOL202/desc_cat_Var.html)
