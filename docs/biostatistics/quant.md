# Quantitative Variables

Quantitaive variables are numerical data: so variables such as height, weight, age. We can describe these variables with the following terms:

!!! example "Explanation of Terms"
    - **Minimum** : smallest value in your variable
    - **Median** : middle value in your variable
    - **Mean** : average value of your variable
    - **Max** : largest value in your variable
    - **Count** : how many values are in your variable
    - **Standard deviation** : measure of the spread of your variable

Let's see how to do this in our code:

```R
library(tidyverse)
meta <- read.table("./data/gbm_cptac_2021/data_clinical_patient.txt",
                     header = T,
                   sep="\t")

height.sum <- meta %>%
  summarise( 
    minimum = min(HEIGHT, na.rm = T),
    median = median(HEIGHT, na.rm = T),
    mean = mean(HEIGHT, na.rm = T),
    max = max(HEIGHT, na.rm = T),
    count = length(HEIGHT[!is.na(HEIGHT)]),
    SD = sd(HEIGHT, na.rm = T))

height.sum
```

```
  minimum median     mean max count       SD
1     150    170 169.6768 196    99 9.875581
```

!!! note
    You'll note here that we explicitly remove `NA` values to calculate these descriptive statistics.

Now that we know how to calculate our descriptive statistics, let's try and visualize our numeric data:

```R
ggplot(meta, aes(x=HEIGHT)) + 
  geom_histogram(fill="lightpink")+
  theme_bw()+
  labs(
    x="Height",
    y="Frequency",
    title=" Histogram of Height"
  )
```

![](images/histogram.png)

## References

- [BIOL202 Tutorials](https://ubco-biology.github.io/BIOL202/desc_cat_Var.html)
