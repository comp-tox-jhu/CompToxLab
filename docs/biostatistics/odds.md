# Odds

We often hear about the odds of something happening, but what does this mean? Well, the odds would be:

$$O = \frac{p}{1-p}$$

!!! example "Explanation of Terms"
    - $O$ the odds of something happening
    - $p$ the probability that thing happens
    
So colloquially, we see that the odds of an event is the probability of that event happening divided by the 
probability of that event not happening. We are going to assess something called lost to follow up. In the setting of a clinical 
trial this term describes losing patients due to various reasons. Let's try using our data to calculate the odds of losing a male patient to 
follow up:

```R
library(tidyverse)
# load meta data
meta <- read.table("./data/gbm_cptac_2021/data_clinical_patient.txt",
                   header = T,
                   sep="\t")

# what are the odds of being a female non-smoker in our dataset?
table = as.data.frame.matrix(
  table(meta$SEX,meta$LOST_TO_FOLLOW_UP)
) %>%
  mutate(row_totals = apply(.,1,sum)) %>%
  rbind(t(data.frame(column_totals=apply(., 2, sum))))

table
```

```
              No Yes row_totals
Female        33  10         43
Male          44  10         54
column_totals 77  20         97
```

Here we can have created what is called a **contingency table** or table that describes the frequency distribution of variables. We see that more patients are **not** lost to follow up. Let's calculate the odds now!

```R
male <- table %>%
  filter(rownames(.) == "Male") %>%
  mutate(
    prob_male_lost = Yes/row_totals,
    prob_male_not_lost = No/row_totals,
    odds_male_lost = prob_male_lost/(1-prob_male_lost),
    odds_male_not_lost = prob_male_not_lost/(1-prob_male_not_lost))

male
```

```
     No Yes row_totals prob_male_lost prob_male_not_lost odds_male_lost odds_male_not_lost
Male 44  10         54      0.1851852          0.8148148      0.2272727                4.4
```

So here we see that the odds of losing a male to follow up are 0.23 to 1. An alternative way of looking at this is that the odds of **not** losing a male to follow up are 4.4 to 1.

## References

1. [BIOL - 202: Analyzing a single categorical variable](https://ubco-biology.github.io/BIOL202/estproportions.html)
2. [Lost To Follow Up](https://en.wikipedia.org/wiki/Lost_to_follow-up)
3. [Contingency Table](https://en.wikipedia.org/wiki/Contingency_table)
