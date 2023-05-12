## Fisher's Exact Test

In the [odds ratio topic note](odds-ratio-risk.md) we noted that the odds ratio could help us determine whether the odds 
were greater in one group versus another. We can test the strength of association of the group and event with **Fisher's Exact Test**. 
Fisher's Exact Test has the following hypotheses:

- $H_0$ or null hypothesis: there is no association between the group and event (Odds ratio = 1)
- $H_a$ or alternative hypothesis: there is an association between the group and event (Odds ratio != 1)

We can calculate the probability given the following contingency table with:

||group 1|group 2|
|-|-|-|
|Event |a |b |
|No Event |c  |d |


$$p = \frac{(a+b)!(c+d)!(a+c)!(b+d)!}{a!b!c!d!(a+b+c+d)!}$$

!!! example "Explanation of Terms"

    - $a$ number of members in group 1 with event
    - $b$ number of members in group 2 with event
    - $c$ number of members in group 1 without event
    - $d$ number of members in group 2 without event
    

Let's assess the relationship between pathiet sex and losing patients to follow up:

```R
library(tidyverse)
# load meta data
meta <- read.table("./data/gbm_cptac_2021/data_clinical_patient.txt",
                   header = T,
                   sep="\t")

# What is the frequency distribution of losing males/females to follow up
table <- as.data.frame.matrix(
  table(meta$SEX,meta$LOST_TO_FOLLOW_UP)
)

table
```

```
       No Yes
Female 33  10
Male   44  10
```

Before we continue, we need to make this table match the contingency table above. With the rows being the event and the columns being the group:
    
```R
# Reorder so that we are assessing the odds ratio of losing patients to follow up
table <- as.data.frame.matrix(
  table(meta$SEX,meta$LOST_TO_FOLLOW_UP)
) %>%
  select(c(Yes,No)) %>%
  t()

table
```

```
    Female Male
Yes     10   10
No      33   44
```

Now let's conduct our hypothesis test:

```R
#apply our test
fisher.test(table)
```

```
	Fisher's Exact Test for Count Data

data:  table
p-value = {==0.6191==}
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 {==0.4395952 4.0274603==}
sample estimates:
odds ratio 
   {==1.32933==} 
```

!!! info "Explanation of Results"

    Here we note:
    
    - our p-value is above 0.05 and thus not strong enough to reject the null (a.k.a the true odds ratio is equal to 1)
    - the 95% confidence interval reveals that our true odds ratio is somewhere between `0.4395952` and `4.0274603`
    - our odds ratio that patient sex is associated with losing the patient to follow up is about `1.3`

## References

1. [BIOL 202](https://ubco-biology.github.io/BIOL202/fishertest.html)
2. [Fisher's Exact Test](https://en.wikipedia.org/wiki/Fisher%27s_exact_test)
