## Risk

We have defined odds as the probability of that event happening divided by the probability of that event not happening. 
In our topic note on [odds](odds.md) we asked what were the odds of losing a male patient to follow up. This is different from **risk**.
Risk can be defined as:


$$Risk = \frac{n_i}{N}$$

!!! example "Explanation of Terms"
    - $n_i$ number of times event happened
    - $N$ total number of events
    
Here we see this is just a proportion or the number of times an event happened divided by the total number of events! Let's assess the risk of 
losing a male patient to follow up:

```R
#risk
library(tidyverse)
# load meta data
meta <- read.table("./data/gbm_cptac_2021/data_clinical_patient.txt",
                   header = T,
                   sep="\t")

# what are the odds of being a female non-smoker in our dataset?
table <- as.data.frame.matrix(
  table(meta$SEX,meta$LOST_TO_FOLLOW_UP)
) %>%
  mutate(row_totals = apply(.,1,sum)) %>%
  rbind(t(data.frame(column_totals=apply(., 2, sum))))

table


male <- table %>%
  filter(rownames(.) == "Male") %>%
  mutate(
    risk_male_lost = Yes/row_totals,
    risk_male_not_lost = No/row_totals)

male
```

```
     No Yes row_totals risk_male_lost risk_male_not_lost
Male 44  10         54      0.1851852          0.8148148
```

We see that the risk of losing a male patient to follow up is about 18%. 

## Relative Risk

You might hear about two risks how relate to each other - also called the **relative risk**. Relative risk can be calculated by:

$$Relative Risk = \frac{Risk_i}{Risk_j}$$

Where :

$$Risk_i = \frac{n_i}{N_i}$$

$$Risk_j = \frac{n_j}{N_j}$$

!!! example "Explanation of Terms"
    - $n_i$ number of times event happens in group i
    - $N_i$ total number of group i members
    - $n_j$ number of times event happens in group j
    - $N_j$ total number of group j members
    
    
We will calculate the relative risk of losing a female patient to follow up versus a male patient. However, before we do so we will visualize 
our data:

```R
# visualize risks:

as.data.frame.matrix(
  table(meta$SEX,meta$LOST_TO_FOLLOW_UP)) %>%
  t() %>%
  reshape2::melt() %>%
  ggplot(.,aes(x=Var2,y=value,fill=Var1))+
  geom_bar(stat = "identity",position="fill") +
  scale_y_continuous(labels = scales::percent)+
  theme_bw()+
  scale_fill_manual(values=c("aquamarine3","lightpink3")) +
  labs(
    x="",
    y="Frequency",
    fill="Lost To Follow Up?"
  )
```

![](images/risk-lost-to-follow.png)

Here we see that females are slighly more prone to being lost to follow up. How does this translate to relative risk?

```R
risks <- table %>%
  filter(rownames(.) != "column_totals") %>%
  mutate(
    risk_lost = Yes/row_totals,
    risk_not_lost = No/row_totals)
    
risks
```

```
       No Yes row_totals risk_lost risk_not_lost
Female 33  10         43 0.2325581     0.7674419
Male   44  10         54 0.1851852     0.8148148
```

Here we can eyeball that the risk of losing a female patient to follow up is greater than losing a  male patient. Let's see what the relative risk is:

```R
relative_risk <- risks$risk_lost["Female"]/risks$risk_lost["Male"]

relative_risk
```

```
  Female 
1.255814 
```

Here we can guage from the relative risk that being a female patient increases the risk of losing the patient to follow up. 

## Odds Ratio

The **odds ratio** of two events is also used as a measure to compare events between groups. However, the odds ratio is the ratio of odds of an event in one group versus another and is defined by:

$$Odds\ Ratio = \frac{Odds_i}{Odds_j}$$

Where:

$$Odds_i = \frac{n_i}{n_{not\ i}}$$

$$Odds_j = \frac{n_j}{n_{not\ j}}$$

!!! example "Explanation of Terms"
    - $Odds_i$ odds event happens in group i
    - $Odds_j$ odds event happens in group j
    - $n_i$ number of times event happens in group i
    - $n_{not\ i}$ number of times event doesn't happen in group i
    - $n_j$ number of times event happens in group j
    - $n_{not\ j}$ number of times event doesn't happen in group j
    
So let's calculate the odds ratio of losing a female patients to follow up, versus male patients to follow up:

```R
odds <- table %>%
  filter(rownames(.) != "column_totals") %>%
  mutate(
    odds = Yes/No)
odds
```

```
       No Yes row_totals odds_lost
Female 33  10         43 0.3030303
Male   44  10         54 0.2272727
```

Here we note that the odds of losing a female patient to follow up are higher than the odds of losing a female. The odds ratio would then be:

```R
odds_ratio <- odds$odds[rownames(odds)=="Female"]/odds$odds[rownames(odds)=="Male"]
odds_ratio
```

```
[1] 1.333333
```

Here we note that the odds of losing a female patient versus a male patient to follow up are 1.3 to 1. We can see that the odds ratio is similar to the relative risk, but the values are nonetheless different. 

!!! tip "Relative Risk/Odds Ratio Interpretation"
    - Odds Ratio/Relative Risk = 1 the factor does not affect the event
    - Odds Ratio/Relative Risk < 1 the factor decreases the risk of the event (protective factor)
    - Odds Ratio/Relative Risk > 1 the factor increases the risk of the event (risk factor)

1. [Relative Risk](https://en.wikipedia.org/wiki/Relative_risk)
2. [Odds Ratio](https://en.wikipedia.org/wiki/Odds_ratio)
3. [BIOL 202](https://ubco-biology.github.io/BIOL202/oddsratio.html)
4. [Common pitfalls in statistical analysis: Odds versus risk](Common pitfalls in statistical analysis: Odds versus risk)
