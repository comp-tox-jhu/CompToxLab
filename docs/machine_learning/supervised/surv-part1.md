## Survival Analysis Part 1

Survival Analysis refers to a statistical framework to assess the time it takes for some event of interest to occur. In medical research this framework
is often used to ask questions about:

- a patient's probability of survival given some time period
- the characteristics that might impact survival
- the differences in survival between groups

## Events and Censoring

The type of event can include:

- Relapse
- Death 
- Progression

Now not all patients in a study will experience the same type of event described above and as such these observations need to be *censored* as to not 
influence our results. Censoring can be caused by:

- a patient not yet experiencing an event of interest (Relapse, Death, Progression, etc.)
- a patient being lost to follow up during the study (meaning they drop out for one reason or another)
- a patient expreiences an alternative event to the event of interest

Given that these are assessed at the end of the study, this is referred to as **right censoring**.

## Kaplan-Meier Survival Probability

The probability that the event of interest **will not happen** in a given time period is referred to as the survival probability. We can calculate this
using the non-parametric Kaplan-Meier method:

$$S(t_i) = S({t_i} - 1) (1 - \frac{d_i}{n_i})$$

$$S(0) = 1$$

$$t_0 = 0$$

!!! example "Explanation of Terms"

    - $S(t_i)$ Survival probability at time $t_i$
    - $S({t_i} - 1)$ Survival probability at time ${t_i} - 1$
    - $n_i$ number of patients that have not experienced event of interest right before time $t_i$
    - $d_i$ number of events of interest at $t_i$
    - $t_0 = 0$ your starting time must be 0
    - $S(0) = 1$ your probability of survival at time 0 is 1 
    
## Pre-Processing

Let's start by loading our meta data and ensuring that we have a censored column:

```R
# load the libraries
.libPaths(c("/cluster/tufts/hpc/tools/R/4.0.0"))
library("survival")
library("survminer")

# load in patient meta data
meta <- read.csv(
  file = "data/gbm_cptac_2021/data_clinical_patient.txt",
  skip=4,
  header = T,
  sep = "\t"
)

# create a censored status column
# 1 being they are censored
# 2 being they are not censored
meta$status <- ifelse(
  (meta$LOST_TO_FOLLOW_UP == "Yes" & meta$VITAL_STATUS == "Living"),
  1,
  2
)

```

## Fitting a Survival Curve

To create our survival curve we will fit a model to determine: What is the survival probability of male and female patients and how do these curves compare?

```R
# fit our survival curve to our data
# time is days
# status is our censor status
# and SEX is our variable of interest
fit <- survfit(
  Surv(PATH_DIAG_TO_DEATH_DAYS, status) ~ SEX,
  data = meta)

# examine our fit
summary(fit)$table
```

```
           records n.max n.start events   *rmean *se(rmean) median 0.95LCL 0.95UCL
SEX=Female      27    27      27     27 308.0741   47.10662    274     147     401
SEX=Male        35    35      35     35 399.2286   39.34208    381     294     511
```

!!! info "What does this table tell us?"

    - First we start off with 27 females and 35 males
    - We then see that we observe the event (death) 27 times in females and 35 times in males
    - the mean survival time in days and it's standard error (`*rmean` and `*se(rmean)`)
    - the median survival time in days (`median`)
    - the confidence interval around our parameter

## Visualizing Survival Curves 

Now if we want to visualize these survival curves we can use the following code:

```R
# plot the survival curves
ggsurvplot(fit,
           pval = TRUE, # add in a p-value
           conf.int = TRUE, # add in a confidence interval
           risk.table = TRUE, # add in a risk table to our plot
           risk.table.col = "strata", # color the risk table by group
           linetype = "strata", # ensure lines are made per group
           surv.median.line = "hv", # add in median survival time line
           ggtheme = theme_bw(), # set theme to black and white
           palette = c("#C06C84","#355C7D")) # change colors of lines
```

![](images/surv-curve-sex.png)

!!! info "What does this graph tell us?"

    - Here we see time against survival probability
    - We also note a risk table down below that lists the number of individuals at risk for the event
    - Given a p-value above 0.05 we do not have enough evidence to rule out the possibility that any difference between the male and female survival curve is due to chance. 
    
    
## Log-Rank Test

We can see above that we asked is there any significant difference between these two curves. The test conducted to determine this is called the log rank-test. This test is non-parametric, meaning it does not make any assumptions about the distribution of survival times. If we were solely interested in testing the difference between these two curves we could use the following code:

```R
survdiff(Surv(PATH_DIAG_TO_DEATH_DAYS, status) ~ SEX, data = meta)
```

```
Call:
survdiff(formula = Surv(PATH_DIAG_TO_DEATH_DAYS, status) ~ SEX, 
    data = meta)

{==n=62==}, {==37==} observations deleted due to missingness.

            N Observed Expected (O-E)^2/E (O-E)^2/V
SEX=Female 27       27     22.6     0.876      1.41
SEX=Male   35       35     39.4     0.501      1.41

 Chisq= {==1.4==}  on 1 degrees of freedom, {==p= 0.2==} 
```

!!! info "What do these results mean?"
    - Here we see that we started with 99 observations but ended up with 37 after removing our censored data points
    - Our chi-square test statistic was 1.4
    - Our p-value is 0.2, again indicating we do not have enough evidence to rule out the possibility that any difference between the male and female survival curve is due to chance. 


## References

1. [Survival Analysis Basics](http://www.sthda.com/english/wiki/survival-analysis-basics)
2. [Survival Analysis Part I: Basic concepts and first analyses](https://www.nature.com/articles/6601118)
3. [Nonparametric Estimation from Incomplete Observations](https://www.jstor.org/stable/2281868#metadata_info_tab_contents)
4. [Survival plots of time-to-event outcomes in clinical trials: good practice and pitfalls](https://www.sciencedirect.com/science/article/pii/S014067360208594X?via%3Dihub)
