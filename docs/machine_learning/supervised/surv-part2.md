## Survival Analysis Part 2

In the [first Survival Analysis topic note](surv-part1.md) we assessed survival probability curves and if two different survival curves were significantly different. Instead of survival probability, here we will discuss **hazard probability**, or the probability that an individual has an event at some time $t$. Survival probability on the other hand, refers to the opposite - that the individual will survive to time $t$. To define the hazard probability we will use the following function:

$$h(t)=h_0(t)×exp(b_{1}x_{1}+b_{2}x_{2}+...+b_{n}x_{n})$$

!!! example "Explanation of Terms"
    - $t$ the survival time
    - $h(t)$ the hazard probability
    - $n$ different number of covariates $x$ (so variables like age, sex, etc.)
    - $b$ are the coefficients of the impact of each of these covariates 
    - $h_0(t)$ would be the baseline hazard

These $b$ coefficients can also be referred to as hazard ratios. And they can be interpreted like so:

!!! info "Hazard Ratio Interpretation"

    - Hazard Ratio < 1: there is a decrease in the hazard
    - HHazard RatioR > 1: there is a increase in the hazard
    - Hazard Ratio = 1: there is no effect on the hazard
    
## Pre-Processing

We will start as we did in the [first Survival Analysis topic note](surv-part1.md) and load our libraries, data, and add a censor status column:

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

## Cox Proportional-Hazards Model

```R
# run the Cox Proportional-Hazards Model on our data 
# with a few covariates and get a summary
cox <- coxph(Surv(PATH_DIAG_TO_DEATH_DAYS, status) ~ SEX + AGE + BMI,
             data = meta)
summary(cox)
```

```
Call:
coxph(formula = Surv(PATH_DIAG_TO_DEATH_DAYS, status) ~ SEX + 
    AGE + BMI, data = meta)

  n= 62, number of events= 62 
   (37 observations deleted due to missingness)

             coef exp(coef)  se(coef)      z Pr(>|z|)
SEXMale {==-0.300913==}  0.740143  0.266021 -1.131    {==0.258==}
AGE     {==-0.003842==}  0.996166  0.011765 -0.327    {==0.744==}
BMI     {==-0.006180==}  0.993839  0.025223 -0.245    {==0.806==}

        exp(coef) exp(-coef) lower .95 upper .95
SEXMale    0.7401      1.351    0.4394     1.247
AGE        0.9962      1.004    0.9735     1.019
BMI        0.9938      1.006    0.9459     1.044

Concordance= 0.564  (se = 0.041 )
Likelihood ratio test= 1.52  on 3 df,   p=0.7
Wald test            = 1.55  on 3 df,   p=0.7
Score (logrank) test = 1.56  on 3 df,   p=0.7
```

!!! info "What Does this mean?"

    - Here we note that all the coefficients for SEX(being male), AGE, and BMI have slight negative coefficients or hazard ratios
    - We also note that no p-value is below 0.05 indicating that any observed effect on our hazard (death) could be due to chance
    - We also note that sex has the most significant p-value and has the largest magnitude of all the coefficients

## Assumptions

Before accepting the results of this model we should discuss the assumptions of the Cox Proportional-Hazards Model

!!! warning "Assumptions of the Cox Proportional-Hazards Model"

    - the hazard curves different groups should be proportional and not cross each other
    - There are no outliers in our data
    - The relationship between the log hazard and the covariates must be linear
    
**Testing Proportionality in Hazards**

```R
# test proportional hazards assumption
cox_test <- cox.zph(cox)
ggcoxzph(cox_test)
```

![](images/test_prop_cox.png)

!!! info "What does this mean?"

    - Here we see that the global test is above 0.05, indicating we have not violated the proportional hazards assumption

**Testing for Outliers**

```R
# test the outlier assumption
ggcoxdiagnostics(cox,
                 type = "dfbeta",
                 linear.predictions = FALSE, 
                 ggtheme = theme_bw())
```

![](images/test_outlier_cox.png)

!!! info "What does this mean?"

    - Here we see that there are a few outliers as the patterns of the residuals, especially for BMI, are not spread equally around 0. 
    
**Testing for Linearity**

NEEDS MORE EXPLANATION HERE

```R
# test for non-linearity for AGE
ggcoxfunctional(Surv(PATH_DIAG_TO_DEATH_DAYS, status) ~ AGE + log(AGE) + sqrt(AGE),
                data = meta)
```

![](images/test_linear_cox.png)


## References

1. [Cox Proportional-Hazards Model](http://www.sthda.com/english/wiki/cox-proportional-hazards-model)
2. [Regression Models and Life-Tables](https://www.jstor.org/stable/2985181#metadata_info_tab_contents)
3. [Survival Analysis Part II: Multivariate data analysis – an introduction to concepts and methods](https://www.nature.com/articles/6601119)
