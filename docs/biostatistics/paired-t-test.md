## Paired T-Test

We have seen how to compare our sample mean to a theory of the underlying population. Here we will ask how do we compare the means of *paired* observations?
We can use the paired t-test to answer any one of the following questions:

- **Does the mean difference, $m$, equal to 0?**
    - $H_0: m = 0$
    - $H_a: m \neq 0$
- **Is the mean difference, $m$, less than 0?**
    - $H_0: m \le 0$
    - $H_a: m > 0$
- **Is the mean difference, $m$, greater than 0?**
    - $H_0: m \ge 0$
    - $H_a: m < 0$

!!! tip
    When we ask if the mean difference is equal to 0 we are conducting a **two-sided test**. When we ask if the sample mean is less than 
    or greater than 0, we are conducting a **one-sided test**.

## Test Statistic

Our paired t-test statistic can be calculated by:

$$t = \frac{m}{\sigma / \sqrt{n}} $$

$$d.f. = n - 1$$

!!! example "Explanation of Terms"
    - $m$ : mean difference
    - $\sigma$ : standard deviation of our sample
    - $n$ : sample size
    - $d.f.$ : degrees of freedom

## Normal Distribution

Just like in the [one-sample t-test topic note](one-t-test.md) we will be comparing our test statistic to the normal distribution with a 
probability density function of:

$$f(x) = \frac{1}{(\sigma\sqrt{2 \pi})} e^{-(\frac{(x - \mu)^2}{2 \sigma^2})}$$

!!! example "Explanation of Terms"
    - $\sigma$ : standard deviation
    - $\mu$ : mean

## Confidence Interval

Similarily to our [one-sample t-test topic note](one-t-test.md) our confidence interval will be defined as:

$$\mu \pm t \frac{\sigma}{\sqrt{n}}$$

!!! example "Explanation of Terms"
    - $\mu$ : sample mean difference
    - $t$ : test statistic for an $\alpha$ of 0.05
    - $\sigma$ : standard deviation of the differences
    - $n$ :  sample size

## Running the Paired T-Test

Seeing as our glioblastoma data does not have a paired observation (essentially a before-after comparison), we are going to need to generate some 
dummy data to apply our paired-test to:

```R
# paired t-test
library(tidyverse)
# create our dummy data
df <- data.frame(
  cholesterol = c(sample(70:120,20),sample(110:220,20)),
  group = c(rep("before",20),rep("after",20))
)

# what is mean of each group
df %>%
  group_by(group) %>%
  summarise(
    n=n(),
    standard_deviation=sd(cholesterol),
    mean=mean(cholesterol)
  )
```

````
# A tibble: 2 × 4
  group      n standard_deviation  mean
  <chr>  <int>              <dbl> <dbl>
1 after     20               32.4 160. 
2 before    20               14.2  94.8
````

Here we note that the before group has much lower cholesterol levels than the after group. Let's now apply the paired t-test to determine the significance
of this difference. We will do this by asking is the difference between groups equal 0?

```R
# run the paired t-test to determine if the 
# difference between groups is equal to 0
t.test(
  cholesterol ~ group, 
  data = df, 
  paired = TRUE)
```

```
	Paired t-test

data:  cholesterol by group
{==t = 8.6603, df = 19, p-value = 5.061e-08==}
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
{== 49.59408 81.20592==}
sample estimates:
mean of the differences 
                   {==65.4==}
```

!!! info "Explanation"
    - our test statistic is `8.6603`
    - The pvalue is ` 5.061e-08`
    - our alternative hypothesis is that the true difference in means is not equal to 0
    - our mean of the differences is `65.4`
    - the 95% confidence interval for our mean of differences is `49.59408` to `81.20592`
    - So we see that we have enough evidence to reject the null hypothesis that the true difference in means is equal to 0
    
## Assumptions

So now that we have conducted our test, we should assess the test's assumptions:

- the values are paired
- the differences of the pairs are normally distributed

The values in our dummy data are intended to be paired. Meaning that the first observation in the before group is the same patient as the first observation
in the after group. Now to check if the differences of the pairs are normally distributed we can use the Shapiro-Wilk test like we did in the 
[one-sample t-test topic note](one-t-test.md):

```R
# differences betweeen the paired observations
difference <- df$cholesterol[df$group == "before"] - df$cholesterol[df$group == "after"]

# Shapiro-Wilk test to determine are they normally distributed
shapiro.test(difference) 
```

```
	Shapiro-Wilk normality test

data:  difference
W = 0.94613, {==p-value = 0.3121==}
```

Given our p-value is **above** 0.05 we do not have enough evidence to reject the null hypothesis of the Shapiro-Wilk Test; that the data are normally distributed. In other words, if the p-value is above 0.05 your data are normally distributed.

## Non-Parametric Alternative

A non-parametric test is often used when either the assumptions about the distribution are not met. Additionally, these tests do not depend on the *parameter* they are assessing. Here, if the assumptions above are not met we can use the non-parametric equivalent, the Wilcoxon signed rank test:

```R
# run the non-parametric alternative to the paired
# t-test the Wilcoxon signed rank test
wilcox.test(cholesterol ~ group, 
            data = df, 
            paired = TRUE,
            alternative = "two.sided")
```

```
	Wilcoxon signed rank exact test

data:  cholesterol by group
V = 210, p-value = 1.907e-06
alternative hypothesis: true location shift is not equal to 0
```

## References

1. [BIOL 202 - Paired T-Test](https://ubco-biology.github.io/BIOL202/paired-t-test.html)
2. [Paired Samples T-test in R](http://www.sthda.com/english/wiki/paired-samples-t-test-in-r)
3. [Normal Distribution](https://en.wikipedia.org/wiki/Normal_distribution)
4. [Paired Samples Wilcoxon Test in R](http://www.sthda.com/english/wiki/paired-samples-wilcoxon-test-in-r)
5. [Nonparametric statistics](https://en.wikipedia.org/wiki/Nonparametric_statistics)
