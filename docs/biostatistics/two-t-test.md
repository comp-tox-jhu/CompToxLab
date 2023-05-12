## Two-Sample T-Test

In the [paired t-test topic note](paired-t-test.md) we examined the difference in means between paired data. To compare the difference in means between
two unpaired groups we will use the two-sample t-test. The two-sample t-test can be used to ask the following:

- **Does the mean of group 1, $\mu_1$, equal the mean of our group 2, $\mu_2$?**
    - $H_0: \mu_1 = \mu_2$
    - $H_a: \mu_1 \neq \mu_2$
- **Is the mean of group 1, $\mu_1$, less than the mean of group 2, $\mu_2$?**
    - $H_0: \mu_1 \le \mu_2$
    - $H_a: \mu_1 > \mu_2$
- **Is the mean of group 1, $\mu_1$, greater than the mean of group 2, $\mu_2$?**
    - $H_0: \mu_1 \ge \mu_2$
    - $H_a: \mu_1 < \mu_2$

!!! tip
    When we ask if mean of group 1 is equal to the mean of group 2 we are conducting a **two-sided test**. When we ask if the mean of group 1 is less than or greater than the mean of group 2, we are conducting a **one-sided test**.

## Test Statistic

Our one-sample t-test statistic can be calculated by:

$$t = \frac{\mu_1 - \mu_2}{\sqrt{\frac{\sigma_1}{n_1} + \frac{\sigma_2}{n_2}}} $$

$$d.f. = n_1 + n_2 - 2$$

!!! example "Explanation of Terms"
    - $\mu_1$ : mean of group 1
    - $\mu_2$ : mean of group 2
    - $\sigma_1$ : standard deviation of group 1
    - $\sigma_2$ : standard deviation of group 2
    - $n_1$ : size of group 1
    - $n_2$ : size of group 2
    - $d.f.$ : degrees of freedom

## Normal Distribution

Just like in the one-sample t-test topic note we will be comparing our test statistic to the normal distribution with a probability density function of:

$$f(x) = \frac{1}{(\sigma\sqrt{2 \pi})} e^{-(\frac{(x - \mu)^2}{2 \sigma^2})}$$

!!! example "Explanation of Terms"
    - $\sigma$ : standard deviation
    - $\mu$ : mean

## Confidence Interval

Just like our proportion tests, we also have a confidence interval around our sample parameter, in this case the sample mean. So for a test statistic, $t$, at an $\alpha$ level of 0.05, our confidence interval would be:

$$(\mu_1 - \mu_2) \pm t * \sqrt{\frac{\sigma_1}{n_1} + \frac{\sigma_2}{n_2}}$$

!!! example "Explanation of Terms"
    - $\mu_1$ : mean of group 1
    - $\mu_2$ : mean of group 2
    - $t$ : test statistic for an $\alpha$ of 0.05
    - $\sigma_1$ : standard deviation of group 1
    - $\sigma_2$ : standard deviation of group 2
    - $n_1$ : size of group 1
    - $n_2$ : size of group 2

## Running the Two-Sample T-Test

Putting this all together let's ask: Is mean age of males equal to the mean age of females in our glioblastoma data? We will start by manually calculating the mean age and standard deviation of males and females:

```R
# what is mean of each group
meta %>%
  group_by(SEX) %>%
  summarise(
    n=n(),
    standard_deviation=sd(AGE),
    mean=mean(AGE)
  )
```

```
# A tibble: 2 × 4
  SEX        n standard_deviation  mean
  <chr>  <int>              <dbl> <dbl>
1 Female    44               13.7  57.9
2 Male      55               11.6  57.9
```

Remarkably we notice that the mean age of males and the mean age of females are the same! We will now use the two-sample t-test to ask if the mean age of males equals the mean age of females:

```R
# run the two-sample t-test to determine if the
# mean age of males is equal to the mean
# age of females
t.test(
  AGE ~ SEX, 
  data = meta, 
  method = "two.sided",
  var.equal = FALSE)
```

```
	Welch Two Sample t-test

data:  AGE by SEX
{==t = 0.014062, df = 84.559, p-value = 0.9888==}
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
{== -5.105621  5.178348 ==}
sample estimates:
mean in group Female   mean in group Male 
{==            57.90909             57.87273 ==}
```

!!! info "Explanation"
    - our test statistic is `0.014062`
    - The pvalue is below `0.9888`
    - our alternative hypothesis is that the true difference in means is not equal to 0
    - the mean female is is `57.90909` and the mean male age is `57.87273`
    - the 95% confidence interval for our difference in means is `-5.105621` to `5.178348`
    - So we see that we do not have enough evidence to reject the null hypothesis that the true difference in means is equal to 0
    
## Assumptions

So now that we have conducted our test, we should assess the test's assumptions:

- the values are independent of one another
- the variances in each group are equal
- the data in each group normally distributed

Here we note that our age values should be independent of one another (the age of one patient should not affect the age of another). Now for the assumption that variances of each group are equal - we avoided this by setting the `var.equal` argument equal to false. If we were to set it to true, we could check our variances with the F-Test:

```R
# use the f test to determine if the variance in 
# group 1 is the same as the variance in group 2
var.test(
  AGE ~ SEX, 
  data = meta
)
```

```
	F test to compare two variances

data:  AGE by SEX
F = 1.3849, num df = 43, denom df = 54, {==p-value = 0.2559==}
alternative hypothesis: true ratio of variances is not equal to 1
95 percent confidence interval:
 0.7879316 2.4804317
sample estimates:
ratio of variances 
          1.384904 
```

Here we will just note that the p-value is above 0.05 and thus there isn't a significant difference in the variance of group 1 and the variance of group 2. Now to check if the data in each group are normally distributed we will use the Shapiro-Wilk test:

```R
# check the normality of male ages and female ages
shapiro.test(meta$AGE[meta$SEX == "Male"])
shapiro.test(meta$AGE[meta$SEX == "Female"])

```

```
	Shapiro-Wilk normality test

data:  meta$AGE[meta$SEX == "Male"]
W = 0.98809, {==p-value = 0.8604==}

	Shapiro-Wilk normality test

data:  meta$AGE[meta$SEX == "Female"]
W = 0.95578, {==p-value = 0.09048==}
```

Given our p-value is **above** 0.05 we do not have enough evidence to reject the null hypothesis of the Shapiro-Wilk Test; that the data are normally distributed. In other words, if the p-value is above 0.05 your data are normally distributed. However, it should be noted that the Shapiro-Wilk test p-value for female ages is rather close to 0.05 so there does seem to be some skey in female ages. 

## Non-Parametric Alternative

A non-parametric test is often used when either the assumptions about the distribution are not met. Additionally, these tests do not depend on the parameter they are assessing. Here, if the assumptions above are not met we can use the non-parametric equivalent, the Wilcoxon signed rank test:

```R
# run the non-parametric alternative to the unpaired
# t-test the Wilcoxon signed rank test
wilcox.test(AGE ~ SEX, 
            data = meta, 
            alternative = "two.sided")
```

```
	Wilcoxon rank sum test with continuity correction

data:  AGE by SEX
W = 1278.5, p-value = 0.6318
alternative hypothesis: true location shift is not equal to 0
```

## References

1. [BIOL 202 - Two-Sample T-Test](https://ubco-biology.github.io/BIOL202/twosamp_ttest.html)
2. [Unpaired Two-Samples T-test in R](http://www.sthda.com/english/wiki/unpaired-two-samples-t-test-in-r)
3. [Normal Distribution](https://en.wikipedia.org/wiki/Normal_distribution)
4. [2 sample test of mean](https://statmagic.info/Content/Help-Content/two-sample-mean.html)
5. [Unpaired Two-Samples Wilcoxon Test in R](http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r)
6. [Nonparametric statistics](https://en.wikipedia.org/wiki/Nonparametric_statistics)
