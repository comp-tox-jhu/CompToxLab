# Statistical Tests

### 1 Dependent Variable/ 0 Independent Variables

| statistical test | description | example in R |
:---------------------|:---------------------|:---------------------|
| **one sample t-test**| test whether a sample mean significantly differs from a hypothesized value (assumes your vairable is normally distributed)| ```t.test(data, mu = meanToTest)``` |
|**one sample median test** | test whether a sample median differs significantly from a hypothesized value (does not assume your variable is normally distributed) | ```wilcox.test(data, mu = medianToTest, alternative = "two.sided")```|
| **binomial test** | test whether the proportion of successes on a two-level categorical dependent variable significantly differs from a hypothesized value - | ```binom.test(numberOfActualSuccesses, numberOfTrialsYouDo, probabilityOfSuccessToTest)``` |
| **chi-square goodness of fit** | test whether the observed proportions for a categorical variable differ from hypothesized proportions | ```chisq.test(vectorOfCounts, hypothesizedProportions)``` |

### 1 Dependent Variable/ 1 Independent Variables with 2 Levels

| statistical test | description | example in R |
:---------------------|:---------------------|:---------------------|
| **two independent samples t-test** | compare the means of a normally distributed dependent variable for two independent groups | ```t.test(data1, data2)``` |
| **wilcoxon-mann-whitney test** | non-parametric analog to the independent samples t-test and can be used when the dependent variable is not normally distributed | ```wilcox.test(data1, data2, alternative = "two.sided")``` |
| **chi-square test** | tests for a relationship between two categorical variables (makes the assumption that each cell has at least 5 when you split by table!)| ```chisq.test(table(variable1,variable2))```|
| **fisher’s exact test** | tests for a relationship between two categorical variables, but can be used when cells have counts less than 5 | ```fisher.test(table(variable1,variable2))``` |

### 1 Dependent Variable/ 1 Independent Variables with 2 or More Levels

| statistical test | description | example in R |
:---------------------|:---------------------|:---------------------|
| **one-way analysis of variance(ANOVA)** | test for differences in the means of the dependent variable broken down by the levels of the independent variable - assumes dependent variable is normally distributed, variances for each of the groups are the same | ```aov(numericDependentVariable ~ categoricalIndependentVariable, data = dataFrameWithBothVariables)```|
| **analysis of co-variance(ANCOVA)** | test for differences in the means of the dependent variable broken down by the levels of two independent variable - assumes dependent variable is normally distributed, variances for each of the groups are the same | ```aov(numericDependentVariable ~ categoricalIndependentVariable1+categoricalIndependentVariable2, data = dataFrameWithBothVariables)```|
| **shapiro-wilk test** | tests for normality of a variable - small p-value = not normally distributed / big p-value = is normally distributed | ```shapiro.test(variableToTest)``` |
| **levene test** | tests for differences in variances among groups - small p-value = there are differences among groups / big p-value = no differences among groups |```levene.test(numericDependentVariable ~ categoricalIndependentVariable, data = dataFrameWithBothVariables)``` |
| **kruskal wallis test**| test for differences between a dependent variable broken down by the levels of the independent variable - non-parametric alternative to one way anova, does not assume normality or equal variances | ```kruskal.test(numericDependentVariable ~ categoricalIndependentVariable, data = dataFrameWithBothVariables)``` |

### 1 Dependent Variable/ 1 Independent Variables with 2 Paired Levels 

| statistical test | description | example in R |
:---------------------|:---------------------|:---------------------|
| **paired t-test** | compare the means of a normally distributed dependent variable for two dependent/related groups | ```t.test(data1, data2, paired = TRUE, alternative = "two.sided")```|
| **wilcoxon signed rank sum test**| non-parametric version of paired t-test - does not assume normality | ```wilcox.test(data1, data2, paired = TRUE, alternative = "two.sided")``` |
| **mcnemar test**| tests for differences in proportions between paired data - like before and after some event | ```mcnemar.test(table(pairedVariable1,pairedVariable2))```|

### 1 Dependent Variable/ 1 Independent Variables with 2 or More Paired Levels

| statistical test | description | example in R |
:---------------------|:---------------------|:---------------------|
| **one-way repeated measures analysis of variance(ANOVA)**| tests for differences between the means of three or more groups where the same subjects show up in each group | ```aov(numericDependentVariable~factor(categoricalVariableThatChanges)+Error(factor(subject)), data = dataFrameWithAllVariables)```|
| **friedman test**|non-parametric alternative to the one-way repeated measures ANOVA | ```friedman.test(y=numericDependentVariable, groups=categoricalVariableThatChanges, blocks=subjects)```|

### 1 Dependent Variable/ 1 Numeric Independent Variables

| statistical test | description | example in R |
:---------------------|:---------------------|:---------------------|
| **correlation**| tests for a relationship between normally distributed variables | ```cor.test(variable1, variable2, method="pearson")``` |
| **non-parametric correlation** | tests for a relationship between non-normally distributed variables | ```cor.test(variable1, variable2, method="spearman")```|

### 2 or More Dependent Variables/ 2 Independent Variables with 2 or More Levels

| statistical test | description | example in R |
:---------------------|:---------------------|:---------------------|
| **one-way multivariate analysis of variance(MANOVA)** | assess how two or more dependent variables are affected by a categorical variable | ```manova(cbind(numericDependentVariable1, numericDependentVariable2) ~ independentCategoricalVariable, data = dataframeWithAllTheVariables)``` |


!!! abstract "References"

    1. [UCLA Statistical Methods and Analysis](https://stats.oarc.ucla.edu/stata/whatstat/what-statistical-analysis-should-i-usestatistical-analyses-using-stata/)
    2. [Statology](https://www.statology.org/)
    3. [STHDA](http://www.sthda.com/)
