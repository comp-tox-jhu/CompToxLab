## Clustering

Clustering, like the name implies, is a way to group data points to find patterns. Clustering is type of unsupervised learning - where patterns are discovered without labeled data. 

## Distance Metrics

Before generating clusters we need to calculate the distance between observations:

![](images/calc_dist_mat.png)

To do so we have a few different options for distance metrics:

**Euclidean Distance**

$$d_{euc}(x,y) = \sqrt{\sum_{i=1}^n{(x_i - y_i)^2}}$$

!!! example "Explanation of Terms"
    - $x$ variable x
    - $y$ variable y
    - $n$ number of observations
    
**Manhattan Distance**

$$d_{man}(x,y) = \sum_{i=1}^n{|(x_i - y_i)|}$$

!!! example "Explanation of Terms"
    - $x$ variable x
    - $y$ variable y
    - $n$ number of observations
  
**Eisen Cosine Correlation Distance**

$$d_{eis}(x,y) = 1 - \frac{|\sum_{i=1}^n{x_iy_i}|}{\sqrt{\sum_{i=1}^n {x_i^2}\sum_{i=1}^n {y_i^2}}}$$

!!! example "Explanation of Terms"
    - $x$ variable x
    - $y$ variable y
    - $n$ number of observations
    
**Pearson Correlation Distance**  

$$d_{pearson}(x,y) = 1 - \frac{\sum_{i=1}^n{(x - \mu_x)(y - \mu_y)}}{\sqrt{\sum_{i=1}^n{(x - \mu_x)^2} \sum_{i=1}^n{(y - \mu_y)^2}}} $$

!!! example "Explanation of Terms"
    - $x$ variable x
    - $y$ variable y
    - $\mu_x$ mean of variable x
    - $\mu_y$ mean of variable y
    - $n$ number of observations

**Spearman Correlation Distance**  

$$d_{spearman}(x,y) = 1 - \frac{\sum{(x\prime - \mu_{x\prime} )(y\prime  - \mu_{y\prime} )}}{\sqrt{\sum{(x\prime  - \mu_{x\prime} )^2} \sum{(y\prime  - \mu_{y\prime} )^2}}}$$

!!! example "Explanation of Terms"
    - $x$ variable x
    - $y$ variable y
    - $\mu_x$ mean of variable x
    - $\mu_y$ mean of variable y
    - $n$ number of observations

## Distance Metrics In R

Let's try creating a distance matrix!

```R
# load the libraries
.libPaths(c("/cluster/tufts/hpc/tools/R/4.0.0"))
library(tidyverse)
library(factoextra)

# load our counts data
counts <- read.csv(
  file="data/gbm_cptac_2021/data_mrna_seq_fpkm.txt",
  header = T,
  sep = "\t")

# make the genes our rownames
rownames(counts) <- make.names(counts$Hugo_Symbol,unique = TRUE)

# remove the gene symbol column
counts <- counts %>%
  select(-c(Hugo_Symbol)) 

# log2 transform our data 
# transpose our data so that our patients are rows
counts <- t(log2(counts + 1))

# Change NA counts to 0
counts[!is.finite(counts)] <- 0

# generate correlation distance matrix
dist <- get_dist(counts,method = "pearson")

# plot correlation distance matrix
fviz_dist(dist) +
  theme(axis.text = element_text(size = 3)) +
  labs(
    title = "Pearson Correlation Distances Between Samples",
    fill = "Pearson Correlation"
  )
```

![](images/sample_corr_mat.png)

## References

1. [Clustering Distance Measures](https://www.datanovia.com/en/lessons/clustering-distance-measures/)
2. [K-Means Clustering in R: Algorithm and Practical Examples](https://www.datanovia.com/en/lessons/k-means-clustering-in-r-algorith-and-practical-examples/)
3. [Agglomerative Hierarchical Clustering](https://www.datanovia.com/en/lessons/agglomerative-hierarchical-clustering/)
4. [Distance Method Formulas](https://www.jmp.com/support/help/14/distance-method-formulas.shtml#177809%C2%A0)
