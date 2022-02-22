
---
title: "Tests for Partially Matched Samples"
author: "Jonathan Gorman"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{PMGorman-vig}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




```r
library(PMGorman)
```


# Setup

The package 'PMGorman' is not currently available on any repositories, so the source should be downloaded and imported into the user's R environment. Installation can be performed as follows:
```
> install.packages('file_absolute_pathname', repos=NULL, type='source')

```
where the absolute pathname is the file location on the user's system after the download finishes. More information about `install.packages()` can be found by entering `?install.packages` into the console.

# Background

To understand partially matched pairs, first consider the conventional matched pairs experimental design. Typically, pairs of similar subjects are exposed to a different treatments. For example, a set of several pairs of twin brothers --  who are nearly identical and whose anatomical traits are likely correlated -- are selected for a matched pairs experiment. One twin brother for each pair is randomly chosen to undergo treatment, while the other is set as a control. Measurements are taken before and after treatment, and the data are organized into a table.

Suppose one twin neglects report their results to the researchers or is unable to report for some tragic reason. This value will be marked as 'NA' in the data table. Consequently, in terms of the table, the data for those twins is no longer paired. If this occurs at least once in the dataset, the total set of pairs is considered \italic{partially matched}. This same phenomenon can occur is the data are lost.

|   1  |   2  |   3  |  ... |  n  |
|:----:|:----:|:----:|:----:|:----:
|  14  |  13  |  NA  |  ... |  5  |
|  NA  |  2   |  7   |  ... |  1  |

From a mathematical standpoint, if there are $n$ pairs in the experiment, then there will ideally be $2n$ observations in the table. In a partially matched situation, there will be fewer than $2n$ observations. We can assign specific labels to the three new groups: $n1$ for the number of complete pairs, $n2$ for the lone observations in the first treatment, and $n3$ for the lone values in the second treatment. Special tests, as highlighted and implemented in this package, can then be used to analyze this partitioned data.

# Functions

The functions in this package are as follows:

* Kim et. al's Modified t-test, `kim_mod_t()`
* Looney's Correct z-test, `looney_z_corr()`
* Liptak's Weighted z-test, `liptak_z_w()`
* Lin and Stiver's MLE Test with Heteroscedasticity, `lin_MLE_hetero()`
* Ekbohm's MLE Test with Homoscedasticity, `ekbohm_MLE_homo()`

They all test if the mean difference between the groups is exactly zero, and are capable of one-sided and two-sided alternative hypotheses. The rationale behind each test can be found in [1].

# Function Input

Every function can accept either two vectors or a matrix with dimensions $ 2 \times k$ or $k \times 2$, where $k$ is the number of partially matched pairs. The parameter list is identical for every function:

```
function(X, y=NULL, alt= c('two.sided', 'less', 'greater'))
```

where `X` is either a vector or matrix. `y` is a vector, and if `y` is present, `X` too must be a vector with equal length to `y`. If `X` is a matrix, `y` should be absent. 

The `alt` parameter holds the alternative hypothesis of the test, and is set to `two.sided` by default.

An example of data input is as follows:

```
X <- c(1, 16, 4, 8, NA, 7, 19, NA, 14, NA, 10, 18, NA, 5)
y <- c(3, NA, 6, 11, 7, NA, 10, 9, NA, 14, 3, 18, 20, 4)
```

which is two vectors. Matrix form would resemble:


```
X <- matrix(c(1, 16, 4, 8, NA, 7, 19, NA, 14, NA, 10, 18, NA, 5, 3, NA, 6, 11, 7, NA, 10, 9, NA, 14, 3, 18, 20, 4), nrow=2, byrow=T)
```


# Reduction Cases

* When data are inputted as vectors, they are combined into a $2 \times k$ matrix. 
* When a $k \times 2$ matrix is inputted, it is converted to a $2 \times k$ matrix.
* Data sets with fewer than three missing entries (NA) in either treatment result in omission of the lone values and their corresponding NA. The function will proceed with a matched pair t-test because all the remaining pairs are complete. A warning message indicates this.
* Data sets with fewer than three complete pairs result in omission of the complete pairs. The function will proceed with a two-sample t-test because all the remaining observations are lone. A warning message indicates this.
* Pairs of NA values (NA in both treatments) will be omitted.
* Looney's corrected z-test reduces to a "...paired sample or two-sample Z-test when $n_2=n_3=0$ or $n_1=0$, respectively [1]."

# Troubleshooting

Possible sources of error are:

* Input vectors of unequal length.
* Input a matrix as `y`.
* Input a matrix as `X` AND a vector as `y`.
* Misspelling the alternative hypothesis input string.
* Covariance calculations may be undefined when there are fewer than three missing entries (NA) in either treatment AND there are fewer than three complete pairs.

# Sample Data and Examples

## Sample Data

The following code simulates a bivariate gaussian with various correlations and variances, as well as two sets of means. It can be used to test the functions, but requires the package `mvtnorm`.


```

#INSTALL PACKAGE 'mvtnorm'

#Generate bivariate gaussian sample data

require(mvtnorm)
set.seed(123) # simulation seed


rho <- c(0.2, 0.7) # correlation test values
n <- c(20, 50) # test sample sizes
mu0 <- c(0, 0) # Equal means
mu1 <- c(1, 4) # Different means
sig1 <- c(1, 3) # heteroscedastic test standard deviations
sig2 <- c(1, 1) # homoscedastic test standard deviations



#Four cases, mu0, ((n=20,rho=0.2), (n=20,rho=0.7), (n=50, rho=0.2), (n=50, rho=0.7)), heteroscedastic

hetero_20_0.2_mu0 <- rmvnorm(n[1], mean=mu0, sigma=matrix(c(sig1[1]^2, rho[1]*sig1[1]*sig1[2], rho[1]*sig1[1]*sig1[2], sig1[2]^2), nrow=2, byrow=T))
hetero_50_0.2_mu0 <- rmvnorm(n[2], mean=mu0, sigma=matrix(c(sig1[1]^2, rho[1]*sig1[1]*sig1[2], rho[1]*sig1[1]*sig1[2], sig1[2]^2), nrow=2, byrow=T))
hetero_20_0.7_mu0 <- rmvnorm(n[1], mean=mu0, sigma=matrix(c(sig1[1]^2, rho[2]*sig1[1]*sig1[2], rho[2]*sig1[1]*sig1[2], sig1[2]^2), nrow=2, byrow=T))
hetero_50_0.7_mu0 <- rmvnorm(n[2], mean=mu0, sigma=matrix(c(sig1[1]^2, rho[2]*sig1[1]*sig1[2], rho[2]*sig1[1]*sig1[2], sig1[2]^2), nrow=2, byrow=T))


#Four cases, mu0, ((n=20,rho=0.2), (n=20,rho=0.7), (n=50, rho=0.2), (n=50, rho=0.7)), homoscedastic

homo_20_0.2_mu0 <- rmvnorm(n[1], mean=mu0, sigma=matrix(c(sig2[1]^2, rho[1]*sig2[1]*sig2[2], rho[1]*sig2[1]*sig2[2], sig2[2]^2), nrow=2, byrow=T))
homo_50_0.2_mu0 <- rmvnorm(n[2], mean=mu0, sigma=matrix(c(sig2[1]^2, rho[1]*sig2[1]*sig2[2], rho[1]*sig2[1]*sig2[2], sig2[2]^2), nrow=2, byrow=T))
homo_20_0.7_mu0 <- rmvnorm(n[1], mean=mu0, sigma=matrix(c(sig2[1]^2, rho[2]*sig2[1]*sig2[2], rho[2]*sig2[1]*sig2[2], sig2[2]^2), nrow=2, byrow=T))
homo_50_0.7_mu0 <- rmvnorm(n[2], mean=mu0, sigma=matrix(c(sig2[1]^2, rho[2]*sig2[1]*sig2[2], rho[2]*sig2[1]*sig2[2], sig2[2]^2), nrow=2, byrow=T))

#Four cases, mu1, ((n=20,rho=0.2), (n=20,rho=0.7), (n=50, rho=0.2), (n=50, rho=0.7)), heteroscedastic

hetero_20_0.2_mu1 <- rmvnorm(n[1], mean=mu1, sigma=matrix(c(sig1[1]^2, rho[1]*sig1[1]*sig1[2], rho[1]*sig1[1]*sig1[2], sig1[2]^2), nrow=2, byrow=T))
hetero_50_0.2_mu1 <- rmvnorm(n[2], mean=mu1, sigma=matrix(c(sig1[1]^2, rho[1]*sig1[1]*sig1[2], rho[1]*sig1[1]*sig1[2], sig1[2]^2), nrow=2, byrow=T))
hetero_20_0.7_mu1 <- rmvnorm(n[1], mean=mu1, sigma=matrix(c(sig1[1]^2, rho[2]*sig1[1]*sig1[2], rho[2]*sig1[1]*sig1[2], sig1[2]^2), nrow=2, byrow=T))
hetero_50_0.7_mu1 <- rmvnorm(n[2], mean=mu1, sigma=matrix(c(sig1[1]^2, rho[2]*sig1[1]*sig1[2], rho[2]*sig1[1]*sig1[2], sig1[2]^2), nrow=2, byrow=T))


#Four cases, mu1, ((n=20,rho=0.2), (n=20,rho=0.7), (n=50, rho=0.2), (n=50, rho=0.7)), homoscedastic

homo_20_0.2_mu1 <- rmvnorm(n[1], mean=mu1, sigma=matrix(c(sig2[1]^2, rho[1]*sig2[1]*sig2[2], rho[1]*sig2[1]*sig2[2], sig2[2]^2), nrow=2, byrow=T))
homo_50_0.2_mu1 <- rmvnorm(n[2], mean=mu1, sigma=matrix(c(sig2[1]^2, rho[1]*sig2[1]*sig2[2], rho[1]*sig2[1]*sig2[2], sig2[2]^2), nrow=2, byrow=T))
homo_20_0.7_mu1 <- rmvnorm(n[1], mean=mu1, sigma=matrix(c(sig2[1]^2, rho[2]*sig2[1]*sig2[2], rho[2]*sig2[1]*sig2[2], sig2[2]^2), nrow=2, byrow=T))
homo_50_0.7_mu1 <- rmvnorm(n[2], mean=mu1, sigma=matrix(c(sig2[1]^2, rho[2]*sig2[1]*sig2[2], rho[2]*sig2[1]*sig2[2], sig2[2]^2), nrow=2, byrow=T))


#Generate indices to make NA

NA_index50 <- sample(1:50, 20) # 20 NA indices generated
NA_index20 <- sample(1:20, 8) # 8 NA indices generated


#Use generated indices to make random entries NA

#mu0

hetero_20_0.2_mu0[,1][NA_index20[1:4]] <- NA
hetero_20_0.2_mu0[,2][NA_index20[5:8]] <- NA
hetero_50_0.2_mu0[,1][NA_index50[1:10]] <- NA
hetero_50_0.2_mu0[,2][NA_index50[11:20]] <- NA
hetero_20_0.7_mu0[,1][NA_index20[1:4]] <- NA
hetero_20_0.7_mu0[,2][NA_index20[5:8]] <- NA
hetero_50_0.7_mu0[,1][NA_index50[1:10]] <- NA
hetero_50_0.7_mu0[,2][NA_index50[11:20]] <- NA

homo_20_0.2_mu0[,1][NA_index20[1:4]] <- NA
homo_20_0.2_mu0[,2][NA_index20[5:8]] <- NA
homo_50_0.2_mu0[,1][NA_index50[1:10]] <- NA
homo_50_0.2_mu0[,2][NA_index50[11:20]] <- NA
homo_20_0.7_mu0[,1][NA_index20[1:4]] <- NA
homo_20_0.7_mu0[,2][NA_index20[5:8]] <- NA
homo_50_0.7_mu0[,1][NA_index50[1:10]] <- NA
homo_50_0.7_mu0[,2][NA_index50[11:20]] <- NA

#mu1

hetero_20_0.2_mu1[,1][NA_index20[1:4]] <- NA
hetero_20_0.2_mu1[,2][NA_index20[5:8]] <- NA
hetero_50_0.2_mu1[,1][NA_index50[1:10]] <- NA
hetero_50_0.2_mu1[,2][NA_index50[11:20]] <- NA
hetero_20_0.7_mu1[,1][NA_index20[1:4]] <- NA
hetero_20_0.7_mu1[,2][NA_index20[5:8]] <- NA
hetero_50_0.7_mu1[,1][NA_index50[1:10]] <- NA
hetero_50_0.7_mu1[,2][NA_index50[11:20]] <- NA

homo_20_0.2_mu1[,1][NA_index20[1:4]] <- NA
homo_20_0.2_mu1[,2][NA_index20[5:8]] <- NA
homo_50_0.2_mu1[,1][NA_index50[1:10]] <- NA
homo_50_0.2_mu1[,2][NA_index50[11:20]] <- NA
homo_20_0.7_mu1[,1][NA_index20[1:4]] <- NA
homo_20_0.7_mu1[,2][NA_index20[5:8]] <- NA
homo_50_0.7_mu1[,1][NA_index50[1:10]] <- NA
homo_50_0.7_mu1[,2][NA_index50[11:20]] <- NA

```

Now 16 different partially matched datasets are available with varied parameters as specified above. They are named `(hetero/homo)_(n)_(rho)_(mu)`. The precise names are as follows:

```
hetero_20_0.2_mu0, hetero_50_0.2_mu0, hetero_20_0.7_mu0, hetero_50_0.7_mu0
homo_20_0.2_mu0, homo_50_0.2_mu0, homo_20_0.7_mu0, homo_50_0.7_mu0

hetero_20_0.2_mu1 hetero_50_0.2_mu1, hetero_20_0.7_mu1, hetero_50_0.7_mu1
homo_20_0.2_mu1, homo_50_0.2_mu1, homo_20_0.7_mu1, homo_50_0.7_mu1
```

# Examples

The data sets above can be split into two vectors, or kept as a matrix. The functions will recognize either.

```
kim_mod_t(homo_20_0.2_mu1, alt='less') # alternative set to 'less'

looney_z_corr(hetero_50_0.7_mu0) # Default alternative 'two.sided'

liptak_z_w(homo_50_0.2_mu1[,1], homo_50_0.2_mu1[,2], alt='less') # Two vectors

lin_MLE_hetero(hetero_20_0.7_mu1) 

ekbohm_MLE_homo(homo_50_0.7_mu0[,1], homo_50_0.7_mu0[,2], alt='greater') # Two vectors

```

# References

[1] Kuan, P.F. and Huang, B. (2013), A simple and robust method for partially matched samples using the p‐values pooling approach. Statist. Med., 32: 3247-3259. https://doi.org/10.1002/sim.5758



