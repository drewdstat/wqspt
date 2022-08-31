
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Weighted Quantile Sum (WQS) Permutation Test

<!-- badges: start -->
<!-- badges: end -->

The goal of wqspt is to implement a permutation test method for the
weighted quantile sum (WQS) regression.

Weighted quantile sum regression is a statistical technique to evaluate
the effect of complex exposure mixtures on an outcome (Carrico et
al. 2015). It is a single-index method which estimates a combined
mixture sum effect as well as weights determining each individual
mixture component’s contributions to the sum effect. However, the model
features a statistical power and Type I error (i.e., false positive)
rate tradeoff, as there is a machine learning step to determine the
weights that optimize the linear model fit. If the full data is used to
estimate both the mixture component weights and the regression
coefficients, there is high power but also a high false positive rate
since coefficient p-values are calculated for a weighted mixture
independent variable calculated using weights that have already been
optimized to find a large effect.

This package provides an alternative method based on a permutation test
that should reliably allow for both high power and low false positive
rate when utilizing the WQSr. The permutation test is a method of
obtaining a p-value by simulating the null distribution through
permutations of the data. The permutation test algorithm is described
more in detail and validated in [Day et al. 2022](https://ehp.niehs.nih.gov/doi/10.1289/EHP10570). 
The version of this permutation test used for a continuous outcome 
variable has been applied in [Loftus et al. 2021](https://www.sciencedirect.com/science/article/pii/S0160412021000337), 
[Day et al. 2021](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9291724/), 
[Wallace et al. 2022](https://www.sciencedirect.com/science/article/pii/S0160412021006644), 
[Barrett et al. 2022](https://www.sciencedirect.com/science/article/pii/S0160412022000034),
and [Freije et al. 2022](https://www.sciencedirect.com/science/article/pii/S0160412022001726).
Another version of the permutation test adapted for logistic WQSr with a binary outcome
variable is applied in [Loftus et al. 2022](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4102800).

## Installation

You can install the development version of wqspt from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("drewdstat/wqspt")
```

## Usage

Here is a [brief tutorial
vignette](http://htmlpreview.github.io/?https://github.com/drewdstat/wqspt/blob/main/vignettes/introduction.html)
on how to use the wqspt package.
