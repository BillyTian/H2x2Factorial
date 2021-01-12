
<!-- README.md is generated from README.Rmd. Please edit that file -->

# H2x2Factorial

<!-- badges: start -->
<!-- badges: end -->

H2x2Factorial encapsulates the sample size methods in “Sample size
calculation in hierarchical 2x2 factorial trials with unequal cluster
sizes”. It contains one function for the power calculation or the sample
size estimation, one function for sample size calculations in a table
format, and one function for plotting illustrative graph of sample size
requirements.

## Installation

The development version of the packages can be installed from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("BillyTian/H2x2Factorial")
```

## Example

Here shows some basic examples:

``` r
library(H2x2Factorial)

H2x2Factorial(n.input=10, delta_x=0.2, delta_z=0.1, test="joint", correction=T, seed.mix=123456, CV=0.38)

table.H2x2Factorial(delta_x=0.2, delta_z=0.1, m_bar=c(10,50,100), CV=c(0, 0.3, 0.5), rho=c(0.01, 0.1), test="cluster")

plot.H2x2Factorial(power=0.9, test="cluster")
```
