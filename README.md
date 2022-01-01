
<!-- README.md is generated from README.Rmd. Please edit that file -->

# H2x2Factorial

<!-- badges: start -->
<!-- badges: end -->

This `H2x2Factorial` package implements the sample size methods for
hierarchical 2x2 factorial trials with unequal cluster sizes. The sample
size calculations support two types of treatment effect estimands and
five types of hypothesis tests based on the two measures. The two
estimands are named as the controlled effect and the natural effect, as
formally defined in Tian et al. (Stat Med accepted); The hypotheses include
(A1) test for the cluster-level controlled effect, (A2) test for the
individual-level controlled effect, (B1) test for the cluster-level
natural effect, (B2) test for the individual-level natural effect, (C)
interaction test for the two treatments, (D1) joint test for the two
controlled treatment effects, (D2) joint test for the two natural
treatment effects, (E1) intersection-union test for the two controlled
treatment effects, (E2) intersection-union test for the two natural
treatment effects. Finite-sample considerations are included for the
tests involving either cluster-level treatment effect, due to the degree
of freedom issues. Three functions are currently contained for
predicting the power or sample size based on given design parameters as
well as delivering illustrative tables or line plots. Specifically, the
`calc.H2x2Factorial` function calculates required number of clusters for
a specific test to achieve a given power, or predicts the actual power
given specified sample size resources, with or without finite-sample
considerations. The `table.H2x2Factorial` function creates a data frame
to show a series of sample size predictions by providing varying mean
cluster sizes, intraclass correlation coefficients, or coefficient of
variations of cluster sizes (CV). The `graph.H2x2Factorial` function
plots sample size requirements under different CV in the form of the
combinations of mean cluster sizes and number of clusters. All of the
hypothesis tests and sample size methodologies are formalized in “Sample
size calculation in hierarchical 2x2 factorial trials with unequal
cluster sizes” (under review).

## Installation

The released version of H2x2Factorial can be installed from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("H2x2Factorial")
```

## Example

This is an example for predicting the required number of clusters based
on fixed design parameters:

``` r
library(H2x2Factorial)
#> Warning: package 'H2x2Factorial' was built under R version 4.0.5
example("calc.H2x2Factorial")
#> 
#> c.H22F> #Predict the actual power of a joint test when the number of clusters is 10
#> c.H22F> joint.power <- calc.H2x2Factorial(n_input=10,
#> c.H22F+                                   delta_x=0.2, delta_z=0.1,
#> c.H22F+                                   rho=0.1, CV=0.38,
#> c.H22F+                                   test="joint", correction=TRUE, seed_mix=123456, verbose=FALSE)
#> 
#> c.H22F> print(joint.power)
#> [1] 0.2131
```

This is an example for displaying a series of sample size predictions in
a table format based on varying design parameters:

``` r
example("table.H2x2Factorial")
#> 
#> t.H22F> #Make a result table by providing three mean cluster sizes, three CV, and three ICC
#> t.H22F> table.cluster <- table.H2x2Factorial(delta_x=0.2, delta_z=0.1,
#> t.H22F+                                      m_bar=c(10,50,100), CV=c(0, 0.3, 0.5), rho=c(0.01, 0.1),
#> t.H22F+                                      test="cluster", verbose=FALSE)
#> 
#> t.H22F> table.cluster
#>    m_bar  rho  CV   n predicted power
#> 1     10 0.01 0.0  86       0.8020410
#> 2     10 0.01 0.3  87       0.8036148
#> 3     10 0.01 0.5  88       0.8027978
#> 4     10 0.10 0.0 150       0.8022800
#> 5     10 0.10 0.3 153       0.8011498
#> 6     10 0.10 0.5 160       0.8023522
#> 7     50 0.01 0.0  24       0.8100115
#> 8     50 0.01 0.3  24       0.8021486
#> 9     50 0.01 0.5  25       0.8036072
#> 10    50 0.10 0.0  93       0.8016170
#> 11    50 0.10 0.3  94       0.8012229
#> 12    50 0.10 0.5  96       0.8011854
#> 13   100 0.01 0.0  16       0.8093656
#> 14   100 0.01 0.3  16       0.8005201
#> 15   100 0.01 0.5  17       0.8078552
#> 16   100 0.10 0.0  86       0.8020410
#> 17   100 0.10 0.3  87       0.8038824
#> 18   100 0.10 0.5  88       0.8035513
```

This is an example for plotting the sample size requirements under
varying coefficients of variation of cluster sizes:

``` r
example("graph.H2x2Factorial")
#> 
#> g.H22F> #Make a plot under the test for marginal cluster-level treatment effect
#> g.H22F> graph.H2x2Factorial(power=0.9, test="cluster", rho=0.1, verbose=FALSE)
```

<img src="man/figures/README-example graph-1.png" width="100%" />
