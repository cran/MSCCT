
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MSCCT

<!-- badges: start -->
<!-- badges: end -->

Multiple Survival Crossing Curves Tests

This package contains tests for comparison or two or more survival
curves when the proportional hazards hypothesis is not verified, in
particular when the survival curves cross each other.

## Installation

You can install the development version of MSCCT from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("https://github.com/HMinP/MSCCT")
```

## Example

``` r
library(MSCCT)
```

This package contains:

- The weighted log-rank test

The log-rank test compares for each group and for each time of event the
expected and the observed number of events. The weighted log-rank adds
weights to each time of event. Some (implemented) exemples are the
Flemming-Harrington test and the Gehan-Wilcoxon test. It is also
possible to chose the weights you want.

``` r
multi_lr(data_under_PH)
#> (Multiple) Weighted log-rank test 
#> 
#> Weighting : Classic log-rank test 
#> Degrees of freedom : 2 
#> 
#>        Statistic p
#> Test 1   142.864 0
```

``` r
multi_lr(data_under_PH, test="fh", rho=1, gamma=0)
#> (Multiple) Weighted log-rank test 
#> 
#> Weighting : Flemming-Harrington test 
#> Parameters : rho = 1 , gamma =  0 
#> Degrees of freedom : 2 
#> 
#>        Statistic p
#> Test 1  112.1108 0
```

- The Restricted Mean Survival Test

The Restricted Mean Survival Time at time $\tau$ is the area under a
survival curve up to time $\tau$. The RMST test compares the areas under
the survival curves and tests the equality to zero of the difference of
RMST.

``` r
multi_rmst(data_under_PH, tau=12, nboot=100, method="bonferroni")
#> (Multiple) test of RMST 
#> Truncation time : 12  
#> Correction : bonferroni 
#> 
#> RMST estimation for each arm 
#>            rmst        sd
#> arm 0 10.232101 0.1777291
#> arm 1  9.380293 0.2147322
#> arm 2  8.179610 0.2163973
#> 
#> Pair-wise comparisons 
#>             dRMST        sd            p   p adjusted
#> 0 VS 1 -0.8518077 0.2787428 2.243925e-03 6.731774e-03
#> 0 VS 2 -2.0524910 0.2800275 2.309264e-13 6.927792e-13
#> 1 VS 2 -1.2006833 0.3048569 8.198752e-05 2.459625e-04
#>  
#> p=6.927792e-13
```

- The Two-stage test

The two-stage test is a combination of two tests. The first one is a
classic log-rank test. When the log-rank test is not significant, this
means that the survival curves are either different or they cross each
other and the log-rank is not powerful enough. In order to differentiate
these cases, a second test is performed. This second test is a weighted
log-rank test with weights that allows to differentiate the two previous
cases.

``` r
multi_ts(data_under_PH, eps=0.1, nboot=100, method="BH")
#> (Multiple) Two-Staged test 
#> Correction : BH  
#> 
#>                  p1   p2            p        adj_p
#> 0 VS 1 5.486672e-09 0.75 8.357007e-08 8.357007e-08
#> 0 VS 2 0.000000e+00 0.05 3.774758e-15 1.132427e-14
#> 1 VS 2 2.401913e-10 0.86 4.813040e-09 7.219560e-09
#>  
#> p=1.132427e-14
```
