
<img width="200" src="man/figures/hubmeta-logo.png?raw=TRUE" alt="hubmeta logo" align="left">

# Hubmeta
## Meta-Analysis and Matrix Completion Toolkit

<!-- badges: start -->
<!-- badges: end -->

Hubmeta is an R toolkit for psychometric meta-analysis and matrix
completion workflows. The current development version includes:

- Hunter-Schmidt style meta-analysis via `meta_analysis()`
- Morris-weight meta-analysis via `morris_weight_analysis()`
- early MICA utilities for diagnosing and deterministically completing
  partial correlation matrices

## Installation

You can install the development version of hubmeta from
[GitHub](https://github.com/hubmeta/R) with:

``` r
install.packages("devtools")
devtools::install_github("hubmeta/R")
```

## Example

This is a basic Stage-1 meta-analysis example:

``` r
library(hubmeta)
data <- meta_analysis(c(.18, .0, .08, .15, .27, .1, .28, .17, .02, .28),
                        c(426, 328, 122, 284, 472, 154, 372, 674, 110, 116),
                        c(.85, .77, .80, .86, .80, .79, .91, .85, .92, .85),
                        c(.63, .63, .62, .39, .24, .85, .89, .48, .68, .84),
                        c(0.95, 0.80)
)
```

This is a small simulated MICA example:

``` r
library(hubmeta)

toy <- matrix(c(
  1.00, 0.32,   NA, 0.21,
  0.32, 1.00, 0.28,   NA,
    NA, 0.28, 1.00, 0.41,
  0.21,   NA, 0.41, 1.00
), nrow = 4, byrow = TRUE)

fit <- mica(toy)
fit
as.matrix(fit)
```
