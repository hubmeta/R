
# hubmeta

<!-- badges: start -->
<!-- badges: end -->

Hubmeta is a meta-analysis toolkit.

## Installation

You can install the released version of hubmeta from [GITHUB](https://github.com/hubmeta/R) with:

``` r
remotes::install_github("psychmeta/psychmeta")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(hubmeta)
data <- meta_analysis(c(.18, .0, .08, .15, .27, .1, .28, .17, .02, .28),
                        c(426, 328, 122, 284, 472, 154, 372, 674, 110, 116),
                        c(.85, .77, .80, .86, .80, .79, .91, .85, .92, .85),
                        c(.63, .63, .62, .39, .24, .85, .89, .48, .68, .84),
                        c(0.95, 0.80)
)
```

