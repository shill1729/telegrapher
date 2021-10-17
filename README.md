
# telegrapher

<!-- badges: start -->
<!-- badges: end -->

The goal of telegrapher is to provide functions for solving the telegrapher PDE and studying the behavior of its solution as well as the associated telegraph or Kac process.

## Installation

You can install the latest version of the package using devtools:

``` r
devtools::install_github("shill1729/telegrapher")
```

## Example

TODO

``` r
library(telegrapher)
# Formal analog of Fokker Planck in one dimensional Minkowski space
# Time intevral and initial point
tt <- 0.5
x0 <- 5
# Default numeric parameters are used.
w <- telegraph_model(x0, tt)
print(w)
```

