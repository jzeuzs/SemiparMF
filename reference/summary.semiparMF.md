# Summary Method for SemiparMF

Summary Method for SemiparMF

## Usage

``` r
# S3 method for class 'semiparMF'
summary(object, ...)
```

## Arguments

- object:

  An object of class `semiparMF`.

- ...:

  Additional arguments.

## Value

A list of class `summary.semiparMF` containing:

- coefficients:

  A list of estimated parameters (beta, gamma, rho).

- last_mspe:

  The Mean Squared Prediction Error from the final iteration.

- residuals_summary:

  Summary statistics of the residuals.

- iterations:

  Total number of iterations performed.

- call:

  The function call.
