# Internal Computational Function for SemiparMF

Implements the iterative backfitting algorithm with Cochrane-Orcutt
updates as described in Malabanan et al. (2022). This function handles
the core estimation loop, separating the non-parametric spline component
from the parametric and temporal components.

## Usage

``` r
.semipar_internal(y, x_high, z_low, w_mat, max_iter = 50, tol = 1e-04, ...)
```

## Arguments

- y:

  Numeric matrix (N x T). The response variable in wide format
  (rows=locations, cols=time).

- x_high:

  Numeric array (N x T x K). The high-frequency covariate.

  - N: Number of spatial locations.

  - T: Number of low-frequency time points (matching `y`).

  - K: Frequency ratio (e.g., 3 for monthly data predicting quarterly
    response).

- z_low:

  Numeric matrix (N x T). The parametric covariate (same frequency as
  `y`).

- w_mat:

  Numeric matrix (N x T). The spatial neighborhood covariate (e.g.,
  spatial lag of Z).

- max_iter:

  Integer. Maximum number of backfitting iterations. Default is 50.

- tol:

  Numeric. Convergence tolerance based on percentage change in Mean
  Squared Prediction Error (MSPE). Default is 1e-4.

- ...:

  Additional arguments passed to
  [`smooth.spline`](https://rdrr.io/r/stats/smooth.spline.html) (e.g.,
  `spar`).

## Value

A list containing:

- coefficients:

  List of scalar estimates: `beta` (covariate effect), `gamma` (spatial
  effect), and `rho` (temporal autocorrelation).

- nonparam:

  List containing the `spline` object and `f_hat` (estimated
  non-parametric component).

- residuals:

  Matrix (N x T) of pure structural residuals \\(Y - \hat{Y})\\.

- fitted.values:

  Matrix (N x T) of fitted values.

- history:

  Vector of MSPE values per iteration.

- iters:

  Number of iterations performed.

## References

Malabanan, V. A., Lansangan, J. R. G., & Barrios, E. B. (2022).
Semiparametric Spatiotemporal Model with Mixed Frequencies: With
Application in Crop Forecasting. *Science & Engineering Journal*, 15(2),
90-107.
