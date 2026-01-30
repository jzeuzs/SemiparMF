# Fit Semiparametric Spatiotemporal Model with Mixed Frequencies

Fits a spatiotemporal model where the response variable is observed at a
lower frequency (e.g., quarterly) than a non-parametric covariate (e.g.,
monthly). The model combines a non-parametric component for the
high-frequency predictor, parametric components for low-frequency
predictors and spatial neighborhood effects, and an autoregressive error
structure.

## Usage

``` r
semiparMF(
  formula,
  data_sf,
  high_freq_data,
  time_col,
  id_col,
  w_matrix = NULL,
  ...
)
```

## Arguments

- formula:

  A `formula` object (e.g., `Y ~ Z`). The left-hand side is the response
  variable (low frequency). The first predictor on the right-hand side
  is the parametric covariate (\$Z\$) measured at the same frequency.

- data_sf:

  An `sf` object containing the panel data in long format. Must contain
  columns for the response, the parametric covariate, the time index,
  and the location ID.

- high_freq_data:

  A numeric array of dimensions (N x T x K), where:

  - `N`: Number of unique spatial locations (must match `data_sf`).

  - `T`: Number of time points (must match `data_sf`).

  - `K`: The frequency ratio (e.g., 3 if predictor is monthly and
    response is quarterly).

- time_col:

  Character string. The name of the column in `data_sf` representing the
  time index.

- id_col:

  Character string. The name of the column in `data_sf` representing the
  location ID.

- w_matrix:

  Optional numeric matrix (N x T). A pre-calculated spatial weight or
  neighborhood variable. If `NULL` (default), a spatial lag of the
  variable \$Z\$ is calculated using Queen Contiguity weights.

- ...:

  Additional arguments passed to the internal backfitting function
  (e.g., `max_iter`, `tol`, `spar`).

## Value

An object of class `semiparMF` containing:

- coefficients:

  A list of estimated parameters: `beta` (parametric covariate effect),
  `gamma` (neighborhood effect), and `rho` (autoregressive parameter).

- nonparam:

  A list containing the fitted smoothing spline and the aggregated
  non-parametric component `f_hat`.

- residuals:

  Matrix (N x T) of model residuals.

- fitted.values:

  Matrix (N x T) of the fitted values (structural part only).

- dims:

  Dimensions of the data (N, T).

- meta:

  Metadata containing location IDs and time indices.

- history:

  Convergence history (MSPE per iteration).

- call:

  The function call.

## References

Malabanan, V. A., Lansangan, J. R. G., & Barrios, E. B. (2022).
Semiparametric Spatiotemporal Model with Mixed Frequencies: With
Application in Crop Forecasting. *Science & Engineering Journal*, 15(2),
90-107.

## Examples

``` r
# Simulate data using the package's included simulation function
sim <- simulate_semipar_data(n_side = 4, t_len = 20, k = 3, rho_error = 0.5)

# Fit the model
fit <- semiparMF(
  formula = Y ~ Z, 
  data_sf = sim$data, 
  high_freq_data = sim$X_high,
  time_col = "time_id", 
  id_col = "location_id"
)
#> Calculating spatial weights (Queen Contiguity)...

# Inspect results
summary(fit)
#> 
#> -- Model Summary --
#> Call:
#> semiparMF(formula = Y ~ Z, data_sf = sim$data, high_freq_data = sim$X_high, 
#>     time_col = "time_id", id_col = "location_id")
#> 
#> Coefficients:
#>   Beta (Parametric):    0.4805 
#>   Gamma (Spatial):      -0.02814 
#>   Rho (Temporal):       0.01043 
#> 
#> Residuals:
#>       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#> -6.9374035 -1.6515231 -0.1171669  0.0000537  1.5207869 10.1848119 
#> 
#> Convergence:
#>   Iterations:  5 
#>   Final MSPE:  6.249991 
```
