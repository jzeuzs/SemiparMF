# Simulate Spatiotemporal Mixed Frequency Data

Generates a synthetic dataset for testing the SemiparMF model. Structure
follows Equation (7) in Malabanan et al. (2022).

## Usage

``` r
simulate_semipar_data(
  n_side = 6,
  t_len = 50,
  k = 3,
  rho_error = 0.5,
  beta = 0.5,
  gamma = 0.3
)
```

## Arguments

- n_side:

  Integer. Grid side length. Total locations N = n_side^2.

- t_len:

  Integer. Length of the time series (T).

- k:

  Integer. Frequency ratio (e.g., 3 months per quarter).

- rho_error:

  Numeric. Autocorrelation coefficient for the error term (default 0.5).

- beta:

  True coefficient for covariate Z (default 0.5).

- gamma:

  True coefficient for neighborhood variable W (default 0.3).

## Value

A list containing:

- data:

  An `sf` object containing Y, Z, W, and geometry.

- X_high:

  An array (N x T x K) of high-frequency covariates.

- true_params:

  List of true parameters used for generation.

## Examples

``` r
# Generate a small dataset
sim <- simulate_semipar_data(n_side = 4, t_len = 10, k = 3)

# Check dimensions
dim(sim$data) # 160 x 5 (16 locations * 10 time points)
#> [1] 160   6
dim(sim$X_high) # 16 x 10 x 3
#> [1] 16 10  3
```
