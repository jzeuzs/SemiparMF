#' Simulate Spatiotemporal Mixed Frequency Data
#'
#' Generates a synthetic dataset for testing the SemiparMF model.
#' Structure follows Equation (7) in Malabanan et al. (2022).
#'
#' @param n_side Integer. Grid side length. Total locations N = n_side^2.
#' @param t_len Integer. Length of the time series (T).
#' @param k Integer. Frequency ratio (e.g., 3 months per quarter).
#' @param rho_error Numeric. Autocorrelation coefficient for the error term (default 0.5).
#' @param beta True coefficient for covariate Z (default 0.5).
#' @param gamma True coefficient for neighborhood variable W (default 0.3).
#'
#' @return A list containing:
#'   \item{data}{An `sf` object containing Y, Z, W, and geometry.}
#'   \item{X_high}{An array (N x T x K) of high-frequency covariates.}
#'   \item{true_params}{List of true parameters used for generation.}
#' @export
simulate_semipar_data <- function(
    n_side = 6,
    t_len = 50,
    k = 3,
    rho_error = 0.5,
    beta = 0.5,
    gamma = 0.3
) {
    n <- n_side * n_side
    ids <- 1:n

    # Uses a regular grid to ensure consistent neighborhood definitions
    grid <- sf::st_make_grid(
        sf::st_as_sfc(sf::st_bbox(c(
            xmin = 0,
            ymin = 0,
            xmax = n_side,
            ymax = n_side
        ))),
        n = c(n_side, n_side)
    )

    # "X_itk ... generated from U(0,1)"
    x_high <- array(stats::runif(n * t_len * k, 0, 1), dim = c(n, t_len, k))

    # "Z_it ... generated from N(100,10)"
    z <- matrix(stats::rnorm(n * t_len, 100, 10), nrow = n, ncol = t_len)

    # "W_it ... generated from Po(lambda)"
    w <- matrix(stats::rpois(n * t_len, 50), nrow = n, ncol = t_len)

    # "epsilon_it = rho * epsilon_i,t-1 + a_it" [cite: 84]
    errors <- matrix(0, nrow = n, ncol = t_len)
    for (i in 1:n) {
        errors[i, ] <- stats::arima.sim(
            list(order = c(1, 0, 0), ar = rho_error),
            n = t_len
        )
    }

    # Model: y = sum(f(x)) + beta*Z + gamma*W + error
    # We use a linear function h(x) = 2x for simulation as implied by "Linear" scenario in Table 1
    y <- matrix(0, nrow = n, ncol = t_len)

    for (i in 1:n) {
        for (t in 1:t_len) {
            # Sum of f(x) over k sub-periods
            # f(x) = 2*x (arbitrary choice for "Linear" scenario)
            f_x_sum <- sum(2 * x_high[i, t, ])

            y[i, t] <- f_x_sum +
                (beta * z[i, t]) +
                (gamma * w[i, t]) +
                errors[i, t]
        }
    }

    df <- expand.grid(location_id = ids, time_id = 1:t_len)
    df <- df[order(df$location_id, df$time_id), ]

    df$Y <- as.vector(t(y))
    df$Z <- as.vector(t(z))
    df$W <- as.vector(t(w))

    geom_df <- data.frame(location_id = ids)
    sf::st_geometry(geom_df) <- grid

    data_sf <- merge(df, geom_df, by = "location_id")
    data_sf <- sf::st_as_sf(data_sf)

    list(
        data = data_sf,
        X_high = x_high,
        true_params = list(beta = beta, gamma = gamma, rho = rho_error)
    )
}
