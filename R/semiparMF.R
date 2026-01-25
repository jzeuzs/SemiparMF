#' Internal Computational Function for SemiparMF
#'
#' Implements the iterative backfitting algorithm with Cochrane-Orcutt updates as described
#' in Malabanan et al. (2022). This function handles the core estimation loop, separating
#' the non-parametric spline component from the parametric and temporal components.
#'
#' @param y Numeric matrix (N x T). The response variable in wide format (rows=locations, cols=time).
#' @param x_high Numeric array (N x T x K). The high-frequency covariate.
#'   - N: Number of spatial locations.
#'   - T: Number of low-frequency time points (matching `y`).
#'   - K: Frequency ratio (e.g., 3 for monthly data predicting quarterly response).
#' @param z_low Numeric matrix (N x T). The parametric covariate (same frequency as `y`).
#' @param w_mat Numeric matrix (N x T). The spatial neighborhood covariate (e.g., spatial lag of Z).
#' @param max_iter Integer. Maximum number of backfitting iterations. Default is 50.
#' @param tol Numeric. Convergence tolerance based on percentage change in Mean Squared Prediction Error (MSPE). Default is 1e-4.
#' @param ... Additional arguments passed to \code{\link[stats]{smooth.spline}} (e.g., \code{spar}).
#'
#' @return A list containing:
#'   \item{coefficients}{List of scalar estimates: \code{beta} (covariate effect), \code{gamma} (spatial effect), and \code{rho} (temporal autocorrelation).}
#'   \item{nonparam}{List containing the \code{spline} object and \code{f_hat} (estimated non-parametric component).}
#'   \item{residuals}{Matrix (N x T) of pure structural residuals \eqn{(Y - \hat{Y})}.}
#'   \item{fitted.values}{Matrix (N x T) of fitted values.}
#'   \item{history}{Vector of MSPE values per iteration.}
#'   \item{iters}{Number of iterations performed.}
#'
#' @references
#' Malabanan, V. A., Lansangan, J. R. G., & Barrios, E. B. (2022). Semiparametric Spatiotemporal Model with Mixed Frequencies: With
#' Application in Crop Forecasting. \emph{Science & Engineering Journal}, 15(2), 90-107.
#'
#' @keywords internal
.semipar_internal <- function(y, x_high, z_low, w_mat, max_iter = 50, tol = 1e-4, ...) {
    # Dimensions
    dims <- dim(x_high)
    n <- dims[1]
    t_len <- dims[2]
    k <- dims[3]

    # Vectorize inputs for spline
    x_vec <- as.vector(x_high)

    # Initialize components to zero as per standard backfitting start
    f_hat <- matrix(0, nrow = n, ncol = t_len)
    beta_hat <- 0
    gamma_hat <- 0
    rho_hat <- 0

    # Initialize residuals for the first Cochrane-Orcutt update
    pure_resid <- matrix(0, nrow = n, ncol = t_len)

    mspe_history <- c()
    old_mspe <- Inf

    for (iter in 1:max_iter) {
        # Update the dependent variable to account for temporal autocorrelation.
        # Y_star = Y_t - rho * e_{t-1}
        resid_lag <- cbind(rep(0, n), pure_resid[, 1:(t_len - 1)])
        y_star <- y - (rho_hat * resid_lag)

        # Estimate f(X) using the whitened response Y_star.
        # Target: Y_star - (Parametric Effects)
        r_np <- y_star - (beta_hat * z_low + gamma_hat * w_mat)

        # Scale residuals to match high-frequency vector length
        r_vec_scaled <- rep(as.vector(r_np) / k, times = k)

        # Fit Smoothing Spline
        spline_args <- list(x = x_vec, y = r_vec_scaled, ...)
        if (!"spar" %in% names(spline_args)) {
            spline_args$spar <- 0.7
        }

        spline_fit <- do.call(stats::smooth.spline, spline_args)

        # Predict and Aggregate f(X)
        f_val_vec <- stats::predict(spline_fit, x_vec)$y
        f_val_arr <- array(f_val_vec, dim = c(n, t_len, k))
        f_hat <- apply(f_val_arr, c(1, 2), sum)

        # Estimate Beta and Gamma using OLS on the partial residuals.
        # Because we use Y_star, this approximates Weighted Least Squares (GLS).
        e_step2 <- y_star - f_hat

        betas <- numeric(n)
        gammas <- numeric(n)

        for (i in 1:n) {
            dat <- data.frame(r = e_step2[i, ], z = z_low[i, ], w = w_mat[i, ])
            mod <- stats::lm(r ~ z + w, data = dat)

            coefs <- stats::coef(mod)
            betas[i] <- coefs["z"]
            gammas[i] <- coefs["w"]
        }
        beta_hat <- mean(betas, na.rm = TRUE)
        gamma_hat <- mean(gammas, na.rm = TRUE)

        # Estimate Rho from the structural residuals of the original model.
        structural_fit <- f_hat + (beta_hat * z_low) + (gamma_hat * w_mat)
        pure_resid <- y - structural_fit

        rhos <- numeric(n)
        for (i in 1:n) {
            resid_i <- pure_resid[i, ]

            # Fit AR(1) to centered residuals
            fit_ar <- try(
                stats::arima(resid_i - mean(resid_i), order = c(1, 0, 0), include.mean = FALSE),
                silent = TRUE
            )

            if (!inherits(fit_ar, "try-error")) {
                rhos[i] <- stats::coef(fit_ar)["ar1"]
            } else {
                rhos[i] <- 0
            }
        }

        rho_hat <- mean(rhos, na.rm = TRUE)

        # Calculate MSPE of the full model (including AR component)
        resid_lag_new <- cbind(rep(0, n), pure_resid[, 1:(t_len - 1)])
        y_pred_full <- structural_fit + (rho_hat * resid_lag_new)

        curr_mspe <- mean((y - y_pred_full)^2)
        mspe_history <- c(mspe_history, curr_mspe)

        pct_change <- abs(curr_mspe - old_mspe) / old_mspe
        if (iter > 1 && pct_change < tol) {
            break
        }

        old_mspe <- curr_mspe
    }

    list(
        coefficients = list(
            beta = as.numeric(beta_hat),
            gamma = as.numeric(gamma_hat),
            rho = as.numeric(rho_hat)
        ),
        nonparam = list(spline = spline_fit, f_hat = f_hat),
        residuals = pure_resid,
        fitted.values = structural_fit,
        history = mspe_history,
        iters = iter
    )
}

#' Fit Semiparametric Spatiotemporal Model with Mixed Frequencies
#'
#' Fits a spatiotemporal model where the response variable is observed at a lower frequency
#' (e.g., quarterly) than a non-parametric covariate (e.g., monthly). The model combines
#' a non-parametric component for the high-frequency predictor, parametric components for
#' low-frequency predictors and spatial neighborhood effects, and an autoregressive error structure.
#'
#' @param formula A \code{formula} object (e.g., \code{Y ~ Z}). The left-hand side is the response
#'   variable (low frequency). The first predictor on the right-hand side is the parametric
#'   covariate ($Z$) measured at the same frequency.
#' @param data_sf An \code{sf} object containing the panel data in long format. Must contain
#'   columns for the response, the parametric covariate, the time index, and the location ID.
#' @param high_freq_data A numeric array of dimensions (N x T x K), where:
#'   \itemize{
#'     \item \code{N}: Number of unique spatial locations (must match \code{data_sf}).
#'     \item \code{T}: Number of time points (must match \code{data_sf}).
#'     \item \code{K}: The frequency ratio (e.g., 3 if predictor is monthly and response is quarterly).
#'   }
#' @param time_col Character string. The name of the column in \code{data_sf} representing the time index.
#' @param id_col Character string. The name of the column in \code{data_sf} representing the location ID.
#' @param w_matrix Optional numeric matrix (N x T). A pre-calculated spatial weight or neighborhood variable.
#'   If \code{NULL} (default), a spatial lag of the variable $Z$ is calculated using Queen Contiguity weights.
#' @param ... Additional arguments passed to the internal backfitting function (e.g., \code{max_iter}, \code{tol}, \code{spar}).
#'
#' @return An object of class \code{semiparMF} containing:
#'   \item{coefficients}{A list of estimated parameters: \code{beta} (parametric covariate effect),
#'     \code{gamma} (neighborhood effect), and \code{rho} (autoregressive parameter).}
#'   \item{nonparam}{A list containing the fitted smoothing spline and the aggregated non-parametric component \code{f_hat}.}
#'   \item{residuals}{Matrix (N x T) of model residuals.}
#'   \item{fitted.values}{Matrix (N x T) of the fitted values (structural part only).}
#'   \item{dims}{Dimensions of the data (N, T).}
#'   \item{meta}{Metadata containing location IDs and time indices.}
#'   \item{history}{Convergence history (MSPE per iteration).}
#'   \item{call}{The function call.}
#'
#' @references
#' Malabanan, V. A., Lansangan, J. R. G., & Barrios, E. B. (2022). Semiparametric Spatiotemporal
#' Model with Mixed Frequencies: With Application in Crop Forecasting. \emph{Science & Engineering Journal},
#' 15(2), 90-107.
#'
#' @examples
#' \dontrun{
#' # Assuming 'corn_data' is an sf object and 'ndvi_array' is the high-freq data
#' fit <- semiparMF(Yield ~ Area, data_sf = corn_data, high_freq_data = ndvi_array,
#'                  time_col = "Quarter", id_col = "Province")
#' summary(fit)
#' plot(fit)
#' }
#'
#' @export
# fmt: skip
semiparMF <- function(formula, data_sf, high_freq_data, time_col, id_col, w_matrix = NULL, ...) { # nolint
    if (!inherits(data_sf, "sf")) {
        stop("Data must be an sf object.")
    }

    # Sort data (Location Primary, Time Secondary) to align with array
    data_sf <- data_sf[order(data_sf[[id_col]], data_sf[[time_col]]), ]

    mf <- stats::model.frame(formula, data = sf::st_drop_geometry(data_sf))
    y_vec <- stats::model.response(mf)
    z_vec <- stats::model.matrix(formula, mf)[, -1] # Exclude intercept

    ids <- unique(data_sf[[id_col]])
    times <- unique(data_sf[[time_col]])
    n <- length(ids)
    t_len <- length(times)

    if (nrow(data_sf) != n * t_len) {
        stop("Panel is unbalanced. Please fill missing values before estimation.")
    }

    if (!all(dim(high_freq_data)[1:2] == c(n, t_len))) {
        stop("High frequency array dimension mismatch. Expected (N x T x K).")
    }

    y_mat <- matrix(y_vec, nrow = n, ncol = t_len, byrow = TRUE)
    z_mat <- matrix(z_vec, nrow = n, ncol = t_len, byrow = TRUE)

    if (is.null(w_matrix)) {
        message("Calculating spatial weights (Queen Contiguity)...")

        geo_unique <- data_sf[!duplicated(data_sf[[id_col]]), ]
        nb <- spdep::poly2nb(geo_unique)
        lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

        w_mat <- matrix(0, nrow = n, ncol = t_len)
        for (t in 1:t_len) {
            # Calculate spatial lag for each time slice
            w_mat[, t] <- spdep::lag.listw(lw, z_mat[, t])
        }
    } else {
        w_mat <- w_matrix
    }

    fit <- .semipar_internal(y_mat, high_freq_data, z_mat, w_mat, ...)

    fit$call <- match.call()
    fit$dims <- list(N = n, T = t_len)
    fit$meta <- list(ids = ids, times = times)
    class(fit) <- "semiparMF"

    fit
}
