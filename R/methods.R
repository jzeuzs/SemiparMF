#' Print Method for SemiparMF
#'
#' @param x An object of class \code{semiparMF}.
#' @param ... Additional arguments.
#' @export
print.semiparMF <- function(x, ...) {
    cat("\nSemiparametric Spatiotemporal Model (Mixed Frequency)\n")
    cat("-----------------------------------------------------\n")
    cat("Dimensions: ", x$dims$N, "Locations x", x$dims$T, "Time Points\n")
    cat("Iterations: ", x$iters, "\n")
    cat("\nCoefficients:\n")
    cat("  Beta (Z): ", format(x$coefficients$beta, digits = 4), "\n")
    cat("  Gamma (W):", format(x$coefficients$gamma, digits = 4), "\n")
    cat("  Rho (AR): ", format(x$coefficients$rho, digits = 4), "\n")
    invisible(x)
}

#' Summary Method for SemiparMF
#'
#' @param object An object of class \code{semiparMF}.
#' @param ... Additional arguments.
#' @export
summary.semiparMF <- function(object, ...) {
    res <- list(
        coefficients = object$coefficients,
        last_mspe = utils::tail(object$history, 1),
        residuals_summary = summary(as.vector(object$residuals)),
        iterations = object$iters,
        call = object$call
    )
    class(res) <- "summary.semiparMF"
    res
}

#' Print Summary for SemiparMF
#'
#' @param x An object of class \code{summary.semiparMF}.
#' @param ... Additional arguments.
#' @export
print.summary.semiparMF <- function(x, ...) {
    cat("\n-- Model Summary --\n")
    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    cat("  Beta (Parametric):   ", round(x$coefficients$beta, 5), "\n")
    cat("  Gamma (Spatial):     ", round(x$coefficients$gamma, 5), "\n")
    cat("  Rho (Temporal):      ", round(x$coefficients$rho, 5), "\n")

    cat("\nResiduals:\n")
    print(x$residuals_summary)

    cat("\nConvergence:\n")
    cat("  Iterations: ", x$iterations, "\n")
    cat("  Final MSPE: ", x$last_mspe, "\n")
}

#' Plot Convergence History
#'
#' Plots the Mean Squared Prediction Error (MSPE) across iterations.
#'
#' @param x An object of class \code{semiparMF}.
#' @param ... Additional graphical parameters.
#' @export
plot.semiparMF <- function(x, ...) {
    if (length(x$history) < 2) {
        warning("Not enough history to plot.")
        return(invisible())
    }
    graphics::plot(
        x$history,
        type = "b",
        col = "blue",
        pch = 16,
        main = "Backfitting Convergence History",
        xlab = "Iteration",
        ylab = "MSPE",
        ...
    )
    graphics::grid()
}

#' Predict Method for SemiparMF
#'
#' Generates predictions for new data.
#' Reference: "Predicted values are calculated by simply adding up the scores... and the linear combination".
#'
#' @param object An object of class \code{semiparMF}.
#' @param new_high_freq Numeric array (N x T_new x K) for the high-frequency covariate.
#' @param new_z Numeric matrix (N x T_new) for the parametric covariate.
#' @param new_w Numeric matrix (N x T_new) for the neighborhood covariate.
#' @param ... Additional arguments.
#' @return A matrix (N x T_new) of predicted values.
#' @export
predict.semiparMF <- function(object, new_high_freq, new_z, new_w, ...) {
    # 1. Nonparametric Prediction
    # We use the spline object stored in the model to predict f(X) for new X data
    x_vec <- as.vector(new_high_freq)

    # Predict using the fitted smoothing spline
    f_vec <- stats::predict(object$nonparam$spline, x_vec)$y

    # Reshape and Sum over K (frequency ratio)
    dims <- dim(new_high_freq)
    if (is.null(dims)) {
        stop("new_high_freq must be an array (N x T x K)")
    }

    f_arr <- array(f_vec, dim = dims)
    f_sum <- apply(f_arr, c(1, 2), sum)

    # 2. Parametric Prediction
    # Y_hat = f(X) + Beta*Z + Gamma*W
    # Note: We do NOT add the autoregressive term (rho * e_{t-1}) for structural prediction
    # as implied by the paper's description of "Predicted values".
    pred <- f_sum +
        (object$coefficients$beta * new_z) +
        (object$coefficients$gamma * new_w)

    pred
}
