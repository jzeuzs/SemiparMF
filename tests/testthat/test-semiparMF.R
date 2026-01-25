library(sf)
library(spdep)

make_small_sim <- function(n_side = 4, t_len = 10, k = 3) {
    simulate_semipar_data(n_side = n_side, t_len = t_len, k = k)
}

test_that("Model returns correct class and components", {
    sim <- make_small_sim()
    model <- semiparMF(Y ~ Z, sim$data, sim$X_high, "time_id", "location_id")

    # Check Class
    expect_s3_class(model, "semiparMF")

    # Check Essential Components
    expected_names <- c("coefficients", "nonparam", "residuals", "fitted.values", "history", "dims")
    expect_true(all(expected_names %in% names(model)))

    # Check Coefficient Structure (Scalar values, no names attribute issues)
    expect_type(model$coefficients$beta, "double")
    expect_type(model$coefficients$gamma, "double")
    expect_type(model$coefficients$rho, "double")

    # Sanity Check: Rho should be a correlation (-1 to 1)
    expect_true(model$coefficients$rho >= -1 && model$coefficients$rho <= 1)

    # Check Dimensions of fitted values
    n <- nrow(sim$data) / 10 # 16 locations
    t <- 10
    expect_equal(dim(model$fitted.values), c(n, t))
})

test_that("Predict method works with new data", {
    sim <- make_small_sim(n_side = 4, t_len = 20)
    train_ids <- 1:16
    train_times <- 1:15

    # Train on first 15 time points
    train_data <- sim$data[sim$data$time_id <= 15, ]
    train_x <- sim$X_high[, 1:15, ]

    model <- semiparMF(Y ~ Z, train_data, train_x, "time_id", "location_id")

    # Predict for the 16th time point (for all locations)
    # Extract "New" Data manually
    test_mask <- sim$data$time_id == 16
    new_z <- matrix(sim$data$Z[test_mask], ncol = 1)
    new_w <- matrix(sim$data$W[test_mask], ncol = 1)
    new_x <- sim$X_high[, 16:16, , drop = FALSE] # Keep array structure 16x1x3

    preds <- predict(model, new_high_freq = new_x, new_z = new_z, new_w = new_w)

    # Validation
    expect_equal(dim(preds), c(16, 1)) # N x T_new
    expect_false(any(is.na(preds)))
    expect_type(preds, "double")
})

test_that("Model catches input errors", {
    sim <- make_small_sim()

    # Case A: High-Freq Data Dimension Mismatch
    # X_high is N x T x K. Let's mess up T.
    bad_x <- array(runif(16 * 100 * 3), dim = c(16, 100, 3))

    expect_error(
        semiparMF(Y ~ Z, sim$data, bad_x, "time_id", "location_id"),
        "High frequency array dimension mismatch"
    )

    # Case B: Unbalanced Panel
    # Remove one observation
    bad_data <- sim$data[-1, ]

    expect_error(
        semiparMF(Y ~ Z, bad_data, sim$X_high, "time_id", "location_id"),
        "Panel is unbalanced"
    )

    # Case C: Not an SF object
    df_data <- sf::st_drop_geometry(sim$data)
    expect_error(
        semiparMF(Y ~ Z, df_data, sim$X_high, "time_id", "location_id"),
        "Data must be an sf object"
    )
})

test_that("Standard S3 methods run without error", {
    sim <- make_small_sim()
    model <- semiparMF(Y ~ Z, sim$data, sim$X_high, "time_id", "location_id")

    # Capture output to avoid cluttering test console
    expect_output(print(model), "Semiparametric Spatiotemporal Model")

    summ <- summary(model)
    expect_s3_class(summ, "summary.semiparMF")
    expect_output(print(summ), "Model Summary")

    # Plot should not error (can't check visual output easily)
    pdf(NULL) # Prevent plot window from popping up
    expect_silent(plot(model))
    dev.off()
})

test_that("Model fits better than a flat line (MSPE check)", {
    sim <- make_small_sim(n_side = 5, t_len = 20)
    model <- semiparMF(Y ~ Z, sim$data, sim$X_high, "time_id", "location_id")

    # Calculate MSPE of the model
    model_mse <- tail(model$history, 1)

    # Calculate MSE of a naive intercept-only model (mean)
    y_vec <- sim$data$Y
    naive_mse <- mean((y_vec - mean(y_vec))^2)

    # The sophisticated model should beat the mean
    expect_lt(model_mse, naive_mse)
})
