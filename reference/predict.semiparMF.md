# Predict Method for SemiparMF

Generates predictions for new data. Reference: "Predicted values are
calculated by simply adding up the scores... and the linear
combination".

## Usage

``` r
# S3 method for class 'semiparMF'
predict(object, new_high_freq, new_z, new_w, ...)
```

## Arguments

- object:

  An object of class `semiparMF`.

- new_high_freq:

  Numeric array (N x T_new x K) for the high-frequency covariate.

- new_z:

  Numeric matrix (N x T_new) for the parametric covariate.

- new_w:

  Numeric matrix (N x T_new) for the neighborhood covariate.

- ...:

  Additional arguments.

## Value

A matrix (N x T_new) of predicted values.
