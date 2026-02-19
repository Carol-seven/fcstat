# Performance Measures

Calculate various measures to evaluate the performance of the estimated
precision matrix.

## Usage

``` r
performance(hatOmega, Omega)
```

## Arguments

- hatOmega:

  The estimated precision matrix.

- Omega:

  The reference precision matrix.

## Value

A list containing the following components:

- Fnorm:

  Frobenius (Hilbert-Schmidt) norm between the true and estimated
  precision matrices.

- KL:

  Kullback-Leibler divergence between the true and estimated precision
  matrices.

- Snorm:

  Spectral (operator) norm of the difference between the true and
  estimated precision matrices.

- precision:

  Precision measure, the ratio of true positives to the total predicted
  positives.

- recall:

  Recall measure, also known as Sensitivity, the ratio of true positives
  to the total actual positives.

- specificity:

  Specificity measure, the ratio of true negatives to the total actual
  negatives.

- F1:

  F1 score, the harmonic mean of Precision and Recall.

- MCC:

  Matthews correlation coefficient, a measure of the quality of binary
  classifications.

- sparsity:

  The proportion of zeros among edges in the estimated precision matrix.
