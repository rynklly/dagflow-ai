source("R/truncated_normal_calibration/trunc_norm_moments.R")
source("R/utils/sample_distribution.R")


set.seed(8242026)

## Oracle check: r2_map_multi at k = 1 must reproduce r2_map to floating point.
## The quadratic form reduces to the scalar map exactly at k = 1 (1x1 algebra),
## so agreement is structural, not statistical. Tolerance 1e-12 (arithmetic
## reordering only; any real defect misses by orders of magnitude).

X  <- rnorm(5000)
Xm <- matrix(X, ncol = 1)

tol <- 1e-12
check <- function(label, beta0, beta1, sigma_eps, a, b, c_val, b_init) {
  r_old <- r2_map(beta0, beta1, sigma_eps, X, a, b)
  r_new <- r2_map_multi(beta0, c_val, sigma_eps, Xm, b_init, a, b)
  ok <- abs(r_old - r_new) < tol
  cat(sprintf("%s  old=%.15f  new=%.15f  %s\n",
              label, r_old, r_new, if (ok) "PASS" else "FAIL"))
}

## beta1 = c_val * b_init in every case
check("symmetric two-sided ", 0.0,  1.0, 0.8, -1.5, 1.5, 1.0,  1.0)
check("asymmetric two-sided", 0.5,  0.7, 1.2, -1.0, 2.5, 0.7,  1.0)
check("one-sided lower     ", 0.0,  1.0, 0.8,  0.0, Inf, 1.0,  1.0)
check("one-sided upper     ", 1.0, -0.6, 0.5, -Inf, 2.0, -0.6, 1.0)
check("untruncated         ", 0.0,  1.0, 0.8, -Inf, Inf, 1.0,  1.0)
check("beta_init carries   ", 0.0, -1.4, 0.9, -2.0, 1.0, 0.7, -2.0)

## Simulation check: r2_map_multi prediction vs realized OLS R^2, k >= 2.
## Known-parameter design (no calibrator): fix (beta0, c, sigma_eps, beta_init),
## predict R^2 via the map, then sample through sample_normal and regress.
## Standard: |map - mean realized| within ~4 * mc_se.

n    <- 5000
reps <- 200

run_condition <- function(label, X, beta0, c_val, b_init, sigma_eps, a, b) {
  pred <- r2_map_multi(beta0, c_val, sigma_eps, X, b_init, a, b)
  beta1 <- c_val * b_init
  r2 <- numeric(reps)
  for (r in seq_len(reps)) {
    y <- sample_normal(n, X = X, beta1 = beta1, beta0 = beta0,
                       sd = sigma_eps, min = a, max = b)
    r2[r] <- summary(stats::lm(y ~ X))$r.squared
  }
  se <- stats::sd(r2) / sqrt(reps)
  dev <- abs(mean(r2) - pred) / se
  cat(sprintf("%-28s map=%.4f  realized=%.4f  mc_se=%.5f  |dev|/se=%.2f  %s\n",
              label, pred, mean(r2), se, dev,
              if (dev < 4) "PASS" else "FAIL"))
}

## Condition 1: k=2 independent parents, symmetric bounds
X1 <- cbind(rnorm(n), rnorm(n))
run_condition("k2 independent sym",  X1, 0, 0.6, c(1, 1),    0.8, -1.5, 1.5)

## Condition 2: k=2 correlated (rho = 0.6) -- overlap correction load-bearing
z1 <- rnorm(n); z2 <- 0.6 * z1 + sqrt(1 - 0.36) * rnorm(n)
X2 <- cbind(z1, z2)
run_condition("k2 rho=.6 sym",       X2, 0, 0.6, c(1, 1),    0.8, -1.5, 1.5)

## Condition 3: same X, mixed-sign beta_init -- signed cancellation
run_condition("k2 rho=.6 mixed sign", X2, 0, 0.6, c(1, -1.5), 0.8, -1.5, 1.5)

## Condition 4: k=3, one uniform parent, asymmetric bounds
X3 <- cbind(rnorm(n), 0.6 * z1 + sqrt(1 - 0.36) * rnorm(n), runif(n, -2, 2))
run_condition("k3 mixed asym",       X3, 0.3, 0.5, c(1, 0.7, -0.8), 1.0, -1, 2.5)

## Condition 5: k=2 unbounded -- untruncated positive control
run_condition("k2 untruncated",      X1, 0, 0.6, c(1, 1),    0.8, -Inf, Inf)