## Calibrator oracle: calibrate_endo_multi at k = 1, beta_init = 1 must
## reproduce calibrate_endo's solutions. Comparison at solver tolerance
## (1e-6), not machine precision: two nleqslv runs may walk different
## paths to the same root. Both termcds must be 1.
source("R/truncated_normal_calibration/trunc_norm_moments.R")
set.seed(8242026)
X  <- rnorm(5000)
Xm <- matrix(X, ncol = 1)

tol <- 1e-6
check_cal <- function(label, mu_t, var_t, r2_t, a, b) {
  old <- calibrate_endo(mu_t, var_t, r2_t, X, a, b)
  new <- calibrate_endo_multi(mu_t, var_t, r2_t, Xm, 1, a, b)
  ok <- old$termcd == 1 && new$termcd == 1 &&
        abs(old$beta0 - new$beta0) < tol &&
        abs(old$beta1 - new$c) < tol &&
        abs(old$sigma_eps - new$sigma_eps) < tol
  cat(sprintf("%s  old=(%.8f, %.8f, %.8f) tc=%d  new=(%.8f, %.8f, %.8f) tc=%d  %s\n",
              label, old$beta0, old$beta1, old$sigma_eps, old$termcd,
              new$beta0, new$c, new$sigma_eps, new$termcd,
              if (ok) "PASS" else "FAIL"))
}

check_cal("symmetric   ", 0.0, 0.5, 0.30, -1.5, 1.5)
check_cal("asymmetric  ", 0.5, 0.4, 0.25, -1.0, 2.5)
check_cal("one-sided   ", 1.0, 0.6, 0.40,  0.0, Inf)
check_cal("untruncated ", 0.0, 1.0, 0.30, -Inf, Inf)

## Round-trip suite, k = 2..5: calibrate -> engine sample -> realized vs declared.
## Single draw per row at n = 50000; reading standard: R2 within ~2 sd
## (~0.007), mean/var proportionally tighter. termcd must be 1 everywhere.

set.seed(8242026)
n <- 50000

round_trip <- function(label, X, beta_init, mu_t, var_t, r2_t, a, b) {
  cal <- calibrate_endo_multi(mu_t, var_t, r2_t, X, beta_init, a, b)
  y   <- sample_normal(nrow(X), X = X, beta1 = cal$c * beta_init,
                       beta0 = cal$beta0, sd = cal$sigma_eps, min = a, max = b)
  fit <- summary(stats::lm(y ~ X))
  cat(sprintf("%-24s tc=%d  mean %.4f (t %.2f)  var %.4f (t %.2f)  R2 %.4f (t %.2f)\n",
              label, cal$termcd, mean(y), mu_t, stats::var(y), var_t,
              fit$r.squared, r2_t))
}

## Shared building blocks
z1 <- rnorm(n); z2 <- 0.6 * z1 + sqrt(1 - 0.36) * rnorm(n)

## 1-3: the original three (unchanged -- continuity with the witnessed run)
round_trip("k2 rho=.6 sym",     cbind(z1, z2), c(1, 1),    0,   0.5, 0.30, -1.5, 1.5)
round_trip("k2 mixed sign",     cbind(z1, z2), c(1, -1.5), 0,   0.5, 0.30, -1.5, 1.5)
round_trip("k3 asym",           cbind(rnorm(n), z1, runif(n, -2, 2)),
                                c(1, 0.7, -0.8), 0.5, 0.4, 0.25, -1, 2.5)

## 4: near-collinear parents (rho = .9) -- stresses the varX inversion
z3 <- 0.9 * z1 + sqrt(1 - 0.81) * rnorm(n)
round_trip("k2 rho=.9",         cbind(z1, z3), c(1, 1),    0,   0.5, 0.30, -1.5, 1.5)

## 5: k=4 with a binary parent, one-sided bounds
X4 <- cbind(rnorm(n), z2, runif(n, 0, 1), rbinom(n, 1, 0.4))
round_trip("k4 binary one-sided", X4, c(1, -0.5, 0.8, 1.2), 1.0, 0.6, 0.35, 0, Inf)

## 6: higher R2, wider bounds
round_trip("k2 r2=.6",          cbind(z1, z2), c(1, 0.5),  0,   0.8, 0.60, -2, 2)

## 7: k=5, factor-correlated parents, mixed signs, tight bounds near the cap
f  <- rnorm(n)
X5 <- sapply(c(0.7, 0.5, 0.3, -0.4, 0.6),
             function(l) l * f + sqrt(1 - l^2) * rnorm(n))
round_trip("k5 tight bounds",   X5, c(1, -1, 0.5, 0.7, -0.3), 0, 0.25, 0.30, -1, 1)