source("R/truncated_normal_calibration/trunc_norm_moments.R")
source("R/utils/sample_distribution.R")
set.seed(8132026)

check <- function(label, got, expected, tol) {
  cat(sprintf("%-52s got %11.6f   expected %11.6f   diff %+9.2e   %s\n",
              label, got, expected, got - expected,
              ifelse(abs(got - expected) < tol, "PASS", "FAIL")))
}

## Tests r2_map against values found in r2_sim.R r2 attenuation test
cat("== r2_map vs Aug 4 simulation table (tol = 4*mc_se ~ 0.006) ==\n")
Xbig <- rnorm(1e6)
sim <- data.frame(
  r2 = rep(c(0.30, 0.60), each = 7),
  b  = rep(c(Inf, 3, 2.5, 2, 1.5, 1.25, 1), 2),
  ref = c(0.3000, 0.2298, 0.1939, 0.1504, 0.1033, 0.0781, 0.0533,
          0.6000, 0.5821, 0.5626, 0.5219, 0.4500, 0.3981, 0.3297))
for (i in seq_len(nrow(sim))) {
  got <- r2_map(0, 1, sqrt((1 - sim$r2[i]) / sim$r2[i]), Xbig, -sim$b[i], sim$b[i])
  check(sprintf("latent %.2f, bounds +/-%-4s", sim$r2[i], format(sim$b[i])),
        got, sim$ref[i], 0.006)
}

## Tests calibrate_endo against sigma_eps identity and an inversion with calculated moments from m_tn and v_tn
cat("\n== calibrate_endo: analytic tests ==\n")
X <- rnorm(5000)
e1 <- calibrate_endo(0, 1, 0.30, X, -50, 50)
check("identity - non-binding bounds: sigma_eps", e1$sigma_eps, sqrt(0.7), 1e-6)   # sigma_eps = sqrt((1-r2_t)*var_t) is sample-independent
cat(sprintf("%-52s beta0 = %.5f  beta1 = %.5f  termcd = %d\n",
            "identity: full solution", e1$beta0, e1$beta1, e1$termcd))
b0 <- 0.3; b1 <- 0.9; s <- 0.7; a <- -1.5; b <- 2
eta <- b0 + b1 * X
m <- m_tn(eta, s, a, b); v <- v_tn(eta, s, a, b)
tgt <- c(mean(m), mean(v) + mean((m - mean(m))^2), r2_map(b0, b1, s, X, a, b))
e2 <- calibrate_endo(tgt[1], tgt[2], tgt[3], X, a, b)
check("inversion: beta0",     e2$beta0,     0.3, 1e-5)
check("inversion: beta1",     e2$beta1,     0.9, 1e-5)
check("inversion: sigma_eps", e2$sigma_eps, 0.7, 1e-5)

## Tests calibrate_endo using sample_normal against chosen mean, var, and r2 for datasets of 5000 and 50000
cat("\n== Round-trip via sample_normal (targets 0, 0.5, 0.30 on +/-1.5) ==\n")
for (n in c(5000, 50000)) {
  Xn  <- rnorm(n)
  cal <- calibrate_endo(0, 0.5, 0.30, Xn, -1.5, 1.5)
  y   <- sample_normal(n, X = matrix(Xn, ncol = 1), beta1 = cal$beta1,
                       beta0 = cal$beta0, sd = cal$sigma_eps,
                       min = -1.5, max = 1.5)
  fit <- summary(lm(y ~ Xn))$r.squared
  cat(sprintf("-- n = %d  (single-dataset, tolerances = 4 x single-dataset sd) --\n", n))
  check("realized mean(y)", mean(y), 0,    4 * sqrt(0.5 / n))
  check("realized var(y)",  var(y),  0.5,  4 * 0.5 * sqrt(2 / n))
  check("realized R^2",     fit,     0.30, 4 * 2*sqrt(0.30)*(1-0.30)/sqrt(n))
}