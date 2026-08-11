source("R/exo_calibration.R")
source("R/r2_map_test.R")
source("R/utils/sample_distribution.R")

calibrate_endo <- function(mu_t, var_t, r2_t, X, a, b) {
    fn <- function(x) {
        beta0 <- x[1]
        beta1 <- x[2]
        sig <- x[3]
        eta <- beta0 + beta1 * X
        m <- m_tn(eta, sig, a, b)
        v <- v_tn(eta, sig, a, b)
        c(mean(m) - mu_t, mean(v) + mean((m-mean(m))^2) - var_t,
          r2_map(beta0, beta1, sig, X, a, b) - r2_t)
    }

    vX <- mean((X - mean(X))^2)
    b1_0 <- sqrt(r2_t * var_t / vX)
    s_0 <- sqrt((1 - r2_t) * var_t)
    b0_0 <- mu_t - b1_0 * mean(X)
    sol <- nleqslv(c(b0_0, b1_0, s_0), fn)
    list(beta0 = sol$x[1], beta1 = sol$x[2], sigma_eps = sol$x[3], termcd = sol$termcd, fvec = sol$fvec)
}

X <- rnorm(50000)

print("cali_endo_1: ")
print(calibrate_endo(0, 1, 0.30, X, -50, 50))

b0 <- 0.3
b1 <- 0.9
s <- 0.7
a <- -1.5
b <- 2
eta <- b0 + b1*X
m <- m_tn(eta, s, a, b)
v <- v_tn(eta, s, a, b)
tgt <- c(mean(m), mean(v) + mean((m-mean(m))^2), r2_map(b0, b1, s, X, a, b))
print("cali_endo_2: ")
print(calibrate_endo(tgt[1], tgt[2], tgt[3], X, a, b))


cal <- calibrate_endo(0, 0.5, 0.30, X, -1.5, 1.5)
y <- sample_normal(length(X), X = matrix(X, ncol=1), beta1 = cal$beta1, 
                    beta0 = cal$beta0, sd = cal$sigma_eps, min = -1.5, max = 1.5)
print("cali_endo_3: ")
print(c(mean(y), var(y), summary(lm(y ~ X))$r.squared))