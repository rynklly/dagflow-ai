source("R/exo_calibration.R")

print(m_tn(c(-1, 0, 1), 0.8, -1.5, 1.5))

r2_map <- function(beta0, beta1, sigma_eps, X, a, b){
    eta <- beta0 + beta1 * X
    m <- m_tn(eta, sigma_eps, a, b)
    v <- v_tn(eta, sigma_eps, a, b)
    covhat <- mean((X - mean(X)) * (m - mean(m)))
    varX <- mean((X - mean(X))^2)
    covhat^2 / (varX * (mean(v) + mean((m - mean(m))^2)))
}

set.seed(8072026)

X <- rnorm(1e6)
for (r2 in c(0.30, 0.60)){
    for (b in c(Inf, 3, 2.5, 2, 1.5, 1.25, 1)) {
        cat(sprintf("latent %.2f bounds +/-%4s r2_map %.4f\n",
        r2, format(b), r2_map(0, 1, sqrt((1-r2)/r2), X, -b, b)))
    }
}