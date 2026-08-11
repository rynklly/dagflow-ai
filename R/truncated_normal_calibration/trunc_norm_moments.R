library(nleqslv)

m_tn <- function(mu, sigma, a, b){
    alpha <- (a-mu)/sigma
    beta <- (b-mu)/sigma
    Z <- pnorm(beta)-pnorm(alpha)
    mu + sigma*((dnorm(alpha)-dnorm(beta))/Z)
}

v_tn <- function(mu, sigma, a, b){
    alpha <- (a-mu)/sigma
    beta <- (b-mu)/sigma
    Z <- pnorm(beta)-pnorm(alpha)
    ta <- ifelse(is.finite(alpha), alpha*dnorm(alpha), 0)
    tb <- ifelse(is.finite(beta), beta*dnorm(beta), 0)
    (sigma^2)*(1+((ta-tb)/Z)-((dnorm(alpha)-dnorm(beta))/Z)^2)
}

r2_map <- function(beta0, beta1, sigma_eps, X, a, b){
    eta <- beta0 + beta1 * X
    m <- m_tn(eta, sigma_eps, a, b)
    v <- v_tn(eta, sigma_eps, a, b)
    covhat <- mean((X - mean(X)) * (m - mean(m)))
    varX <- mean((X - mean(X))^2)
    covhat^2 / (varX * (mean(v) + mean((m - mean(m))^2)))
}

calibrate_exo <- function(mu_t, var_t, a, b){
    fn <- function(x){
        c(m_tn(x[1], x[2], a, b) - mu_t,
          v_tn(x[1], x[2], a, b) - var_t)
    }
    x0 <- c(mu_t, sqrt(var_t))
    sol <- nleqslv(x0, fn)
    list(mu = sol$x[1], sigma = sol$x[2],
         termcd = sol$termcd, fvec = sol$fvec)
}

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

