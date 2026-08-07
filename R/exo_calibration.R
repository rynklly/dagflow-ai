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

print(m_tn(0, 1, -1.5, 1.5))
print(v_tn(0, 1, -1.5, 1.5))
print(m_tn(0, 1, 1, Inf))
print(v_tn(3, 2, -Inf, Inf))
print(m_tn(1, 0.5, 0, 3))
print(v_tn(1, .5, 0, 3))

library(nleqslv)

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

print(calibrate_exo(3, 4  -50, 50))
print(calibrate_exo(0, 0.5515244, -1.5, 1.5))
print(calibrate_exo(1.027556, 0.2214789, 0, 3))

cal <- calibrate_exo(1, 0.5, 0, 3)
n <- 1e6
u <- runif(n, pnorm(0, cal$mu, cal$sigma), pnorm(3, cal$mu, cal$sigma))
y <- qnorm(u, cal$mu, cal$sigma)

print(c(mean(y), var(y)))
print(c(cal$mu, cal$sigma))

source("R/utils/sample_distribution.R")
set.seed(8072026)

y2 <- sample_normal(1e6, mean = cal$mu, sd = cal$sigma, min = 0, max = 3)
print(c(mean(y2), var(y2)))