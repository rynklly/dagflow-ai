library(nleqslv)

## Returns realized mean of TruncNormal(mu, sigma^2, a, b)
m_tn <- function(mu, sigma, a, b){                                          
    alpha <- (a-mu)/sigma
    beta <- (b-mu)/sigma
    Z <- pnorm(beta)-pnorm(alpha)

    mu + sigma*((dnorm(alpha)-dnorm(beta))/Z)
}


## Returns realized variance of TruncNormal(mu, sigma^2, a, b)
v_tn <- function(mu, sigma, a, b){                                          
    alpha <- (a-mu)/sigma
    beta <- (b-mu)/sigma
    Z <- pnorm(beta)-pnorm(alpha)
    ta <- ifelse(is.finite(alpha), alpha*dnorm(alpha), 0)                   #if unbounded, ta and tb -> 0 since exponential decay beats the polynomial factor
    tb <- ifelse(is.finite(beta), beta*dnorm(beta), 0)
    (sigma^2)*(1+((ta-tb)/Z)-((dnorm(alpha)-dnorm(beta))/Z)^2)
}

## Returns realized OLS r2 found as a function of parameters known at calibration time for a single parent.
r2_map <- function(beta0, beta1, sigma_eps, X, a, b){                       
    eta <- beta0 + beta1 * X
    m <- m_tn(eta, sigma_eps, a, b)
    v <- v_tn(eta, sigma_eps, a, b)
    covhat <- mean((X - mean(X)) * (m - mean(m)))
    varX <- mean((X - mean(X))^2)
    covhat^2 / (varX * (mean(v) + mean((m - mean(m))^2)))
}

## Exogenous calibration. Solves mu_tn(mu, sigma)=mu_t, v_tn(mu, sigma)=v_t via 2-D root finding (nleqslv).
## Starts from the untruncated solution, returns the parameters mu, sigma to achieve targets after truncation.
## Also returns termination code (termcd) to ensure convergence (need termcd=1) and vector of function values (fvec).
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

## Endogenous calibration, single parent. Solves for beta0, beta1, sigma_eps using known targets mu_t, var_t, r2_t via 3-D root finding (nleqslv).
## Starts with values solved for in the untruncated system, returns the beta0, beta1, sigma_eps needed to achieve targets after truncation.
## Also returns termination code (termcd) to ensure convergence and vector of function values (fvec).
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


r2_map_multi <- function(beta0, c, sigma_eps, X, beta_init, a, b){
    eta <- beta0 + c * drop(X %*% beta_init)
    m <- m_tn(eta, sigma_eps, a, b)
    v <- v_tn(eta, sigma_eps, a, b)

    n <- ncol(X)                                                                    ## Store number of columns of X
    covhat <- numeric(n)                                                            ## Create covhat vector
    for (j in 1:n) {                                                                ## Loop through all X columns
        xj <- X[, j]                                                                ## Store jth column of X
        covhat[j] <- mean((xj - mean(xj)) * (m - mean(m)))                          ## Compute covariance of parent j with conditional mean
    }
    
    varX <- matrix(0, n, n)
    for (j in 1:n) for (k in 1:n) {
        xj <- X[, j]
        xk <- X[, k]
        varX[j, k] <- mean((xj - mean(xj)) * (xk - mean(xk)))
    }
    varX_inv <- solve(varX)                                                         ## Invert varX
    drop(t(covhat) %*% varX_inv %*% covhat) / (mean(v) + mean((m - mean(m))^2))     ## Assemble R^2 map
}

calibrate_endo_multi <- function(mu_t, var_t, r2_t, X, beta_init, a, b) {
    fn <- function(x) {
        beta0 <- x[1]
        c <- x[2]
        sig <- x[3]
        eta <- beta0 + c * drop(X %*% beta_init)
        m <- m_tn(eta, sig, a, b)
        v <- v_tn(eta, sig, a, b)
        c(mean(m) - mu_t, mean(v) + mean((m-mean(m))^2) - var_t,
            r2_map_multi(beta0, c, sig, X, beta_init, a, b) - r2_t)
    }

    vX <- mean((drop(X %*% beta_init) - mean(drop(X %*% beta_init)))^2)
    c_0 <- sqrt(r2_t * var_t / vX)
    s_0 <- sqrt((1-r2_t) * var_t)
    b0_0 <- mu_t - c_0 * mean(drop(X %*% beta_init))
    sol <- nleqslv(c(b0_0, c_0, s_0), fn)
    list(beta0=sol$x[1], c=sol$x[2], sigma_eps=sol$x[3], termcd = sol$termcd, fvec = sol$fvec)

}
