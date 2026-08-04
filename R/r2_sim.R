source("R/utils/sample_distribution.R")
set.seed(7312026)

n <- 5000; reps <-60
latent_r2 <- c(0.30, 0.60)
bounds <- c(Inf, 3, 2.5, 2, 1.5, 1.25, 1)

cat(sprintf("%-10s %10s | %9s %9s\n", "latent_R2", "bounds", "mean_R2", "(mc_se)"))
cat(strrep("-", 45), "\n")
for(r2 in latent_r2){
    sig <- sqrt((1-r2)/r2)
    for (b in bounds) {
        r2s <- numeric(reps)
        for (r in seq_len(reps)){
            X <- rnorm(n)
            y <- sample_normal(n, X=X, beta1=1, beta0=0, sd=sig, min = -b, max = b)
            r2s[r] <- summary(lm(y ~ X))$r.squared
        }
        cat(sprintf("%-10.2f %10s | %9.4f (%6.4f)\n",
            r2, ifelse(is.finite(b), sprintf("+/-%.2f", b), "none"), mean(r2s), sd(r2s)/sqrt(reps)))
    }
}