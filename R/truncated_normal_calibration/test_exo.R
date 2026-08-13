source("R/truncated_normal_calibration/trunc_norm_moments.R")
source("R/utils/sample_distribution.R")
set.seed(8132026)

check <- function(label, got, expected, tol) {
  cat(sprintf("%-52s got %11.6f   expected %11.6f   diff %+9.2e   %s\n",
              label, got, expected, got - expected,
              ifelse(abs(got - expected) < tol, "PASS", "FAIL")))
}

## Tests m_tn and v_tn against calculations done using the same formulas in python and wolfram alpha
cat("== Moment functions vs independent reference values ==\n")
check("m_tn(0,1,-1.5,1.5)   [symmetric]",     m_tn(0, 1, -1.5, 1.5),   0,        1e-12)
check("v_tn(0,1,-1.5,1.5)   [symmetric]",  v_tn(0, 1, -1.5, 1.5),   0.551524, 1e-6)
check("m_tn(0,1,1,Inf)      [one-sided]", m_tn(0, 1, 1, Inf),     1.525135, 1e-6)
check("v_tn(3,2,-Inf,Inf)   [identity]",     v_tn(3, 2, -Inf, Inf),   4,        1e-12)
check("m_tn(1,0.5,0,3)      [asymmetric]",   m_tn(1, 0.5, 0, 3),      1.027556, 1e-6)
check("v_tn(1,0.5,0,3)      [asymmetric]",   v_tn(1, 0.5, 0, 3),      0.221479, 1e-6)

## Tests calibrate_exo against identity and recovering chosen parameters from computed moments. (Inverts the moment values verified above)
cat("\n== calibrate_exo: analytic tests ==\n")
s1 <- calibrate_exo(3, 4, -50, 50)
check("identity: mu    (non-binding bounds)", s1$mu,    3, 1e-8)
check("identity: sigma",                      s1$sigma, 2, 1e-8)
s2 <- calibrate_exo(0, 0.5515244, -1.5, 1.5)
check("inversion - symmetric: mu   (std normal +/-1.5)", s2$mu,    0, 1e-6)
check("inversion - symmetric: sigma",                    s2$sigma, 1, 1e-6)
s3 <- calibrate_exo(1.027556, 0.2214789, 0, 3)
check("inversion - asymmetric: mu    (fwd moments of 1,0.5)", s3$mu,    1,   1e-5)
check("inversion - asymmetric: sigma",                        s3$sigma, 0.5, 1e-5)
cat(sprintf("%-52s termcd = %d %d %d   (1 = converged)\n",
            "solver exit codes", s1$termcd, s2$termcd, s3$termcd))


## Tests calibrate_exo parameters using direct inverse-CDF sampling as well as DagFlow's sample_normal
cat("\n== Round-trip: solve, sample, measure (targets 1, 0.5 on [0,3]) ==\n")
cal <- calibrate_exo(1, 0.5, 0, 3)
cat(sprintf("%-52s mu_p = %.4f   sigma_p = %.4f\n",
            "solved parameters", cal$mu, cal$sigma))
n <- 1e6                                                                        # MC tolerance with se ~ 0.0007, tol = 4*se
u  <- runif(n, pnorm(0, cal$mu, cal$sigma), pnorm(3, cal$mu, cal$sigma))
y1 <- qnorm(u, cal$mu, cal$sigma)
check("inverse-CDF sampler: mean(y)", mean(y1), 1,   0.003)
check("inverse-CDF sampler: var(y)",  var(y1),  0.5, 0.003)
y2 <- sample_normal(n, mean = cal$mu, sd = cal$sigma, min = 0, max = 3)
check("DagFlow sample_normal: mean(y)",      mean(y2), 1,   0.003)
check("DagFlow sample_normal: var(y)",       var(y2),  0.5, 0.003)