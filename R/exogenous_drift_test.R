source("R/sample_distribution.R")
set.seed(7292026)

n <- 1e6

conds <- list(
    list("unbounded", -Inf, Inf, 0, 1),
    list("sym +/- 3 SD", -3, 3, 0, 1),
    list("sym +/- 2.5 SD", -2.5, 2.5, 0, 1),
    list("sym +/- 2 SD", -2, 2, 0, 1),
    list("sym +/- 1.5 SD", -1.5, 1.5, 0, 1),
    list("sym +/- 1 SD", -1, 1, 0, 1),
    list("asym [0, Inf]", 0, Inf, 1, 1),
    list("asym [-Inf, 2]", -Inf, 2, 1, 1),
    list("asym [-1,2]", -1, 2, 0, 1),
    list("asym [-2,1]", -2, 1, 0, 1)
)

cat(sprintf("%-20s | %11s %11s | %9s %9s\n",
    "condition", "target mean", "target sd", "real mean", "real sd"))
cat(strrep("-", 70), "\n")
for (cc in conds) {
  cond <- cc[[1]]; a <- cc[[2]]; b <- cc[[3]]; tm <- cc[[4]]; ts <- cc[[5]]
  x <- sample_normal(n, mean = tm, sd = ts, min = a, max = b)
  cat(sprintf("%-20s | %11.2f %11.2f | %9.4f %9.4f\n",
      cond, tm, ts, mean(x), sd(x)))
}