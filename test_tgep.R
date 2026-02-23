library(metafor)
source("C:/Models/TGEP_Development/R/TGEP.R")

set.seed(42)
yi <- rnorm(10, 0.3, 0.2)
vi <- rep(0.05, 10)

cat("=== Test 1: tgep_meta basic ===\n")
res <- tgep_meta(yi, vi, n_boot = 0)
cat(sprintf("Estimate: %.4f, SE: %.4f\n", res$estimate, res$se))
cat(sprintf("CI: [%.4f, %.4f]\n", res$ci_lb, res$ci_ub))
cat(sprintf("Weights: GRMA=%.3f, WRD=%.3f, SWA=%.3f\n", res$weights[1], res$weights[2], res$weights[3]))
stopifnot(abs(res$ci_ub - res$ci_lb) > 0)
cat("PASS\n\n")

cat("=== Test 2: conf.level parameter ===\n")
res90 <- tgep_meta(yi, vi, n_boot = 0, conf.level = 0.90)
res99 <- tgep_meta(yi, vi, n_boot = 0, conf.level = 0.99)
# 90% CI should be narrower than 95% which is narrower than 99%
width90 <- res90$ci_ub - res90$ci_lb
width95 <- res$ci_ub - res$ci_lb
width99 <- res99$ci_ub - res99$ci_lb
cat(sprintf("Width 90%%: %.4f, 95%%: %.4f, 99%%: %.4f\n", width90, width95, width99))
stopifnot(width90 < width95)
stopifnot(width95 < width99)
cat("PASS\n\n")

cat("=== Test 3: compare_tgep ===\n")
comp <- compare_tgep(yi, vi, n_boot = 0)
print(comp)
stopifnot(nrow(comp) == 3)
stopifnot(all(c("TGEP", "REML", "HKSJ") %in% comp$Method))
cat("PASS\n\n")

cat("=== Test 4: bootstrap SE ===\n")
set.seed(42)
res_boot <- tgep_meta(yi, vi, n_boot = 200)
cat(sprintf("Bootstrap SE: %.4f vs Analytical SE: %.4f\n", res_boot$se, res$se))
stopifnot(res_boot$se > 0)
cat("PASS\n\n")

cat("=== Test 5: small k (k=3) ===\n")
res_small <- tgep_meta(yi[1:3], vi[1:3], n_boot = 0)
cat(sprintf("k=3 Estimate: %.4f\n", res_small$estimate))
stopifnot(!is.na(res_small$estimate))
cat("PASS\n\n")

cat("=== Test 6: k=2 edge case ===\n")
res_k2 <- tgep_meta(yi[1:2], vi[1:2], n_boot = 0)
cat(sprintf("k=2 Estimate: %.4f\n", res_k2$estimate))
stopifnot(!is.na(res_k2$estimate))
cat("PASS\n\n")

cat("=== Test 7: publication bias scenario ===\n")
set.seed(123)
k <- 20
theta <- rnorm(k, 0.3, sqrt(0.05))
vi_pb <- rexp(k, 10) + 0.01
yi_pb <- rnorm(k, theta, sqrt(vi_pb))
p <- 2 * (1 - pnorm(abs(yi_pb / sqrt(vi_pb))))
keep <- (p < 0.05) | (runif(k) < 0.2)
if (sum(keep) < 3) keep <- rank(p) <= 3
yi_pb <- yi_pb[keep]; vi_pb <- vi_pb[keep]
res_pb <- tgep_meta(yi_pb, vi_pb, n_boot = 0)
reml_pb <- rma(yi_pb, vi_pb, method = "REML")
cat(sprintf("TGEP: %.4f, REML: %.4f, True: 0.3\n", res_pb$estimate, as.numeric(reml_pb$beta)))
cat("PASS\n\n")

cat("=== All 7 tests PASSED ===\n")
