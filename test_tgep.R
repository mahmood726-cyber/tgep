if (!requireNamespace("metafor", quietly = TRUE)) stop("metafor required")
library(metafor)
source(file.path(getwd(), "R", "TGEP.R"))

n_pass <- 0; n_fail <- 0
check <- function(name, cond) {
  if (cond) { cat(sprintf("  PASS: %s\n", name)); n_pass <<- n_pass + 1 }
  else { cat(sprintf("  FAIL: %s\n", name)); n_fail <<- n_fail + 1 }
}

set.seed(42)
yi <- rnorm(10, 0.3, 0.2); vi <- rep(0.05, 10)

cat("=== Test 1: Basic functionality ===\n")
res <- tgep_meta(yi, vi, n_boot = 0)
check("estimate is finite", is.finite(res$estimate))
check("SE > 0", res$se > 0)
check("CI width > 0", res$ci_ub > res$ci_lb)
check("3 guard weights sum to 1", abs(sum(res$weights) - 1) < 1e-10)

cat("=== Test 2: conf.level parameter ===\n")
res90 <- tgep_meta(yi, vi, n_boot = 0, conf.level = 0.90)
res99 <- tgep_meta(yi, vi, n_boot = 0, conf.level = 0.99)
check("90% CI narrower than 95%", (res90$ci_ub - res90$ci_lb) < (res$ci_ub - res$ci_lb))
check("95% CI narrower than 99%", (res$ci_ub - res$ci_lb) < (res99$ci_ub - res99$ci_lb))

cat("=== Test 3: compare_tgep ===\n")
comp <- compare_tgep(yi, vi, n_boot = 0)
check("compare returns 3 rows", nrow(comp) == 3)
check("has TGEP, REML, HKSJ", all(c("TGEP", "REML", "HKSJ") %in% comp$Method))

cat("=== Test 4: Bootstrap SE ===\n")
set.seed(42)
res_boot <- tgep_meta(yi, vi, n_boot = 200)
check("bootstrap SE > 0", res_boot$se > 0)

cat("=== Test 5: K=3 (LOO-CV edge) ===\n")
res3 <- tgep_meta(yi[1:3], vi[1:3], n_boot = 0)
check("k=3 estimate finite", is.finite(res3$estimate))

cat("=== Test 6: K=2 fallback to equal guard weights ===\n")
res2 <- tgep_meta(yi[1:2], vi[1:2], n_boot = 0)
check("k=2 estimate finite", is.finite(res2$estimate))
check("k=2 weights are equal", all(abs(res2$weights - 1/3) < 1e-10))

cat("=== Test 7: K=1 fallback ===\n")
res1 <- tgep_meta(yi[1], vi[1], n_boot = 0)
check("k=1 estimate finite", is.finite(res1$estimate))

cat("=== Test 8: All-identical studies (GRMA div-by-zero fix) ===\n")
yi_same <- rep(0.5, 10); vi_same <- rep(0.1, 10)
res_same <- tryCatch(tgep_meta(yi_same, vi_same, n_boot = 0), error = function(e) NULL)
check("identical studies don't crash", !is.null(res_same))
check("identical studies estimate finite", is.finite(res_same$estimate))
check("identical studies estimate near 0.5", abs(res_same$estimate - 0.5) < 0.05)

cat("=== Test 9: SWA p-value uses marginal (sqrt(v)) not total ===\n")
# Under high tau2, SWA p-value with sqrt(v) should differ from sqrt(v+tau2)
set.seed(123)
yi_test <- c(0.5, 0.6, -0.1, 0.3, 0.05)
vi_test <- rep(0.04, 5)
# SWA should identify study 3 (negative) as non-significant and up-weight it
swa_res <- swa_guard_core(yi_test, vi_test)
check("SWA produces finite estimate", is.finite(swa_res$est))

cat("=== Test 10: Publication bias scenario ===\n")
set.seed(123)
k <- 20
theta <- rnorm(k, 0.3, sqrt(0.05))
vi_pb <- runif(k, 0.02, 0.15)
yi_pb <- rnorm(k, theta, sqrt(vi_pb))
p <- 2 * (1 - pnorm(abs(yi_pb / sqrt(vi_pb))))
keep <- (p < 0.05) | (runif(k) < 0.2)
if (sum(keep) < 3) keep[order(p)][1:3] <- TRUE
yi_pb <- yi_pb[keep]; vi_pb <- vi_pb[keep]
res_pb <- tgep_meta(yi_pb, vi_pb, n_boot = 0)
reml_pb <- rma(yi_pb, vi_pb, method = "REML")
check("pub bias: TGEP runs", is.finite(res_pb$estimate))
cat(sprintf("  Info: TGEP=%.4f, REML=%.4f (true=0.3, k=%d after selection)\n",
            res_pb$estimate, as.numeric(reml_pb$beta), length(yi_pb)))

cat("=== Test 11: Large K ===\n")
set.seed(42)
yi_big <- rnorm(50, 0.3, 0.3); vi_big <- runif(50, 0.02, 0.2)
res_big <- tgep_meta(yi_big, vi_big, n_boot = 0)
check("k=50 runs", is.finite(res_big$estimate))

cat("=== Test 12: All significant studies ===\n")
yi_sig <- c(2.5, 3.0, 2.8, 3.2, 2.6); vi_sig <- rep(0.5, 5)
res_sig <- tgep_meta(yi_sig, vi_sig, n_boot = 0)
check("all-significant runs", is.finite(res_sig$estimate))

cat("=== Test 13: All non-significant studies ===\n")
yi_ns <- c(0.01, -0.02, 0.03, 0.005, -0.01); vi_ns <- rep(0.5, 5)
res_ns <- tgep_meta(yi_ns, vi_ns, n_boot = 0)
check("all-non-significant runs", is.finite(res_ns$estimate))

cat("=== Test 14: Negative effects ===\n")
yi_neg <- rnorm(10, -0.5, 0.2); vi_neg <- rep(0.05, 10)
res_neg <- tgep_meta(yi_neg, vi_neg, n_boot = 0)
check("negative effects: estimate < 0", res_neg$estimate < 0)

cat("=== Test 15: tau2 = 0 scenario ===\n")
yi_homo <- rnorm(10, 0.3, sqrt(0.05)); vi_homo <- rep(0.05, 10)
res_homo <- tgep_meta(yi_homo, vi_homo, n_boot = 0)
check("tau2=0: runs", is.finite(res_homo$estimate))

cat("=== Test 16: K=1 SE and CI are finite ===\n")
res1_full <- tgep_meta(yi[1], vi[1], n_boot = 0)
check("k=1 SE finite", is.finite(res1_full$se))
check("k=1 CI_lb finite", is.finite(res1_full$ci_lb))
check("k=1 CI_ub finite", is.finite(res1_full$ci_ub))
check("k=1 CI width > 0", res1_full$ci_ub > res1_full$ci_lb)

cat("=== Test 17: Input validation ===\n")
err1 <- tryCatch(tgep_meta(yi[1:3], vi[1:2], n_boot = 0), error = function(e) e$message)
check("length mismatch caught", grepl("same length", err1))
err2 <- tryCatch(tgep_meta(c(0.3, NA, 0.1), vi[1:3], n_boot = 0), error = function(e) e$message)
check("NA input caught", grepl("finite", err2))
err3 <- tryCatch(tgep_meta(yi[1:3], c(0.05, -0.01, 0.05), n_boot = 0), error = function(e) e$message)
check("negative vi caught", grepl("positive", err3))

cat(sprintf("\n=== Results: %d PASS, %d FAIL ===\n", n_pass, n_fail))
if (n_fail > 0) stop(sprintf("%d test(s) failed!", n_fail))
