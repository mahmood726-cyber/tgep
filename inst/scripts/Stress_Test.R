library(metafor)
source("C:/Models/TGEP_Development/R/TGEP.R")

cat("=== TGEP STRESS TEST: EXTREME HETEROGENEITY ===\n")

# 1. Simulate High Heterogeneity Data (tau2 = 0.8)
set.seed(2026)
k <- 20
true_effect <- 0.5
tau2 <- 0.8
vi <- rexp(k, 10) + 0.05
theta <- rnorm(k, true_effect, sqrt(tau2))
yi <- rnorm(k, theta, sqrt(vi))

# 2. Run TGEP
res <- tgep_meta(yi, vi, n_boot = 0) # Fast run for diagnostic
cat(sprintf("True Effect: %.2f | TGEP Est: %.4f | SE: %.4f\n", true_effect, res$estimate, res$se))

# 3. Generate Diagnostic Plots (Save to file since no UI)
png("C:/Models/TGEP_Development/Stress_Test_Diagnostics.png", width = 800, height = 400)
plot.tgep(res)
dev.off()

png("C:/Models/TGEP_Development/Sparsity_Plot.png", width = 600, height = 400)
plot_sparsity(yi, vi)
dev.off()

cat("Stress test complete. Plots saved to C:/Models/TGEP_Development/\n")

# 4. Comparative Table
reml <- rma(yi, vi, method="REML")
cat("\nCOMPARISON:\n")
cat(sprintf("REML: %.4f (95%% CI: %.4f, %.4f)\n", as.numeric(reml$beta), reml$ci.lb, reml$ci.ub))
cat(sprintf("TGEP: %.4f (95%% CI: %.4f, %.4f)\n", res$estimate, res$ci_lb, res$ci_ub))
