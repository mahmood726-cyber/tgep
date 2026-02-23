#!/usr/bin/env Rscript
# TGEP Comprehensive Simulation Study (ADEMP Framework)
# Author: Mahmood Ul Hassan
# Date: 2026-02-23
#
# ADEMP:
#   Aims: Compare TGEP vs REML coverage, bias, RMSE under heterogeneity & pub bias
#   Data-generating mechanisms: Normal random-effects with/without pub bias
#   Estimands: Pooled mean effect (true_mu)
#   Methods: REML, HKSJ, TGEP (T=1.0, n_boot=0)
#   Performance: Bias, RMSE, Coverage, MCSE

library(metafor)
source("C:/Models/TGEP_Development/R/TGEP.R")

set.seed(20260223)
N_SIM <- 500
t_global <- Sys.time()

# Data-generating mechanism
generate_meta_data <- function(k, true_mu, tau2, bias_type = "none", bias_strength = 0.8) {
  theta_i <- rnorm(k, true_mu, sqrt(tau2))
  vi <- 1 / rgamma(k, shape = 5, rate = 5 * 20)  # ~N(20) per study -> se^2 ~ 1/20
  yi <- rnorm(k, theta_i, sqrt(vi))

  if (bias_type == "significance") {
    p <- 2 * (1 - pnorm(abs(yi / sqrt(vi))))
    prob_pub <- ifelse(p < 0.05, 1.0, 1 - bias_strength)
    keep <- runif(k) < prob_pub
    if (sum(keep) < 3) keep[order(p)][1:3] <- TRUE
    yi <- yi[keep]; vi <- vi[keep]
  } else if (bias_type == "onesided") {
    # One-sided: only publish if effect in expected direction
    prob_pub <- ifelse(yi > 0, 1.0, 1 - bias_strength)
    keep <- runif(k) < prob_pub
    if (sum(keep) < 3) keep[order(-yi)][1:3] <- TRUE
    yi <- yi[keep]; vi <- vi[keep]
  }

  list(yi = yi, vi = vi, k_final = length(yi))
}

# Define scenarios
scenarios <- list(
  # Type I / Null scenarios
  list(name = "Null_k10_tau0.05", type = "null", k = 10, true_mu = 0.0, tau2 = 0.05, bias_type = "none"),
  list(name = "Null_k10_tau0.20", type = "null", k = 10, true_mu = 0.0, tau2 = 0.20, bias_type = "none"),

  # Baseline (no bias) — varying k
  list(name = "Base_k5",  type = "base", k = 5,  true_mu = 0.3, tau2 = 0.05, bias_type = "none"),
  list(name = "Base_k10", type = "base", k = 10, true_mu = 0.3, tau2 = 0.05, bias_type = "none"),
  list(name = "Base_k20", type = "base", k = 20, true_mu = 0.3, tau2 = 0.05, bias_type = "none"),
  list(name = "Base_k30", type = "base", k = 30, true_mu = 0.3, tau2 = 0.05, bias_type = "none"),

  # Varying heterogeneity
  list(name = "Hetero_tau0.01", type = "hetero", k = 15, true_mu = 0.3, tau2 = 0.01, bias_type = "none"),
  list(name = "Hetero_tau0.10", type = "hetero", k = 15, true_mu = 0.3, tau2 = 0.10, bias_type = "none"),
  list(name = "Hetero_tau0.50", type = "hetero", k = 15, true_mu = 0.3, tau2 = 0.50, bias_type = "none"),
  list(name = "Hetero_tau1.00", type = "hetero", k = 15, true_mu = 0.3, tau2 = 1.00, bias_type = "none"),

  # Publication bias — significance-based
  list(name = "PubBias_k10_mild",   type = "pubbias", k = 10, true_mu = 0.3, tau2 = 0.05, bias_type = "significance", bias_strength = 0.5),
  list(name = "PubBias_k10_strong", type = "pubbias", k = 10, true_mu = 0.3, tau2 = 0.05, bias_type = "significance", bias_strength = 0.8),
  list(name = "PubBias_k20_strong", type = "pubbias", k = 20, true_mu = 0.3, tau2 = 0.05, bias_type = "significance", bias_strength = 0.8),
  list(name = "PubBias_k20_severe", type = "pubbias", k = 20, true_mu = 0.3, tau2 = 0.05, bias_type = "significance", bias_strength = 0.95),

  # One-sided publication bias
  list(name = "OneSided_k10", type = "onesided", k = 10, true_mu = 0.3, tau2 = 0.05, bias_type = "onesided", bias_strength = 0.7),
  list(name = "OneSided_k20", type = "onesided", k = 20, true_mu = 0.3, tau2 = 0.10, bias_type = "onesided", bias_strength = 0.7),

  # Combined: high heterogeneity + publication bias
  list(name = "Combined_hetero_bias", type = "combined", k = 20, true_mu = 0.3, tau2 = 0.30, bias_type = "significance", bias_strength = 0.8),

  # Null effect + publication bias (worst case)
  list(name = "Null_PubBias", type = "null_bias", k = 20, true_mu = 0.0, tau2 = 0.05, bias_type = "significance", bias_strength = 0.8)
)

cat(sprintf("Running %d scenarios x %d reps = %d total iterations\n", length(scenarios), N_SIM, length(scenarios) * N_SIM))

# Run simulation
all_results <- list()
for (si in seq_along(scenarios)) {
  sc <- scenarios[[si]]
  cat(sprintf("[%d/%d] %s (k=%d, tau2=%.2f, bias=%s)\n", si, length(scenarios), sc$name, sc$k, sc$tau2, sc$bias_type))
  t0 <- Sys.time()

  for (r in 1:N_SIM) {
    dat <- generate_meta_data(sc$k, sc$true_mu, sc$tau2, sc$bias_type,
                              ifelse(is.null(sc$bias_strength), 0, sc$bias_strength))

    # REML
    reml_res <- tryCatch({
      fit <- rma(dat$yi, dat$vi, method = "REML")
      list(est = as.numeric(fit$beta), se = fit$se, ci_lb = fit$ci.lb, ci_ub = fit$ci.ub)
    }, error = function(e) NULL)

    # HKSJ
    hksj_res <- tryCatch({
      fit <- rma(dat$yi, dat$vi, method = "REML", test = "knha")
      list(est = as.numeric(fit$beta), se = fit$se, ci_lb = fit$ci.lb, ci_ub = fit$ci.ub)
    }, error = function(e) NULL)

    # TGEP
    tgep_res <- tryCatch({
      fit <- tgep_meta(dat$yi, dat$vi, n_boot = 0)
      list(est = fit$estimate, se = fit$se, ci_lb = fit$ci_lb, ci_ub = fit$ci_ub)
    }, error = function(e) NULL)

    for (method_info in list(
      list(name = "REML", res = reml_res),
      list(name = "HKSJ", res = hksj_res),
      list(name = "TGEP", res = tgep_res)
    )) {
      if (!is.null(method_info$res)) {
        all_results[[length(all_results) + 1]] <- data.frame(
          Scenario = sc$name, Type = sc$type, k = sc$k, tau2 = sc$tau2,
          bias_type = sc$bias_type, true_mu = sc$true_mu,
          Method = method_info$name,
          Estimate = method_info$res$est,
          SE = method_info$res$se,
          CI_Lower = method_info$res$ci_lb,
          CI_Upper = method_info$res$ci_ub,
          Covered = (sc$true_mu >= method_info$res$ci_lb & sc$true_mu <= method_info$res$ci_ub),
          k_final = dat$k_final,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  %.0fs\n", dt))
}

results_df <- do.call(rbind, all_results)

# Compute summary statistics
summary_df <- aggregate(
  cbind(Estimate, Covered) ~ Scenario + Type + k + tau2 + bias_type + true_mu + Method,
  data = results_df,
  FUN = function(x) c(mean = mean(x))
)
# Fix column names after aggregate
colnames(summary_df)[8:9] <- c("Mean_Est", "Coverage")

# Compute bias, RMSE, MCSE separately
compute_metrics <- function(df) {
  do.call(rbind, lapply(split(df, paste(df$Scenario, df$Method)), function(sub) {
    bias <- mean(sub$Estimate - sub$true_mu)
    rmse <- sqrt(mean((sub$Estimate - sub$true_mu)^2))
    coverage <- mean(sub$Covered)
    mcse_cov <- sqrt(coverage * (1 - coverage) / nrow(sub))
    mean_se <- mean(sub$SE)
    emp_se <- sd(sub$Estimate)
    mean_k <- mean(sub$k_final)
    data.frame(
      Scenario = sub$Scenario[1], Type = sub$Type[1], k = sub$k[1], tau2 = sub$tau2[1],
      bias_type = sub$bias_type[1], true_mu = sub$true_mu[1], Method = sub$Method[1],
      Mean_Bias = bias, RMSE = rmse, Coverage = coverage, MCSE_Coverage = mcse_cov,
      Mean_SE = mean_se, Empirical_SE = emp_se, Mean_k_final = mean_k,
      N_sim = nrow(sub),
      stringsAsFactors = FALSE
    )
  }))
}

metrics <- compute_metrics(results_df)
rownames(metrics) <- NULL

# Save
dir.create("C:/Models/TGEP_Development/output", showWarnings = FALSE)
write.csv(metrics, "C:/Models/TGEP_Development/output/simulation_results.csv", row.names = FALSE)
write.csv(results_df, "C:/Models/TGEP_Development/output/simulation_raw.csv", row.names = FALSE)

# Print summary
cat("\n=== SIMULATION SUMMARY ===\n")
for (sc_name in unique(metrics$Scenario)) {
  cat(sprintf("\n--- %s ---\n", sc_name))
  sub <- metrics[metrics$Scenario == sc_name, ]
  for (i in 1:nrow(sub)) {
    cat(sprintf("  %s: Bias=%.4f RMSE=%.4f Coverage=%.3f (MCSE=%.4f) SE_ratio=%.3f\n",
                sub$Method[i], sub$Mean_Bias[i], sub$RMSE[i], sub$Coverage[i],
                sub$MCSE_Coverage[i], sub$Mean_SE[i] / sub$Empirical_SE[i]))
  }
}

cat(sprintf("\nTotal time: %.1f min\n", as.numeric(difftime(Sys.time(), t_global, units = "mins"))))
cat("Results saved to output/simulation_results.csv\n")
