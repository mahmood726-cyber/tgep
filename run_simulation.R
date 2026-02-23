#!/usr/bin/env Rscript
# TGEP Comprehensive Simulation Study (ADEMP Framework)
# Author: Mahmood Ul Hassan
# Date: 2026-02-23
# Version: 2.0 — fixed DGM, 500 reps, added trim-and-fill & PET-PEESE
#
# ADEMP:
#   Aims: Compare TGEP vs REML, HKSJ, Trim-and-Fill, PET-PEESE under bias
#   Data-generating mechanisms: Normal RE with realistic variances
#   Estimands: Unconditional population mean mu (pre-selection)
#   Methods: REML, HKSJ, TGEP (T=1.0, n_boot=0), Trim-Fill, PET-PEESE
#   Performance: Bias, RMSE, Coverage, MCSE, SE ratio

library(metafor)
source(file.path(getwd(), "R", "TGEP.R"))

set.seed(20260223)
N_SIM <- 500
t_global <- Sys.time()

# Data-generating mechanism with REALISTIC variances
# For SMD with n_per_arm ~ U(20,100): vi ~ 2/n + d^2/(2*2n) ~ 0.02-0.1
generate_meta_data <- function(k, true_mu, tau2, bias_type = "none", bias_strength = 0.8) {
  n_per_arm <- sample(20:100, k, replace = TRUE)
  vi <- 2 / n_per_arm + true_mu^2 / (2 * 2 * n_per_arm)  # SMD variance approximation
  theta_i <- rnorm(k, true_mu, sqrt(tau2))
  yi <- rnorm(k, theta_i, sqrt(vi))

  if (bias_type == "significance") {
    p <- 2 * (1 - pnorm(abs(yi / sqrt(vi))))
    prob_pub <- ifelse(p < 0.05, 1.0, 1 - bias_strength)
    keep <- runif(k) < prob_pub
    if (sum(keep) < 3) keep[order(p)][1:3] <- TRUE
    yi <- yi[keep]; vi <- vi[keep]
  } else if (bias_type == "onesided") {
    prob_pub <- ifelse(yi > 0, 1.0, 1 - bias_strength)
    keep <- runif(k) < prob_pub
    if (sum(keep) < 3) keep[order(-yi)][1:3] <- TRUE
    yi <- yi[keep]; vi <- vi[keep]
  }

  list(yi = yi, vi = vi, k_final = length(yi))
}

# PET-PEESE helper
pet_peese <- function(yi, vi) {
  se <- sqrt(vi)
  # PET: regress on SE
  pet_fit <- tryCatch(rma(yi, vi, mods = ~ se, method = "REML"), error = function(e) NULL)
  if (is.null(pet_fit)) return(NULL)
  pet_est <- as.numeric(pet_fit$beta[1])
  pet_pval <- pet_fit$pval[1]
  # If PET significant, use PEESE (regress on vi)
  if (pet_pval < 0.10) {
    peese_fit <- tryCatch(rma(yi, vi, mods = ~ vi, method = "REML"), error = function(e) NULL)
    if (!is.null(peese_fit)) {
      return(list(est = as.numeric(peese_fit$beta[1]), se = peese_fit$se[1],
                  ci_lb = peese_fit$ci.lb[1], ci_ub = peese_fit$ci.ub[1]))
    }
  }
  list(est = pet_est, se = pet_fit$se[1],
       ci_lb = pet_fit$ci.lb[1], ci_ub = pet_fit$ci.ub[1])
}

# Define scenarios
scenarios <- list(
  # Null scenarios
  list(name = "Null_k10_tau0.05", type = "null", k = 10, true_mu = 0.0, tau2 = 0.05, bias_type = "none"),
  list(name = "Null_k10_tau0.20", type = "null", k = 10, true_mu = 0.0, tau2 = 0.20, bias_type = "none"),
  # Baseline (no bias) — varying k
  list(name = "Base_k5",  type = "base", k = 5,  true_mu = 0.3, tau2 = 0.05, bias_type = "none"),
  list(name = "Base_k10", type = "base", k = 10, true_mu = 0.3, tau2 = 0.05, bias_type = "none"),
  list(name = "Base_k20", type = "base", k = 20, true_mu = 0.3, tau2 = 0.05, bias_type = "none"),
  # Base_k30 removed for computational efficiency; k=5,10,20 cover the range
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
  # Small k + bias
  list(name = "PubBias_k5_strong",  type = "pubbias", k = 5,  true_mu = 0.3, tau2 = 0.05, bias_type = "significance", bias_strength = 0.8),
  # One-sided publication bias
  list(name = "OneSided_k10", type = "onesided", k = 10, true_mu = 0.3, tau2 = 0.05, bias_type = "onesided", bias_strength = 0.7),
  list(name = "OneSided_k20", type = "onesided", k = 20, true_mu = 0.3, tau2 = 0.10, bias_type = "onesided", bias_strength = 0.7),
  # Combined: high heterogeneity + publication bias
  list(name = "Combined_hetero_bias", type = "combined", k = 20, true_mu = 0.3, tau2 = 0.30, bias_type = "significance", bias_strength = 0.8),
  # Null effect + publication bias
  list(name = "Null_PubBias", type = "null_bias", k = 20, true_mu = 0.0, tau2 = 0.05, bias_type = "significance", bias_strength = 0.8)
)

cat(sprintf("Running %d scenarios x %d reps = %d iterations\n", length(scenarios), N_SIM, length(scenarios) * N_SIM))

methods_list <- c("REML", "HKSJ", "TGEP", "TrimFill", "PETPEESE")

all_results <- list()
for (si in seq_along(scenarios)) {
  sc <- scenarios[[si]]
  cat(sprintf("[%d/%d] %s (k=%d, tau2=%.2f, bias=%s)\n", si, length(scenarios), sc$name, sc$k, sc$tau2, sc$bias_type))
  t0 <- Sys.time()
  n_discard <- 0

  for (r in 1:N_SIM) {
    dat <- generate_meta_data(sc$k, sc$true_mu, sc$tau2, sc$bias_type,
                              ifelse(is.null(sc$bias_strength), 0, sc$bias_strength))

    # Run all methods; if any fail, discard entire iteration
    reml_res <- tryCatch({
      fit <- rma(dat$yi, dat$vi, method = "REML")
      list(est = as.numeric(fit$beta), se = fit$se, ci_lb = fit$ci.lb, ci_ub = fit$ci.ub)
    }, error = function(e) NULL)

    hksj_res <- tryCatch({
      fit <- rma(dat$yi, dat$vi, method = "REML", test = "knha")
      list(est = as.numeric(fit$beta), se = fit$se, ci_lb = fit$ci.lb, ci_ub = fit$ci.ub)
    }, error = function(e) NULL)

    tgep_res <- tryCatch({
      fit <- tgep_meta(dat$yi, dat$vi, n_boot = 0)
      list(est = fit$estimate, se = fit$se, ci_lb = fit$ci_lb, ci_ub = fit$ci_ub)
    }, error = function(e) NULL)

    tf_res <- tryCatch({
      fit <- rma(dat$yi, dat$vi, method = "REML")
      tf <- trimfill(fit)
      list(est = as.numeric(tf$beta), se = tf$se, ci_lb = tf$ci.lb, ci_ub = tf$ci.ub)
    }, error = function(e) NULL)

    pp_res <- tryCatch(pet_peese(dat$yi, dat$vi), error = function(e) NULL)

    # Discard if any core method fails
    if (is.null(reml_res) || is.null(hksj_res) || is.null(tgep_res)) {
      n_discard <- n_discard + 1
      next
    }

    for (method_info in list(
      list(name = "REML", res = reml_res),
      list(name = "HKSJ", res = hksj_res),
      list(name = "TGEP", res = tgep_res),
      list(name = "TrimFill", res = tf_res),
      list(name = "PETPEESE", res = pp_res)
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
  cat(sprintf("  %.0fs (discarded %d)\n", dt, n_discard))
}

results_df <- do.call(rbind, all_results)

# Compute metrics
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
      Mean_SE = mean_se, Empirical_SE = emp_se, SE_Ratio = mean_se / emp_se,
      Mean_k_final = mean_k, N_sim = nrow(sub),
      stringsAsFactors = FALSE
    )
  }))
}

metrics <- compute_metrics(results_df)
rownames(metrics) <- NULL

# Save
dir.create(file.path(getwd(), "output"), showWarnings = FALSE)
write.csv(metrics, file.path(getwd(), "output", "simulation_results.csv"), row.names = FALSE)

# Print summary
cat("\n=== SIMULATION SUMMARY ===\n")
for (sc_name in unique(metrics$Scenario)) {
  cat(sprintf("\n--- %s ---\n", sc_name))
  sub <- metrics[metrics$Scenario == sc_name, ]
  for (i in 1:nrow(sub)) {
    cat(sprintf("  %s: Bias=%.4f RMSE=%.4f Cov=%.3f (MCSE=%.4f) SE_r=%.3f N=%d\n",
                sub$Method[i], sub$Mean_Bias[i], sub$RMSE[i], sub$Coverage[i],
                sub$MCSE_Coverage[i], sub$SE_Ratio[i], sub$N_sim[i]))
  }
}

cat(sprintf("\nTotal time: %.1f min\n", as.numeric(difftime(Sys.time(), t_global, units = "mins"))))
cat("Results saved to output/simulation_results.csv\n")
