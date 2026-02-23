#' TGEP: Comprehensive Methodology Evaluation (Updated)
library(metafor)
library(data.table)

# Source core TGEP logic from local directory
source("C:/Models/TGEP_Development/TGEP.R")

# Simulation Engine
generate_data <- function(k = 10, true_effect = 0.3, tau2 = 0.05, bias = FALSE) {
  theta <- rnorm(k, true_effect, sqrt(tau2))
  vi <- rexp(k, 10) + 0.01
  yi <- rnorm(k, theta, sqrt(vi))
  if (bias) {
    p <- 2 * (1 - pnorm(abs(yi / sqrt(vi))))
    keep <- (p < 0.05) | (runif(k) < 0.2)
    if(sum(keep) < 3) keep <- rank(p) <= 3
    yi <- yi[keep]; vi <- vi[keep]
  }
  list(yi = yi, vi = vi, k = length(yi))
}

run_sim_scenarios <- function(n_reps = 50) {
  scenarios <- list(
    list(name = "Small k, Low Hetero", k = 5, tau2 = 0.02, bias = FALSE),
    list(name = "Med k, High Hetero", k = 20, tau2 = 0.20, bias = FALSE),
    list(name = "Med k, Pub Bias", k = 20, tau2 = 0.05, bias = TRUE)
  )
  methods <- list(
    REML = function(y, v) { 
      f <- rma(y, v, method="REML");
      list(est = as.numeric(f$beta), ci_lb = f$ci.lb, ci_ub = f$ci.ub) 
    },
    TGEP = function(y, v) {
      # Use the main tgep_meta function from TGEP.R
      f <- tgep_meta(y, v, n_boot = 0) # n_boot=0 for fast simulation
      list(est = f$estimate, ci_lb = f$ci_lb, ci_ub = f$ci_ub)
    }
  )
  all_results <- list()
  for (sc in scenarios) {
    cat("Running Scenario:", sc$name, "\n")
    reps_data <- list()
    for (r in 1:n_reps) {
      dat <- generate_data(k = sc$k, true_effect = 0.3, tau2 = sc$tau2, bias = sc$bias)
      for (m in names(methods)) {
        res <- tryCatch(methods[[m]](dat$yi, dat$vi), error = function(e) NULL)
        if(!is.null(res)) {
          reps_data[[length(reps_data)+1]] <- data.table(
            scenario = sc$name, method = m, 
            estimate = res$est, 
            covered = (0.3 >= res$ci_lb & 0.3 <= res$ci_ub)
          )
        }
      }
    }
    all_results[[sc$name]] <- rbindlist(reps_data)
  }
  rbindlist(all_results)
}

sink("C:/Models/TGEP_Development/Full_Methodology_Report_Updated.txt")
cat("TGEP COMPREHENSIVE EVALUATION REPORT (UPDATED)\n")
cat("==============================================\n\n")

# Simulation Study
cat("SIMULATION STUDY (Bias, RMSE, Coverage)\n")
cat("------------------------------------------\n")
sim_res <- run_sim_scenarios(n_reps = 100)
summary_stats <- sim_res[, .(
  Mean_Bias = mean(estimate - 0.3),
  RMSE = sqrt(mean((estimate - 0.3)^2)),
  Coverage = mean(covered)
), by = .(scenario, method)]
print(summary_stats)
sink()

cat("Analysis complete. Report saved to C:/Models/TGEP_Development/Full_Methodology_Report_Updated.txt\n")
