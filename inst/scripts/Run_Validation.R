# Run Validation for TGEP Methodology (Updated)
library(metafor)
source("C:/Models/TGEP_Development/TGEP.R")

cat("=== TGEP Validation Run (Fully Adaptive) ===\n")

# 1. Load Data
# Using local path if possible, or same absolute if required by data location
load("C:/Users/user/OneDrive - NHS/Documents/Pairwise70/data/CD000028_pub4_data.rda")
df <- CD000028_pub4_data

# 2. Extract a specific meta-analysis
ma_data <- df[df$Analysis.number == 1, ]
cat(sprintf("Selected Analysis: %s (k = %d)\n", ma_data$Analysis.name[1], nrow(ma_data)))

# 3. Calculate effect sizes (Log OR)
dat <- escalc(measure="OR", 
              ai=Experimental.cases, n1i=Experimental.N,
              ci=Control.cases, n2i=Control.N, 
              data=ma_data)

dat <- dat[!is.na(dat$yi) & !is.na(dat$vi), ][1:5, ]
cat(sprintf("Studies after filtering: %d\n", nrow(dat)))

# 4. Run Comparison (with bootstrap)
tgep_res <- tgep_meta(dat$yi, dat$vi, n_boot = 0)
results <- compare_tgep(dat$yi, dat$vi)

# 5. Print Results
print(results)

sink("C:/Models/TGEP_Development/Validation_Results_Updated.txt")
cat("TGEP VALIDATION REPORT (UPDATED)\n")
cat("================================\n\n")
cat(sprintf("Dataset: CD000028_pub4_data\n"))
cat(sprintf("Analysis: %s\n", ma_data$Analysis.name[1]))
cat(sprintf("Studies: %d\n\n", nrow(dat)))

cat("ENSEMBLE SUMMARY:\n")
cat("----------------\n")
cat(sprintf("Final TGEP Estimate: %.4f (SE: %.4f)\n", tgep_res$estimate, tgep_res$se))
cat(sprintf("95%% CI: [%.4f, %.4f]\n", tgep_res$ci_lb, tgep_res$ci_ub))
cat(sprintf("P-value: %.4f\n\n", tgep_res$pvalue))

cat("GUARD PERFORMANCE:\n")
cat("------------------\n")
diag_df <- data.frame(
  Guard = names(tgep_res$guard_estimates),
  Estimate = as.numeric(tgep_res$guard_estimates),
  Weight = as.numeric(tgep_res$ensemble_weights),
  CV_Error = as.numeric(tgep_res$cv_scores)
)
print(diag_df)

cat("\nMETHOD COMPARISON:\n")
cat("------------------\n")
print(results)
sink()

cat("\nValidation complete. Results saved to C:/Models/TGEP_Development/Validation_Results_Updated.txt\n")
