# Run Validation for TGEP Methodology (Updated)
library(metafor)

# Use system.file() when installed as package, or source from project root
if (requireNamespace("TGEP", quietly = TRUE)) {
  library(TGEP)
} else {
  source(file.path(getwd(), "R", "TGEP.R"))
}

cat("=== TGEP Validation Run (Fully Adaptive) ===\n")

# 1. Load Data
# Expects CD000028_pub4_data.rda in the working directory or a user-specified path.
# To use: place the .rda file in getwd(), or set DATA_PATH before sourcing this script.
data_path <- Sys.getenv("TGEP_DATA_PATH", unset = file.path(getwd(), "CD000028_pub4_data.rda"))
if (!file.exists(data_path)) stop("Data file not found: ", data_path,
                                   "\nSet TGEP_DATA_PATH environment variable or place file in working directory.")
load(data_path)
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

output_path <- file.path(getwd(), "output", "Validation_Results_Updated.txt")
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
sink(output_path)
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
  Guard = names(tgep_res$guards),
  Estimate = as.numeric(tgep_res$guards),
  Weight = as.numeric(tgep_res$weights)
)
print(diag_df)

cat("\nMETHOD COMPARISON:\n")
cat("------------------\n")
print(results)
sink()

cat(sprintf("\nValidation complete. Results saved to %s\n", output_path))
