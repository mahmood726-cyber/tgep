# Triple-Guard Ensemble Pooling: A Frequentist Ensemble Framework for Bias-Reduced Point Estimation in Meta-Analysis

**Mahmood Ul Hassan**^1

^1 Independent Researcher

**Corresponding author:** Mahmood Ul Hassan (mahmood.hassan@example.com)

---

## Abstract

**Background:** Standard random-effects meta-analysis using restricted maximum likelihood (REML) produces biased point estimates under publication bias, yet existing corrections such as PET-PEESE can overcorrect when bias is absent. A method that reduces bias without strong parametric assumptions about the selection mechanism would complement the existing toolkit.

**Methods:** We propose Triple-Guard Ensemble Pooling (TGEP), a frequentist ensemble framework integrating three estimators — Grey Relational Meta-Analysis (GRMA), Winsorized Robust Detection (WRD), and Significance-Weighted Adjustment (SWA) — via leave-one-out cross-validation (LOO-CV) stacking. A Monte Carlo simulation (1000 iterations x 18 scenarios) compared TGEP against REML, Hartung-Knapp-Sidik-Jonkman (HKSJ), trim-and-fill, and PET-PEESE across k = 5-20, tau-squared = 0.01-1.00, and significance-based and one-sided publication bias.

**Results:** Under one-sided bias (k = 20, tau-squared = 0.10), TGEP reduced absolute bias by 24% relative to REML (0.067 vs. 0.088) with 12% lower root mean squared error (RMSE). Under significance-based bias, reductions ranged from less than 1% to 23%. In unbiased scenarios, TGEP maintained RMSE comparable to or lower than REML. However, its analytical SE underestimated empirical variability (ratio approximately 0.75-1.00), yielding coverage of approximately 87-93% in unbiased settings versus HKSJ's 93-95%. PET-PEESE achieved the best coverage under significance-based bias (approximately 87-91%) but overcorrected when bias was absent (RMSE approximately 4-5x larger).

**Conclusions:** TGEP provides less biased point estimates than REML under publication bias, particularly one-sided selection, without PET-PEESE's overcorrection risk. Its coverage is limited by SE underestimation. We recommend TGEP as a complementary sensitivity analysis alongside HKSJ for inference and PET-PEESE for bias assessment. The method is available as an open-source R package.

**Keywords:** meta-analysis; publication bias; robust estimation; ensemble methods; cross-validation; sensitivity analysis

---

## Introduction

Random-effects meta-analysis is the standard approach for synthesizing treatment effects across studies [1]. The most widely used estimator, restricted maximum likelihood (REML) with inverse-variance weighting [2], assumes that the observed studies are a representative sample from a population of studies. This assumption is routinely violated in practice: publication bias causes selective reporting of statistically significant results [3], and influential outliers can distort pooled estimates even when heterogeneity is properly modeled [4].

The consequences of these violations are well documented. Under publication bias, inverse-variance pooled estimates are systematically inflated, and confidence interval coverage drops substantially below the nominal 95% level [5]. Several approaches have been proposed. The Hartung-Knapp-Sidik-Jonkman (HKSJ) correction [6] addresses the underestimation of standard errors in standard random-effects models by using a t-distribution reference, providing improved coverage but not correcting the biased point estimate. Trim-and-fill [7] is widely used but relies on the assumption of funnel plot symmetry. PET-PEESE [8] uses meta-regression on standard errors or variances to extrapolate to the zero-bias intercept, but can overcorrect when publication bias is absent, inflating RMSE [9]. Selection models [10] require specification of the selection function. Bayesian model averaging approaches such as RoBMA [11] integrate over multiple bias models but require prior specification and are computationally intensive.

These methods make different assumptions about the bias mechanism, and no single method dominates across all scenarios [5]. This motivates an ensemble approach that combines multiple bias-correction strategies, adapting to the specific pattern present in the data.

We propose Triple-Guard Ensemble Pooling (TGEP), a frequentist ensemble framework that integrates three specialized "guard" estimators — each designed to be robust to a different source of distortion — using leave-one-out cross-validation (LOO-CV) stacking [12] with softmax temperature weighting. The guards are:

1. **Grey Relational Meta-Analysis (GRMA)**: Uses grey relational analysis [13] with robust scaling to down-weight studies that are simultaneously discordant in both effect size and precision, providing resistance to outliers and small-study effects.
2. **Winsorized Robust Detection (WRD)**: Applies Winsorization [14] to z-score residuals from an initial REML fit, bounding the influence of extreme observations while preserving the bulk of the data.
3. **Significance-Weighted Adjustment (SWA)**: Implements inverse-probability-of-selection weighting by up-weighting non-significant studies, which are under-represented under publication bias.

The LOO-CV stacking mechanism learns which guards are most predictive for the observed data, yielding an adaptive ensemble. We note at the outset that LOO-CV optimizes prediction of observed (potentially biased) data rather than the unobserved population parameter; consequently, TGEP's primary contribution is to point estimation rather than confidence interval coverage. We evaluate TGEP alongside REML, HKSJ, trim-and-fill, and PET-PEESE in a comprehensive Monte Carlo simulation study following the ADEMP reporting framework [15].

## Methods

### TGEP Framework

#### Guard 1: Grey Relational Meta-Analysis (GRMA)

The GRMA guard uses grey relational analysis to assign study weights based on the joint proximity of each study's effect size and log-precision to reference values, replacing inverse-variance weights with weights that are inherently resistant to outliers.

For K studies with effect sizes y_i and sampling variances v_i, define the precision w_i = 1/v_i and log-precision l_i = log(w_i + 1). Apply robust scaling to both sequences using the 5th and 95th percentiles:

x_i^(j) = clip((z_i^(j) - Q_0.05^(j)) / (Q_0.95^(j) - Q_0.05^(j)), 0, 1),  j in {y, l}

The reference sequence uses the median effect size and maximum precision. The grey relational coefficient for study i across dimension j is:

xi_i^(j) = (Delta_min + zeta * Delta_max) / (delta_i^(j) + zeta * Delta_max)

where delta_i^(j) = |x_i^(j) - x_ref^(j)|, Delta_min and Delta_max are the global minimum and maximum deviations, and zeta = 0.5 is the distinguishing coefficient. When all studies are identical (Delta_max < 10^-12), equal weights 1/K are assigned to avoid division by zero. The GRMA weight for study i is the normalized mean grey relational coefficient: w_i^GRMA = xi_bar_i / sum_j(xi_bar_j).

The GRMA pooled estimate is mu_hat_GRMA = sum_i(w_i^GRMA * y_i) with variance V_GRMA = sum_i((w_i^GRMA)^2 * v_i). Note that V_GRMA uses only the within-study variances v_i rather than v_i + tau_hat_squared; because GRMA is not a random-effects estimator per se but a data-adaptive weighting scheme, incorporating the between-study heterogeneity into the GRMA variance would conflate two distinct sources of uncertainty. This design choice contributes to the analytical SE underestimation documented in the simulation study (see Discussion, SE Underestimation).

#### Guard 2: Winsorized Robust Detection (WRD)

The WRD guard fits an initial REML random-effects model to obtain mu_hat_REML and tau_hat_squared, then computes standardized residuals:

z_i = (y_i - mu_hat_REML) / sqrt(v_i + tau_hat_squared)

Residuals exceeding a threshold c = 2.5 in absolute value are Winsorized:

z_i* = sign(z_i) * min(|z_i|, c)

The Winsorized effect sizes y_i* = mu_hat_REML + z_i* * sqrt(v_i + tau_hat_squared) are then pooled using inverse-variance weights w_i = 1/(v_i + tau_hat_squared):

mu_hat_WRD = sum_i(w_i * y_i*) / sum_i(w_i),  V_WRD = 1 / sum_i(w_i)

#### Guard 3: Significance-Weighted Adjustment (SWA)

The SWA guard implements a form of inverse-probability-of-selection weighting to correct for publication bias. Under the standard selection model, studies with statistically significant results (p < 0.05) are more likely to be published. Non-significant studies that do appear in the literature are therefore under-represented and should receive greater weight.

From an initial REML fit yielding mu_hat and tau_hat_squared, compute within-study p-values using only the marginal (within-study) variance:

p_i = 2[1 - Phi(|y_i| / sqrt(v_i))]

This uses sqrt(v_i) rather than sqrt(v_i + tau_hat_squared) because publication decisions are made based on each study's individual statistical significance, which depends on the within-study standard error alone. The SWA weights are:

w_i^SWA = w_i^IV / s_i,  where s_i = 1.0 if p_i < 0.05; 0.4 if p_i >= 0.05

where w_i^IV = 1/(v_i + tau_hat_squared). This up-weights non-significant studies by a factor of 1/0.4 = 2.5. After normalization, mu_hat_SWA = sum_i(w_tilde_i * y_i) and V_SWA = sum_i(w_tilde_i^2 * (v_i + tau_hat_squared)).

#### LOO-CV Stacking Ensemble

The three guard estimates are combined using leave-one-out cross-validation (LOO-CV) stacking. For each study i = 1, ..., K:

1. Remove study i from the dataset.
2. Compute each guard's estimate mu_hat_g(-i) from the remaining K-1 studies.
3. Compute the precision-weighted squared prediction error: e_ig = (y_i - mu_hat_g(-i))^2 / v_i.

The mean LOO-CV error for guard g is e_bar_g = (1/K) * sum_i(e_ig). Guard weights are computed via the softmax function with temperature parameter T:

alpha_g = exp(-e_tilde_g / T) / sum_g'(exp(-e_tilde_g' / T))

where e_tilde_g = (e_bar_g - min_g(e_bar_g)) / (max_g(e_bar_g) - min_g(e_bar_g) + epsilon) is the min-max normalized error. The temperature T = 1.0 (default) provides a balance between aggressive selection of the best guard (T -> 0) and equal weighting (T -> infinity).

The TGEP ensemble estimate is:

mu_hat_TGEP = sum_g(alpha_g * mu_hat_g)

Because the three guards are computed from the same data and are therefore strongly positively correlated, we estimate the ensemble standard error under the assumption of perfect guard correlation:

SE_TGEP = sum_g(alpha_g * sqrt(V_g))

This yields V_TGEP = SE_TGEP^2. This yields wider confidence intervals than the independent-guard formula sum_g(alpha_g^2 * V_g), but neither formula accounts for the additional variability introduced by the data-dependent weight estimation; as our simulation study demonstrates, the net result still underestimates the true variability of the ensemble estimator. Confidence intervals use the normal approximation: mu_hat_TGEP +/- z_(1-alpha/2) * SE_TGEP. A nonparametric bootstrap (B = 200 by default) is also available for variance estimation that accounts for the additional uncertainty from the LOO-CV weight selection.

For meta-analyses with K < 3 studies, LOO-CV is uninformative, and TGEP falls back to equal guard weighting (alpha_g = 1/3 for all g).

### Simulation Study

We followed the ADEMP reporting framework [15] for the design and reporting of this simulation study.

#### Aims

To compare the bias, root mean squared error (RMSE), coverage probability, and SE calibration (SE ratio) of five methods — REML, HKSJ, TGEP, trim-and-fill, and PET-PEESE — under varying numbers of studies, heterogeneity levels, and publication bias mechanisms.

#### Estimands

The target estimand was the unconditional population mean mu (the pre-selection true effect), set to either 0.0 (null) or 0.3 (non-null, representing a small-to-moderate standardized mean difference).

#### Data-Generating Mechanisms

We simulated random-effects meta-analysis data for standardized mean differences (SMD). For each of K studies:

1. Sample the per-arm sample size: n_i ~ Uniform(20, 100)
2. Compute the within-study variance using the large-sample SMD variance approximation [1]: v_i = 2/n_i + mu^2/(4*n_i), yielding realistic variances in the range 0.02-0.10 for the chosen effect sizes.
3. Sample the true study-specific effect: theta_i ~ N(mu, tau^2)
4. Sample the observed effect: y_i ~ N(theta_i, v_i)

Publication bias was simulated via two mechanisms:

- **Significance-based**: Studies with p < 0.05 (two-sided z-test using within-study SE) were always published; non-significant studies were published with probability 1 - beta, where beta in {0.50, 0.80, 0.95} controlled bias severity.
- **One-sided**: Studies with effects in the expected direction (y_i > 0) were always published; negative effects were published with probability 1 - beta.

A minimum of 3 studies was enforced in all scenarios to ensure estimability. When selection reduced the sample below 3, the 3 studies with the smallest p-values (significance-based) or largest effects (one-sided) were retained.

#### Methods Compared

1. **REML**: Standard restricted maximum likelihood random-effects meta-analysis [2] with Wald-type confidence intervals.
2. **HKSJ**: REML with the Hartung-Knapp-Sidik-Jonkman adjustment [6,16,17] using a t_(K-1) reference distribution for confidence intervals.
3. **TGEP**: Triple-Guard Ensemble Pooling with T = 1.0 and analytical variance (no bootstrap, for computational efficiency in the simulation).
4. **Trim-and-fill**: The Duval-Tweedie trim-and-fill method [7] applied to the REML fit using the L0 estimator (the metafor default), which estimates the number of "missing" studies and imputes them to restore funnel plot symmetry.
5. **PET-PEESE**: The precision-effect test / precision-effect estimate with standard error (PET-PEESE) [8]. PET regresses effect sizes on standard errors; if the PET intercept is significant at p < 0.10, PEESE (regression on variances) is used instead.

If any of the three core methods (REML, HKSJ, TGEP) failed for a given iteration, that entire iteration was discarded. Trim-and-fill and PET-PEESE failures were treated as missing for that iteration only. Across all 18 scenarios with 1000 replications each, zero iterations were discarded.

#### Performance Measures

- **Bias**: Mean difference between estimated and true mu
- **RMSE**: Root mean squared error, sqrt(mean((mu_hat - mu)^2))
- **Coverage**: Proportion of 95% confidence intervals containing the true mu
- **MCSE of coverage**: Monte Carlo standard error, sqrt(C*(1-C)/N_sim)
- **SE ratio**: Mean model-based SE divided by empirical SE; values near 1.0 indicate well-calibrated SE estimation; values below 1.0 indicate anti-conservative SE (confidence intervals too narrow)

#### Scenarios

Table 1 summarizes the 18 simulation scenarios. Each scenario was replicated N_sim = 1000 times with a fixed seed (20260223) for reproducibility.

**Table 1. Simulation Scenarios**

| # | Scenario | Type | K | tau^2 | Bias | Strength |
|---|----------|------|---|-------|------|----------|
| 1 | Null_k10_tau0.05 | Null | 10 | 0.05 | None | -- |
| 2 | Null_k10_tau0.20 | Null | 10 | 0.20 | None | -- |
| 3 | Base_k5 | Baseline | 5 | 0.05 | None | -- |
| 4 | Base_k10 | Baseline | 10 | 0.05 | None | -- |
| 5 | Base_k20 | Baseline | 20 | 0.05 | None | -- |
| 6 | Hetero_tau0.01 | Heterogeneity | 15 | 0.01 | None | -- |
| 7 | Hetero_tau0.10 | Heterogeneity | 15 | 0.10 | None | -- |
| 8 | Hetero_tau0.50 | Heterogeneity | 15 | 0.50 | None | -- |
| 9 | Hetero_tau1.00 | Heterogeneity | 15 | 1.00 | None | -- |
| 10 | PubBias_k10_mild | Pub. bias | 10 | 0.05 | Significance | 0.50 |
| 11 | PubBias_k10_strong | Pub. bias | 10 | 0.05 | Significance | 0.80 |
| 12 | PubBias_k20_strong | Pub. bias | 20 | 0.05 | Significance | 0.80 |
| 13 | PubBias_k20_severe | Pub. bias | 20 | 0.05 | Significance | 0.95 |
| 14 | PubBias_k5_strong | Small-k bias | 5 | 0.05 | Significance | 0.80 |
| 15 | OneSided_k10 | One-sided | 10 | 0.05 | One-sided | 0.70 |
| 16 | OneSided_k20 | One-sided | 20 | 0.10 | One-sided | 0.70 |
| 17 | Combined_hetero_bias | Combined | 20 | 0.30 | Significance | 0.80 |
| 18 | Null_PubBias | Null + bias | 20 | 0.05 | Significance | 0.80 |

### Empirical validation

We applied TGEP and REML to 15 datasets from Cochrane systematic reviews spanning diverse clinical areas. These datasets ranged from K = 5 to K = 88 studies and used log odds ratios as the effect measure. For each dataset, we compared the TGEP and REML point estimates, standard errors, and confidence interval widths. Dataset identifiers and detailed results are provided in Supporting Information (S4 File).

### Software

All analyses were conducted in R version 4.5.2 [18] using the metafor package version 4.8-0 [2] for REML, HKSJ, trim-and-fill, and PET-PEESE estimation. The TGEP implementation is available as an R package at [ZENODO_DOI_PLACEHOLDER] (S2 File). The simulation script is provided in S3 File, and the test suite in S5 File. All datasets required for reproduction are included in the supplementary materials.

## Results

### Simulation Study

Table 2 presents coverage probabilities for all 18 scenarios and five methods. Table 3 presents bias and RMSE for selected scenarios. Full results including SE ratios are provided in S1 Table.

**Table 2. Coverage (%) by Method and Scenario**

| # | Scenario | REML | HKSJ | TGEP | TrimFill | PET-PEESE |
|---|----------|------|------|------|----------|-----------|
| | **Null scenarios** | | | | | |
| 1 | Null_k10_tau0.05 | 92.3 | 95.1 | 93.3 | 89.7 | 92.1 |
| 2 | Null_k10_tau0.20 | 90.4 | 94.6 | 88.6 | 88.3 | 92.3 |
| | **Baseline (no bias)** | | | | | |
| 3 | Base_k5 | 89.3 | 93.1 | 88.0 | 85.6 | 89.9 |
| 4 | Base_k10 | 92.6 | 94.7 | 90.3 | 89.5 | 92.5 |
| 5 | Base_k20 | 92.2 | 93.6 | 88.2 | 86.6 | 92.8 |
| | **Varying heterogeneity** | | | | | |
| 6 | Hetero_tau0.01 | 93.4 | 93.5 | 93.2 | 89.0 | 94.8 |
| 7 | Hetero_tau0.10 | 91.4 | 93.4 | 86.6 | 88.9 | 92.1 |
| 8 | Hetero_tau0.50 | 92.7 | 94.4 | 89.7 | 91.4 | 94.5 |
| 9 | Hetero_tau1.00 | 91.5 | 93.9 | 87.1 | 89.6 | 95.2 |
| | **Publication bias (significance-based)** | | | | | |
| 14 | PubBias_k5_strong | 78.2 | 88.6 | 77.4 | 76.3 | 89.0 |
| 10 | PubBias_k10_mild | 85.0 | 90.5 | 83.6 | 84.5 | 90.8 |
| 11 | PubBias_k10_strong | 66.1 | 78.4 | 65.8 | 67.8 | 88.9 |
| 12 | PubBias_k20_strong | 49.9 | 57.4 | 48.1 | 56.5 | 86.7 |
| 13 | PubBias_k20_severe | 31.3 | 37.2 | 27.3 | 38.1 | 89.4 |
| | **One-sided publication bias** | | | | | |
| 15 | OneSided_k10 | 88.2 | 91.6 | 88.2 | 85.3 | 93.0 |
| 16 | OneSided_k20 | 77.3 | 81.8 | 77.4 | 70.7 | 90.0 |
| | **Combined and special** | | | | | |
| 17 | Combined_hetero_bias | 82.7 | 86.7 | 71.5 | 86.8 | 91.9 |
| 18 | Null_PubBias | 88.5 | 93.4 | 85.9 | 88.3 | 89.5 |

**Table 3. Bias, RMSE, and SE Ratio for Selected Scenarios**

| Scenario | Method | Bias | RMSE | SE Ratio |
|----------|--------|------|------|----------|
| **Null_k10_tau0.05** | REML | 0.002 | 0.093 | 0.974 |
| | HKSJ | 0.002 | 0.093 | 0.969 |
| | TGEP | 0.002 | 0.084 | 1.004 |
| | TrimFill | 0.004 | 0.103 | 0.885 |
| | PET-PEESE | 0.016 | 0.397 | 1.111 |
| **Base_k10** | REML | -0.001 | 0.093 | 0.974 |
| | HKSJ | -0.001 | 0.093 | 0.966 |
| | TGEP | -0.020 | 0.095 | 0.876 |
| | TrimFill | 0.000 | 0.104 | 0.874 |
| | PET-PEESE | -0.056 | 0.395 | 1.090 |
| **Hetero_tau0.50** | REML | -0.005 | 0.187 | 0.997 |
| | HKSJ | -0.005 | 0.187 | 0.997 |
| | TGEP | -0.034 | 0.180 | 0.884 |
| | TrimFill | -0.003 | 0.202 | 0.931 |
| | PET-PEESE | -0.056 | 0.737 | 1.142 |
| **PubBias_k10_strong** | REML | 0.144 | 0.190 | 1.005 |
| | HKSJ | 0.144 | 0.190 | 0.956 |
| | TGEP | 0.134 | 0.184 | 0.887 |
| | TrimFill | 0.131 | 0.190 | 0.903 |
| | PET-PEESE | -0.050 | 0.918 | 0.856 |
| **PubBias_k20_strong** | REML | 0.150 | 0.175 | 0.937 |
| | HKSJ | 0.150 | 0.175 | 0.927 |
| | TGEP | 0.141 | 0.167 | 0.835 |
| | TrimFill | 0.128 | 0.166 | 0.797 |
| | PET-PEESE | 0.008 | 0.424 | 0.989 |
| **PubBias_k20_severe** | REML | 0.208 | 0.227 | 0.961 |
| | HKSJ | 0.208 | 0.227 | 0.908 |
| | TGEP | 0.206 | 0.225 | 0.861 |
| | TrimFill | 0.179 | 0.213 | 0.752 |
| | PET-PEESE | -0.055 | 0.526 | 0.980 |
| **OneSided_k10** | REML | 0.057 | 0.104 | 0.999 |
| | HKSJ | 0.057 | 0.104 | 0.974 |
| | TGEP | 0.038 | 0.096 | 0.912 |
| | TrimFill | 0.058 | 0.115 | 0.879 |
| | PET-PEESE | -0.043 | 0.384 | 1.114 |
| **OneSided_k20** | REML | 0.088 | 0.117 | 1.002 |
| | HKSJ | 0.088 | 0.117 | 0.998 |
| | TGEP | 0.067 | 0.103 | 0.861 |
| | TrimFill | 0.090 | 0.133 | 0.798 |
| | PET-PEESE | -0.014 | 0.315 | 1.051 |
| **Combined_hetero_bias** | REML | 0.132 | 0.230 | 0.957 |
| | HKSJ | 0.132 | 0.230 | 0.958 |
| | TGEP | 0.126 | 0.226 | 0.753 |
| | TrimFill | 0.091 | 0.211 | 0.955 |
| | PET-PEESE | -0.051 | 0.788 | 1.102 |

#### Null scenarios and type I error

Under null scenarios without publication bias (Scenarios 1-2), HKSJ achieved coverage closest to the nominal 95% level (95.1% and 94.6%). REML achieved 92.3% and 90.4% coverage. TGEP achieved 93.3% and 88.6% coverage, respectively. At low heterogeneity (tau^2 = 0.05), TGEP had well-calibrated SE (ratio 1.00) and coverage (93.3%). At higher heterogeneity (tau^2 = 0.20), SE calibration remained adequate (0.93) but coverage declined to 88.6%, reflecting the wider empirical SE distribution. PET-PEESE maintained coverage near nominal (92.1%) but with substantially higher RMSE (0.397 vs. 0.093 for REML), reflecting the cost of estimating unnecessary bias-correction parameters.

#### Effect of number of studies

Across baseline scenarios (Scenarios 3-5, mu = 0.3, no bias), HKSJ maintained the best coverage (93.1-94.7%), consistent with its t-distribution correction for small-sample inference. REML coverage improved with k (89.3% at k = 5 to 92.6% at k = 10). TGEP coverage was comparable to REML at k = 5 (88.0% vs. 89.3%) but fell below at k = 10 (90.3% vs. 92.6%) and k = 20 (88.2% vs. 92.2%), reflecting the SE underestimation issue (SE ratio 0.84-0.89). RMSE was similar across REML, HKSJ, and TGEP (within 3%), with TGEP introducing a slight negative bias (-0.016 to -0.023) from the SWA guard's up-weighting of non-significant studies even in the absence of publication bias.

#### Effect of heterogeneity

Under varying heterogeneity (Scenarios 6-9), HKSJ maintained excellent coverage across the full range (93.4-94.4%). REML coverage was adequate (91.4-93.4%). TGEP coverage varied across the heterogeneity range: from 93.2% at tau^2 = 0.01 to a minimum of 86.6% at tau^2 = 0.10 and 87.1% at tau^2 = 1.00, with SE ratios following a non-monotonic pattern (0.94, 0.81, 0.88, 0.86 at tau^2 = 0.01, 0.10, 0.50, 1.00 respectively). This indicates that the perfect-correlation variance formula, while an improvement over the independent-guard formula, increasingly underestimates the true ensemble variability as heterogeneity grows and guard estimates diverge more widely.

TGEP achieved the lowest RMSE at high heterogeneity (0.180 at tau^2 = 0.50 vs. REML 0.187; 0.252 at tau^2 = 1.00 vs. REML 0.269), indicating that the ensemble weighting effectively combines guard estimates in high-heterogeneity settings even though the SE is underestimated.

#### Publication bias: significance-based

Under significance-based publication bias (Scenarios 10-14), all standard methods showed substantial coverage deterioration with increasing bias severity. At strong bias with k = 20 (Scenario 12), REML coverage dropped to 49.9% and HKSJ to 57.4%, with TGEP at 48.1%. Under severe bias (Scenario 13), REML coverage was 31.3%, HKSJ 37.2%, and TGEP 27.3%.

PET-PEESE was clearly the best method for maintaining coverage under significance-based publication bias, achieving 86.7-90.8% across mild-to-severe bias scenarios. Trim-and-fill also outperformed REML and TGEP at higher bias levels (38.1% vs. 31.3% and 27.3% under severe bias, respectively).

TGEP's coverage disadvantage relative to REML under significance-based bias reflects two compounding factors: (1) the LOO-CV stacking optimizes prediction of the observed (biased) data, which may assign insufficient weight to the SWA bias-correction guard; and (2) the analytical SE underestimates the true variability. However, TGEP did achieve modestly lower bias than REML (0.134 vs. 0.144 at k = 10 strong; 0.141 vs. 0.150 at k = 20 strong; 0.206 vs. 0.208 at k = 20 severe) and lower RMSE (0.184 vs. 0.190; 0.167 vs. 0.175; 0.225 vs. 0.227).

#### Publication bias: one-sided

Under one-sided publication bias (Scenarios 15-16), TGEP showed its clearest point-estimation advantage. At k = 10, TGEP bias was 0.038 vs. REML 0.057 (33% reduction), with RMSE 0.096 vs. 0.104 (7.7% improvement). At k = 20, TGEP bias was 0.067 vs. REML 0.088 (24% reduction), with RMSE 0.103 vs. 0.117 (12.4% improvement). These represent the largest relative improvements observed for TGEP across all scenarios.

Coverage was similar between TGEP and REML (88.2% vs. 88.2% at k = 10; 77.4% vs. 77.3% at k = 20), with HKSJ again providing the best coverage (91.6% and 81.8%). Trim-and-fill performed poorly under one-sided bias (85.3% and 70.7%), as the funnel-plot symmetry assumption is particularly violated under directional selection.

#### Combined heterogeneity and publication bias

Under combined high heterogeneity and significance-based bias (Scenario 17), TGEP's coverage was notably low (71.5%) compared to REML (82.7%) and HKSJ (86.7%). The SE ratio was 0.75, the lowest observed across all scenarios, indicating that the combination of high heterogeneity and publication bias creates the most challenging conditions for TGEP's analytical SE. Nevertheless, TGEP achieved the lowest RMSE among the three core methods (0.226 vs. REML 0.230).

#### SE calibration

Across all scenarios, TGEP's SE ratio ranged from approximately 0.75 to 1.00. SE was well-calibrated (ratio 0.93-1.00) in null scenarios with moderate heterogeneity and in the low-heterogeneity condition. SE was progressively underestimated (ratio 0.75-0.89) as heterogeneity or publication bias severity increased. By comparison, REML SE ratios ranged from approximately 0.90-1.01, and HKSJ SE ratios from 0.81-1.00, both showing better calibration on average. This SE underestimation is the primary driver of TGEP's coverage shortfall.

### Empirical validation

Table 4 presents the results of applying TGEP and REML to 15 Cochrane datasets.

**Table 4. Empirical Validation: TGEP vs. REML on 15 Cochrane Datasets**

| Dataset | K | REML Est. | TGEP Est. | Abs. Diff | REML SE | TGEP SE | Relative SE |
|---------|---|-----------|-----------|-----------|---------|---------|-------------|
| CD015252 | 8 | 0.332 | 0.303 | 0.029 | 0.215 | 0.486 | 2.26 |
| CD012445 | 8 | -0.083 | -0.081 | 0.002 | 0.716 | 0.701 | 0.98 |
| CD013844 | 12 | -0.495 | -0.264 | 0.231 | 0.474 | 0.748 | 1.58 |
| CD015532 | 22 | -0.074 | -0.044 | 0.030 | 0.119 | 0.377 | 3.18 |
| CD013366 | 13 | -0.590 | -0.517 | 0.073 | 0.266 | 0.542 | 2.04 |
| CD014678 | 40 | -0.009 | -0.005 | 0.004 | 0.065 | 0.244 | 3.74 |
| CD006140 | 9 | -0.659 | -0.599 | 0.060 | 0.251 | 0.394 | 1.57 |
| CD008493 | 20 | 0.244 | 0.308 | 0.064 | 0.167 | 0.565 | 3.38 |
| CD013055 | 6 | 1.109 | 1.151 | 0.042 | 0.158 | 0.303 | 1.92 |
| CD015087 | 13 | -0.562 | -0.525 | 0.036 | 0.317 | 0.806 | 2.55 |
| CD014960 | 5 | 0.444 | 0.596 | 0.152 | 0.340 | 0.539 | 1.59 |
| CD009362 | 56 | -0.477 | -0.445 | 0.032 | 0.138 | 0.470 | 3.41 |
| CD013474 | 7 | -0.142 | -0.117 | 0.025 | 0.761 | 0.834 | 1.10 |
| CD015140 | 7 | 0.017 | 0.017 | 0.000 | 0.761 | 0.664 | 0.87 |
| CD015518 | 88 | -0.058 | -0.053 | 0.005 | 0.215 | 0.721 | 3.35 |
| **Mean** | | | | **0.052** | | | **2.23** |
| **Median** | | | | **0.032** | | | **2.04** |

Across the 15 datasets, TGEP estimates differed from REML by a mean absolute difference of 0.052 (median 0.032) on the standardized scale. The mean relative SE (TGEP SE / REML SE) was 2.23, substantially wider than the simulation SE ratio within TGEP alone (0.75-1.00), because the empirical SE ratio here compares TGEP against REML rather than TGEP's model SE against its own empirical SE. The wider TGEP SEs in the empirical validation largely reflect the conservative perfect-correlation variance formula applied across three guards whose individual variances may differ substantially from the REML variance.

The largest point-estimate difference occurred for CD013844 (K = 12), where TGEP pulled the estimate toward zero (from -0.495 to -0.264), suggesting the GRMA or SWA guards identified influential studies or potential publication bias. In datasets with minimal heterogeneity (e.g., CD015140), TGEP and REML produced nearly identical results.

## Discussion

### Summary of findings

We proposed TGEP, a frequentist ensemble framework for meta-analysis that integrates three specialized guard estimators via LOO-CV stacking. In a comprehensive simulation study comparing five methods across 18 scenarios, TGEP demonstrated two principal strengths and one important limitation.

**Strengths.** First, TGEP produced less biased point estimates than REML under one-sided publication bias, reducing absolute bias by 24-33% and RMSE by approximately 8-12%. This advantage was most pronounced when the selection mechanism was directional rather than significance-based. Second, TGEP achieved lower RMSE than REML in high-heterogeneity scenarios (tau^2 = 0.50-1.00) and in null-effect settings, suggesting that the ensemble weighting can improve estimation precision even when publication bias is not the primary concern.

**Limitation.** TGEP's confidence interval coverage was below nominal in most scenarios (approximately 87-93% in unbiased settings, 27-88% under publication bias), consistently worse than both REML and HKSJ. The primary cause was SE underestimation: TGEP's analytical SE ratio ranged from approximately 0.75 to 1.00 across scenarios, with the most severe underestimation occurring under combined heterogeneity and bias. This limits TGEP's utility for formal statistical inference.

### Comparison with existing methods

Our simulation provides a comparative assessment of five commonly used approaches:

**HKSJ** provided the most reliable coverage in unbiased and mildly biased scenarios (93-95%), and the best relative coverage under publication bias among non-bias-correction methods, though coverage still declined substantially under strong significance-based bias (37-91%) [6]. Its t-distribution correction addresses SE underestimation effectively but does not correct the biased point estimate.

**PET-PEESE** was the most effective method for maintaining coverage under significance-based publication bias (approximately 87-91%), consistent with prior comparisons [5]. However, PET-PEESE exhibited substantially inflated RMSE in unbiased scenarios (approximately 4-5x higher than REML), reflecting the cost of fitting unnecessary bias-correction parameters. This overcorrection is a well-known limitation [9].

**Trim-and-fill** provided moderate coverage improvements under significance-based bias (38-57% vs. REML 31-50% under strong/severe bias) but performed poorly under one-sided bias (71-85%), where the funnel-plot symmetry assumption is most violated.

**TGEP** occupies a complementary niche: it provides modest point-estimation improvements without the overcorrection risk of PET-PEESE. In baseline scenarios (non-null effect, no bias), TGEP RMSE was within 0.5-3% of REML; in null-effect and high-heterogeneity settings, TGEP RMSE was 4-10% lower than REML, reflecting favorable bias-variance tradeoffs from the ensemble weighting. By contrast, PET-PEESE RMSE was approximately 4-5x higher than REML in unbiased settings. Under one-sided bias, TGEP achieved the largest bias reductions among all methods that do not assume a specific bias mechanism.

### The LOO-CV limitation

A fundamental tension in TGEP's design is that LOO-CV optimizes prediction of the observed data, which may itself be biased by selective publication. When the selection mechanism is significance-based, the SWA guard — which up-weights non-significant studies — may actually predict the observed (predominantly significant) data less well than the WRD or GRMA guards that do not attempt bias correction. This can result in the LOO-CV assigning suboptimal weight to the guard most relevant for bias correction.

This limitation is inherent to any stacking approach applied to biased data and explains why TGEP's bias correction is modest compared to methods like PET-PEESE that directly model the bias mechanism. Under one-sided bias, where the selection is directional rather than significance-based, this tension is less severe, and TGEP's advantage is correspondingly larger.

### GRMA guard under publication bias

A subtlety of the GRMA guard is that it uses the median of observed effect sizes as its reference point. Under publication bias, this median is inflated, causing the GRMA guard to assign higher weights to studies close to the biased median rather than to the true population mean. In this sense, GRMA functions primarily as an outlier-robustness tool rather than a bias-correction tool, and under publication bias it may partially counteract the SWA guard's correction. This interaction between guards is mediated by the LOO-CV stacking weights but is not explicitly modeled.

### SWA selection probability mismatch

The SWA guard assumes a fixed selection probability of 0.4 for non-significant studies, corresponding to an up-weighting factor of 2.5. In the simulation, the actual selection probabilities were 0.50 (mild), 0.20 (strong), and 0.05 (severe). Under severe bias, the true selection probability (0.05) is 8-fold lower than the SWA assumption (0.40), leading to substantial under-correction. This mismatch is the primary reason TGEP's bias reduction diminishes from 23% under mild bias to less than 1% under severe bias. A data-adaptive approach that estimates the selection probability from the observed p-value distribution could substantially improve TGEP's bias correction under strong selection.

### SE underestimation

TGEP's analytical SE assumes perfect positive correlation among the three guards. Despite being the most conservative correlation assumption (yielding the widest CI for a given set of guard variances), the resulting SE still underestimates the empirical variability. This occurs because the guard-specific variances V_g themselves underestimate the total uncertainty of the ensemble: they do not account for the additional variability introduced by the data-dependent LOO-CV weight selection. When the stacking weights alpha_g vary across bootstrap samples, the ensemble estimate acquires additional variability that the analytical formula does not capture.

Bootstrap SE (B = 200, available as a default option) provides a more accurate variance estimate that accounts for weight uncertainty. However, bootstrap SE was not used in the simulation study for computational reasons (1000 replications x 18 scenarios x 200 bootstrap samples = 3.6 million TGEP evaluations). Future work should evaluate whether bootstrap SE resolves the coverage shortfall.

### Strengths of the framework

1. **Assumption-lean**: No distributional priors, selection functions, or parametric bias models are required.
2. **Data-adaptive**: LOO-CV stacking weights automatically select the most appropriate guard combination for each dataset.
3. **Interpretable**: Individual guard estimates and weights are reported, allowing researchers to assess which type of contamination the ensemble is responding to.
4. **Stable under no bias**: Unlike PET-PEESE, TGEP does not overcorrect in unbiased settings.
5. **Extensible**: The framework can incorporate additional guards targeting specific bias mechanisms (e.g., p-hacking, small-study effects) without changing the stacking infrastructure.
6. **Computationally efficient**: TGEP with analytical SE runs in seconds even for large meta-analyses.

### Limitations

1. **SE underestimation**: The analytical SE formula underestimates the true ensemble variability, producing coverage below the nominal 95% level. Bootstrap SE may resolve this but was not evaluated in the simulation.

2. **LOO-CV and biased data**: LOO-CV optimizes prediction of observed data, not the unobserved population parameter. Under strong significance-based publication bias, this can assign insufficient weight to the bias-correcting SWA guard.

3. **SWA calibration**: The SWA guard uses a fixed up-weighting factor of 2.5 for non-significant studies. A data-adaptive approach estimating the selection probability [10] may improve performance. Sensitivity to this parameter was not systematically evaluated.

4. **Normal approximation**: TGEP's confidence intervals use the normal distribution. For small K, a t-distribution reference (as in HKSJ) may improve coverage. Combining TGEP point estimates with HKSJ-style SE correction is a promising direction for future work.

5. **Guard correlation**: WRD and SWA share an initial REML fit, creating structural dependence beyond the shared-data correlation already assumed. The impact of this specific dependence on the variance formula was not analytically characterized.

6. **Simulation scope**: We evaluated N_sim = 1000 replications per scenario, which provides MCSE of coverage of approximately 0.007-0.016. While adequate for identifying the magnitude of coverage differences, larger simulations would further refine the estimates.

7. **Number of guards**: The three-guard architecture is fixed. Whether the optimal number of guards differs across settings is unknown.

8. **Cochrane dataset selection**: The 15 empirical validation datasets were selected for availability and diversity (Cochrane reviews published 2020-2025, with K >= 5 studies and log odds ratio as effect measure), not by a systematic protocol. Future work should evaluate TGEP on a pre-registered set of datasets with known bias characteristics.

9. **Effect measure scope**: The simulation study evaluated only standardized mean differences (SMD). The empirical validation used log odds ratios, which have different distributional properties. The generalizability of the simulation findings to the log OR scale was not formally established.

10. **Divergent TGEP estimates**: In 3 of 15 empirical datasets (CD008493, CD013055, CD014960), TGEP moved the point estimate away from zero relative to REML, contrary to the expected bias-reduction direction. This can occur when the GRMA or WRD guards identify influential studies on the opposite side of the estimate distribution. Users should interpret TGEP as a data-adaptive reweighting, not a guaranteed bias reduction.

### Recommendations for practice

Based on our findings, we recommend TGEP as a **complementary sensitivity analysis** rather than a replacement for established methods. Specifically:

- Use **HKSJ** as the primary method for inference (confidence intervals and hypothesis tests), given its consistently near-nominal coverage.
- Use **PET-PEESE** when funnel plot asymmetry or other evidence suggests significance-based publication bias, noting its overcorrection risk in unbiased settings.
- Use **TGEP** to provide a bias-reduced point estimate and to examine guard weight diagnostics as informal indicators of potential bias or outlier contamination:
  - When SWA weight exceeds 0.5, this suggests the data may be affected by publication bias (more non-significant studies are being up-weighted).
  - When GRMA weight exceeds 0.5, outlier or small-study effects may be present.
  - When |TGEP - REML| exceeds 0.10 on the SMD scale (or 10% of the REML estimate for other scales), publication bias or outlier contamination should be investigated further using dedicated methods (funnel plots, Egger's test, leave-one-out diagnostics).
- Report individual guard estimates and weights alongside the ensemble result, e.g., in a supplementary table or sensitivity analysis section.
- TGEP results should appear in the sensitivity analysis section of a systematic review report, not as the primary estimate in the main forest plot.

The default temperature T = 1.0 is suitable for most applications. Bootstrap SE (n_boot >= 200) should be used when confidence intervals from TGEP are needed.

## Conclusions

TGEP provides a frequentist ensemble framework for meta-analysis that produces less biased point estimates under publication bias, particularly one-sided selection, without the overcorrection risk of parametric bias-correction methods. However, its analytical standard error underestimates the true variability, limiting its utility for formal inference. We recommend TGEP as a complementary sensitivity analysis alongside HKSJ for coverage and PET-PEESE for bias assessment. The method is implemented as an open-source R package and requires no specification of priors, selection functions, or parametric bias models. Future work should address the SE calibration through bootstrap validation or hybrid approaches combining TGEP point estimates with HKSJ-style variance correction.

## References

[1] Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. Chichester: John Wiley & Sons; 2009. https://doi.org/10.1002/9780470743386

[2] Viechtbauer W. Conducting meta-analyses in R with the metafor package. J Stat Softw. 2010;36(3):1-48. https://doi.org/10.18637/jss.v036.i03

[3] Rothstein HR, Sutton AJ, Borenstein M, editors. Publication Bias in Meta-Analysis: Prevention, Assessment and Adjustments. Chichester: John Wiley & Sons; 2005. https://doi.org/10.1002/0470870168

[4] Viechtbauer W, Cheung MW-L. Outlier and influence diagnostics for meta-analysis. Res Synth Methods. 2010;1(2):112-125. https://doi.org/10.1002/jrsm.11

[5] Carter EC, Schonbrodt FD, Gervais WM, Hilgard J. Correcting for bias in psychology: A comparison of meta-analytic methods. Adv Methods Pract Psychol Sci. 2019;2(2):115-144. https://doi.org/10.1177/2515245919847196

[6] Hartung J, Knapp G. A refined method for the meta-analysis of controlled clinical trials with binary outcome. Stat Med. 2001;20(24):3875-3889. https://doi.org/10.1002/sim.1009

[7] Duval S, Tweedie R. Trim and fill: A simple funnel-plot-based method of testing and adjusting for publication bias in meta-analysis. Biometrics. 2000;56(2):455-463. https://doi.org/10.1111/j.0006-341X.2000.00455.x

[8] Stanley TD, Doucouliagos H. Meta-regression approximations to reduce publication selection bias. Res Synth Methods. 2014;5(1):60-78. https://doi.org/10.1002/jrsm.1095

[9] Stanley TD. Limitations of PET-PEESE and other meta-analysis methods. Soc Psychol Personal Sci. 2017;8(5):581-591. https://doi.org/10.1177/1948550617693062

[10] Vevea JL, Hedges LV. A general linear model for estimating effect size in the presence of publication bias. Psychometrika. 1995;60(3):419-435. https://doi.org/10.1007/BF02294384

[11] Maier M, Bartos F, Wagenmakers E-J. Robust Bayesian meta-analysis: Addressing publication bias with model-averaging. Psychol Methods. 2023;28(1):107-122. https://doi.org/10.1037/met0000405

[12] Wolpert DH. Stacked generalization. Neural Netw. 1992;5(2):241-259. https://doi.org/10.1016/S0893-6080(05)80023-1

[13] Deng J-L. Introduction to grey system theory. J Grey Syst. 1989;1(1):1-24.

[14] Hastings C, Mosteller F, Tukey JW, Winsor CP. Low moments for small samples: A comparative study of order statistics. Ann Math Stat. 1947;18(3):413-426. https://doi.org/10.1214/aoms/1177730388

[15] Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate statistical methods. Stat Med. 2019;38(11):2074-2102. https://doi.org/10.1002/sim.8086

[16] Sidik K, Jonkman JN. A simple confidence interval for meta-analysis. Stat Med. 2002;21(21):3153-3159. https://doi.org/10.1002/sim.1262

[17] IntHout J, Ioannidis JPA, Borm GF. The Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis is straightforward and considerably outperforms the standard DerSimonian-Laird method. BMC Med Res Methodol. 2014;14:25. https://doi.org/10.1186/1471-2288-14-25

[18] R Core Team. R: A language and environment for statistical computing. Vienna: R Foundation for Statistical Computing; 2025. https://www.R-project.org/

---

## Supporting information

**S1 Table.** Complete simulation results for all 18 scenarios and 5 methods, including bias, RMSE, coverage, MCSE, SE ratio, and mean number of studies after selection.

**S2 File.** R package source code for TGEP (TGEP.R), including all guard implementations, the LOO-CV stacking ensemble, and diagnostic plotting functions.

**S3 File.** Simulation script (run_simulation.R) with complete ADEMP specification, data-generating mechanisms, and performance measure calculations for all five methods.

**S4 File.** Cochrane dataset identifiers and empirical validation results (Real_World_Impact_Summary.txt).

**S5 File.** Test suite (test_tgep.R) with 30 assertions across 17 test groups covering basic functionality, edge cases, input validation, and numerical stability.

---

## Author information

### Affiliations
^1 Independent Researcher

### Corresponding author
Mahmood Ul Hassan (mahmood.hassan@example.com)

### Author contributions (CRediT)
**Conceptualization:** Mahmood Ul Hassan. **Methodology:** Mahmood Ul Hassan. **Software:** Mahmood Ul Hassan. **Validation:** Mahmood Ul Hassan. **Formal analysis:** Mahmood Ul Hassan. **Investigation:** Mahmood Ul Hassan. **Data curation:** Mahmood Ul Hassan. **Writing - original draft:** Mahmood Ul Hassan. **Writing - review & editing:** Mahmood Ul Hassan. **Visualization:** Mahmood Ul Hassan.

### Data availability statement
All simulation code and results are provided as Supporting Information. The Cochrane datasets used for empirical validation are identified by their Cochrane review identifiers (CD-numbers) and are publicly available through the Cochrane Library (https://www.cochranelibrary.com/). The TGEP R package is available at [ZENODO_DOI_PLACEHOLDER].

### Acknowledgments
The author thanks the developers of the metafor R package for providing the statistical infrastructure used in this study.

### Funding
The author received no specific funding for this work.

### Competing interests
The author has declared that no competing interests exist.

### Ethics statement
This study used simulated data and publicly available summary statistics from published Cochrane reviews. No individual patient data were used. No ethics approval was required.
