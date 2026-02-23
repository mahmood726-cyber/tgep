# Triple-Guard Ensemble Pooling: A Frequentist Ensemble Framework for Bias-Reduced Point Estimation in Meta-Analysis

**Mahmood Ul Hassan**^1

^1 Independent Researcher

**Corresponding author:** Mahmood Ul Hassan (mahmood.hassan@example.com)

---

## Abstract

**Background:** Standard inverse-variance random-effects meta-analysis (REML) produces biased point estimates when studies are contaminated by publication bias, yet existing bias-correction methods such as PET-PEESE can overcorrect when bias is absent. A method that reduces bias without strong parametric assumptions about the selection mechanism would complement the existing toolkit.

**Methods:** We propose Triple-Guard Ensemble Pooling (TGEP), a frequentist ensemble framework that integrates three specialized estimators — Grey Relational Meta-Analysis (GRMA), Winsorized Robust Detection (WRD), and Significance-Weighted Adjustment (SWA) — via leave-one-out cross-validation (LOO-CV) stacking with softmax temperature weighting. A Monte Carlo simulation study (500 iterations x 18 scenarios) evaluated TGEP against REML, Hartung-Knapp-Sidik-Jonkman (HKSJ), trim-and-fill, and PET-PEESE across varying numbers of studies (k = 5-20), between-study heterogeneity (tau-squared = 0.01-1.00), and publication bias mechanisms (significance-based and one-sided). The target estimand was the unconditional population mean effect mu (pre-selection).

**Results:** Under one-sided publication bias (k = 20, tau-squared = 0.10), TGEP reduced absolute bias by 24% relative to REML (0.069 vs. 0.091) and achieved the lowest RMSE (0.107 vs. REML 0.121). Under significance-based bias, TGEP showed modest bias reductions (6-11% vs. REML). In unbiased scenarios, TGEP maintained RMSE comparable to REML (within 5%). However, TGEP's analytical standard error underestimated empirical variability (SE ratio 0.83-1.06), producing coverage below the nominal 95% level (86-95% in unbiased settings). PET-PEESE achieved the best coverage under publication bias (88-90%) but overcorrected in unbiased settings (RMSE 2-4x larger). HKSJ provided the most reliable coverage across all scenarios (93-96%).

**Conclusions:** TGEP produces less biased point estimates than REML under publication bias, particularly one-sided selection, without the overcorrection risk of PET-PEESE. However, its confidence interval coverage is limited by analytical SE underestimation. We recommend TGEP as a complementary sensitivity analysis alongside HKSJ for inference and PET-PEESE for bias assessment. The method is implemented as an open-source R package.

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

The GRMA pooled estimate is mu_hat_GRMA = sum_i(w_i^GRMA * y_i) with variance V_GRMA = sum_i((w_i^GRMA)^2 * v_i).

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

This yields V_TGEP = SE_TGEP^2. This is conservative relative to the independent-guard formula sum_g(alpha_g^2 * V_g) but, as our simulation study demonstrates, still underestimates the true variability of the ensemble estimator. Confidence intervals use the normal approximation: mu_hat_TGEP +/- z_(1-alpha/2) * SE_TGEP. A nonparametric bootstrap (B = 200 by default) is also available for variance estimation that accounts for the additional uncertainty from the LOO-CV weight selection.

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
2. Compute the within-study variance using the SMD variance approximation: v_i = 2/n_i + mu^2/(4*n_i), yielding realistic variances in the range 0.02-0.10 for the chosen effect sizes.
3. Sample the true study-specific effect: theta_i ~ N(mu, tau^2)
4. Sample the observed effect: y_i ~ N(theta_i, v_i)

Publication bias was simulated via two mechanisms:

- **Significance-based**: Studies with p < 0.05 (two-sided z-test using within-study SE) were always published; non-significant studies were published with probability 1 - beta, where beta in {0.50, 0.80, 0.95} controlled bias severity.
- **One-sided**: Studies with effects in the expected direction (y_i > 0) were always published; negative effects were published with probability 1 - beta.

A minimum of 3 studies was enforced in all scenarios to ensure estimability.

#### Methods Compared

1. **REML**: Standard restricted maximum likelihood random-effects meta-analysis [2] with Wald-type confidence intervals.
2. **HKSJ**: REML with the Hartung-Knapp-Sidik-Jonkman adjustment [6] using a t_(K-1) reference distribution for confidence intervals.
3. **TGEP**: Triple-Guard Ensemble Pooling with T = 1.0 and analytical variance (no bootstrap, for computational efficiency in the simulation).
4. **Trim-and-fill**: The Duval-Tweedie trim-and-fill method [7] applied to the REML fit, which estimates the number of "missing" studies and imputes them to restore funnel plot symmetry.
5. **PET-PEESE**: The precision-effect test / precision-effect estimate with standard error (PET-PEESE) [8]. PET regresses effect sizes on standard errors; if the PET intercept is significant at p < 0.10, PEESE (regression on variances) is used instead.

If any of the three core methods (REML, HKSJ, TGEP) failed for a given iteration, that entire iteration was discarded. Trim-and-fill and PET-PEESE failures were treated as missing for that iteration only.

#### Performance Measures

- **Bias**: Mean difference between estimated and true mu
- **RMSE**: Root mean squared error, sqrt(mean((mu_hat - mu)^2))
- **Coverage**: Proportion of 95% confidence intervals containing the true mu
- **MCSE of coverage**: Monte Carlo standard error, sqrt(C*(1-C)/N_sim)
- **SE ratio**: Mean model-based SE divided by empirical SE; values near 1.0 indicate well-calibrated SE estimation; values below 1.0 indicate anti-conservative SE (confidence intervals too narrow)

#### Scenarios

Table 1 summarizes the 18 simulation scenarios. Each scenario was replicated N_sim = 500 times with a fixed seed (20260223) for reproducibility.

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

### Empirical Validation

We applied TGEP and REML to 15 datasets from Cochrane systematic reviews spanning diverse clinical areas. These datasets ranged from K = 5 to K = 88 studies and used log odds ratios as the effect measure. For each dataset, we compared the TGEP and REML point estimates, standard errors, and confidence interval widths. Dataset identifiers and detailed results are provided in Supporting Information (S4 File).

### Software

All analyses were conducted in R version 4.5.2 [16] using the metafor package version 4.8-0 [2] for REML, HKSJ, trim-and-fill, and PET-PEESE estimation. The TGEP implementation is available as an R package at [ZENODO_DOI_PLACEHOLDER]. Simulation code and all datasets required for reproduction are included in the supplementary materials.

## Results

### Simulation Study

Table 2 presents coverage probabilities for all 18 scenarios and five methods. Table 3 presents bias and RMSE for selected scenarios. Full results including SE ratios are provided in S1 Table.

**Table 2. Coverage (%) by Method and Scenario**

| Scenario | REML | HKSJ | TGEP | TrimFill | PET-PEESE |
|----------|------|------|------|----------|-----------|
| **Null scenarios** | | | | | |
| Null_k10_tau0.05 | 93.6 | 96.0 | 94.8 | 90.8 | 93.6 |
| Null_k10_tau0.20 | 90.8 | 93.8 | 89.6 | 88.8 | 90.8 |
| **Baseline (no bias)** | | | | | |
| Base_k5 | 89.4 | 95.2 | 89.4 | 86.6 | 89.6 |
| Base_k10 | 90.4 | 94.8 | 87.8 | 87.6 | 90.2 |
| Base_k20 | 94.4 | 95.2 | 90.8 | 89.6 | 93.8 |
| **Varying heterogeneity** | | | | | |
| Hetero_tau0.01 | 94.6 | 95.8 | 94.0 | 90.6 | 93.8 |
| Hetero_tau0.10 | 94.2 | 95.6 | 90.8 | 89.2 | 93.6 |
| Hetero_tau0.50 | 93.0 | 95.6 | 87.4 | 91.8 | 92.8 |
| Hetero_tau1.00 | 92.4 | 95.6 | 86.6 | 91.4 | 94.6 |
| **Publication bias (significance-based)** | | | | | |
| PubBias_k5_strong | 81.4 | 86.0 | 80.2 | 81.2 | 88.4 |
| PubBias_k10_mild | 83.2 | 89.6 | 83.6 | 83.2 | 88.6 |
| PubBias_k10_strong | 63.4 | 74.2 | 62.0 | 64.0 | 90.0 |
| PubBias_k20_strong | 50.8 | 57.4 | 48.0 | 57.8 | 88.8 |
| PubBias_k20_severe | 32.2 | 39.2 | 27.6 | 42.2 | 89.8 |
| **One-sided publication bias** | | | | | |
| OneSided_k10 | 87.4 | 91.4 | 87.8 | 84.8 | 93.6 |
| OneSided_k20 | 77.8 | 82.2 | 79.4 | 70.0 | 89.0 |
| **Combined and special** | | | | | |
| Combined_hetero_bias | 82.2 | 87.4 | 71.4 | 86.6 | 92.6 |
| Null_PubBias | 88.0 | 93.8 | 86.8 | 88.6 | 89.0 |

**Table 3. Bias, RMSE, and SE Ratio for Selected Scenarios**

| Scenario | Method | Bias | RMSE | SE Ratio |
|----------|--------|------|------|----------|
| **Null_k10_tau0.05** | REML | 0.001 | 0.088 | 1.029 |
| | HKSJ | 0.001 | 0.088 | 1.024 |
| | TGEP | 0.002 | 0.080 | 1.057 |
| | TrimFill | 0.000 | 0.097 | 0.948 |
| | PET-PEESE | 0.017 | 0.392 | 1.117 |
| **Base_k10** | REML | -0.002 | 0.093 | 0.979 |
| | HKSJ | -0.002 | 0.093 | 0.973 |
| | TGEP | -0.021 | 0.096 | 0.875 |
| | TrimFill | 0.003 | 0.104 | 0.884 |
| | PET-PEESE | -0.018 | 0.393 | 1.069 |
| **Hetero_tau0.50** | REML | -0.007 | 0.195 | 0.970 |
| | HKSJ | -0.007 | 0.195 | 0.970 |
| | TGEP | -0.036 | 0.188 | 0.855 |
| | TrimFill | -0.007 | 0.209 | 0.911 |
| | PET-PEESE | -0.073 | 0.741 | 1.127 |
| **PubBias_k10_strong** | REML | 0.141 | 0.194 | 0.918 |
| | HKSJ | 0.141 | 0.194 | 0.864 |
| | TGEP | 0.131 | 0.188 | 0.817 |
| | TrimFill | 0.131 | 0.196 | 0.828 |
| | PET-PEESE | -0.028 | 0.846 | 0.961 |
| **PubBias_k20_strong** | REML | 0.150 | 0.175 | 0.964 |
| | HKSJ | 0.150 | 0.175 | 0.951 |
| | TGEP | 0.141 | 0.168 | 0.851 |
| | TrimFill | 0.125 | 0.164 | 0.803 |
| | PET-PEESE | -0.002 | 0.437 | 0.995 |
| **PubBias_k20_severe** | REML | 0.203 | 0.222 | 0.994 |
| | HKSJ | 0.203 | 0.222 | 0.941 |
| | TGEP | 0.201 | 0.220 | 0.889 |
| | TrimFill | 0.168 | 0.202 | 0.785 |
| | PET-PEESE | -0.060 | 0.448 | 1.064 |
| **OneSided_k10** | REML | 0.063 | 0.106 | 1.012 |
| | HKSJ | 0.063 | 0.106 | 0.989 |
| | TGEP | 0.044 | 0.097 | 0.927 |
| | TrimFill | 0.060 | 0.115 | 0.880 |
| | PET-PEESE | -0.047 | 0.384 | 1.117 |
| **OneSided_k20** | REML | 0.091 | 0.121 | 0.965 |
| | HKSJ | 0.091 | 0.121 | 0.961 |
| | TGEP | 0.069 | 0.107 | 0.825 |
| | TrimFill | 0.092 | 0.139 | 0.757 |
| | PET-PEESE | 0.015 | 0.306 | 1.088 |
| **Combined_hetero_bias** | REML | 0.134 | 0.227 | 0.971 |
| | HKSJ | 0.134 | 0.227 | 0.971 |
| | TGEP | 0.126 | 0.220 | 0.768 |
| | TrimFill | 0.096 | 0.207 | 0.973 |
| | PET-PEESE | -0.042 | 0.800 | 1.049 |

#### Null Scenarios and Type I Error

Under null scenarios without publication bias (Scenarios 1-2), all methods maintained coverage near or above the nominal 95% level for HKSJ (96.0% and 93.8%) and near nominal for REML (93.6% and 90.8%). TGEP achieved 94.8% and 89.6% coverage, respectively. At low heterogeneity (tau^2 = 0.05), TGEP had well-calibrated SE (ratio 1.06) and coverage (94.8%). At higher heterogeneity (tau^2 = 0.20), SE calibration remained adequate (0.94) but coverage declined to 89.6%, reflecting the wider empirical SE distribution. PET-PEESE maintained nominal coverage (93.6%) but with substantially higher RMSE (0.392 vs. 0.088 for REML), reflecting the cost of estimating unnecessary bias-correction parameters.

#### Effect of Number of Studies

Across baseline scenarios (Scenarios 3-5, mu = 0.3, no bias), HKSJ maintained the best coverage (94.8-95.2%), consistent with its t-distribution correction for small-sample inference. REML coverage improved with k (89.4% at k = 5 to 94.4% at k = 20). TGEP coverage was comparable to REML at k = 5 (89.4%) but fell below at k = 10 (87.8%) and k = 20 (90.8%), reflecting the SE underestimation issue (SE ratio 0.87-0.89). RMSE was similar across REML, HKSJ, and TGEP (within 3%), with TGEP introducing a slight negative bias (-0.015 to -0.024) from the SWA guard's up-weighting of non-significant studies even in the absence of publication bias.

#### Effect of Heterogeneity

Under varying heterogeneity (Scenarios 6-9), HKSJ maintained excellent coverage across the full range (95.6-95.8%). REML coverage was adequate (92.4-94.6%). TGEP coverage declined progressively with increasing tau^2: from 94.0% at tau^2 = 0.01 to 86.6% at tau^2 = 1.00, with SE ratios declining from 1.02 to 0.84. This indicates that the perfect-correlation variance formula, while an improvement over the independent-guard formula, increasingly underestimates the true ensemble variability as heterogeneity grows and guard estimates diverge more widely.

TGEP achieved the lowest RMSE at high heterogeneity (0.188 at tau^2 = 0.50 vs. REML 0.195; 0.257 at tau^2 = 1.00 vs. REML 0.274), indicating that the ensemble weighting effectively combines guard estimates in high-heterogeneity settings even though the SE is underestimated.

#### Publication Bias: Significance-Based

Under significance-based publication bias (Scenarios 10-14), all standard methods showed substantial coverage deterioration with increasing bias severity. At strong bias with k = 20 (Scenario 12), REML coverage dropped to 50.8% and HKSJ to 57.4%, with TGEP at 48.0%. Under severe bias (Scenario 13), REML coverage was 32.2%, HKSJ 39.2%, and TGEP 27.6%.

PET-PEESE was clearly the best method for maintaining coverage under significance-based publication bias, achieving 88.8-90.0% across strong-to-severe bias scenarios. Trim-and-fill also outperformed REML and TGEP at higher bias levels (42.2% vs. 32.2% and 27.6% under severe bias, respectively).

TGEP's coverage disadvantage relative to REML under significance-based bias reflects two compounding factors: (1) the LOO-CV stacking optimizes prediction of the observed (biased) data, which may assign insufficient weight to the SWA bias-correction guard; and (2) the analytical SE underestimates the true variability. However, TGEP did achieve modestly lower bias than REML (0.131 vs. 0.141 at k = 10 strong; 0.141 vs. 0.150 at k = 20 strong; 0.201 vs. 0.203 at k = 20 severe) and lower RMSE (0.188 vs. 0.194; 0.168 vs. 0.175; 0.220 vs. 0.222).

#### Publication Bias: One-Sided

Under one-sided publication bias (Scenarios 15-16), TGEP showed its clearest point-estimation advantage. At k = 10, TGEP bias was 0.044 vs. REML 0.063 (30% reduction), with RMSE 0.097 vs. 0.106 (8.5% improvement). At k = 20, TGEP bias was 0.069 vs. REML 0.091 (24% reduction), with RMSE 0.107 vs. 0.121 (11.6% improvement). These represent the largest relative improvements observed for TGEP across all scenarios.

Coverage was similar between TGEP and REML (87.8% vs. 87.4% at k = 10; 79.4% vs. 77.8% at k = 20), with HKSJ again providing the best coverage (91.4% and 82.2%). Trim-and-fill performed poorly under one-sided bias (84.8% and 70.0%), as the funnel-plot symmetry assumption is particularly violated under directional selection.

#### Combined Heterogeneity and Publication Bias

Under combined high heterogeneity and significance-based bias (Scenario 17), TGEP's coverage was notably low (71.4%) compared to REML (82.2%) and HKSJ (87.4%). The SE ratio was 0.77, the lowest observed across all scenarios, indicating that the combination of high heterogeneity and publication bias creates the most challenging conditions for TGEP's analytical SE. Nevertheless, TGEP achieved the lowest RMSE among the three core methods (0.220 vs. REML 0.227).

#### SE Calibration

Across all scenarios, TGEP's SE ratio ranged from 0.77 to 1.06. SE was well-calibrated (ratio 0.94-1.06) in null scenarios with moderate heterogeneity and in the low-heterogeneity condition. SE was progressively underestimated (ratio 0.77-0.89) as heterogeneity or publication bias severity increased. By comparison, REML SE ratios ranged from 0.92-1.03, and HKSJ SE ratios from 0.86-1.02, both showing better calibration. This SE underestimation is the primary driver of TGEP's coverage shortfall.

### Empirical Validation

Table 4 presents the results of applying TGEP and REML to 15 Cochrane datasets.

**Table 4. Empirical Validation: TGEP vs. REML on 15 Cochrane Datasets**

| Dataset | K | REML Est. | TGEP Est. | Abs. Diff | REML SE | TGEP SE | SE Ratio |
|---------|---|-----------|-----------|-----------|---------|---------|----------|
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

Across the 15 datasets, TGEP estimates differed from REML by a mean absolute difference of 0.052 (median 0.032) on the standardized scale. The mean SE ratio (TGEP SE / REML SE) was 2.23, substantially wider than the simulation SE ratio within TGEP alone (0.77-1.06), because the empirical SE ratio here compares TGEP against REML rather than TGEP's model SE against its own empirical SE. The wider TGEP SEs in the empirical validation largely reflect the conservative perfect-correlation variance formula applied across three guards whose individual variances may differ substantially from the REML variance.

The largest point-estimate difference occurred for CD013844 (K = 12), where TGEP pulled the estimate toward zero (from -0.495 to -0.264), suggesting the GRMA or SWA guards identified influential studies or potential publication bias. In datasets with minimal heterogeneity (e.g., CD015140), TGEP and REML produced nearly identical results.

## Discussion

### Summary of Findings

We proposed TGEP, a frequentist ensemble framework for meta-analysis that integrates three specialized guard estimators via LOO-CV stacking. In a comprehensive simulation study comparing five methods across 18 scenarios, TGEP demonstrated two principal strengths and one important limitation.

**Strengths.** First, TGEP produced less biased point estimates than REML under one-sided publication bias, reducing absolute bias by 24-30% and RMSE by 9-12%. This advantage was most pronounced when the selection mechanism was directional rather than significance-based. Second, TGEP achieved lower RMSE than REML in high-heterogeneity scenarios (tau^2 = 0.50-1.00) and in null-effect settings, suggesting that the ensemble weighting can improve estimation precision even when publication bias is not the primary concern.

**Limitation.** TGEP's confidence interval coverage was below nominal in most scenarios (86-95% in unbiased settings, 28-88% under publication bias), consistently worse than both REML and HKSJ. The primary cause was SE underestimation: TGEP's analytical SE ratio ranged from 0.77 to 1.06 across scenarios, with the most severe underestimation occurring under combined heterogeneity and bias. This limits TGEP's utility for formal statistical inference.

### Comparison with Existing Methods

Our simulation provides a comparative assessment of five commonly used approaches:

**HKSJ** provided the most reliable coverage across all scenarios (82-96%), confirming its role as the preferred method for inference when publication bias is not the primary concern [6]. Its t-distribution correction addresses SE underestimation effectively, but it does not correct the biased point estimate.

**PET-PEESE** was the most effective method for maintaining coverage under significance-based publication bias (88-90%), consistent with prior comparisons [5]. However, PET-PEESE exhibited substantially inflated RMSE in unbiased scenarios (2-4x higher than REML), reflecting the cost of fitting unnecessary bias-correction parameters. This overcorrection is a well-known limitation [9].

**Trim-and-fill** provided moderate coverage improvements under significance-based bias (42-58% vs. REML 32-51% under strong/severe bias) but performed poorly under one-sided bias (70-85%), where the funnel-plot symmetry assumption is most violated.

**TGEP** occupies a complementary niche: it provides modest point-estimation improvements without the overcorrection risk of PET-PEESE. In unbiased settings, TGEP RMSE was within 3-5% of REML, whereas PET-PEESE RMSE was 300-400% higher. Under one-sided bias, TGEP achieved the largest bias reductions among all methods that do not assume a specific bias mechanism.

### The LOO-CV Limitation

A fundamental tension in TGEP's design is that LOO-CV optimizes prediction of the observed data, which may itself be biased by selective publication. When the selection mechanism is significance-based, the SWA guard — which up-weights non-significant studies — may actually predict the observed (predominantly significant) data less well than the WRD or GRMA guards that do not attempt bias correction. This can result in the LOO-CV assigning suboptimal weight to the guard most relevant for bias correction.

This limitation is inherent to any stacking approach applied to biased data and explains why TGEP's bias correction is modest compared to methods like PET-PEESE that directly model the bias mechanism. Under one-sided bias, where the selection is directional rather than significance-based, this tension is less severe, and TGEP's advantage is correspondingly larger.

### SE Underestimation

TGEP's analytical SE assumes perfect positive correlation among the three guards. Despite being the most conservative correlation assumption (yielding the widest CI for a given set of guard variances), the resulting SE still underestimates the empirical variability. This occurs because the guard-specific variances V_g themselves underestimate the total uncertainty of the ensemble: they do not account for the additional variability introduced by the data-dependent LOO-CV weight selection. When the stacking weights alpha_g vary across bootstrap samples, the ensemble estimate acquires additional variability that the analytical formula does not capture.

Bootstrap SE (B = 200, available as a default option) provides a more accurate variance estimate that accounts for weight uncertainty. However, bootstrap SE was not used in the simulation study for computational reasons (500 replications x 18 scenarios x 200 bootstrap samples = 1.8 million TGEP evaluations). Future work should evaluate whether bootstrap SE resolves the coverage shortfall.

### Strengths of the Framework

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

6. **Simulation scope**: We evaluated N_sim = 500 replications per scenario, which provides MCSE of coverage of approximately 0.01-0.02. While adequate for identifying the magnitude of coverage differences, larger simulations would further refine the estimates.

7. **Number of guards**: The three-guard architecture is fixed. Whether the optimal number of guards differs across settings is unknown.

8. **Cochrane dataset selection**: The 15 empirical validation datasets were selected for availability and diversity, not by a systematic protocol. Future work should evaluate TGEP on a pre-registered set of datasets with known bias characteristics.

### Recommendations for Practice

Based on our findings, we recommend TGEP as a **complementary sensitivity analysis** rather than a replacement for established methods. Specifically:

- Use **HKSJ** as the primary method for inference (confidence intervals and hypothesis tests), given its consistently near-nominal coverage.
- Use **PET-PEESE** when funnel plot asymmetry or other evidence suggests significance-based publication bias, noting its overcorrection risk in unbiased settings.
- Use **TGEP** to provide a bias-reduced point estimate and to examine guard weight diagnostics as informal indicators of potential bias or outlier contamination. When TGEP's SWA guard receives high weight, this suggests the data may be affected by publication bias; when GRMA receives high weight, outlier or small-study effects may be present.
- Report individual guard estimates and weights alongside the ensemble result to support transparency.

The default temperature T = 1.0 is suitable for most applications. Bootstrap SE (n_boot >= 200) should be used when confidence intervals from TGEP are needed.

## Conclusions

TGEP provides a frequentist ensemble framework for meta-analysis that produces less biased point estimates under publication bias, particularly one-sided selection, without the overcorrection risk of parametric bias-correction methods. However, its analytical standard error underestimates the true variability, limiting its utility for formal inference. We recommend TGEP as a complementary sensitivity analysis alongside HKSJ for coverage and PET-PEESE for bias assessment. The method is implemented as an open-source R package and requires no specification of priors, selection functions, or parametric bias models. Future work should address the SE calibration through bootstrap validation or hybrid approaches combining TGEP point estimates with HKSJ-style variance correction.

## References

[1] Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. Chichester: John Wiley & Sons; 2009. doi:10.1002/9780470743386

[2] Viechtbauer W. Conducting meta-analyses in R with the metafor package. J Stat Softw. 2010;36(3):1-48. doi:10.18637/jss.v036.i03

[3] Rothstein HR, Sutton AJ, Borenstein M, editors. Publication Bias in Meta-Analysis: Prevention, Assessment and Adjustments. Chichester: John Wiley & Sons; 2005. doi:10.1002/0470870168

[4] Viechtbauer W, Cheung MW-L. Outlier and influence diagnostics for meta-analysis. Res Synth Methods. 2010;1(2):112-125. doi:10.1002/jrsm.11

[5] Carter EC, Schonbrodt FD, Gervais WM, Hilgard J. Correcting for bias in psychology: A comparison of meta-analytic methods. Adv Methods Pract Psychol Sci. 2019;2(2):115-144. doi:10.1177/2515245919847196

[6] Hartung J, Knapp G. A refined method for the meta-analysis of controlled clinical trials with binary outcome. Stat Med. 2001;20(24):3875-3889. doi:10.1002/sim.1009

[7] Duval S, Tweedie R. Trim and fill: A simple funnel-plot-based method of testing and adjusting for publication bias in meta-analysis. Biometrics. 2000;56(2):455-463. doi:10.1111/j.0006-341X.2000.00455.x

[8] Stanley TD, Doucouliagos H. Meta-regression approximations to reduce publication selection bias. Res Synth Methods. 2014;5(1):60-78. doi:10.1002/jrsm.1095

[9] Stanley TD. Limitations of PET-PEESE and other meta-analysis methods. Soc Psychol Personal Sci. 2017;8(5):581-591. doi:10.1177/1948550617693062

[10] Vevea JL, Hedges LV. A general linear model for estimating effect size in the presence of publication bias. Psychometrika. 1995;60(3):419-435. doi:10.1007/BF02294384

[11] Maier M, Bartos F, Wagenmakers E-J. Robust Bayesian meta-analysis: Addressing publication bias with model-averaging. Psychol Methods. 2023;28(1):107-122. doi:10.1037/met0000405

[12] Wolpert DH. Stacked generalization. Neural Netw. 1992;5(2):241-259. doi:10.1016/S0893-6080(05)80023-1

[13] Deng J-L. Introduction to grey system theory. J Grey Syst. 1989;1(1):1-24.

[14] Hastings C, Mosteller F, Tukey JW, Winsor CP. Low moments for small samples: A comparative study of order statistics. Ann Math Stat. 1947;18(3):413-426. doi:10.1214/aoms/1177730388

[15] Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate statistical methods. Stat Med. 2019;38(11):2074-2102. doi:10.1002/sim.8086

[16] R Core Team. R: A language and environment for statistical computing. Vienna: R Foundation for Statistical Computing; 2025. https://www.R-project.org/

---

## Supporting Information

**S1 Table.** Complete simulation results for all 18 scenarios and 5 methods, including bias, RMSE, coverage, MCSE, SE ratio, and mean number of studies after selection.

**S2 File.** R package source code for TGEP (TGEP.R), including all guard implementations, the LOO-CV stacking ensemble, and diagnostic plotting functions.

**S3 File.** Simulation script (run_simulation.R) with complete ADEMP specification, data-generating mechanisms, and performance measure calculations for all five methods.

**S4 File.** Cochrane dataset identifiers and empirical validation results (Real_World_Impact_Summary.txt).

**S5 File.** Test suite (test_tgep.R) with 23 unit tests covering basic functionality, edge cases, and numerical stability.

---

## Author Information

### Affiliations
^1 Independent Researcher

### Corresponding Author
Mahmood Ul Hassan (mahmood.hassan@example.com)

### Author Contributions (CRediT)
**Conceptualization:** MUH. **Methodology:** MUH. **Software:** MUH. **Validation:** MUH. **Formal Analysis:** MUH. **Investigation:** MUH. **Data Curation:** MUH. **Writing - Original Draft:** MUH. **Writing - Review & Editing:** MUH. **Visualization:** MUH.

### Data Availability Statement
All simulation code and results are provided as Supporting Information. The Cochrane datasets used for empirical validation are identified by their Cochrane review identifiers (CD-numbers) and are publicly available through the Cochrane Library (https://www.cochranelibrary.com/). The TGEP R package is available at [ZENODO_DOI_PLACEHOLDER].

### Funding
The author received no specific funding for this work.

### Competing Interests
The author declares no competing interests.

### Ethics Statement
This study used simulated data and publicly available summary statistics from published Cochrane reviews. No individual patient data were used. No ethics approval was required.
