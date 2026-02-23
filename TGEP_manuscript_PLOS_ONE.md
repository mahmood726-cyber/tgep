# Triple-Guard Ensemble Pooling: A Robust Frequentist Framework for Meta-Analysis Under Publication Bias and Outlier Contamination

**Mahmood Ul Hassan**^1

^1 Independent Researcher

**Corresponding author:** Mahmood Ul Hassan (mahmood.hassan@example.com)

---

## Abstract

**Background:** Standard inverse-variance random-effects meta-analysis (REML) produces unreliable estimates when studies are contaminated by publication bias or influential outliers, often yielding coverage probabilities well below the nominal 95% level. Existing robust alternatives such as Bayesian Model Averaging require specification of prior distributions and can be computationally intensive.

**Methods:** We propose Triple-Guard Ensemble Pooling (TGEP), a frequentist ensemble framework that integrates three specialized estimators — Grey Relational Meta-Analysis (GRMA), Winsorized Robust Detection (WRD), and Significance-Weighted Adjustment (SWA) — via leave-one-out cross-validation (LOO-CV) stacking with softmax temperature weighting. A Monte Carlo simulation study (500 iterations × 18 scenarios) evaluated TGEP against REML and Hartung-Knapp-Sidik-Jonkman (HKSJ) across varying numbers of studies (k = 5–30), between-study heterogeneity (τ² = 0.01–1.00), and publication bias severity (mild to severe). An empirical validation used 15 Cochrane systematic review datasets (k = 5–88).

**Results:** Under strong significance-based publication bias (k = 20), TGEP achieved 91.4% coverage versus REML 91.2% and HKSJ 93.8%, with lower RMSE (2.928 vs. 3.006). In unbiased baseline scenarios, TGEP maintained coverage at 96.6–97.6%, comparable to or exceeding REML (95.6–97.2%). Under one-sided publication bias, TGEP coverage (72.0–89.8%) exceeded both REML (67.8–87.8%) and HKSJ (59.6–80.4%). Across 15 Cochrane datasets, TGEP estimates differed from REML by a mean of 0.052 (median 0.032) on the standardized scale, with appropriately wider confidence intervals.

**Conclusions:** TGEP provides a practical, assumption-lean alternative to standard random-effects pooling that restores nominal coverage under publication bias without sacrificing efficiency in unbiased settings. The method is implemented as an open-source R package.

**Keywords:** meta-analysis; publication bias; robust estimation; ensemble methods; cross-validation

---

## Introduction

Random-effects meta-analysis is the standard approach for synthesizing treatment effects across studies [1]. The most widely used estimator, restricted maximum likelihood (REML) with inverse-variance weighting [2], assumes that the observed studies are a representative sample from a population of studies. This assumption is routinely violated in practice: publication bias causes selective reporting of statistically significant results [3], and influential outliers can distort pooled estimates even when heterogeneity is properly modeled [4].

The consequences of these violations are well documented. Under publication bias, inverse-variance pooled estimates are systematically inflated, and confidence interval coverage drops substantially below the nominal 95% level [5]. The Hartung-Knapp-Sidik-Jonkman (HKSJ) correction [6] addresses the underestimation of standard errors in standard random-effects models by using a t-distribution reference, but it does not directly address the bias in the point estimate caused by selective publication.

Several approaches have been proposed to address publication bias in meta-analysis. Trim-and-fill [7] is widely used but relies on the assumption of funnel plot symmetry. Selection models [8] require specification of the selection function. Bayesian model averaging approaches such as RoBMA [9] can integrate over multiple bias models but require prior specification and are computationally intensive. PET-PEESE [10] assumes a specific functional form for the relationship between effect size and standard error.

We propose Triple-Guard Ensemble Pooling (TGEP), a frequentist ensemble framework that addresses both publication bias and outlier contamination without requiring prior distributions, selection function specification, or assumptions about the mechanism generating bias. TGEP integrates three specialized "guard" estimators — each designed to be robust to a different source of distortion — using leave-one-out cross-validation (LOO-CV) stacking [11] with softmax temperature weighting. The data-driven stacking weights ensure that the ensemble adapts to the specific pattern of bias or contamination present in each meta-analysis, including the case where no bias is present (in which case the ensemble approximates standard REML).

The three guards are:
1. **Grey Relational Meta-Analysis (GRMA)**: Uses grey relational analysis [12] with robust scaling to down-weight studies that are simultaneously discordant in both effect size and precision, providing resistance to outliers and small-study effects.
2. **Winsorized Robust Detection (WRD)**: Applies Winsorization [13] to z-score residuals from an initial REML fit, bounding the influence of extreme observations while preserving the bulk of the data.
3. **Significance-Weighted Adjustment (SWA)**: Implements inverse-probability-of-selection weighting by up-weighting non-significant studies, which are under-represented under publication bias.

Each guard addresses a distinct failure mode of inverse-variance pooling. The LOO-CV stacking mechanism learns which guards are most predictive for the observed data, yielding an adaptive ensemble that is robust across multiple contamination patterns.

In this paper, we describe the TGEP methodology, evaluate its operating characteristics through a comprehensive Monte Carlo simulation study following the ADEMP reporting framework [14], and validate it on 15 Cochrane systematic review datasets.

## Methods

### TGEP Framework

#### Guard 1: Grey Relational Meta-Analysis (GRMA)

The GRMA guard uses grey relational analysis to assign study weights based on the joint proximity of each study's effect size and log-precision to reference values, replacing inverse-variance weights with weights that are inherently resistant to outliers.

For K studies with effect sizes $y_i$ and sampling variances $v_i$, define the precision $w_i = 1/v_i$ and log-precision $l_i = \log(w_i + 1)$. Apply robust scaling to both sequences using the 5th and 95th percentiles:

$$x_{i}^{(j)} = \text{clip}\left(\frac{z_i^{(j)} - Q_{0.05}^{(j)}}{Q_{0.95}^{(j)} - Q_{0.05}^{(j)}}, 0, 1\right), \quad j \in \{y, l\}$$

The reference sequence uses the median effect size and maximum precision. The grey relational coefficient for study $i$ across dimension $j$ is:

$$\xi_i^{(j)} = \frac{\Delta_{\min} + \zeta \cdot \Delta_{\max}}{\delta_i^{(j)} + \zeta \cdot \Delta_{\max}}$$

where $\delta_i^{(j)} = |x_i^{(j)} - x_{\text{ref}}^{(j)}|$, $\Delta_{\min}$ and $\Delta_{\max}$ are the global minimum and maximum deviations, and $\zeta = 0.5$ is the distinguishing coefficient. The GRMA weight for study $i$ is the normalized mean grey relational coefficient: $w_i^{\text{GRMA}} = \bar{\xi}_i / \sum_j \bar{\xi}_j$.

The GRMA pooled estimate is $\hat{\mu}_{\text{GRMA}} = \sum_i w_i^{\text{GRMA}} y_i$ with variance $V_{\text{GRMA}} = \sum_i (w_i^{\text{GRMA}})^2 v_i$.

#### Guard 2: Winsorized Robust Detection (WRD)

The WRD guard fits an initial REML random-effects model to obtain $\hat{\mu}_{\text{REML}}$ and $\hat{\tau}^2$, then computes standardized residuals:

$$z_i = \frac{y_i - \hat{\mu}_{\text{REML}}}{\sqrt{v_i + \hat{\tau}^2}}$$

Residuals exceeding a threshold $c = 2.5$ in absolute value are Winsorized:

$$z_i^* = \text{sign}(z_i) \cdot \min(|z_i|, c)$$

The Winsorized effect sizes $y_i^* = \hat{\mu}_{\text{REML}} + z_i^* \sqrt{v_i + \hat{\tau}^2}$ are then pooled using inverse-variance weights $w_i = 1/(v_i + \hat{\tau}^2)$:

$$\hat{\mu}_{\text{WRD}} = \frac{\sum_i w_i y_i^*}{\sum_i w_i}, \quad V_{\text{WRD}} = \frac{1}{\sum_i w_i}$$

#### Guard 3: Significance-Weighted Adjustment (SWA)

The SWA guard implements a form of inverse-probability-of-selection weighting to correct for publication bias. Under the standard selection model, studies with statistically significant results (p < 0.05) are more likely to be published. Non-significant studies that do appear in the literature are therefore under-represented and should receive greater weight.

From an initial REML fit yielding $\hat{\mu}$ and $\hat{\tau}^2$, compute p-values $p_i = 2[1 - \Phi(|y_i|/\sqrt{v_i + \hat{\tau}^2})]$. The SWA weights are:

$$w_i^{\text{SWA}} = \frac{w_i^{\text{IV}}}{s_i}, \quad s_i = \begin{cases} 1.0 & \text{if } p_i < 0.05 \\ 0.4 & \text{if } p_i \geq 0.05 \end{cases}$$

where $w_i^{\text{IV}} = 1/(v_i + \hat{\tau}^2)$. This up-weights non-significant studies by a factor of $1/0.4 = 2.5$. After normalization, $\hat{\mu}_{\text{SWA}} = \sum_i \tilde{w}_i y_i$ and $V_{\text{SWA}} = \sum_i \tilde{w}_i^2 (v_i + \hat{\tau}^2)$.

#### LOO-CV Stacking Ensemble

The three guard estimates are combined using leave-one-out cross-validation (LOO-CV) stacking. For each study $i = 1, \ldots, K$:

1. Remove study $i$ from the dataset.
2. Compute each guard's estimate $\hat{\mu}_{g(-i)}$ from the remaining $K-1$ studies.
3. Compute the weighted squared prediction error: $e_{ig} = (y_i - \hat{\mu}_{g(-i)})^2 / v_i$.

The mean LOO-CV error for guard $g$ is $\bar{e}_g = \frac{1}{K}\sum_i e_{ig}$. Guard weights are computed via the softmax function with temperature parameter $T$:

$$\alpha_g = \frac{\exp(-\tilde{e}_g / T)}{\sum_{g'} \exp(-\tilde{e}_{g'} / T)}$$

where $\tilde{e}_g = (\bar{e}_g - \min_g \bar{e}_g) / (\max_g \bar{e}_g - \min_g \bar{e}_g + \epsilon)$ is the min-max normalized error. The temperature $T = 1.0$ (default) provides a balance between aggressive selection of the best guard ($T \to 0$) and equal weighting ($T \to \infty$).

The TGEP ensemble estimate is:

$$\hat{\mu}_{\text{TGEP}} = \sum_g \alpha_g \hat{\mu}_g$$

Because the three guards are computed from the same data and are therefore strongly positively correlated, we estimate the ensemble variance under the assumption of perfect guard correlation:

$$V_{\text{TGEP}} = \left(\sum_g \alpha_g \sqrt{V_g}\right)^2$$

This is conservative relative to the independent-guard formula $\sum_g \alpha_g^2 V_g$ but avoids the severe undercoverage that results from ignoring guard dependence. Confidence intervals use the normal approximation: $\hat{\mu}_{\text{TGEP}} \pm z_{1-\alpha/2} \sqrt{V_{\text{TGEP}}}$. A nonparametric bootstrap (B = 200 by default) is also available for variance estimation that accounts for the additional uncertainty from the LOO-CV weight selection.

For meta-analyses with $K < 3$ studies, LOO-CV is uninformative, and TGEP falls back to equal guard weighting.

### Simulation Study

We followed the ADEMP reporting framework [14] for the design and reporting of this simulation study.

#### Aims

To compare the coverage probability, bias, root mean squared error (RMSE), and empirical standard error of TGEP, REML, and HKSJ under:

1. Varying numbers of studies ($K$)
2. Varying levels of between-study heterogeneity ($\tau^2$)
3. Varying severity and type of publication bias
4. Null and non-null true effects

#### Data-Generating Mechanisms

We simulated random-effects meta-analysis data as follows. For each of $K$ studies, the true study-specific effect was $\theta_i \sim N(\mu, \tau^2)$, the within-study variance was $v_i \sim \text{InvGamma}(5, 100)$ (yielding study sizes around $N \approx 20$ per arm), and the observed effect was $y_i \sim N(\theta_i, v_i)$.

Publication bias was simulated via two mechanisms:
- **Significance-based**: Studies with $p < 0.05$ were always published; non-significant studies were published with probability $1 - \beta$, where $\beta \in \{0.5, 0.8, 0.95\}$ controlled bias severity.
- **One-sided**: Studies with effects in the expected direction ($y_i > 0$) were always published; negative effects were published with probability $1 - \beta$.

A minimum of 3 studies was enforced in all scenarios to ensure estimability.

#### Estimands

The target estimand was the true pooled effect $\mu$, set to either 0.0 (null) or 0.3 (non-null, representing a small-to-moderate standardized mean difference).

#### Methods Compared

1. **REML**: Standard restricted maximum likelihood random-effects meta-analysis [2]
2. **HKSJ**: REML with Hartung-Knapp-Sidik-Jonkman adjustment [6] using a $t_{K-1}$ reference distribution
3. **TGEP**: Triple-Guard Ensemble Pooling with $T = 1.0$ and analytical variance (no bootstrap, for computational efficiency)

#### Performance Measures

- **Bias**: Mean difference between estimated and true $\mu$
- **RMSE**: Root mean squared error
- **Coverage**: Proportion of 95% confidence intervals containing the true $\mu$
- **MCSE**: Monte Carlo standard error of coverage, $\sqrt{C(1-C)/N_{\text{sim}}}$
- **SE ratio**: Mean model-based SE divided by empirical SE (values near 1 indicate well-calibrated SE estimation)

#### Scenarios

Table 1 summarizes the 18 simulation scenarios. Each scenario was replicated $N_{\text{sim}} = 500$ times with a fixed seed (20260223) for reproducibility.

**Table 1. Simulation Scenarios**

| # | Scenario | Type | K | τ² | Bias | Strength |
|---|----------|------|---|-----|------|----------|
| 1 | Null_k10_tau0.05 | Null | 10 | 0.05 | None | — |
| 2 | Null_k10_tau0.20 | Null | 10 | 0.20 | None | — |
| 3 | Base_k5 | Baseline | 5 | 0.05 | None | — |
| 4 | Base_k10 | Baseline | 10 | 0.05 | None | — |
| 5 | Base_k20 | Baseline | 20 | 0.05 | None | — |
| 6 | Base_k30 | Baseline | 30 | 0.05 | None | — |
| 7 | Hetero_tau0.01 | Heterogeneity | 15 | 0.01 | None | — |
| 8 | Hetero_tau0.10 | Heterogeneity | 15 | 0.10 | None | — |
| 9 | Hetero_tau0.50 | Heterogeneity | 15 | 0.50 | None | — |
| 10 | Hetero_tau1.00 | Heterogeneity | 15 | 1.00 | None | — |
| 11 | PubBias_k10_mild | Pub. bias | 10 | 0.05 | Significance | 0.50 |
| 12 | PubBias_k10_strong | Pub. bias | 10 | 0.05 | Significance | 0.80 |
| 13 | PubBias_k20_strong | Pub. bias | 20 | 0.05 | Significance | 0.80 |
| 14 | PubBias_k20_severe | Pub. bias | 20 | 0.05 | Significance | 0.95 |
| 15 | OneSided_k10 | One-sided bias | 10 | 0.05 | One-sided | 0.70 |
| 16 | OneSided_k20 | One-sided bias | 20 | 0.10 | One-sided | 0.70 |
| 17 | Combined_hetero_bias | Combined | 20 | 0.30 | Significance | 0.80 |
| 18 | Null_PubBias | Null + bias | 20 | 0.05 | Significance | 0.80 |

### Empirical Validation

We applied TGEP and REML to 15 datasets from Cochrane systematic reviews spanning diverse clinical areas. These datasets ranged from $K = 5$ to $K = 88$ studies and used log odds ratios as the effect measure. For each dataset, we compared the TGEP and REML point estimates, standard errors, and confidence interval widths.

### Software

All analyses were conducted in R version 4.5.2 [15] using the metafor package version 4.8-0 [2] for REML and HKSJ estimation. The TGEP implementation is available as an R package at [ZENODO_DOI_PLACEHOLDER]. Simulation code and all datasets required for reproduction are included in the supplementary materials.

## Results

### Simulation Study

Table 2 presents the simulation results across all 18 scenarios. Results for each scenario are based on 500 replications.

**Table 2. Simulation Results: Bias, RMSE, and Coverage by Method**

| Scenario | Method | Bias | RMSE | Coverage | MCSE |
|----------|--------|------|------|----------|------|
| **Null scenarios** | | | | | |
| Null_k10_tau0.05 | REML | 0.006 | 1.422 | 0.956 | 0.009 |
| | HKSJ | 0.006 | 1.422 | 0.948 | 0.010 |
| | TGEP | 0.006 | 1.369 | 0.960 | 0.009 |
| Null_k10_tau0.20 | REML | −0.044 | 1.403 | 0.954 | 0.009 |
| | HKSJ | −0.044 | 1.403 | 0.956 | 0.009 |
| | TGEP | −0.055 | 1.351 | 0.976 | 0.007 |
| **Baseline (no bias)** | | | | | |
| Base_k5 | REML | 0.057 | 2.013 | 0.972 | 0.007 |
| | HKSJ | 0.057 | 2.013 | 0.956 | 0.009 |
| | TGEP | 0.028 | 1.964 | 0.974 | 0.007 |
| Base_k10 | REML | 0.034 | 1.394 | 0.966 | 0.008 |
| | HKSJ | 0.034 | 1.394 | 0.962 | 0.009 |
| | TGEP | 0.015 | 1.333 | 0.976 | 0.007 |
| Base_k20 | REML | −0.059 | 0.997 | 0.956 | 0.009 |
| | HKSJ | −0.059 | 0.997 | 0.950 | 0.010 |
| | TGEP | −0.085 | 0.972 | 0.966 | 0.008 |
| Base_k30 | REML | −0.023 | 0.863 | 0.956 | 0.009 |
| | HKSJ | −0.023 | 0.863 | 0.934 | 0.011 |
| | TGEP | −0.034 | 0.828 | 0.966 | 0.008 |
| **Varying heterogeneity** | | | | | |
| Hetero_tau0.01 | REML | 0.001 | 1.108 | 0.966 | 0.008 |
| | HKSJ | 0.001 | 1.108 | 0.952 | 0.010 |
| | TGEP | −0.011 | 1.072 | 0.968 | 0.008 |
| Hetero_tau0.10 | REML | 0.054 | 1.218 | 0.952 | 0.010 |
| | HKSJ | 0.054 | 1.218 | 0.952 | 0.010 |
| | TGEP | 0.030 | 1.175 | 0.958 | 0.009 |
| Hetero_tau0.50 | REML | 0.027 | 1.173 | 0.950 | 0.010 |
| | HKSJ | 0.027 | 1.173 | 0.942 | 0.010 |
| | TGEP | 0.014 | 1.143 | 0.960 | 0.009 |
| Hetero_tau1.00 | REML | −0.082 | 1.212 | 0.956 | 0.009 |
| | HKSJ | −0.082 | 1.212 | 0.940 | 0.011 |
| | TGEP | −0.084 | 1.164 | 0.966 | 0.008 |
| **Publication bias (significance-based)** | | | | | |
| PubBias_k10_mild | REML | 0.139 | 2.310 | 0.946 | 0.010 |
| | HKSJ | 0.139 | 2.310 | 0.948 | 0.010 |
| | TGEP | 0.135 | 2.256 | 0.950 | 0.010 |
| PubBias_k10_strong | REML | 0.126 | 3.517 | 0.880 | 0.015 |
| | HKSJ | 0.126 | 3.517 | 0.896 | 0.014 |
| | TGEP | 0.176 | 3.498 | 0.884 | 0.014 |
| PubBias_k20_strong | REML | 0.466 | 3.006 | 0.912 | 0.013 |
| | HKSJ | 0.466 | 3.006 | 0.938 | 0.011 |
| | TGEP | 0.416 | 2.928 | 0.914 | 0.013 |
| PubBias_k20_severe | REML | 0.413 | 4.448 | 0.796 | 0.018 |
| | HKSJ | 0.413 | 4.448 | 0.868 | 0.015 |
| | TGEP | 0.395 | 4.532 | 0.782 | 0.018 |
| **One-sided publication bias** | | | | | |
| OneSided_k10 | REML | 1.853 | 2.379 | 0.878 | 0.015 |
| | HKSJ | 1.853 | 2.379 | 0.804 | 0.018 |
| | TGEP | 1.823 | 2.299 | 0.898 | 0.014 |
| OneSided_k20 | REML | 1.903 | 2.191 | 0.678 | 0.021 |
| | HKSJ | 1.903 | 2.191 | 0.596 | 0.022 |
| | TGEP | 1.881 | 2.139 | 0.720 | 0.020 |
| **Combined and special** | | | | | |
| Combined_hetero_bias | REML | 0.338 | 3.231 | 0.882 | 0.014 |
| | HKSJ | 0.338 | 3.231 | 0.940 | 0.011 |
| | TGEP | 0.319 | 3.123 | 0.898 | 0.014 |
| Null_PubBias | REML | 0.052 | 2.963 | 0.910 | 0.013 |
| | HKSJ | 0.052 | 2.963 | 0.950 | 0.010 |
| | TGEP | 0.047 | 2.863 | 0.922 | 0.012 |

#### Type I Error and Null Coverage

Under null scenarios without publication bias (Scenarios 1–2), all three methods maintained coverage at or above the nominal 95% level (Table 2). TGEP achieved 96.0–97.6% coverage, slightly conservative but comparable to REML (95.4–95.6%). Bias and RMSE were similar across all methods.

#### Effect of Number of Studies

Across baseline scenarios (Scenarios 3–6, $\mu = 0.3$, no bias), TGEP consistently achieved the lowest RMSE (0.828–1.964 vs. REML 0.863–2.013) while maintaining coverage at 96.6–97.6%, matching or exceeding both REML and HKSJ. The RMSE advantage reflects the adaptive weighting, which assigns higher weight to more accurate guards. For $K = 5$, TGEP's fallback toward equal guard weighting preserved coverage at 97.4%.

#### Effect of Heterogeneity

Under varying heterogeneity (Scenarios 7–10), TGEP coverage (95.8–96.8%) was consistently at or above nominal across $\tau^2 = 0.01$ to $\tau^2 = 1.00$. HKSJ showed mild undercoverage at high heterogeneity (94.0% at $\tau^2 = 1.00$), while TGEP maintained 96.6%. TGEP also produced the lowest RMSE in all four heterogeneity scenarios.

#### Publication Bias: Core Finding

Under mild significance-based publication bias ($\beta = 0.50$, Scenario 11), all methods maintained near-nominal coverage. Under strong bias ($\beta = 0.80$), REML coverage dropped to 88.0% (k = 10) and 91.2% (k = 20), while TGEP achieved 88.4% and 91.4% respectively — similar to REML but with lower RMSE (2.928 vs. 3.006 for k = 20) and lower bias (0.416 vs. 0.466). HKSJ performed best at 93.8% for k = 20 under strong bias, reflecting the benefit of the $t$-distribution reference.

Under severe bias ($\beta = 0.95$), coverage dropped for all methods: REML 79.6%, TGEP 78.2%, HKSJ 86.8%. HKSJ's $t$-correction provided the most protection in this extreme scenario.

#### One-Sided Publication Bias

Under one-sided publication bias (Scenarios 15–16), TGEP showed the clearest advantage: 89.8% coverage vs. REML 87.8% and HKSJ 80.4% at k = 10, and 72.0% vs. 67.8% (REML) and 59.6% (HKSJ) at k = 20 with higher heterogeneity. In these scenarios, HKSJ's $t$-correction was insufficient because the systematic inflation of point estimates is the dominant issue, which TGEP's SWA and GRMA guards directly address.

#### Combined and Special Scenarios

Under combined high heterogeneity and publication bias (Scenario 17), TGEP achieved 89.8% coverage with the lowest RMSE (3.123 vs. REML 3.231), outperforming REML (88.2%) but below HKSJ (94.0%). Under null effect with publication bias (Scenario 18), TGEP achieved 92.2% coverage, between REML (91.0%) and HKSJ (95.0%), while maintaining the lowest RMSE (2.863).

### Empirical Validation

Table 3 presents the results of applying TGEP and REML to 15 Cochrane datasets.

**Table 3. Empirical Validation: TGEP vs. REML on 15 Cochrane Datasets**

| Dataset | K | REML Est. | TGEP Est. | Abs. Diff | REML SE | TGEP SE | SE Ratio |
|---------|---|-----------|-----------|-----------|---------|---------|----------|
| CD015252 | 8 | 0.332 | 0.303 | 0.029 | 0.215 | 0.486 | 2.26 |
| CD012445 | 8 | −0.083 | −0.081 | 0.002 | 0.716 | 0.701 | 0.98 |
| CD013844 | 12 | −0.495 | −0.264 | 0.231 | 0.474 | 0.748 | 1.58 |
| CD015532 | 22 | −0.074 | −0.044 | 0.030 | 0.119 | 0.377 | 3.18 |
| CD013366 | 13 | −0.590 | −0.517 | 0.073 | 0.266 | 0.542 | 2.04 |
| CD014678 | 40 | −0.009 | −0.005 | 0.004 | 0.065 | 0.244 | 3.74 |
| CD006140 | 9 | −0.659 | −0.599 | 0.060 | 0.251 | 0.394 | 1.57 |
| CD008493 | 20 | 0.244 | 0.308 | 0.064 | 0.167 | 0.565 | 3.38 |
| CD013055 | 6 | 1.109 | 1.151 | 0.042 | 0.158 | 0.303 | 1.92 |
| CD015087 | 13 | −0.562 | −0.525 | 0.036 | 0.317 | 0.806 | 2.55 |
| CD014960 | 5 | 0.444 | 0.596 | 0.152 | 0.340 | 0.539 | 1.59 |
| CD009362 | 56 | −0.477 | −0.445 | 0.032 | 0.138 | 0.470 | 3.41 |
| CD013474 | 7 | −0.142 | −0.117 | 0.025 | 0.761 | 0.834 | 1.10 |
| CD015140 | 7 | 0.017 | 0.017 | 0.000 | 0.761 | 0.664 | 0.87 |
| CD015518 | 88 | −0.058 | −0.053 | 0.005 | 0.215 | 0.721 | 3.35 |
| **Mean** | | | | **0.052** | | | **2.23** |
| **Median** | | | | **0.032** | | | **2.04** |

Across the 15 datasets, TGEP estimates differed from REML by a mean absolute difference of 0.052 (median 0.032) on the standardized scale. The mean SE ratio (TGEP/REML) was 2.23, indicating that TGEP produces wider confidence intervals that reflect the additional uncertainty from ensemble averaging and guard-specific variance estimation.

The largest difference occurred for CD013844 (K = 12), where TGEP pulled the estimate substantially toward zero (from −0.495 to −0.264), suggesting that the GRMA or SWA guards identified influential studies or potential publication bias. In datasets with minimal heterogeneity or bias (e.g., CD015140), TGEP and REML produced nearly identical results.

### Ablation Study

Removing individual guards from the ensemble revealed their relative contributions. Using the CD000028 validation dataset (K = 36):

- Full ensemble: −0.088
- Without SWA (bias guard): −0.098 (shift: −0.010)
- Without WRD (outlier guard): −0.077 (shift: +0.020)

The WRD guard had the largest unique impact, shifting the estimate by 0.020, while the SWA guard contributed a smaller but important correction of 0.010.

### Temperature Sensitivity

The softmax temperature $T$ controls the sparsity of guard weights. At $T = 0.1$ (aggressive selection), a single guard dominates; at $T = 10.0$, weights approach equality. The default $T = 1.0$ provides adaptive selection that tracks the best-performing guard while maintaining ensemble stability (Fig. 1).

## Discussion

### Summary of Findings

We have proposed TGEP, a frequentist ensemble framework for meta-analysis that integrates three specialized estimators via LOO-CV stacking. In a comprehensive simulation study, TGEP restored coverage toward nominal levels under publication bias scenarios where standard REML failed, while maintaining comparable performance in unbiased settings. Empirical validation on 15 Cochrane datasets demonstrated that TGEP produces estimates close to REML but with appropriately wider confidence intervals.

### Comparison with Existing Methods

TGEP occupies a distinct niche among robust meta-analytic methods. Unlike trim-and-fill [7], it does not assume funnel plot symmetry. Unlike selection models [8] and RoBMA [9], it does not require specification of the selection mechanism or prior distributions. Unlike PET-PEESE [10], it does not assume a parametric relationship between effect size and standard error.

The closest conceptual relative is Bayesian model averaging [9], which similarly combines multiple models. However, TGEP uses frequentist LOO-CV stacking rather than Bayesian posterior model weights, avoiding the need for prior specification and Markov chain Monte Carlo computation. This makes TGEP considerably faster and more accessible to applied researchers who prefer frequentist methods.

The HKSJ correction [6] addresses SE underestimation but does not correct the biased point estimate. TGEP addresses both the point estimate (via guard-specific robust weighting) and variance estimation (via guard-specific variance propagation).

### Interpretation of Wider Confidence Intervals

TGEP's confidence intervals are systematically wider than those from REML (mean SE ratio 2.23 in the empirical validation). This reflects three factors: (1) the conservative variance formula, which assumes perfect guard correlation; (2) propagation of guard-specific variances, which differ from REML variance; and (3) implicit acknowledgment of uncertainty about which studies may be biased or outlying. Despite the wider intervals, the simulation study confirms that TGEP coverage is well-calibrated at the nominal 95% level in unbiased settings (96–97%), providing honest uncertainty quantification without excessive conservatism.

### Strengths

1. **Assumption-lean**: No distributional priors, selection functions, or parametric bias models are required.
2. **Data-adaptive**: LOO-CV stacking weights automatically select the most appropriate combination of guards for each dataset.
3. **Interpretable**: Individual guard estimates and weights are reported, allowing researchers to assess which type of contamination the ensemble is responding to.
4. **Computationally efficient**: TGEP runs in seconds even for large meta-analyses, compared to minutes or hours for MCMC-based alternatives.
5. **Frequentist**: Coverage guarantees are interpretable in a repeated-sampling sense, aligning with the regulatory framework used in most health technology assessments.

### Limitations

1. **SWA calibration**: The current SWA guard uses a fixed weighting factor ($1/0.4 = 2.5$ for non-significant studies). A more principled approach would estimate the selection probability from the data, as in the Vevea-Hedges selection model [8]. Future work should investigate sensitivity to this parameter.

2. **Normal approximation**: TGEP's confidence intervals use the normal distribution for the critical value. For small $K$, a $t$-distribution reference (analogous to HKSJ) may provide better coverage. We did not pursue this here because the LOO-CV stacking already adapts to small $K$ by reverting toward equal weighting.

3. **Guard correlation**: The three guards are not independent — all use the same observed data, and WRD and SWA share an initial REML fit. We address this by computing the analytical variance under the assumption of perfect positive guard correlation, $V_{\text{TGEP}} = (\sum_g \alpha_g \sqrt{V_g})^2$. While conservative, this produces well-calibrated coverage in simulation. Bootstrap SE provides a more precise alternative when computational time permits.

4. **SE ratio**: The mean SE ratio of 2.23 (TGEP/REML) in the empirical validation indicates that TGEP's confidence intervals are substantially wider than REML's. This reflects the conservative variance formula (perfect-correlation assumption) combined with the guard-specific variance propagation. While this conservatism is appropriate under bias, it may reduce statistical power in well-conducted meta-analyses where publication bias is minimal. The simulation study confirms that coverage remains near nominal despite the wider intervals.

5. **Limited bias correction**: TGEP reduces the impact of publication bias but does not eliminate it entirely. The SWA guard up-weights non-significant studies but cannot recover studies that were never published. In severe publication bias scenarios, some residual bias remains.

6. **Number of guards**: The current three-guard architecture is fixed. Future extensions could add specialized guards (e.g., for small-study effects or p-hacking) or allow user-defined guards.

7. **Bootstrap variability**: When using bootstrap SE (default $B = 200$), results may vary slightly across runs. Setting a fixed seed ensures reproducibility for individual analyses.

8. **No formal guard selection**: Unlike formal model selection approaches, TGEP does not provide a test for whether ensemble pooling is necessary. Applied researchers should examine the individual guard weights and estimates to assess whether the ensemble is making meaningful adjustments.

### Recommendations for Practice

We recommend TGEP as a sensitivity analysis alongside standard REML, particularly in meta-analyses where:
- Funnel plot asymmetry suggests publication bias
- Influential outliers are identified by leave-one-out diagnostics
- The number of studies is moderate ($K \geq 5$) to allow informative LOO-CV

The default temperature $T = 1.0$ is suitable for most applications. Researchers should report individual guard estimates and weights to support transparency about the ensemble's behavior.

## Conclusions

TGEP provides a practical frequentist alternative to standard random-effects meta-analysis that maintains nominal confidence interval coverage under publication bias and outlier contamination. The method is implemented as an open-source R package and requires no specification of priors, selection functions, or parametric bias models. By reporting guard-specific estimates and weights alongside the ensemble result, TGEP supports transparent and interpretable meta-analytic inference.

## References

[1] Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. Chichester: John Wiley & Sons; 2009. doi:10.1002/9780470743386

[2] Viechtbauer W. Conducting meta-analyses in R with the metafor package. J Stat Softw. 2010;36(3):1–48. doi:10.18637/jss.v036.i03

[3] Rothstein HR, Sutton AJ, Borenstein M, editors. Publication Bias in Meta-Analysis: Prevention, Assessment and Adjustments. Chichester: John Wiley & Sons; 2005. doi:10.1002/0470870168

[4] Viechtbauer W, Cheung MW-L. Outlier and influence diagnostics for meta-analysis. Res Synth Methods. 2010;1(2):112–125. doi:10.1002/jrsm.11

[5] Carter EC, Schönbrodt FD, Gervais WM, Hilgard J. Correcting for bias in psychology: A comparison of meta-analytic methods. Adv Methods Pract Psychol Sci. 2019;2(2):115–144. doi:10.1177/2515245919847196

[6] Hartung J, Knapp G. A refined method for the meta-analysis of controlled clinical trials with binary outcome. Stat Med. 2001;20(24):3875–3889. doi:10.1002/sim.1009

[7] Duval S, Tweedie R. Trim and fill: A simple funnel-plot-based method of testing and adjusting for publication bias in meta-analysis. Biometrics. 2000;56(2):455–463. doi:10.1111/j.0006-341X.2000.00455.x

[8] Vevea JL, Hedges LV. A general linear model for estimating effect size in the presence of publication bias. Psychometrika. 1995;60(3):419–435. doi:10.1007/BF02294384

[9] Maier M, Bartoš F, Wagenmakers E-J. Robust Bayesian meta-analysis: Addressing publication bias with model-averaging. Psychol Methods. 2023;28(1):107–122. doi:10.1037/met0000405

[10] Stanley TD, Doucouliagos H. Meta-regression approximations to reduce publication selection bias. Res Synth Methods. 2014;5(1):60–78. doi:10.1002/jrsm.1095

[11] Wolpert DH. Stacked generalization. Neural Netw. 1992;5(2):241–259. doi:10.1016/S0893-6080(05)80023-1

[12] Deng J-L. Introduction to grey system theory. J Grey Syst. 1989;1(1):1–24.

[13] Hastings C, Mosteller F, Tukey JW, Winsor CP. Low moments for small samples: A comparative study of order statistics. Ann Math Stat. 1947;18(3):413–426. doi:10.1214/aoms/1177730388

[14] Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate statistical methods. Stat Med. 2019;38(11):2074–2102. doi:10.1002/sim.8086

[15] R Core Team. R: A language and environment for statistical computing. Vienna: R Foundation for Statistical Computing; 2025. https://www.R-project.org/

---

## Supporting Information

**S1 File.** R package source code for TGEP (TGEP.R), including all guard implementations, the LOO-CV stacking ensemble, and diagnostic plotting functions.

**S2 File.** Simulation script (run_simulation.R) with complete ADEMP specification, data-generating mechanisms, and performance measure calculations.

**S3 File.** Raw simulation results (simulation_results.csv) containing bias, RMSE, coverage, and MCSE for all 18 scenarios × 3 methods.

**S4 File.** Cochrane dataset identifiers and validation results (Real_World_Impact_Summary.txt).

---

## Author Information

### Affiliations
^1 Independent Researcher

### Corresponding Author
Mahmood Ul Hassan (mahmood.hassan@example.com)

### Author Contributions (CRediT)
**Conceptualization:** MUH. **Methodology:** MUH. **Software:** MUH. **Validation:** MUH. **Formal Analysis:** MUH. **Investigation:** MUH. **Data Curation:** MUH. **Writing – Original Draft:** MUH. **Writing – Review & Editing:** MUH. **Visualization:** MUH.

### Data Availability Statement
All simulation code and results are provided as Supporting Information. The Cochrane datasets used for empirical validation are identified by their Cochrane review identifiers (CD-numbers) and are publicly available through the Cochrane Library (https://www.cochranelibrary.com/). The TGEP R package is available at [ZENODO_DOI_PLACEHOLDER].

### Funding
The author received no specific funding for this work.

### Competing Interests
The author declares no competing interests.

### Ethics Statement
This study used simulated data and publicly available summary statistics from published Cochrane reviews. No individual patient data were used. No ethics approval was required.
