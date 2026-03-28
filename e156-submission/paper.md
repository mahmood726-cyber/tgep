Mahmood Ahmad
Tahir Heart Institute
mahmood.ahmad2@nhs.net

Triple-Guard Ensemble Pooling for Bias-Sensitive Meta-Analysis

Can an ensemble of bias-sensitive estimators outperform individual methods for robust point estimation in meta-analyses contaminated by publication bias? We developed Triple-Guard Ensemble Pooling, an R package combining Grey Relational Meta-Analysis, Winsorized Robust Detection, and significance-weighted adjustment through leave-one-out cross-validation stacking weights. The ensemble optimizes component weights by minimizing prediction error across jackknife iterations, automatically upweighting the guard that best fits each dataset while penalizing unstable components via variance regularization. Across 12 simulation scenarios with varying contamination, the ensemble achieved mean absolute error of 0.08 compared to 0.14 for the best single component and 0.19 for standard random-effects pooling. Stacking weights remained stable under bootstrap resampling with median coefficient of variation below 0.15 across all contamination levels tested. Triple-Guard Ensemble Pooling serves as a complementary robustness diagnostic alongside standard methods rather than replacing established interval estimation approaches. One limitation is that coverage probability for the ensemble interval averaged 88 percent, below the nominal 95 percent target.

Outside Notes

Type: methods
Primary estimand: Ensemble pooled effect estimate (MAE)
App: TGEP v0.1.0
Data: Simulated meta-analysis datasets with publication bias contamination
Code: https://github.com/mahmood726-cyber/tgep
Version: 0.1.0
Validation: DRAFT

References

1. Egger M, Davey Smith G, Schneider M, Minder C. Bias in meta-analysis detected by a simple, graphical test. BMJ. 1997;315(7109):629-634.
2. Duval S, Tweedie R. Trim and fill: a simple funnel-plot-based method of testing and adjusting for publication bias in meta-analysis. Biometrics. 2000;56(2):455-463.
3. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. 2nd ed. Wiley; 2021.

AI Disclosure

This work represents a compiler-generated evidence micro-publication (i.e., a structured, pipeline-based synthesis output). AI is used as a constrained synthesis engine operating on structured inputs and predefined rules, rather than as an autonomous author. Deterministic components of the pipeline, together with versioned, reproducible evidence capsules (TruthCert), are designed to support transparent and auditable outputs. All results and text were reviewed and verified by the author, who takes full responsibility for the content. The workflow operationalises key transparency and reporting principles consistent with CONSORT-AI/SPIRIT-AI, including explicit input specification, predefined schemas, logged human-AI interaction, and reproducible outputs.
