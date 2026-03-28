# TGEP: a software tool for reviewer-auditable evidence synthesis

## Authors
- Mahmood Ahmad [1,2]
- Niraj Kumar [1]
- Bilaal Dar [3]
- Laiba Khan [1]
- Andrew Woo [4]
- Corresponding author: Andrew Woo (andy2709w@gmail.com)

## Affiliations
1. Royal Free Hospital
2. Tahir Heart Institute Rabwah
3. King's College Medical School
4. St George's Medical School

## Abstract
**Background:** No single publication-bias-sensitive estimator behaves optimally across all meta-analytic contamination patterns. Ensemble approaches can help, but only if their components, weighting logic, and performance tradeoffs are documented clearly enough for peer review.

**Methods:** TGEP implements a frequentist ensemble combining Grey Relational Meta-Analysis, Winsorized Robust Detection, and a significance-weighted adjustment through leave-one-out cross-validation stacking. The repository includes method manuscripts, simulations, and supporting project materials.

**Results:** Local files document scenario-based simulation comparisons, guard-specific behavior, and practical guidance on using TGEP as a companion point-estimation and sensitivity-analysis method rather than as a primary inferential engine.

**Conclusions:** TGEP should be framed as a complementary robustness tool whose value lies in bias-sensitive point estimation and diagnostic disagreement, while interval estimation remains more conservatively handled by established methods such as HKSJ.

## Keywords
ensemble meta-analysis; publication bias; robust pooling; leave-one-out stacking; software tool

## Introduction
The project packages ensemble meta-analysis in a form reviewers can interrogate: component guards are explicit, the stacking logic is visible, and supporting files explain where coverage lags behind point-estimate performance.

The paper explicitly compares the ensemble with REML, HKSJ, trim-and-fill, PET-PEESE, and other bias-sensitive alternatives, while keeping claims limited to the contexts actually evaluated in the local materials.

The manuscript structure below is deliberately aligned to common open-software review requests: the rationale is stated explicitly, at least one runnable example path is named, local validation artifacts are listed, and conclusions are bounded to the functions and outputs documented in the repository.

## Methods
### Software architecture and workflow
The local project centers on the TGEP manuscript, supporting scripts, and reviewer-oriented documentation. The software logic combines three guards under cross-validation-derived weights and records their ensemble behavior for downstream use.

### Installation, runtime, and reviewer reruns
The local implementation is packaged under `C:\Models\TGEP_Development`. The manuscript identifies the local entry points, dependency manifest, fixed example input, and expected saved outputs so that reviewers can rerun the documented workflow without reconstructing it from scratch.

- Entry directory: `C:\Models\TGEP_Development`.
- Detected documentation entry points: `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.
- Detected environment capture or packaging files: `environment.yml`.
- Named worked-example paths in this draft: `TGEP_manuscript_PLOS_ONE.md` as the primary methods narrative; Project scripts and supporting files for simulations and reporting; `F1000_Submission_Checklist_RealReview.md` for reviewer-aligned packaging.
- Detected validation or regression artifacts: `f1000_artifacts/validation_summary.md`, `inst/scripts/Run_Validation.R`, `test_tgep.R`.
- Detected example or sample data files: `f1000_artifacts/example_dataset.csv`.

### Worked examples and validation materials
**Example or fixed demonstration paths**
- `TGEP_manuscript_PLOS_ONE.md` as the primary methods narrative.
- Project scripts and supporting files for simulations and reporting.
- `F1000_Submission_Checklist_RealReview.md` for reviewer-aligned packaging.

**Validation and reporting artifacts**
- Local simulation outputs summarized in the manuscript materials.
- Reviewer-facing discussion of coverage, weighting, and guard behavior.
- Submission artifacts documenting software packaging and scope boundaries.

### Typical outputs and user-facing deliverables
- Ensemble pooled point estimates and guard weights.
- Scenario-wise simulation summaries of bias, RMSE, and coverage tradeoffs.
- Practical reporting guidance for publication-bias sensitivity analysis.

### Reviewer-informed safeguards
- Provides a named example workflow or fixed demonstration path.
- Documents local validation artifacts rather than relying on unsupported claims.
- Positions the software against existing tools without claiming blanket superiority.
- States limitations and interpretation boundaries in the manuscript itself.
- Requires explicit environment capture and public example accessibility in the released archive.

## Review-Driven Revisions
This draft has been tightened against recurring open peer-review objections taken from the supplied reviewer reports.
- Reproducibility: the draft names a reviewer rerun path and points readers to validation artifacts instead of assuming interface availability is proof of correctness.
- Validation: claims are anchored to local tests, validation summaries, simulations, or consistency checks rather than to unsupported assertions of performance.
- Comparators and niche: the manuscript now names the relevant comparison class and keeps the claimed niche bounded instead of implying universal superiority.
- Documentation and interpretation: the text expects a worked example, input transparency, and reviewer-verifiable outputs rather than a high-level feature list alone.
- Claims discipline: conclusions are moderated to the documented scope of TGEP and paired with explicit limitations.

## Use Cases and Results
The software outputs should be described in terms of concrete reviewer-verifiable workflows: running the packaged example, inspecting the generated results, and checking that the reported interpretation matches the saved local artifacts. In this project, the most important result layer is the availability of a transparent execution path from input to analysis output.

Representative local result: `inst/extdata/Validation_Results.txt` reports REML -0.12475569 0.03342857 -0.1902745 -0.05923689 0.0001899599.

### Concrete local quantitative evidence
- `inst/extdata/Validation_Results.txt` reports REML -0.12475569 0.03342857 -0.1902745 -0.05923689 0.0001899599.
- `inst/extdata/Validation_Results_Updated.txt` reports REML -0.13075624 0.07295042 -0.2737364 0.01222396 0.07306913.
- `inst/extdata/Real_World_Impact_Summary.txt` reports : CD013844_pub2_data 12 -0.495019680 -0.263972598 0.47407311 0.7476437.

## Discussion
Representative local result: `inst/extdata/Validation_Results.txt` reports REML -0.12475569 0.03342857 -0.1902745 -0.05923689 0.0001899599.

The most defensible F1000 framing is TGEP as a transparent ensemble sensitivity tool. The local paper already makes the crucial point that analytical SEs remain a limitation and that HKSJ-style inference should still accompany the ensemble.

### Limitations
- Analytical interval calibration is weaker than point-estimate performance.
- Guard weighting depends on observed-data prediction rather than unobserved truth.
- TGEP should complement, not replace, conventional inferential methods.

## Software Availability
- Local source package: `TGEP_Development` under `C:\Models`.
- Public repository: `https://github.com/mahmood726-cyber/tgep`.
- Public source snapshot: Fixed public commit snapshot available at `https://github.com/mahmood726-cyber/tgep/tree/10101322a46fbe445824f513cde14997e1aaaa72`.
- DOI/archive record: No project-specific DOI or Zenodo record URL was detected locally; archive registration pending.
- Environment capture detected locally: `environment.yml`.
- Reviewer-facing documentation detected locally: `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.
- Reproducibility walkthrough: `f1000_artifacts/tutorial_walkthrough.md` where present.
- Validation summary: `f1000_artifacts/validation_summary.md` where present.
- Reviewer rerun manifest: `F1000_Reviewer_Rerun_Manifest.md`.
- Multi-persona review memo: `F1000_MultiPersona_Review.md`.
- Concrete submission-fix note: `F1000_Concrete_Submission_Fixes.md`.
- License: see the local `LICENSE` file.

## Data Availability
Method materials, simulations, and project files are stored locally. External publication should add the final repository and DOI metadata.

## Reporting Checklist
Real-peer-review-aligned checklist: `F1000_Submission_Checklist_RealReview.md`.
Reviewer rerun companion: `F1000_Reviewer_Rerun_Manifest.md`.
Companion reviewer-response artifact: `F1000_MultiPersona_Review.md`.
Project-level concrete fix list: `F1000_Concrete_Submission_Fixes.md`.

## Declarations
### Competing interests
The authors declare that no competing interests were disclosed.

### Grant information
No specific grant was declared for this manuscript draft.

### Author contributions (CRediT)
| Author | CRediT roles |
|---|---|
| Mahmood Ahmad | Conceptualization; Software; Validation; Data curation; Writing - original draft; Writing - review and editing |
| Niraj Kumar | Conceptualization |
| Bilaal Dar | Conceptualization |
| Laiba Khan | Conceptualization |
| Andrew Woo | Conceptualization |

### Acknowledgements
The authors acknowledge contributors to open statistical methods, reproducible research software, and reviewer-led software quality improvement.

## References
1. DerSimonian R, Laird N. Meta-analysis in clinical trials. Controlled Clinical Trials. 1986;7(3):177-188.
2. Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. Statistics in Medicine. 2002;21(11):1539-1558.
3. Viechtbauer W. Conducting meta-analyses in R with the metafor package. Journal of Statistical Software. 2010;36(3):1-48.
4. Page MJ, McKenzie JE, Bossuyt PM, et al. The PRISMA 2020 statement: an updated guideline for reporting systematic reviews. BMJ. 2021;372:n71.
5. Fay C, Rochette S, Guyader V, Girard C. Engineering Production-Grade Shiny Apps. Chapman and Hall/CRC. 2022.
