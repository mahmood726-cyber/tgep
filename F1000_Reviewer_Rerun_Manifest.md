# TGEP: reviewer rerun manifest

This manifest is the shortest reviewer-facing rerun path for the local software package. It lists the files that should be sufficient to recreate one worked example, inspect saved outputs, and verify that the manuscript claims remain bounded to what the repository actually demonstrates.

## Reviewer Entry Points
- Project directory: `C:\Models\TGEP_Development`.
- Preferred documentation start points: `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.
- Detected public repository root: `https://github.com/mahmood726-cyber/tgep`.
- Detected public source snapshot: Fixed public commit snapshot available at `https://github.com/mahmood726-cyber/tgep/tree/10101322a46fbe445824f513cde14997e1aaaa72`.
- Detected public archive record: No project-specific DOI or Zenodo record URL was detected locally; archive registration pending.
- Environment capture files: `environment.yml`.
- Validation/test artifacts: `f1000_artifacts/validation_summary.md`, `inst/scripts/Run_Validation.R`, `test_tgep.R`.

## Worked Example Inputs
- Manuscript-named example paths: `TGEP_manuscript_PLOS_ONE.md` as the primary methods narrative; Project scripts and supporting files for simulations and reporting; `F1000_Submission_Checklist_RealReview.md` for reviewer-aligned packaging; f1000_artifacts/example_dataset.csv.
- Auto-detected sample/example files: `f1000_artifacts/example_dataset.csv`.

## Expected Outputs To Inspect
- Ensemble pooled point estimates and guard weights.
- Scenario-wise simulation summaries of bias, RMSE, and coverage tradeoffs.
- Practical reporting guidance for publication-bias sensitivity analysis.

## Minimal Reviewer Rerun Sequence
- Start with the README/tutorial files listed below and keep the manuscript paths synchronized with the public archive.
- Create the local runtime from the detected environment capture files if available: `environment.yml`.
- Run at least one named example path from the manuscript and confirm that the generated outputs match the saved validation materials.
- Quote one concrete numeric result from the local validation snippets below when preparing the final software paper.

## Local Numeric Evidence Available
- `inst/extdata/Validation_Results.txt` reports REML -0.12475569 0.03342857 -0.1902745 -0.05923689 0.0001899599.
- `inst/extdata/Validation_Results_Updated.txt` reports REML -0.13075624 0.07295042 -0.2737364 0.01222396 0.07306913.
- `inst/extdata/Real_World_Impact_Summary.txt` reports : CD013844_pub2_data 12 -0.495019680 -0.263972598 0.47407311 0.7476437.
