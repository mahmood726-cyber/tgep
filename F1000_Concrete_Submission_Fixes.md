# TGEP: concrete submission fixes

This file converts the multi-persona review into repository-side actions that should be checked before external submission of the F1000 software paper for `TGEP_Development`.

## Detectable Local State
- Documentation files detected: `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.
- Environment lock or container files detected: `environment.yml`.
- Package manifests detected: `DESCRIPTION`.
- Example data files detected: `f1000_artifacts/example_dataset.csv`.
- Validation artifacts detected: `f1000_artifacts/validation_summary.md`, `inst/scripts/Run_Validation.R`, `test_tgep.R`.
- Detected public repository root: `https://github.com/mahmood726-cyber/tgep`.
- Detected public source snapshot: Fixed public commit snapshot available at `https://github.com/mahmood726-cyber/tgep/tree/10101322a46fbe445824f513cde14997e1aaaa72`.
- Detected public archive record: No project-specific DOI or Zenodo record URL was detected locally; archive registration pending.

## High-Priority Fixes
- Check that the manuscript's named example paths exist in the public archive and can be run without repository archaeology.
- Confirm that the cited repository root (`https://github.com/mahmood726-cyber/tgep`) resolves to the same fixed public source snapshot used for submission.
- Archive the tagged release and insert the Zenodo DOI or record URL once it has been minted; no project-specific archive DOI was detected locally.
- Reconfirm the quoted benchmark or validation sentence after the final rerun so the narrative text matches the shipped artifacts.

## Numeric Evidence Available To Quote
- `inst/extdata/Validation_Results.txt` reports REML -0.12475569 0.03342857 -0.1902745 -0.05923689 0.0001899599.
- `inst/extdata/Validation_Results_Updated.txt` reports REML -0.13075624 0.07295042 -0.2737364 0.01222396 0.07306913.
- `inst/extdata/Real_World_Impact_Summary.txt` reports : CD013844_pub2_data 12 -0.495019680 -0.263972598 0.47407311 0.7476437.

## Manuscript Files To Keep In Sync
- `F1000_Software_Tool_Article.md`
- `F1000_Reviewer_Rerun_Manifest.md`
- `F1000_MultiPersona_Review.md`
- `F1000_Submission_Checklist_RealReview.md` where present
- README/tutorial files and the public repository release metadata
