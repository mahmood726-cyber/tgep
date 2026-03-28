# TGEP: multi-persona peer review

This memo applies the recurring concerns in the supplied peer-review document to the current F1000 draft for this project (`TGEP_Development`). It distinguishes changes already made in the draft from repository-side items that still need to hold in the released repository and manuscript bundle.

## Detected Local Evidence
- Detected documentation files: `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.
- Detected environment capture or packaging files: `environment.yml`.
- Detected validation/test artifacts: `f1000_artifacts/validation_summary.md`, `inst/scripts/Run_Validation.R`, `test_tgep.R`.
- Detected browser deliverables: no HTML file detected.
- Detected public repository root: `https://github.com/mahmood726-cyber/tgep`.
- Detected public source snapshot: Fixed public commit snapshot available at `https://github.com/mahmood726-cyber/tgep/tree/10101322a46fbe445824f513cde14997e1aaaa72`.
- Detected public archive record: No project-specific DOI or Zenodo record URL was detected locally; archive registration pending.

## Reviewer Rerun Companion
- `F1000_Reviewer_Rerun_Manifest.md` consolidates the shortest reviewer-facing rerun path, named example files, environment capture, and validation checkpoints.

## Detected Quantitative Evidence
- `inst/extdata/Validation_Results.txt` reports REML -0.12475569 0.03342857 -0.1902745 -0.05923689 0.0001899599.
- `inst/extdata/Validation_Results_Updated.txt` reports REML -0.13075624 0.07295042 -0.2737364 0.01222396 0.07306913.
- `inst/extdata/Real_World_Impact_Summary.txt` reports : CD013844_pub2_data 12 -0.495019680 -0.263972598 0.47407311 0.7476437.

## Current Draft Strengths
- States the project rationale and niche explicitly: No single publication-bias-sensitive estimator behaves optimally across all meta-analytic contamination patterns. Ensemble approaches can help, but only if their components, weighting logic, and performance tradeoffs are documented clearly enough for peer review.
- Names concrete worked-example paths: `TGEP_manuscript_PLOS_ONE.md` as the primary methods narrative; Project scripts and supporting files for simulations and reporting; `F1000_Submission_Checklist_RealReview.md` for reviewer-aligned packaging.
- Points reviewers to local validation materials: Local simulation outputs summarized in the manuscript materials; Reviewer-facing discussion of coverage, weighting, and guard behavior; Submission artifacts documenting software packaging and scope boundaries.
- Moderates conclusions and lists explicit limitations for TGEP.

## Remaining High-Priority Fixes
- Keep one minimal worked example public and ensure the manuscript paths match the released files.
- Ensure README/tutorial text, software availability metadata, and public runtime instructions stay synchronized with the manuscript.
- Confirm that the cited repository root resolves to the same fixed public source snapshot used for the submission package.
- Mint and cite a Zenodo DOI or record URL for the tagged release; none was detected locally.
- Reconfirm the quoted benchmark or validation sentence after the final rerun so the narrative text stays synchronized with the shipped artifacts.

## Persona Reviews

### Reproducibility Auditor
- Review question: Looks for a frozen computational environment, a fixed example input, and an end-to-end rerun path with saved outputs.
- What the revised draft now provides: The revised draft names concrete rerun assets such as `TGEP_manuscript_PLOS_ONE.md` as the primary methods narrative; Project scripts and supporting files for simulations and reporting and ties them to validation files such as Local simulation outputs summarized in the manuscript materials; Reviewer-facing discussion of coverage, weighting, and guard behavior.
- What still needs confirmation before submission: Before submission, freeze the public runtime with `environment.yml` and keep at least one minimal example input accessible in the external archive.

### Validation and Benchmarking Statistician
- Review question: Checks whether the paper shows evidence that outputs are accurate, reproducible, and compared against known references or stress tests.
- What the revised draft now provides: The manuscript now cites concrete validation evidence including Local simulation outputs summarized in the manuscript materials; Reviewer-facing discussion of coverage, weighting, and guard behavior; Submission artifacts documenting software packaging and scope boundaries and frames conclusions as being supported by those materials rather than by interface availability alone.
- What still needs confirmation before submission: Concrete numeric evidence detected locally is now available for quotation: `inst/extdata/Validation_Results.txt` reports REML -0.12475569 0.03342857 -0.1902745 -0.05923689 0.0001899599; `inst/extdata/Validation_Results_Updated.txt` reports REML -0.13075624 0.07295042 -0.2737364 0.01222396 0.07306913.

### Methods-Rigor Reviewer
- Review question: Examines modeling assumptions, scope conditions, and whether method-specific caveats are stated instead of implied.
- What the revised draft now provides: The architecture and discussion sections now state the method scope explicitly and keep caveats visible through limitations such as Analytical interval calibration is weaker than point-estimate performance; Guard weighting depends on observed-data prediction rather than unobserved truth.
- What still needs confirmation before submission: Retain method-specific caveats in the final Results and Discussion and avoid collapsing exploratory thresholds or heuristics into universal recommendations.

### Comparator and Positioning Reviewer
- Review question: Asks what gap the tool fills relative to existing software and whether the manuscript avoids unsupported superiority claims.
- What the revised draft now provides: The introduction now positions the software against an explicit comparator class: The paper explicitly compares the ensemble with REML, HKSJ, trim-and-fill, PET-PEESE, and other bias-sensitive alternatives, while keeping claims limited to the contexts actually evaluated in the local materials.
- What still needs confirmation before submission: Keep the comparator discussion citation-backed in the final submission and avoid phrasing that implies blanket superiority over better-established tools.

### Documentation and Usability Reviewer
- Review question: Looks for a README, tutorial, worked example, input-schema clarity, and short interpretation guidance for outputs.
- What the revised draft now provides: The revised draft points readers to concrete walkthrough materials such as `TGEP_manuscript_PLOS_ONE.md` as the primary methods narrative; Project scripts and supporting files for simulations and reporting; `F1000_Submission_Checklist_RealReview.md` for reviewer-aligned packaging and spells out expected outputs in the Methods section.
- What still needs confirmation before submission: Make sure the public archive exposes a readable README/tutorial bundle: currently detected files include `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.

### Software Engineering Hygiene Reviewer
- Review question: Checks for evidence of testing, deployment hygiene, browser/runtime verification, secret handling, and removal of obvious development leftovers.
- What the revised draft now provides: The draft now foregrounds regression and validation evidence via `f1000_artifacts/validation_summary.md`, `inst/scripts/Run_Validation.R`, `test_tgep.R`, and browser-facing projects are described as self-validating where applicable.
- What still needs confirmation before submission: Before submission, remove any dead links, exposed secrets, or development-stage text from the public repo and ensure the runtime path described in the manuscript matches the shipped code.

### Claims-and-Limitations Editor
- Review question: Verifies that conclusions are bounded to what the repository actually demonstrates and that limitations are explicit.
- What the revised draft now provides: The abstract and discussion now moderate claims and pair them with explicit limitations, including Analytical interval calibration is weaker than point-estimate performance; Guard weighting depends on observed-data prediction rather than unobserved truth; TGEP should complement, not replace, conventional inferential methods.
- What still needs confirmation before submission: Keep the conclusion tied to documented functions and artifacts only; avoid adding impact claims that are not directly backed by validation, benchmarking, or user-study evidence.

### F1000 and Editorial Compliance Reviewer
- Review question: Checks for manuscript completeness, software/data availability clarity, references, and reviewer-facing support files.
- What the revised draft now provides: The revised draft is more complete structurally and now points reviewers to software availability, data availability, and reviewer-facing support files.
- What still needs confirmation before submission: Confirm repository/archive metadata, figure/export requirements, and supporting-file synchronization before release.
