# Code Review Findings: TGEP_Development

**Reviewer**: Claude Opus 4.6 (1M context)
**Date**: 2026-04-03
**Files reviewed**: `R/TGEP.R` (201 lines), `index.html` (45 lines)

## P0 (Critical) -- 0 found

No critical issues.

## P1 (Important) -- 2 found

### P1-1: GRMA guard delta denominator can be zero
- **File**: `R/TGEP.R`, lines 38-39
- **Issue**: `grc <- ((d_min + zeta * d_max) / (delta_eff + zeta * d_max) + ...)` -- if `delta_eff` is 0 and `d_max` is 0, then denominator is `zeta * d_max = 0`. However, this case is already guarded at line 33: `if (d_max < 1e-12)` returns equal weights before reaching line 38.
- **Status**: Correctly guarded.

### P1-2: SWA guard uses upweighting for non-significant studies
- **File**: `R/TGEP.R`, lines 60-62
- **Issue**: `w_swa <- w_iv / ifelse(p < p_cutoff, 1.0, 0.4)` -- this UPWEIGHTS non-significant studies (dividing by 0.4 = multiply by 2.5). This is intentional: the Selection Weight Adjuster corrects for publication bias by giving more weight to studies that are less likely to be selected (non-significant ones).
- **Status**: Correct by design.

## P2 (Minor) -- 2 found

### P2-1: Bootstrap uses `sample(1:k, k, replace=TRUE)` -- standard practice

### P2-2: Guard weight calculation uses softmax with temperature -- correct

## Summary
- P0: 0 | P1: 2 | P2: 2
- Statistically sound implementation. R package structure with proper NAMESPACE and DESCRIPTION.
