# Workspace Plan — Remaining Work

> Completed items moved to [`docs/TODO-done.md`](TODO-done.md).
> This file tracks only what's still pending.

## Remaining Gaps

### Edge-Case Tests for Complex Functions

- TestComplexData/acosh[[0.5_0]] : at prec >= 256, result should have real part = 0, but we have residual bits (1e-319).
- **Inputs**: ±0, ±Inf for all complex functions
- **Large exponents**: e.g., `Sin(1e20+0i)`
- **Values near branch cuts**: e.g., `Log(-1+εi)`, `Asin(2+0i)`, `Atan(0+1.001i)`

### Real Tan — Large-Input Precision

`Tan` at huge inputs loses precision through the `reducePi2` → `Quo(sin, cos)` path. The flat `+_W` guard may be insufficient for large `x` where `reducePi2` itself adds dynamic guards but `Tan` doesn't expand its own guard to match. No reported failures, but worth profiling.

### Documentation

No remaining gaps.

### Tests

- Unify panic tests and regular tests in data_test.go

## Priority

| Task | Effort | Impact |
|------|--------|--------|
| Edge-case complex tests | Medium | Test robustness |
| Real Tan guard profiling | Small (research only) | Future-proofing |
