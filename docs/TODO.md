# Workspace Plan — Remaining Work

> Completed items moved to [`docs/TODO-done.md`](TODO-done.md).
> This file tracks only what's still pending.

## Remaining Gaps

### Edge-Case Tests for Complex Functions

- **Inputs**: ±0, ±Inf for all complex functions
- **Large exponents**: e.g., `Sin(1e20+0i)`
- **Values near branch cuts**: e.g., `Log(-1+εi)`, `Asin(2+0i)`, `Atan(0+1.001i)`

### Documentation

No remaining gaps.

## Priority

| Task | Effort | Impact |
|------|--------|--------|
| Edge-case complex tests | Medium | Test robustness |
