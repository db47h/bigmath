# Workspace Plan — Remaining Work

> Completed items moved to [`docs/TODO-done.md`](TODO-done.md).
> This file tracks only what's still pending.

## Remaining Gaps

### Branch Cut Documentation on Complex Inverse Functions

`Complex.Log` documents its branch cut (negative real axis). The following complex methods lack branch cut / special-case doc comments:

| Method | File line | Missing |
|--------|-----------|---------|
| `Complex.Asin` | 281 | Branch cut, special cases |
| `Complex.Acos` | 308 | Branch cut, special cases |
| `Complex.Atan` | 146 | Special cases (branch cut: imaginary axis beyond ±i) |
| `Complex.Asinh` | 336 | Branch cut, special cases |
| `Complex.Acosh` | 354 | Branch cut, special cases |
| `Complex.Atanh` | 372 | Branch cut, special cases |

Standard branch cuts follow `math/cmplx` conventions:
- `Asin`: two cuts on the real axis, outside [-1, +1]
- `Acos`: two cuts on the real axis, outside [-1, +1]
- `Atan`: two cuts on the imaginary axis, outside [-i, +i]
- `Asinh`: two cuts on the imaginary axis, outside [-i, +i]
- `Acosh`: a cut on the real axis, x < 1
- `Atanh`: two cuts on the real axis, outside [-1, +1]

### Edge-Case Tests for Complex Functions

Current complex test coverage (`TestComplex_AgainstStd`) tests a single input `(0.5, 0.7)` at 53-bit. Missing:

- **Inputs**: ±0, ±Inf for all complex functions
- **Large exponents**: e.g., `Sin(1e20+0i)`
- **Values near branch cuts**: e.g., `Log(-1+εi)`, `Asin(2+0i)`, `Atan(0+1.001i)`
- **Parameterized tests at higher precisions** (128-bit, 256-bit) against `math/cmplx`

### Real Tan — Large-Input Precision

`Tan` at huge inputs loses precision through the `reducePi2` → `Quo(sin, cos)` path. The flat `+_W` guard may be insufficient for large `x` where `reducePi2` itself adds dynamic guards but `Tan` doesn't expand its own guard to match. No reported failures, but worth profiling.

### Documentation

No remaining gaps.

## Priority

| Task | Effort | Impact |
|------|--------|--------|
| Branch cut doc comments | Small | Surface-level completeness |
| Edge-case complex tests | Medium | Test robustness |
| Real Tan guard profiling | Small (research only) | Future-proofing |
