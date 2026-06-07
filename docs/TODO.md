# Workspace Plan — Remaining Work

> Completed items moved to [`docs/TODO-done.md`](TODO-done.md).
> This file tracks only what's still pending.

## Missing tests

Float.Acosh(-Inf)
Float.Acoth(-Inf)
Test Float to string: test integer path for fexp in initHuge.

## Aliasing doc review

docs/aliasing.md details which functions of the big.Float API do not cause hidden
allocations when called with arguments aliasing the receiver. We need to extend
that to the full bigmath.Float API. Aliasing is fully supported but:

- some functions do not have any specific handling for this (because it's not needed)
- others have defensive aliasing handling (copy x, just in case...)

This should be reviewed and documented. Functions with defensive handling of
aliasing may need updating to make copies of x only if needed.

- [ ] Finish updating the aliasing doc.

## Gamma/Lgamma review

- [ ] Review Complex.Gamma, Complex.Lgamma

## Dynamic precision guards

- [ ] Review the dynamic precision guards (they&#39;re overkill) and update edge cases for these

## Test gen warnings

- [ ] Fix warnings in test gen about Complex.Gamma tests flagged as `!panic` but not producing NaN

## Missing Public API Functions

Functions ordered by priority (highest first). All items apply to both Float and Complex unless noted.

### P3 — Special Functions

Not in the original package description but expected by the intended audience (arbitrary-precision math consumers).

| # | Function | Float | Complex | Notes |
|---|----------|-------|---------|-------|
| 15 | **`Erf`** | ❌ | ❌ | Gauss error function — series/asymptotic expansion |
| 16 | **`Erfc`** | ❌ | ❌ | Complementary error function: `1 - Erf(x)` (with numerical care for large x) |
| 17 | **`Gamma`** | ✅ | ❌ | Γ(x) — Stirling series with Bernoulli numbers |
| 18 | **`Lgamma`** | ✅ | ❌ | Log-Gamma: `(log|Γ|, sign)` with reflection formula |

### P4 — Gamma/Lgamma Deferred Optimizations

| # | Task | Effort |
|---|------|--------|
| 19 | Lanczos fast path for Gamma at ≤128-bit | Small |
| 20 | Dynamic precision optimization for high-order Stirling terms | Small |
| 21 | Rising factorial loop: use product of reciprocals instead of per-term Log | Small |

### P6 — Test Robustness

| # | Task | Effort |
|---|------|--------|
| 22 | Edge-case tests for complex functions (±0, ±Inf, large exponents, branch cuts) | Medium |
