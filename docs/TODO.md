# Workspace Plan — Remaining Work

> Completed items moved to [`docs/TODO-done.md`](TODO-done.md).
> This file tracks only what's still pending.

## Missing Public API Functions

Functions ordered by priority (highest first). All items apply to both Float and Complex unless noted.

### P3 — Special Functions

Not in the original package description but expected by the intended audience (arbitrary-precision math consumers).

| # | Function | Float | Complex | Notes |
|---|----------|-------|---------|-------|
| 15 | **`Erf`** | ❌ | ❌ | Gauss error function — series/asymptotic expansion |
| 16 | **`Erfc`** | ❌ | ❌ | Complementary error function: `1 - Erf(x)` (with numerical care for large x) |
| 17 | **`Gamma`** | ❌ | ❌ | Γ(x) — Stirling/Lanczos approximation |
| 18 | **`Lgamma`** | ❌ | ❌ | Log-Gamma: `Log(Gamma(x))` |

### P5 — Test Robustness

| # | Task | Effort |
|---|------|--------|
| 20 | Edge-case tests for complex functions (±0, ±Inf, large exponents, branch cuts) | Medium |
