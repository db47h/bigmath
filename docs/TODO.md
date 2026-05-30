# Workspace Plan — Remaining Work

> Completed items moved to [`docs/TODO-done.md`](TODO-done.md).
> This file tracks only what's still pending.

## Missing Public API Functions

Functions ordered by priority (highest first). All items apply to both Float and Complex unless noted.

### P0 — High Impact, Low Effort

These have supporting infrastructure already in place (`ln10`, `ln2` cached) or are fundamental operations expected in any math library.

| # | Function | Float | Complex | Via |
|---|----------|-------|---------|-----|
| 3 | **`Cbrt`** | ❌ | ❌ | Float: `Pow(x, ⅓)` with Newton refinement; Complex: `Pow(x, ⅓)` |

### P1 — Trivial One-Liners (Reciprocal Trig / Hyperbolic)

Reciprocal functions are expressible as `Inv(parentFunction)` but expected in a complete API.

| # | Function | Float | Complex | Via |
|---|----------|-------|---------|-----|
| 4 | **`Sec`** | ❌ | ❌ | `Inv(Cos(x))` |
| 5 | **`Csc`** | ❌ | ❌ | `Inv(Sin(x))` |
| 6 | **`Coth`** | ❌ | ❌ | `Inv(Tanh(x))` or `Cosh(x)/Sinh(x)` |
| 7 | **`Sech`** | ❌ | ❌ | `Inv(Cosh(x))` |
| 8 | **`Csch`** | ❌ | ❌ | `Inv(Sinh(x))` |

### P2 — Inverse Reciprocal Trig / Hyperbolic

Standard identities, less commonly used.

| # | Function | Float | Complex | Via |
|---|----------|-------|---------|-----|
| 9 | **`Acot`** | ❌ | ❌ | `Atan2(one, x)` or `π/2 - Atan(x)` |
| 10 | **`Asec`** | ❌ | ❌ | `Acos(Inv(x))` |
| 11 | **`Acsc`** | ❌ | ❌ | `Asin(Inv(x))` |
| 12 | **`Acoth`** | ❌ | ❌ | `Atanh(Inv(x))` (domain \|x\| > 1) |
| 13 | **`Asech`** | ❌ | ❌ | Log identity (domain (0, 1]) |
| 14 | **`Acsch`** | ❌ | ❌ | Log identity (domain x ≠ 0) |

### P3 — Special Functions

Not in the original package description but expected by the intended audience (arbitrary-precision math consumers).

| # | Function | Float | Complex | Notes |
|---|----------|-------|---------|-------|
| 15 | **`Erf`** | ❌ | ❌ | Gauss error function — series/asymptotic expansion |
| 16 | **`Erfc`** | ❌ | ❌ | Complementary error function: `1 - Erf(x)` (with numerical care for large x) |
| 17 | **`Gamma`** | ❌ | ❌ | Γ(x) — Stirling/Lanczos approximation |
| 18 | **`Lgamma`** | ❌ | ❌ | Log-Gamma: `Log(Gamma(x))` |

### P4 — Test Robustness

| # | Task | Effort |
|---|------|--------|
| 19 | Edge-case tests for complex functions (±0, ±Inf, large exponents, branch cuts) | Medium |
