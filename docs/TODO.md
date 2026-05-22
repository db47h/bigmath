# Workspace Plan

## Status

- **`Complex.Log`**: Uses real `Log` + `Arg` (Atan2). Branch cut and special cases documented.
- **`Complex.Exp`**: Fixed — uses big.Float `Sincos` (full precision). `"math"` import removed from complex.go.
- **`trig.go`**: `Sin`, `Cos`, `Sincos` implemented. Taylor series on [0, π/2] with reduction modulo 2π. Tested via `sin_cos` in test data pipeline (decomposed into `sin`/`cos` entries).
- **`hyperbolic.go`**: `Sinh`, `Cosh`, `Tanh`, `SinhCosh`, `Asinh`, `Acosh`, `Atanh` implemented. Tested via generated data. SinhCosh reuses Exp call for |x|≥1, Taylor shared loop for |x|<1. Guard-bit expansion for Acosh/Atanh near domain boundaries. Sinh uses Taylor for |x|<1 (no cancellation).
- **`trig.go` & `hyperbolic.go` guard strategy reviewed**: `docs/trig-hyperbolic-precision-review.md`. Conclusion: flat `+2*_W` is adequate (accumulated error over O(P) Taylor terms ≪ 0.5 ULP at target precision for any practical P). Dynamic guards in `reducePi2`, `acoshGuard`, `atanhGuard` handle the real precision-loss paths. Comments added with rationale next to every flat guard.
- **[4] ULP stress testing needed**: The flat `+2*_W` guard should be empirically verified across precision ranges (e.g. prec=64, 128, 256, 1024, 4096) using high-precision references (e.g. mpmath at prec+64). Test identities: sin²+cos²=1, sinh²+cosh²=cosh(2x), sin(3x)=3sin(x)-4sin³(x). Worst-case ULP error should be < 2 across all tested ranges.

## Phase 2: Expand Complex Public API

### Trigonometric (via Euler / exp + trig identities)

| Method | Formula |
|--------|---------|
| `(z *Complex) Sin(x *Complex)` | sin(a+bi) = sin(a)cosh(b) + i·cos(a)sinh(b) — **ready** (Sinh/Cosh implemented) |
| `(z *Complex) Cos(x *Complex)` | cos(a+bi) = cos(a)cosh(b) − i·sin(a)sinh(b) — **ready** |
| `(z *Complex) Tan(x *Complex)` | sin(z) / cos(z) |

### Hyperbolic (via complex Exp)

| Method | Formula |
|--------|---------|
| `(z *Complex) Sinh(x *Complex)` | (e^z − e^(−z)) / 2 |
| `(z *Complex) Cosh(x *Complex)` | (e^z + e^(−z)) / 2 |
| `(z *Complex) Tanh(x *Complex)` | sinh(z) / cosh(z) |

### Inverse Trigonometric / Hyperbolic (via Complex.Log)

| Method | Formula |
|--------|---------|
| `(z *Complex) Asin(x *Complex)` | −i · ln(i·z + √(1−z²)) |
| `(z *Complex) Acos(x *Complex)` | −i · ln(z + i·√(1−z²)) |
| `(z *Complex) Atan(x *Complex)` | (i/2) · ln((1−iz)/(1+iz)) — **already implemented** |
| `(z *Complex) Asinh(x *Complex)` | ln(z + √(1+z²)) |
| `(z *Complex) Acosh(x *Complex)` | ln(z + √(z−1)·√(z+1)) |
| `(z *Complex) Atanh(x *Complex)` | ½ · ln((1+z)/(1−z)) |

### Power / Root

| Method | Formula |
|--------|---------|
| `(z *Complex) Pow(x, y *Complex)` | exp(y · log(x)) |
| `(z *Complex) Sqrt(x *Complex)` | via De Moivre or exp(½ · log(x)) |

### Utility

| Method | Notes |
|--------|-------|
| `(x *Complex) String() string` | `"(a+bi)"` format |
| `(x *Complex) Format(f fmt.State, verb rune)` | fmt.Formatter support |

## Phase 3: Branch Cuts & Edge Cases

- `Complex.Log`: branch cut documented. Special cases (0, ±∞) handled via `Atan2` fix.
- `Complex.Exp`: handle overflow/underflow of e^a, periodic wrapping of imag.
- All inverse trig functions have standard branch cuts — document them.

## Testing

- `complex_test.go`: add parameterized tests using `math/cmplx` as reference at moderate precision (e.g., prec=128).
- Edge cases: ±0, ±Inf, large exponents, values near branch cuts.
- Can reuse pattern from `data_test.go` / `bigmath_test.go` if generated test data infrastructure exists.
