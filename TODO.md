# Workspace Plan

## Status

- **`Complex.Log`**: Uses real `Log` + `Arg` (Atan2). Branch cut and special cases documented.
- **`Complex.Exp`**: Fixed — uses big.Float `Sincos` (full precision). `"math"` import removed from complex.go.
- **`trig.go`**: `Sin`, `Cos`, `Sincos` implemented. Taylor series on [0, π/2] with reduction modulo 2π. Tested via `sin_cos` in test data pipeline (decomposed into `sin`/`cos` entries).
- **`trig.go` review needed**: `workPrec = z.Prec() + 2*_W` is a simplified guard-bit heuristic. Should be revised to follow the proportional/profiled approach used in atan.go (`atanExtraBits`) or exp.go (`prec * 0.15` + fixed guard) for better precision-to-performance trade-off across all precision ranges.

## Phase 2: Expand Complex Public API

### Trigonometric (via Euler / exp + trig identities)

| Method | Formula |
|--------|---------|
| `(z *Complex) Sin(x *Complex)` | sin(a+bi) = sin(a)cosh(b) + i·cos(a)sinh(b) |
| `(z *Complex) Cos(x *Complex)` | cos(a+bi) = cos(a)cosh(b) − i·sin(a)sinh(b) |
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
