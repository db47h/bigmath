# Workspace Plan

## Status

- **`Complex.Log`**: Uses real `Log` + `Arg` (Atan2). Branch cut and special cases documented.
- **`Complex.Exp`**: Fixed — uses big.Float `Sincos` (full precision). `"math"` import removed from complex.go.
- **`trig.go`**: `Sin`, `Cos`, `Sincos` implemented. Taylor series on [0, π/2] with reduction modulo 2π. Tested via `sin_cos` in test data pipeline (decomposed into `sin`/`cos` entries).
- **`hyperbolic.go`**: `Sinh`, `Cosh`, `Tanh`, `SinhCosh`, `Asinh`, `Acosh`, `Atanh` implemented. Tested via generated data. SinhCosh reuses Exp call for |x|≥1, Taylor shared loop for |x|<1. Guard-bit expansion for Acosh/Atanh near domain boundaries. Sinh uses Taylor for |x|<1 (no cancellation).
- **`trig.go` & `hyperbolic.go` guard strategy reviewed**: `docs/trig-hyperbolic-precision-review.md`. Conclusion: flat `+2*_W` is adequate (accumulated error over O(P) Taylor terms ≪ 0.5 ULP at target precision for any practical P). Dynamic guards in `reducePi2`, `acoshGuard`, `atanhGuard` handle the real precision-loss paths. Comments added with rationale next to every flat guard.
- **[4] ULP stress testing — DONE**: Verified across prec=64, 128, 256, 1024. All tests pass at `maxULP=0.5`.
- **Test data generation fixed**: `REF_MARGIN` changed from 64 to 0 so that `x` is parsed at the same precision in both Go and gmpy2.
- **`reducePi2` guard calculation fixed**: Simplified to `workPrec = prec + _W + xExp` (when `xExp>0`), avoiding under-guarding for large inputs.
- **All identity tests rewritten** with guard bits (`prec+64` or `prec+128`) on polynomial evaluations to prevent precision loss in multiplications.
- **1024-bit precision added** to test suite.

## Known Bugs & Discrepancies

### 🐛 Cos(3π/2) 256-bit: Go vs gmpy2 differ after bit 63

Go's `Cos(3π/2)` at 256-bit precision disagrees with gmpy2 (Python) when both
receive the **same x** at the **same precision**. The first 63 bits match;
everything after bit 63 (0-indexed from MSB) diverges — 104 differing bits in
total. At 64-bit and 128-bit precision, Go matches gmpy2 exactly.

**How to reproduce:**

```go
// In Go:
prec := uint(256)
x := new(big.Float).SetPrec(prec)
x.Parse("0x4.b65f1fccc8748d3c9ca64f450528ace6f60dd4333e6ecab80c4677e56275a2dp+0", 0)
z := new(big.Float).SetPrec(prec)
bigmath.Cos(z, x)
fmt.Println(z.Text('p', -1))
// Result: 0x.8610f349aab1f8b2f7ca8cd9e69d218d9a03298cc33a71b019af9b932d5b1a41p-254
```

```python
# In Python (same x, same precision):
import gmpy2
gmpy2.get_context().precision = 256
x = gmpy2.mpfr("0x4.b65f1fccc8748d3c9ca64f450528ace6f60dd4333e6ecab80c4677e56275a2dp+0")
print(format(gmpy2.cos(x), 'a'))
# Result: 0x2.1843cd26aac7e2cc628165c930a26d5cdefdc16c51c586b420b8bf6f7015725cp-256
```

These are the same value expressed in different hex-float conventions
(Go's `0x.xxx...pE` vs gmpy2's `0xM.FFF...pE`). Converting:

```
Go:    0x.8610f349aab1f8b2f7ca8cd9e69d218d9a03298cc33a71b019af9b932d5b1a41p-254
     = 0x21843cd26aac7e2cbdf2a33679a748636680ca6330ce9c6c066be6e4cb56c6904 * 2^(-512)

gmpy2: 0x2.1843cd26aac7e2cc628165c930a26d5cdefdc16c51c586b420b8bf6f7015725cp-256
     = 0x21843cd26aac7e2cc628165c930a26d5cdefdc16c51c586b420b8bf6f7015725c * 2^(-512)
```

First differing bit at position 63 (0-indexed from MSB of the 257-bit unsigned
integer). The ULP error relative to gmpy2 is ~2.8e-20 (nearly zero), but the
result only carries ~63 correct bits out of 256.

**Hint:** Go's internal `workPrec` for Cos at 256-bit is `prec + _W = 256 + 64 = 320`,
and `reducePi2` adds `xExp=2` → workPrec=322. `cosCore` uses `prec+2*_W = 384`.
These should be sufficient. Check whether the issue is in argument reduction
(`reducePi2`), Taylor series loop termination, or a subtle precision leak in
the core computation at high precisions.

### 🐛 TestSinTripleAngle: two issues

✅ **Fixed** — see `ulp_test.go`

1. **RHS computed without guard bits:** The identity test computes
   `3·sin(x) − 4·sin³(x)` at `prec` throughout. The `sin³` term loses
   precision through three consecutive multiplications at `prec`, and the
   final linear combination amplifies that loss. For small `sin(x)` (e.g.
   near zeros), `sin³(x)` underflows to zero or carries far fewer than `prec`
   bits, making the RHS significantly less accurate than the LHS (`sin3x`
   from `bigmath.Sin(3x)`). This can cause false-positive failures or,
   worse, mask real accuracy problems by comparing two equally wrong values.

2. **Ignores reference data:** The JSON test data includes a `sin_3x_ref`
   field computed by gmpy2 at `ref_prec`. The test never compares `sin3x`
   against this reference — it only checks the LHS-vs-RHS identity. This
   means a bug where both `Sin(3x)` and `Sin(x)` are wrong in a correlated
   way (e.g. a sign error in both) would pass undetected.

   **Fix applied:**
   - Added primary ULP comparison against `pt.Sin3xRef` (like `TestULPErrorDirect`
     does for `sin_ref`/`cos_ref`), with `maxULP = 0.5`.
   - Kept the identity check as a secondary sanity test with guard bits
     (`prec + 64`) on the RHS computation to prevent precision loss in the
     polynomial evaluation. Identity errors >2 ULPs are logged but do not
     cause test failure.

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
