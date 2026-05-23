# Completed Work

> Everything listed here has been implemented, tested, and is passing.
> Moved from `docs/TODO.md` on 2026-05-24.

## Core Functions

### Real Trigonometry (`trig.go`)
- **`Sin`**, **`Cos`**, **`Sincos`** — Taylor series on [0, π/2] with reduction modulo 2π.
  - Tested via `sin`/`cos` in test data pipeline (decomposed from `sin_cos`).
- **`Sin`/`Cos`/`Sincos` dispatch fixed** — now select `sinCore`/`cosCore` based on quadrant (previously always used the same core, only correct for quadrant 0).

### Real Hyperbolic (`hyperbolic.go`)
- **`Sinh`**, **`Cosh`**, **`Tanh`**, **`SinhCosh`**, **`Asinh`**, **`Acosh`**, **`Atanh`** implemented.
- SinhCosh reuses `Exp` call for |x|≥1, Taylor shared loop for |x|<1.
- Sinh uses Taylor for |x|<1 (no cancellation).
- Guard-bit expansion for Acosh/Atanh near domain boundaries.
- All special cases (±0, ±Inf) documented.

### Real Inverse Trig (`asin_acos.go`, `atan.go`)
- **`Asin`**, **`Acos`** — guard-bit aware (`asinGuard` near x=±1). All special cases documented.
- **`Atan`**, **`Atan2`** — all special cases documented.

### Real Exp / Log / Pow (`exp.go`, `log.go`, `pow.go`)
- **`Exp`** — Taylor series with argument reduction (Brent-Zimmermann method). Overflow guard at exp > 31. All special cases (+0, -Inf, +Inf) documented.
- **`Log`** — uses `computeLn` (artanh series), handles ±0, ±Inf. Panics on negative inputs.
- **`Pow`** — `x^y = exp(y·log(x))`. Overflow/underflow detection for huge exponents. Special cases documented.

### Payne-Hanek Argument Reduction
- Replaced old `reducePi2` (which had a precision mismatch comparing x at target_prec against π at workPrec) with a proper Payne-Hanek approach:
  - Multiplies x by precomputed `2/π`, extracts integer part via precision clamping + `ToZero`, computes remainder via `fma(x, n, π/2)`.
- **`halfPi`** and **`twoOverPi`** constants added (`const.go`).

### Guard-Bit Strategy Review
- **`docs/trig-hyperbolic-precision-review.md`** — comprehensive analysis.
- Conclusion: flat `+2*_W` is adequate for Taylor loops; dynamic guards in `reducePi2`, `acoshGuard`, `atanhGuard` handle the real precision-loss paths.
- Comments added next to every flat guard.

## Complex Public API (`complex.go`)

All methods implemented and supporting `z == x` aliasing:

### Arithmetic
| Method | Notes |
|--------|-------|
| **`Add`** | |
| **`Sub`** | |
| **`Mul`** | True FMA semantics for `ac-bd`, `ad+bc` |
| **`Quo`** | Full-precision denominator via FMA |
| **`Neg`** | |
| **`Conj`** | |
| **`Abs`** | Delegates to real `Hypot` |
| **`Arg`** | Delegates to real `Atan2` |

### Transcendental
| Method | Formula |
|--------|---------|
| **`Log`** | ½·ln(x²+y²) + i·Arg(x,y) — branch cut along negative real axis. Special cases (0, ±∞) documented. |
| **`Exp`** | e^a·cos(b) + i·e^a·sin(b) via real `Exp` + `Sincos` |
| **`Sin`** | sin(a)cosh(b) + i·cos(a)sinh(b) |
| **`Cos`** | cos(a)cosh(b) − i·sin(a)sinh(b) |
| **`Tan`** | sin(z) / cos(z) via Quo |
| **`Sinh`** | sinh(a)cos(b) + i·cosh(a)sin(b) |
| **`Cosh`** | cosh(a)cos(b) + i·sinh(a)sin(b) |
| **`Tanh`** | sinh(z) / cosh(z) via Quo |
| **`Asin`** | −i·ln(i·z + √(1−z²)) |
| **`Acos`** | −i·ln(z + i·√(1−z²)) |
| **`Atan`** | (i/2)·ln((1−iz)/(1+iz)) — all special cases documented |
| **`Asinh`** | ln(z + √(1+z²)) |
| **`Acosh`** | ln(z + √(z−1)·√(z+1)) |
| **`Atanh`** | ½·ln((1+z)/(1−z)) |
| **`Sqrt`** | exp(½·log(x)) |
| **`Pow`** | exp(y·log(x)) |

### Utility
| Method | Notes |
|--------|-------|
| **`String()`** | Delegates to `Format` |
| **`Format()`** | `(a+bi)` format. ∞ rendered as `±∞`. Special cases: `0 → "0"`, `0+bi → "bi"`, `a+0i → "a"`, `±1i → "±i"`. |
| **`Copy`** | Deep copy |
| **`Equals`** | Real and imag exact comparison |
| **`IsReal`** | Imaginary part is zero |
| **`IsZero`** | Both parts zero |
| **`Prec`** / **`SetPrec`** | Precision management |
| **`setPrec`** / **`setPrec2`** | Internal precision normalization |

## Testing

| Test | Coverage |
|------|----------|
| **ULP stress testing** | prec=64, 128, 256, 1024. All pass at `maxULP=0.5` (including a one-off at 10240 bits). |
| **Data-driven tests** (`TestBigMath`) | Generated via Python/gmpy2: sin, cos, tan, asin, acos, atan, atan2, sinh, cosh, tanh, asinh, acosh, atanh, exp, log, pow, π constant. |
| **Real identity tests** | `sin²+cos²=1`, `sin(3x)=3sin(x)−4sin³(x)`, `cosh²−sinh²=1`. All at 64/128/256/1024. |
| **Complex identity tests** | `sin²(z)+cos²(z)=1` at 256-bit. |
| **Complex vs math/cmplx** | All 13 complex functions tested against Go stdlib at 53-bit for input `(0.5, 0.7)`. |
| **Complex Format** | Full 5×5 grid of reals×imags spanning {−∞, neg, 0, pos, +∞} plus 20+ special cases (`0`, `i`, `-i`, `1+i`, `1-i`, `+∞+3i`, etc.) and `%.2f` formatting. |
| **Complex aliasing** | `Sin(z)`, `Mul(z,z)`, `Quo(z,z)` — all pass with `z==x`. |
| **Error handling** | `Log(-1)` etc. panic with `ErrNaN`. |
| **Test data generation fixed** | `REF_MARGIN` changed from 64 to 0 so `x` is parsed at the same precision in Go and gmpy2. |
| **Identity tests with guard bits** | All rewritten to use `prec+64` or `prec+128` on polynomial evaluations. |
| **1024-bit precision added** | to test suite. |

## Documentation

- **`docs/TODO.md`** — project plan
- **`docs/Payne-Hanek.md`** — Payne-Hanek argument reduction technique write-up.
- **`docs/trig-hyperbolic-precision-review.md`** — guard-bit precision strategy review with analysis and conclusion.
- **Code doc comments** on all exported functions with special cases listed.
- Branch cut documented on `Complex.Log`.

## Configuration / Infrastructure

- `go.mod` / `go.sum` — module setup.
- Test data pipeline — Python generator (`testdata/gen_go_tests.py`, `testdata/data.txt`).
- Benchmark test (`bench_test.go`).
- Constants (`const.go`) — `pi`, `ln2`, `ln10`, `sqrt2`, `halfPi`, `twoOverPi` via thread-safe caching.
