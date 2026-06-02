# bigmath

[![Go Reference](https://pkg.go.dev/badge/github.com/db47h/bigmath.svg)](https://pkg.go.dev/github.com/db47h/bigmath)
[![Go Report Card](https://goreportcard.com/badge/github.com/db47h/bigmath)](https://goreportcard.com/report/github.com/db47h/bigmath)

**bigmath** extends Go's [`big.Float`](https://pkg.go.dev/math/big#Float) with a
complete set of elementary, trigonometric, hyperbolic, power, rounding, and
special functions, plus a `Complex` type with full arithmetic and transcendental
support. It follows the same rounding and precision semantics as `math/big`.

```go
z := new(bigmath.Float)
z.SetPrec(128)
z.Sin(z)  // same rounding/precision contract as big.Float
```

## Features

### Float

| Category | Functions |
|----------|-----------|
| **Elementary** | `Exp`, `Log`, `Log2`, `Log10`, `Pow`, `Cbrt`, `Hypot`, `FMA`, `FMS` |
| **Trigonometric** | `Sin`, `Cos`, `Sincos`, `Tan`, `Cot`, `Sec`, `Csc` |
| **Inverse Trig** | `Asin`, `Acos`, `Atan`, `Atan2` |
| **Hyperbolic** | `Sinh`, `Cosh`, `Tanh`, `SinhCosh`, `Coth`, `Sech`, `Csch` |
| **Inverse Hyperbolic** | `Asinh`, `Acosh`, `Atanh` |
| **Rounding** | `Floor`, `Ceil` |
| **Special** | `Pi` (constant), `Inv` (multiplicative inverse) |
| **Constants** | π, √2, ln 2, ln 10 (thread-safe cache) |

### Complex

`Complex` stores its real and imaginary parts as `Float` values.

| Category | Methods |
|----------|---------|
| **Arithmetic** | `Add`, `Sub`, `Mul`, `Quo`, `Neg`, `Conj`, `Inv` |
| **Transcendental** | `Exp`, `Log`, `Pow`, `Sqrt` |
| **Trigonometric** | `Sin`, `Cos`, `Tan`, `Cot`, `Sec`, `Csc` |
| **Inverse Trig** | `Asin`, `Acos`, `Atan` |
| **Hyperbolic** | `Sinh`, `Cosh`, `Tanh`, `Coth`, `Sech`, `Csch` |
| **Inverse Hyperbolic** | `Asinh`, `Acosh`, `Atanh` |
| **Utilities** | `Abs`, `Arg`, `Copy`, `Set`, `Equals`, `IsReal`, `IsZero`, `Prec`, `SetPrec`, `String`, `Format` |

`Complex` supports full aliasing (`z == x` is safe in all operations).

## Usage

```go
package main

import (
    "fmt"
    "github.com/db47h/bigmath"
)

func main() {
    z := new(bigmath.Float)
    z.SetPrec(256)

    // π
    z.Pi()
    fmt.Println(z.Text('x', -1))

    // sin(π/2) — π·2⁻¹ after z.Pi() set z = π
    z.Sin(new(bigmath.Float).SetMantExp(z, -1))
    fmt.Println(z.Text('x', -1))

    // Floor
    f := new(bigmath.Float)
    f.SetPrec(128)
    f.Parse("3.141592653589793", 0)
    f.Floor(f) // → 3

    // Complex
    a := &bigmath.Complex{Real: *bigmath.NewFloat(1), Imag: *bigmath.NewFloat(2)}
    b := &bigmath.Complex{Real: *bigmath.NewFloat(3), Imag: *bigmath.NewFloat(4)}
    c := new(bigmath.Complex)
    c.Mul(a, b) // → (-5+10i)
    c.Exp(a)    // → e¹(cos2 + i·sin2)
    fmt.Println(c.String())
}
```

## Precision & Rounding

All functions follow the same contract as [`big.Float`](https://pkg.go.dev/math/big#Float):

- **Precision propagation**: `z.Prec()` determines the result precision. If `z.Prec() == 0`, it inherits `x.Prec()` (or the maximum operand precision).
- **Rounding**: Results are rounded to `z.Prec()` bits using `z`'s rounding mode. Functions never round internally to a lower precision than the receiver.
- **Guard bits**: Each function uses a profiled guard strategy (typically `+2*_W` for Taylor loops, dynamic precision for argument reduction) to ensure 0.5-ULP accuracy.
- **Aliasing**: All functions support `z == x`. Internal allocations avoid hidden copies when possible — see `CONTEXT.md` for details.

## Status

- **Float**: All listed functions implemented and tested.
- **Complex**: All listed functions implemented and tested — branch cuts documented per ISO C standard.
- **Test coverage**: Golden-comparison tests against MPFR/gmpy2 at 128-bit, plus ULP stress tests at precisions 64–1024, plus identity tests (sin²+cos²=1, etc.) and cross-validation against `math/cmplx`.
- **Upcoming**: `Acot`, `Asec`, `Acsc`, `Acoth`, `Asech`, `Acsch`, `Erf`, `Erfc`, `Gamma`, `Lgamma`, engineering-notation string conversion.

## Performance

The package prioritises correctness and robustness over peak performance. Key
design choices:

- **True FMA semantics**: multiply-add rounds once (not twice), matching IEEE 754.
- **Hybrid argument reduction** (trig): a Payne-Hanek–inspired precision-scaled quotient estimation followed by a Ziv-style dynamic-precision remainder loop — accurate for arbitrarily large inputs without a precomputed `2/π` word table.
- **Thread-safe constant cache**: constants are computed once and reused across calls.
- **No intra-loop allocation**: Taylor/series loops recycle temps via pointer swapping.

## Requirements

- **Runtime**: Go 1.22+ — the library source is compatible with Go 1.22 and later.
- **`go generate`**: Go 1.26+ — the code generator in `internal/gen/` uses `go/types` APIs
  added in Go 1.26. The generated output (`float_gen.go`) is checked in, so end users
  at any Go >= 1.22 are unaffected.

## License

MIT — see [`LICENSE`](LICENSE) or the SPDX header in individual source files.
