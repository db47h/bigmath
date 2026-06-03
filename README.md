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
| **Inverse Trig** | `Asin`, `Acos`, `Atan`, `Atan2`, `Acot`, `Asec`, `Acsc` |
| **Hyperbolic** | `Sinh`, `Cosh`, `Tanh`, `SinhCosh`, `Coth`, `Sech`, `Csch` |
| **Inverse Hyperbolic** | `Asinh`, `Acosh`, `Atanh`, `Acoth`, `Asech`, `Acsch` |
| **Rounding** | `Floor`, `Ceil` |
| **Special** | `Pi` (constant), `Inv` (multiplicative inverse), `AbsCmp` (compare absolute values) |
| **Constants** | π, √2, ln 2, ln 10 (thread-safe cache) |
| **Formatting** | `Text`, `Format` (`fmt.Formatter`), `String`, `Append`, `AppendText`, `MarshalText` — supports `'e'`, `'E'`, `'f'`, `'g'`, `'G'`, `'b'`, `'p'`, `'x'`, plus **engineering notation** `'n'`/`'N'` |

### Complex

`Complex` stores its real and imaginary parts as `Float` values.

| Category | Methods |
|----------|---------|
| **Arithmetic** | `Add`, `Sub`, `Mul`, `Quo`, `Neg`, `Conj`, `Inv` |
| **Transcendental** | `Exp`, `Log`, `Pow`, `Sqrt` |
| **Trigonometric** | `Sin`, `Cos`, `Tan`, `Cot`, `Sec`, `Csc` |
| **Inverse Trig** | `Asin`, `Acos`, `Atan`, `Acot`, `Asec`, `Acsc` |
| **Hyperbolic** | `Sinh`, `Cosh`, `Tanh`, `Coth`, `Sech`, `Csch` |
| **Inverse Hyperbolic** | `Asinh`, `Acosh`, `Atanh`, `Acoth`, `Asech`, `Acsch` |
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

## Requirements

- **Runtime**: Go 1.22+ — the library source is compatible with Go 1.22 and later.
- **`go generate`**: Go 1.26+ — the code generator in `internal/gen/` uses `go/types` APIs
  added in Go 1.26. The generated output (`float_gen.go`) is checked in, so end users
  at any Go >= 1.22 are unaffected.

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
  
  The golden test data is checked in as JSON files under `testdata/`. After adding entries to
  [`testdata/data.txt`](testdata/data.txt) or [`testdata/cplx_data.txt`](testdata/cplx_data.txt),
  regenerate with:
  
  ```bash
  go generate ./...
  ```
  
  This runs the [`go:generate`](https://go.dev/wiki/Gen) directives in
  [`data_test.go`](data_test.go) and [`cplx_data_test.go`](cplx_data_test.go), which invoke
  [`testdata/gen_go_tests.py`](testdata/gen_go_tests.py) at 128-bit precision.
  
  **Prerequisites**: [Python 3](https://www.python.org/) with
  [**gmpy2**](https://gmpy2.readthedocs.io/) (the MPFR/MPC wrapper). Install with:
  
  ```bash
  pip install gmpy2
  ```
  
  To generate at a different precision, use the script directly:
  
  ```bash
  python testdata/gen_go_tests.py testdata/data.txt -o testdata/data_tests.json -p 256
  python testdata/gen_go_tests.py -m cplx testdata/cplx_data.txt -o testdata/cplx_data_tests.json -p 256
  ```
- **Upcoming**: `Erf`, `Erfc`, `Gamma`, `Lgamma`.

## Performance

The package prioritises correctness and robustness over peak performance. Key
design choices:

- **True FMA semantics**: multiply-add rounds once (not twice), matching IEEE 754.
- **Hybrid argument reduction** (trig): a Payne-Hanek–inspired precision-scaled quotient estimation followed by a Ziv-style dynamic-precision remainder loop — accurate for arbitrarily large inputs without a precomputed `2/π` word table.
- **Thread-safe constant cache**: constants are computed once and reused across calls.
- **No intra-loop allocation**: Taylor/series loops recycle temps via pointer swapping.

## Formatting

`Float` implements the full `fmt.Formatter` interface and supports:

| Verb | Description |
|------|-------------|
| `'e'`/`'E'` | Scientific notation (`-d.dddde±dd`) |
| `'f'`/`'F'` | Fixed-point (`-ddddd.dddd`) |
| `'g'`/`'G'` | General format (like `'e'` for large/small exponents, `'f'` otherwise) |
| `'n'`/`'N'` | **Engineering notation** — exponent is always a multiple of 3, mantissa has 1–3 leading digits. Lowercase uses `'e'`, uppercase uses `'E'`. (non-standard, specific to this package) |
| `'b'` | Decimal mantissa with binary exponent (non-standard) |
| `'p'` | Hexadecimal mantissa with binary exponent (non-standard) |
| `'x'` | Hexadecimal mantissa with decimal power-of-two exponent |

### Engineering notation examples

| Value | `Text('e', 6)` | `Text('n', 6)` | `Text('n', -1)` (shortest) |
|-------|----------------|----------------|----------------------------|
| `12345` | `1.234500e+04` | `12.345000e+03` | `12.345e+03` |
| `0.000123` | `1.230000e-04` | `123.000000e-06` | `123e-06` |
| `1e-13` | `1.000000e-13` | `100.000000e-15` | `100e-15` |
| `1e100` | `1.000000e+100` | `10.000000e+99` | `10e+99` |
| `1` | `1.000000e+00` | `1.000000e+00` | `1e+00` |
| `0` | `0.000000e+00` | `0.000000e+00` | `0e+00` |

### Workaround for Go issue [#11068](https://github.com/golang/go/issues/11068)

The standard Go `big.Float.Text()` method suffers from **O(exp × prec²) time
complexity** when formatting numbers with very large exponents (e.g. `1e1000000`),
making it unusable for exponent magnitudes above a few million. This package's
`Append()` implementation uses a custom fast path for the `'e'`, `'E'`, `'f'`,
`'g'`, `'G'`, `'n'`, and `'N'` formats that avoids the expensive decimal
conversion and scales linearly with the exponent magnitude. Formats `'b'`, `'p'`,
and `'x'` delegate to `math/big` (already fast for those forms).

```go
// bigmath.Float handles this instantly:
f := new(bigmath.Float).SetPrec(128)
f.SetString("1e10000000")
fmt.Println(f.Text('g', -1))  // fast, not O(exp × prec²)
```

## License

MIT — see [`LICENSE`](LICENSE). The big.Float forwarders and the formatting code duplicate parts of the Go standard library, covered by the [`Go LICENSE`](LICENSE-go).
