# bigmath Context — Domain & Coding Conventions

## Project

`bigmath` provides arbitrary-precision mathematical functions for `[big.Float]` (Go's `math/big`). It extends `big.Float` with trigonometry, exponentiation, logarithms, hyperbolic functions, and complex-number support — all with the same rounding and precision semantics as the standard library.

See [TODO.md](./TODO.md) for the current workspace plan and status.

## Key Design Decisions

### Precision Propagation
All public functions follow the pattern `func F(z, x *big.Float) *big.Float`:
- If `z.Prec() == 0`, inherit precision from `x.Prec()`
- `z` is returned as the result
- Internal computation uses a higher `workPrec` with guard bits; result is rounded down to `z.Prec()` on return

### Fused Multiply-Add (fma / FMA in bigmath.go)
**`fma(z, x, y, t, temp)` is correct. Do not question or attempt to optimize it.**

It provides genuine FMA semantics (one rounding for `x*y + t`):

- `big.Float.Mul` always computes the full mantissa product (`O(n×m)` words)
  internally, then rounds to the target precision. The full computation
  happens regardless of what precision you set on the result.
- By sizing `temp` at `x.Prec() + y.Prec()`, we tell `Mul` to keep every bit
  of that product — no intermediate rounding.
- The addition `z.Add(temp.Mul(x, y), t)` rounds once to `z.Prec()`.
- **Result: one rounding for the entire expression.**

The sum-of-precisions temp size (`x.Prec() + y.Prec()`) is NOT wasteful:
the full product was computed internally by `Mul` anyway; the cost is only
storing it unrounded. Shrinking `temp` would introduce intermediate rounding
and break the FMA guarantee.

`Complex.Mul` follows the same principle: `bd` and `bc` use sum-of-precisions
so they feed full-precision products into `fma`, giving genuine FMA for the
complex product formula `(ac−bd) + i·(ad+bc)`.

### Aliasing Guarantee
**All `bigmath` API functions — including internal helpers — support aliasing.** If `z == x`, the function computes correctly. This mirrors `big.Float`'s own aliasing guarantee.

The reason to **avoid** passing the receiver as an argument to `big.Float` methods in internal code is **not** semantic correctness (that's always fine), but **hidden allocations**: `big.Float`'s Add, Sub, Mul, Quo, etc. must internally copy the receiver when it aliases an argument, which defeats the performance goal of recycling temps. Avoid aliasing in internal `big.Float` calls unless you know the specific operation is allocation-free when aliasing (e.g. `Neg`, `Set`, `SetMantExp`, `SetPrec`, read-only accessors like `Sign`, `Signbit`, `MantExp(nil)`, `IsInf`, `Prec`, `Cmp`).

Some `bigmath` internal helpers (like `reducePi2`) are explicitly designed to handle aliasing without hidden allocations by copying early.

### Constant Cache (`const.go`)
Mathematical constants (`π`, `ln2`, `ln10`, `√2`) are computed once per precision via a thread-safe `cache()` wrapper that returns the internal pointer directly — no copy. Callers MUST treat returned `*big.Float` as **read-only** (modifying them corrupts the cache).

Pre-allocated singletons: `zero`, `one`, `two`, `ten`, `minusOne`. These are also read-only — copy if mutation is needed.

### Guard-Bit Heuristics
Each function uses a profiled guard-bit strategy rather than a flat `+ 2*_W`:
- `exp.go` Taylor loop: `prec += uint(math.Log(float64(prec))) + 1` plus `2 * uint(squarings)` for squaring steps (see code for full logic)
- `exp.go` Arg reduction: `-k < exp` causes `prec += 2 * uint(squarings)`
- `atan.go` Atan reduction loop: `prec + 8 + uint(max(0, min(x.MantExp(nil), 2)+u))` where `u = max(int(math.Sqrt(float64(prec))) / 4, 2)`
- `atan.go` atanCore: `prec + 2*_W` (flat, core series only)
- `log.go`: `prec + _W`
- `trig.go` sinCore/cosCore/sincosCore: `prec + 2*_W` — noted with a TODO for revision toward a proportional/profiled approach
- `hyperbolic.go` sinhCore/sinhcoshCore: `prec + 2*_W`

The convention: core Taylor-series loops use a flat guard; argument-reduction layers add extra precision proportional to the reduction depth.

## Coding Style Rules

These rules MUST be followed in all new code and edits. The workflow prompt instructs every subagent to read this file.

### 1. SPDX License Header
Every `.go` file must begin with:
```go
// SPDX-License-Identifier: MIT

package bigmath
```

### 2. Avoid Aliasing in Internal `big.Float` Calls
- All bigmath API functions and helpers support aliasing correctly — the concern is **performance**, not correctness.
- In internal code, avoid calling `big.Float` methods with the receiver aliasing an argument, because `big.Float` will allocate a hidden copy internally.
- Exception: operations that are allocation-free when aliasing — `Neg`, `Set`, `SetMantExp`, `SetPrec`, `SetUint64`, `SetInt64`, and read-only methods (`Sign`, `Signbit`, `MantExp(nil)`, `IsInf`, `Prec`, `Cmp`).
- Exception: bigmath's own `reducePi2` handles aliasing without hidden allocations (copies early).
- When aliasing is unavoidable in a bigmath function, copy inputs to temps at the top and defer `z.Set(result)` until the end.

### 3. Minimize Temporary Allocations
- Allocate all scratch `*big.Float` variables at the **top** of the function (or loop preamble) with the appropriate `workPrec`.
- Use `newFloat(workPrec)` for allocation — never `new(big.Float)` directly (it leaves precision at 0, causing reallocation on every operation).
- **Never allocate inside a loop body.** Use pointer swapping (`sum, t0 = t0, sum`) to recycle temps.
- Pattern:
  ```go
  t0 := newFloat(workPrec)
  t1 := newFloat(workPrec)
  // ... use t0, t1 throughout, swap pointers instead of re-allocating
  ```

### 4. Use Project Constants, Don't Recompute
- Use `pi(prec)`, `ln2(prec)`, `ln10(prec)`, `sqrt2(prec)` from `const.go` — they're cached and precision-adaptive.
- Use `one`, `two`, `ten`, `zero`, `minusOne` for trivial constants.
- **Do not** compute π, ln2, √2 inline — always use the cache.
- When you need `π/2`, set a temp from `pi(workPrec)` then halve it (e.g. `t0.Set(pi(workPrec)).SetMantExp(t0, -1)`).

### 5. Internal Helper Signature Pattern
Internal functions take a `temp` or work-precision parameter explicitly when the caller controls allocation:
```go
func fma(z, x, y, t, temp *big.Float) *big.Float
func sinCore(z, x *big.Float) *big.Float        // allocates temps internally
```

### 6. Error Handling
- Panic with `ErrNaN("descriptive message")` for NaN conditions (per `big.Float` semantics).
- Do NOT use `big.ErrNaN` — this package has its own `ErrNaN` type defined in `bigmath.go`.
- Document all special cases in doc comments (e.g. `Sin(±0) = ±0`, `Sin(±Inf) = panic(ErrNaN)`).

### 7. Loop Exit Condition
Taylor-series loops exit when `term.Sign() == 0 || term.MantExp(nil) < ULPExponent(sum)`. This avoids unnecessary iterations once the term is below the ULP of the accumulated sum.

### 9. Precision

Always assume that the precision can be anything, not just 128 or 256.
