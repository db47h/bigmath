# bigmath Context — Domain & Coding Conventions

## Project

Arbitrary-precision math extending `big.Float` with trig, exp, log, hyperbolic, and complex functions — same rounding/precision semantics as `math/big`.

## Key Design Decisions

### Precision Propagation
Every function follows: `prec := z.Prec(); if prec == 0 { prec = x.Prec() }`. Compute at `workPrec := prec + guard`, round down to `z.Prec()` on return.

### Fused Multiply-Add
**`fma(z, x, y, t, temp)` is correct. Do not alter.**

True FMA semantics (one rounding for `x*y + t`):
- `big.Float.Mul` always computes the full mantissa product internally — summing the operand precisions for `temp` just stores the unrounded result.
- The `Add` in `fma` rounds once to `z.Prec()`.
- `Complex.Mul` follows the same principle for its `bd`/`bc` intermediates.

### Aliasing Guarantee
All functions support `z == x` correctly. Avoid passing the receiver as an argument to `big.Float` methods (Add, Sub, Mul, Quo) because they allocate internally when receiver aliases an argument. Allocation-free when aliasing: `Neg`, `Set`, `SetMantExp`, `SetPrec`, `SetUint64`, `SetInt64`, and read-only methods (`Sign`, `Signbit`, `MantExp(nil)`, `IsInf`, `Prec`, `Cmp`).

### Constant Cache
Computed constants via thread-safe `constCache`: `pi.get(prec)`, `ln2.get(prec)`, `ln10.get(prec)`, `sqrt2.get(prec)`. Returned pointer is **read-only** — mutation corrupts the cache.

 ** To make a mutable copy of computed constants, use**:
```go
mutablePi := newFloat(prec).setConst(pi)
```

### Guard-Bit Heuristics
Each function uses a profiled guard strategy. Convention: core Taylor loops use a flat guard (`prec + _W` typically); argument-reduction layers add precision proportional to reduction depth. See individual `.go` files for specifics.

## Coding Style Rules

Rules MUST be followed in all new code. The workflow prompt instructs subagents to read this file.

### 1. SPDX License Header
```
// SPDX-License-Identifier: MIT

package bigmath
```

### 2. Avoid Aliasing in Internal `big.Float` Calls
- API functions support aliasing correctly — the concern is performance, not correctness.
- Avoid calling `big.Float` methods with the receiver aliasing an argument, which forces a hidden copy internally.
- Exceptions (allocation-free when aliasing): `Neg`, `Set`, `SetMantExp`, `SetPrec`, `SetUint64`, `SetInt64`, and read-only methods (`Sign`, `Signbit`, `MantExp(nil)`, `IsInf`, `Prec`, `Cmp`).
- When aliasing is unavoidable, copy inputs to temps at the top and defer `z.Set(result)` to the end.

### 3. Temporary Allocation
- **Never allocate inside a loop body.** Use pointer swapping (`sum, t0 = t0, sum`) to recycle temps.
- Outside loops, inline allocations are fine — escape analysis keeps `new` on the stack when the result doesn't escape.
- Default allocation: `newFloat(workPrec)`. This is correct for most scratch values.
- `new(big.Float)` is correct when precision is handled elsewhere:
  - Temp passed to `fma` — `fma` resets precision internally.
  - `*new(big.Float).Neg(op)` in `Complex` i·value construction — precision is inherited from `op`.
- Pattern for loops:
  ```go
  t0 := newFloat(workPrec)
  t1 := newFloat(workPrec)
  for ... {
      // use t0, t1, swap pointers instead of re-allocating
      t0, t1 = t1, t0
  }
  ```

### 4. Use Project Constants, Don't Recompute
See [Constant Cache](#constant-cache) above for API. Use cached constants instead of recomputing values like π.

### 5. Internal Helper Signature
Caller-controlled allocation: pass explicit `temp` or work-precision param. Self-contained helpers allocate internally.

### 6. Error Handling
- Panic with `ErrNaN("message")` for NaN conditions.
- Document special cases in doc comments (e.g. `Sin(±0) = ±0`).

### 7. Loop Exit Condition
Taylor-series loops exit when `term.Sign() == 0 || term.MantExp(nil) < ULPExponent(sum)`. No extra iterations once the term is below the accumulated sum's ULP.

### 8. Precision
Assume arbitrary precision (not just 128 or 256).
