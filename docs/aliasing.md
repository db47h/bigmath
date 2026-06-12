# big.Float Aliasing and Temporary Allocations

A reference for which methods of `*big.Float` support zero-allocation aliasing
(where the receiver and a `*Float` argument point to the same value).

Source: `/usr/lib/go/src/math/big/float.go` (Go standard library).

---

## Legend

- **✅ Safe** — aliasing the receiver with the `*Float` argument does **not** allocate
  temporary mantissa storage. The method operates in-place.
- **❌ Allocates** — aliasing triggers one or more temporary allocations for the
  mantissa (or a shifted copy of it).

---

## Methods where aliasing is allocation-free

These methods, when called as `z.Method(z)` (or `z.Method(z, …)` for two-arg
setters), complete without allocating any new mantissa backing storage.

### `Neg(x *Float) *Float`

```go
func (z *Float) Neg(x *Float) *Float {
    z.Set(x)       // no-op when z == x
    z.neg = !z.neg // flip the sign bit
    return z
}
```

`z.Set(x)` has an explicit `if z != x` guard — when the receiver and argument
are the same, the entire copy body is skipped and only `z.acc = Exact` is set
(which is benign because the value is unchanged). Then `z.neg` is toggled.
**No mantissa bytes are read or written.**

### `Abs(x *Float) *Float`

```go
func (z *Float) Abs(x *Float) *Float {
    z.Set(x)       // no-op when z == x
    z.neg = false  // clear the sign bit
    return z
}
```

Identical pattern to `Neg`.

### `Set(x *Float) *Float`

```go
func (z *Float) Set(x *Float) *Float {
    z.acc = Exact
    if z != x {                    // <-- aliasing guard
        z.form = x.form
        z.neg = x.neg
        if x.form == finite {
            z.exp = x.exp
            z.mant = z.mant.set(x.mant)
        }
        if z.prec == 0 {
            z.prec = x.prec
        } else if z.prec < x.prec {
            z.round(0)
        }
    }
    return z
}
```

When `z == x` the method becomes a no-op (only `z.acc = Exact` runs). When `z
!= x`, `z.mant = z.mant.set(x.mant)` *does* allocate — that's expected because
we're copying into a different value.

### `Copy(x *Float) *Float`

```go
func (z *Float) Copy(x *Float) *Float {
    if z != x {                    // <-- aliasing guard
        z.prec = x.prec
        z.mode = x.mode
        z.acc = x.acc
        z.form = x.form
        z.neg = x.neg
        if z.form == finite {
            z.mant = z.mant.set(x.mant)
            z.exp = x.exp
        }
    }
    return z
}
```

Same guard pattern as `Set`.

### `SetMantExp(mant *Float, exp int) *Float`

The doc comment on the method explicitly says:

> *z and mant may be the same in which case z's exponent is set to exp.*

```go
func (z *Float) SetMantExp(mant *Float, exp int) *Float {
    z.Copy(mant)        // no-op when z == mant
    if z.form == finite {
        z.setExpAndRound(int64(z.exp)+int64(exp), 0)
    }
    return z
}
```

`z.Copy(z)` is a no-op, then `setExpAndRound` adjusts the exponent and calls
`round`, which only shortens or bit-fiddles `z.mant` in-place. **No allocation.**

### `MantExp(mant *Float) (exp int)`

```go
func (x *Float) MantExp(mant *Float) (exp int) {
    if x.form == finite {
        exp = int(x.exp)
    }
    if mant != nil {
        mant.Copy(x)     // no-op when mant == x
        if mant.form == finite {
            mant.exp = 0 // just clears the exponent
        }
    }
    return
}
```

This method has an explicit aliasing test in the standard library:
`TestFloatMantExpAliasing` in `float_test.go`.

### `bigmath.Float` methods

For `bigmath.Float` methods that are allocation-free under aliasing, there
are two categories:

- **Argument transformation** — the method applies an absolute value,
  argument reduction, or another mathematical transform to `x` before
  computing the result. Examples: `Acos`, `Sin`, `Exp`, `Log`, `Pow`.
- **In-place mantissa modification** — `IntRound`, `Floor`, `Ceil`, `Trunc`
  directly truncate the mantissa of `x`.

---

## Methods where aliasing **does** allocate

These methods allocate temporary `nat` storage when the receiver aliases an
argument. In many cases they would produce wrong results without the allocation
(because reading and writing the same slice mid-computation would corrupt
inputs).

### `Add(x, y *Float) *Float` ❌

The internal `uadd` (and `usub` for opposite-sign cases) checks:

```go
al := alias(z.mant, x.mant) || alias(z.mant, y.mant)
```

When `al` is true it allocates a temporary shifted copy:

```go
if al {
    t := nat(nil).lsh(y.mant, uint(ey-ex))  // extra allocation
    z.mant = z.mant.add(x.mant, t)
}
```

The non-aliased path reuses `z.mant` for the shift and then adds into it
directly, avoiding the `t` allocation.

### `Sub(x, y *Float) *Float` ❌

Same internal machinery as `Add` — `usub` has the identical aliasing check and
temporary shift allocation.

### `Mul(x, y *Float) *Float` ❌

Calls `umul`, which invokes `nat.mul` (or `nat.sqr` when `x == y`):

```go
// natmul.go
func (z nat) mul(stk *stack, x, y nat) nat {
    // ...
    if alias(z, x) || alias(z, y) {
        z = nil       // discard receiver — forces fresh allocation
    }
    z = z.make(m + n) // always allocates new backing store
    // ...
}
```

When the receiver's mantissa aliases `x.mant` or `y.mant`, the existing storage
is discarded and entirely new storage is allocated for the result.

### `Quo(x, y *Float) *Float` ❌

Calls `uquo`, which uses `nat.div` → `nat.divLarge`:

```go
// natdiv.go
func (z nat) divLarge(stk *stack, u, uIn, vIn nat) (q, r nat) {
    // ...
    if alias(z, u) {
        z = nil       // alias guard
    }
    q = z.make(m + 1) // fresh allocation
    // ...
}
```

Additionally, `uquo` may allocate `xadj` when the dividend needs extra words
for sufficient result precision.

### `Inv(x *Float) *Float` ❌

`Inv` delegates to `Quo`:

```go
func (z *Float) Inv(x *Float) *Float {
    return z.Quo(one, x)   // one is a package-level *Float with value 1
}
```

When `z == x`, this becomes `z.Quo(one, z)`, hitting the same `nat.div`
aliasing guard described for `Quo` above and forcing a fresh allocation for
the result mantissa.

### `Sqrt(x *Float) *Float` ❌

Internally calls `z.sqrtInverse(z)` when aliased. That function creates fresh
working `Float` values (`u`, `v`, `sqi` — already allocated regardless of
aliasing), but the final step is `z.Mul(x, sqi)`. When `z == x`, this calls
`Mul(z, sqi)`, hitting the same `nat.mul` aliasing guard described above and
forcing a fresh allocation for the result mantissa.

---

## Quick reference table

| Method | Aliasing allocates? |
|---|---|
| `Abs(x)`       | **No** ✅¹ |
| `Acos(x)`      | **No** ✅ |
| `Acosh(x)`     | **No** ✅ |
| `Acot(x)`      | **No** ✅ |
| `Acoth(x)`     | **No** ✅ |
| `Acsc(x)`      | **No** ✅ |
| `Acsch(x)`     | **No** ✅ |
| `Add(x, y)`    | **Yes** ❌² |
| `Asec(x)`      | **No** ✅ |
| `Asech(x)`     | **No** ✅ |
| `Asin(x)`      | **No** ✅ |
| `Asinh(x)`     | **No** ✅ |
| `Atan(x)`      | **No** ✅ |
| `Atan2(y, x)`  | **No** ✅ |
| `Atanh(x)`     | **No** ✅ |
| `Cbrt(x)`      | **No** ✅ |
| `Ceil(x)`      | **No** ✅³ |
| `Copy(x)`      | **No** ✅⁴ |
| `Cos(x)`       | **No** ✅ |
| `Cosh(x)`      | **No** ✅ |
| `Cot(x)`       | **No** ✅ |
| `Coth(x)`      | **No** ✅ |
| `Csc(x)`       | **No** ✅ |
| `Csch(x)`      | **No** ✅ |
| `Exp(x)`       | **No** ✅ |
| `Floor(x)`     | **No** ✅³ |
| `FMA(x, y, t)` | **No** ✅⁵ |
| `FMS(x, y, t)` | **No** ✅⁵ |
| `Gamma(x)`     | **No** ✅ |
| `Hypot(x, y)`  | **No** ✅ |
| `IntRound(x)`  | **No** ✅³ |
| `Inv(x)`       | **Yes** ❌ |
| `Lgamma(x)`    | **No** ✅ |
| `Log(x)`       | **No** ✅ |
| `Log10(x)`     | **No** ✅ |
| `Log2(x)`      | **No** ✅ |
| `MantExp(mant)` | **No** ✅ |
| `Mul(x, y)`    | **Yes** ❌⁶ |
| `Neg(x)`       | **No** ✅¹ |
| `Pow(x, y)`    | **No** ✅ |
| `Quo(x, y)`    | **Yes** ❌⁶ |
| `Sec(x)`       | **No** ✅ |
| `Sech(x)`      | **No** ✅ |
| `Set(x)`       | **No** ✅⁴ |
| `SetMantExp(mant, exp)` | **No** ✅¹ |
| `Sin(x)`       | **No** ✅ |
| `Sinh(x)`      | **No** ✅ |
| `Sqrt(x)`      | **Yes** ❌⁶ |
| `Sub(x, y)`    | **Yes** ❌² |
| `Tan(x)`       | **No** ✅ |
| `Tanh(x)`      | **No** ✅ |
| `Trunc(x)`     | **No** ✅³ |

¹ Internal fields (sign, exponent) are updated in-place without touching the mantissa.
² Allocates a temporary shifted `nat` via `nat(nil).lsh`.
³ Modifies the mantissa in-place (no argument transformation needed).
⁴ No-op when aliased (the `if z != x` guard skips all work).
⁵ FMA/FMS always allocate a temp for the product x*y.
⁶ Forces a fresh `nat.make` allocation because `nat.mul`/`nat.sqr`/`nat.div`
  discard the receiver slice when aliasing is detected.

---

## Why `Neg` and `Abs` can be allocation-free while binary ops can't

The critical difference is that `Neg` and `Abs` only modify the **sign bit**
(`z.neg`) — they never touch the mantissa. The `Set` call at the top is a
no-op when the receiver and argument are the same, so the mantissa stays
untouched.

Binary operations (`Add`, `Sub`, `Mul`, `Quo`, `Sqrt`) must compute a new
mantissa whose size differs from (or overlaps with) the inputs. The
underlying `nat` arithmetic methods detect aliasing between their receiver
and arguments as a correctness requirement — writing into a slice that's
still being read as input would produce corrupt results mid-computation.
