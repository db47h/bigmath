# bigmath.Float functions causing hidden allocs when aliasing

## Legend

- **✅ Safe** — aliasing the receiver with the `*Float` argument does **not** allocate
  temporary mantissa storage. The method operates in-place.
- **❌ Allocates** — aliasing triggers one or more temporary allocations for the
  mantissa (or a shifted copy of it).
- **?** - not checked yet

## Detail

[✅] func (z *Float) Acos(x *Float) *Float
[✅] func (z *Float) Acosh(x *Float) *Float
[✅] func (z *Float) Acot(x *Float) *Float
[✅] func (z *Float) Acoth(x *Float) *Float
[✅] func (z *Float) Acsc(x *Float) *Float
[✅] func (z *Float) Acsch(x *Float) *Float
[✅] func (z *Float) Asec(x *Float) *Float
[✅] func (z *Float) Asech(x *Float) *Float
[✅] func (z *Float) Asin(x *Float) *Float
[✅] func (z *Float) Asinh(x *Float) *Float
[✅] func (z *Float) Atan(x *Float) *Float
[✅] func (z *Float) Atan2(y, x *Float) *Float
[✅] func (z *Float) Atanh(x *Float) *Float
[✅] func (z *Float) Cbrt(x *Float) *Float
[✅] func (z *Float) Ceil(x *Float) *Float
[✅] func (z *Float) Cos(x *Float) *Float
[✅] func (z *Float) Cosh(x *Float) *Float
[✅] func (z *Float) Cot(x *Float) *Float
[✅] func (z *Float) Coth(x *Float) *Float
[✅] func (z *Float) Csc(x *Float) *Float
[✅] func (z *Float) Csch(x *Float) *Float
[✅] func (z *Float) Exp(x *Float) *Float
[✅] func (z *Float) FMA(x, y, t *Float)¹ *Float
[✅] func (z *Float) FMS(x, y, t *Float)¹ *Float
[✅] func (z *Float) Floor(x *Float) *Float
[✅] func (z *Float) Gamma(x *Float) *Float
[✅] func (z *Float) Hypot(x, y *Float) *Float
[✅] func (z *Float) IntRound(x *Float) *Float
[❌] func (z *Float) Inv(x *Float) *Float
[✅] func (z *Float) Lgamma(x *Float) (*Float, int)
[✅] func (z *Float) Log(x *Float) *Float
[✅] func (z *Float) Log10(x *Float) *Float
[✅] func (z *Float) Log2(x *Float) *Float
[?] func (z *Float) Pow(x, y *Float) *Float
[?] func (z *Float) Sec(x *Float) *Float
[?] func (z *Float) Sech(x *Float) *Float
[?] func (z *Float) Sin(x *Float) *Float
[?] func (z *Float) Sinh(x *Float) *Float
[?] func (z *Float) Tan(x *Float) *Float
[?] func (z *Float) Tanh(x *Float) *Float
[✅] func (z *Float) Trunc(x *Float) *Float

¹: FMA/FMS always allocate a temp for the product x*y.
