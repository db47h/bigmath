// SPDX-License-Identifier: MIT

package bigmath

// Cbrt sets z to the cube root of x and returns z.
//
// If z's precision is 0, it is changed to x's precision before the operation.
// Rounding is performed according to z's precision and rounding mode.
//
// Special cases:
//
//	Cbrt(±0) = ±0
//	Cbrt(+Inf) = +Inf
//	Cbrt(-Inf) = -Inf
func (z *Float) Cbrt(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.Sign() == 0 {
		return z.Set(x)
	}
	if x.IsInf() {
		return z.SetInf(x.Signbit())
	}

	// Guard bits for intermediate calculations
	workPrec := prec + _W

	neg := x.Signbit()
	// allocate t0 with workPrec bits but use only prec for xAbs
	t0 := newFloat(workPrec).SetPrec(prec)
	t0.Abs(x)

	// Compute cbrt(|x|) = exp(ln(|x|) / 3)
	t1 := newFloat(workPrec).Log(t0)
	t0.SetPrec(workPrec).Quo(t1, three)
	z.Exp(t0)

	if neg {
		z.Neg(z)
	}
	return z
}
