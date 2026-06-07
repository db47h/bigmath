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
		// we use z as a temp, delay z.SetPrec
	}

	if x.Sign() == 0 {
		return z.SetPrec(prec).Set(x)
	}
	if x.IsInf() {
		return z.SetPrec(prec).SetInf(x.Signbit())
	}

	// Guard bits for intermediate calculations
	workPrec := prec + _W

	neg := x.Signbit()
	z.SetPrec(x.Prec()).Abs(x)

	// Compute cbrt(|x|) = exp(ln(|x|) / 3)
	t0 := newFloat(workPrec).Log(z)
	t1 := newFloat(workPrec).Quo(t0, three)
	z.SetPrec(prec).Exp(t1)

	if neg {
		z.Neg(z)
	}
	return z
}
