// SPDX-License-Identifier: MIT

package bigmath

// Sin sets z to the sine of x and returns z.
//
// Formula: sin(a+bi) = sin(a)cosh(b) + i·cos(a)sinh(b)
func (z *Complex) Sin(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	s, c := Sincos(newFloat(workPrec), newFloat(workPrec), &x.Real)
	sh, ch := SinhCosh(newFloat(workPrec), newFloat(workPrec), &x.Imag)

	z.Real.Mul(s, ch)
	z.Imag.Mul(c, sh)

	return z
}

// Cos sets z to the cosine of x and returns z.
//
// Formula: cos(a+bi) = cos(a)cosh(b) − i·sin(a)sinh(b)
func (z *Complex) Cos(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	s, c := Sincos(newFloat(workPrec), newFloat(workPrec), &x.Real)
	sh, ch := SinhCosh(newFloat(workPrec), newFloat(workPrec), &x.Imag)

	z.Real.Mul(c, ch)
	z.Imag.Mul(s, sh).Neg(&z.Imag)

	return z
}

// Tan sets z to the tangent of x and returns z.
func (z *Complex) Tan(x *Complex) *Complex {
	// TODO: take this as an example for future review of guard bits strategy.
	// for example, Sin/Cos already add guard bits, and here we need at least 2 more for Quo.
	workPrec := z.setPrec(x) + _W

	// tan(z) = sin(z) / cos(z)
	s := newComplex(workPrec).Sin(x)
	c := newComplex(workPrec).Cos(x)
	return z.Quo(s, c)
}
