// SPDX-License-Identifier: MIT

package bigmath

// Exp sets z to e^x, the base-e exponential of x, and returns z.
//
// Special cases are:
//
//	Exp(-Inf + i·y) = 0                         (negative infinite real part)
//	Exp(x + i·±Inf) = panic                     (infinite imaginary part)
//	Exp(x + i·0)    = Exp(x)                    (for any x)
//	Exp(+Inf + i·y) = +Inf·(cos(y) + i·sin(y))  (for any finite y)
func (z *Complex) Exp(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	if x.Imag.IsInf() {
		if x.Real.IsInf() && x.Real.Signbit() {
			// -Inf + i·±Inf, continuity with Exp(-Inf + i·y) with y finite.
			z.Real.Set(zero)
			z.Imag.Set(zero)
			return z
		}
		panic(ErrNaN("complex exponential of infinite imaginary part"))
	}
	expA := Exp(newFloat(workPrec), &x.Real)
	if x.Imag.Sign() == 0 {
		z.Real.Set(expA)
		z.Imag.Set(&x.Imag)
		return z
	}
	s, c := Sincos(newFloat(workPrec), newFloat(workPrec), &x.Imag)
	if expA.IsInf() {
		switch {
		case expA.Signbit():
			z.Real.Set(zero)
			z.Imag.Set(zero)
			return z
		case c.Sign() == 0:
			// in this unlikely case, choose geometric continuity over panicking
			z.Real.SetInf(false)
			z.Imag.SetInf(s.Signbit())
			return z
		}
	}

	z.Real.Mul(expA, c)
	z.Imag.Mul(expA, s)

	return z
}

// Log sets z to the principal logarithm of x and returns z.
//
// The branch cut is along the negative real axis. The imaginary part of
// the result lies in the interval [-π, π].
//
// Special cases:
//
//	Log(0) = -Inf + i·0
//	Log(+Inf + i·y) = +Inf + i·0
//	Log(-Inf + i·y) = +Inf + i·π for finite y
//	Log(x + i·±Inf) = +Inf + i·π/2 for finite x
func (z *Complex) Log(x *Complex) *Complex {
	prec := z.setPrec(x)
	workPrec := prec + _W

	// ln(x+iy) = ln|x+iy| + i*arg(x+iy)
	// ln|x+iy| = 0.5 * ln(x^2 + y^2)
	// We use Abs and then real Log to avoid precision loss.
	abs := x.Abs(newFloat(workPrec))
	x.Arg(&z.Imag)
	Log(&z.Real, abs)

	return z
}
