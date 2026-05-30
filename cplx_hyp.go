// SPDX-License-Identifier: MIT

package bigmath

// Sinh sets z to the hyperbolic sine of x and returns z.
//
// Formula: sinh(a+bi) = sinh(a)cos(b) + i·cosh(a)sin(b)
func (z *Complex) Sinh(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	if x.IsReal() {
		z.Real.Sinh(&x.Real)
		z.Imag.Set(&x.Imag)
		return z
	}

	s, c := Sincos(newFloat(workPrec), newFloat(workPrec), &x.Imag)
	sh, ch := SinhCosh(newFloat(workPrec), newFloat(workPrec), &x.Real)

	z.Real.Mul(sh, c)
	z.Imag.Mul(ch, s)

	return z
}

// Cosh sets z to the hyperbolic cosine of x and returns z.
//
// Formula: cosh(a+bi) = cosh(a)cos(b) + i·sinh(a)sin(b)
func (z *Complex) Cosh(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	if x.IsReal() {
		z.Real.Cosh(&x.Real)
		z.Imag.Set(&x.Imag)
		return z
	}

	s, c := Sincos(newFloat(workPrec), newFloat(workPrec), &x.Imag)
	sh, ch := SinhCosh(newFloat(workPrec), newFloat(workPrec), &x.Real)

	z.Real.Mul(ch, c)
	z.Imag.Mul(sh, s)

	return z
}

// Tanh sets z to the hyperbolic tangent of x and returns z.
func (z *Complex) Tanh(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	if x.IsReal() {
		z.Real.Tanh(&x.Real)
		z.Imag.Set(&x.Imag)
		return z
	}

	// tanh(z) = sinh(z) / cosh(z)
	sh := newComplex(workPrec).Sinh(x)
	ch := newComplex(workPrec).Cosh(x)
	return z.Quo(sh, ch)
}

// Sech sets z to the hyperbolic secant of x, sech(x) = 1/cosh(x), and returns z.
func (z *Complex) Sech(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W
	if x.IsReal() {
		z.Real.Sech(&x.Real)
		z.Imag.Set(&x.Imag)
		return z
	}
	t := newComplex(workPrec).Cosh(x)
	return z.Inv(t)
}

// Csch sets z to the hyperbolic cosecant of x, csch(x) = 1/sinh(x), and returns z.
func (z *Complex) Csch(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W
	if x.IsReal() {
		z.Real.Csch(&x.Real)
		z.Imag.Set(&x.Imag)
		return z
	}
	t := newComplex(workPrec).Sinh(x)
	return z.Inv(t)
}
