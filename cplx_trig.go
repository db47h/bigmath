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

// Cot sets z to the cotangent of x and returns z.
//
// Formula: cot(z) = cos(z) / sin(z)
func (z *Complex) Cot(x *Complex) *Complex {
	// Extra guard word: Cos/Sin internally compute at prec+_W and store the
	// result at prec+_W precision (setPrec only increases, never decreases).
	// Adding a second guard word here ensures Quo operates on Cos/Sin results
	// with 2*_W guard bits, keeping the final rounding to z.Prec() within
	// 0.5 ULP of gmpy2 reference.
	workPrec := z.setPrec(x) + 2*_W

	c := newComplex(workPrec).Cos(x)
	s := newComplex(workPrec).Sin(x)
	return z.Quo(c, s)
}

// Sec sets z to the secant of x, sec(x) = 1/cos(x), and returns z.
func (z *Complex) Sec(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W
	if x.IsReal() {
		z.Real.Sec(&x.Real)
		z.Imag.Set(&x.Imag)
		return z
	}
	t := newComplex(workPrec).Cos(x)
	return z.Inv(t)
}

// Csc sets z to the cosecant of x, csc(x) = 1/sin(x), and returns z.
func (z *Complex) Csc(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W
	if x.IsReal() {
		z.Real.Csc(&x.Real)
		z.Imag.Set(&x.Imag)
		return z
	}
	t := newComplex(workPrec).Sin(x)
	return z.Inv(t)
}

// Coth sets z to the hyperbolic cotangent of x, coth(x) = 1/tanh(x), and returns z.
func (z *Complex) Coth(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W
	if x.IsReal() {
		z.Real.Coth(&x.Real)
		z.Imag.Set(&x.Imag)
		return z
	}
	t := newComplex(workPrec).Tanh(x)
	return z.Inv(t)
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

// Acot sets z to the inverse cotangent of x, acot(z) = π/2 - atan(z), and returns z.
//
// The branch cut is along the imaginary axis, outside the interval [-i, +i].
// The real part of the result lies in the interval [0, π].
func (z *Complex) Acot(x *Complex) *Complex {
	// Atan internally computes at prec+_W, then result stored at workPrec.
	prec := z.setPrec(x)
	workPrec := prec + _W

	// x.IsReal shortcut: Float.Acot (Atan2(one, x)) handles all real x including 0.
	if x.IsReal() {
		z.Real.Acot(&x.Real)
		z.Imag.Set(&x.Imag)
		return z
	}
	t := newComplex(workPrec).Atan(x)
	// acot(z) = π/2 - atan(z)
	z.Real.Sub(halfPi.get(workPrec), &t.Real)
	z.Imag.Neg(&t.Imag)
	return z
}

// Asec sets z to the inverse secant of x, asec(z) = acos(1/z), and returns z.
//
// The branch cut is along the real axis, in the interval [-1, +1].
// The real part of the result lies in the interval [0, π].
func (z *Complex) Asec(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W
	// Only use Float shortcut when |x.Real| >= 1 (Float.Asec panics for |x| < 1;
	// the complex formula handles all inputs via Acos(Inv(z))).
	if x.IsReal() && x.Real.absCmpOne() >= 0 {
		z.Real.Asec(&x.Real)
		z.Imag.Set(&x.Imag)
		return z
	}
	return z.Acos(newComplex(workPrec).Inv(x))
}

// Acsc sets z to the inverse cosecant of x, acsc(z) = asin(1/z), and returns z.
//
// The branch cut is along the real axis, in the interval [-1, +1].
// The real part of the result lies in the interval [-π/2, π/2].
func (z *Complex) Acsc(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W
	// Only use Float shortcut when |x.Real| >= 1 (Float.Acsc panics for |x| < 1;
	// the complex formula handles all inputs via Asin(Inv(z))).
	if x.IsReal() && x.Real.absCmpOne() >= 0 {
		z.Real.Acsc(&x.Real)
		z.Imag.Set(&x.Imag)
		return z
	}
	return z.Asin(newComplex(workPrec).Inv(x))
}
