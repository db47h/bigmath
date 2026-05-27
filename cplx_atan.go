// SPDX-License-Identifier: MIT

package bigmath

// Atan sets z to the inverse tangent of x and returns z.
//
// Formula: atan(x) = (i/2) · ln((1−ix)/(1+ix))
//
// The branch cut is along the imaginary axis, outside the interval [-i, +i].
// The real part of the result lies in the interval [-π/2, π/2].
//
// Special cases:
//
//	Atan(0 + i·0) = 0 + i·0
func (z *Complex) Atan(x *Complex) *Complex {
	// atan(z) = (i/2) * ln((1-iz)/(1+iz))
	// let w = (1-iz)/(1+iz)
	// atan(z) = (i/2) * Log(w)
	prec := z.setPrec(x)
	workPrec := prec + _W

	if x.IsReal() {
		z.Real.Atan(&x.Real)
		z.Imag.Set(&x.Imag)
		return z
	}
	if x.Real.Sign() == 0 {
		if x.Imag.absCmpOne() <= 0 {
			z.Real.Set(&x.Real)
			z.Imag.Atanh(&x.Imag)
			return z
		}
		z.Real.Set(halfPi(prec))
		z.Imag.Atanh(newFloat(workPrec).Quo(one, &x.Imag))
		return z
	}
	if x.Real.IsInf() || x.Imag.IsInf() {
		sgn := x.Real.Signbit()
		sgnI := x.Imag.Signbit()
		z.Real.Set(halfPi(prec))
		if sgn {
			z.Real.Neg(&z.Real)
		}
		z.Imag.Set(zero)
		if sgnI {
			z.Imag.Neg(&z.Imag)
		}
		return z
	}

	ix := &Complex{Real: *new(Float).Neg(&x.Imag), Imag: x.Real}

	// num = 1 - ix
	num := newComplex(workPrec)
	num.Real.Sub(one, &ix.Real)
	num.Imag.Neg(&ix.Imag)
	// den = 1 + ix
	den := newComplex(workPrec)
	den.Real.Add(one, &ix.Real)
	den.Imag.Set(&ix.Imag)

	w := newComplex(workPrec).Quo(num, den)
	lw := newComplex(workPrec).Log(w)

	// z = (i/2) * lw = (-lw.Imag/2) + i*(lw.Real/2)
	z.Real.Neg(lw.Imag.SetMantExp(&lw.Imag, -1))
	z.Imag.Set(lw.Real.SetMantExp(&lw.Real, -1))

	return z
}

// Asin sets z to the inverse sine of x and returns z.
//
// Formula: asin(x) = −i · ln(i·x + √(1−x²))
//
// The branch cut is along the real axis, outside the interval [-1, +1].
// The real part of the result lies in the interval [-π/2, π/2].
//
// Special cases:
//
//	Asin(0 + i·0) = 0 + i·0
func (z *Complex) Asin(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	switch {
	case x.IsReal() && x.Real.absCmpOne() <= 0:
		z.Real.Asin(&x.Real)
		z.Imag.Set(&x.Imag)
		return z
	case x.Real.Sign() == 0 && x.Imag.absCmpOne() <= 0:
		z.Real.Set(&x.Real)
		z.Imag.Asinh(&x.Imag)
		return z
	}
	// t0 = x^2
	t0 := newComplex(workPrec).Mul(x, x)
	// t1 = 1 - x^2
	t1 := newComplex(workPrec)
	t1.Real.Sub(one, &t0.Real)
	t1.Imag.Neg(&t0.Imag)
	// t0 = sqrt(1 - x^2)
	t0.Sqrt(t1)

	// t1 = i * x
	ix := &Complex{Real: *new(Float).Neg(&x.Imag), Imag: x.Real}

	// t1 = i·x + √(1−x²)
	t1.Add(ix, t0)

	// res = -i * ln(t2)
	t0.Log(t1)
	z.Real.Set(&t0.Imag)
	z.Imag.Neg(&t0.Real)

	return z
}

// Acos sets z to the inverse cosine of x and returns z.
//
// Formula: acos(x) = −i · ln(x + i·√(1−x²))
//
// The branch cut is along the real axis, outside the interval [-1, +1].
// The real part of the result lies in the interval [0, π].
//
// Special cases:
//
//	Acos(0 + i·0) = π/2 + i·0
func (z *Complex) Acos(x *Complex) *Complex {
	prec := z.setPrec(x)
	workPrec := prec + _W

	switch {
	case x.IsReal() && x.Real.absCmpOne() <= 0:
		z.Real.Acos(&x.Real)
		z.Imag.Set(&x.Imag)
		return z
	case x.Real.Sign() == 0 && x.Imag.absCmpOne() <= 0:
		z.Real.Set(halfPi(prec))
		z.Imag.Asinh(&x.Imag)
		z.Imag.Neg(&z.Imag)
		return z
	}

	// t0 = x^2
	t0 := newComplex(workPrec).Mul(x, x)
	// t1 = 1 - z^2
	t1 := newComplex(workPrec)
	t1.Real.Sub(one, &t0.Real)
	t1.Imag.Neg(&t0.Imag)
	// t0 = sqrt(1 - z^2)
	t0.Sqrt(t1)

	// it = i * t0
	t0.Imag.Neg(&t0.Imag)
	it := &Complex{t0.Imag, t0.Real}

	// t1 = x + it
	t1.Add(x, it)

	// res = -i * ln(t1)
	t0.Log(t1)
	z.Real.Set(&t0.Imag)
	z.Imag.Neg(&t0.Real)

	return z
}

// Asinh sets z to the inverse hyperbolic sine of x and returns z.
//
// Formula: asinh(x) = ln(x + √(x²+1))
//
// The branch cut is along the imaginary axis, outside the interval [-i, +i].
// The imaginary part of the result lies in the interval [-π/2, π/2].
//
// Special cases:
//
//	Asinh(0 + i·0) = 0 + i·0
func (z *Complex) Asinh(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	// t0 = x^2
	t0 := newComplex(workPrec).Mul(x, x)
	// t1 = x^2 + 1
	t1 := newComplex(workPrec) //.Add(t0, oneC)
	t1.Real.Add(&t0.Real, one)
	t1.Imag.Set(&t0.Imag)

	// t0 = sqrt(x^2 + 1)
	t0.Sqrt(t1)
	// t1 = x + t0
	t1.Add(x, t0)

	return z.Log(t1)
}

// Acosh sets z to the inverse hyperbolic cosine of x and returns z.
//
// Formula: acosh(x) = ln(x + √(x−1)·√(x+1))
//
// The branch cut is along the real axis, for x < 1.
// The imaginary part of the result lies in the interval [0, π].
//
// Special cases:
//
//	Acosh(0 + i·0) = 0 + i·π/2
func (z *Complex) Acosh(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	// t0 = sqrt(x - 1)
	t1 := newComplex(workPrec)
	t1.Real.Sub(&x.Real, one)
	t1.Imag.Set(&x.Imag)
	t0 := newComplex(workPrec).Sqrt(t1)
	// t1 = x + 1
	t2 := newComplex(workPrec)
	t2.Real.Add(&x.Real, one)
	t2.Imag.Set(&x.Imag)
	t1.Sqrt(t2)
	// t2 = sqrt(x-1) * sqrt(x+1)
	t2.Mul(t0, t1)
	// log(x + t2)
	return z.Log(t0.Add(x, t2))
}

// Atanh sets z to the inverse hyperbolic tangent of x and returns z.
//
// Formula: atanh(x) = ½ · ln((1+x)/(1−x))
//
// The branch cut is along the real axis, outside the interval [-1, +1].
// The imaginary part of the result lies in the interval [-π/2, π/2].
//
// Special cases:
//
//	Atanh(0 + i·0) = 0 + i·0
func (z *Complex) Atanh(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	// num = 1 + z
	num := newComplex(workPrec)
	num.Real.Add(one, &x.Real)
	num.Imag.Set(&x.Imag)
	// den = 1 - z
	den := newComplex(workPrec)
	den.Real.Sub(one, &x.Real)
	den.Imag.Neg(&x.Imag)

	w := newComplex(workPrec).Quo(num, den)
	// directly use z as target. The z /= 2 is lossless
	z.Log(w)

	// z /= 2
	z.Real.SetMantExp(&z.Real, -1)
	z.Imag.SetMantExp(&z.Imag, -1)

	return z
}
