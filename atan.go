// SPDX-License-Identifier: MIT

package bigmath

import (
	"math/big"
)

// atanCore computes arctan(1/n) for n > 0 using the Taylor series
// arctan(1/n) = Σ (-1)^k / ((2k+1) * n^(2k+1))
// Using the identity: arctan(1/n) = 1/n * Σ (-1)^k / ((2k+1) * n^(2k))
func atanCore(z, n *big.Float) *big.Float {
	prec := z.Prec()
	// Guard bits to ensure precision
	workPrec := prec + _W

	// temps
	t := newFloat(workPrec)
	t2 := newFloat(workPrec)

	// term = 1/n
	term := newFloat(workPrec).Quo(one, n)

	// sum = term
	sum := newFloat(workPrec).Set(term)

	// v = 1/n^2
	v := newFloat(workPrec).Mul(n, n)
	v.Quo(one, v)

	for i := uint64(1); ; i++ {
		// term = (term * v * (2i-1)) / (2i+1)
		t.Mul(term, v)

		t2.SetUint64(2*i - 1)
		t.Mul(t, t2)

		t2.SetUint64(2*i + 1)
		term.Quo(t, t2)

		if term.Sign() == 0 || term.MantExp(nil) < ULPExponent(sum) {
			break
		}

		if i%2 != 0 {
			sum.Sub(sum, term)
		} else {
			sum.Add(sum, term)
		}
	}
	return z.Set(sum)
}

// Atan sets z to the rounded value of arctan(x) and returns z.
//
// Special cases:
//
//	Atan(±0) = ±0
//	Atan(±Inf) = ±π/2
func Atan(z, x *big.Float) *big.Float {
	if x.IsInf() {
		Pi(z)
		z.SetMantExp(z, -1)
		if x.Signbit() {
			z.Neg(z)
		}
		return z.SetPrec(z.Prec()) // Ensure rounded to precision
	}
	if x.Sign() == 0 {
		return z.Set(zero)
	}

	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	workPrec := prec + _W

	xVal := new(big.Float).SetPrec(workPrec).Copy(x)
	var neg bool
	if xVal.Signbit() {
		neg = true
		xVal.Neg(xVal)
	}

	// Reduction
	nReductions := 0
	limit := new(big.Float).SetFloat64(0.1)
	for xVal.Cmp(limit) > 0 {
		// x = x / (1 + sqrt(1+x^2))
		t := newFloat(workPrec).Mul(xVal, xVal)
		t.Add(t, one)
		t.Sqrt(t)
		t.Add(t, one)
		xVal.Quo(xVal, t)
		nReductions++
	}

	// Now x is small, use atanCore(1/n) where n = 1/x
	n := newFloat(workPrec).Quo(one, xVal)
	atanCore(z, n)

	// Undo the double angle reductions: Atan(x) = 2^n * Atan(x_reduced)
	// z is calculated at workPrec, we scale it then round to prec.
	if nReductions > 0 {
		z.SetMantExp(z, nReductions)
	}

	if neg {
		z.Neg(z)
	}
	return z.SetPrec(prec)
}

// Hypot sets z to sqrt(x*x + y*y) and returns z.
//
// Special cases are:
//
//	Hypot(±Inf, y) = +Inf
//	Hypot(x, ±Inf) = +Inf
//	Hypot(NaN, y) = NaN
//	Hypot(x, NaN) = NaN
func Hypot(z, x, y *big.Float) *big.Float {
	if x.IsInf() || y.IsInf() {
		return z.SetInf(false)
	}

	if x.Sign() == 0 {
		return z.Abs(y)
	}
	if y.Sign() == 0 {
		return z.Abs(x)
	}

	prec := z.Prec()
	if prec == 0 {
		prec = max(x.Prec(), y.Prec())
		z.SetPrec(prec)
	}
	workPrec := prec + _W

	// Use FMA for better precision: sqrt(x*x + y*y)
	t := newFloat(workPrec)
	t2 := newFloat(workPrec).Mul(y, y)
	FMA(t, x, x, t2)
	t.Sqrt(t)

	return z.Set(t).SetPrec(prec)
}

// Atan2 sets z to the arc tangent of y/x, using the signs of the
// arguments to determine the quadrant of the result.
//
// Special cases:
//
//	Atan2(+0, x>=0) = +0
//	Atan2(-0, x>=0) = -0
//	Atan2(+0, x<0) = +π
//	Atan2(-0, x<0) = -π
//	Atan2(y>0, 0) = +π/2
//	Atan2(y<0, 0) = -π/2
//	Atan2(+Inf, x) = +π/2
//	Atan2(-Inf, x) = -π/2
//	Atan2(y, +Inf) = 0
//	Atan2(y, -Inf) = ±π
func Atan2(z, y, x *big.Float) *big.Float {
	prec := z.Prec()
	if prec == 0 {
		prec = max(y.Prec(), x.Prec())
		z.SetPrec(prec)
	}

	if y.Sign() == 0 {
		if x.Signbit() {
			Pi(z)
			if y.Signbit() {
				z.Neg(z)
			}
			return z.SetPrec(prec)
		}
		return z.Set(y).SetPrec(prec)
	}

	if x.Sign() == 0 {
		Pi(z)
		z.SetMantExp(z, -1)
		if y.Signbit() {
			z.Neg(z)
		}
		return z.SetPrec(prec)
	}

	if x.IsInf() {
		if x.Signbit() {
			// Atan2(y, -Inf) = ±π
			Pi(z)
			if y.Signbit() {
				z.Neg(z)
			}
			return z.SetPrec(prec)
		}
		// Atan2(y, +Inf) = 0
		return z.Set(zero).SetPrec(prec)
	}

	if y.IsInf() {
		// Atan2(±Inf, x) = ±π/2
		Pi(z)
		z.SetMantExp(z, -1)
		if y.Signbit() {
			z.Neg(z)
		}
		return z.SetPrec(prec)
	}

	// Atan2(y, x) = Atan(y/x) + quadrant adjustment
	workPrec := prec + _W
	q := newFloat(workPrec).Quo(y, x)
	Atan(z, q)

	if x.Signbit() {
		p := newFloat(workPrec)
		Pi(p)
		if y.Signbit() {
			z.Sub(z, p)
		} else {
			z.Add(z, p)
		}
	}

	return z.SetPrec(prec)
}
