// SPDX-License-Identifier: MIT

package bigmath

import (
	"math"
	"math/big"
)

// atanCore computes arctan(1/n) for n > 0 using the Taylor series
// arctan(1/n) = Σ (-1)^k / ((2k+1) * n^(2k+1))
// Using the identity: arctan(1/n) = 1/n * Σ (-1)^k / ((2k+1) * n^(2k))
func atanCore(z, n *big.Float) *big.Float {
	prec := z.Prec()
	// Guard bits to ensure precision
	workPrec := prec + 2*_W

	// temps
	t0 := newFloat(workPrec)
	t1 := newFloat(workPrec)
	t2 := newFloat(workPrec)

	// term = 1/n
	term := newFloat(workPrec).Quo(one, n)

	// sum = term
	sum := newFloat(workPrec).Set(term)

	// v = 1/n^2
	t0.Mul(n, n)
	v := new(big.Float).Quo(one, t0)

	for i := uint64(1); ; i++ {
		// term = (term * v * (2i-1)) / (2i+1)
		t0.Mul(term, v)

		t1.SetUint64(2*i - 1)
		t2.Mul(t0, t1)

		t1.SetUint64(2*i + 1)
		term.Quo(t2, t1)

		if term.Sign() == 0 || term.MantExp(nil) < ULPExponent(sum) {
			break
		}

		if i%2 != 0 {
			t0.Sub(sum, term)
		} else {
			t0.Add(sum, term)
		}
		sum, t0 = t0, sum
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
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.IsInf() {
		z.Set(pi(z.Prec()))
		z.SetMantExp(z, -1)
		if x.Signbit() {
			z.Neg(z)
		}
		return z
	}

	if x.Sign() == 0 {
		return z.Set(zero)
	}

	// u is the reduction threshold where we'll reduce x until x < 2^(-√(prec)/4)
	// with a minimum of 2^-2 (x < 0.25) to ensure convergence.
	u := max(int(math.Sqrt(float64(prec)))/4, 2)
	// We need 8 bits of added precision plus one bit per reduction step
	// with at least one more reduction step for x.exp > 1
	extraPrec := 8 + max(0, min(x.MantExp(nil), 2)+u)
	workPrec := prec + uint(extraPrec)

	xVal := newFloat(workPrec).Set(x)
	var neg bool
	if xVal.Signbit() {
		neg = true
		xVal.Neg(xVal)
	}

	// Reduction
	nReductions := 0
	t0 := newFloat(workPrec)
	t1 := newFloat(workPrec)
	for xVal.MantExp(nil) > -u {
		// x = x / (1 + sqrt(1+x^2))
		FMA(t0, xVal, xVal, one)
		t1.Sqrt(t0)
		t0.Add(t1, one)
		xVal.Quo(xVal, t0)
		nReductions++
	}

	// Now x is small, use atanCore(1/n) where n = 1/x
	t0.Quo(one, xVal)
	atanCore(t1, t0)

	// Undo the double angle reductions: Atan(x) = 2^n * Atan(x_reduced)
	// z is calculated at workPrec, we scale it then round to prec.
	if nReductions > 0 {
		t1.SetMantExp(t1, nReductions)
	}

	if neg {
		t1.Neg(t1)
	}
	return z.Set(t1)
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

	return z.Set(t)
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
			return z
		}
		return z.Set(y)
	}

	if x.Sign() == 0 {
		Pi(z)
		z.SetMantExp(z, -1)
		if y.Signbit() {
			z.Neg(z)
		}
		return z
	}

	if x.IsInf() {
		if x.Signbit() {
			// Atan2(y, -Inf) = ±π
			Pi(z)
			if y.Signbit() {
				z.Neg(z)
			}
			return z
		}
		// Atan2(y, +Inf) = 0
		return z.Set(zero)
	}

	if y.IsInf() {
		// Atan2(±Inf, x) = ±π/2
		Pi(z)
		z.SetMantExp(z, -1)
		if y.Signbit() {
			z.Neg(z)
		}
		return z
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

	return z
}
