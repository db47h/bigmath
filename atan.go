// SPDX-License-Identifier: MIT

package bigmath

import (
	"math"
	"math/big"
)

// atanCore computes arctan(x) for |x| < 1 using the Taylor series:
// arctan(x) = Σ (-1)^k * x^(2k+1) / (2k+1)
func atanCore(z, x *big.Float) *big.Float {
	prec := z.Prec()
	// Guard bits to ensure precision
	workPrec := prec + 2*_W

	// temps
	t0 := newFloat(workPrec)
	t1 := newFloat(workPrec)
	t2 := newFloat(workPrec)

	// term = x
	term := newFloat(workPrec).Set(x)

	// sum = term
	sum := newFloat(workPrec).Set(term)

	// v = x^2
	v := newFloat(workPrec).Mul(x, x)

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

// atanReciprocal computes arctan(1/n) for n > 0 using the Taylor series.
// This version is optimized for small integer n, avoiding full multiplications
// in the loop by using divisions by n^2.
func atanReciprocal(z *big.Float, n uint64) *big.Float {
	prec := z.Prec()
	workPrec := prec + 2*_W

	n2 := n * n
	bigN2 := newFloat(workPrec).SetUint64(n2)

	// term = 1/n
	term := newFloat(workPrec).Quo(one, newFloat(workPrec).SetUint64(n))
	sum := newFloat(workPrec).Set(term)

	t0 := newFloat(workPrec)
	t1 := newFloat(workPrec)
	t2 := newFloat(workPrec)

	for i := uint64(1); ; i++ {
		// term = (term * (2i-1)) / (n^2 * (2i+1))
		t0.SetUint64(2*i - 1)
		t1.Mul(term, t0) // Fast O(P)

		t0.SetUint64(2*i + 1)
		t2.Mul(t0, bigN2) // Fast O(1)
		term.Quo(t1, t2)  // O(P^2) but avoids full Mul

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
		Pi(z)
		z.SetMantExp(z, -1)
		if x.Signbit() {
			z.Neg(z)
		}
		return z
	}

	if x.Sign() == 0 {
		return z.Set(x)
	}

	// u is the reduction threshold where we'll reduce x until x < 2^-u.
	u := max(int(math.Sqrt(float64(prec)))/4, 2)
	// We need some added precision plus one bit per reduction step.
	workPrec := prec + 8 + uint(max(0, min(x.MantExp(nil), 2)+u))

	xVal := newFloat(workPrec).Set(x)
	var neg bool
	if xVal.Signbit() {
		neg = true
		xVal.Neg(xVal)
	}

	// Reduction identity: arctan(x) = 2 arctan(x / (1 + sqrt(1+x^2)))
	nReductions := 0
	t0 := newFloat(workPrec)
	t1 := newFloat(workPrec)
	for xVal.MantExp(nil) > -u {
		FMA(t0, xVal, xVal, one)
		t1.Sqrt(t0)
		t0.Add(t1, one)
		xVal.Quo(xVal, t0)
		nReductions++
	}

	atanCore(t0, xVal)

	// Undo the double angle reductions: Atan(x) = 2^n * Atan(x_reduced)
	if nReductions > 0 {
		t0.SetMantExp(t0, nReductions)
	}

	if neg {
		t0.Neg(t0)
	}
	return z.Set(t0)
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
	prec := z.Prec()
	if prec == 0 {
		prec = max(x.Prec(), y.Prec())
		z.SetPrec(prec)
	}

	if x.IsInf() || y.IsInf() {
		return z.SetInf(false)
	}

	if x.Sign() == 0 {
		return z.Abs(y)
	}
	if y.Sign() == 0 {
		return z.Abs(x)
	}

	workPrec := prec + _W

	// Use FMA for better precision: sqrt(x*x + y*y)
	t := FMA(newFloat(workPrec), x, x, newFloat(2*y.Prec()).Mul(y, y))

	return z.Sqrt(t)
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
//	Atan2(+Inf, +Inf) = +π/4
//	Atan2(-Inf, +Inf) = -π/4
//	Atan2(+Inf, -Inf) = +3π/4
//	Atan2(-Inf, -Inf) = -3π/4
//	Atan2(y, +Inf) = 0
//	Atan2(y, -Inf) = ±π
//	Atan2(+Inf, x) = +π/2
//	Atan2(-Inf, x) = -π/2
func Atan2(z, y, x *big.Float) *big.Float {
	prec := z.Prec()
	if prec == 0 {
		prec = max(y.Prec(), x.Prec())
		z.SetPrec(prec)
	}

	// check y.sign first so that if y == z, we don't lose the sign when setting z = π
	yNeg := y.Signbit()
	if y.Sign() == 0 {
		if x.Signbit() {
			Pi(z)
			if yNeg {
				z.Neg(z)
			}
			return z
		}
		return z.Set(y)
	}

	if x.Sign() == 0 {
		Pi(z)
		z.SetMantExp(z, -1)
		if yNeg {
			z.Neg(z)
		}
		return z
	}

	if x.IsInf() || y.IsInf() {
		if x.IsInf() && y.IsInf() {
			if x.Signbit() {
				prec := z.Prec() + _W
				p := newFloat(prec)
				q := newFloat(prec)
				Pi(p)
				Pi(q)
				q.SetMantExp(q, -2) // π/4 at same precision as p
				z.Sub(p, q)
			} else {
				Pi(z)
				z.SetMantExp(z, -2) // π/4
			}
			if yNeg {
				z.Neg(z)
			}
			return z
		}
		if y.IsInf() {
			// Atan2(±Inf, x) = ±π/2
			Pi(z)
			z.SetMantExp(z, -1)
			if yNeg {
				z.Neg(z)
			}
			return z
		}
		// x is Inf, y is finite
		if x.Signbit() {
			// Atan2(y, -Inf) = ±π
			Pi(z)
			if yNeg {
				z.Neg(z)
			}
			return z
		}
		// Atan2(y, +Inf) = 0
		return z.Set(zero)
	}

	// Atan2(y, x) = Atan(y/x) + quadrant adjustment
	workPrec := prec + _W
	q := newFloat(workPrec).Quo(y, x)
	res := Atan(newFloat(workPrec), q)

	if x.Signbit() {
		p := pi(workPrec)
		if yNeg {
			q.Sub(res, p)
		} else {
			q.Add(res, p)
		}
		res, q = q, res
		_ = q // quiet warnings about unused q.
	}

	return z.Set(res)
}
