// SPDX-License-Identifier: MIT

package bigmath

import (
	"math"
)

// atanCore computes arctan(x) for |x| < 1 using the Taylor series:
// arctan(x) = Σ (-1)^k * x^(2k+1) / (2k+1)
func atanCore(z, x *Float) *Float {
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

		if term.Sign() == 0 || term.MantExp(nil) < sum.ULPExponent() {
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
func atanReciprocal(z *Float, n uint64) *Float {
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

		if term.Sign() == 0 || term.MantExp(nil) < sum.ULPExponent() {
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
func (z *Float) Atan(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.IsInf() {
		z.Pi()
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
		t0.FMA(xVal, xVal, one)
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
func (z *Float) Atan2(y, x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = max(y.Prec(), x.Prec())
		z.SetPrec(prec)
	}

	// check y.sign first so that if y == z, we don't lose the sign when setting z = π
	yNeg := y.Signbit()
	if y.Sign() == 0 {
		if x.Signbit() {
			z.Pi()
			if yNeg {
				z.Neg(z)
			}
			return z
		}
		return z.Set(y)
	}

	if x.Sign() == 0 {
		z.Pi()
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
				p := newFloat(prec).Pi()
				q := newFloat(prec).Pi()
				q.SetMantExp(q, -2) // π/4 at same precision as p
				z.Sub(p, q)
			} else {
				z.Pi()
				z.SetMantExp(z, -2) // π/4
			}
			if yNeg {
				z.Neg(z)
			}
			return z
		}
		if y.IsInf() {
			// Atan2(±Inf, x) = ±π/2
			z.Pi()
			z.SetMantExp(z, -1)
			if yNeg {
				z.Neg(z)
			}
			return z
		}
		// x is Inf, y is finite
		if x.Signbit() {
			// Atan2(y, -Inf) = ±π
			z.Pi()
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
	res := newFloat(workPrec).Atan(q)

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

// asinGuard computes the working precision needed for asin(x) near x=±1,
// where 1-x² loses significant bits. Returns the required working precision.
func asinGuard(x *Float, prec uint) uint {
	// Work with |x|
	xAbs := newFloat(prec + 2*_W).Abs(x)

	// For |x| ≤ 0.5, 1-x² ≥ 0.75, so no catastrophic cancellation.
	if xAbs.Cmp(newFloat(0).SetFloat64(0.5)) <= 0 {
		return prec + 2*_W
	}

	// Near ±1: compute 1 - |x| to estimate bit loss from cancellation.
	// Since 1 - x² = (1-x)(1+x), the precision loss from computing 1 - x²
	// directly (as 1 - xVal*xVal) is bounded by the loss from 1 - |x|.
	oneMinus := newFloat(prec+2*_W).Sub(one, xAbs)
	subExp := oneMinus.MantExp(nil)
	if -subExp <= 2 {
		return prec + 2*_W
	}
	need := max(prec+uint(-subExp)+3*_W, prec+2*_W)
	return need
}

// Asin sets z to the rounded value of arc sine of x and returns z.
//
// Special cases:
//
//	Asin(±0) = ±0
//	Asin(±1) = ±π/2
//	Asin(|x| > 1) = panic(ErrNaN)
func (z *Float) Asin(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.Sign() == 0 {
		return z.Set(x)
	}

	// Domain check: |x| ≤ 1
	switch absCmpOne(x) {
	case 1:
		panic(ErrNaN("asin of x outside [-1, 1]"))

	case 0:
		// |x| == 1 → ±π/2
		z.Pi()              // z = π at z's precision
		z.SetMantExp(z, -1) // z = π/2, keeps z's precision
		if x.Signbit() {
			z.Neg(z)
		}
		return z
	}

	workPrec := asinGuard(x, prec)

	xVal := newFloat(workPrec).Set(x)
	neg := xVal.Signbit()
	if neg {
		xVal.Neg(xVal)
	}

	// Asin(x) = Atan(x / sqrt(1 - x²))
	x2 := newFloat(workPrec).Mul(xVal, xVal)
	oneMinusX2 := newFloat(workPrec).Sub(one, x2)
	sqrt := newFloat(workPrec).Sqrt(oneMinusX2)
	ratio := newFloat(workPrec).Quo(xVal, sqrt)
	result := newFloat(workPrec).Atan(ratio)

	if neg {
		result.Neg(result)
	}
	return z.Set(result)
}

// Acos sets z to the rounded value of arc cosine of x and returns z.
//
// Special cases:
//
//	Acos(1) = 0
//	Acos(0) = π/2
//	Acos(-1) = π
//	Acos(|x| > 1) = panic(ErrNaN)
func (z *Float) Acos(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	// Domain check: |x| ≤ 1
	switch absCmpOne(x) {
	case 1:
		panic(ErrNaN("acos of x outside [-1, 1]"))
	case 0:
		if x.Signbit() {
			return z.Pi()
		}
		return z.Set(zero)
	}

	// Acos(x) = π/2 - Asin(x)
	workPrec := asinGuard(x, prec)

	tmp := newFloat(workPrec).Asin(x)
	result := newFloat(workPrec)
	result.Pi()                   // result = π at workPrec
	result.SetMantExp(result, -1) // result = π/2, keeps workPrec
	result.Sub(result, tmp)

	return z.Set(result)
}
