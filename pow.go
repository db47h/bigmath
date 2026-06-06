// SPDX-License-Identifier: MIT

package bigmath

import (
	"math/big"
)

// Pow sets z to the rounded value of x^y and returns z.
//
// If z's precision is 0, it is changed to x's precision before the operation.
// Rounding is performed according to z's precision and rounding mode.
//
// Special cases:
//
//	Pow(x, ±0) = 1 for any x
//	Pow(1, y) = 1 for any y
//	Pow(x, 1) = x for any x
//	Pow(±0, y) where y is an odd integer:
//	    Pow(±0, y) = ±Inf for y < 0
//	    Pow(±0, y) = ±0 for y > 0
//	Pow(±0, y) where y is not an odd integer:
//	    Pow(±0, y) = +Inf for y < 0
//	    Pow(±0, y) = +0 for y > 0
//	Pow(-1, ±Inf) = 1
//	Pow(x, +Inf) for |x| > 1 is +Inf
//	Pow(x, +Inf) for |x| < 1 is +0
//	Pow(x, -Inf) for |x| > 1 is +0
//	Pow(x, -Inf) for |x| < 1 is +Inf
//	Pow(+Inf, y) for y > 0 is +Inf
//	Pow(+Inf, y) for y < 0 is +0
//	Pow(-Inf, y) where y is an odd integer:
//	    Pow(-Inf, y) = -Inf for y > 0
//	    Pow(-Inf, y) = -0 for y < 0
//	Pow(-Inf, y) where y is not an odd integer:
//	    Pow(-Inf, y) = +Inf for y > 0
//	    Pow(-Inf, y) = +0 for y < 0
//
// For finite x < 0 and finite non-integer y, Pow panics with ErrNaN.
func (z *Float) Pow(x, y *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = max(x.Prec(), y.Prec())
		z.SetPrec(prec)
	}

	if y.Sign() == 0 {
		return z.Set(one)
	}
	if x.Cmp(one) == 0 {
		return z.Set(one)
	}
	if y.Cmp(one) == 0 {
		return z.Set(x)
	}

	if x.Sign() == 0 {
		if y.Sign() < 0 {
			if y.isOdd() {
				return z.SetInf(x.Signbit())
			}
			return z.SetInf(false)
		}
		if y.isOdd() {
			return z.Set(x)
		}
		return z.Set(zero)
	}

	if y.IsInf() {
		cmp := x.absCmpOne()
		if cmp == 0 {
			return z.Set(one)
		}
		if (y.Sign() > 0 && cmp > 0) || (y.Sign() < 0 && cmp < 0) {
			return z.SetInf(false)
		}
		return z.Set(zero)
	}

	if x.IsInf() {
		if y.Sign() > 0 {
			if x.Signbit() && y.isOdd() {
				return z.SetInf(true)
			}
			return z.SetInf(false)
		}
		// y < 0
		if x.Signbit() && y.isOdd() {
			return z.Set(zero).Neg(z)
		}
		return z.Set(zero)
	}

	if x.Signbit() && !y.IsInt() {
		panic(ErrNaN("Pow(x, y) where x < 0 and y is not an integer"))
	}

	// For integer exponents, use optimizations.
	if y.IsInt() {
		// Detect huge exponents that will surely overflow or underflow.
		// yExp > 64 implies |y| >= 2^64 (as f*2^e with 0.5 <= f < 1).
		yExp := y.MantExp(nil)
		if yExp > 64 {
			cmp := x.absCmpOne()
			if (y.Sign() > 0 && cmp > 0) || (y.Sign() < 0 && cmp < 0) {
				if x.Signbit() && y.isOdd() {
					return z.SetInf(true)
				}
				return z.SetInf(false)
			}
			// underflow
			if x.Signbit() && y.isOdd() {
				return z.Set(zero).Neg(z)
			}
			return z.Set(zero)
		}

		// Since yExp <= 64, |y| fits in a big.Int with BitLen() <= 64.
		// Binary exponentiation is efficient for this range.
		n := new(big.Int)
		y.Int(n)
		return z.powInt(x, n)
	}

	// General case: x^y = exp(y * ln(x))
	// Working precision: add at least one word of guard bits.
	prec += _W
	l := newFloat(prec).Log(x)
	t := newFloat(prec).Mul(l, y)
	return z.Set(l.Exp(t))
}

// powInt computes x^n using exponentiation by squaring.
func (z *Float) powInt(x *Float, n *big.Int) *Float {
	prec := z.Prec()
	workPrec := prec + _W

	neg := n.Sign() < 0
	absN := new(big.Int).Abs(n)

	res := newFloat(workPrec).Set(one)
	temp := newFloat(workPrec)
	z.SetPrec(x.Prec()).Abs(x)

	for i := absN.BitLen() - 1; i >= 0; i-- {
		temp.Mul(res, res)
		res, temp = temp, res

		if absN.Bit(i) != 0 {
			temp.Mul(res, z)
			res, temp = temp, res
		}
		if res.IsInf() {
			break
		}
	}

	z.SetPrec(prec)
	if neg {
		return z.Inv(res)
	}

	if x.Sign() < 0 && absN.Bit(0) != 0 {
		return z.Neg(res)
	}
	return z.Set(res)
}
