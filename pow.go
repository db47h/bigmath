// SPDX-License-Identifier: MIT

package bigmath

import (
	"math/big"
)

func isOdd(x *big.Float) bool {
	if !x.IsInt() {
		return false
	}
	var i big.Int
	x.Int(&i)
	return i.Bit(0) != 0
}

func absCmpOne(x *big.Float) int {
	exp := x.MantExp(nil)
	if exp > 1 {
		return 1
	}
	if exp < 1 {
		return -1
	}
	// exp == 1 => 1 <= |x| < 2
	if x.Signbit() {
		return -x.Cmp(minusOne)
	}
	return x.Cmp(one)
}

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
//	Pow(NaN, y) = NaN (panics with ErrNaN)
//	Pow(x, NaN) = NaN (panics with ErrNaN)
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
func Pow(z, x, y *big.Float) *big.Float {
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
			if isOdd(y) {
				return z.SetInf(x.Signbit())
			}
			return z.SetInf(false)
		}
		if isOdd(y) {
			return z.Set(x)
		}
		return z.Set(zero)
	}

	if y.IsInf() {
		cmp := absCmpOne(x)
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
			if x.Signbit() && isOdd(y) {
				return z.SetInf(true)
			}
			return z.SetInf(false)
		}
		// y < 0
		if x.Signbit() && isOdd(y) {
			return z.Set(zero).Neg(z)
		}
		return z.Set(zero)
	}

	if x.Signbit() && !y.IsInt() {
		panic(ErrNaN("Pow(x, y) where x < 0 and y is not an integer"))
	}

	// For integer exponents, use optimizations.
	if y.IsInt() {
		// Detect huge exponents that will surely overflow or underflow
		yExp := y.MantExp(nil)
		if yExp > 64 {
			// y is huge.
			cmp := absCmpOne(x)
			if (y.Sign() > 0 && cmp > 0) || (y.Sign() < 0 && cmp < 0) {
				if x.Signbit() && isOdd(y) {
					return z.SetInf(true)
				}
				return z.SetInf(false)
			}
			// underflow
			if x.Signbit() && isOdd(y) {
				return z.Set(zero).Neg(z)
			}
			return z.Set(zero)
		}

		n := new(big.Int)
		y.Int(n)

		// Binary exponentiation is O(log(n)) multiplications.
		// Exp(y * Log(x)) is O(log(prec)) multiplications.
		if n.BitLen() <= 128 {
			return powInt(z, x, n)
		}

		if x.Signbit() {
			absX := new(big.Float).Abs(x)
			odd := n.Bit(0) != 0
			z = Pow(z, absX, y)
			if odd {
				z.Neg(z)
			}
			return z
		}
	}

	// General case: x^y = exp(y * ln(x))
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
	}

	// Working precision
	workPrec := addPrec(prec, 64)

	l := newFloat(workPrec)
	Log(l, x)
	
	temp := newFloat(workPrec)
	temp.Mul(l, y)
	l, temp = temp, l

	// Set z precision if it was 0
	if z.Prec() == 0 {
		z.SetPrec(prec)
	}

	return Exp(z, l)
}

// powInt computes x^n using exponentiation by squaring.
func powInt(z, x *big.Float, n *big.Int) *big.Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
	}
	workPrec := addPrec(prec, 64)

	neg := n.Sign() < 0
	absN := new(big.Int).Abs(n)

	res := newFloat(workPrec).SetUint64(1)
	temp := newFloat(workPrec)
	base := newFloat(workPrec).Set(x)

	for i := absN.BitLen() - 1; i >= 0; i-- {
		temp.Mul(res, res)
		res, temp = temp, res
		if res.IsInf() {
			break
		}
		if absN.Bit(i) != 0 {
			temp.Mul(res, base)
			res, temp = temp, res
		}
		if res.IsInf() {
			break
		}
	}

	if neg {
		res.Quo(one, res)
	}

	return z.Set(res)
}
