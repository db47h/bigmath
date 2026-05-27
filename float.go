// SPDX-License-Identifier: MIT

package bigmath

import (
	"math/big"
)

// newFloat allocates a new *big.Float with the given precision.
func newFloat(prec uint) *big.Float {
	return new(big.Float).SetPrec(prec)
}

// ULPExponent returns the exponent of the Unit in the Last Place (ULP) of x,
// i.e., the value of the least significant mantissa bit.
func ULPExponent(x *big.Float) int {
	return x.MantExp(nil) - int(x.Prec())
}

// isOdd returns true if x is a non-zero odd integer.
func isOdd(x *big.Float) bool {
	return x.IsInt() && x.Sign() != 0 && x.MantExp(nil) == int(x.MinPrec())
}

// absCmpOne compares |x| to 1.
// Returns -1 if |x| < 1, 0 if |x| == 1, 1 if |x| > 1.
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
