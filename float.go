// SPDX-License-Identifier: MIT

package bigmath

import (
	"math/big"
)

//go:generate go run internal/gen/main.go

// Float is a drop-in replacement for [big.Float] that adds transcendentals,
// power, and trigonometric functions via the [bigmath] package.
//
// It has the same memory layout as [big.Float] and forwards all [big.Float]
// methods via generated code (see float_gen.go), so it can be used as a
// drop-in replacement. Conversion between *Float and *big.Float is zero-cost:
//
//	f := (*Float)(bf)
//	bf := (*big.Float)(f)
//
// Functions like [Sin], [Cos], [Exp], [Log], [Pow], and [Pi] accept
// *Float arguments and follow the same precision and rounding semantics
// as [big.Float].
type Float big.Float

// newFloat allocates a new *Float with the given precision.
func newFloat(prec uint) *Float {
	return new(Float).SetPrec(prec)
}

// ULPExponent returns the exponent of the Unit in the Last Place (ULP) of x,
// i.e., the value of the least significant mantissa bit.
func ULPExponent(x *Float) int {
	return x.MantExp(nil) - int(x.Prec())
}

// isOdd returns true if x is a non-zero odd integer.
func isOdd(x *Float) bool {
	return x.IsInt() && x.Sign() != 0 && x.MantExp(nil) == int(x.MinPrec())
}

// absCmpOne compares |x| to 1.
// Returns -1 if |x| < 1, 0 if |x| == 1, 1 if |x| > 1.
func absCmpOne(x *Float) int {
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
func Hypot(z, x, y *Float) *Float {
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
