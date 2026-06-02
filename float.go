// SPDX-License-Identifier: MIT

package bigmath

import (
	"math/big"
)

//go:generate go run internal/gen/main.go
// NOTE: go generate requires Go 1.26+ (the generator uses go/types APIs
// added in 1.26). The generated float_gen.go is checked in, so end users
// at any Go >= 1.22 are unaffected.

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
func (x *Float) ULPExponent() int {
	return x.MantExp(nil) - int(x.Prec())
}

// isOdd returns true if x is a non-zero odd integer.
func (x *Float) isOdd() bool {
	return x.IsInt() && x.Sign() != 0 && x.MantExp(nil) == int(x.MinPrec())
}

// absCmpOne compares |x| to 1.
// Returns -1 if |x| < 1, 0 if |x| == 1, 1 if |x| > 1.
func (x *Float) absCmpOne() int {
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
func (z *Float) Hypot(x, y *Float) *Float {
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
	t := newFloat(workPrec).FMA(x, x, newFloat(2*y.Prec()).Mul(y, y))

	return z.Sqrt(t)
}

// fma implements fused multiply-add: z = x*y + t with a single rounding.
//
// big.Float.Mul always computes the full mantissa product (O(n×m) words)
// before rounding to the target precision. By sizing temp at
// x.Prec() + y.Prec(), we capture the full product without intermediate
// rounding. Then z.Add rounds once to z.Prec(). Result: one rounding
// for the entire x*y + t expression — genuine FMA semantics.
//
// The sum-of-precisions temp size is NOT wasteful. The full product would
// have been computed internally by Mul regardless; sizing temp this way
// just preserves it. Shrinking temp would introduce intermediate rounding
// and break the FMA guarantee.
//
// z may alias x or y without extra allocations.
// temp receives the full-precision product and must NOT alias any argument.
func (z *Float) fma(x, y, t, temp *Float) *Float {
	// Size temp to hold the full product: Mul computes the full mantissa
	// product internally, so setting temp's precision to x.Prec() + y.Prec()
	// prevents it from rounding the product away. SetPrec(0) first to
	// free any previous mantissa, ensuring a fresh allocation of the
	// correct size.
	temp.SetPrec(0).SetPrec(x.Prec() + y.Prec())

	if z.Prec() == 0 {
		z.SetPrec(max(x.Prec(), y.Prec(), t.Prec()))
	}
	return z.Add(temp.Mul(x, y), t)
}

// fma implements fused multiply-sub: z = x*y - t with a single rounding.
// See fma.
func (z *Float) fms(x, y, t, temp *Float) *Float {
	temp.SetPrec(0).SetPrec(x.Prec() + y.Prec())

	if z.Prec() == 0 {
		z.SetPrec(max(x.Prec(), y.Prec(), t.Prec()))
	}
	return z.Sub(temp.Mul(x, y), t)
}

// FMA sets z to x*y + t with a single rounding (fused multiply-add) and
// returns z. The product x*y is computed without intermediate rounding,
// then added to t and rounded once to z's precision.
func (z *Float) FMA(x, y, t *Float) *Float {
	return z.fma(x, y, t, new(Float))
}

// FMS sets z to x*y - t with a single rounding (fused multiply-sub) and
// returns z. The product x*y is computed without intermediate rounding,
// then added to t and rounded once to z's precision.
func (z *Float) FMS(x, y, t *Float) *Float {
	return z.fms(x, y, t, new(Float))
}

// Inv sets z to the rounded quotient 1/y and returns z.
// Precision, rounding, and accuracy reporting are as for [Float.Add].
func (z *Float) Inv(x *Float) *Float {
	return z.Quo(one, x)
}
