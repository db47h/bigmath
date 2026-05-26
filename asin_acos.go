// SPDX-License-Identifier: MIT

package bigmath

import (
	"math/big"
)

// asinGuard computes the working precision needed for asin(x) near x=±1,
// where 1-x² loses significant bits. Returns the required working precision.
func asinGuard(x *big.Float, prec uint) uint {
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
func Asin(z, x *big.Float) *big.Float {
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
		Pi(z)               // z = π at z's precision
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
	result := Atan(newFloat(workPrec), ratio)

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
func Acos(z, x *big.Float) *big.Float {
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
			return Pi(z)
		}
		return z.Set(zero)
	}

	// Acos(x) = π/2 - Asin(x)
	workPrec := asinGuard(x, prec)

	tmp := Asin(newFloat(workPrec), x)
	result := newFloat(workPrec)
	Pi(result)                    // result = π at workPrec
	result.SetMantExp(result, -1) // result = π/2, keeps workPrec
	result.Sub(result, tmp)

	return z.Set(result)
}
