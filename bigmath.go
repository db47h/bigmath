// SPDX-License-Identifier: MIT

// Package bigmath provides arbitrary precision mathematical functions for [big.Float].
// The functions in this package follow the same rounding and precision semantics
// as the standard [big.Float] operations.
package bigmath

import (
	"math/big"
)

// An ErrNaN panic is raised by a [big.Float] operation that would lead to
// a NaN under IEEE 754 rules. An ErrNaN implements the error interface.
type ErrNaN string

// TODO: We cannot reuse big.ErrNaN (since its implementation is private) nor
// easily embed it without using panic/recover. Also, the bigmath code should
// not trigger any big.ErrNaN. Need to review the code to make sure this does
// not happen. The possible sources for big.ErrNaN are:
//
// - NewFloat panics with [ErrNaN] if x is a NaN.
// - SetFloat64 panics with [ErrNaN] if x is a NaN.
// - Add panics with [ErrNaN] if x and y are infinities with opposite signs
// - Sub panics with [ErrNaN] if x and y are infinities with equal signs
// - Mul panics with [ErrNaN] if one operand is zero and the other is infinity
// - Quo panics with [ErrNaN] if both operands are zero or infinities

func (err ErrNaN) Error() string {
	return string(err)
}

func newFloat(prec uint) *big.Float {
	return new(big.Float).SetPrec(prec)
}

func ULPExponent(x *big.Float) int {
	return x.MantExp(nil) - int(x.Prec())
}

// fma sets z to x * y + t and returns z.
// z may be an alias of x or y without causing any extra memory allocations.
// temp is a scratch variable that will hold the temp result of the multiplication with added precision.
// temp should not be an alias of any other argument.
func fma(z, x, y, t, temp *big.Float) *big.Float {
	// Use full precision for the product. Mul computes the product with full
	// precision before rounding. As a result, setting temp's precision to
	// x.prec + z.prec does not cause any extra allocations, even if
	// x.MinPrec() < x.Prec().
	// Since z's precision may change and z could be an alias for x or y, set
	//  temp's precision early.
	temp.SetPrec(0).SetPrec(x.Prec() + y.Prec())

	if z.Prec() == 0 {
		z.SetPrec(max(x.Prec(), y.Prec(), t.Prec()))
	}
	return z.Add(temp.Mul(x, y), t)
}

// FMA sets z to x * y + t and returns z.
// The operation is performed with extra precision to minimize rounding errors.
func FMA(z, x, y, t *big.Float) *big.Float {
	return fma(z, x, y, t, new(big.Float))
}
