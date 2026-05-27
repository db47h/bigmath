// SPDX-License-Identifier: MIT

// Package bigmath builds upon [big.Float] and provides [Float], a drop-in
// replacement with a complete set of methods for transcendentals, power,
// trigonometric and hyperbolic functions, and a [Complex] with full arithmetic
// (add, sub, mul, quo), transcendental and trigonometric support.
//
// The functions in this package follow the same rounding and precision semantics
// as the standard [big.Float] operations. Use them exactly like [big.Float]:
//
//	var z bigmath.Float
//	z.SetPrec(128)
//	bigmath.Sin(&z, &z)
package bigmath

// An ErrNaN panic is raised by a [Float] operation that would lead to
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
