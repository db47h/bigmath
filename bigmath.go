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
func fma(z, x, y, t, temp *Float) *Float {
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

// FMA sets z to x*y + t with a single rounding (fused multiply-add) and
// returns z. The product x*y is computed without intermediate rounding,
// then added to t and rounded once to z's precision.
//
// This provides genuine FMA semantics: one rounding for the entire
// expression, not two. See fma for the implementation details.
func FMA(z, x, y, t *Float) *Float {
	return fma(z, x, y, t, new(Float))
}
