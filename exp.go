// SPDX-License-Identifier: MIT

package bigmath

import (
	"math"
)

// Exp sets z to the rounded value of e^x, and returns z.
//
// If z's precision is 0, it is changed to x's precision before the operation.
// Rounding is performed according to z's precision and rounding mode.
//
// The operation uses the Taylor series e^x = ∑(x^n/n!) for n ≥ 0.
func (z *Float) Exp(x *Float) *Float {
	sgn := x.Sign()
	if sgn == 0 {
		return z.Set(one)
	}

	// exp(x) is finite if 0.5 × 2^big.MinExp ≤ exp(x) < 1 × 2^big.MaxExp
	//   ⇒ log(2) × (big.MinExp-1) ≤ x < log(2) × big.MaxExp
	// While this function properly handles values of x outside of this range,
	// exit early on extreme values to prevent long running times and simplify the
	// bounds check to x.exp-1 < log2(big.MaxExp)
	exp := x.MantExp(nil)
	if x.IsInf() || exp > 31 {
		if sgn < 0 {
			return z.Set(zero)
		}
		return z.SetInf(false)
	}

	// make a modifyable copy of x. Also covers the case where z == x.
	x = new(Float).Copy(x)

	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	// The following is based on R. P. Brent, P. Zimmermann, Modern Computer
	// Arithmetic, Cambridge Monographs on Computational and Applied Mathematics
	// (No. 18), Cambridge University Press
	// https://members.loria.fr/PZimmermann/mca/pub226.html
	//
	// Argument reduction: bring x in the range [0.5, 1)×2^-k for faster
	// convergence. This also brings extreme values of x for which exp(x) is
	// 0 or +Inf into a computable range (i.e. for z=x^-k, ∑(z^n/n!) is finite).

	// For x < 0, compute exp(-x) = 1/exp(x).
	// This is to prevent alternating signs in the power series terms and avoid
	// cancellation in the summation, as well as keeping the summation in a
	// known range after argument reduction (1 <= ∑(x^n/n!) < 1+2^(-k+1)).
	var invert bool
	if x.Signbit() {
		invert = true
		x.Neg(x)
	}
	// §4.3.1 & §4.4.2
	// 1 ≤ k ≤ ceil(√MaxPrec)
	k := int(math.Ceil(math.Sqrt(float64(prec))))
	// Working precision (§4.4)
	prec += uint(math.Log(float64(prec))) + 1
	squarings := 0
	if -k < exp {
		// -46341 ≤ -k ≤ -1 < exp ≤ 31
		// 1 < exp + k ≤ 46371
		squarings = exp + k
		// 0 ≤ k-1 < exp (condition needed to undo argument reduction)
		x.SetMantExp(x, -squarings)
		// 2 bits of added precision per multiplication when undoing argument reduction.
		prec += 2 * uint(squarings)
	}

	// temp vars
	n := new(Float)
	t0 := newFloat(prec)
	t1 := newFloat(prec)
	term := newFloat(prec).Set(one)
	sum := newFloat(prec).Set(one)

	// The exit condition of the loop is: term = 0 || term.exp < ULPExponent(sum)
	// Because of argument reduction, 1 ≤ sum < 1+2^(-k+1) ⇒ sum.exp == 1
	// So we can pre-calculate ULPExponent(sum)
	minExp := 1 - int(prec)

	// term(n) = term(n-1) × x/n is faster than term(n) = x^n / n! (saves one .Mul)
	for i := uint64(1); ; i++ {
		t0.Quo(x, n.SetUint64(i))
		// term.Mul(term, t) and sum.Add(sum, term) require a temp Float for the
		// result. Manage that ourselves by using our own temps t0, t1, then swap the
		// pointers.
		t1.Mul(term, t0)
		t1, term = term, t1
		if term.Sign() == 0 || term.MantExp(nil) < minExp {
			break
		}
		t1.Add(sum, term)
		t1, sum = sum, t1
	}

	// Undo argument reduction if exp > 0
	for range squarings {
		// Prevent temp allocations using the same trick as above
		t0.Mul(sum, sum)
		t0, sum = sum, t0
	}

	if invert {
		return z.Quo(one, sum)
	}
	return z.Set(sum)
}
