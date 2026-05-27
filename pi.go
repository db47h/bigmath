// SPDX-License-Identifier: MIT

package bigmath

// Pi sets z to the rounded value of PI and returns z.
func (z *Float) Pi() *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = 53
	}
	return z.Set(pi(prec))
}

// computePi computes PI using Machin's formula: PI/4 = 4*arctan(1/5) - arctan(1/239)
func computePi(prec uint) *Float {
	workPrec := prec + _W

	p1 := newFloat(workPrec)
	p2 := newFloat(workPrec)

	// 4*arctan(1/5)
	// Use the optimized atanReciprocal for reciprocal integers.
	atanReciprocal(p1, 5)
	p1.SetMantExp(p1, 2)

	// arctan(1/239)
	atanReciprocal(p2, 239)

	// pi/4 = 4*arctan(1/5) - arctan(1/239)
	res := newFloat(workPrec).Sub(p1, p2)

	// pi = 4 * (pi/4)
	return res.SetMantExp(res, 2).SetPrec(prec)
}
