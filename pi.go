// SPDX-License-Identifier: MIT

package bigmath

import (
	"math/big"
)

// Pi sets z to the rounded value of PI and returns z.
func Pi(z *big.Float) *big.Float {
	prec := z.Prec()
	if prec == 0 {
		prec = 53
	}
	return z.Set(pi(prec))
}

// computePi computes PI using Machin's formula: PI/4 = 4*arctan(1/5) - arctan(1/239)
func computePi(prec uint) *big.Float {
	workPrec := prec + 4

	t := newFloat(workPrec)
	p1 := newFloat(workPrec)
	p2 := newFloat(workPrec)

	// 4*arctan(1/5)
	atanCore(p1, five)
	p1.SetMantExp(p1, 2)

	// arctan(1/239)
	atanCore(p2, twoHundredThirtyNine)

	// pi/4 = 4*arctan(1/5) - arctan(1/239)
	t.Sub(p1, p2)

	// pi = 4 * (pi/4)
	return t.SetMantExp(t, 2).SetPrec(prec)
}
