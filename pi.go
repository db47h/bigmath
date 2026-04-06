// SPDX-License-Identifier: MIT

package bigmath

import (
	"math/big"
)

// atan1over computes arctan(1/n) for n > 0 using the Taylor series
// arctan(1/n) = Σ (-1)^k / ((2k+1) * n^(2k+1))
// Using the identity: arctan(1/n) = 1/n * Σ (-1)^k / ((2k+1) * n^(2k))
func atan1over(n *big.Float, prec uint) *big.Float {
	// Guard bits to ensure precision
	workPrec := addPrec(prec, 2)

	// temps
	t := newFloat(workPrec)
	t2 := newFloat(workPrec)

	// term = 1/n
	term := newFloat(workPrec).Quo(one, n)

	// sum = term
	sum := newFloat(workPrec).Set(term)

	// v = 1/n^2
	v := newFloat(workPrec).Quo(one, t.Mul(n, n))

	for i := uint64(1); ; i++ {
		// term = term * v * (2i-1) / (2i+1)
		t.Mul(term, v)

		t2.SetUint64(2*i + 1)
		term.Quo(t, t2)

		t2.SetUint64(2*i - 1)
		t.Mul(term, t2)

		term, t = t, term

		if term.Sign() == 0 || term.MantExp(nil) < ULPExponent(sum) {
			break
		}

		if i%2 != 0 {
			sum.Sub(sum, term)
		} else {
			sum.Add(sum, term)
		}
	}
	return sum
}

// computePi computes PI using Machin's formula: PI/4 = 4*arctan(1/5) - arctan(1/239)
func computePi(prec uint) *big.Float {
	workPrec := addPrec(prec, 4)

	t := newFloat(workPrec)

	// 4*arctan(1/5)
	p1 := atan1over(five, workPrec)
	t.Mul(p1, four)

	// arctan(1/239)
	p2 := atan1over(twoHundredThirtyNine, workPrec)

	// pi/4 = 4*arctan(1/5) - arctan(1/239)
	p1.Sub(t, p2)

	// pi = 4 * (pi/4)
	return t.Mul(p1, four).SetPrec(prec)
}

// Pi sets z to the rounded value of PI and returns z.
func Pi(z *big.Float) *big.Float {
	prec := z.Prec()
	if prec == 0 {
		prec = 53
	}
	return z.Set(computePi(prec))
}
