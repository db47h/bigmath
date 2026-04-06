// SPDX-License-Identifier: MIT

package bigmath

import (
	"math/big"
)

// atanCore computes arctan(1/n) for n > 0 using the Taylor series
// arctan(1/n) = Σ (-1)^k / ((2k+1) * n^(2k+1))
// Using the identity: arctan(1/n) = 1/n * Σ (-1)^k / ((2k+1) * n^(2k))
func atanCore(z, n *big.Float) *big.Float {
	prec := z.Prec()
	// Guard bits to ensure precision
	workPrec := prec + 4

	// temps
	t := newFloat(workPrec)
	t2 := newFloat(workPrec)

	// term = 1/n
	term := newFloat(workPrec).Quo(one, n)

	// sum = term
	sum := newFloat(workPrec).Set(term)

	// v = 1/n^2
	v := newFloat(workPrec).Mul(n, n)
	v.Quo(one, v)

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
			t.Sub(sum, term)
		} else {
			t.Add(sum, term)
		}
		sum, t = t, sum
	}
	return z.Set(sum)
}

// Atan sets z to the rounded value of arctan(x) and returns z.
//
// Special cases:
//
//	Atan(±0) = ±0
//	Atan(±Inf) = ±π/2
func Atan(z, x *big.Float) *big.Float {
	if x.IsInf() {
		Pi(z)
		z.Quo(z, two)
		if x.Signbit() {
			z.Neg(z)
		}
		return z
	}
	if x.Sign() == 0 {
		return z.Set(zero)
	}

	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	prec += 4

	x = new(big.Float).SetPrec(prec).Copy(x)
	var neg bool
	if x.Signbit() {
		neg = true
		x.Neg(x)
	}

	// Reduction
	nReductions := 0
	for x.Cmp(new(big.Float).SetFloat64(0.1)) > 0 {
		// x = x / (1 + sqrt(1+x^2))
		t := newFloat(prec).Mul(x, x)
		t.Add(t, one)
		t.Sqrt(t)
		t.Add(t, one)
		x.Quo(x, t)
		nReductions++
	}

	// Now x is small, use atanCore(1/n) where n = 1/x
	n := newFloat(prec).Quo(one, x)
	atanCore(z, n)

	// Undo the double angle reductions
	z.SetMantExp(z, nReductions)

	if neg {
		z.Neg(z)
	}
	return z
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

// Pi sets z to the rounded value of PI and returns z.
func Pi(z *big.Float) *big.Float {
	prec := z.Prec()
	if prec == 0 {
		prec = 53
	}
	return z.Set(pi(prec))
}
