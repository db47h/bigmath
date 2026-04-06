// SPDX-License-Identifier: MIT

package bigmath

import (
	"math/big"
)

// atan1over computes arctan(1/n) for n > 0 using the Taylor series
// arctan(1/n) = sum_{k=0}^{\infty} (-1)^k / ((2k+1) * n^{2k+1})
// Using the identity: arctan(1/n) = 1/n * sum_{k=0}^{\infty} (-1)^k / ((2k+1) * n^{2k})
func atan1over(n uint64, prec uint) *big.Float {
	// Guard bits to ensure precision
	workPrec := addPrec(prec, 2)
	
	n2 := newFloat(workPrec).SetUint64(n * n)
	
	// term = 1/n
	term := newFloat(workPrec).SetUint64(1)
	term.Quo(term, newFloat(workPrec).SetUint64(n))
	
	// sum = term
	sum := newFloat(workPrec).Set(term)
	
	// v = 1/n^2
	v := newFloat(workPrec).SetUint64(1)
	v.Quo(v, n2)
	
	for i := uint64(1); ; i++ {
		// term = term * v * (2i-1) / (2i+1)
		term.Mul(term, v)
		
		t := newFloat(workPrec).SetUint64(2*i + 1)
		term.Quo(term, t)
		
		// term *= (2i-1)
		term.Mul(term, newFloat(workPrec).SetUint64(2*i-1))
		
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
	
	// 4*arctan(1/5)
	p1 := atan1over(5, workPrec)
	p1.Mul(p1, newFloat(workPrec).SetUint64(4))
	
	// arctan(1/239)
	p2 := atan1over(239, workPrec)
	
	// pi/4 = 4*arctan(1/5) - arctan(1/239)
	res := p1.Sub(p1, p2)
	
	// pi = 4 * (pi/4)
	res.Mul(res, newFloat(workPrec).SetUint64(4))
	
	return res.SetPrec(prec)
}

// Pi sets z to the rounded value of PI and returns z.
func Pi(z *big.Float) *big.Float {
	prec := z.Prec()
	if prec == 0 {
		prec = 53
	}
	return z.Set(computePi(prec))
}
