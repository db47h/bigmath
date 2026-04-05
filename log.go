// SPDX-License-Identifier: MIT

package bigmath

import (
	"math/big"
)

// computeLn calculates ln(x) using the artanh series.
// Optimized for constants where x is a small integer near 1.
// prec is the desired precision for the result, however, the result is not
// guaranteed to be correct in the last word of the mantissa.
func computeLn(x *big.Float, prec uint) *big.Float {
	// target = (x-1)/(x+1)
	t0 := newFloat(prec).Sub(x, one)
	t1 := newFloat(prec).Add(x, one)
	term := newFloat(prec).Quo(t0, t1)
	sum := newFloat(prec).Set(term)
	v2 := newFloat(prec).Mul(term, term)

	for i := uint64(1); ; i++ {
		// term = term * v^2
		t0.Mul(term, v2)
		term, t0 = t0, term

		// t0 = term / (2i + 1)
		t0.Quo(term, t1.SetUint64(2*i+1))

		if t0.Sign() == 0 || t0.MantExp(nil) < ULPExponent(sum) {
			break
		}
		t0.Add(sum, t0)
		sum, t0 = t0, sum
	}

	return sum.SetMantExp(sum, 1)
}

func Log(z, x *big.Float) *big.Float {
	if x.Sign() <= 0 {
		if x.Sign() == 0 {
			return z.SetInf(true).Neg(z) // ln(0) = -Inf
		}
		return z.SetInf(false) // ln(neg) = NaN (big.Float uses Inf for simplicity or you can handle differently)
	}
	if x.IsInf() {
		return z.Set(x)
	}

	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
	}
	z.SetPrec(0).SetPrec(prec)

	// Guard bits for intermediate calculations
	prec = addPrec(prec, 64)

	// 1. Primary Reduction: x = m * 2^exp
	m := new(big.Float).Copy(x)
	exp := m.MantExp(m) // m is now [0.5, 1)

	// 2. Secondary Reduction: Center m around 1.0
	// If m < sqrt(0.5), it's closer to 0.5. Shift it to [0.707, 1.414]
	if m.Cmp(sqrt2(prec)) < 0 {
		m.SetMantExp(m, 1)
		exp--
	}

	// 3. Compute ln(m) using the artanh series
	// ln(m) = 2 * artanh((m-1)/(m+1))
	res := computeLn(m, prec)

	// 4. Combine: res = ln(m) + exp × ln(2)
	if exp != 0 {
		m.SetPrec(0).SetInt64(int64(exp))
		termExp := newFloat(prec).Mul(m, ln2(prec))
		res.Add(res, termExp)
	}

	return z.Set(res)
}
