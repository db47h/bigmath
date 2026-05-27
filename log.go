// SPDX-License-Identifier: MIT

package bigmath

// computeLn calculates z = ln(x) using the artanh series using z's precision
// and returns z. Optimized for constants where x is a small integer near 1. The
// result is not guaranteed to be correct in the last word of the mantissa.
func computeLn(z, x *Float) *Float {
	prec := z.Prec()
	// target = (x-1)/(x+1)
	t0 := newFloat(prec).Sub(x, one)
	t1 := newFloat(prec).Add(x, one)
	term := newFloat(prec).Quo(t0, t1)
	z.Set(term)
	v2 := newFloat(prec).Mul(term, term)

	for i := uint64(1); ; i++ {
		// term = term * v^2
		t0.Mul(term, v2)
		term, t0 = t0, term

		// t0 = term / (2i + 1)
		t0.Quo(term, t1.SetUint64(2*i+1))

		if t0.Sign() == 0 || t0.MantExp(nil) < ULPExponent(z) {
			break
		}
		t0.Add(z, t0)
		z, t0 = t0, z
	}

	return z.SetMantExp(z, 1)
}

func Log(z, x *Float) *Float {
	if x.Sign() <= 0 {
		if x.Sign() == 0 {
			return z.SetInf(true) // ln(0) = -Inf
		}
		panic(ErrNaN("logarithm of negative number"))
	}
	if x.IsInf() {
		return z.Set(x)
	}

	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	// Guard bits for intermediate calculations
	prec += _W

	// 1. Primary Reduction: x = m * 2^exp
	m := new(Float).Copy(x)
	exp := m.MantExp(m) // m is now [0.5, 1)

	// 2. Secondary Reduction: Center m around 1.0
	// If m < sqrt(0.5), it's closer to 0.5. Shift it to [0.707, 1.414]
	if m.Cmp(sqrt2(prec)) < 0 {
		m.SetMantExp(m, 1)
		exp--
	}

	// 3. Compute ln(m) using the artanh series
	// ln(m) = 2 * artanh((m-1)/(m+1))
	lnM := computeLn(newFloat(prec), m)

	// 4. Combine: ln(x) = ln(m) + exp × ln(2)
	if exp != 0 {
		// FMA is cheap here since the internal precision will be prec+64
		// and it will handle temps nicely.
		return FMA(z, m.SetPrec(0).SetInt64(int64(exp)), ln2(prec), lnM)
	}

	return z.Set(lnM)
}
