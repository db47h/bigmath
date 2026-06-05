// SPDX-License-Identifier: MIT

package bigmath

// lnCore calculates z = ln(x) using the artanh series using z's precision
// and returns z. Optimized for constants where x is a small integer near 1. The
// result is not guaranteed to be correct in the last word of the mantissa.
func (z *Float) lnCore(x *Float) *Float {
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

		if t0.Sign() == 0 || t0.MantExp(nil) < z.ULPExponent() {
			break
		}
		t0.Add(z, t0)
		z, t0 = t0, z
	}

	return z.SetMantExp(z, 1)
}

func (z *Float) Log(x *Float) *Float {
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
	if m.Cmp(sqrt2.get(prec)) < 0 {
		m.SetMantExp(m, 1)
		exp--
	}

	// 3. Compute ln(m) using the artanh series
	// ln(m) = 2 * artanh((m-1)/(m+1))
	lnM := newFloat(prec).lnCore(m)

	// 4. Combine: ln(x) = ln(m) + exp × ln(2)
	if exp != 0 {
		// FMA is cheap here since the internal precision will be prec+64
		// and it will handle temps nicely.
		return z.FMA(m.SetPrec(0).SetInt64(int64(exp)), ln2.get(prec), lnM)
	}

	return z.Set(lnM)
}

// Log10 sets z to the base-10 logarithm of x and returns z.
func (z *Float) Log10(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}
	// Log10(x) = Log(x) / ln(10)
	workPrec := prec + _W
	lnX := newFloat(workPrec).Log(x)
	return z.Quo(lnX, ln10.get(workPrec))
}

// Log2 sets z to the base-2 logarithm of x and returns z.
func (z *Float) Log2(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}
	// Log2(x) = Log(x) / ln(2)
	workPrec := prec + _W
	lnX := newFloat(workPrec).Log(x)
	return z.Quo(lnX, ln2.get(workPrec))
}
