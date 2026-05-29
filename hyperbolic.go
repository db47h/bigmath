// SPDX-License-Identifier: MIT

package bigmath

// sinhCore computes sinh(x) for x in [0, 1] using the Taylor series.
// sinh(x) = x + x³/3! + x⁵/5! + x⁷/7! + ...
//
// workPrec uses a flat +2*_W guard. Unlike sinCore, this series has no
// alternating signs (all terms are positive), so there is no subtractive
// cancellation. A smaller guard (+_W) would likely suffice, but +2*_W is
// harmless — the performance difference is negligible at any precision.
func (z *Float) sinhCore(x *Float) *Float {
	prec := z.Prec()
	workPrec := prec + 2*_W

	x2 := newFloat(workPrec).Mul(x, x)
	term := newFloat(workPrec).Set(x)
	sum := newFloat(workPrec).Set(x)

	t0 := newFloat(workPrec)
	t1 := newFloat(workPrec)

	for i := uint64(1); ; i++ {
		t0.Mul(term, x2)
		t1.SetUint64(2 * i * (2*i + 1))
		term.Quo(t0, t1)

		if term.Sign() == 0 || term.MantExp(nil) < sum.ULPExponent() {
			break
		}

		t0.Add(sum, term)
		sum, t0 = t0, sum
	}
	return z.Set(sum)
}

// sinhcoshCore computes both sinh(x) and cosh(x) for x in [0, 1]
// using a single Taylor series loop that shares the computation of x²
// and the factorial denominator between both series.
//
// Same flat +2*_W guard as sinhCore. Both series are all-positive (no
// alternating signs), so the guard is conservative but adequate.
func sinhcoshCore(zs, zc, x *Float) (*Float, *Float) {
	prec := zs.Prec()
	workPrec := prec + 2*_W

	x2 := newFloat(workPrec).Mul(x, x)
	sinTerm := newFloat(workPrec).Set(x)
	cosTerm := newFloat(workPrec).Set(one)
	sinSum := newFloat(workPrec).Set(x)
	cosSum := newFloat(workPrec).Set(one)

	t0 := newFloat(workPrec)
	t1 := newFloat(workPrec)

	for i := uint64(1); ; i++ {
		t0.Mul(cosTerm, x2)
		t1.SetUint64((2*i - 1) * (2 * i))
		cosTerm.Quo(t0, t1)

		t0.Mul(sinTerm, x2)
		t1.SetUint64(2 * i * (2*i + 1))
		sinTerm.Quo(t0, t1)

		sinDone := sinTerm.Sign() == 0 || sinTerm.MantExp(nil) < sinSum.ULPExponent()
		cosDone := cosTerm.Sign() == 0 || cosTerm.MantExp(nil) < cosSum.ULPExponent()
		if sinDone && cosDone {
			break
		}

		t0.Add(sinSum, sinTerm)
		t1.Add(cosSum, cosTerm)
		sinSum, t0 = t0, sinSum
		cosSum, t1 = t1, cosSum
	}

	s := newFloat(prec).Set(sinSum)
	c := newFloat(prec).Set(cosSum)
	zs.Set(s)
	zc.Set(c)
	return zs, zc
}

// Sinh sets z to the hyperbolic sine of x and returns z.
//
// Special cases:
//
//	Sinh(±0) = ±0
//	Sinh(±Inf) = ±Inf
func (z *Float) Sinh(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.IsInf() {
		return z.SetInf(x.Signbit())
	}
	if x.Sign() == 0 {
		return z.Set(x)
	}

	// Flat +_W guard for the full computation path. For |x|<1 the Taylor
	// series sinhCore uses +2*_W internally; for |x|≥1 the Exp call brings
	// its own proportional guard. The outer guard covers sign and branch
	// handling around both paths.
	workPrec := prec + _W

	t0 := new(Float).Abs(x)
	neg := x.Signbit()

	if t0.Cmp(one) < 0 {
		z.sinhCore(t0)
		if neg {
			z.Neg(z)
		}
		return z
	}

	// |x| >= 1: use (eˣ − e⁻ˣ) / 2
	t1 := newFloat(workPrec).Exp(t0)
	t0.SetPrec(0).SetPrec(workPrec).Inv(t1)

	z.Sub(t1, t0)
	z.SetMantExp(z, -1)
	if neg {
		z.Neg(z)
	}
	return z
}

// Cosh sets z to the hyperbolic cosine of x and returns z.
//
// Special cases:
//
//	Cosh(±0) = 1
//	Cosh(±Inf) = +Inf
func (z *Float) Cosh(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.IsInf() {
		return z.SetInf(false)
	}
	if x.Sign() == 0 {
		return z.Set(one)
	}

	// Flat +_W guard. Always uses Exp internally (no Taylor path), so
	// Exp's own proportional guard plus this outer margin covers it.
	workPrec := prec + _W

	t0 := new(Float).Abs(x)

	t1 := newFloat(workPrec).Exp(t0)
	t0.SetPrec(0).SetPrec(workPrec).Inv(t1)

	z.Add(t1, t0)
	z.SetMantExp(z, -1)

	return z
}

// SinhCosh sets zs to sinh(x) and zc to cosh(x) and returns both.
// The precision of zs determines the working precision; if zs's precision is 0,
// it is set from x's precision. zc's precision is set to match zs.
//
// Special cases:
//
//	SinhCosh(±0, zc) = ±0, 1
//	SinhCosh(±Inf, zc) = ±Inf, +Inf
func SinhCosh(zs, zc, x *Float) (*Float, *Float) {
	prec := zs.Prec()
	if prec == 0 {
		prec = x.Prec()
		zs.SetPrec(prec)
	}
	zc.SetPrec(prec)

	if x.IsInf() {
		zs.SetInf(x.Signbit())
		zc.SetInf(false)
		return zs, zc
	}
	if x.Sign() == 0 {
		zs.Set(x)
		zc.Set(one)
		return zs, zc
	}

	// Flat +_W guard. Dispatches to sinhcoshCore (|x|<1) or shared Exp
	// (|x|≥1); in both paths the internal temps use matching precision.
	workPrec := prec + _W

	t0 := new(Float).Abs(x)
	neg := x.Signbit()

	if t0.Cmp(one) < 0 {
		sinhcoshCore(zs, zc, t0)
		if neg {
			zs.Neg(zs)
		}
		return zs, zc
	}

	// |x| >= 1: share the Exp call between sinh and cosh.
	t1 := newFloat(workPrec).Exp(t0)
	t0.SetPrec(0).SetPrec(workPrec).Inv(t1)

	zs.Sub(t1, t0)
	zs.SetMantExp(zs, -1)
	zc.Add(t1, t0)
	zc.SetMantExp(zc, -1)

	if neg {
		zs.Neg(zs)
	}
	return zs, zc
}

// Tanh sets z to the hyperbolic tangent of x and returns z.
//
// Special cases:
//
//	Tanh(±0) = ±0
//	Tanh(±Inf) = ±1
func (z *Float) Tanh(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.Sign() == 0 {
		return z.Set(x)
	}
	if x.IsInf() {
		if x.Signbit() {
			z.Set(minusOne)
		} else {
			z.Set(one)
		}
		return z
	}

	// Flat +2*_W guard. Delegates to SinhCosh internally, which in turn
	// provides sufficient precision for both sinh and cosh paths.
	workPrec := prec + _W

	s, c := SinhCosh(newFloat(workPrec), newFloat(workPrec), x)

	if s.IsInf() && c.IsInf() {
		neg := x.Signbit()
		z.Set(one)
		if neg {
			z.Neg(z)
		}
		return z
	}

	return z.Quo(s, c)
}

// Asinh sets z to the inverse hyperbolic sine of x and returns z.
//
// Special cases:
//
//	Asinh(±0) = ±0
//	Asinh(±Inf) = ±Inf
func (z *Float) Asinh(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.Sign() == 0 {
		return z.Set(x)
	}
	if x.IsInf() {
		return z.SetInf(x.Signbit())
	}

	// Flat +2*_W guard. The computation involves a sqrt, an add, and a Log
	// call (which brings its own +_W guard internally). The outer margin
	// covers intermediate rounding in sqrt and addition.
	workPrec := prec + 2*_W

	neg := x.Signbit()
	xVal := new(Float).Abs(x)

	t0 := newFloat(workPrec).Mul(xVal, xVal)
	if t0.IsInf() {
		// For extremely large x, asinh(x) ≈ ln(x) + ln(2)
		t0.Log(xVal)
		z.Add(t0, ln2(workPrec))
		if neg {
			z.Neg(z)
		}
		return z
	}

	t1 := newFloat(workPrec).Add(t0, one)
	t0.Sqrt(t1)
	t1.Add(xVal, t0)
	z.Log(t1)
	if neg {
		z.Neg(z)
	}
	return z
}

// acoshGuard computes the working precision needed for acosh(x) near x=1,
// where x-1 loses significant bits. Returns the required working precision.
func (x *Float) acoshGuard(prec uint) uint {
	xExp := x.MantExp(nil)
	if xExp > 1 {
		return prec + 2*_W
	}
	// x is in [1, 2). Compute x-1 at the default working precision to
	// determine the guard bits needed.
	workPrec := prec + 2*_W
	t := newFloat(workPrec).Sub(x, one)
	subExp := t.MantExp(nil)
	if -subExp <= 2 {
		return prec + 2*_W
	}
	need := max(prec+uint(-subExp)+3*_W, workPrec)
	return need
}

// Acosh sets z to the inverse hyperbolic cosine of x and returns z.
//
// Special cases:
//
//	Acosh(1) = 0
//	Acosh(x < 1) = panic(ErrNaN)
//	Acosh(+Inf) = +Inf
func (z *Float) Acosh(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.IsInf() {
		return z.SetInf(false)
	}
	switch x.Cmp(one) {
	case -1:
		panic(ErrNaN("acosh of x < 1"))
	case 0:
		return z.Set(zero)
	}

	workPrec := x.acoshGuard(prec)

	t0 := newFloat(workPrec).Sub(x, one)
	t1 := newFloat(workPrec).Add(x, one)

	// √(x-1) * √(x+1) for numerical stability (avoids computing x² directly).
	t0.Sqrt(t0)
	t1.Sqrt(t1)
	t0.Mul(t0, t1)
	t0.Add(x, t0)
	t0.Log(t0)

	return z.Set(t0)
}

// atanhGuard computes the working precision needed for atanh(x) near x=±1,
// where 1-x (or 1+x) loses significant bits. Returns the required precision.
func (x *Float) atanhGuard(prec uint) uint {
	// Work with |x| to always check 1-|x|.
	xAbs := newFloat(prec + 2*_W).Abs(x)
	oneMinus := newFloat(prec+2*_W).Sub(one, xAbs)
	subExp := oneMinus.MantExp(nil)
	if -subExp <= 2 {
		return prec + 2*_W
	}
	need := max(prec+uint(-subExp)+3*_W, prec+2*_W)
	return need
}

// Atanh sets z to the inverse hyperbolic tangent of x and returns z.
//
// Special cases:
//
//	Atanh(±0) = ±0
//	Atanh(1) = +Inf
//	Atanh(-1) = -Inf
//	Atanh(|x| > 1) = panic(ErrNaN)
func (z *Float) Atanh(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.Sign() == 0 {
		return z.Set(x)
	}

	// Domain check: |x| < 1
	// |x| == 1: atanh diverges to ±∞
	// |x| > 1: mathematically undefined for real atanh
	switch x.absCmpOne() {
	case 0:
		// atanh(1) = +Inf, atanh(-1) = -Inf
		return z.SetInf(x.Signbit())
	case 1:
		panic(ErrNaN("atanh of |x| > 1"))
	}

	workPrec := x.atanhGuard(prec)

	neg := x.Signbit()
	t0 := new(Float).Abs(x)

	t1 := newFloat(workPrec).Sub(one, t0)
	t2 := newFloat(workPrec).Add(one, t0)
	t0.SetPrec(0).SetPrec(workPrec).Quo(t2, t1)
	z.Log(t0)
	z.SetMantExp(z, -1)

	if neg {
		z.Neg(z)
	}
	return z
}
