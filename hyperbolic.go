// SPDX-License-Identifier: MIT

package bigmath

import (
	"math/big"
)

// sinhCore computes sinh(x) for x in [0, 1] using the Taylor series.
// sinh(x) = x + x³/3! + x⁵/5! + x⁷/7! + ...
func sinhCore(z, x *big.Float) *big.Float {
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

		if term.Sign() == 0 || term.MantExp(nil) < ULPExponent(sum) {
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
func sinhcoshCore(zs, zc, x *big.Float) (*big.Float, *big.Float) {
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

		sinDone := sinTerm.Sign() == 0 || sinTerm.MantExp(nil) < ULPExponent(sinSum)
		cosDone := cosTerm.Sign() == 0 || cosTerm.MantExp(nil) < ULPExponent(cosSum)
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
func Sinh(z, x *big.Float) *big.Float {
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

	workPrec := prec + 2*_W

	xVal := newFloat(workPrec).Set(x)
	neg := xVal.Signbit()
	if neg {
		xVal.Neg(xVal)
	}

	if xVal.Cmp(one) < 0 {
		t := sinhCore(newFloat(workPrec), xVal)
		if neg {
			t.Neg(t)
		}
		return z.Set(t)
	}

	// |x| >= 1: use (eˣ − e⁻ˣ) / 2
	ep := Exp(newFloat(workPrec), xVal)
	em := newFloat(workPrec).Quo(one, ep)

	t0 := newFloat(workPrec).Sub(ep, em)
	t0.SetMantExp(t0, -1)

	if neg {
		t0.Neg(t0)
	}
	return z.Set(t0)
}

// Cosh sets z to the hyperbolic cosine of x and returns z.
//
// Special cases:
//
//	Cosh(±0) = 1
//	Cosh(±Inf) = +Inf
func Cosh(z, x *big.Float) *big.Float {
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

	workPrec := prec + 2*_W

	xVal := newFloat(workPrec).Set(x)
	if xVal.Signbit() {
		xVal.Neg(xVal)
	}

	ep := Exp(newFloat(workPrec), xVal)
	em := newFloat(workPrec).Quo(one, ep)

	t0 := newFloat(workPrec).Add(ep, em)
	t0.SetMantExp(t0, -1)

	return z.Set(t0)
}

// SinhCosh sets zs to sinh(x) and zc to cosh(x) and returns both.
// The precision of zs determines the working precision; if zs's precision is 0,
// it is set from x's precision. zc's precision is set to match zs.
//
// Special cases:
//
//	SinhCosh(±0, zc) = ±0, 1
//	SinhCosh(±Inf, zc) = ±Inf, +Inf
func SinhCosh(zs, zc, x *big.Float) (*big.Float, *big.Float) {
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

	workPrec := prec + 2*_W

	xVal := newFloat(workPrec).Set(x)
	neg := xVal.Signbit()
	if neg {
		xVal.Neg(xVal)
	}

	if xVal.Cmp(one) < 0 {
		return sinhcoshCore(zs, zc, xVal)
	}

	// |x| >= 1: share the Exp call between sinh and cosh.
	ep := Exp(newFloat(workPrec), xVal)
	em := newFloat(workPrec).Quo(one, ep)

	s := newFloat(workPrec).Sub(ep, em)
	s.SetMantExp(s, -1)
	c := newFloat(workPrec).Add(ep, em)
	c.SetMantExp(c, -1)

	if neg {
		s.Neg(s)
	}
	zs.Set(s)
	zc.Set(c)
	return zs, zc
}

// Tanh sets z to the hyperbolic tangent of x and returns z.
//
// Special cases:
//
//	Tanh(±0) = ±0
//	Tanh(±Inf) = ±1
func Tanh(z, x *big.Float) *big.Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.Sign() == 0 {
		return z.Set(x)
	}
	if x.IsInf() {
		z.SetUint64(1)
		if x.Signbit() {
			z.Neg(z)
		}
		return z
	}

	workPrec := prec + 2*_W

	s := newFloat(workPrec)
	c := newFloat(workPrec)
	SinhCosh(s, c, x)

	if s.IsInf() && c.IsInf() {
		z.SetUint64(1)
		if x.Signbit() {
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
func Asinh(z, x *big.Float) *big.Float {
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

	workPrec := prec + 2*_W

	xVal := newFloat(workPrec).Set(x)
	neg := xVal.Signbit()
	if neg {
		xVal.Neg(xVal)
	}

	x2 := newFloat(workPrec).Mul(xVal, xVal)
	if x2.IsInf() {
		t := Log(newFloat(workPrec), xVal)
		ln2 := Log(newFloat(workPrec), newFloat(workPrec).SetUint64(2))
		t.Add(t, ln2)
		if neg {
			t.Neg(t)
		}
		return z.Set(t)
	}

	t0 := newFloat(workPrec).Add(x2, one)
	t0.Sqrt(t0)
	t0.Add(xVal, t0)
	Log(t0, t0)
	if neg {
		t0.Neg(t0)
	}
	return z.Set(t0)
}

// acoshGuard computes the working precision needed for acosh(x) near x=1,
// where x-1 loses significant bits. Returns the required working precision.
func acoshGuard(x *big.Float, prec uint) uint {
	xExp := x.MantExp(nil)
	if xExp > 1 {
		return prec + 2*_W
	}
	// x is in [1, 4). Compute x-1 at the default working precision to
	// determine the guard bits needed.
	workPrec := prec + 2*_W
	t := newFloat(workPrec).Sub(x, one)
	subExp := t.MantExp(nil)
	if -subExp <= 2 {
		return prec + 2*_W
	}
	need := prec + uint(-subExp) + 3*_W
	if need < workPrec {
		need = workPrec
	}
	return need
}

// Acosh sets z to the inverse hyperbolic cosine of x and returns z.
//
// Special cases:
//
//	Acosh(1) = 0
//	Acosh(x < 1) = panic(ErrNaN)
//	Acosh(+Inf) = +Inf
func Acosh(z, x *big.Float) *big.Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.IsInf() {
		return z.SetInf(false)
	}
	if x.Cmp(one) < 0 {
		panic(ErrNaN("acosh of x < 1"))
	}
	if x.Cmp(one) == 0 {
		return z.Set(zero)
	}

	workPrec := acoshGuard(x, prec)

	xVal := newFloat(workPrec).Set(x)
	t0 := newFloat(workPrec).Sub(xVal, one)
	t1 := newFloat(workPrec).Add(xVal, one)

	// √(x-1) * √(x+1) for numerical stability (avoids computing x² directly).
	t0.Sqrt(t0)
	t1.Sqrt(t1)
	t0.Mul(t0, t1)
	t0.Add(xVal, t0)
	Log(t0, t0)

	return z.Set(t0)
}

// atanhGuard computes the working precision needed for atanh(x) near x=±1,
// where 1-x (or 1+x) loses significant bits. Returns the required precision.
func atanhGuard(x *big.Float, prec uint) uint {
	// Work with |x| to always check 1-|x|.
	xAbs := newFloat(prec + 2*_W).Abs(x)
	oneMinus := newFloat(prec + 2*_W).Sub(one, xAbs)
	subExp := oneMinus.MantExp(nil)
	if -subExp <= 2 {
		return prec + 2*_W
	}
	need := prec + uint(-subExp) + 3*_W
	if need < prec+2*_W {
		need = prec + 2*_W
	}
	return need
}

// Atanh sets z to the inverse hyperbolic tangent of x and returns z.
//
// Special cases:
//
//	Atanh(±0) = ±0
//	Atanh(|x| = 1) = panic(ErrNaN)
//	Atanh(|x| > 1) = panic(ErrNaN)
func Atanh(z, x *big.Float) *big.Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.Sign() == 0 {
		return z.Set(x)
	}

	// Domain check: |x| < 1
	xAbs := newFloat(prec).Abs(x)
	if xAbs.Cmp(one) >= 0 {
		panic(ErrNaN("atanh of |x| >= 1"))
	}

	workPrec := atanhGuard(x, prec)

	xVal := newFloat(workPrec).Set(x)
	neg := xVal.Signbit()
	if neg {
		xVal.Neg(xVal)
	}

	t0 := newFloat(workPrec).Sub(one, xVal)
	t1 := newFloat(workPrec).Add(one, xVal)
	t0.Quo(t1, t0)
	Log(t0, t0)
	t0.SetMantExp(t0, -1)

	if neg {
		t0.Neg(t0)
	}
	return z.Set(t0)
}
