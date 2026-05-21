// SPDX-License-Identifier: MIT

// TODO: review workPrec heuristic. Currently uses z.Prec() + 2*_W for guard
// bits in sinCore/cosCore/sincosCore. The other Taylor-series functions in
// this package use proportional/profiled heuristics: exp.go scales with
// prec*0.15 + 2*_W, atan.go uses atanExtraBits(x, prec), log.go uses +_W.
// A similar approach would improve the precision-to-performance trade-off.

package bigmath

import (
	"math/big"
)

// reducePi2 reduces x modulo 2π and maps the result to [0, π/2).
// It returns the quadrant (0–3) of the original reduced value.
// r may alias x. r's precision determines the target precision of the result.
func reducePi2(r, x *big.Float) int {
	if x.Sign() == 0 {
		return 0
	}

	targetPrec := r.Prec()

	xAbs := newFloat(targetPrec).Set(x)
	if xAbs.Signbit() {
		xAbs.Neg(xAbs)
	}

	xExp := xAbs.MantExp(nil)
	guard := uint(0)
	if xExp > 4 {
		guard = uint(xExp-2) + _W
	}
	workPrec := targetPrec + guard

	twoPi := newFloat(workPrec).Set(pi(workPrec))
	twoPi.SetMantExp(twoPi, 1)

	nFloat := newFloat(workPrec).Quo(xAbs, twoPi)
	nInt := new(big.Int)
	nFloat.Int(nInt)

	nFloat.SetInt(nInt)
	rTmp := newFloat(workPrec).Sub(xAbs, nFloat.Mul(nFloat, twoPi))

	if rTmp.Sign() < 0 {
		rTmp.Add(rTmp, twoPi)
	} else if rTmp.Cmp(twoPi) == 0 {
		rTmp.Set(zero)
	}

	pWork := newFloat(workPrec).Set(pi(workPrec))
	halfPi := newFloat(workPrec).SetMantExp(pWork, -1)
	pVal := newFloat(workPrec).SetMantExp(halfPi, 1)

	quad := 0
	switch {
	case rTmp.Cmp(halfPi) < 0:
	case rTmp.Cmp(pVal) < 0:
		quad = 1
		rTmp.Sub(pVal, rTmp)
	default:
		threeHalfPi := newFloat(workPrec).Add(pVal, halfPi)
		if rTmp.Cmp(threeHalfPi) < 0 {
			quad = 2
			rTmp.Sub(rTmp, pVal)
		} else {
			quad = 3
			twoPiHF := newFloat(workPrec).SetMantExp(pVal, 1)
			rTmp.Sub(twoPiHF, rTmp)
		}
	}

	r.Set(rTmp)
	return quad
}

// sinCore computes sin(x) for x in [0, π/2] using the Taylor series.
// sin(x) = Σ (−1)ⁿ · x^(2n+1) / (2n+1)!
func sinCore(z, x *big.Float) *big.Float {
	prec := z.Prec()
	workPrec := prec + 2*_W

	v := newFloat(workPrec).Mul(x, x)
	term := newFloat(workPrec).Set(x)
	sum := newFloat(workPrec).Set(x)

	t0 := newFloat(workPrec)
	t1 := newFloat(workPrec)

	for i := uint64(1); ; i++ {
		t0.Mul(term, v)
		t1.SetUint64(2*i * (2*i + 1))
		term.Quo(t0, t1)

		if term.Sign() == 0 || term.MantExp(nil) < ULPExponent(sum) {
			break
		}

		if i%2 != 0 {
			t0.Sub(sum, term)
		} else {
			t0.Add(sum, term)
		}
		sum, t0 = t0, sum
	}
	return z.Set(sum)
}

// cosCore computes cos(x) for x in [0, π/2] using the Taylor series.
// cos(x) = Σ (−1)ⁿ · x^(2n) / (2n)!
func cosCore(z, x *big.Float) *big.Float {
	prec := z.Prec()
	workPrec := prec + 2*_W

	v := newFloat(workPrec).Mul(x, x)
	term := newFloat(workPrec).Set(one)
	sum := newFloat(workPrec).Set(one)

	t0 := newFloat(workPrec)
	t1 := newFloat(workPrec)

	for i := uint64(1); ; i++ {
		t0.Mul(term, v)
		t1.SetUint64((2*i - 1) * (2 * i))
		term.Quo(t0, t1)

		if term.Sign() == 0 || term.MantExp(nil) < ULPExponent(sum) {
			break
		}

		if i%2 != 0 {
			t0.Sub(sum, term)
		} else {
			t0.Add(sum, term)
		}
		sum, t0 = t0, sum
	}
	return z.Set(sum)
}

// sincosCore computes both sin(x) and cos(x) for x in [0, π/2]
// using a single Taylor series loop that shares the computation of x²
// and the factorial denominator between both series.
func sincosCore(zs, zc, x *big.Float) (*big.Float, *big.Float) {
	prec := zs.Prec()
	workPrec := prec + 2*_W

	v := newFloat(workPrec).Mul(x, x)
	sinTerm := newFloat(workPrec).Set(x)
	cosTerm := newFloat(workPrec).Set(one)
	sinSum := newFloat(workPrec).Set(x)
	cosSum := newFloat(workPrec).Set(one)

	t0 := newFloat(workPrec)
	t1 := newFloat(workPrec)

	for i := uint64(1); ; i++ {
		t0.Mul(cosTerm, v)
		t1.SetUint64((2*i - 1) * (2 * i))
		cosTerm.Quo(t0, t1)

		t0.Mul(sinTerm, v)
		t1.SetUint64(2 * i * (2*i + 1))
		sinTerm.Quo(t0, t1)

		sinDone := sinTerm.Sign() == 0 || sinTerm.MantExp(nil) < ULPExponent(sinSum)
		cosDone := cosTerm.Sign() == 0 || cosTerm.MantExp(nil) < ULPExponent(cosSum)
		if sinDone && cosDone {
			break
		}

		if i%2 != 0 {
			t0.Sub(sinSum, sinTerm)
			t1.Sub(cosSum, cosTerm)
		} else {
			t0.Add(sinSum, sinTerm)
			t1.Add(cosSum, cosTerm)
		}
		sinSum, t0 = t0, sinSum
		cosSum, t1 = t1, cosSum
	}

	// Handle zs == zc aliasing: copy to temps before Set
	s := newFloat(prec).Set(sinSum)
	c := newFloat(prec).Set(cosSum)
	zs.Set(s)
	zc.Set(c)
	return zs, zc
}

// Sin sets z to the sine of x and returns z.
//
// Special cases:
//
//	Sin(±0) = ±0
//	Sin(±Inf) = panic(ErrNaN)
func Sin(z, x *big.Float) *big.Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.IsInf() {
		panic(ErrNaN("sin of infinity"))
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

	r := newFloat(workPrec)
	quad := reducePi2(r, xVal)

	sinCore(z, r)

	if quad >= 2 {
		z.Neg(z)
	}
	if neg {
		z.Neg(z)
	}

	return z
}

// Cos sets z to the cosine of x and returns z.
//
// Special cases:
//
//	Cos(±0) = 1
//	Cos(±Inf) = panic(ErrNaN)
func Cos(z, x *big.Float) *big.Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.IsInf() {
		panic(ErrNaN("cos of infinity"))
	}
	if x.Sign() == 0 {
		return z.Set(one)
	}

	workPrec := prec + 2*_W

	xVal := newFloat(workPrec).Set(x)
	if xVal.Signbit() {
		xVal.Neg(xVal)
	}

	r := newFloat(workPrec)
	quad := reducePi2(r, xVal)

	cosCore(z, r)

	if quad == 1 || quad == 2 {
		z.Neg(z)
	}

	return z
}

// Sincos sets zs to sin(x) and zc to cos(x) and returns both.
// The precision of zs determines the working precision; if zs's precision is 0,
// it is set from x's precision. zc's precision is set to match zs.
//
// Special cases:
//
//	Sincos(±0, zc) = ±0, 1
//	Sincos(±Inf, zc) = panic(ErrNaN)
func Sincos(zs, zc, x *big.Float) (*big.Float, *big.Float) {
	prec := zs.Prec()
	if prec == 0 {
		prec = x.Prec()
		zs.SetPrec(prec)
	}
	zc.SetPrec(prec)

	if x.IsInf() {
		panic(ErrNaN("sincos of infinity"))
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

	r := newFloat(workPrec)
	quad := reducePi2(r, xVal)

	sincosCore(zs, zc, r)

	if quad >= 2 {
		zs.Neg(zs)
	}
	if quad == 1 || quad == 2 {
		zc.Neg(zc)
	}
	if neg {
		zs.Neg(zs)
	}

	return zs, zc
}
