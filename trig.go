// SPDX-License-Identifier: MIT

// Guard-bit strategy: core Taylor-series loops (sinCore, cosCore, sincosCore)
// use a flat +2*_W guard. This is adequate because the accumulated rounding
// error over N iterations is bounded by N × 2^-(prec+2*_W), which stays well
// below the 0.5-ULP rounding threshold at target precision for any practical
// iteration count. (The number of terms needed for convergence at precision P
// is O(P), and 2*_W = 128 on 64-bit — far more than log₂(max_terms).)
// The outer Sin/Cos/Sincos functions apply the same flat guard for the
// entire computation path (including argument reduction and sign handling).
// This is conservative but harmless — the important adaptive guard is inside
// reducePi2, which scales with input magnitude for large arguments.
// See docs/trig-hyperbolic-precision-review.md for the full analysis.

package bigmath

import (
	"math/big"
)

// reducePi2 reduces x modulo 2π and maps the result to [0, π/2).
// It returns the quadrant (0–3) of the original reduced value.
// z may alias x. z's precision determines the target precision of the result.
func reducePi2(z, x *big.Float) int {
	if x.Sign() == 0 {
		return 0
	}

	prec := z.Prec()

	xAbs := newFloat(prec).Set(x)
	if xAbs.Signbit() {
		xAbs.Neg(xAbs)
	}

	xExp := xAbs.MantExp(nil)
	workPrec := prec + _W
	if xExp > 0 {
		workPrec += uint(xExp)
	}

	twoPi := newFloat(workPrec).Set(pi(workPrec))
	twoPi.SetMantExp(twoPi, 1)

	t0 := newFloat(workPrec).Quo(xAbs, twoPi)
	nInt := new(big.Int)
	t0.Int(nInt)

	t0.SetInt(nInt)
	t1 := newFloat(workPrec).Mul(t0, twoPi)
	rTmp := newFloat(workPrec).Sub(xAbs, t1)

	if rTmp.Sign() < 0 {
		t0.Add(rTmp, twoPi)
		t0, rTmp = rTmp, t0
	} else if rTmp.Cmp(twoPi) == 0 {
		rTmp.Set(zero)
	}

	pVal := pi(workPrec)
	// t1 is free after Sub(xAbs, t1) above; reuse as π/2.
	// In the default/else branch below, t1 is overwritten to 2π.
	t1.SetMantExp(pVal, -1)

	quad := 0
	switch {
	case rTmp.Cmp(t1) < 0:
	case rTmp.Cmp(pVal) < 0:
		quad = 1
		t0.Sub(pVal, rTmp)
		t0, rTmp = rTmp, t0
	default:
		t0.Add(pVal, t1) // pVal + π/2 = 3π/2
		if rTmp.Cmp(t0) < 0 {
			quad = 2
			t0.Sub(rTmp, pVal)
			t0, rTmp = rTmp, t0
		} else {
			quad = 3
			t1.SetMantExp(pVal, 1)
			t0.Sub(t1, rTmp)
			t0, rTmp = rTmp, t0
		}
	}

	z.Set(rTmp)
	return quad
}

// sinCore computes sin(x) for x in [0, π/2] using the Taylor series.
// sin(x) = Σ (−1)ⁿ · x^(2n+1) / (2n+1)!
//
// workPrec uses a flat +2*_W guard. The alternating series introduces
// subtractive cancellation as terms approach the sum, so +2*_W provides
// a comfortable margin. At target precision P, O(P) terms are needed;
// the accumulated error (N × 2^-(P+2*_W)) is negligible for any N reachable
// in practice (see file header for the full rationale).
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
		t1.SetUint64(2 * i * (2*i + 1))
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
//
// Same flat +2*_W guard as sinCore — same alternating-series cancellation
// characteristics. See sinCore doc for the rationale.
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
//
// Same flat +2*_W guard as sinCore/cosCore. The shared loop handles both
// alternating series; the guard covers the combined rounding error.
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

	// Flat +_W guard for the full computation path:
	// xVal at this precision feeds into reducePi2 (which adds its own
	// dynamic guard internally) and then into sinCore (whose temps are
	// also at prec+2*_W). The guard covers sign handling, reduction
	// output rounding, and the Taylor series accumulation.
	workPrec := prec + _W

	xVal := newFloat(workPrec).Set(x)
	neg := xVal.Signbit()
	if neg {
		xVal.Neg(xVal)
	}

	// reducePi2 handles z == x aliasing; reuse xVal as both input and output.
	quad := reducePi2(xVal, xVal)

	sinCore(z, xVal)

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

	// Flat +_W guard — same rationale as Sin. The computation path
	// (reducePi2 → cosCore) is identical in structure.
	workPrec := prec + _W

	xVal := newFloat(workPrec).Set(x)
	if xVal.Signbit() {
		xVal.Neg(xVal)
	}

	// reducePi2 handles z == x aliasing; reuse xVal as both input and output.
	quad := reducePi2(xVal, xVal)

	cosCore(z, xVal)

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

	// Flat +_W guard — same rationale as Sin. The computation path
	// (reducePi2 → sincosCore) is identical in structure.
	workPrec := prec + _W

	xVal := newFloat(workPrec).Set(x)
	neg := xVal.Signbit()
	if neg {
		xVal.Neg(xVal)
	}

	// reducePi2 handles z == x aliasing; reuse xVal as both input and output.
	quad := reducePi2(xVal, xVal)

	sincosCore(zs, zc, xVal)

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
