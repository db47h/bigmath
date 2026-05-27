// SPDX-License-Identifier: MIT

// Guard-bit strategy: core Taylor-series loops (sinCore, cosCore, sincosCore)
// use a flat +2*_W guard. Tan shares the same cores and strategy. This is adequate because the accumulated rounding
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
// x must be non-negative; panics otherwise.
func (z *Float) reducePi2(x *Float) int {
	if x.Sign() == 0 {
		return 0
	}

	prec := z.Prec()

	xAbs := newFloat(prec).Set(x)
	if xAbs.Signbit() {
		panic("reducePi2: x must be >= 0")
	}

	xExp := xAbs.MantExp(nil)
	workPrec := prec + _W
	if xExp > 0 {
		workPrec += uint(xExp)
	}

	// t0 = xAbs * (2/π)
	t0 := newFloat(workPrec).Mul(xAbs, twoOverPi(workPrec))

	// Get the number of bits required to represent the integer part
	E := t0.MantExp(nil)
	t1 := newFloat(workPrec)
	var quadrant int

	if E <= 0 {
		// x * (2/π) < 1, so the integer part n = 0
		z.Set(xAbs)
		return 0
	}

	// 1. Truncate t0 to an exact integer (Floor)
	// By setting the precision exactly to the exponent with ToZero rounding,
	// we keep all integer bits and drop all fractional bits in place.
	t0.SetMode(big.ToZero)
	t0.SetPrec(uint(E))

	// 2. Fast Modulo 4
	if E <= 2 {
		// E is 1 or 2, so the integer is <= 3. Modulo 4 is just the number itself.
		q64, _ := t0.Int64()
		quadrant = int(q64)
	} else {
		// Create a copy rounded to E-2 bits. This keeps bits down to the 2^2 (4) place,
		// effectively dropping the lowest 2 bits (equivalent to calculating t0 - (t0 % 4)).
		t1.SetPrec(uint(E - 2)).SetMode(big.ToZero).Set(t0)

		// The difference is exactly the modulo 4 remainder (0, 1, 2, or 3).
		rem := newFloat(workPrec).Sub(t0, t1)
		q64, _ := rem.Int64()
		quadrant = int(q64)

		// reset t1
		t1.SetMode(0).SetPrec(0).SetPrec(workPrec)
	}

	// -x + n * π/2
	xAbs.Neg(xAbs)
	rTmp := newFloat(workPrec).fma(t0, halfPi(workPrec), xAbs, t1)

	// x - n * π/2
	rTmp.Neg(rTmp)

	z.Set(rTmp)
	return quadrant
}

// sinCore computes sin(x) for x in [0, π/2) using the Taylor series.
// sin(x) = Σ (−1)ⁿ · x^(2n+1) / (2n+1)!
//
// workPrec uses a flat +2*_W guard. The alternating series introduces
// subtractive cancellation as terms approach the sum, so +2*_W provides
// a comfortable margin. At target precision P, O(P) terms are needed;
// the accumulated error (N × 2^-(P+2*_W)) is negligible for any N reachable
// in practice (see file header for the full rationale).
func (z *Float) sinCore(x *Float) *Float {
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

		if term.Sign() == 0 || term.MantExp(nil) < sum.ULPExponent() {
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

// cosCore computes cos(x) for x in [0, π/2) using the Taylor series.
// cos(x) = Σ (−1)ⁿ · x^(2n) / (2n)!
//
// Same flat +2*_W guard as sinCore — same alternating-series cancellation
// characteristics. See sinCore doc for the rationale.
func (z *Float) cosCore(x *Float) *Float {
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

		if term.Sign() == 0 || term.MantExp(nil) < sum.ULPExponent() {
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

// sincosCore computes both sin(x) and cos(x) for x in [0, π/2)
// using a single Taylor series loop that shares the computation of x²
// and the factorial denominator between both series.
//
// Same flat +2*_W guard as sinCore/cosCore. The shared loop handles both
// alternating series; the guard covers the combined rounding error.
func sincosCore(zs, zc, x *Float) (*Float, *Float) {
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

		sinDone := sinTerm.Sign() == 0 || sinTerm.MantExp(nil) < sinSum.ULPExponent()
		cosDone := cosTerm.Sign() == 0 || cosTerm.MantExp(nil) < cosSum.ULPExponent()
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

	return zs.Set(sinSum), zc.Set(cosSum)
}

// Sin sets z to the sine of x and returns z.
//
// Special cases:
//
//	Sin(±0) = ±0
//	Sin(±Inf) = panic(ErrNaN)
func (z *Float) Sin(x *Float) *Float {
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
	quad := xVal.reducePi2(xVal)

	if quad == 0 || quad == 2 {
		z.sinCore(xVal)
	} else {
		z.cosCore(xVal)
	}

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
func (z *Float) Cos(x *Float) *Float {
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
	quad := xVal.reducePi2(xVal)

	if quad == 0 || quad == 2 {
		z.cosCore(xVal)
	} else {
		z.sinCore(xVal)
	}

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
func Sincos(zs, zc, x *Float) (*Float, *Float) {
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
	quad := xVal.reducePi2(xVal)

	sincosCore(zs, zc, xVal)
	if quad == 1 || quad == 3 {
		*zs, *zc = *zc, *zs
	}
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

// Tan sets z to the tangent of x and returns z.
//
// Special cases:
//
//	Tan(±0) = ±0
//	Tan(±Inf) = panic(ErrNaN)
//	Tan(π/2 + nπ) = ±Inf
func (z *Float) Tan(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.IsInf() {
		panic(ErrNaN("tan of infinity"))
	}
	if x.Sign() == 0 {
		return z.Set(x)
	}

	// Flat +_W guard — same rationale as Sin/Cos. The computation path
	// (reducePi2 → sinCore/cosCore → Quo) is identical in structure.
	workPrec := prec + _W

	xVal := newFloat(workPrec).Set(x)
	neg := xVal.Signbit()
	if neg {
		xVal.Neg(xVal)
	}

	// reducePi2 handles z == x aliasing; reuse xVal as both input and output.
	quad := xVal.reducePi2(xVal)

	s := newFloat(workPrec).sinCore(xVal)
	c := newFloat(workPrec).cosCore(xVal)

	// tan(x) after reduction to [0, π/2):
	//   Q0: tan = sinR / cosR   → s / c
	//   Q1: tan = cosR / -sinR  → -c / s
	//   Q2: tan = -sinR / -cosR → s / c
	//   Q3: tan = -cosR / sinR  → -c / s
	// The formula choice determines the sign; no extra quadrant sign flip needed.
	if quad&1 == 0 {
		// Q0, Q2: tan = s / c
		if c.Sign() == 0 {
			z.SetInf(false)
			if neg {
				z.Neg(z)
			}
			return z
		}
		z.Quo(s, c)
	} else {
		// Q1, Q3: tan = -c / s
		if s.Sign() == 0 {
			z.SetInf(true)
			if neg {
				z.Neg(z)
			}
			return z
		}
		c.Neg(c)
		z.Quo(c, s)
	}

	if neg {
		z.Neg(z)
	}
	return z
}
