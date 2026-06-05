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
// modPi2, which scales with input magnitude for large arguments.
// See docs/trig-hyperbolic-precision-review.md for the full analysis.

package bigmath

import (
	"math/big"
)

// modPi2 reduces x modulo π/2 and maps the result strictly to [-π/4, π/4).
// Returns the reduced value in z and the quadrant mapping:
//
//	0: [-π/4, π/4)
//	1: [π/4, 3π/4)
//	2: [3π/4, 5π/4)
//	3: [5π/4, 7π/4)
func (z *Float) modPi2(x *Float) (*Float, int) {
	// x == 0 || |x| < 2^-1
	if x.Sign() == 0 || x.MantExp(nil) < 0 {
		return x, 0
	}
	prec := z.Prec()

	// 1. Compute q = x * (2/π) + 0.5
	// Scale precision by the magnitude of x so the integer part of x * (2/π)
	// is exactly representable. Without this, |x| > 2^prec makes q + 0.5 a no-op
	// for large x.
	xExp := x.MantExp(nil)
	multPrec, _ := addPrec(prec, uint(xExp))
	q := newFloat(multPrec).setConst(twoOverPi)
	t := newFloat(multPrec).Mul(x, q)
	q.Add(t, half)

	// 2. q = Floor(q)
	q.SetMode(ToNegativeInf)
	E := q.MantExp(nil)
	if E > 0 {
		q.SetPrec(uint(E))
	} else if q.Sign() < 0 {
		q.Copy(minusOne)
	} else {
		return z.Set(x), 0
	}

	// 4. Extract Octant Mapping (0, 1, 2, 3)
	qInt, _ := q.Int(new(big.Int)) // Int() strictly truncates toward zero
	var quad int
	if qInt.Sign() < 0 {
		qInt.Neg(qInt)
		quad = (4 - int(qInt.Bit(1)<<1|qInt.Bit(0))) % 4
	} else {
		quad = int(qInt.Bit(1)<<1 | qInt.Bit(0))
	}

	// Dynamic precision loop (Ziv's strategy)
	var r Float
	workPrec, _ := addPrec(prec, 2*_W)
	for {
		// 3. Compute r = x - q * (π/2)
		// Force a copy: pi() can return a value with a much higher precision
		// and Mul uses the full precision of its arguments.
		pi := newFloat(workPrec).Pi()
		qPi2 := newFloat(workPrec).Mul(q, pi.SetMantExp(pi, -1))
		r.SetPrec(workPrec).Sub(x, qPi2)

		// 4. Measure bit loss from cancellation
		if r.Sign() == 0 {
			break
		}
		loss := max(x.MantExp(nil)-r.MantExp(nil), 0)

		// 5. Minimum precision required to guarantee targetPrec bits in the result
		requiredPrec, overflow := addPrec(prec, uint(loss))
		if overflow || workPrec >= requiredPrec {
			break
		}

		// 6. Threshold breached: escalate precision.
		workPrec, overflow = addPrec(requiredPrec, _W)
		if overflow {
			break
		}
	}

	return z.Set(&r), quad
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
	if x.MantExp(nil) < -int(prec/2) {
		return z.Set(x)
	}
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
	if x.MantExp(nil) < -int(prec/2) {
		return z.Set(one)
	}
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
	// x at this precision feeds into modPi2 and then into sinCore (whose temps
	// are at prec+2*_W). The guard covers sign handling, reduction output rounding,
	// and the Taylor series accumulation.
	workPrec := prec + _W

	// modPi2 handles z == x aliasing; reuse xVal as both input and output.
	xr, quad := newFloat(workPrec).modPi2(x)

	if quad == 0 || quad == 2 {
		z.sinCore(xr)
	} else {
		z.cosCore(xr)
	}
	if quad >= 2 {
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
	// (modPi2 → cosCore) is identical in structure.
	workPrec := prec + _W

	// modPi2 handles z == x aliasing; reuse xVal as both input and output.
	xr, quad := newFloat(workPrec).modPi2(x)

	if quad == 0 || quad == 2 {
		z.cosCore(xr)
	} else {
		z.sinCore(xr)
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
	// (modPi2 → sincosCore) is identical in structure.
	workPrec := prec + _W

	// modPi2 handles z == x aliasing; reuse xr as both input and output.
	xr, quad := newFloat(workPrec).modPi2(x)

	sincosCore(zs, zc, xr)
	if quad == 1 || quad == 3 {
		*zs, *zc = *zc, *zs
	}
	if quad >= 2 {
		zs.Neg(zs)
	}
	if quad == 1 || quad == 2 {
		zc.Neg(zc)
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

	// Flat +_W guard, reduced from +2*_W after empirical testing up to
	// 2048 bits (see TestTanULP). The +_W guard in sincosCore is sufficient,
	// and the Quo of two prec+_W operands into prec produces a correctly
	// rounded result for all tested inputs.
	workPrec := prec + _W

	// modPi2 handles z == x aliasing; reuse xVal as both input and output.
	xr, quad := newFloat(workPrec).modPi2(x)
	s, c := sincosCore(newFloat(workPrec), newFloat(workPrec), xr)

	// tan(x) after reduction to [0, π/2):
	//   Q0: tan = sinR / cosR   → s / c
	//   Q1: tan = cosR / -sinR  → -c / s
	//   Q2: tan = -sinR / -cosR → s / c
	//   Q3: tan = -cosR / sinR  → -c / s
	// The formula choice determines the sign; no extra quadrant sign flip needed.
	if quad&1 == 0 {
		// Q0, Q2: tan = s / c
		z.Quo(s, c)
	} else {
		// Q1, Q3: tan = -c / s
		c.Neg(c)
		z.Quo(c, s)
	}
	return z
}

// Cot sets z to the cotangent of x and returns z.
//
// Special cases:
//
//	Cot(±0) = ±Inf
//	Cot(±Inf) = panic(ErrNaN)
func (z *Float) Cot(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.IsInf() {
		panic(ErrNaN("cot of infinity"))
	}

	workPrec := prec + _W

	xr, quad := newFloat(workPrec).modPi2(x)

	s, c := sincosCore(newFloat(workPrec), newFloat(workPrec), xr)

	// cot(x) after reduction to [0, π/2):
	//   Q0: cot = cosR / sinR   → c / s
	//   Q1: cot = -sinR / cosR  → -s / c
	//   Q2: cot = -cosR / -sinR → c / s
	//   Q3: cot = sinR / -cosR  → -s / c
	if quad&1 == 0 {
		z.Quo(c, s)
	} else {
		s.Neg(s)
		z.Quo(s, c)
	}
	return z
}

// Sec sets z to the secant of x, sec(x) = 1/cos(x), and returns z.
//
// Special cases:
//
//	Sec(±0) = 1
//	Sec(±Inf) = panic(ErrNaN)
func (z *Float) Sec(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}
	if x.IsInf() {
		panic(ErrNaN("sec of infinity"))
	}
	if x.Sign() == 0 {
		return z.Set(one)
	}
	t := newFloat(prec + _W).Cos(x)
	return z.Inv(t)
}

// Csc sets z to the cosecant of x, csc(x) = 1/sin(x), and returns z.
//
// Special cases:
//
//	Csc(±0) = ±Inf
//	Csc(±Inf) = panic(ErrNaN)
func (z *Float) Csc(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}
	if x.IsInf() {
		panic(ErrNaN("csc of infinity"))
	}
	t := newFloat(prec + _W).Sin(x)
	return z.Inv(t)
}

// Acot sets z to the inverse cotangent of x, acot(x) = atan2(1, x), and returns z.
//
// Special cases:
//
//	Acot(±0) = ±π/2
//	Acot(±Inf) = ±0
func (z *Float) Acot(x *Float) *Float {
	return z.Atan2(one, x)
}

// Asec sets z to the inverse secant of x, asec(x) = acos(1/x), and returns z.
//
// Special cases:
//
//	Asec(|x| < 1) = panic(ErrNaN)
//	Asec(1) = 0
//	Asec(-1) = π
//	Asec(±Inf) = π/2
func (z *Float) Asec(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}
	if x.IsInf() {
		z.Pi()
		z.SetMantExp(z, -1)
		return z
	}
	switch x.absCmpOne() {
	case -1:
		panic(ErrNaN("asec of |x| < 1"))
	case 0:
		if x.Signbit() {
			return z.Pi()
		}
		return z.Set(zero)
	}
	return z.Acos(newFloat(prec + _W).Inv(x))
}

// Acsc sets z to the inverse cosecant of x, acsc(x) = asin(1/x), and returns z.
//
// Special cases:
//
//	Acsc(|x| < 1) = panic(ErrNaN)
//	Acsc(±1) = ±π/2
//	Acsc(±Inf) = ±0
func (z *Float) Acsc(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}
	if x.IsInf() {
		// return Asin(1/±Inf) = ±0
		sgn := x.Signbit()
		z.Set(zero)
		if sgn {
			z.Neg(z)
		}
		return z
	}
	if x.absCmpOne() < 0 {
		panic(ErrNaN("acsc of |x| < 1"))
	}
	return z.Asin(newFloat(prec + _W).Inv(x))
}

// addPrec returns prec + extra saturated to big.MaxPrec, and a boolean
// indicating whether the addition exceeded MaxPrec (or wrapped on 32-bit).
func addPrec(prec, extra uint) (uint, bool) {
	sum := prec + extra
	if sum < prec || sum > MaxPrec {
		return MaxPrec, true
	}
	return sum, false
}
