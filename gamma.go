// SPDX-License-Identifier: MIT

package bigmath

import (
	"math"
	"math/big"
	"sync"
)

// Gamma-related constants for argument reduction (Arb convention).
const gammaBeta = 0.2 // β: threshold coefficient for Stirling series convergence

// bernoulliCache caches Bernoulli numbers B₂, B₄, ..., B_{2n} as exact
// math/big.Rat fractions, converted to *Float on demand.
type bernoulliCache struct {
	mu   sync.Mutex
	vals []*big.Rat // vals[k] = B_{2(k+1)}
}

var bernoulli = new(bernoulliCache)

// ensure computes Bernoulli numbers up to B_{2n} if not already cached.
func (bc *bernoulliCache) ensure(n int) {
	bc.mu.Lock()
	defer bc.mu.Unlock()

	if len(bc.vals) >= n {
		return
	}

	m := 2 * n // highest index needed
	B := make([]*big.Rat, m+1)
	B[0] = new(big.Rat).SetFrac64(1, 1) // B₀ = 1

	// Sequential recurrence: B_i = -1/(i+1) · Σ_{k=0}^{i-1} C(i+1, k) · B_k
	binom := new(big.Int)
	term := new(big.Rat)
	for i := 1; i <= m; i++ {
		// For odd i > 1, B_i = 0
		if i > 1 && i%2 == 1 {
			B[i] = new(big.Rat)
			continue
		}

		sum := new(big.Rat)
		for k := 0; k < i; k++ {
			binom.Binomial(int64(i+1), int64(k))
			term.SetFrac(binom, big.NewInt(1))
			term.Mul(term, B[k])
			sum.Add(sum, term)
		}

		// B_i = -sum / (i+1)
		B[i] = new(big.Rat).Neg(sum)
		B[i].Mul(B[i], new(big.Rat).SetFrac64(1, int64(i+1)))
	}

	// Extract even-index Bernoulli numbers B₂, B₄, ..., B_{2n}
	// vals[0] = B₂, vals[1] = B₄, ...
	bc.vals = make([]*big.Rat, n)
	for i := range bc.vals {
		bc.vals[i] = B[2*(i+1)]
	}
}

// floats returns Bernoulli numbers B₂, B₄, ..., B_{2n} as *Float values
// at the given precision.
func (bc *bernoulliCache) floats(n int, prec uint) []*Float {
	bc.ensure(n)
	r := make([]*Float, n)
	for i, v := range bc.vals[:n] {
		r[i] = newFloat(prec).SetRat(v)
	}
	return r
}

// lgammaStirling computes log|Γ(x)| for x ≥ ½ using the Stirling series
// with Bernoulli numbers. x must be positive. The result is written to z.
//
// The Stirling series:
//
//	log Γ(z) = (z - ½)·log(z) - z + log(2π)/2 + Σ_{k=1}^{N-1} B_{2k} / (2k·(2k-1)·z^{2k-1})
//
// Argument reduction via rising factorial shifts small x up to the
// convergent range (|x| ≥ β·prec bits).
func (z *Float) lgammaStirling(x *Float) *Float {
	prec := z.Prec()
	workPrec := prec + _W

	// --- Argument reduction ---
	// If |x| < β·prec, shift x up by r via rising factorial:
	//   Γ(x) = Γ(x+r) / (∏_{j=0}^{r-1} (x+j))
	//   log|Γ(x)| = log|Γ(x+r)| - Σ_{j=0}^{r-1} log(x+j)

	var shiftSum *Float // non-nil when argument reduction is applied
	betaPrec := gammaBeta * float64(prec)
	exp := x.MantExp(nil) // log₂(|x|) — approximate magnitude

	if exp < int(math.Floor(math.Log2(betaPrec))) {
		// |x| < β·prec: need to shift
		xMag := math.Exp2(float64(exp)) // rough |x| in [0.5, 1)·2^exp
		r := int(math.Ceil(betaPrec - xMag))
		if r < 1 {
			r = 1
		}

		// Compute shiftSum = Σ_{j=0}^{r-1} log(x + j)
		// Use a separate temp for the log to avoid mutating the value being accumulated.
		shiftSum = newFloat(workPrec)
		logTerm := newFloat(workPrec)
		t := newFloat(workPrec).Set(x)
		for j := 0; j < r; j++ {
			logTerm.Log(t) // logTerm = log(x + j)
			shiftSum.Add(shiftSum, logTerm)
			t.Add(t, one) // t = x + j + 1 for next iteration
		}

		// Set x = x + r as the new argument for Stirling
		x = newFloat(workPrec).Set(x)
		x.Add(x, newFloat(workPrec).SetInt64(int64(r)))
	}

	// --- Stirling series evaluation ---
	// log Γ(z) = (z - ½)·log(z) - z + log(2π)/2 + S
	// where S = Σ B_{2k} / (2k·(2k-1)·z^{2k-1})

	// Determine number of terms: approximately prec / (4π·β) ≈ prec / 2.5,
	// but at least 8 terms for low precision.
	maxTerms := int(math.Ceil(float64(prec) / (4 * math.Pi * gammaBeta)))
	if maxTerms < 8 {
		maxTerms = 8
	}

	// Get Bernoulli numbers B₂ ... B_{2·maxTerms}
	bnums := bernoulli.floats(maxTerms, workPrec)

	// Pre-compute z² and 1/z
	zSquared := newFloat(workPrec).Mul(x, x) // z²
	zInv := newFloat(workPrec).Quo(one, x)   // 1/z

	// zInvPow tracks z^{-(2k-1)}: starts at z^{-1}, then z^{-3}, z^{-5}, ...
	zInvPow := newFloat(workPrec).Set(zInv)

	sum := newFloat(workPrec).SetPrec(workPrec) // series accumulator
	t0 := newFloat(workPrec)
	t1 := newFloat(workPrec)
	t2 := newFloat(workPrec)

	for k := 1; k <= maxTerms; k++ {
		// term = B_{2k} · z^{-(2k-1)} / (2k·(2k-1))
		t0.Mul(bnums[k-1], zInvPow)           // B_{2k} · z^{-(2k-1)}
		t2.SetInt64(int64(2 * k * (2*k - 1))) // (2k)(2k-1)
		t1.Quo(t0, t2)                        // B_{2k} / (2k(2k-1) · z^{2k-1})

		// Stop if term underflows or series starts diverging
		if t1.Sign() == 0 || t1.MantExp(nil) < sum.ULPExponent() {
			break
		}

		// Accumulate: sum += term
		t0.Add(sum, t1)
		sum, t0 = t0, sum // pointer swap

		// Update zInvPow for next iteration: z^{-(2k+1)} = z^{-(2k-1)} / z²
		t0.Quo(zInvPow, zSquared)
		zInvPow, t0 = t0, zInvPow // pointer swap
	}

	// --- Combine terms ---
	// log(z) term
	logZ := newFloat(workPrec).Log(x) // log(z)

	// (z - ½)·log(z)
	t0.Sub(x, half)
	t0.Mul(t0, logZ) // (z - ½)·log(z)

	// -z
	t0.Sub(t0, x) // (z - ½)·log(z) - z

	// + log(2π)/2 = (ln(2) + ln(π))/2
	ln2Val := newFloat(workPrec).Set(ln2(workPrec))
	piCopy := newFloat(workPrec).Set(pi(workPrec))
	logPi := newFloat(workPrec).Log(piCopy)
	t1.Add(ln2Val, logPi) // ln(2π)
	t1.Mul(t1, half)      // ln(2π)/2
	t0.Add(t0, t1)        // (z-½)·log(z) - z + log(2π)/2

	// + series sum
	t0.Add(t0, sum) // log|Γ(z)| via Stirling

	// Subtract rising factorial if argument reduction was applied
	if shiftSum != nil {
		t0.Sub(t0, shiftSum)
	}

	return z.SetPrec(prec).Set(t0)
}

// lgammaReflect computes log|Γ(x)| and the sign for x < 0 (non-integer)
// via the reflection formula:
//
//	log|Γ(x)| = log(π) - log|sin(πx)| - log|Γ(1-x)|
//	sign = sign(sin(πx))
func (z *Float) lgammaReflect(x *Float) (*Float, int) {
	prec := z.Prec()
	workPrec := prec + 2*_W // extra guard bits for sin(πx) near integers

	// 1. Compute 1 - x
	oneMinusX := newFloat(workPrec).Sub(one, x)
	logGamma1mX := newFloat(workPrec).SetPrec(workPrec).lgammaStirling(oneMinusX)

	// 2. Compute sin(πx)
	piCopy := newFloat(workPrec).Set(pi(workPrec))
	piX := newFloat(workPrec).Mul(piCopy, x)
	sinVal := newFloat(workPrec + _W).Sin(piX) // Sin uses workPrec+_W internally

	// Determine sign from sin(πx)
	sign := 1
	if sinVal.Sign() < 0 {
		sign = -1
	}

	// 3. log|sin(πx)|
	absSin := newFloat(workPrec).Abs(sinVal)
	if absSin.Sign() == 0 {
		// Pole at integer — should not happen (caught by caller)
		return z.SetInf(false), 1
	}
	logSin := newFloat(workPrec).Log(absSin)

	// 4. log(π)
	piCopy2 := newFloat(workPrec).Set(pi(workPrec))
	logPiVal := newFloat(workPrec).Log(piCopy2)

	// 5. log|Γ(x)| = log(π) - log|sin(πx)| - log|Γ(1-x)|
	z.SetPrec(prec)
	z.Sub(logPiVal, logSin)
	z.Sub(z, logGamma1mX)

	return z, sign
}

// Lgamma sets z to the natural logarithm of |Γ(x)| and returns (z, sign).
// sign is +1 when Γ(x) > 0 and -1 when Γ(x) < 0.
// At poles (zero and negative integers), Lgamma returns (+Inf, +1).
//
// The precision of z determines the working precision; if z's precision is 0,
// it is set from x's precision.
//
// Special cases:
//
//	Lgamma(+Inf) = (+Inf, +1)
//	Lgamma(-Inf) = (+Inf, -1)
//	Lgamma(+0)   = (+Inf, +1)
//	Lgamma(-0)   = (+Inf, -1)
//	Lgamma(1)    = (0, +1)
//	Lgamma(2)    = (0, +1)
//	Lgamma(n)    for negative integer n: (+Inf, +1)
func (z *Float) Lgamma(x *Float) (*Float, int) {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	// --- Infinity ---
	if x.IsInf() {
		if x.Signbit() {
			return z.SetInf(false), -1 // -Inf → (+Inf, -1), matching gmpy2 convention
		}
		return z.SetInf(false), 1 // +Inf → (+Inf, +1)
	}

	s := x.Sign()

	// --- Zero ---
	if s == 0 {
		if x.Signbit() {
			return z.SetInf(false), -1 // -0 → (+Inf, -1)
		}
		return z.SetInf(false), 1 // +0 → (+Inf, +1)
	}

	// --- Negative values ---
	if x.Signbit() {
		// x < 0
		if x.IsInt() {
			return z.SetInf(false), 1 // pole at negative integer
		}
		// Non-integer negative: use reflection
		workPrec := prec + 2*_W
		return newFloat(workPrec).SetPrec(workPrec).lgammaReflect(x)
	}

	// x > 0 from here

	// --- x == 1 or x == 2 ---
	if x.Cmp(one) == 0 || x.Cmp(two) == 0 {
		return z.Set(zero), 1
	}

	// --- General case: Stirling series ---
	workPrec := prec + _W
	result := newFloat(workPrec).SetPrec(workPrec).lgammaStirling(x)
	return z.Set(result), 1
}

// Gamma sets z to Γ(x) and returns z.
//
// The precision of z determines the working precision; if z's precision is 0,
// it is set from x's precision.
//
// Special cases:
//
//	Gamma(+Inf) = +Inf
//	Gamma(-Inf) = NaN (panics with ErrNaN)
//	Gamma(+0)   = +Inf
//	Gamma(-0)   = -Inf
//	Gamma(n)    for negative integer n: panics with ErrNaN
//	Gamma(1)    = 1
//	Gamma(2)    = 1
func (z *Float) Gamma(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	// --- Infinity ---
	if x.IsInf() {
		if x.Signbit() {
			panic(ErrNaN("gamma of -Inf"))
		}
		return z.SetInf(false)
	}

	s := x.Sign()

	// --- Zero ---
	if s == 0 {
		if x.Signbit() {
			return z.SetInf(true) // -0 → -Inf
		}
		return z.SetInf(false) // +0 → +Inf
	}

	// --- Negative integer (pole) ---
	if s < 0 && x.IsInt() {
		panic(ErrNaN("gamma of negative integer"))
	}

	// --- Exact values ---
	if x.Cmp(one) == 0 || x.Cmp(two) == 0 {
		return z.Set(one)
	}

	// --- General case: Gamma(x) = exp(Lgamma(x)), then apply sign ---
	workPrec := prec + _W
	logGamma, sign := newFloat(workPrec).SetPrec(workPrec).Lgamma(x)

	// Exp handles overflow (returns ±Inf when exponent exceeds MaxExp).
	result := newFloat(workPrec).Exp(logGamma)
	if sign < 0 {
		result.Neg(result)
	}

	return z.Set(result)
}
