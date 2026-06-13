// SPDX-License-Identifier: MIT

package bigmath

import (
	"math"
	"sync"
)

// Gamma-related constants for argument reduction (Arb convention).
const gammaBeta = 0.2 // β: threshold coefficient for Stirling series convergence

// bernoulliCache caches Bernoulli numbers B₂, B₄, ..., B_{2n} as *Float
// values computed at cachedPrec bits of precision.
type bernoulliCache struct {
	mu         sync.Mutex
	vals       []*Float // vals[k] = B_{2(k+1)} at cachedPrec bits
	cachedPrec uint     // precision at which vals were stored
}

var bernoulli = new(bernoulliCache)

// ensure computes Bernoulli numbers up to B_{2n} at the requested
// precision, stored at a higher guard-bit precision. Cache hit when
// len(vals) >= n AND cachedPrec >= prec.
func (bc *bernoulliCache) ensure(n int, prec uint) {
	bc.mu.Lock()
	defer bc.mu.Unlock()

	if len(bc.vals) >= n && bc.cachedPrec >= prec {
		return
	}

	// Guard bits to absorb rounding in the float recurrence.
	// Per-operation loss: 2 bits/Mul, 1 bit/Add, 2 bits/Quo, plus
	// ⌈log₂(n)⌉ for alternating-sign cancellation.
	guardBits := uint(1.5*float64(n) + 2 + math.Ceil(math.Log2(float64(n))))
	workPrec := prec + guardBits + _W

	bc.vals = computeBernoulliFloats(n, workPrec)
	bc.cachedPrec = prec
}

// floats returns Bernoulli numbers B₂, B₄, ..., B_{2n} as *Float values
// at the given precision.
func (bc *bernoulliCache) floats(n int, prec uint) []*Float {
	bc.ensure(n, prec)
	r := make([]*Float, n)
	for i, v := range bc.vals[:n] {
		if v.Prec() == prec {
			r[i] = v
		} else {
			r[i] = newFloat(prec).Set(v)
		}
	}
	return r
}

// computeBernoulliFloats computes Bernoulli numbers B₂, B₄, ..., B_{2n}
// directly as *Float values using the recurrence:
//
//	B_i = -1/(i+1) · Σ_{k=0}^{i-1} C(i+1, k) · B_k
//
// The binomial coefficient C(i+1, k) is maintained iteratively as a *Float
// within the inner loop.
func computeBernoulliFloats(n int, prec uint) []*Float {
	m := 2 * n
	B := make([]*Float, m+1)

	// B₀ = 1
	B[0] = one

	// B₁ = -½ (exact in binary)
	if m >= 1 {
		B[1] = new(Float).Neg(half)
	}

	// Reusable temporaries for the inner loop
	sum := newFloat(prec)
	term := newFloat(prec)
	s0 := newFloat(prec)
	denom := newFloat(prec)
	binom := newFloat(prec)
	ratio := newFloat(prec)

	for i := 2; i <= m; i++ {
		// For odd i > 1, B_i = 0
		if i%2 == 1 {
			continue
		}
		sum.SetInt64(0)
		binom.SetUint64(1) // C(i+1, 0) = 1

		for k := 0; k < i; {
			// binom = C(i+1, k) — invariant: holds at k=0, maintained by advance below

			if k%2 != 1 || k <= 1 {
				// B_k is non-zero (B₀, B₁, or even k ≥ 2)
				term.Mul(binom, B[k])
				s0.Add(sum, term)
				sum, s0 = s0, sum
			}

			k++
			if k >= i {
				break
			}

			// Advance binom: C(i+1, k) = C(i+1, k-1) · (i+2-k) / k
			ratio.SetInt64(int64(i + 2 - k))
			s0.Mul(binom, ratio) // s0 = binom · (i+2-k)
			ratio.SetInt64(int64(k))
			binom.Quo(s0, ratio) // binom = s0 / k
		}

		// B[i] = -(sum) / (i+1)
		s0.Neg(sum)
		denom.SetInt64(int64(i + 1))
		B[i] = newFloat(prec).Quo(s0, denom)
	}

	// Extract even-index Bernoulli numbers B₂, B₄, ..., B_{2n}
	result := make([]*Float, n)
	for k := range result {
		result[k] = B[2*(k+1)]
	}
	return result
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
	t0 := newFloat(workPrec)
	t1 := newFloat(workPrec)
	betaPrec := gammaBeta * float64(workPrec)
	exp := x.MantExp(nil) // log₂(|x|) — approximate magnitude

	if exp < int(math.Floor(math.Log2(betaPrec))) {
		// |x| < β·prec: need to shift
		xMag := math.Exp2(float64(exp)) // rough |x| in [0.5, 1)·2^exp
		r := max(int(math.Ceil(betaPrec-xMag)), 1)

		// Compute shiftSum = Σ_{j=0}^{r-1} log(x + j)
		// Use a separate temp for the log to avoid mutating the value being accumulated.
		shiftSum = newFloat(workPrec)
		t0.Set(x)
		for range r {
			t1.Log(t0) // logTerm = log(x + j)
			shiftSum.Add(shiftSum, t1)
			t1.Add(t0, one) // t = x + j + 1 for next iteration
			t0, t1 = t1, t0
		}

		// Set x = x + r as the new argument for Stirling
		x = newFloat(workPrec).Add(x, newFloat(workPrec).SetInt64(int64(r)))
	}

	// --- Stirling series evaluation ---
	// log Γ(z) = (z - ½)·log(z) - z + log(2π)/2 + S
	// where S = Σ B_{2k} / (2k·(2k-1)·z^{2k-1})

	// Determine number of terms: approximately prec / (4π·β) ≈ prec / 2.5,
	// but at least 8 terms for low precision.
	maxTerms := max(int(math.Ceil(float64(workPrec)/(4*math.Pi*gammaBeta))), 8)

	// Get Bernoulli numbers B₂ ... B_{2·maxTerms}
	bnums := bernoulli.floats(maxTerms, workPrec)

	// t0 = 1/z²
	t0.Inv(t1.Mul(x, x))
	// t1 tracks z^{-(2k-1)}: starts at z^{-1}, then z^{-3}, z^{-5}, ...
	t1.Inv(x)
	sum := newFloat(workPrec) // series accumulator

	t2 := newFloat(workPrec)
	t3 := newFloat(workPrec)
	t4 := newFloat(workPrec)

	for k := 1; k <= maxTerms; k++ {
		// term = B_{2k} · z^{-(2k-1)} / (2k·(2k-1))
		t2.Mul(bnums[k-1], t1)                // B_{2k} · z^{-(2k-1)}
		t4.SetInt64(int64(2 * k * (2*k - 1))) // (2k)(2k-1)
		t3.Quo(t2, t4)                        // B_{2k} / (2k(2k-1) · z^{2k-1})

		// Stop if term underflows or series starts diverging
		if t3.Sign() == 0 || t3.MantExp(nil) < sum.ULPExponent() {
			break
		}

		// Accumulate: sum += term
		t2.Add(sum, t3)
		sum, t2 = t2, sum

		// Update zInvPow for next iteration: z^{-(2k+1)} = z^{-(2k-1)} / z²
		t2.Mul(t1, t0)
		t2, t1 = t1, t2
	}

	// --- Combine terms ---

	// (z - ½)·log(z)
	t1.Sub(x, half)
	t2.Mul(t1, t0.Log(x)) // (z - ½)·log(z)
	// -z
	t0.Sub(t2, x) // (z - ½)·log(z) - z

	t2.SetMantExp(log2Pi.get(workPrec), -1) // log(2π)/2
	t1.Add(t0, t2)                          // (z-½)·log(z) - z + log(2π)/2

	// Subtract rising factorial if argument reduction was applied
	if shiftSum != nil {
		// + series sum
		t2.Add(t1, sum) // log|Γ(z)| via Stirling
		return z.Sub(t2, shiftSum)
	}

	// + series sum
	return z.Add(t1, sum) // log|Γ(z)| via Stirling
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
	t0 := newFloat(workPrec).Sub(one, x)
	logGamma1mX := newFloat(workPrec).lgammaStirling(t0)

	// 2. Compute sin(πx)
	t2 := newFloat(workPrec).setConst(pi)
	t1 := newFloat(workPrec).Mul(t2, x)
	t0.Sin(t1)

	// Determine sign from sin(πx)
	sign := 1
	if t0.Sign() < 0 {
		sign = -1
	}

	// 3. log|sin(πx)|
	t0.Abs(t0)
	if t0.Sign() == 0 {
		// Pole at integer — should not happen (caught by caller)
		return z.SetInf(false), 1
	}
	t1.Log(t0)

	// 4. log(π)
	t0.Log(t2)

	// 5. log|Γ(x)| = log(π) - log|sin(πx)| - log|Γ(1-x)|
	t2.Sub(t0, t1)
	z.Sub(t2, logGamma1mX)

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
		return z.lgammaReflect(x)
	}

	// x > 0 from here

	// --- x == 1 or x == 2 ---
	if x.Cmp(one) == 0 || x.Cmp(two) == 0 {
		return z.Set(zero), 1
	}

	// --- General case: Stirling series ---
	return z.lgammaStirling(x), 1
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
	logGamma, sign := newFloat(prec + _W).Lgamma(x)

	// Exp handles overflow (returns ±Inf when exponent exceeds MaxExp).
	z.Exp(logGamma)
	if sign < 0 {
		z.Neg(z)
	}

	return z
}
