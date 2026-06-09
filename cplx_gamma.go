// SPDX-License-Identifier: MIT

package bigmath

import "math"

// Lgamma sets z to the principal value of log Γ(x) and returns z.
//
// The branch cut is along the negative real axis. The imaginary part of
// the result is continuous from above across the cut.
//
// Special cases:
//
//	Lgamma(+Inf + i·y)  = +Inf + 0i
//	Lgamma(-Inf + i·y)  = +Inf + 0i
//	Lgamma(x + i·±Inf)  = +Inf + 0i
//	Lgamma(0 + 0i)      = +Inf + 0i
//	Lgamma(-0 + 0i)     = +Inf - iπ
//	Lgamma(n + 0i)      for negative integer n: +Inf + 0i
//	Lgamma(1 + 0i)      = 0 + 0i
//	Lgamma(2 + 0i)      = 0 + 0i
func (z *Complex) Lgamma(x *Complex) *Complex {
	prec := z.setPrec(x)
	workPrec := prec + 2*_W

	// Infinity (any component infinite)
	if x.Real.IsInf() || x.Imag.IsInf() {
		z.Real.SetInf(false) // +Inf
		z.Imag.Set(zero)
		return z
	}

	if x.IsReal() {
		// Short-circuit to Float.Lgamma for purely positive real inputs
		xr := &x.Real
		if !xr.Signbit() {
			// z >= 0: Lgamma is purely real
			lg, _ := newFloat(workPrec).Lgamma(xr)
			z.Real.Set(lg)
			z.Imag.Set(zero)
			return z
		}

		// z < 0
		// Principal branch of log Γ(z) for real z < 0 (non-integer):
		//   log Γ(z) = log|Γ(z)| - i·π   when Γ(z) < 0 (sign = -1)
		//   log Γ(z) = log|Γ(z)|          when Γ(z) > 0 (sign = +1)
		// The -π convention comes from the Euler reflection formula:
		//   log Γ(z) = log(π) - log(sin(πz)) - log Γ(1-z)
		// and the principal Complex.Log: Log(negative + 0i) = log|negative| + iπ,
		// so -Log(sin(πz)) contributes -iπ when sin(πz) < 0.
		_, sign := z.Real.Lgamma(xr)
		if sign < 0 {
			z.Imag.setConst(pi)
			z.Imag.Neg(&z.Imag) // -π
		} else {
			z.Imag.Set(zero)
		}
		return z
	}

	// Reflection for Re(z) < 0.5
	if x.Real.Cmp(half) < 0 {
		return z.lgammaReflect(x, workPrec)
	}

	// General case: Stirling series (with optional rising factorial shift)
	return z.lgammaStirling(x, workPrec)
}

// lgammaReflect computes log Γ(z) for Re(z) < 0.5 via Euler's reflection
// formula: log Γ(z) = log π - log sin(πz) - log Γ(1-z)
//
// The recursion Lgamma(1-z) is safe because Re(1-z) >= 0.5, so no
// reentrancy occurs.
func (z *Complex) lgammaReflect(x *Complex, workPrec uint) *Complex {
	// 1. Compute 1 - z
	w := newComplex(workPrec)
	w.Real.Sub(one, &x.Real)
	w.Imag.Neg(&x.Imag)

	// 2. log Γ(1-z) — recursive call; Re(w) >= 0.5 guarantees termination.
	logGamma1mZ := w.Lgamma(w) // re-use w

	// 3. log sin(πz)
	logSin := newComplex(workPrec)
	logSin.Real.setConst(pi)
	logSin.Mul(logSin, x)
	logSin.Log(logSin.Sin(logSin))

	// 4. result = log π - log sin(πz) - log Γ(1-z)
	//           = log π - (log sin(πz) + log Γ(1-z))
	t := newFloat(workPrec)
	t.Sub(logPi.get(workPrec), &logSin.Real)
	z.Real.Sub(t, &logGamma1mZ.Real)
	// imaginary part of log π is 0, just negate (logSin.Imag+logGamma1mZ.Imag)
	z.Imag.Add(&logSin.Imag, &logGamma1mZ.Imag)
	z.Imag.Neg(&z.Imag)

	return z
}

// lgammaStirling computes log Γ(z) for Re(z) >= 0.5 using the Stirling
// series with Bernoulli numbers. If |z| is small, the argument is shifted
// up via the rising factorial (matching the real lgammaStirling pattern).
//
// The Stirling series:
//
//	log Γ(z) = (z - ½)·log(z) - z + ½·log(2π) + Σ_{k=1}^{N-1} B_{2k} / (2k·(2k-1)·z^{2k-1})
func (z *Complex) lgammaStirling(x *Complex, workPrec uint) *Complex {
	// Argument reduction via rising factorial
	// If |z| < γ·prec, shift z up by r so that |z+r| >= γ·prec.
	tempFloat := x.Abs(newFloat(workPrec))

	betaPrec := gammaBeta * float64(workPrec)

	var shiftSum *Complex // non-nil when argument reduction applied
	t0 := newComplex(workPrec)
	t1 := newComplex(workPrec)

	exp := tempFloat.MantExp(nil) // log₂(|z|) — approximate magnitude
	if exp < int(math.Floor(math.Log2(betaPrec))) {
		xMag := math.Exp2(float64(exp)) // rough |z| in [0.5, 1)·2^exp
		r := max(int(math.Ceil(betaPrec-xMag)), 1)

		shiftSum = newComplex(workPrec)
		t0.Set(x)
		for range r {
			t1.Log(t0) // log(z + j)
			shiftSum.Add(shiftSum, t1)
			tempFloat.Add(&t0.Real, one)              // z + j + 1
			t0.Real, *tempFloat = *tempFloat, t0.Real // swap by value
		}

		// Set x = x + r for Stirling
		t := newComplex(workPrec)
		t.Real.Add(&x.Real, tempFloat.SetInt64(int64(r)))
		t.Imag.Set(&x.Imag)
		x = t
	}

	// Stirling series evaluation
	maxTerms := max(int(math.Ceil(float64(workPrec)/(4*math.Pi*gammaBeta))), 8)
	bnums := bernoulli.floats(maxTerms, workPrec)

	sum := newComplex(workPrec) // series accumulator

	t1.Inv(x)                                  // z^{-(2k-1)}, starts at z^{-1}
	zInvSq := newComplex(workPrec).Mul(t1, t1) // zInvSq = 1/x²
	term := newComplex(workPrec)               // current term
	denom := new(Float)                        // (2k)(2k-1) as float

	for k := 1; k <= maxTerms; k++ {
		// term = B_{2k} · z^{-(2k-1)} / (2k·(2k-1))
		// B_{2k} is real, so multiply components separately.
		term.Real.Mul(bnums[k-1], &t1.Real)
		term.Imag.Mul(bnums[k-1], &t1.Imag)

		denom.SetInt64(int64(2 * k * (2*k - 1)))
		t0.Real.Quo(&term.Real, denom)
		t0.Imag.Quo(&term.Imag, denom)
		term, t0 = t0, term

		// Check convergence on both components.
		if term.Real.Sign() == 0 && term.Imag.Sign() == 0 {
			break
		}
		if term.Real.MantExp(nil) < sum.Real.ULPExponent() &&
			term.Imag.MantExp(nil) < sum.Imag.ULPExponent() {
			break
		}

		// Accumulate: sum += term (pointer swap)
		t0.Add(sum, term)
		sum, t0 = t0, sum

		// Update zPow for next iteration:
		// z^{-(2k+1)} = z^{-(2k-1)} · z^{-2}
		t1.Mul(t1, zInvSq)
	}

	// Combine terms: (z-½)·log(z) - z + ½·log(2π) + Σ
	tempFloat.setConst(log2Pi)
	tempFloat.SetMantExp(tempFloat, -1) // (log(2π)/2, 0)

	// (z - ½)·log(z)
	t0.Real.Sub(&x.Real, half) // z - 0.5
	t0.Imag.Set(&x.Imag)
	t1.Log(x)      // log(z)
	t0.Mul(t0, t1) // (z-0.5)·log(z)

	// - z
	t1.Sub(t0, x) // (z-0.5)·log(z) - z

	// + log(2π)/2
	t0.Real.Add(&t1.Real, tempFloat) // (z-0.5)·log(z) - z + log(2π)/2
	t0.Imag.Set(&t1.Imag)

	// + Σ (series sum)
	t1.Add(t0, sum)

	// Subtract rising factorial if argument reduction was applied
	if shiftSum != nil {
		return z.Sub(t1, shiftSum)
	}
	return z.Set(t1)
}

// Gamma sets z to Γ(x) and returns z.
//
// Special cases:
//
//	Gamma(+Inf + i·y)  = +Inf + 0i
//	Gamma(-Inf + i·y)  = NaN (panics with ErrNaN)
//	Gamma(0 + 0i)      = +Inf + 0i
//	Gamma(-0 + 0i)     = -Inf + 0i
//	Gamma(n + 0i)      for negative integer n: panics with ErrNaN
//	Gamma(1 + 0i)      = 1 + 0i
//	Gamma(2 + 0i)      = 1 + 0i
func (z *Complex) Gamma(x *Complex) *Complex {
	prec := z.setPrec(x)
	workPrec := prec + _W

	// Infinity
	if x.Real.IsInf() || x.Imag.IsInf() {
		if x.Real.IsInf() && x.Real.Signbit() {
			// -Inf in real part → NaN
			panic(ErrNaN("gamma of (-Inf + i·y)"))
		}
		z.Real.SetInf(false) // +Inf
		z.Imag.Set(zero)
		return z
	}

	// Zero
	if x.IsZero() {
		z.Real.SetInf(x.Real.Signbit()) // sign mirrors real component
		z.Imag.Set(zero)
		return z
	}

	if x.IsReal() {
		// Negative real integer (pole)
		if x.Real.Signbit() && x.Real.IsInt() {
			panic(ErrNaN("gamma of negative integer"))
		}

		// Exact values
		// Only exact for purely real 1 or 2; gamma(non-real) ≠ 1 even with Re(z)=1.
		// Float.Gamma short-circuit above handles the real case, so this is for
		// non-real inputs only — but we guard with IsReal to avoid returning (1,0)
		// for e.g. gamma(1+i).
		if x.Real.Cmp(one) == 0 || x.Real.Cmp(two) == 0 {
			z.Real.Set(one)
			z.Imag.Set(zero)
			return z
		}

		// Short-circuit to Float.Gamma for purely real inputs
		// Avoids numerical noise in the imaginary part from the complex Exp path.
		z.Real.Gamma(&x.Real)
		z.Imag.Set(zero)
		return z

	}

	// General case: Gamma(z) = exp(Lgamma(z))
	lgr := newComplex(workPrec).Lgamma(x)
	return z.Exp(lgr)
}
