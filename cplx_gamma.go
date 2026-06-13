// SPDX-License-Identifier: MIT

package bigmath

import "math"

// Lgamma sets z to the principal value of log Γ(x) and returns z.
//
// The principal branch of the log-gamma function has a single branch cut
// along the negative real axis. It is continuous from above across the
// cut: for negative real x (non-integer),
//
//	log Γ(x) = log|Γ(x)| - (⌊|x|⌋ + 1)·π·i
//
// For x just below the cut (Im(x) < 0, Re(x) < 0) the imaginary part
// jumps to +(⌊|x|⌋ + 1)·π·i.
//
// The real part matches ln|Γ(z)| for all z where defined; the imaginary
// part differs from the actual ln(Γ(z)) by an integer multiple of 2π
// (the branch choice).
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
	workPrec := prec + _W

	// Infinity (any component infinite)
	if x.Real.IsInf() || x.Imag.IsInf() {
		z.Real.SetInf(false) // +Inf
		z.Imag.Set(zero)
		return z
	}

	if x.IsReal() {
		xr := &x.Real
		if !xr.Signbit() {
			// Positive real — use Float.Lgamma directly
			lg, _ := newFloat(workPrec).Lgamma(xr)
			z.Real.Set(lg)
			z.Imag.Set(zero)
			return z
		}

		// Negative real — check for integer pole
		if xr.IsInt() {
			// Non-zero negative integer: pole
			if xr.Sign() != 0 {
				z.Imag.Set(zero)
			} else {
				// -0+0i — branch cut: limit from above gives +Inf - iπ
				z.Imag.Neg(pi.get(workPrec)) // -π
			}
			z.Real.SetInf(false) // +Inf
			return z
		}

		// Principal branch for non-integer negative real z.
		// Compute imag part FIRST to avoid aliasing: z.Real.Lgamma(xr)
		// would overwrite xr when z == x.
		// For negative non-integer x: |⌊x⌋| = ⌊|x|⌋ + 1.
		floorX := newFloat(workPrec).Floor(xr) // ⌊x⌋ (negative)
		floorX.Neg(floorX)                     // |⌊x⌋|
		z.Imag.Mul(floorX, pi.get(workPrec))   // |⌊x⌋|·π
		z.Imag.Neg(&z.Imag)                    // -|⌊x⌋|·π = -(⌊|x|⌋+1)·π

		z.Real.Lgamma(xr) // log|Γ(x)|
		return z
	}

	// For non-real z:
	if x.Real.Signbit() {
		return z.lgammaNegateReflect(x, workPrec)
	}

	// General case: Stirling series (with optional rising factorial shift)
	return z.lgammaStirling(x, workPrec)
}

// lgammaNegateReflect computes log Γ(z) for Re(z) < 0 via
// negation to the right half-plane + Wolfram reflection formula
func (z *Complex) lgammaNegateReflect(x *Complex, workPrec uint) *Complex {
	t0 := newComplex(workPrec).Neg(x)
	t1 := newComplex(workPrec).Log(t0) // Log(-z)
	t2 := newComplex(workPrec)
	pi := pi.get(workPrec) // cache pi@workPrec

	// 1. Compute LogGamma[-z] via Lgamma (Re(-z) ≥ 0, goes to lgammaStirling)
	t0.Lgamma(t0)

	// 2. result = -LogGamma[-z] - Log[-z]
	t2.Neg(t2.Add(t0, t1))

	// 3. Branch correction: Sign(Im(z)) · ⌊Re(z)⌋ · π · i
	f0 := newFloat(workPrec).Floor(&x.Real) // ⌊Re(z)⌋
	f1 := newFloat(workPrec).Mul(f0, pi)

	if x.Imag.Signbit() {
		f1.Neg(f1) // flip sign for negative imag
	}
	t0.Imag.Add(&t2.Imag, f1)

	// 4. + Log(π)
	t0.Real.Add(&t2.Real, logPi.get(workPrec))

	// 5. − Log(Sin(π·(z − ⌊Re(z)⌋)))
	t1.Real.Sub(&x.Real, f0)
	t1.Imag.Set(&x.Imag)

	// π·frac
	t2.Real.Mul(&t1.Real, pi)
	t2.Imag.Mul(&t1.Imag, pi)

	// Log(Sin(π·frac))
	t2.Log(t2.Sin(t2))
	return z.Sub(t0, t2)
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
	bernoulli.ensure(maxTerms, workPrec)

	sum := newComplex(workPrec) // series accumulator

	t1.Inv(x)                                  // z^{-(2k-1)}, starts at z^{-1}
	zInvSq := newComplex(workPrec).Mul(t1, t1) // zInvSq = 1/x²
	term := newComplex(workPrec)               // current term
	denom := new(Float)                        // (2k)(2k-1) as float

	for k := 1; k <= maxTerms; k++ {
		// term = B_{2k} · z^{-(2k-1)} / (2k·(2k-1))
		// B_{2k} is real, so multiply components separately.
		bk := bernoulli.get(k-1, workPrec)
		term.Real.Mul(bk, &t1.Real)
		term.Imag.Mul(bk, &t1.Imag)

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
// Γ(z) is meromorphic with simple poles at z = 0, -1, -2, ... and is
// otherwise analytic everywhere. It is computed via Γ(z) = exp(Lgamma(z)),
// which eliminates the Lgamma branch cut — the function has no branch cuts.
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
