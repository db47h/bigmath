// SPDX-License-Identifier: MIT

package bigmath

// Cbrt sets z to the principal cube root of x and returns z.
//
// The principal cube root is computed using polar decomposition:
//
//	z = |x|^(1/3) * e^(i * arg(x) / 3)
func (z *Complex) Cbrt(x *Complex) *Complex {
	prec := z.setPrec(x)
	workPrec := prec + _W

	if x.IsZero() {
		return z.Set(&Complex{})
	}

	// Fast path for purely real inputs: avoid polar decomposition.
	if x.IsReal() {
		if x.Real.Sign() >= 0 {
			z.Real.Cbrt(&x.Real)
			z.Imag.Set(zero)
		} else {
			t0 := newFloat(workPrec).Neg(&x.Real)
			t1 := newFloat(workPrec).Cbrt(t0)
			z.Imag.Mul(t1, t0.SetMantExp(sqrt3(workPrec), -1))
			z.Real.SetMantExp(t1.SetPrec(prec), -1)
		}
		return z
	}

	// Handle infinite inputs.
	// For positive-real infinity, the polar form gives θ=0 and
	// ∞ * sin(0) = NaN, so special-case it.
	if x.Real.IsInf() || x.Imag.IsInf() {
		if x.Imag.Sign() == 0 && x.Real.Sign() > 0 {
			z.Real.SetInf(false)
			z.Imag.Set(zero)
			return z
		}
		z.Real.SetInf(false)
		z.Imag.SetInf(false)
		return z
	}

	// Polar decomposition
	r := x.Abs(newFloat(workPrec))
	th := x.Arg(newFloat(workPrec))

	// ρ = cbrt(|x|)
	rho := newFloat(workPrec).Cbrt(r)

	// θ/3
	th3 := newFloat(workPrec).Quo(th, three)

	// z.Real = ρ * cos(θ/3)
	// z.Imag = ρ * sin(θ/3)
	s, c := Sincos(newFloat(workPrec), newFloat(workPrec), th3)
	z.Real.Mul(rho, c)
	z.Imag.Mul(rho, s)

	return z
}
