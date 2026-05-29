// SPDX-License-Identifier: MIT

package bigmath

import "math/big"

// Pow sets z to x^y and returns z.
//
//	Pow(0, ±0) returns 1+0i
//	Pow(0, c) for real(c)<0 returns Inf+0i if imag(c) is zero, otherwise Inf+Inf i.
func (z *Complex) Pow(x, y *Complex) *Complex {
	workPrec := z.setPrec2(x, y) + _W

	t0 := x.Abs(newFloat(workPrec))
	if t0.Sign() == 0 {
		switch sgn := y.Real.Sign(); {
		case sgn == 0:
			z.Real.Set(one)
			z.Imag.Set(zero)
		case sgn < 0:
			z.Real.SetInf(y.IsReal())
			z.Imag.Set(zero)
		default:
			z.Set(&Complex{})
		}
		return z
	}

	// Special case: y real
	if y.IsReal() {
		if x.IsReal() {
			// If x.Real: delegate to real-domain Pow.
			z.Real.Pow(&x.Real, &y.Real)
			z.Imag.Set(zero)
			return z
		} else if n, acc := y.Real.Int64(); acc == big.Exact {
			// x Complex: use Complex multiplication.
			return z.powInt(x, n)
		}
	}

	t1 := x.Arg(newFloat(workPrec))
	t2 := newFloat(workPrec)
	r := newFloat(workPrec).Pow(t0, &y.Real)
	theta := newFloat(workPrec).Mul(&y.Real, t1)
	if y.Imag.Sign() != 0 {
		t2.Add(theta, t0.Mul(&y.Imag, t2.Log(t0)))
		theta, t2 = t2, theta
		t1.Exp(t0.Neg(t0.Mul(&y.Imag, t1))) // t1 = e^(-d⋅arg(x))
		t0.Mul(r, t1)
		t0, r = r, t0
	}
	Sincos(t0, t1, theta)
	z.Real.Mul(r, t1)
	z.Imag.Mul(r, t0)
	return z
}

// powInt computes z^n using exact complex multiplications
func (z *Complex) powInt(x *Complex, n int64) *Complex {
	workPrec := z.setPrec(x) + _W

	// Handle negative exponents: x^(-n) = (1/x)^n
	isNegative := n < 0
	base := newComplex(workPrec).Set(x)
	t0 := newComplex(workPrec)
	if isNegative {
		n = -n
		t0.Inv(base) // Assuming you have a reciprocal/inverse function: 1 / base
		t0, base = base, t0
	}

	// Initialize result to 1 + 0i
	res := newComplex(workPrec)
	res.Real.Set(one)
	res.Imag.Set(zero)

	// Power-by-squaring loop
	for n > 0 {
		if n&1 == 1 {
			t0.Mul(res, base)
			t0, res = res, t0
		}
		t0.Mul(base, base)
		t0, base = base, t0
		n >>= 1
	}

	z.Set(res)
	return z
}

// Sqrt sets z to the square root of x and returns z.
func (z *Complex) Sqrt(x *Complex) *Complex {
	// Algebraic square root: sqrt(a+bi) = sqrt((r+a)/2) + i·sgn(b)*sqrt((r-a)/2)
	workPrec := z.setPrec(x) + _W

	if x.IsReal() {
		switch x.Real.Sign() {
		case -1:
			z.Imag.Sqrt(new(Float).Neg(&x.Real))
			if x.Imag.Signbit() {
				z.Imag.Neg(&z.Imag)
			}
			z.Real.Set(zero)
		case 0:
			z.Real.Set(zero)
			z.Imag.Set(&x.Imag)
		case 1:
			z.Real.Sqrt(&x.Real)
			z.Imag.Set(&x.Imag)
		}
		return z
	}
	if x.Imag.IsInf() {
		z.Real.SetInf(false)
		z.Imag.Set(&x.Imag)
		return z
	}
	if x.Real.Sign() == 0 {
		if x.Imag.Sign() < 0 {
			r := new(Float).SetMantExp(&x.Imag, -1)
			z.Real.Sqrt(r.Neg(r))
			z.Imag.Neg(&z.Real)
			return z
		}
		z.Real.Sqrt(new(Float).SetMantExp(&x.Imag, -1))
		z.Imag.Set(&z.Real)
		return z
	}

	// r = |x|
	r := x.Abs(newFloat(workPrec))
	si := x.Imag.Sign()
	a := &x.Real
	if x == z {
		a = new(Float).Copy(a)
	}

	// real part = sqrt((r + a) / 2)
	t := newFloat(workPrec).Add(r, a)
	t.SetMantExp(t, -1)
	z.Real.Sqrt(t)

	// imag part = sqrt((r - a) / 2)
	t.Sub(r, a)
	t.SetMantExp(t, -1)
	z.Imag.Sqrt(t)
	if si < 0 {
		z.Imag.Neg(&z.Imag)
	}

	return z
}
