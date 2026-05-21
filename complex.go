// SPDX-License-Identifier: MIT

package bigmath

import (
	"math"
	"math/big"
)

// Complex represents a complex number with real and imaginary parts as *big.Float.
type Complex struct {
	Real big.Float
	Imag big.Float
}

func (z *Complex) Copy(x *Complex) *Complex {
	if z == x {
		return z
	}
	z.Real.Copy(&x.Real)
	z.Imag.Copy(&x.Imag)
	return z
}

// Add sets z to x+y and returns z.
func (z *Complex) Add(x, y *Complex) *Complex {
	// check for aliasing
	if z == x {
		x = new(Complex).Copy(x)
	}
	if z == y {
		y = new(Complex).Copy(y)
	}
	z.Real.Add(&x.Real, &y.Real)
	z.Imag.Add(&x.Imag, &y.Imag)
	return z
}

// Sub sets z to x-y and returns z.
func (z *Complex) Sub(x, y *Complex) *Complex {
	// check for aliasing
	if z == x {
		x = new(Complex).Copy(x)
	}
	if z == y {
		y = new(Complex).Copy(y)
	}
	z.Real.Sub(&x.Real, &y.Real)
	z.Imag.Sub(&x.Imag, &y.Imag)
	return z
}

func resultPrec(x, y *Complex) uint {
	return max(x.Real.Prec(), x.Imag.Prec(), y.Real.Prec(), y.Imag.Prec())
}

// Mul sets z to x*y and returns z.
func (z *Complex) Mul(x, y *Complex) *Complex {
	// check for aliasing
	// TODO: handle the case where z = x = y
	if z == x {
		x = new(Complex).Copy(x)
	}
	if z == y {
		y = new(Complex).Copy(y)
	}

	if z.Real.Prec() == 0 {
		z.Real.SetPrec(resultPrec(x, y))
	}
	if z.Imag.Prec() == 0 {
		z.Imag.SetPrec(resultPrec(x, y))
	}

	t := new(big.Float)

	// (a+bi)(c+di) = (ac-bd) + (ad+bc)i
	// Real: ac - bd
	bd := newFloat(x.Imag.Prec()+y.Imag.Prec()).Mul(&x.Imag, &y.Imag)
	fma(&z.Real, &x.Real, &y.Real, bd.Neg(bd), t)

	// Imag: ad + bc
	bc := newFloat(x.Imag.Prec()+y.Real.Prec()).Mul(&x.Imag, &y.Real)
	fma(&z.Imag, &x.Real, &y.Imag, bc, t)

	return z
}

// Quo sets z to x/y and returns z.
func (z *Complex) Quo(x, y *Complex) *Complex {
	// check for aliasing
	// TODO: handle the case where z = x = y
	if z == x {
		x = new(Complex).Copy(x)
	}
	if z == y {
		y = new(Complex).Copy(y)
	}

	t := new(big.Float)
	prec := max(resultPrec(x, y))
	if z.Real.Prec() == 0 {
		z.Real.SetPrec(prec)
	}
	if z.Imag.Prec() == 0 {
		z.Imag.SetPrec(prec)
	}
	workPrec := prec + 2

	// (a+bi)/(c+di) = ((ac+bd) + (bc-ad)i) / (c^2+d^2)
	// denom = c^2 + d^2.
	c2 := newFloat(y.Real.Prec()*2).Mul(&y.Real, &y.Real)
	denom := fma(newFloat(workPrec), &y.Imag, &y.Imag, c2, t)

	// ac + bd
	ac := newFloat(x.Real.Prec()+y.Real.Prec()).Mul(&x.Real, &y.Real)
	re := fma(newFloat(workPrec), &x.Imag, &y.Imag, ac, t)
	z.Real.Quo(re, denom)

	// bc - ad
	ad := newFloat(x.Real.Prec()+y.Imag.Prec()).Mul(&x.Real, &y.Imag)
	im := fma(newFloat(workPrec), &x.Imag, &y.Real, ad.Neg(ad), t)
	z.Imag.Quo(im, denom)

	return z
}

// Neg sets z to -x and returns z.
func (z *Complex) Neg(x *Complex) *Complex {
	z.Real.Neg(&x.Real)
	z.Imag.Neg(&x.Imag)
	return z
}

// Conj sets z to the complex conjugate of x and returns z.
func (z *Complex) Conj(x *Complex) *Complex {
	z.Real.Copy(&x.Real)
	z.Imag.Neg(&x.Imag)
	return z
}

// Abs sets res to the rounded value of |x| and returns res.
func (z *Complex) Abs(res *big.Float, x *Complex) *big.Float {
	return Hypot(res, &x.Real, &x.Imag)
}

// Arg sets res to the rounded value of arg(x) and returns res.
func (z *Complex) Arg(res *big.Float, x *Complex) *big.Float {
	return Atan2(res, &x.Imag, &x.Real)
}

// Log sets z to the principal logarithm of x and returns z.
//
// The branch cut is along the negative real axis. The imaginary part of
// the result lies in the interval [-π, π].
//
// Special cases:
//
//	Log(0) = -Inf + i·0
//	Log(+Inf + i·y) = +Inf + i·0
//	Log(-Inf + i·y) = +Inf + i·π for finite y
//	Log(x + i·±Inf) = +Inf + i·π/2 for finite x
func (z *Complex) Log(x *Complex) *Complex {
	// check for aliasing
	if z == x {
		x = new(Complex).Copy(x)
	}

	prec := resultPrec(x, x)
	if z.Real.Prec() == 0 {
		z.Real.SetPrec(prec)
	}
	if z.Imag.Prec() == 0 {
		z.Imag.SetPrec(prec)
	}

	// ln(x+iy) = ln|x+iy| + i*arg(x+iy)
	// ln|x+iy| = 0.5 * ln(x^2 + y^2)
	// We use Abs and then real Log to avoid precision loss.
	t := newFloat(prec + _W)
	x.Abs(t, x)
	Log(&z.Real, t)

	x.Arg(&z.Imag, x)

	return z
}

// Atan sets z to the rounded value of arctan(x) and returns z.
func (z *Complex) Atan(x *Complex) *Complex {
	// atan(z) = (i/2) * ln((1-iz)/(1+iz))
	// let w = (1-iz)/(1+iz)
	// atan(z) = (i/2) * Log(w)

	// check for aliasing
	if z == x {
		x = new(Complex).Copy(x)
	}

	prec := resultPrec(x, x)
	workPrec := prec + _W
	oneC := &Complex{Real: *one, Imag: *zero}
	iz := &Complex{Real: *newFloat(workPrec).Neg(&x.Imag), Imag: *newFloat(workPrec).Copy(&x.Real)}

	// num = 1 - iz
	num := new(Complex).Sub(oneC, iz)
	// den = 1 + iz
	den := new(Complex).Add(oneC, iz)

	w := new(Complex).Quo(num, den)
	lw := new(Complex).Log(w)

	// z = (i/2) * lw = (-lw.Imag/2) + i*(lw.Real/2)
	z.Real.SetMantExp(&lw.Imag, -1).Neg(&z.Real)
	z.Imag.SetMantExp(&lw.Real, -1)

	return z
}

// Equals checks if x and y are equal.
func (x *Complex) Equals(y *Complex) bool {
	return x.Real.Cmp(&y.Real) == 0 && x.Imag.Cmp(&y.Imag) == 0
}

// IsReal checks if the complex number is a real number (imaginary part is zero).
func (x *Complex) IsReal() bool {
	return x.Imag.Sign() == 0
}

// IsZero checks if the complex number is zero.
func (x *Complex) IsZero() bool {
	return x.Real.Sign() == 0 && x.Imag.Sign() == 0
}

// Exp sets z to e^x and returns z.
// Uses float64-precision sin/cos for the imaginary part.
// TODO: use big.Float trig for full precision.
func (z *Complex) Exp(x *Complex) *Complex {
	if z == x {
		x = new(Complex).Copy(x)
	}

	prec := resultPrec(x, x)
	if z.Real.Prec() == 0 {
		z.Real.SetPrec(prec)
	}
	if z.Imag.Prec() == 0 {
		z.Imag.SetPrec(prec)
	}

	var expA big.Float
	Exp(&expA, &x.Real)

	b64, _ := x.Imag.Float64()
	s, c := math.Sincos(b64)

	z.Real.Mul(&expA, new(big.Float).SetFloat64(c))
	z.Imag.Mul(&expA, new(big.Float).SetFloat64(s))

	return z
}
