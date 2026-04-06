// SPDX-License-Identifier: MIT

package bigmath

import (
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
