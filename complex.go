// SPDX-License-Identifier: MIT

package bigmath

import (
	"fmt"
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
	prec := resultPrec(x, y)
	if z.Real.Prec() == 0 {
		z.Real.SetPrec(prec)
	}
	if z.Imag.Prec() == 0 {
		z.Imag.SetPrec(prec)
	}

	workPrec := prec + _W
	t := newFloat(workPrec)
	re := newFloat(workPrec)
	im := newFloat(workPrec)

	// (a+bi)(c+di) = (ac-bd) + (ad+bc)i
	// Real: ac - bd
	bd := newFloat(x.Imag.Prec()+y.Imag.Prec()).Mul(&x.Imag, &y.Imag)
	fma(re, &x.Real, &y.Real, bd.Neg(bd), t)

	// Imag: ad + bc
	bc := newFloat(x.Imag.Prec()+y.Real.Prec()).Mul(&x.Imag, &y.Real)
	fma(im, &x.Real, &y.Imag, bc, t)

	z.Real.Set(re)
	z.Imag.Set(im)
	return z
}

// Quo sets z to x/y and returns z.
func (z *Complex) Quo(x, y *Complex) *Complex {
	prec := resultPrec(x, y)
	if z.Real.Prec() == 0 {
		z.Real.SetPrec(prec)
	}
	if z.Imag.Prec() == 0 {
		z.Imag.SetPrec(prec)
	}

	workPrec := prec + _W
	t := newFloat(workPrec)
	re := newFloat(workPrec)
	im := newFloat(workPrec)

	// (a+bi)/(c+di) = ((ac+bd) + (bc-ad)i) / (c^2+d^2)
	// denom = c^2 + d^2.
	c2 := newFloat(y.Real.Prec()*2).Mul(&y.Real, &y.Real)
	denom := fma(newFloat(workPrec), &y.Imag, &y.Imag, c2, t)

	// ac + bd
	ac := newFloat(x.Real.Prec()+y.Real.Prec()).Mul(&x.Real, &y.Real)
	fma(re, &x.Imag, &y.Imag, ac, t)

	// bc - ad
	ad := newFloat(x.Real.Prec()+y.Imag.Prec()).Mul(&x.Real, &y.Imag)
	fma(im, &x.Imag, &y.Real, ad.Neg(ad), t)

	z.Real.Quo(re, denom)
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

// Abs sets z to the rounded value of |x| and returns z.
func (x *Complex) Abs(z *big.Float) *big.Float {
	return Hypot(z, &x.Real, &x.Imag)
}

// Arg sets z to the rounded value of arg(x) and returns z.
func (x *Complex) Arg(z *big.Float) *big.Float {
	return Atan2(z, &x.Imag, &x.Real)
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
	x.Abs(t)
	Log(&z.Real, t)

	x.Arg(&z.Imag)

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

	prec := max(x.Real.Prec(), x.Imag.Prec())
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

	var s, c big.Float
	Sincos(&s, &c, &x.Imag)

	z.Real.Mul(&expA, &c)
	z.Imag.Mul(&expA, &s)

	return z
}

// Sin sets z to the sine of x and returns z.
//
// Formula: sin(a+bi) = sin(a)cosh(b) + i·cos(a)sinh(b)
func (z *Complex) Sin(x *Complex) *Complex {
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

	workPrec := prec + _W
	var s, c, sh, ch big.Float
	s.SetPrec(workPrec)
	c.SetPrec(workPrec)
	sh.SetPrec(workPrec)
	ch.SetPrec(workPrec)

	Sincos(&s, &c, &x.Real)
	SinhCosh(&sh, &ch, &x.Imag)

	z.Real.Mul(&s, &ch)
	z.Imag.Mul(&c, &sh)

	return z
}

// Cos sets z to the cosine of x and returns z.
//
// Formula: cos(a+bi) = cos(a)cosh(b) − i·sin(a)sinh(b)
func (z *Complex) Cos(x *Complex) *Complex {
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

	workPrec := prec + _W
	var s, c, sh, ch big.Float
	s.SetPrec(workPrec)
	c.SetPrec(workPrec)
	sh.SetPrec(workPrec)
	ch.SetPrec(workPrec)

	Sincos(&s, &c, &x.Real)
	SinhCosh(&sh, &ch, &x.Imag)

	z.Real.Mul(&c, &ch)
	z.Imag.Mul(&s, &sh).Neg(&z.Imag)

	return z
}

// Tan sets z to the tangent of x and returns z.
func (z *Complex) Tan(x *Complex) *Complex {
	// tan(z) = sin(z) / cos(z)
	s := new(Complex).Sin(x)
	c := new(Complex).Cos(x)
	return z.Quo(s, c)
}

// Sinh sets z to the hyperbolic sine of x and returns z.
//
// Formula: sinh(a+bi) = sinh(a)cos(b) + i·cosh(a)sin(b)
func (z *Complex) Sinh(x *Complex) *Complex {
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

	workPrec := prec + _W
	var s, c, sh, ch big.Float
	s.SetPrec(workPrec)
	c.SetPrec(workPrec)
	sh.SetPrec(workPrec)
	ch.SetPrec(workPrec)

	Sincos(&s, &c, &x.Imag)
	SinhCosh(&sh, &ch, &x.Real)

	z.Real.Mul(&sh, &c)
	z.Imag.Mul(&ch, &s)

	return z
}

// Cosh sets z to the hyperbolic cosine of x and returns z.
//
// Formula: cosh(a+bi) = cosh(a)cos(b) + i·sinh(a)sin(b)
func (z *Complex) Cosh(x *Complex) *Complex {
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

	workPrec := prec + _W
	var s, c, sh, ch big.Float
	s.SetPrec(workPrec)
	c.SetPrec(workPrec)
	sh.SetPrec(workPrec)
	ch.SetPrec(workPrec)

	Sincos(&s, &c, &x.Imag)
	SinhCosh(&sh, &ch, &x.Real)

	z.Real.Mul(&ch, &c)
	z.Imag.Mul(&sh, &s)

	return z
}

// Tanh sets z to the hyperbolic tangent of x and returns z.
func (z *Complex) Tanh(x *Complex) *Complex {
	// tanh(z) = sinh(z) / cosh(z)
	sh := new(Complex).Sinh(x)
	ch := new(Complex).Cosh(x)
	return z.Quo(sh, ch)
}

// Asin sets z to the inverse sine of x and returns z.
//
// Formula: asin(z) = −i · ln(i·z + √(1−z²))
func (z *Complex) Asin(x *Complex) *Complex {
	if z == x {
		x = new(Complex).Copy(x)
	}

	prec := resultPrec(x, x)
	workPrec := prec + _W
	oneC := &Complex{Real: *one, Imag: *zero}

	// z2 = z^2
	z2 := new(Complex).Mul(x, x)
	// t0 = 1 - z^2
	t0 := new(Complex).Sub(oneC, z2)
	// t1 = sqrt(1 - z^2)
	t1 := new(Complex).Sqrt(t0)

	// iz = i * z
	iz := &Complex{Real: *newFloat(workPrec).Neg(&x.Imag), Imag: *newFloat(workPrec).Copy(&x.Real)}

	// t2 = iz + t1
	t2 := new(Complex).Add(iz, t1)

	// res = -i * ln(t2)
	ln := new(Complex).Log(t2)
	z.Real.Set(&ln.Imag)
	z.Imag.Set(&ln.Real).Neg(&z.Imag)

	return z
}

// Acos sets z to the inverse cosine of x and returns z.
//
// Formula: acos(z) = −i · ln(z + i·√(1−z²))
func (z *Complex) Acos(x *Complex) *Complex {
	if z == x {
		x = new(Complex).Copy(x)
	}

	prec := resultPrec(x, x)
	workPrec := prec + _W
	oneC := &Complex{Real: *one, Imag: *zero}

	// z2 = z^2
	z2 := new(Complex).Mul(x, x)
	// t0 = 1 - z^2
	t0 := new(Complex).Sub(oneC, z2)
	// t1 = sqrt(1 - z^2)
	t1 := new(Complex).Sqrt(t0)

	// it1 = i * t1
	it1 := &Complex{Real: *newFloat(workPrec).Neg(&t1.Imag), Imag: *newFloat(workPrec).Copy(&t1.Real)}

	// t2 = z + it1
	t2 := new(Complex).Add(x, it1)

	// res = -i * ln(t2)
	ln := new(Complex).Log(t2)
	z.Real.Set(&ln.Imag)
	z.Imag.Set(&ln.Real).Neg(&z.Imag)

	return z
}

// Asinh sets z to the inverse hyperbolic sine of x and returns z.
//
// Formula: asinh(z) = ln(z + √(z²+1))
func (z *Complex) Asinh(x *Complex) *Complex {
	if z == x {
		x = new(Complex).Copy(x)
	}

	oneC := &Complex{Real: *one, Imag: *zero}

	// z2 = z^2
	z2 := new(Complex).Mul(x, x)
	// t0 = z^2 + 1
	t0 := new(Complex).Add(z2, oneC)
	// t1 = sqrt(z^2 + 1)
	t1 := new(Complex).Sqrt(t0)
	// t2 = z + t1
	t2 := new(Complex).Add(x, t1)

	return z.Log(t2)
}

// Acosh sets z to the inverse hyperbolic cosine of x and returns z.
//
// Formula: acosh(z) = ln(z + √(z−1)·√(z+1))
func (z *Complex) Acosh(x *Complex) *Complex {
	if z == x {
		x = new(Complex).Copy(x)
	}

	oneC := &Complex{Real: *one, Imag: *zero}

	// t0 = z - 1
	t0 := new(Complex).Sub(x, oneC)
	// t1 = z + 1
	t1 := new(Complex).Add(x, oneC)
	// t2 = sqrt(z-1) * sqrt(z+1)
	t2 := new(Complex).Mul(new(Complex).Sqrt(t0), new(Complex).Sqrt(t1))
	// t3 = z + t2
	t3 := new(Complex).Add(x, t2)

	return z.Log(t3)
}

// Atanh sets z to the inverse hyperbolic tangent of x and returns z.
//
// Formula: atanh(z) = ½ · ln((1+z)/(1−z))
func (z *Complex) Atanh(x *Complex) *Complex {
	if z == x {
		x = new(Complex).Copy(x)
	}

	oneC := &Complex{Real: *one, Imag: *zero}

	// num = 1 + z
	num := new(Complex).Add(oneC, x)
	// den = 1 - z
	den := new(Complex).Sub(oneC, x)

	w := new(Complex).Quo(num, den)
	lw := new(Complex).Log(w)

	// z = lw / 2
	z.Real.SetMantExp(&lw.Real, -1)
	z.Imag.SetMantExp(&lw.Imag, -1)

	return z
}

// Sqrt sets z to the square root of x and returns z.
func (z *Complex) Sqrt(x *Complex) *Complex {
	// sqrt(z) = exp(0.5 * log(z))
	l := new(Complex).Log(x)
	l.Real.SetMantExp(&l.Real, -1)
	l.Imag.SetMantExp(&l.Imag, -1)
	return z.Exp(l)
}

// Pow sets z to x^y and returns z.
func (z *Complex) Pow(x, y *Complex) *Complex {
	// x^y = exp(y * log(x))
	l := new(Complex).Log(x)
	return z.Exp(z.Mul(y, l))
}

func (x *Complex) String() string {
	return fmt.Sprint(x)
}

func (x *Complex) Format(s fmt.State, verb rune) {
	x.Real.Format(s, verb)
	if x.Imag.Sign() >= 0 && !x.Imag.IsInf() {
		fmt.Fprintf(s, "+")
	}
	x.Imag.Format(s, verb)
	fmt.Fprint(s, "i")
}

func (x *Complex) Prec() uint {
	return max(x.Real.Prec(), x.Imag.Prec())
}
