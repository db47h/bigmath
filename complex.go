// SPDX-License-Identifier: MIT

package bigmath

import (
	"fmt"
	"math/big"
)

var oneC = &Complex{Real: *one}

// Complex represents a complex number with real and imaginary parts as big.Float.
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
	z.setPrec2(x, y)
	z.Real.Add(&x.Real, &y.Real)
	z.Imag.Add(&x.Imag, &y.Imag)
	return z
}

// Sub sets z to x-y and returns z.
func (z *Complex) Sub(x, y *Complex) *Complex {
	z.setPrec2(x, y)
	z.Real.Sub(&x.Real, &y.Real)
	z.Imag.Sub(&x.Imag, &y.Imag)
	return z
}

// Mul sets z to the rounded product x*y and returns z.
// If z's precision is 0, it is changed to the larger of x's or y's precision
// before the operation.
func (z *Complex) Mul(x, y *Complex) *Complex {
	if (x.Real.Sign() == 0 && y.Real.IsInf()) || (x.Imag.Sign() == 0 && y.Imag.IsInf()) {
		// Define result as ±Inf according to sign rules, avoid panic.
	}
	workPrec := z.setPrec2(x, y) + _W

	temp := new(big.Float) // temp for fma. prec will be handled by fma()
	re := newFloat(workPrec)
	im := newFloat(workPrec)

	// (a+bi)(c+di) = (ac-bd) + (ad+bc)i
	// Real: ac - bd
	bd := newFloat(x.Imag.Prec()+y.Imag.Prec()).Mul(&x.Imag, &y.Imag)
	fma(re, &x.Real, &y.Real, bd.Neg(bd), temp)

	// Imag: ad + bc
	bc := newFloat(x.Imag.Prec()+y.Real.Prec()).Mul(&x.Imag, &y.Real)
	fma(im, &x.Real, &y.Imag, bc, temp)

	z.Real.Set(re)
	z.Imag.Set(im)
	return z
}

// Quo sets z to x/y and returns z.
func (z *Complex) Quo(x, y *Complex) *Complex {
	workPrec := z.setPrec2(x, y) + _W

	temp := new(big.Float) // temp for fma. prec will be handled by fma()
	re := newFloat(workPrec)
	im := newFloat(workPrec)

	// (a+bi)/(c+di) = ((ac+bd) + (bc-ad)i) / (c^2+d^2)
	// denom = c^2 + d^2.
	c2 := newFloat(y.Real.Prec()*2).Mul(&y.Real, &y.Real)
	denom := fma(newFloat(workPrec), &y.Imag, &y.Imag, c2, temp)

	// ac + bd
	ac := newFloat(x.Real.Prec()+y.Real.Prec()).Mul(&x.Real, &y.Real)
	fma(re, &x.Imag, &y.Imag, ac, temp)

	// bc - ad
	ad := newFloat(x.Real.Prec()+y.Imag.Prec()).Mul(&x.Real, &y.Imag)
	fma(im, &x.Imag, &y.Real, ad.Neg(ad), temp)

	z.Real.Quo(re, denom)
	z.Imag.Quo(im, denom)

	return z
}

// Neg sets z to -x and returns z.
func (z *Complex) Neg(x *Complex) *Complex {
	z.setPrec(x)
	z.Real.Neg(&x.Real)
	z.Imag.Neg(&x.Imag)
	return z
}

// Conj sets z to the complex conjugate of x and returns z.
func (z *Complex) Conj(x *Complex) *Complex {
	z.setPrec(x)
	z.Real.Set(&x.Real)
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
	prec := z.setPrec(x)
	workPrec := prec + _W

	// ln(x+iy) = ln|x+iy| + i*arg(x+iy)
	// ln|x+iy| = 0.5 * ln(x^2 + y^2)
	// We use Abs and then real Log to avoid precision loss.
	abs := x.Abs(newFloat(workPrec))
	x.Arg(&z.Imag)
	Log(&z.Real, abs)

	return z
}

// Atan sets z to the inverse tangent of x and returns z.
//
// Formula: atan(x) = (i/2) · ln((1−ix)/(1+ix))
//
// The branch cut is along the imaginary axis, outside the interval [-i, +i].
// The real part of the result lies in the interval [-π/2, π/2].
//
// Special cases:
//
//	Atan(0 + i·0) = 0 + i·0
func (z *Complex) Atan(x *Complex) *Complex {
	// atan(z) = (i/2) * ln((1-iz)/(1+iz))
	// let w = (1-iz)/(1+iz)
	// atan(z) = (i/2) * Log(w)
	workPrec := z.setPrec(x) + _W

	iz := &Complex{Real: *new(big.Float).Neg(&x.Imag), Imag: x.Real}

	// num = 1 - iz
	num := newComplex(workPrec).Sub(oneC, iz)
	// den = 1 + iz
	den := newComplex(workPrec).Add(oneC, iz)

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

// Exp sets z to e^x, the base-e exponential of x, and returns z.
//
// Special cases are:
//
//	Exp(x + i·±Inf) = panic                     (infinite imaginary part)
//	Exp(x + i·0)    = Exp(x)                    (for any x)
//	Exp(-Inf + i·y) = 0                         (exact underflow, for finite y)
//	Exp(+Inf + i·y) = +Inf·(cos(y) + i·sin(y))  (for any finite y)
func (z *Complex) Exp(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	if x.Imag.IsInf() {
		panic(ErrNaN("complex exponential of infinite imaginary part"))
	}

	expA := Exp(newFloat(workPrec), &x.Real)
	if x.Imag.Sign() == 0 {
		z.Real.Set(expA)
		z.Imag.Set(&x.Imag)
		return z
	}
	s, c := Sincos(newFloat(workPrec), newFloat(workPrec), &x.Imag)
	if expA.IsInf() {
		switch {
		case expA.Signbit():
			z.Real.Set(zero)
			z.Imag.Set(zero)
			return z
		case c.Sign() == 0:
			// in this unlikely case, choose geometric continuity over panicking
			z.Real.SetInf(false)
			z.Imag.SetInf(s.Signbit())
			return z
		}
	}

	z.Real.Mul(expA, c)
	z.Imag.Mul(expA, s)

	return z
}

// Sin sets z to the sine of x and returns z.
//
// Formula: sin(a+bi) = sin(a)cosh(b) + i·cos(a)sinh(b)
func (z *Complex) Sin(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	s, c := Sincos(newFloat(workPrec), newFloat(workPrec), &x.Real)
	sh, ch := SinhCosh(newFloat(workPrec), newFloat(workPrec), &x.Imag)

	z.Real.Mul(s, ch)
	z.Imag.Mul(c, sh)

	return z
}

// Cos sets z to the cosine of x and returns z.
//
// Formula: cos(a+bi) = cos(a)cosh(b) − i·sin(a)sinh(b)
func (z *Complex) Cos(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	s, c := Sincos(newFloat(workPrec), newFloat(workPrec), &x.Real)
	sh, ch := SinhCosh(newFloat(workPrec), newFloat(workPrec), &x.Imag)

	z.Real.Mul(c, ch)
	z.Imag.Mul(s, sh).Neg(&z.Imag)

	return z
}

// Tan sets z to the tangent of x and returns z.
func (z *Complex) Tan(x *Complex) *Complex {
	// TODO: take this as an example for future review of guard bits strategy.
	// for example, Sin/Cos already add guard bits, and here we need at least 2 more for Quo.
	workPrec := z.setPrec(x) + _W

	// tan(z) = sin(z) / cos(z)
	s := newComplex(workPrec).Sin(x)
	c := newComplex(workPrec).Cos(x)
	return z.Quo(s, c)
}

// Sinh sets z to the hyperbolic sine of x and returns z.
//
// Formula: sinh(a+bi) = sinh(a)cos(b) + i·cosh(a)sin(b)
func (z *Complex) Sinh(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	s, c := Sincos(newFloat(workPrec), newFloat(workPrec), &x.Imag)
	sh, ch := SinhCosh(newFloat(workPrec), newFloat(workPrec), &x.Real)

	z.Real.Mul(sh, c)
	z.Imag.Mul(ch, s)

	return z
}

// Cosh sets z to the hyperbolic cosine of x and returns z.
//
// Formula: cosh(a+bi) = cosh(a)cos(b) + i·sinh(a)sin(b)
func (z *Complex) Cosh(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	s, c := Sincos(newFloat(workPrec), newFloat(workPrec), &x.Imag)
	sh, ch := SinhCosh(newFloat(workPrec), newFloat(workPrec), &x.Real)

	z.Real.Mul(ch, c)
	z.Imag.Mul(sh, s)

	return z
}

// Tanh sets z to the hyperbolic tangent of x and returns z.
func (z *Complex) Tanh(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W
	// tanh(z) = sinh(z) / cosh(z)
	sh := newComplex(workPrec).Sinh(x)
	ch := newComplex(workPrec).Cosh(x)
	return z.Quo(sh, ch)
}

// Asin sets z to the inverse sine of x and returns z.
//
// Formula: asin(x) = −i · ln(i·x + √(1−x²))
//
// The branch cut is along the real axis, outside the interval [-1, +1].
// The real part of the result lies in the interval [-π/2, π/2].
//
// Special cases:
//
//	Asin(0 + i·0) = 0 + i·0
func (z *Complex) Asin(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	// t0 = x^2
	t0 := newComplex(workPrec).Mul(x, x)
	// t1 = 1 - x^2
	t1 := newComplex(workPrec).Sub(oneC, t0)
	// t0 = sqrt(1 - x^2)
	t0.Sqrt(t1)

	// t1 = i * x
	ix := &Complex{Real: *new(big.Float).Neg(&x.Imag), Imag: x.Real}

	// t1 = i·x + √(1−x²)
	t1.Add(ix, t0)

	// res = -i * ln(t2)
	t0.Log(t1)
	z.Real.Set(&t0.Imag)
	z.Imag.Set(&t0.Real).Neg(&z.Imag)

	return z
}

// Acos sets z to the inverse cosine of x and returns z.
//
// Formula: acos(x) = −i · ln(x + i·√(1−x²))
//
// The branch cut is along the real axis, outside the interval [-1, +1].
// The real part of the result lies in the interval [0, π].
//
// Special cases:
//
//	Acos(0 + i·0) = π/2 + i·0
func (z *Complex) Acos(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	// t0 = x^2
	t0 := newComplex(workPrec).Mul(x, x)
	// t1 = 1 - z^2
	t1 := newComplex(workPrec).Sub(oneC, t0)
	// t0 = sqrt(1 - z^2)
	t0.Sqrt(t1)

	// it = i * t0
	t0.Imag.Neg(&t0.Imag)
	it := &Complex{t0.Imag, t0.Real}

	// t1 = x + it
	t1.Add(x, it)

	// res = -i * ln(t0)
	t0.Log(t1)
	z.Real.Set(&t0.Imag)
	z.Imag.Set(&t0.Real).Neg(&z.Imag)

	return z
}

// Asinh sets z to the inverse hyperbolic sine of x and returns z.
//
// Formula: asinh(x) = ln(x + √(x²+1))
//
// The branch cut is along the imaginary axis, outside the interval [-i, +i].
// The imaginary part of the result lies in the interval [-π/2, π/2].
//
// Special cases:
//
//	Asinh(0 + i·0) = 0 + i·0
func (z *Complex) Asinh(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	// t0 = x^2
	t0 := newComplex(workPrec).Mul(x, x)
	// t1 = x^2 + 1
	t1 := newComplex(workPrec).Add(t0, oneC)
	// t0 = sqrt(x^2 + 1)
	t0.Sqrt(t1)
	// t1 = x + t0
	t1.Add(x, t0)

	return z.Log(t1)
}

// Acosh sets z to the inverse hyperbolic cosine of x and returns z.
//
// Formula: acosh(x) = ln(x + √(x−1)·√(x+1))
//
// The branch cut is along the real axis, for x < 1.
// The imaginary part of the result lies in the interval [0, π].
//
// Special cases:
//
//	Acosh(0 + i·0) = 0 + i·π/2
func (z *Complex) Acosh(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	// t0 = sqrt(x - 1)
	t1 := newComplex(workPrec).Sub(x, oneC)
	t0 := newComplex(workPrec).Sqrt(t1)
	// t1 = x + 1
	t2 := newComplex(workPrec).Add(x, oneC)
	t1.Sqrt(t2)
	// t2 = sqrt(x-1) * sqrt(x+1)
	t2.Mul(t0, t1)
	// log(x + t2)
	return z.Log(t0.Add(x, t2))
}

// Atanh sets z to the inverse hyperbolic tangent of x and returns z.
//
// Formula: atanh(x) = ½ · ln((1+x)/(1−x))
//
// The branch cut is along the real axis, outside the interval [-1, +1].
// The imaginary part of the result lies in the interval [-π/2, π/2].
//
// Special cases:
//
//	Atanh(0 + i·0) = 0 + i·0
func (z *Complex) Atanh(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W

	// num = 1 + z
	num := newComplex(workPrec).Add(oneC, x)
	// den = 1 - z
	den := newComplex(workPrec).Sub(oneC, x)

	w := newComplex(workPrec).Quo(num, den)
	// directly use z as target. The z /= 2 is lossless
	z.Log(w)

	// z /= 2
	z.Real.SetMantExp(&z.Real, -1)
	z.Imag.SetMantExp(&z.Imag, -1)

	return z
}

// Sqrt sets z to the square root of x and returns z.
func (z *Complex) Sqrt(x *Complex) *Complex {
	// Algebraic square root: sqrt(a+bi) = sqrt((r+a)/2) + i·sgn(b)*sqrt((r-a)/2)
	workPrec := z.setPrec(x) + _W

	// r = |x|
	r := x.Abs(newFloat(workPrec))
	si := x.Imag.Sign()
	a := &x.Real
	if x == z {
		a = new(big.Float).Copy(a)
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

// Pow sets z to x^y and returns z.
func (z *Complex) Pow(x, y *Complex) *Complex {
	workPrec := z.setPrec2(x, y) + _W

	// Special case: x and y are real, meaning imaginary parts are 0
	if x.Imag.Sign() == 0 && y.Imag.Sign() == 0 {
		// If x.Real >= 0, or y.Real is an integer, the result is purely real.
		// (For negative x and non-integer y, it will panic with ErrNaN in Pow,
		// which matches the real-domain Pow behavior).
		if x.Real.Sign() >= 0 || y.Real.IsInt() {
			Pow(&z.Real, &x.Real, &y.Real)
			z.Imag.Set(zero)
			return z
		}
	}

	// x^y = exp(y * log(x))
	l := newComplex(workPrec).Log(x)
	prod := newComplex(workPrec).Mul(y, l)
	return z.Exp(prod)
}

func (x *Complex) String() string {
	return fmt.Sprint(x)
}

// formatFloat writes f to s, replacing "Inf" with "∞" for readability.
// If addSign is true, positive infinity is printed as "+∞".
func formatFloat(s fmt.State, verb rune, f *big.Float, addSign bool) {
	if f.IsInf() {
		if f.Signbit() {
			fmt.Fprint(s, "-∞")
		} else if addSign {
			fmt.Fprint(s, "+∞")
		} else {
			fmt.Fprint(s, "∞")
		}
		return
	}
	f.Format(s, verb)
}

// formatImag writes the imaginary part f to s, replacing "Inf" with "∞".
// When reZero is true and f is +Inf, the value is printed as "+∞i".
func formatImag(s fmt.State, verb rune, f *big.Float, reZero bool) {
	if f.IsInf() {
		if f.Signbit() {
			fmt.Fprint(s, "-∞i")
		} else if reZero {
			fmt.Fprint(s, "+∞i")
		} else {
			fmt.Fprint(s, "∞i")
		}
		return
	}
	f.Format(s, verb)
	fmt.Fprint(s, "i")
}

func (x *Complex) Format(s fmt.State, verb rune) {
	reZero := x.Real.Sign() == 0
	imZero := x.Imag.Sign() == 0

	if reZero && imZero {
		formatFloat(s, verb, &x.Real, false)
		return
	}

	if !reZero {
		formatFloat(s, verb, &x.Real, true)
	}

	if imZero {
		return
	}

	if !reZero && !x.Imag.Signbit() {
		fmt.Fprint(s, "+")
	}

	if x.Imag.Cmp(one) == 0 {
		fmt.Fprint(s, "i")
	} else if x.Imag.Cmp(minusOne) == 0 {
		fmt.Fprint(s, "-i")
	} else {
		formatImag(s, verb, &x.Imag, reZero)
	}
}

func (x *Complex) Prec() uint {
	return max(x.Real.Prec(), x.Imag.Prec())
}

func (x *Complex) SetPrec(prec uint) *Complex {
	x.Real.SetPrec(prec)
	x.Imag.SetPrec(prec)
	return x
}

// setPrec normalizes z's precision and returns the working precision.
// If z has zero precision, inherits x's precision.
// Always sets both parts to the same precision, resolving Real≠Imag mismatches.
func (z *Complex) setPrec(x *Complex) (prec uint) {
	if prec = z.Prec(); prec == 0 {
		prec = x.Prec()
	}
	z.SetPrec(prec)
	return
}

// setPrec2 normalizes z's precision and returns the working precision.
// If z has zero precision, inherits from the operands' maximum.
// Always sets both parts to the same precision, resolving Real≠Imag mismatches.
func (z *Complex) setPrec2(x, y *Complex) (prec uint) {
	if prec = z.Prec(); prec == 0 {
		prec = argsPrec(x, y)
	}
	z.SetPrec(prec)
	return
}

func argsPrec(x, y *Complex) uint {
	return max(x.Real.Prec(), x.Imag.Prec(), y.Real.Prec(), y.Imag.Prec())
}

func newComplex(prec uint) *Complex {
	return new(Complex).SetPrec(prec)
}
