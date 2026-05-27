// SPDX-License-Identifier: MIT

package bigmath

import (
	"fmt"
)

// Complex represents a complex number with real and imaginary parts as Float.
type Complex struct {
	Real Float
	Imag Float
}

func (z *Complex) Copy(x *Complex) *Complex {
	if z == x {
		return z
	}
	z.Real.Copy(&x.Real)
	z.Imag.Copy(&x.Imag)
	return z
}

func (z *Complex) Set(x *Complex) *Complex {
	if z == x {
		return z
	}
	z.Real.Set(&x.Real)
	z.Imag.Set(&x.Imag)
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
	workPrec := z.setPrec2(x, y) + _W

	temp := new(Float) // temp for fma. prec will be handled by fma()
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

	temp := new(Float) // temp for fma. prec will be handled by fma()
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

// Inv sets z to the multiplicative inverse of x (1/x) and returns z.
func (z *Complex) Inv(x *Complex) *Complex {
	workPrec := z.setPrec(x) + _W
	// denom = real² + imag²
	a2 := newFloat(workPrec).Mul(&x.Real, &x.Real)
	b2 := newFloat(workPrec).Mul(&x.Imag, &x.Imag)
	denom := newFloat(workPrec).Add(a2, b2)

	// Compute real and imag parts
	z.Real.Quo(&x.Real, denom)
	z.Imag.Quo(&x.Imag, denom)
	z.Imag.Neg(&z.Imag) // imag = -imag
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
func (x *Complex) Abs(z *Float) *Float {
	return z.Hypot(&x.Real, &x.Imag)
}

// Arg sets z to the rounded value of arg(x) and returns z.
func (x *Complex) Arg(z *Float) *Float {
	return z.Atan2(&x.Imag, &x.Real)
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

func (x *Complex) String() string {
	return fmt.Sprint(x)
}

// formatFloat writes f to s, replacing "Inf" with "∞" for readability.
// If addSign is true, positive infinity is printed as "+∞".
func formatFloat(s fmt.State, verb rune, f *Float, addSign bool) {
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
func formatImag(s fmt.State, verb rune, f *Float, reZero bool) {
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
