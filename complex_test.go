// SPDX-License-Identifier: MIT

package bigmath_test

import (
	"fmt"
	"math/big"
	"math/cmplx"
	"testing"

	"github.com/db47h/bigmath"
)

func newComplex(re, im float64, prec uint) *bigmath.Complex {
	c := &bigmath.Complex{}
	c.Real.SetPrec(prec).SetFloat64(re)
	c.Imag.SetPrec(prec).SetFloat64(im)
	return c
}

func TestComplex_Add(t *testing.T) {
	prec := uint(64)
	x := newComplex(1, 2, prec)
	y := newComplex(3, 4, prec)
	z := newComplex(0, 0, prec)
	want := newComplex(4, 6, prec)
	z.Add(x, y)
	if !z.Equals(want) {
		t.Errorf("Add failed: got %v + %vi, want 4 + 6i", &z.Real, &z.Imag)
	}
}

func TestComplex_Mul(t *testing.T) {
	prec := uint(64)
	x := newComplex(1, 2, prec)
	y := newComplex(3, 4, prec)
	z := newComplex(0, 0, prec)
	want := newComplex(-5, 10, prec)
	z.Mul(x, y)
	if !z.Equals(want) {
		t.Errorf("Mul failed: got %v + %vi, want -5 + 10i", &z.Real, &z.Imag)
	}
}

func TestComplex_Identities(t *testing.T) {
	prec := uint(256)
	x := newComplex(0.5, 0.7, prec)
	oneC := newComplex(1, 0, prec)

	// sin^2(z) + cos^2(z) = 1
	s := new(bigmath.Complex).Sin(x)
	c := new(bigmath.Complex).Cos(x)
	s2 := new(bigmath.Complex).Mul(s, s)
	c2 := new(bigmath.Complex).Mul(c, c)
	res := new(bigmath.Complex).Add(s2, c2)

	// Check if the difference is small.
	diff := new(big.Float).SetPrec(prec).Sub(&res.Real, &oneC.Real)
	diff.Abs(diff)

	limit := new(big.Float).SetPrec(prec).SetUint64(1)
	limit.SetMantExp(limit, 1-int(prec))

	if diff.Cmp(limit) > 0 {
		t.Errorf("sin^2(z) + cos^2(z) != 1: got %v, diff %v", res.Real.Text('g', 10), diff.Text('g', 10))
	}
	if res.Imag.Sign() != 0 && res.Imag.MantExp(nil) > 1-int(prec) {
		t.Errorf("sin^2(z) + cos^2(z) has non-zero Imag: got %v", res.Imag.Text('g', 10))
	}
}

func TestComplex_AgainstStd(t *testing.T) {
	prec := uint(53) // matches float64 precision
	re, im := 0.5, 0.7
	x := newComplex(re, im, prec)
	xc := complex(re, im)

	tests := []struct {
		name string
		f    func(z, x *bigmath.Complex) *bigmath.Complex
		std  func(x complex128) complex128
	}{
		{"Sin", (*bigmath.Complex).Sin, cmplx.Sin},
		{"Cos", (*bigmath.Complex).Cos, cmplx.Cos},
		{"Tan", (*bigmath.Complex).Tan, cmplx.Tan},
		{"Sinh", (*bigmath.Complex).Sinh, cmplx.Sinh},
		{"Cosh", (*bigmath.Complex).Cosh, cmplx.Cosh},
		{"Tanh", (*bigmath.Complex).Tanh, cmplx.Tanh},
		{"Asin", (*bigmath.Complex).Asin, cmplx.Asin},
		{"Acos", (*bigmath.Complex).Acos, cmplx.Acos},
		{"Atan", (*bigmath.Complex).Atan, cmplx.Atan},
		{"Asinh", (*bigmath.Complex).Asinh, cmplx.Asinh},
		{"Acosh", (*bigmath.Complex).Acosh, cmplx.Acosh},
		{"Atanh", (*bigmath.Complex).Atanh, cmplx.Atanh},
		{"Sqrt", (*bigmath.Complex).Sqrt, cmplx.Sqrt},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := new(bigmath.Complex)
			tt.f(got, x)

			want := tt.std(xc)

			gotRe, _ := got.Real.Float64()
			gotIm, _ := got.Imag.Float64()

			if cmplx.Abs(complex(gotRe, gotIm)-want) > 1e-15 {
				t.Errorf("%s failed: got %v+%vi, want %v", tt.name, gotRe, gotIm, want)
			}
		})
	}
}

func TestComplex_Format(t *testing.T) {
	c := newComplex(1.23, -4.56, 64)
	got := fmt.Sprintf("%.2f", c)
	want := "1.23-4.56i"
	if got != want {
		t.Errorf("Format failed: got %q, want %q", got, want)
	}

	c2 := newComplex(1.23, 4.56, 64)
	got2 := fmt.Sprintf("%.2f", c2)
	want2 := "1.23+4.56i"
	if got2 != want2 {
		t.Errorf("Format failed: got %q, want %q", got2, want2)
	}
}

func TestComplex_Aliasing(t *testing.T) {
	prec := uint(64)
	x := newComplex(1, 2, prec)
	y := newComplex(3, 4, prec)

	_ = x // x used only for aliasing tests below

	// Mul
	z := newComplex(1, 2, prec)
	z.Mul(z, y)
	want := newComplex(-5, 10, prec)
	if !z.Equals(want) {
		t.Errorf("Mul aliasing failed: got %v, want %v", z, want)
	}

	// Quo
	z = newComplex(1, 2, prec)
	z.Quo(z, z)
	want = newComplex(1, 0, prec)
	if !z.Equals(want) {
		t.Errorf("Quo aliasing (z=x=y) failed: got %v, want %v", z, want)
	}

	// Sin
	z = newComplex(0.5, 0.7, prec)
	wantZ := new(bigmath.Complex).Sin(z)
	z.Sin(z)
	if !z.Equals(wantZ) {
		t.Errorf("Sin aliasing failed: got %v, want %v", z, wantZ)
	}
}
