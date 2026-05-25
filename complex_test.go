// SPDX-License-Identifier: MIT

package bigmath_test

import (
	"fmt"
	"math"
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
				t.Errorf("%s failed: got %v, want %v", tt.name, got, want)
			}
		})
	}
}

func TestComplex_Format(t *testing.T) {
	prec := uint(64)

	type testCase struct {
		name   string
		re, im float64
		format string
		want   string
	}

	// full grid: real × imag = {-Inf, -5, 0, 3, +Inf}
	reVals := []struct {
		v     float64
		label string
	}{
		{math.Inf(-1), "-Inf"},
		{-5, "neg"},
		{0, "zero"},
		{3, "pos"},
		{math.Inf(1), "+Inf"},
	}
	imVals := []struct {
		v     float64
		label string
	}{
		{math.Inf(-1), "-Inf"},
		{-7, "neg"},
		{0, "zero"},
		{4, "pos"},
		{math.Inf(1), "+Inf"},
	}

	// expected formats for %v (computed manually)
	// real row -> re\im col
	expectedGrid := [5][5]string{
		// im: -Inf      neg      zero    pos       +Inf
		{"-∞-∞i", "-∞-7i", "-∞", "-∞+4i", "-∞+∞i"}, // re: -Inf
		{"-5-∞i", "-5-7i", "-5", "-5+4i", "-5+∞i"}, // re: -5
		{"-∞i", "-7i", "0", "4i", "+∞i"},           // re: 0
		{"3-∞i", "3-7i", "3", "3+4i", "3+∞i"},      // re: 3
		{"+∞-∞i", "+∞-7i", "+∞", "+∞+4i", "+∞+∞i"}, // re: +Inf
	}

	var tests []testCase
	for ri, rv := range reVals {
		for ii, iv := range imVals {
			name := rv.label + "_" + iv.label
			tests = append(tests, testCase{
				name: name,
				re:   rv.v,
				im:   iv.v,
				want: expectedGrid[ri][ii],
			})
		}
	}

	// Special cases from user specification
	specialCases := []testCase{
		// User-provided examples
		{name: "neg2", re: -2, im: 0, want: "-2"},
		{name: "pos54", re: 54, im: 0, want: "54"},
		{name: "only_imag_42", re: 0, im: 42, want: "42i"},
		{name: "only_imag_neg1", re: 0, im: -1, want: "-i"},
		{name: "only_imag_1", re: 0, im: 1, want: "i"},
		{name: "pos1_pos1", re: 1, im: 1, want: "1+i"},
		{name: "pos1_neg1", re: 1, im: -1, want: "1-i"},
		{name: "pos1_pos2", re: 1, im: 2, want: "1+2i"},
		{name: "pos1_neg2", re: 1, im: -2, want: "1-2i"},

		// More boundary cases
		{name: "neg1_zero", re: -1, im: 0, want: "-1"},
		{name: "neg1_pos1", re: -1, im: 1, want: "-1+i"},
		{name: "neg1_neg1", re: -1, im: -1, want: "-1-i"},
		{name: "pos2_pos1", re: 2, im: 1, want: "2+i"},
		{name: "neg2_pos1", re: -2, im: 1, want: "-2+i"},
		{name: "neg2_neg1", re: -2, im: -1, want: "-2-i"},
		{name: "zero_neg2", re: 0, im: -2, want: "-2i"},

		// User-specified inf formatting examples
		{name: "inf_plus_3i", re: math.Inf(1), im: 3, want: "+∞+3i"},
		{name: "1_plus_inf_i", re: 1, im: math.Inf(1), want: "1+∞i"},
	}
	tests = append(tests, specialCases...)

	// %.2f format tests (legacy format test)
	f2Tests := []testCase{
		{name: "fmt_2f_neg_imag", re: 1.23, im: -4.56, format: "%.2f", want: "1.23-4.56i"},
		{name: "fmt_2f_pos_imag", re: 1.23, im: 4.56, format: "%.2f", want: "1.23+4.56i"},
		{name: "fmt_2f_only_imag", re: 0, im: 4.56, format: "%.2f", want: "4.56i"},
		{name: "fmt_2f_only_real", re: 1.23, im: 0, format: "%.2f", want: "1.23"},
		{name: "fmt_2f_neg_imag_only", re: 0, im: -4.56, format: "%.2f", want: "-4.56i"},
		{name: "fmt_2f_one_imag", re: 1, im: 1, format: "%.2f", want: "1.00+i"},
		{name: "fmt_2f_neg_one_imag", re: 1, im: -1, format: "%.2f", want: "1.00-i"},
	}
	tests = append(tests, f2Tests...)

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			c := newComplex(tt.re, tt.im, prec)
			var got string
			if tt.format != "" {
				got = fmt.Sprintf(tt.format, c)
			} else {
				got = fmt.Sprint(c)
			}
			if got != tt.want {
				t.Errorf("Format(%s,%v,%v): got %q, want %q", tt.format, tt.re, tt.im, got, tt.want)
			}
		})
	}
}

func TestComplex_Exp_EdgeCases(t *testing.T) {
	prec := uint(53)

	t.Run("large positive real triggers overflow", func(t *testing.T) {
		z := newComplex(0, 0, prec)
		x := newComplex(1e10, 0, prec)
		got := z.Exp(x)
		if !got.Real.IsInf() || got.Real.Signbit() {
			t.Errorf("Real: expected +Inf, got %v", &got.Real)
		}
		if got.Imag.Sign() != 0 {
			t.Errorf("Imag: expected +Inf, got %v", &got.Imag)
		}
	})

	t.Run("large positive real with imag triggers overflow", func(t *testing.T) {
		z := newComplex(0, 0, prec)
		x := newComplex(1e10, 1, prec)
		got := z.Exp(x)
		if !got.Real.IsInf() || got.Real.Signbit() {
			t.Errorf("Real: expected +Inf, got %v", &got.Real)
		}
		if !got.Imag.IsInf() || got.Imag.Signbit() {
			t.Errorf("Imag: expected +Inf, got %v", &got.Imag)
		}
	})

	t.Run("large negative real triggers underflow", func(t *testing.T) {
		z := newComplex(0, 0, prec)
		x := newComplex(-1e10, 0, prec)
		got := z.Exp(x)
		if got.Real.Sign() != 0 {
			t.Errorf("Real: expected 0, got %v", &got.Real)
		}
		if got.Imag.Sign() != 0 {
			t.Errorf("Imag: expected 0, got %v", &got.Imag)
		}
	})

	t.Run("large negative real with imag triggers underflow", func(t *testing.T) {
		z := newComplex(0, 0, prec)
		x := newComplex(-1e10, 1.5, prec)
		got := z.Exp(x)
		if got.Real.Sign() != 0 {
			t.Errorf("Real: expected 0, got %v", &got.Real)
		}
		if got.Imag.Sign() != 0 {
			t.Errorf("Imag: expected 0, got %v", &got.Imag)
		}
	})

	t.Run("real +Inf", func(t *testing.T) {
		z := newComplex(0, 0, prec)
		x := newComplex(math.Inf(1), 0, prec)
		got := z.Exp(x)
		if !got.Real.IsInf() || got.Real.Signbit() {
			t.Errorf("Real: expected +Inf, got %v", &got.Real)
		}
		if z.Imag.Sign() != 0 {
			t.Errorf("Imag: expected +Inf, got %v", &got.Imag)
		}
	})

	t.Run("real +Inf with imag", func(t *testing.T) {
		z := newComplex(0, 0, prec)
		x := newComplex(math.Inf(1), 1, prec)
		got := z.Exp(x)
		if !got.Real.IsInf() || got.Real.Signbit() {
			t.Errorf("Real: expected +Inf, got %v", &got.Real)
		}
		if !got.Imag.IsInf() || got.Imag.Signbit() {
			t.Errorf("Imag: expected +Inf, got %v", &got.Imag)
		}
	})

	t.Run("real -Inf", func(t *testing.T) {
		z := newComplex(0, 0, prec)
		x := newComplex(math.Inf(-1), 0, prec)
		got := z.Exp(x)
		if got.Real.Sign() != 0 {
			t.Errorf("Real: expected 0, got %v", &got.Real)
		}
		if got.Imag.Sign() != 0 {
			t.Errorf("Imag: expected 0, got %v", &got.Imag)
		}
	})

	t.Run("real -Inf with imag", func(t *testing.T) {
		z := newComplex(0, 0, prec)
		x := newComplex(math.Inf(-1), 1.5, prec)
		got := z.Exp(x)
		if got.Real.Sign() != 0 {
			t.Errorf("Real: expected 0, got %v", &got.Real)
		}
		if got.Imag.Sign() != 0 {
			t.Errorf("Imag: expected 0, got %v", &got.Imag)
		}
	})

	t.Run("zero", func(t *testing.T) {
		z := newComplex(0, 0, prec)
		x := newComplex(0, 0, prec)
		got := z.Exp(x)
		// e^0 = 1
		float64RE, _ := got.Real.Float64()
		float64IM, _ := got.Imag.Float64()
		want := cmplx.Exp(complex(0, 0))
		if cmplx.Abs(complex(float64RE, float64IM)-want) > 1e-15 {
			t.Errorf("Exp(0): got %v, want %v", got, want)
		}
	})

	t.Run("normal point", func(t *testing.T) {
		z := newComplex(0, 0, prec)
		x := newComplex(0.5, 0.7, prec)
		got := z.Exp(x)
		float64RE, _ := got.Real.Float64()
		float64IM, _ := got.Imag.Float64()
		want := cmplx.Exp(complex(0.5, 0.7))
		if cmplx.Abs(complex(float64RE, float64IM)-want) > 1e-15 {
			t.Errorf("Exp(0.5+0.7i): got %v, want %v", got, want)
		}
	})

	// --- Panic tests for infinite imaginary part ---
	t.Run("imag +Inf panics", func(t *testing.T) {
		defer func() {
			r := recover()
			if r == nil {
				t.Error("expected panic, got none")
			} else if _, ok := r.(bigmath.ErrNaN); !ok {
				t.Errorf("expected ErrNaN panic, got %T(%v)", r, r)
			}
		}()
		z := newComplex(0, 0, prec)
		x := newComplex(0, math.Inf(1), prec)
		z.Exp(x)
	})

	t.Run("imag -Inf panics", func(t *testing.T) {
		defer func() {
			r := recover()
			if r == nil {
				t.Error("expected panic, got none")
			} else if _, ok := r.(bigmath.ErrNaN); !ok {
				t.Errorf("expected ErrNaN panic, got %T(%v)", r, r)
			}
		}()
		z := newComplex(0, 0, prec)
		x := newComplex(0, math.Inf(-1), prec)
		z.Exp(x)
	})

	t.Run("real + imag +Inf panics", func(t *testing.T) {
		defer func() {
			r := recover()
			if r == nil {
				t.Error("expected panic, got none")
			} else if _, ok := r.(bigmath.ErrNaN); !ok {
				t.Errorf("expected ErrNaN panic, got %T(%v)", r, r)
			}
		}()
		z := newComplex(0, 0, prec)
		x := newComplex(1, math.Inf(1), prec)
		z.Exp(x)
	})

	// --- 256-bit sanity tests for large finite imaginary ---
	t.Run("large imag 1e6 at 256-bit", func(t *testing.T) {
		prec := uint(256)
		z := newComplex(0, 0, prec)
		x := newComplex(0, 1e6, prec)
		got := z.Exp(x)
		if got.Real.IsInf() {
			t.Errorf("Real: expected finite, got Inf")
		}
		if got.Imag.IsInf() {
			t.Errorf("Imag: expected finite, got Inf")
		}
	})

	t.Run("large imag 1e20 at 256-bit", func(t *testing.T) {
		prec := uint(256)
		z := newComplex(0, 0, prec)
		x := newComplex(0, 1e20, prec)
		got := z.Exp(x)
		if got.Real.IsInf() {
			t.Errorf("Real: expected finite, got Inf")
		}
		if got.Imag.IsInf() {
			t.Errorf("Imag: expected finite, got Inf")
		}
	})

	t.Run("medium real + large imag at 256-bit", func(t *testing.T) {
		prec := uint(256)
		z := newComplex(0, 0, prec)
		x := newComplex(1, 1e6, prec)
		got := z.Exp(x)
		if got.Real.IsInf() {
			t.Errorf("Real: expected finite, got Inf")
		}
		if got.Imag.IsInf() {
			t.Errorf("Imag: expected finite, got Inf")
		}
	})
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
