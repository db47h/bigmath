// SPDX-License-Identifier: MIT

package bigmath_test

import (
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
		t.Errorf("Add failed: got %v + %vi, want 4 + 6i", z.Real, z.Imag)
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
		t.Errorf("Mul failed: got %v + %vi, want -5 + 10i", z.Real, z.Imag)
	}
}
