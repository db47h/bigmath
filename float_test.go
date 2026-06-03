// SPDX-License-Identifier: MIT

package bigmath_test

import (
	"math"
	"testing"

	"github.com/db47h/bigmath"
)

// Helper functions for constructing test Float values.

func floatVal(f float64) *bigmath.Float {
	return new(bigmath.Float).SetFloat64(f)
}

func neg(f *bigmath.Float) *bigmath.Float { return new(bigmath.Float).Neg(f) }

var (
	negInf  = new(bigmath.Float).SetInf(true)
	inf     = new(bigmath.Float).SetInf(false)
	zero    = new(bigmath.Float)
	negZero = new(bigmath.Float).SetFloat64(math.Copysign(0, -1))
	one     = floatVal(1)
	three   = floatVal(3)
	four    = floatVal(4)
	five    = floatVal(5)
)

// parseFloat parses a decimal string into a 128-bit precision Float.
func parseFloat(s string) *bigmath.Float {
	f := new(bigmath.Float).SetPrec(128)
	_, _, err := f.Parse(s, 0)
	if err != nil {
		panic("parseFloat: " + err.Error())
	}
	return f
}

func TestFloat_AbsCmp(t *testing.T) {
	tests := []struct {
		name string
		x, y *bigmath.Float
		want int
	}{
		// ── Infinity (paths 1-2) ──
		{"Inf_Inf", inf, inf, 0},
		{"Inf_NegInf", inf, negInf, 0},
		{"NegInf_NegInf", negInf, negInf, 0},
		{"Inf_1", inf, one, 1},
		{"1_Inf", one, inf, -1},
		{"negInf_1", inf, one, 1},
		{"1_negInf", one, inf, -1},

		// ── Zero (paths 3-4) ──
		{"0_0", zero, zero, 0},
		{"Neg0_0", negZero, zero, 0},
		{"0_1", zero, one, -1},
		{"1_0", one, zero, 1},

		// ── Different exponents (path 5: ex > ey / ex < ey) ──
		{"1e100_1", parseFloat("1e100"), one, 1},
		{"1_1e100", one, parseFloat("1e100"), -1},
		{"Neg1e100_1", parseFloat("-1e100"), one, 1},
		{"1_Neg1e100", one, parseFloat("-1e100"), -1},
		// 3 (exp=2) vs 5 (exp=3): different exponents.
		{"3_5", three, five, -1},
		{"5_3", five, three, 1},

		// ── Same exponent, same sign (path 6: sx == sy) ──
		// 4 and 5 both have MantExp exponent 3 (both in [4,8)).
		{"4_5", four, five, -1},
		{"5_4", five, four, 1},
		{"4_4", four, four, 0},
		// Both negative, same exponent → y.Cmp(x) branch.
		{"Neg4_Neg5", neg(four), neg(five), -1},
		{"Neg5_Neg4", neg(five), neg(four), 1},

		// ── Same exponent, opposite signs (path 7: temporary Float + Cmp) ──
		// -3 and 3: both have exponent 2. sx<0 → t.Neg(x).Cmp(y).
		{"Neg3_3", neg(three), three, 0},
		// -5 and 4: both have exponent 3. sx<0 → t.Neg(x).Cmp(y).
		{"Neg5_4", neg(five), four, 1},
		// 4 and -5: both have exponent 3. sx>0 → x.Cmp(t.Neg(y)).
		{"4_Neg5", four, neg(five), -1},

		// ── Edge cases: Inf/zero combinations ──
		{"0_Inf", zero, inf, -1},
		{"0_NegInf", zero, negInf, -1},
		{"Inf_0", inf, zero, 1},
		{"NegInf_0", negInf, zero, 1},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := tt.x.AbsCmp(tt.y)
			if got != tt.want {
				t.Errorf("AbsCmp(%v, %v) = %d, want %d", tt.x, tt.y, got, tt.want)
			}
		})
	}
}
