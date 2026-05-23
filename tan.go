// SPDX-License-Identifier: MIT

package bigmath

import (
	"math/big"
)

// Tan sets z to the tangent of x and returns z.
//
// Special cases:
//
//	Tan(±0) = ±0
//	Tan(±Inf) = panic(ErrNaN)
//	Tan(π/2 + nπ) = ±Inf
func Tan(z, x *big.Float) *big.Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.IsInf() {
		panic(ErrNaN("tan of infinity"))
	}
	if x.Sign() == 0 {
		return z.Set(x)
	}

	// Flat +_W guard — same rationale as Sin/Cos. The computation path
	// (reducePi2 → sinCore/cosCore → Quo) is identical in structure.
	workPrec := prec + _W

	xVal := newFloat(workPrec).Set(x)
	neg := xVal.Signbit()
	if neg {
		xVal.Neg(xVal)
	}

	// reducePi2 handles z == x aliasing; reuse xVal as both input and output.
	quad := reducePi2(xVal, xVal)

	s := newFloat(workPrec)
	c := newFloat(workPrec)
	sinCore(s, xVal)
	cosCore(c, xVal)

	// tan(x) after reduction to [0, π/2):
	//   Q0: tan = sinR / cosR   → s / c
	//   Q1: tan = cosR / -sinR  → -c / s
	//   Q2: tan = -sinR / -cosR → s / c
	//   Q3: tan = -cosR / sinR  → -c / s
	// The formula choice determines the sign; no extra quadrant sign flip needed.
	if quad&1 == 0 {
		// Q0, Q2: tan = s / c
		if c.Sign() == 0 {
			z.SetInf(false)
			if neg {
				z.Neg(z)
			}
			return z
		}
		z.Quo(s, c)
	} else {
		// Q1, Q3: tan = -c / s
		if s.Sign() == 0 {
			z.SetInf(true)
			if neg {
				z.Neg(z)
			}
			return z
		}
		c.Neg(c)
		z.Quo(c, s)
	}

	if neg {
		z.Neg(z)
	}
	return z
}
