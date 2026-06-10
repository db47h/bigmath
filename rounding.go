// SPDX-License-Identifier: MIT

package bigmath

// Floor sets z to the greatest integer value less than or equal to x and
// returns z.
//
// Special cases:
//
//	Floor(±Inf) = ±Inf
//	Floor(±0)   = ±0
func (z *Float) Floor(x *Float) *Float {
	m := z.Mode()
	return z.SetMode(ToNegativeInf).IntRound(x).SetMode(m)
}

// Ceil sets z to the least integer value greater than or equal to x and
// returns z.
//
// Special cases:
//
//	Ceil(±Inf) = ±Inf
//	Ceil(±0)   = ±0
func (z *Float) Ceil(x *Float) *Float {
	m := z.Mode()
	return z.SetMode(ToPositiveInf).IntRound(x).SetMode(m)
}

// Trunc sets z to the integer value nearest to x, rounding toward zero, and
// returns z.
//
// Special cases:
//
//	Trunc(±Inf) = ±Inf
//	Trunc(±0)   = ±0
func (z *Float) Trunc(x *Float) *Float {
	m := z.Mode()
	return z.SetMode(ToZero).IntRound(x).SetMode(m)
}

// IntRound sets z to the integer value nearest to x, selecting the rounding
// direction according to z's rounding mode, and returns z.
//
// If z has precision 0, it is set to x's precision.
//
// In mode-specific terms:
//
//	ToNearestEven  — round to nearest, ties to even
//	ToNearestAway  — round to nearest, ties away from zero
//	ToZero         — truncation toward zero (standard Trunc)
//	ToNegativeInf  — floor (greatest integer ≤ x)
//	ToPositiveInf  — ceil (least integer ≥ x)
//	AwayFromZero   — round away from zero
//
// IntRound is allocation-free when z == x.
//
// Special cases:
//
//	IntRound(±Inf) = ±Inf
//	IntRound(±0)   = ±0
func (z *Float) IntRound(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	z.Set(x)

	if z.IsInf() || z.Sign() == 0 {
		return z
	}

	E := z.MantExp(nil)
	if E <= 0 {
		return z.intRoundZero(E)
	}

	return z.SetPrec(uint(E)).SetPrec(prec)
}

func (z *Float) intRoundZero(exp int) *Float {
	sign := z.Signbit()
	switch z.Mode() {
	case ToNearestAway:
		if exp < 0 {
			// x < 0.5, truncate to 0
			break
		}
		fallthrough
	case AwayFromZero:
		// ± 1
		z.Set(one)
		if sign {
			z.Neg(z)
		}
		return z
	case ToNegativeInf:
		if sign {
			return z.Set(minusOne)
		}
	case ToPositiveInf:
		if !sign {
			return z.Set(one)
		}
	}
	// general case and fallthrough: truncate to 0
	return z.Set(zero)
}
