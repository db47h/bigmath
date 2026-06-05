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
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.Sign() == 0 || x.IsInf() {
		return z.Set(x)
	}

	E := x.MantExp(nil)
	if E <= 0 {
		if x.Sign() < 0 {
			return z.Set(minusOne)
		}
		return z.Set(zero)
	}

	// Ensure work precision covers E bits.
	if x.IsInt() {
		return z.Set(x)
	}

	m := z.Mode()
	z.SetMode(ToNegativeInf).SetPrec(uint(E))
	z.Set(x)
	return z.SetMode(m).SetPrec(prec)
}

// Ceil sets z to the least integer value greater than or equal to x and
// returns z.
//
// Special cases:
//
//	Ceil(±Inf) = ±Inf
//	Ceil(±0)   = ±0
func (z *Float) Ceil(x *Float) *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = x.Prec()
		z.SetPrec(prec)
	}

	if x.Sign() == 0 || x.IsInf() {
		return z.Set(x)
	}

	E := x.MantExp(nil)
	if E <= 0 {
		if x.Sign() < 0 {
			return z.Set(zero)
		}
		return z.Set(one)
	}

	if x.IsInt() {
		return z.Set(x)
	}

	m := z.Mode()
	z.SetMode(ToPositiveInf).SetPrec(uint(E))
	z.Set(x)

	return z.SetMode(m).SetPrec(prec)
}
