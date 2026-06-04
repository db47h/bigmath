// SPDX-License-Identifier: MIT

package bigmath

import (
	"math/bits"
	"sync"
)

// word size for the mantissa of a big.Float (big.Word)
const _W = bits.UintSize

// constProvider is an internal function type used to generate a mathematical
// constant at a specific precision.
//
// Implementation Requirement: A constProvider must return a *Float that
// has at least the requested precision. To ensure the integrity of the
// cache, a constProvider should return a new *Float instance on each
// call, or at least must never return the same pointer for different
// precision values unless that pointer already satisfies the highest
// requested precision.
type constProvider func(prec uint) *Float

// The following variables provide access to common mathematical constants.
//
// Performance Note: These use an internal caching mechanism to avoid
// recomputation. To eliminate allocation overhead, the functions return
// a direct pointer to the cached *Float.
//
// IMMUTABILITY CONTRACT: Callers MUST treat the returned *Float as
// read-only. Modifying the value, precision, or rounding mode of a
// returned constant will corrupt the cache and lead to undefined
// behavior in subsequent calculations across the entire program.
var (
	// static constants not managed by the cache
	minusOne = new(Float).SetInt64(-1)
	zero     = new(Float)
	one      = new(Float).SetUint64(1)
	two      = new(Float).SetUint64(2)
	three    = new(Float).SetUint64(3)
	ten      = new(Float).SetUint64(10)
	half     = NewFloat(0.5)

	// cached constants
	pi        = cache(computePi)
	halfPi    = cache(func(prec uint) *Float { return newFloat(prec).SetMantExp(pi(prec), -1) })
	twoOverPi = cache(func(prec uint) *Float { return newFloat(prec).Quo(two, pi(prec+2)) })
	sqrt2     = cache(func(prec uint) *Float { return newFloat(prec).Sqrt(two) })
	sqrt3     = cache(func(prec uint) *Float { return newFloat(prec).Sqrt(three) })
	ln2       = cache(func(prec uint) *Float { return newFloat(prec + _W).lnCore(two).SetPrec(prec) })
	ln10      = cache(func(prec uint) *Float { return newFloat(prec).Log(ten) })
	ln10Of2   = cache(func(prec uint) *Float { return newFloat(prec).Quo(ln2(prec+2), ln10(prec+2)) })
	log2Pi    = cache(func(prec uint) *Float {
		twoPi := newFloat(prec+_W).SetMantExp(pi(prec+_W), 1)
		return newFloat(prec).Log(twoPi)
	})
)

// cache wraps a constProvider with thread-safe memoization.
//
// It only invokes the provider when the requested precision exceeds the
// precision of the currently cached value. It returns the internal pointer
// directly to avoid the cost of copying.
func cache(fn constProvider) constProvider {
	var (
		m sync.Mutex
		v *Float
	)
	return func(prec uint) *Float {
		m.Lock()
		defer m.Unlock()
		if v != nil && v.Prec() >= prec {
			return v
		}
		v = fn(prec)
		return v
	}
}

// Pi sets z to the rounded value of π and returns z.
//
// The precision of the returned value is determined by z's precision
// (z.Prec()). If z has the default precision (0), a precision of 53 bits
// (IEEE 754 double precision) is used.
func (z *Float) Pi() *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = 53
	}
	return z.Set(pi(prec))
}

// computePi computes PI using Machin's formula: PI/4 = 4*arctan(1/5) - arctan(1/239)
func computePi(prec uint) *Float {
	workPrec := prec + _W

	p1 := newFloat(workPrec)
	p2 := newFloat(workPrec)

	// 4*arctan(1/5)
	// Use the optimized atanReciprocal for reciprocal integers.
	p1.atanReciprocal(5)
	p1.SetMantExp(p1, 2)

	// arctan(1/239)
	p2.atanReciprocal(239)

	// pi/4 = 4*arctan(1/5) - arctan(1/239)
	res := newFloat(workPrec).Sub(p1, p2)

	// pi = 4 * (pi/4)
	return res.SetMantExp(res, 2).SetPrec(prec)
}

// Ln2 sets z to the rounded value of ln(2) and returns z.
//
// The precision of the returned value is determined by z's precision
// (z.Prec()). If z has the default precision (0), a precision of 53 bits
// (IEEE 754 double precision) is used.
func (z *Float) Ln2() *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = 53
	}
	return z.Set(ln2(prec))
}

// Ln10 sets z to the rounded value of ln(10) and returns z.
//
// The precision of the returned value is determined by z's precision
// (z.Prec()). If z has the default precision (0), a precision of 53 bits
// (IEEE 754 double precision) is used.
func (z *Float) Ln10() *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = 53
	}
	return z.Set(ln10(prec))
}

// Sqrt2 sets z to the rounded value of √2 and returns z.
//
// The precision of the returned value is determined by z's precision
// (z.Prec()). If z has the default precision (0), a precision of 53 bits
// (IEEE 754 double precision) is used.
func (z *Float) Sqrt2() *Float {
	prec := z.Prec()
	if prec == 0 {
		prec = 53
	}
	return z.Set(sqrt2(prec))
}
