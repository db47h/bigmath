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
	pi        = &constCache{fn: computePi}
	halfPi    = &constCache{fn: func(prec uint) *Float { return newFloat(prec).SetMantExp(pi.get(prec), -1) }}
	twoOverPi = &constCache{fn: func(prec uint) *Float { return newFloat(prec).Quo(two, pi.get(prec)) }}
	sqrt2     = &constCache{fn: func(prec uint) *Float { return newFloat(prec).Sqrt(two) }}
	sqrt3     = &constCache{fn: func(prec uint) *Float { return newFloat(prec).Sqrt(three) }}
	ln2       = &constCache{fn: func(prec uint) *Float { return newFloat(prec + _W).lnCore(two).SetPrec(prec) }}
	ln10      = &constCache{fn: func(prec uint) *Float { return newFloat(prec).Log(ten) }}
	ln10Of2   = &constCache{fn: func(prec uint) *Float { return newFloat(prec).Quo(ln2.get(prec+2), ln10.get(prec+2)) }}
	logPi     = &constCache{fn: func(prec uint) *Float { return newFloat(prec).Log(pi.get(prec + _W)) }}
	log2Pi    = &constCache{fn: func(prec uint) *Float {
		twoPi := newFloat(prec+_W).SetMantExp(pi.get(prec+_W), 1)
		return newFloat(prec).Log(twoPi)
	}}
)

// constCache is thread-safe a [*Float] cache.
//
// It only invokes the constProvider function when the requested precision exceeds the
// precision of the currently cached value.
type constCache struct {
	m  sync.Mutex
	fn constProvider
	v  *Float
}

// get returns the cached constant rounded to prec bits. The returned *Float MUST NOT be modified.
// If the precision of the cached value matches exactly the requested precision,
// it returns a pointer to the cached value without allocations.
func (c *constCache) get(prec uint) *Float {
	c.m.Lock()
	defer c.m.Unlock()
	if c.v != nil {
		switch pv := c.v.Prec(); {
		case pv > prec:
			return newFloat(prec).Set(c.v)
		case pv == prec:
			return c.v
		}
	}
	// c.v == nil || c.v.Prec() < prec
	c.v = c.fn(prec)
	return c.v
}

// get performs the operation z.Set(c.get(z.Prec())) optimizing
// out any intermediate allocation and rounding.
func (z *Float) setConst(c *constCache) *Float {
	prec := z.Prec()
	c.m.Lock()
	defer c.m.Unlock()
	if c.v != nil && c.v.Prec() >= prec {
		return z.Set(c.v)
	}
	c.v = c.fn(prec)
	return z.Set(c.v)
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
	return z.setConst(pi)
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
	return z.setConst(ln2)
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
	return z.setConst(ln10)
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
	return z.setConst(sqrt2)
}
