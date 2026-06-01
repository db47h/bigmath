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
	sqrt2 = cache(func(prec uint) *Float { return newFloat(prec).Sqrt(two) })
	sqrt3 = cache(func(prec uint) *Float { return newFloat(prec).Sqrt(three) })
	ln2   = cache(func(prec uint) *Float {
		return newFloat(prec + _W).lnCore(two).SetPrec(prec)
	})
	ln10      = cache(func(prec uint) *Float { return newFloat(prec).Log(ten) })
	pi        = cache(computePi)
	halfPi    = cache(func(prec uint) *Float { return newFloat(prec).SetMantExp(pi(prec), -1) })
	twoOverPi = cache(func(prec uint) *Float { return newFloat(prec).Quo(two, pi(prec+2)) })
	ln10Of2   = cache(func(prec uint) *Float { return newFloat(prec).Quo(ln2(prec+2), ln10(prec+2)) })
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
