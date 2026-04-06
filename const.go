// SPDX-License-Identifier: MIT

package bigmath

import (
	"math/big"
	"math/bits"
	"sync"
)

// word size for the mantissa of a big.Float (big.Word)
const _W = bits.UintSize

// constProvider is an internal function type used to generate a mathematical
// constant at a specific precision.
//
// Implementation Requirement: A constProvider must return a *big.Float that
// has at least the requested precision. To ensure the integrity of the
// cache, a constProvider should return a new *big.Float instance on each
// call, or at least must never return the same pointer for different
// precision values unless that pointer already satisfies the highest
// requested precision.
type constProvider func(prec uint) *big.Float

// The following variables provide access to common mathematical constants.
//
// Performance Note: These use an internal caching mechanism to avoid
// recomputation. To eliminate allocation overhead, the functions return
// a direct pointer to the cached *big.Float.
//
// IMMUTABILITY CONTRACT: Callers MUST treat the returned *big.Float as
// read-only. Modifying the value, precision, or rounding mode of a
// returned constant will corrupt the cache and lead to undefined
// behavior in subsequent calculations across the entire program.
var (
	// static constants not managed by the cache
	zero                 = new(big.Float)
	one                  = new(big.Float).SetUint64(1)
	two                  = new(big.Float).SetUint64(2)
	five                 = new(big.Float).SetUint64(5)
	ten                  = new(big.Float).SetUint64(10)
	twoHundredThirtyNine = new(big.Float).SetUint64(239)

	minusOne = new(big.Float).SetInt64(-1)

	// cached constants
	sqrt2 = cache(func(prec uint) *big.Float { return newFloat(prec).Sqrt(two) })
	ln2   = cache(func(prec uint) *big.Float { return computeLn(newFloat(prec+_W), two).SetPrec(prec) })
	ln10  = cache(func(prec uint) *big.Float { return Log(newFloat(prec), ten) })
	pi    = cache(computePi)
)

// cache wraps a constProvider with thread-safe memoization.
//
// It only invokes the provider when the requested precision exceeds the
// precision of the currently cached value. It returns the internal pointer
// directly to avoid the cost of copying.
func cache(fn constProvider) constProvider {
	var (
		m sync.Mutex
		v *big.Float
	)
	return func(prec uint) *big.Float {
		m.Lock()
		defer m.Unlock()
		if v != nil && v.Prec() >= prec {
			return v
		}
		v = fn(prec)
		return v
	}
}
