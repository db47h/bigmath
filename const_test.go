// SPDX-License-Identifier: MIT

package bigmath

import (
	"math/big"
	"testing"
)

// TestCacheBehavior verifies that the cache only recomputes when requested
// precision increases.
func TestCacheBehavior(t *testing.T) {
	var count int
	cached := cache(func(prec uint) *big.Float {
		count++
		return new(big.Float).SetPrec(prec).SetFloat64(3.14)
	})

	// Initial call
	v1 := cached(64)
	if count != 1 {
		t.Errorf("expected 1 call, got %d", count)
	}

	// Call with lower precision, should not recompute
	v2 := cached(32)
	if count != 1 {
		t.Errorf("expected 1 call, got %d", count)
	}
	if v1 != v2 {
		t.Errorf("expected same pointer for lower precision")
	}

	// Call with higher precision, should recompute
	v3 := cached(128)
	if count != 2 {
		t.Errorf("expected 2 calls, got %d", count)
	}
	if v1 == v3 {
		t.Errorf("expected different pointer for higher precision")
	}
}
