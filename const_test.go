// SPDX-License-Identifier: MIT

package bigmath

import (
	"testing"
)

// TestCacheBehavior verifies that the cache only recomputes when requested
// precision increases.
func TestCacheBehavior(t *testing.T) {
	var count int
	cached := &constCache{fn: func(prec uint) *Float {
		count++
		return new(Float).SetPrec(prec).SetFloat64(3.14)
	}}

	// Initial call
	v1 := cached.get(64)
	if count != 1 {
		t.Errorf("expected 1 call, got %d", count)
	}

	// Call with lower precision, should not recompute
	v2 := cached.get(32)
	if count != 1 {
		t.Errorf("expected 1 call, got %d", count)
	}
	if v1 == v2 {
		t.Errorf("expected different pointer for lower precision")
	}

	// Call with higher precision, should recompute
	v3 := cached.get(128)
	if count != 2 {
		t.Errorf("expected 2 calls, got %d", count)
	}
	if v1 == v3 {
		t.Errorf("expected different pointer for higher precision")
	}
	// Call with higher precision, should recompute
	v4 := cached.get(128)
	if count != 2 {
		t.Errorf("expected 2 calls, got %d", count)
	}
	if v3 != v4 {
		t.Errorf("expected same pointer for same precision")
	}
}
