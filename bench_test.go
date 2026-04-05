// SPDX-License-Identifier: MIT

package bigmath_test

import (
	"crypto/rand"
	"fmt"
	"math/big"
	"testing"

	"github.com/db47h/bigmath"
)

func benchmarkPowInt(b *testing.B, x, y, z *big.Float) {
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		bigmath.Pow(z, x, y)
	}
}

func BenchmarkPow(b *testing.B) {
	precisions := []uint{256, 1024, 2048}
	exponents := []int64{10, 100, 1000, 10000}

	for _, prec := range precisions {
		// Pre-generate x once for each precision.
		// Use a random prime to ensure it's not a power of 2 and has full precision.
		xInt, err := rand.Prime(rand.Reader, int(prec))
		if err != nil {
			b.Fatalf("Failed to generate random prime: %v", err)
		}
		x := new(big.Float).SetPrec(prec).SetInt(xInt)

		for _, exp := range exponents {
			b.Run(fmt.Sprintf("Prec%d/Exp%d", prec, exp), func(b *testing.B) {
				y := new(big.Float).SetPrec(prec).SetInt64(exp)
				z := new(big.Float).SetPrec(prec)
				benchmarkPowInt(b, x, y, z)
			})
		}
	}
}
