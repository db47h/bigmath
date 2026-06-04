// SPDX-License-Identifier: MIT

package bigmath_test

import (
	"crypto/rand"
	"fmt"
	"testing"

	"github.com/db47h/bigmath"
)

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
		x := new(bigmath.Float).SetPrec(prec).SetInt(xInt)

		for _, exp := range exponents {
			b.Run(fmt.Sprintf("Prec%d/Exp%d", prec, exp), func(b *testing.B) {
				y := new(bigmath.Float).SetPrec(prec).SetInt64(exp)
				z := new(bigmath.Float).SetPrec(prec)
				benchmarkPowInt(b, x, y, z)
			})
		}
	}
}

func benchmarkPowInt(b *testing.B, x, y, z *bigmath.Float) {
	b.Helper()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		z.Pow(x, y)
	}
}

func BenchmarkFloat_LGamma_stirling(b *testing.B) {
	prec := uint(256)
	x := new(bigmath.Float).SetFloat64(10.5)
	z := new(bigmath.Float).SetPrec(prec)
	z.Lgamma(x)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		z.Lgamma(x)
	}
}

func BenchmarkFloat_LGamma_reflect(b *testing.B) {
	prec := uint(256)
	x := new(bigmath.Float).SetFloat64(-10.5)
	z := new(bigmath.Float).SetPrec(prec)
	z.Lgamma(x)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		z.Lgamma(x)
	}
}
