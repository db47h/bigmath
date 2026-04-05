// SPDX-License-Identifier: MIT

package bigmath_test

import (
	"fmt"
	"math/big"
	"math/rand/v2"
	"testing"

	"github.com/db47h/bigmath"
)

func benchmarkPowInt(b *testing.B, x, y, z *big.Float) {
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		bigmath.Pow(z, x, y)
	}
}

func randomBigInt(rng *rand.Rand, prec uint) *big.Int {
	bytes := make([]byte, (prec+7)/8)
	for i := range bytes {
		bytes[i] = byte(rng.Uint32())
	}
	n := new(big.Int).SetBytes(bytes)
	// Ensure it has exactly 'prec' bits
	n.SetBit(n, int(prec-1), 1)
	// Ensure it's not a power of 2 by making it odd
	n.SetBit(n, 0, 1)
	return n
}

func BenchmarkPow(b *testing.B) {
	precisions := []uint{256, 1024, 2048}
	exponents := []int64{10, 100, 1000, 10000}

	rng := rand.New(rand.NewPCG(42, 42))

	for _, prec := range precisions {
		// Pre-generate x once for each precision
		xInt := randomBigInt(rng, prec)
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
