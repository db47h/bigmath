// SPDX-License-Identifier: MIT

package bigmath_test

import (
	"fmt"
	"math/big"
	"testing"

	"github.com/db47h/bigmath"
)

func benchmarkPowInt(b *testing.B, prec uint, yInt int64) {
	x := new(big.Float).SetPrec(prec).SetUint64(2)
	y := new(big.Float).SetPrec(prec).SetInt64(yInt)
	z := new(big.Float).SetPrec(prec)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		bigmath.Pow(z, x, y)
	}
}

func BenchmarkPow(b *testing.B) {
	precisions := []uint{256, 1024, 2048}
	exponents := []int64{10, 100, 1000, 10000}

	for _, prec := range precisions {
		for _, exp := range exponents {
			b.Run(fmt.Sprintf("Prec%d/Exp%d", prec, exp), func(b *testing.B) {
				benchmarkPowInt(b, prec, exp)
			})
		}
	}
}
