package bigmath

import (
	"fmt"
	"math"
	"testing"
)

// BenchmarkBernoulliOnly measures just Bernoulli ensure() on a cold cache.
// Before each iteration we clear the global cache so each call pays the
// full recurrence cost.
func BenchmarkBernoulliOnly(b *testing.B) {
	cases := []struct {
		prec uint
	}{
		{128},
		{256},
		{512},
		{1024},
	}
	for _, c := range cases {
		n := max(int(math.Ceil(float64(c.prec)/(4*math.Pi*0.2))), 8)
		b.Run(fmt.Sprintf("N%d", n), func(b *testing.B) {
			b.ReportAllocs()
			for i := 0; i < b.N; i++ {
				bernoulli.mu.Lock()
				bernoulli.vals = nil
				bernoulli.mu.Unlock()
				bernoulli.ensure(n)
			}
		})
	}
}

// BenchmarkGammaLoopOnly measures just the gamma loop iteration cost
// with warm Bernoulli cache
func BenchmarkGammaLoopOnly(b *testing.B) {
	cases := []struct {
		prec uint
	}{
		{128},
		{256},
		{512},
		{1024},
	}
	// Warm global cache for the worst-case precision we benchmark
	warmPrec := uint(0)
	for _, c := range cases {
		if c.prec > warmPrec {
			warmPrec = c.prec
		}
	}
	{
		bernoulli.mu.Lock()
		bernoulli.vals = nil
		bernoulli.mu.Unlock()
		warmN := max(int(math.Ceil(float64(warmPrec)/(4*math.Pi*0.2))), 8)
		bernoulli.ensure(warmN)
	}
	for _, c := range cases {
		n := max(int(math.Ceil(float64(c.prec)/(4*math.Pi*0.2))), 8)
		b.Run(fmt.Sprintf("N%d", n), func(b *testing.B) {
			b.ReportAllocs()
			x := newFloat(c.prec).SetFloat64(1e10)
			z := new(Float).SetPrec(c.prec)
			b.ResetTimer()
			for i := 0; i < b.N; i++ {
				z.Lgamma(x)
			}
		})
	}
}
