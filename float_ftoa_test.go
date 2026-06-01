package bigmath

import (
	"math/big"
	"testing"
)

func TestQuickCompare(t *testing.T) {
	f := new(Float).SetPrec(53)
	bf := new(big.Float).SetPrec(53)

	for _, val := range []string{"0", "-0", "1.0", "3.14", "-3.14", "Inf", "-Inf", "1e100", "1e-100"} {
		f.SetString(val)
		bf.SetString(val)

		t.Run(val, func(t *testing.T) {
			for _, fmtChar := range []byte{'e', 'E', 'f', 'g', 'G'} {
				for _, prec := range []int{-1, 0, 5, 10} {
					got := f.Text(fmtChar, prec)
					want := bf.Text(fmtChar, prec)
					if got != want {
						t.Errorf("%s Text('%c', %d)\n  got:  %q\n  want: %q", val, fmtChar, prec, got, want)
					}
				}
			}
			// String()
			got := f.String()
			want := bf.String()
			if got != want {
				t.Errorf("%s String()\n  got:  %q\n  want: %q", val, got, want)
			}
		})
	}
}

func TestQuickLarge(t *testing.T) {
	f := new(Float).SetPrec(53)
	bf := new(big.Float).SetPrec(53)

	for _, val := range []string{"1e100", "1e1000", "1e10000", "1e-100", "1e-1000", "1e-10000"} {
		f.SetString(val)
		bf.SetString(val)

		for _, fmtChar := range []byte{'e', 'E', 'g', 'G'} {
			for _, prec := range []int{-1, 0, 5, 10} {
				got := f.Text(fmtChar, prec)
				want := bf.Text(fmtChar, prec)
				if got != want {
					t.Errorf("%s Text('%c', %d)\n  got:  %q\n  want: %q", val, fmtChar, prec, got, want)
				}
			}
		}
	}
}

var switchExp = 8000

func BenchmarkBigFloatString(b *testing.B) {
	t := new(big.Float).SetPrec(128).SetFloat64(0.5)
	t.SetMantExp(t, switchExp)
	for b.Loop() {
		_ = t.Text('g', -1)
	}
}

func BenchmarkFloatString(b *testing.B) {
	t := new(Float).SetPrec(128).SetFloat64(0.5)
	t.SetMantExp(t, switchExp)
	for b.Loop() {
		_ = t.Text('g', -1)
	}
}
