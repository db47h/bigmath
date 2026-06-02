package bigmath

import (
	"fmt"
	"math/big"
	"strings"
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
	for i := 0; i < b.N; i++ {
		_ = t.Text('g', -1)
	}
}

func BenchmarkFloatString(b *testing.B) {
	t := new(Float).SetPrec(128).SetFloat64(0.5)
	t.SetMantExp(t, switchExp)
	for i := 0; i < b.N; i++ {
		_ = t.Text('g', -1)
	}
}

// engineeringNotationVerifier checks that a string produced by 'n'/'N' format
// satisfies the engineering notation invariants.
func engineeringNotationVerifier(t *testing.T, s string, isUpper bool) {
	// +Inf / -Inf / Inf are valid outputs for any format
	if s == "+Inf" || s == "-Inf" || s == "Inf" {
		return
	}

	// Strip optional leading sign for parsing
	body := s
	if len(body) > 0 && (body[0] == '+' || body[0] == '-') {
		body = body[1:]
	}

	// Must contain 'e' or 'E'
	sep := byte('e')
	if isUpper {
		sep = 'E'
	}
	eIdx := strings.IndexByte(body, sep)
	if eIdx < 0 {
		t.Errorf("output %q has no exponent separator %c", s, sep)
		return
	}

	// Check mantissa format: must have 1-3 digits before decimal point
	mant := body[:eIdx]
	dotIdx := strings.IndexByte(mant, '.')
	if dotIdx >= 0 {
		intPart := mant[:dotIdx]
		if len(intPart) < 1 || len(intPart) > 3 {
			t.Errorf("output %q: integer part has %d digits, want 1-3", s, len(intPart))
		}
		for _, ch := range intPart {
			if ch < '0' || ch > '9' {
				t.Errorf("output %q: non-digit in integer part", s)
			}
		}
	} else {
		// No decimal point means prec=0, integer-only mantissa
		if len(mant) < 1 || len(mant) > 3 {
			t.Errorf("output %q: mantissa has %d digits without decimal, want 1-3", s, len(mant))
		}
	}

	// Exponent must be multiple of 3
	expStr := body[eIdx+1:]
	if len(expStr) < 2 || (expStr[0] != '+' && expStr[0] != '-') {
		t.Errorf("output %q: invalid exponent format", s)
		return
	}
	var exp int64
	if _, err := fmt.Sscanf(expStr, "%d", &exp); err != nil {
		t.Errorf("output %q: cannot parse exponent: %v", s, err)
		return
	}
	if exp%3 != 0 {
		t.Errorf("output %q: exponent %d is not a multiple of 3", s, exp)
	}
}

func TestEngineeringNotation(t *testing.T) {
	f := new(Float).SetPrec(53)

	// -------------------------------------------------------------------
	// 1. Exact output for canonical examples at default precision (prec=6)
	// -------------------------------------------------------------------
	t.Run("canonical", func(t *testing.T) {
		tests := []struct {
			val     string
			want    string // expected with 'n' (prec=6)
			wantCap string // expected with 'N' (prec=6)
		}{
			// Plan's canonical table — verified by construction
			{"12345", "12.345000e+03", "12.345000E+03"},
			{"0.00123", "1.230000e-03", "1.230000E-03"},
			{"0.000123", "123.000000e-06", "123.000000E-06"},
			{"1e-12", "1.000000e-12", "1.000000E-12"},
			{"1e-13", "100.000000e-15", "100.000000E-15"},
			{"0", "0.000000e+00", "0.000000E+00"},
			{"Inf", "+Inf", "+Inf"},
			{"-Inf", "-Inf", "-Inf"},
		}
		for _, tc := range tests {
			f.SetString(tc.val)
			got := f.Text('n', 6)
			if got != tc.want {
				t.Errorf("Text('n',6) for %q: got %q, want %q", tc.val, got, tc.want)
			}
			gotCap := f.Text('N', 6)
			if gotCap != tc.wantCap {
				t.Errorf("Text('N',6) for %q: got %q, want %q", tc.val, gotCap, tc.wantCap)
			}
		}
	})

	// -------------------------
	// 2. Zero precision (prec=0)
	// -------------------------
	t.Run("prec_zero", func(t *testing.T) {
		tests := []struct {
			val  string
			want string
		}{
			{"12345", "12e+03"},
			{"0.000123", "123e-06"},
			{"0", "0e+00"},
			{"1000", "1e+03"},
		}
		for _, tc := range tests {
			f.SetString(tc.val)
			got := f.Text('n', 0)
			if got != tc.want {
				t.Errorf("Text('n',0) for %q: got %q, want %q", tc.val, got, tc.want)
			}
		}
	})

	// --------------------------
	// 3. "n" vs "N" case difference
	// --------------------------
	t.Run("case_diff", func(t *testing.T) {
		for _, val := range []string{"12345", "0.00123", "0.000123", "1e100", "1e-13"} {
			f.SetString(val)
			lower := f.Text('n', 6)
			upper := f.Text('N', 6)
			// Same length and structure, only 'e' vs 'E' differs
			if len(lower) != len(upper) {
				t.Errorf("Text('n',6) and Text('N',6) have different lengths for %q: %q vs %q", val, lower, upper)
			}
			// Replace 'e' with 'E' and compare
			normalized := strings.ReplaceAll(lower, "e", "E")
			if normalized != upper {
				t.Errorf("Text('n',6)=%q and Text('N',6)=%q differ (normalized: %q)", lower, upper, normalized)
			}
		}
	})

	// -------------------------------------
	// 4. Structural invariants for all cases
	// -------------------------------------
	t.Run("invariants", func(t *testing.T) {
		for _, val := range []string{
			"0", "-0", "1", "12345", "3.14", "-3.14",
			"0.00123", "0.000123", "1e-12", "1e-13", "1e100", "1e1000",
			"1e-100", "1e-1000", "Inf", "-Inf",
		} {
			f.SetString(val)
			for _, fmtChar := range []byte{'n', 'N'} {
				for _, prec := range []int{-1, 0, 2, 6} {
					got := f.Text(fmtChar, prec)
					isUpper := fmtChar == 'N'
					engineeringNotationVerifier(t, got, isUpper)
				}
			}
		}
	})

	// -----------------------------------------
	// 5. Semantic equivalence: %n ≈ %e (with sufficient precision)
	// Note: at prec=0, %e produces 1 significant digit while %n produces
	// 1-3 significant digits, so direct value comparison is not meaningful.
	// We compare at prec >= 2 where both have at least 3+prec total digits.
	// Even then, %n has (shift+1) leading digits vs %e's 1 leading digit,
	// so %n is more precise. The tolerance must account for %e's rounding
	// at 1+prec significant digits.
	// -----------------------------------------
	t.Run("semantic_equivalence", func(t *testing.T) {
		for _, val := range []string{
			"0", "1", "12345", "3.14", "-3.14",
			"0.00123", "0.000123", "1e-12", "1e-13", "1e100", "1e-100",
		} {
			f.SetString(val)
			for _, prec := range []int{2, 6} {
				eStr := f.Text('e', prec)
				nStr := f.Text('n', prec)

				// Parse both back and compare values
				var eVal, nVal Float
				if _, _, err := eVal.Parse(eStr, 0); err != nil {
					t.Fatalf("cannot parse %%e output %q: %v", eStr, err)
				}
				if _, _, err := nVal.Parse(nStr, 0); err != nil {
					t.Fatalf("cannot parse %%n output %q: %v", nStr, err)
				}

				// Compute absolute difference
				diff := new(Float).Sub(&eVal, &nVal)
				diff.Abs(diff)

				// Tolerance: %e at prec P has 1+P significant digits,
				// so the ULP is 10^(exp - P) where exp is the %e exponent.
				// The %n value should be within 1 ULP of %e.
				// Parse the %e exponent from the string.
				eExp := int64(0)
				if eIdx := strings.IndexAny(eStr, "eE"); eIdx >= 0 {
					eExpStr := eStr[eIdx+1:]
					fmt.Sscanf(eExpStr, "%d", &eExp)
				}
				// %e at precision prec has 1+prec significant digits.
				// The ULP of %e output is 10^(exp - prec).
				// %n has (shift+1)+prec significant digits (shift ∈ {0,1,2}),
				// so %n is at least as precise as %e.
				// Tolerance = 10^(eExp - prec + 1) (one order above ULP to
				// account for rounding boundary cases).
				tol := new(Float)
				pow := eExp - int64(prec) + 1
				if pow >= 0 {
					tol.SetFloat64(10.0).Pow(tol, new(Float).SetInt64(pow))
				} else {
					// negative power: compute reciprocal
					tol.SetFloat64(10.0).Pow(tol, new(Float).SetInt64(-pow))
					onef := new(Float).SetFloat64(1.0)
					tol.Quo(onef, tol)
				}
				if diff.Cmp(tol) > 0 {
					t.Errorf("%%e=%q and %%n=%q for %q differ by %v > %v (eExp=%d, prec=%d)",
						eStr, nStr, val, diff, tol, eExp, prec)
				}
			}
		}
	})

	// -------------------------
	// 6. Shortest mode: no trailing zeros
	// -------------------------
	t.Run("shortest_no_trailing_zeros", func(t *testing.T) {
		for _, val := range []string{
			"0", "1", "12345", "3.14", "0.00123",
			"0.000123", "1e-12", "1e-13", "1e100", "1e-100",
		} {
			f.SetString(val)
			nStr := f.Text('n', -1)
			eIdx := strings.IndexAny(nStr, "eE")
			if eIdx < 0 {
				if nStr == "+Inf" || nStr == "-Inf" || nStr == "Inf" {
					continue
				}
				t.Errorf("shortest %%n output for %q = %q has no exponent", val, nStr)
				continue
			}
			mant := nStr[:eIdx]
			dotIdx := strings.IndexByte(mant, '.')
			if dotIdx >= 0 {
				fracPart := mant[dotIdx+1:]
				if len(fracPart) > 0 && fracPart[len(fracPart)-1] == '0' {
					t.Errorf("shortest %%n output for %q = %q has trailing zero in fraction", val, nStr)
				}
			}
		}
	})

	// -----------------------------------------
	// 7. Format flag integration with fmt.Sprintf
	// -----------------------------------------
	t.Run("format_flags", func(t *testing.T) {
		f.SetString("12345")

		// %+n — always show sign
		got := fmt.Sprintf("%+n", f)
		if !strings.HasPrefix(got, "+") && !strings.HasPrefix(got, "-") {
			t.Errorf("%%+n should include sign, got %q", got)
		}

		// % n — space for sign
		got2 := fmt.Sprintf("% n", f)
		if !strings.HasPrefix(got2, " ") && !strings.HasPrefix(got2, "-") {
			t.Errorf("%% n should include space sign, got %q", got2)
		}

		// %-12.4n — left justification
		got3 := fmt.Sprintf("%-12.4n", f)
		if len(got3) < 12 {
			t.Errorf("%%-12.4n should pad to width 12, got %q (len=%d)", got3, len(got3))
		}

		// %+N — uppercase E
		got4 := fmt.Sprintf("%+N", f)
		if !strings.Contains(got4, "E") {
			t.Errorf("%%+N should use uppercase E, got %q", got4)
		}
	})

	// -----------------------------------------
	// 8. Negative values preserve sign
	// -----------------------------------------
	t.Run("negative_sign", func(t *testing.T) {
		for _, val := range []string{"-12345", "-0.000123", "-1e-13"} {
			f.SetString(val)
			nStr := f.Text('n', 6)
			if nStr[0] != '-' {
				t.Errorf("Text('n',6) for %q should start with '-', got %q", val, nStr)
			}
			eStr := f.Text('e', 6)
			if eStr[0] != '-' {
				t.Errorf("Text('e',6) for %q should start with '-', got %q", val, eStr)
			}
		}
	})

	// -----------------------------------------
	// 9. Unknown format verb fallback preserved
	// -----------------------------------------
	t.Run("unknown_verb", func(t *testing.T) {
		f.SetString("12345")
		got := f.Text('z', 6)
		// Text returns "%" + format byte for unknown formats
		if got != "%z" {
			t.Errorf("unknown verb should produce '%%z', got %q", got)
		}
	})

	// -----------------------------------------
	// 10. Large exponent path (initHuge)
	// -----------------------------------------
	t.Run("initHuge", func(t *testing.T) {
		f.SetPrec(128)
		for _, val := range []string{"1e10000", "1e-10000"} {
			f.SetString(val)
			nStr := f.Text('n', 6)
			eStr := f.Text('e', 6)
			t.Logf("%%e(%s)=%q  %%n(%s)=%q", val, eStr, val, nStr)
			engineeringNotationVerifier(t, nStr, false)
		}
	})
}
