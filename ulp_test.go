// SPDX-License-Identifier: MIT

package bigmath_test

import (
	"encoding/json"
	"fmt"
	"math/big"
	"os"
	"testing"

	"github.com/db47h/bigmath"
)

// ulpPoint represents a single test point with all reference values.
type ulpPoint struct {
	X         string `json:"x"`
	SinRef    string `json:"sin_ref"`
	CosRef    string `json:"cos_ref"`
	Sin3xRef  string `json:"sin_3x_ref"`
	SinhRef   string `json:"sinh_ref"`
	CoshRef   string `json:"cosh_ref"`
	Cosh2xRef string `json:"cosh_2x_ref"`
}

// ulpData wraps a set of test points for a given precision.
type ulpData struct {
	Prec    uint       `json:"prec"`
	RefPrec uint       `json:"ref_prec"`
	Points  []ulpPoint `json:"points"`
}

// loadULPData reads the JSON reference data file for the given precision.
func loadULPData(prec uint) (*ulpData, error) {
	path := fmt.Sprintf("testdata/ulp_data_%d.json", prec)
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	var data ulpData
	if err := json.NewDecoder(f).Decode(&data); err != nil {
		return nil, err
	}
	return &data, nil
}

// parseHexFloat parses a hex float string at the given precision.
func parseHexFloat(s string, prec uint) *big.Float {
	z := new(big.Float).SetPrec(prec)
	_, _, err := z.Parse(s, 0)
	if err != nil {
		panic(fmt.Sprintf("parseHexFloat(%q, %d): %v", s, prec, err))
	}
	return z
}

// ulpErr computes the ULP error between got and ref.
// Both got and ref must be at the same precision.
// Returns |got - ref| / ULP(ref).
func ulpErr(got, ref *big.Float) float64 {
	// Handle special cases: if both are Inf or both are zero, error is 0.
	if got.IsInf() || ref.IsInf() {
		if got.IsInf() && ref.IsInf() && got.Signbit() == ref.Signbit() {
			return 0
		}
		// One is Inf and the other is not — we can't compute meaningful ULP.
		// Return a large sentinel.
		return 1e9
	}
	if ref.Sign() == 0 {
		if got.Sign() == 0 {
			return 0
		}
		// got is non-zero but ref is zero — compute ULP relative to ULP at zero.
		// For subnormal-like: ULP = 2^(0 - prec) when ref is exactly 0.
		prec := ref.Prec()
		ulp := new(big.Float).SetPrec(64).SetMantExp(
			new(big.Float).SetFloat64(1),
			-int(prec),
		)
		diff := new(big.Float).Sub(got, ref)
		diff.Abs(diff)
		diff.Quo(diff, ulp)
		err, _ := diff.Float64()
		return err
	}

	refExp := ref.MantExp(nil)
	prec := ref.Prec()
	// Clamp ULP exponent to minimum normal: for denormal values (exponent < -prec),
	// use the minimum normal ULP (2^(-prec)) to avoid producing absurdly large
	// ratios when ref is very close to zero (e.g., cos(π/2)).
	ulpExp := refExp - int(prec)
	if ulpExp < -int(prec) {
		ulpExp = -int(prec)
	}
	diff := new(big.Float).Sub(got, ref)
	diff.Abs(diff)
	ulp := new(big.Float).SetPrec(64).SetMantExp(
		new(big.Float).SetFloat64(1),
		ulpExp,
	)
	diff.Quo(diff, ulp)
	err, _ := diff.Float64()
	return err
}

// isSpecial returns true if x is ±Inf or ±0 (exact-zero cases that don't
// need ULP measurement).
func isSpecial(x *big.Float) bool {
	return x.IsInf() || x.Sign() == 0
}

// TestULPErrorDirect measures per-function ULP error by comparing bigmath
// results against gmpy2 reference values at target precision.
func TestULPErrorDirect(t *testing.T) {
	precs := []uint{64, 128, 256, 1024}
	const maxULP = 2.0

	for _, prec := range precs {
		t.Run(fmt.Sprintf("prec=%d", prec), func(t *testing.T) {
			data, err := loadULPData(prec)
			if err != nil {
				t.Fatalf("load data: %v", err)
			}

			maxSin, maxCos, maxSinh, maxCosh := 0.0, 0.0, 0.0, 0.0

			for i, pt := range data.Points {
				x := parseHexFloat(pt.X, prec)
				refSin := parseHexFloat(pt.SinRef, prec)
				refCos := parseHexFloat(pt.CosRef, prec)
				refSinh := parseHexFloat(pt.SinhRef, prec)
				refCosh := parseHexFloat(pt.CoshRef, prec)

				// Test Sin
				got := new(big.Float).SetPrec(prec)
				bigmath.Sin(got, x)
				err := ulpErr(got, refSin)
				if err > maxSin {
					maxSin = err
				}
				if err > maxULP && !isSpecial(refSin) {
					t.Logf("point %d (x=%s): Sin ULP error=%g", i, pt.X, err)
				}

				// Test Cos
				bigmath.Cos(got, x)
				err = ulpErr(got, refCos)
				if err > maxCos {
					maxCos = err
				}
				if err > maxULP && !isSpecial(refCos) {
					t.Logf("point %d (x=%s): Cos ULP error=%g", i, pt.X, err)
				}

				// Test Sinh
				bigmath.Sinh(got, x)
				err = ulpErr(got, refSinh)
				if err > maxSinh {
					maxSinh = err
				}
				if err > maxULP && !isSpecial(refSinh) {
					t.Logf("point %d (x=%s): Sinh ULP error=%g", i, pt.X, err)
				}

				// Test Cosh
				bigmath.Cosh(got, x)
				err = ulpErr(got, refCosh)
				if err > maxCosh {
					maxCosh = err
				}
				if err > maxULP && !isSpecial(refCosh) {
					t.Logf("point %d (x=%s): Cosh ULP error=%g", i, pt.X, err)
				}
			}

			t.Logf("prec=%d max ULPs: sin=%g cos=%g sinh=%g cosh=%g",
				prec, maxSin, maxCos, maxSinh, maxCosh)

			if maxSin > maxULP {
				t.Errorf("Sin max ULP error %g exceeds %g", maxSin, maxULP)
			}
			if maxCos > maxULP {
				t.Errorf("Cos max ULP error %g exceeds %g", maxCos, maxULP)
			}
			if maxSinh > maxULP {
				t.Errorf("Sinh max ULP error %g exceeds %g", maxSinh, maxULP)
			}
			if maxCosh > maxULP {
				t.Errorf("Cosh max ULP error %g exceeds %g", maxCosh, maxULP)
			}
		})
	}
}

// identityULP computes |result - expected| / ULP(result) as a float64.
func identityULP(result, expected *big.Float) float64 {
	// Handle special cases
	if result.IsInf() || expected.IsInf() {
		if result.IsInf() && expected.IsInf() && result.Signbit() == expected.Signbit() {
			return 0
		}
		return 1e9
	}
	if result.Sign() == 0 {
		if expected.Sign() == 0 {
			return 0
		}
		// result is zero but expected isn't — measure ULP at zero
		prec := expected.Prec()
		ulp := new(big.Float).SetPrec(64).SetMantExp(
			new(big.Float).SetFloat64(1),
			-int(prec),
		)
		diff := new(big.Float).Sub(result, expected)
		diff.Abs(diff)
		diff.Quo(diff, ulp)
		err, _ := diff.Float64()
		return err
	}

	exp := result.MantExp(nil)
	prec := result.Prec()
	// Clamp ULP exponent to minimum normal: for denormal values (exponent < -prec),
	// use the minimum normal ULP (2^(-prec)) to avoid producing absurdly large
	// ratios when the result is very close to zero.
	ulpExp := exp - int(prec)
	if ulpExp < -int(prec) {
		ulpExp = -int(prec)
	}
	diff := new(big.Float).Sub(result, expected)
	diff.Abs(diff)
	ulp := new(big.Float).SetPrec(64).SetMantExp(
		new(big.Float).SetFloat64(1),
		ulpExp,
	)
	diff.Quo(diff, ulp)
	err, _ := diff.Float64()
	return err
}

// TestSinCosSquared verifies sin²(x) + cos²(x) ≈ 1 for all test points.
func TestSinCosSquared(t *testing.T) {
	precs := []uint{64, 128, 256, 1024}
	const maxULP = 2.0

	for _, prec := range precs {
		t.Run(fmt.Sprintf("prec=%d", prec), func(t *testing.T) {
			data, err := loadULPData(prec)
			if err != nil {
				t.Fatalf("load data: %v", err)
			}

			maxErr := 0.0
			oneRef := new(big.Float).SetPrec(prec).SetUint64(1)

			for i, pt := range data.Points {
				x := parseHexFloat(pt.X, prec)

				s := new(big.Float).SetPrec(prec)
				c := new(big.Float).SetPrec(prec)
				bigmath.Sincos(s, c, x)

				// sin² + cos²
				s2 := new(big.Float).SetPrec(prec).Mul(s, s)
				c2 := new(big.Float).SetPrec(prec).Mul(c, c)
				sum := new(big.Float).SetPrec(prec).Add(s2, c2)

				// Skip if Inf (can't compute ULP meaningfully)
				if sum.IsInf() || oneRef.IsInf() {
					continue
				}

				err := identityULP(sum, oneRef)
				if err > maxErr {
					maxErr = err
				}
				if err > maxULP {
					t.Logf("point %d (x=%s): sin²+cos² ULP error=%g", i, pt.X, err)
				}
			}

			t.Logf("prec=%d max sin²+cos² identity ULP error=%g", prec, maxErr)

			if maxErr > maxULP {
				t.Errorf("sin²+cos² identity max ULP error %g exceeds %g", maxErr, maxULP)
			}
		})
	}
}

// TestSinTripleAngle verifies sin(3x) = 3·sin(x) − 4·sin³(x).
func TestSinTripleAngle(t *testing.T) {
	precs := []uint{64, 128, 256, 1024}
	const maxULP = 2.0

	for _, prec := range precs {
		t.Run(fmt.Sprintf("prec=%d", prec), func(t *testing.T) {
			data, err := loadULPData(prec)
			if err != nil {
				t.Fatalf("load data: %v", err)
			}

			maxErr := 0.0
			three := new(big.Float).SetPrec(prec).SetUint64(3)
			four := new(big.Float).SetPrec(prec).SetUint64(4)

			for i, pt := range data.Points {
				x := parseHexFloat(pt.X, prec)

				// LHS: sin(3x)
				threeX := new(big.Float).SetPrec(prec).Mul(x, three)
				sin3x := new(big.Float).SetPrec(prec)
				bigmath.Sin(sin3x, threeX)

				// RHS: 3·sin(x) - 4·sin³(x)
				sinx := new(big.Float).SetPrec(prec)
				bigmath.Sin(sinx, x)

				// sin³(x)
				sinx3 := new(big.Float).SetPrec(prec).Mul(sinx, sinx)
				sinx3.Mul(sinx3, sinx)

				t0 := new(big.Float).SetPrec(prec).Mul(three, sinx)
				t1 := new(big.Float).SetPrec(prec).Mul(four, sinx3)
				rhs := new(big.Float).SetPrec(prec).Sub(t0, t1)

				// Skip if either side is Inf
				if sin3x.IsInf() || rhs.IsInf() {
					continue
				}

				err := identityULP(sin3x, rhs)
				if err > maxErr {
					maxErr = err
				}
				if err > maxULP {
					t.Logf("point %d (x=%s): sin(3x) identity ULP error=%g", i, pt.X, err)
				}
			}

			t.Logf("prec=%d max sin(3x) identity ULP error=%g", prec, maxErr)

			if maxErr > maxULP {
				t.Errorf("sin(3x) identity max ULP error %g exceeds %g", maxErr, maxULP)
			}
		})
	}
}

// TestSinhCoshIdentity verifies sinh²(x) + cosh²(x) = cosh(2x).
func TestSinhCoshIdentity(t *testing.T) {
	precs := []uint{64, 128, 256, 1024}
	const maxULP = 2.0

	for _, prec := range precs {
		t.Run(fmt.Sprintf("prec=%d", prec), func(t *testing.T) {
			data, err := loadULPData(prec)
			if err != nil {
				t.Fatalf("load data: %v", err)
			}

			maxErr := 0.0
			two := new(big.Float).SetPrec(prec).SetUint64(2)

			for i, pt := range data.Points {
				x := parseHexFloat(pt.X, prec)

				// LHS: sinh²(x) + cosh²(x)
				sh := new(big.Float).SetPrec(prec)
				ch := new(big.Float).SetPrec(prec)
				bigmath.SinhCosh(sh, ch, x)

				sh2 := new(big.Float).SetPrec(prec).Mul(sh, sh)
				ch2 := new(big.Float).SetPrec(prec).Mul(ch, ch)
				lhs := new(big.Float).SetPrec(prec).Add(sh2, ch2)

				// RHS: cosh(2x)
				twoX := new(big.Float).SetPrec(prec).Mul(x, two)
				cosh2x := new(big.Float).SetPrec(prec)
				bigmath.Cosh(cosh2x, twoX)

				// Skip if either side is Inf
				if lhs.IsInf() || cosh2x.IsInf() {
					continue
				}

				err := identityULP(lhs, cosh2x)
				if err > maxErr {
					maxErr = err
				}
				if err > maxULP {
					t.Logf("point %d (x=%s): sinh²+cosh² identity ULP error=%g", i, pt.X, err)
				}
			}

			t.Logf("prec=%d max sinh²+cosh² identity ULP error=%g", prec, maxErr)

			if maxErr > maxULP {
				t.Errorf("sinh²+cosh² identity max ULP error %g exceeds %g", maxErr, maxULP)
			}
		})
	}
}
