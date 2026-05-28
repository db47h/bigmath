// SPDX-License-Identifier: MIT

//go:generate python testdata/gen_go_tests.py -m cplx testdata/cplx_data.txt -o testdata/cplx_data_tests.json -p 128

package bigmath_test

import (
	"encoding/json"
	"fmt"
	"math/big"
	"os"
	"testing"

	"github.com/db47h/bigmath"
)

// cplxTestCase represents one complex test case.
type cplxTestCase struct {
	Fn     string      `json:"fn"`
	Args   [][2]string `json:"args"` // each element is [re, im]
	ResRe  string      `json:"res_re,omitempty"`
	ResIm  string      `json:"res_im,omitempty"`
	Panics bool        `json:"panics,omitempty"`
	Reason string      `json:"reason,omitempty"`
}

type cplxTestData struct {
	Prec  uint           `json:"prec"`
	Mode  string         `json:"mode"`
	Cases []cplxTestCase `json:"cases"`
}

// loadCplxTestData reads and parses the complex JSON test file.
func loadCplxTestData(path string) (*cplxTestData, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	var data cplxTestData
	if err := json.NewDecoder(f).Decode(&data); err != nil {
		return nil, err
	}
	return &data, nil
}

// parseCplxHex parses a hex float string at the given precision.
func parseCplxHex(x *bigmath.Float, s string) {
	_, _, err := x.Parse(s, 0)
	if err != nil {
		panic(fmt.Sprintf("parseCplxHex(%q): %v", s, err))
	}
}

// parseComplexArg builds a *bigmath.Complex from a JSON [re,im] pair at prec.
func parseComplexArg(c *bigmath.Complex, pair [2]string) {
	_, _, err := c.Real.Parse(pair[0], 0)
	if err != nil {
		panic(fmt.Sprintf("parseComplexArg real(%q): %v", pair[0], err))
	}
	_, _, err = c.Imag.Parse(pair[1], 0)
	if err != nil {
		panic(fmt.Sprintf("parseComplexArg imag(%q): %v", pair[1], err))
	}
}

// callComplexFn dispatches to the correct *Complex method by function name.
// z is the result pointer, x is the first arg, y is the second (for binary functions).
func callComplexFn(fn string, z, x, y *bigmath.Complex) {
	switch fn {
	case "add":
		z.Add(x, y)
	case "sub":
		z.Sub(x, y)
	case "mul":
		z.Mul(x, y)
	case "quo":
		z.Quo(x, y)
	case "neg":
		z.Neg(x)
	case "conj":
		z.Conj(x)
	case "exp":
		z.Exp(x)
	case "log":
		z.Log(x)
	case "sin":
		z.Sin(x)
	case "cos":
		z.Cos(x)
	case "cot":
		z.Cot(x)
	case "tan":
		z.Tan(x)
	case "sinh":
		z.Sinh(x)
	case "cosh":
		z.Cosh(x)
	case "tanh":
		z.Tanh(x)
	case "asin":
		z.Asin(x)
	case "acos":
		z.Acos(x)
	case "atan":
		z.Atan(x)
	case "asinh":
		z.Asinh(x)
	case "acosh":
		z.Acosh(x)
	case "atanh":
		z.Atanh(x)
	case "sqrt":
		z.Sqrt(x)
	case "pow":
		z.Pow(x, y)
	default:
		panic(fmt.Sprintf("unknown complex function: %s", fn))
	}
}

// isBinaryFn returns true if the function takes two complex arguments.
func isBinaryFn(fn string) bool {
	switch fn {
	case "add", "sub", "mul", "quo", "pow":
		return true
	default:
		return false
	}
}

// TestComplexData runs golden comparison tests for complex functions.
// Only processes cases where panics is false (or unset).
func TestComplexData(t *testing.T) {
	data, err := loadCplxTestData("testdata/cplx_data_tests.json")
	if err != nil {
		t.Fatalf("Failed to load test data: %v", err)
	}

	if data.Mode != "cplx" {
		t.Fatalf("Expected mode 'cplx', got %q", data.Mode)
	}

	x := new(bigmath.Complex).SetPrec(data.Prec)
	y := new(bigmath.Complex).SetPrec(data.Prec)
	got := new(bigmath.Complex).SetPrec(data.Prec)
	wantRe := new(bigmath.Float).SetPrec(data.Prec)
	wantIm := new(bigmath.Float).SetPrec(data.Prec)

	for _, d := range data.Cases {
		t.Run(fmt.Sprintf("%s%v", d.Fn, d.Args), func(t *testing.T) {
			defer func() {
				r := recover()
				if r == nil {
					if d.Panics {
						t.Errorf("%s(%v): expected ErrNaN panic, got none", d.Fn, d.Args)
					}
					return
				}
				e, ok := r.(error)
				if !ok {
					panic(r)
				}
				if !d.Panics {
					t.Errorf("%s(%v): expected (%s, %s), got error: %v", d.Fn, d.Args,
						d.ResRe, d.ResIm, e)
					return
				}
				switch e.(type) {
				case bigmath.ErrNaN:
				case big.ErrNaN:
				default:
					t.Errorf("%s(%v): expected ErrNaN, got: %T %v", d.Fn, d.Args, e, e)
				}
			}()

			// Parse arguments
			parseComplexArg(x, d.Args[0])
			if isBinaryFn(d.Fn) {
				parseComplexArg(y, d.Args[1])
			}

			// Dispatch
			got.Set(x) // test aliasing on x at the same time
			callComplexFn(d.Fn, got, got, y)

			// Parse expected results
			parseCplxHex(wantRe, d.ResRe)
			parseCplxHex(wantIm, d.ResIm)

			if got.Real.Cmp(wantRe) != 0 || got.Imag.Cmp(wantIm) != 0 {
				t.Fatalf("%s(%v): got (%s, %s), want (%s, %s)",
					d.Fn, d.Args,
					got.Real.Text('x', -1), got.Imag.Text('x', -1),
					wantRe.Text('x', -1), wantIm.Text('x', -1))
			}
		})
	}
}
