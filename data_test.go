// SPDX-License-Identifier: MIT

//go:generate python testdata/gen_go_tests.py testdata/data.txt -o testdata/data_tests.json -p 128

package bigmath_test

import (
	"encoding/json"
	"fmt"
	"math/big"
	"os"
	"reflect"
	"strings"
	"testing"

	"github.com/db47h/bigmath"
)

// floatTestCase represents one test case from data_tests.json.
type floatTestCase struct {
	Fn     string   `json:"fn"`
	Args   []string `json:"args"`
	Res    string   `json:"res,omitempty"`
	Res2   string   `json:"res2,omitempty"`
	Panics bool     `json:"panics,omitempty"`
	Reason string   `json:"reason,omitempty"`
}

// floatTestData wraps the top-level JSON structure.
type floatTestData struct {
	Prec  uint            `json:"prec"`
	Cases []floatTestCase `json:"cases"`
}

// loadFloatTestData reads and parses the JSON test file.
func loadFloatTestData(path string) (*floatTestData, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	var data floatTestData
	if err := json.NewDecoder(f).Decode(&data); err != nil {
		return nil, err
	}
	return &data, nil
}

// parseHex parses a hex float string at the given precision.
func parseHex(s string, prec uint) *bigmath.Float {
	x, _, err := new(bigmath.Float).SetPrec(prec).Parse(s, 0)
	if err != nil {
		panic(fmt.Sprintf("parseHex(%q, %d): %v (this is usually means that gmpy2 returned a NaN)", s, prec, err))
	}
	return x
}

// fnMap maps function names to actual bigmath functions.
var fnMap = map[string]any{
	"cbrt":      (*bigmath.Float).Cbrt,
	"ceil":      (*bigmath.Float).Ceil,
	"exp":       (*bigmath.Float).Exp,
	"log":       (*bigmath.Float).Log,
	"log10":     (*bigmath.Float).Log10,
	"log2":      (*bigmath.Float).Log2,
	"pow":       (*bigmath.Float).Pow,
	"atan":      (*bigmath.Float).Atan,
	"atan2":     (*bigmath.Float).Atan2,
	"const_pi":  testPiConst,
	"sin":       (*bigmath.Float).Sin,
	"cos":       (*bigmath.Float).Cos,
	"sin_cos":   bigmath.Sincos,
	"sinh":      (*bigmath.Float).Sinh,
	"cosh":      (*bigmath.Float).Cosh,
	"sinh_cosh": bigmath.SinhCosh,
	"tanh":      (*bigmath.Float).Tanh,
	"asinh":     (*bigmath.Float).Asinh,
	"acosh":     (*bigmath.Float).Acosh,
	"atanh":     (*bigmath.Float).Atanh,
	"asin":      (*bigmath.Float).Asin,
	"acos":      (*bigmath.Float).Acos,
	"tan":       (*bigmath.Float).Tan,
	"cot":       (*bigmath.Float).Cot,
	"sec":       (*bigmath.Float).Sec,
	"csc":       (*bigmath.Float).Csc,
	"coth":      (*bigmath.Float).Coth,
	"sech":      (*bigmath.Float).Sech,
	"csch":      (*bigmath.Float).Csch,
	"hypot":     (*bigmath.Float).Hypot,
	"floor":     (*bigmath.Float).Floor,
	"fma":       (*bigmath.Float).FMA,
}

func testPiConst(z *bigmath.Float) *bigmath.Float {
	return z.Pi()
}

// makeReflectArgs builds a reflect.Value slice for calling fn.
// gots contains one *big.Float per result pointer (1 for single-result,
// 2 for tuple). strArgs are parsed at prec and appended after the results.
func makeReflectArgs(gots []*bigmath.Float, strArgs []string, prec uint) []reflect.Value {
	args := make([]reflect.Value, len(gots))
	for i, g := range gots {
		args[i] = reflect.ValueOf(g)
	}
	for _, arg := range strArgs {
		x := parseHex(arg, prec)
		args = append(args, reflect.ValueOf(x))
	}
	return args
}

// buildReflectArgs prepares the reflect.Value arguments for a single-result function call.
// First arg is the result *big.Float, remaining args are parsed from strings.
func buildReflectArgs(got *bigmath.Float, strArgs []string, prec uint) []reflect.Value {
	return makeReflectArgs([]*bigmath.Float{got}, strArgs, prec)
}

// TestFloatData runs golden comparison tests for float functions.
// Handles both success cases (panics=false) and expected panic cases (panics=true).
func TestFloatData(t *testing.T) {
	data, err := loadFloatTestData("testdata/data_tests.json")
	if err != nil {
		t.Fatalf("Failed to load test data: %v", err)
	}

	for _, d := range data.Cases {
		t.Run(fmt.Sprintf("%s%v", d.Fn, d.Args), func(t *testing.T) {
			defer func() {
				r := recover()
				if r == nil {
					if d.Panics {
						t.Errorf("%s(%v): expected ErrNaN panic, got none",
							d.Fn, d.Args)
					}
					return
				}
				e, ok := r.(error)
				if !ok {
					panic(r)
				}
				if !d.Panics {
					t.Errorf("%s(%v): expected %s, got error: %v",
						d.Fn, d.Args, d.Res, e)
					return
				}
				switch e.(type) {
				case bigmath.ErrNaN:
				case big.ErrNaN:
				default:
					t.Errorf("%s(%v): expected ErrNaN, got: %T %v",
						d.Fn, d.Args, e, e)
				}
			}()

			fn := fnMap[strings.ToLower(d.Fn)]
			if fn == nil {
				t.Fatalf("unknown function %v", d.Fn)
			}

			if d.Res2 != "" {
				// Tuple function: two results
				got1 := new(bigmath.Float)
				got2 := new(bigmath.Float)
				args := makeReflectArgs([]*bigmath.Float{got1, got2}, d.Args, data.Prec)
				reflect.ValueOf(fn).Call(args)

				want1 := parseHex(d.Res, data.Prec)
				if got1.Cmp(want1) != 0 {
					t.Fatalf("%s(%v)[0]: got %s, want %s",
						d.Fn, d.Args, got1.Text('x', -1), want1.Text('x', -1))
				}
				want2 := parseHex(d.Res2, data.Prec)
				if got2.Cmp(want2) != 0 {
					t.Fatalf("%s(%v)[1]: got %s, want %s",
						d.Fn, d.Args, got2.Text('x', -1), want2.Text('x', -1))
				}
			} else {
				got := new(bigmath.Float)
				if d.Fn == "const_pi" {
					got.SetPrec(data.Prec)
				}
				// Single-result function (existing logic)
				args := buildReflectArgs(got, d.Args, data.Prec)
				reflect.ValueOf(fn).Call(args)

				// Parse expected result
				want := parseHex(d.Res, data.Prec)
				if got.Cmp(want) != 0 {
					t.Fatalf("%s(%v): got %s, want %s",
						d.Fn, d.Args, got.Text('x', -1), want.Text('x', -1))
				}
			}
		})
	}
}
