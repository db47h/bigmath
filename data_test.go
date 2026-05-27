// SPDX-License-Identifier: MIT

//go:generate python testdata/gen_go_tests.py testdata/data.txt -o testdata/data_tests.json -p 128

package bigmath_test

import (
	"encoding/json"
	"fmt"
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
		panic(fmt.Sprintf("parseHex(%q, %d): %v", s, prec, err))
	}
	return x
}

// fnMap maps function names to actual bigmath functions.
var fnMap = map[string]any{
	"exp":       bigmath.Exp,
	"log":       bigmath.Log,
	"pow":       bigmath.Pow,
	"atan":      bigmath.Atan,
	"atan2":     bigmath.Atan2,
	"const_pi":  testPiConst,
	"sin":       bigmath.Sin,
	"cos":       bigmath.Cos,
	"sin_cos":   bigmath.Sincos,
	"sinh":      bigmath.Sinh,
	"cosh":      bigmath.Cosh,
	"sinh_cosh": bigmath.SinhCosh,
	"tanh":      bigmath.Tanh,
	"asinh":     bigmath.Asinh,
	"acosh":     bigmath.Acosh,
	"atanh":     bigmath.Atanh,
	"asin":      bigmath.Asin,
	"acos":      bigmath.Acos,
	"tan":       bigmath.Tan,
}

func testPiConst(z *bigmath.Float) *bigmath.Float {
	return bigmath.Pi(z)
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

// TestFloatData runs golden comparison tests.
// Only processes cases where panics is false (or unset).
func TestFloatData(t *testing.T) {
	data, err := loadFloatTestData("testdata/data_tests.json")
	if err != nil {
		t.Fatalf("Failed to load test data: %v", err)
	}

	got := new(bigmath.Float).SetPrec(data.Prec)
	for _, d := range data.Cases {
		if d.Panics {
			continue // skip panic tests
		}
		t.Run(fmt.Sprintf("%s%v", d.Fn, d.Args), func(t *testing.T) {
			fn := fnMap[strings.ToLower(d.Fn)]
			if fn == nil {
				t.Fatalf("unknown function %v", d.Fn)
			}

			if d.Res2 != "" {
				// Tuple function: two results
				got1 := new(bigmath.Float).SetPrec(data.Prec)
				got2 := new(bigmath.Float).SetPrec(data.Prec)
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

// TestFloatPanics tests that certain inputs panic with ErrNaN.
// Only processes cases where panics is true.
func TestFloatPanics(t *testing.T) {
	data, err := loadFloatTestData("testdata/data_tests.json")
	if err != nil {
		t.Fatalf("Failed to load test data: %v", err)
	}

	for _, d := range data.Cases {
		if !d.Panics {
			continue // skip golden tests
		}
		t.Run(fmt.Sprintf("%s%v", d.Fn, d.Args), func(t *testing.T) {
			defer func() {
				r := recover()
				if r == nil {
					t.Errorf("%s(%v): expected ErrNaN panic, got none",
						d.Fn, d.Args)
					return
				}
				if _, ok := r.(bigmath.ErrNaN); ok {
					return // expected type
				}
				// ErrNaN satisfies error with Error() == "ErrNaN"
				if errStr, ok := r.(error); ok && errStr.Error() == "ErrNaN" {
					return
				}
				t.Errorf("%s(%v): expected ErrNaN, got %T(%v)",
					d.Fn, d.Args, r, r)
			}()

			if d.Res2 != "" {
				got1 := new(bigmath.Float).SetPrec(data.Prec)
				got2 := new(bigmath.Float).SetPrec(data.Prec)
				args := makeReflectArgs([]*bigmath.Float{got1, got2}, d.Args, data.Prec)
				reflect.ValueOf(fnMap[strings.ToLower(d.Fn)]).Call(args)
			} else {
				got := new(bigmath.Float).SetPrec(data.Prec)
				args := buildReflectArgs(got, d.Args, data.Prec)
				reflect.ValueOf(fnMap[strings.ToLower(d.Fn)]).Call(args)
			}
		})
	}
}
