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
func parseHex(s string, prec uint) *big.Float {
	x, _, err := new(big.Float).SetPrec(prec).Parse(s, 0)
	if err != nil {
		panic(fmt.Sprintf("parseHex(%q, %d): %v", s, prec, err))
	}
	return x
}

// fnMap maps function names to actual bigmath functions.
var fnMap = map[string]any{
	"exp":      bigmath.Exp,
	"log":      bigmath.Log,
	"pow":      bigmath.Pow,
	"atan":     bigmath.Atan,
	"atan2":    bigmath.Atan2,
	"const_pi": testPiConst,
	"sin":      bigmath.Sin,
	"cos":      bigmath.Cos,
	"sinh":     bigmath.Sinh,
	"cosh":     bigmath.Cosh,
	"tanh":     bigmath.Tanh,
	"asinh":    bigmath.Asinh,
	"acosh":    bigmath.Acosh,
	"atanh":    bigmath.Atanh,
	"asin":     bigmath.Asin,
	"acos":     bigmath.Acos,
	"tan":      bigmath.Tan,
}

func testPiConst(z *big.Float) *big.Float {
	return bigmath.Pi(z)
}

// buildReflectArgs prepares the reflect.Value arguments for a function call.
// First arg is the result *big.Float, remaining args are parsed from strings.
func buildReflectArgs(got *big.Float, strArgs []string, prec uint) []reflect.Value {
	args := []reflect.Value{reflect.ValueOf(got)}
	for _, arg := range strArgs {
		x := parseHex(arg, prec)
		args = append(args, reflect.ValueOf(x))
	}
	return args
}

// TestFloatData runs golden comparison tests.
// Only processes cases where panics is false (or unset).
func TestFloatData(t *testing.T) {
	data, err := loadFloatTestData("testdata/data_tests.json")
	if err != nil {
		t.Fatalf("Failed to load test data: %v", err)
	}

	got := new(big.Float).SetPrec(data.Prec)
	for _, d := range data.Cases {
		if d.Panics {
			continue // skip panic tests
		}
		t.Run(fmt.Sprintf("%s%v", d.Fn, d.Args), func(t *testing.T) {
			fn := fnMap[strings.ToLower(d.Fn)]
			if fn == nil {
				t.Fatalf("unknown function %v", d.Fn)
			}

			// Build reflect args
			args := buildReflectArgs(got, d.Args, data.Prec)
			reflect.ValueOf(fn).Call(args)

			// Parse expected result
			want := parseHex(d.Res, data.Prec)
			if got.Cmp(want) != 0 {
				t.Fatalf("%s(%v): got %s, want %s",
					d.Fn, d.Args, got.Text('x', -1), want.Text('x', -1))
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

			got := new(big.Float).SetPrec(data.Prec)
			args := buildReflectArgs(got, d.Args, data.Prec)
			reflect.ValueOf(fnMap[strings.ToLower(d.Fn)]).Call(args)
		})
	}
}
