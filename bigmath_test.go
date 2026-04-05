// SPDX-License-Identifier: MIT

//go:generate python testdata/gen_go_tests.py testdata/data.txt -o data_test.go -p 2048

package bigmath_test

import (
	"fmt"
	"math/big"
	"reflect"
	"strings"
	"testing"

	"github.com/db47h/bigmath"
)

var fnMap = map[string]any{
	"exp": bigmath.Exp,
	"log": bigmath.Log,
	"pow": bigmath.Pow,
}

// func testExp(z *big.Float, args ...big.Float) *big.Float {
// 	return nil
// }

func TestBigMath(t *testing.T) {
	// Create a big.Float with the same precision used in Python
	got := new(big.Float).SetPrec(dataPrec)
	for _, d := range data {
		t.Run(fmt.Sprintf("%s%v", d.fn, d.args), func(t *testing.T) {
			// lookup function
			fn := fnMap[strings.ToLower(d.fn)]
			if fn == nil {
				t.Fatalf("unknown function %v", d.fn)
			}
			fnv := reflect.ValueOf(fn)
			if fnv.Kind() != reflect.Func {
				t.Fatalf("invalid function mapping for %v, not a function", fn)
			}
			// parse args
			var args []reflect.Value
			args = append(args, reflect.ValueOf(got))
			for _, arg := range d.args {
				x, _, err := new(big.Float).SetPrec(dataPrec).Parse(arg, 0)
				if err != nil {
					t.Fatalf("Failed to parse argument: %v", err)
				}
				args = append(args, reflect.ValueOf(x))
			}
			fnv.Call(args)
			// Create the reference from the Python string
			want, _, err := new(big.Float).SetPrec(dataPrec).Parse(d.res, 0)
			if err != nil {
				t.Fatalf("Failed to parse reference result: %v", err)
			}

			if got.Cmp(want) != 0 {
				t.Fatalf("%s(%v): got %s, want %s", d.fn, d.args, got.Text('x', -1), want.Text('x', -1))
			}
		})
	}
}

func TestLogErrNaN(t *testing.T) {
	tests := []string{"-1", "-10", "-inf"}
	for _, tt := range tests {
		t.Run(fmt.Sprintf("Log(%s)", tt), func(t *testing.T) {
			defer func() {
				if r := recover(); r == nil {
					t.Errorf("Log(%s) did not panic", tt)
				} else if _, ok := r.(bigmath.ErrNaN); !ok {
					t.Errorf("Log(%s) panicked with %v, want ErrNaN", tt, r)
				}
			}()
			x, _, _ := new(big.Float).SetPrec(dataPrec).Parse(tt, 0)
			bigmath.Log(new(big.Float).SetPrec(dataPrec), x)
		})
	}
}

func TestPowErrNaN(t *testing.T) {
	tests := []struct {
		x, y string
	}{
		{"-2", "0.5"},
		{"-1", "1.5"},
	}
	for _, tt := range tests {
		t.Run(fmt.Sprintf("Pow(%s,%s)", tt.x, tt.y), func(t *testing.T) {
			defer func() {
				if r := recover(); r == nil {
					t.Errorf("Pow(%s, %s) did not panic", tt.x, tt.y)
				} else if _, ok := r.(bigmath.ErrNaN); !ok {
					t.Errorf("Pow(%s, %s) panicked with %v, want ErrNaN", tt.x, tt.y, r)
				}
			}()
			x, _, _ := new(big.Float).SetPrec(dataPrec).Parse(tt.x, 0)
			y, _, _ := new(big.Float).SetPrec(dataPrec).Parse(tt.y, 0)
			bigmath.Pow(new(big.Float).SetPrec(dataPrec), x, y)
		})
	}
}
