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
