// Copyright 2014 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE-go file.

// This file contains duplicates of all big.Float text formatting functions that
// end up calling Float.Append.

package bigmath

import (
	"fmt"
)

// String formats x like x.Text('g', 10).
// (String must be called explicitly, [Float.Format] does not support %s verb.)
func (x *Float) String() string {
	return x.Text('g', 10)
}

// Text converts the floating-point number x to a string according
// to the given format and precision prec. The format is one of:
//
//	'e'	-d.dddde±dd, decimal exponent, at least two (possibly 0) exponent digits
//	'E'	-d.ddddE±dd, decimal exponent, at least two (possibly 0) exponent digits
//	'f'	-ddddd.dddd, no exponent
//	'g'	like 'e' for large exponents, like 'f' otherwise
//	'G'	like 'E' for large exponents, like 'f' otherwise
//	'x'	-0xd.dddddp±dd, hexadecimal mantissa, decimal power of two exponent
//	'p'	-0x.dddp±dd, hexadecimal mantissa, decimal power of two exponent (non-standard)
//	'b'	-ddddddp±dd, decimal mantissa, decimal power of two exponent (non-standard)
//
// For the power-of-two exponent formats, the mantissa is printed in normalized form:
//
//	'x'	hexadecimal mantissa in [1, 2), or 0
//	'p'	hexadecimal mantissa in [½, 1), or 0
//	'b'	decimal integer mantissa using x.Prec() bits, or 0
//
// Note that the 'x' form is the one used by most other languages and libraries.
//
// If format is a different character, Text returns a "%" followed by the
// unrecognized format character.
//
// The precision prec controls the number of digits (excluding the exponent)
// printed by the 'e', 'E', 'f', 'g', 'G', and 'x' formats.
// For 'e', 'E', 'f', and 'x', it is the number of digits after the decimal point.
// For 'g' and 'G' it is the total number of digits. A negative precision selects
// the smallest number of decimal digits necessary to identify the value x uniquely
// using x.Prec() mantissa bits.
// The prec value is ignored for the 'b' and 'p' formats.
func (x *Float) Text(format byte, prec int) string {
	buf := x.Append(nil, format, prec)
	return string(buf)
}

// Format implements the [fmt.Formatter] interface.
// It supports the verbs 'b', 'e', 'E', 'f', 'F', 'g', 'G', 'p', 'v',
// and 'x' for floating-point numbers. 'v' is handled like 'g',
// and 'F' is handled like 'f'. The precision specifier sets the number
// of digits, the output field width, as well as the format flags
// '+' and ' ' for sign control, '0' for space or zero padding,
// and '-' for left or right justification. See the fmt package
// for details.
func (x *Float) Format(s fmt.State, verb rune) {
	prec, hasPrec := s.Precision()
	if !hasPrec {
		prec = 6 // default precision for 'e', 'f'
	}

	switch verb {
	case 'e', 'E', 'f', 'b', 'p', 'x':
		// nothing to do
	case 'F':
		// (*Float).Text doesn't support 'F'; handle like 'f'
		verb = 'f'
	case 'v':
		// handle like 'g'
		verb = 'g'
		fallthrough
	case 'g', 'G':
		if !hasPrec {
			prec = -1 // default precision for 'g', 'G'
		}
	default:
		fmt.Fprintf(s, "%%!%c(*big.Float=%s)", verb, x.String())
		return
	}
	var buf []byte
	buf = x.Append(buf, byte(verb), prec)
	if len(buf) == 0 {
		buf = []byte("?") // should never happen, but don't crash
	}
	// len(buf) > 0

	var sign string
	switch {
	case buf[0] == '-':
		sign = "-"
		buf = buf[1:]
	case buf[0] == '+':
		// +Inf
		sign = "+"
		if s.Flag(' ') {
			sign = " "
		}
		buf = buf[1:]
	case s.Flag('+'):
		sign = "+"
	case s.Flag(' '):
		sign = " "
	}

	var padding int
	if width, hasWidth := s.Width(); hasWidth && width > len(sign)+len(buf) {
		padding = width - len(sign) - len(buf)
	}

	switch {
	case s.Flag('0') && !x.IsInf():
		// 0-padding on left
		writeMultiple(s, sign, 1)
		writeMultiple(s, "0", padding)
		s.Write(buf)
	case s.Flag('-'):
		// padding on right
		writeMultiple(s, sign, 1)
		s.Write(buf)
		writeMultiple(s, " ", padding)
	default:
		// padding on left
		writeMultiple(s, " ", padding)
		writeMultiple(s, sign, 1)
		s.Write(buf)
	}
}

// AppendText implements the [encoding.TextAppender] interface.
// Only the [Float] value is marshaled (in full precision), other
// attributes such as precision or accuracy are ignored.
func (x *Float) AppendText(b []byte) ([]byte, error) {
	if x == nil {
		return append(b, "<nil>"...), nil
	}
	return x.Append(b, 'g', -1), nil
}

// MarshalText implements the [encoding.TextMarshaler] interface.
// Only the [Float] value is marshaled (in full precision), other
// attributes such as precision or accuracy are ignored.
func (x *Float) MarshalText() (text []byte, err error) {
	return x.AppendText(nil)
}

// writeMultiple writes s count times.
func writeMultiple(s fmt.State, str string, count int) {
	for range count {
		fmt.Fprint(s, str)
	}
}
