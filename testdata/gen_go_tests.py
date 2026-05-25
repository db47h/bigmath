#!/usr/bin/env python3

# SPDX-License-Identifier: MIT

import sys
import json
import argparse
import gmpy2
import builtins


def format_or_nan(val):
    """Format an mpfr value as hex float string.
    Returns (hex_string, is_nan) tuple.
    is_nan is True if the value is NaN."""
    if gmpy2.is_nan(val):
        return None, True
    return format(val, 'a'), False


def generate_json_tests(input_file, output_file, precision):
    ctx = gmpy2.get_context()
    ctx.precision = precision
    # Attempt to match Go's math/big exponent limits as closely as MPFR allows
    # Go: -2^31 to 2^31-1. MPFR is often limited to 2^30-1 or 2^62-1.
    GO_EMAX = 2147483647
    GO_EMIN = -2147483648
    ctx.emax = min(GO_EMAX, gmpy2.get_emax_max())
    ctx.emin = max(GO_EMIN, gmpy2.get_emin_min())

    cases = []

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            func_name = parts[0]
            raw_args = parts[1:]

            # Detect !panic annotation
            is_panic = "!panic" in raw_args
            str_args = [a for a in raw_args if a != "!panic"]

            # For constant functions, ignore passed args (they are placeholders)
            if func_name.startswith("const_"):
                str_args = []

            # Parse arguments
            mpfr_args = [gmpy2.mpfr(arg) for arg in str_args]

            # Lookup function: gmpy2 -> globals -> builtins
            if func_name.startswith("const_"):
                func = getattr(gmpy2, func_name, None)
                if func is None:
                    raise AttributeError(f"constant function '{func_name}' not found in gmpy2")
                result = func()
            else:
                func = getattr(gmpy2, func_name, None)
                if func is None:
                    func = globals().get(func_name)
                if func is None:
                    func = getattr(builtins, func_name, None)
                if func is None:
                    raise AttributeError(f"function '{func_name}' not found")
                result = func(*mpfr_args)

            # Handle tuple results (sin_cos, sinh_cosh)
            if isinstance(result, tuple):
                entry = {"fn": func_name, "args": str_args}
                hex_val0, nan0 = format_or_nan(result[0])
                hex_val1, nan1 = format_or_nan(result[1])
                if nan0 or nan1 or is_panic:
                    entry["panics"] = True
                    if nan0 or nan1:
                        entry["reason"] = "gmpy2 returned NaN"
                    if is_panic and not (nan0 or nan1):
                        print(
                            f"Warning: {func_name} {str_args} tagged !panic "
                            f"but gmpy2 did not produce NaN",
                            file=sys.stderr,
                        )
                else:
                    entry["res"] = hex_val0
                    entry["res2"] = hex_val1
                cases.append(entry)
                continue

            # Single result
            hex_val, is_nan = format_or_nan(result)

            entry = {"fn": func_name, "args": str_args}
            if is_nan or is_panic:
                entry["panics"] = True
                if is_nan:
                    entry["reason"] = "gmpy2 returned NaN"
                if is_panic and not is_nan:
                    print(
                        f"Warning: {func_name} {str_args} tagged !panic "
                        f"but gmpy2 did not produce NaN",
                        file=sys.stderr,
                    )
            else:
                entry["res"] = hex_val

            cases.append(entry)

    output = {"prec": precision, "cases": cases}

    out = open(output_file, 'w') if output_file else sys.stdout
    try:
        json.dump(output, out, indent=2)
        out.write('\n')
    finally:
        if output_file:
            out.close()


def generate_cplx_tests(input_file, output_file, precision):
    """Generate complex golden test data from cplx_data.txt format.

    Each input line:
        <func_name> <re1> <im1> [<re2> <im2> ...] [!panic]

    Parses real/imag pairs into gmpy2.mpc values, dispatches via gmpy2,
    and writes res_re/res_im for non-panic cases.
    """
    ctx = gmpy2.get_context()
    ctx.precision = precision
    GO_EMAX = 2147483647
    GO_EMIN = -2147483648
    ctx.emax = min(GO_EMAX, gmpy2.get_emax_max())
    ctx.emin = max(GO_EMIN, gmpy2.get_emin_min())

    cases = []

    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            func_name = parts[0]
            raw_pairs = parts[1:]

            # Detect !panic annotation
            is_panic = "!panic" in raw_pairs
            tokens = [t for t in raw_pairs if t != "!panic"]

            # Pair-wise parse into mpc values
            # tokens = [re1, im1, re2, im2, ...]
            mpc_args = []
            for i in range(0, len(tokens), 2):
                re = gmpy2.mpfr(tokens[i])
                im = gmpy2.mpfr(tokens[i + 1])
                mpc_args.append(gmpy2.mpc(re, im))

            # Build args pairs as strings for JSON output
            str_args = [[tokens[i], tokens[i + 1]] for i in range(0, len(tokens), 2)]

            # Lookup function and call with mpc_args
            # Some Go function names differ from gmpy2 names, handled here
            if func_name == "quo":
                result = gmpy2.div(*mpc_args)
            elif func_name == "neg":
                result = -mpc_args[0]
            elif func_name == "conj":
                result = mpc_args[0].conjugate()
            elif func_name == "pow":
                result = builtins.pow(*mpc_args)
            elif func_name.startswith("const_"):
                func = getattr(gmpy2, func_name, None)
                if func is None:
                    raise AttributeError(f"constant function '{func_name}' not found in gmpy2")
                result = func()
            else:
                func = getattr(gmpy2, func_name, None)
                if func is None:
                    raise AttributeError(f"complex function '{func_name}' not found in gmpy2")
                result = func(*mpc_args)

            # Extract real and imaginary parts
            hex_re, nan_re = format_or_nan(result.real)
            hex_im, nan_im = format_or_nan(result.imag)

            entry = {"fn": func_name, "args": str_args}

            if nan_re or nan_im or is_panic:
                entry["panics"] = True
                if nan_re or nan_im:
                    entry["reason"] = "gmpy2 returned NaN"
                if is_panic and not (nan_re or nan_im):
                    print(
                        f"Warning: {func_name} {tokens} tagged !panic "
                        f"but gmpy2 did not produce NaN",
                        file=sys.stderr,
                    )
            else:
                entry["res_re"] = hex_re
                entry["res_im"] = hex_im

            cases.append(entry)

    output = {"prec": precision, "mode": "cplx", "cases": cases}

    out = open(output_file, 'w') if output_file else sys.stdout
    try:
        json.dump(output, out, indent=2)
        out.write('\n')
    finally:
        if output_file:
            out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate test data in JSON format.")
    parser.add_argument("input", help="Input text file")
    parser.add_argument("-o", "--output", help="Output JSON file")
    parser.add_argument(
        "-p", "--precision", type=int, default=256,
        help="Precision in bits (default: 256)",
    )
    parser.add_argument(
        "-m", "--mode", choices=["float", "cplx"], default="float",
        help="Test data mode (default: float)",
    )

    args = parser.parse_args()
    if args.mode == "cplx":
        generate_cplx_tests(args.input, args.output, args.precision)
    else:
        generate_json_tests(args.input, args.output, args.precision)
